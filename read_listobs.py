#!/usr/bin/env python
"""
read_listobs.py

Pre-flight sanity check for a radio interferometer measurement set.
Pulls the same information you'd otherwise have to squint at in a
listobs dump and turns it into a short, searchable report: antennas,
frequency setup, field classification, scan timing, current flag
fraction (including per-antenna and per-channel breakdowns to help
spot dead antennas or RFI), and a checklist.

Most of the report is generic CASA/msmetadata output and works on any
MS. A thin per-telescope "profile" supplies the instrument-specific
knowledge: antenna naming convention, receiver band edges, known RFI
ranges, and standard flux calibrator names. Two profiles are included
now: GMRT/uGMRT and VLA/EVLA.

TELESCOPE SELECTION
    The script tries to auto-detect the telescope from the MS's own
    OBSERVATION table. You can override this explicitly:

        casa -c read_listobs.py --telescope gmrt yourfile.ms
        casa -c read_listobs.py --telescope vla   yourfile.ms

    or set the TELESCOPE variable below and use execfile() from the
    CASA prompt.

FITS INPUT (before you've converted to an MS)
    You can also point this at a raw UVFITS file (.fits/.FITS) instead
    of an MS:

        casa -c read_listobs.py --telescope gmrt yourfile.FITS

    This mode reads the FITS file's own header and binary tables (AN,
    FQ, SU) directly, no MS is created, nothing is written to disk
    except the report. It needs the 'astropy' package (pip install
    astropy) since CASA's own msmetadata/table tools only work on an
    MS, not raw FITS.

    IMPORTANT: a raw FITS file doesn't have "scans" or "flags" yet,
    those only get built during MS import (CASA's importgmrt task, the
    same one CAPTURE calls internally). So FITS mode cannot report scan
    timing, current flagging status, elevation during a specific
    field's track, or the scan/elevation-dependent self-cal risk
    factors, the report prints an explicit notice at the top listing
    exactly what's skipped and why. Everything else (antennas,
    frequency/band setup, suggested cell size/FOV/image size, field
    classification and positions, polarization products, the risk
    factors that don't need scan timing) is read directly from the
    FITS file itself.

ADDING A NEW TELESCOPE
    Add one entry to TELESCOPE_PROFILES with the same keys as the
    existing ones. Nothing else in the script needs to change. See
    the comments inside TELESCOPE_PROFILES for what each key does and
    which fields are safe to leave empty/None if you don't have
    verified data yet (don't guess RFI ranges or band edges, an empty
    list/None just means that check is skipped and noted as such).

Output is printed to the terminal and also written to
read_listobs_report.txt in the current directory.
"""

import os
import sys
import re
import numpy as np

# ---------------------------------------------------------------------
# CASA tools/tasks: available as globals inside casashell (casa -c ...).
# Imported explicitly as a fallback so the script also runs under a
# plain modular-CASA python interpreter.
# ---------------------------------------------------------------------
try:
    msmetadata
    table
    flagdata
    quanta
    measures
except NameError:
    from casatools import msmetadata, table, quanta, measures
    from casatasks import flagdata

# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------
MSNAME = None       # set explicitly, or pass the ms path as a script arg
TELESCOPE = None    # 'gmrt' or 'vla', or leave None to auto-detect

# The per-channel RFI breakdown is the one part of this script that reads
# actual visibility data (the FLAG column) rather than lightweight metadata
# tables, so it's by far the slowest part on a large MS (scales with
# scans x baselines x channels x polarizations). Set to False to skip it
# and get a much faster run when you don't need the per-channel view,
# e.g. once you've already checked it and are just re-running for
# everything else. Can also override per-run with --nochanflags.
CHECK_PERCHANNEL_FLAGS = True

# Shared calibrator catalog file. This is the VLA calibrator manual
# list, used both to identify VLA phase calibrators directly and by
# GMRT/uGMRT observers (who borrow the same catalog for cm-wave phase
# calibrator selection). Same file, same role, for either telescope.
PHASE_CAL_FILE = 'vla-cals.list'

# a field named like "0834+555" or "J0834+5533" looks calibrator-shaped
CALLIKE_NAME = re.compile(r'^(J)?\d{4}[+-]\d{2,4}$')

# Per-channel spike detection (below, in the flagging section) marks a
# channel as a statistical "outlier" if it's more than RFI_SPIKE_SIGMA
# standard deviations above that spw's median flagged fraction, with a
# floor of RFI_SPIKE_FLOOR_PCT percentage points so a very flat/quiet spw
# doesn't flag ordinary noise as a spike. This is a heuristic this script
# made up, not sourced from CAPTURE or any external RFI-detection
# standard, tune it if it's missing real RFI or over-flagging. It's only
# used to decide which channels get an "elevated" marker, every channel
# above its spw's median is still printed either way, so nothing that
# looks unusual is hidden behind this cutoff.
RFI_SPIKE_SIGMA = 3.0
RFI_SPIKE_FLOOR_PCT = 10.0


# ---------------------------------------------------------------------
# TELESCOPE PROFILES
#
# match_names     : substrings (case-insensitive) matched against the
#                   MS OBSERVATION table's telescope name, for
#                   auto-detection.
# antenna_pattern : compiled regex an antenna name is expected to
#                   match. Used only to flag oddities, not a hard
#                   validator of which specific antennas exist.
# band_ranges     : list of (label, lo_Hz, hi_Hz) receiver bands, used
#                   to guess which band your data is in.
# rfi_bands       : list of (lo_Hz, hi_Hz) known persistent-RFI ranges.
#                   Leave as [] if you don't have a vetted list, the
#                   report will say so explicitly rather than silently
#                   claiming "no RFI found".
# std_flux_cals   : either a flat list of exact field-name strings, or
#                   a dict of {canonical_name: [known aliases]} if the
#                   same calibrator shows up under multiple naming
#                   conventions in different datasets.
# default_quack_s : a typical quack/edge-flagging interval in seconds,
#                   used only for the "is this scan too short" check.
#                   Set to None if there's no standard convention for
#                   this telescope/pipeline, that check is then skipped.
# ---------------------------------------------------------------------
TELESCOPE_PROFILES = {
    'gmrt': {
        'label': 'GMRT / uGMRT',
        'match_names': ['GMRT'],
        # name to look up in CASA's own Observatories table (me.observatory())
        # for elevation/day-night calculations, so no coordinates are hardcoded here
        'measures_site_name': 'GMRT',
        # Central square C00-C14, East arm E02-E06, South arm S01-S06,
        # West arm W01-W06 (~30 antennas total).
        'antenna_pattern': re.compile(r'^(C0[0-9]|C1[0-4]|E0[2-6]|S0[1-6]|W0[1-6])$', re.IGNORECASE),
        'band_ranges': [
            # uGMRT wideband receivers (GMRT_specs.pdf, NCRA)
            ('uGMRT Band 2', 120e6, 250e6),
            ('uGMRT Band 3', 250e6, 500e6),
            ('uGMRT Band 4', 550e6, 850e6),
            ('uGMRT Band 5', 1000e6, 1460e6),
            # legacy GMRT narrowband receivers (GMRT observer's manual)
            ('Legacy 150 MHz', 130e6, 170e6),
            ('Legacy 235 MHz', 225e6, 245e6),
            ('Legacy 325 MHz', 300e6, 360e6),
            ('Legacy 610 MHz', 580e6, 660e6),
            ('Legacy 1420 MHz', 1000e6, 1450e6),
        ],
        # persistent RFI ranges hardcoded in CAPTURE's flagbadfreq step
        'rfi_bands': [
            (0.36e9, 0.3796e9),
            (0.486e9, 0.49355e9),
            (0.7646e9, 0.769092e9),
            (0.8808e9, 0.885596e9),
        ],
        'std_flux_cals': ['3C48', '3C147', '3C286', '0542+498', '1331+305', '0137+331'],
        'default_quack_s': 10.0,  # matches CAPTURE's config_capture.ini default
        # Known physical dish diameter, used only in FITS mode (below) as a
        # fallback for the field-of-view calculation, since a raw UVFITS
        # AN table doesn't carry a dish-diameter column the way an MS's
        # ANTENNA table does. MS mode always reads the real value from the
        # MS itself instead of this constant. GMRT/uGMRT dishes are 45m.
        'dish_diameter_m': 45.0,
    },
    'vla': {
        'label': 'VLA / EVLA',
        'match_names': ['VLA', 'EVLA'],
        'measures_site_name': 'VLA',
        # standard EVLA antenna naming is ea01-ea28; older pre-upgrade
        # data may use other conventions, treat mismatches as informational.
        'antenna_pattern': re.compile(r'^ea\d{2}$', re.IGNORECASE),
        'band_ranges': [
            # NRAO VLA Observational Status Summary / frequency bands page
            ('4-band', 54e6, 86e6),
            ('P-band', 200e6, 500e6),
            ('L-band', 1.0e9, 2.0e9),
            ('S-band', 2.0e9, 4.0e9),
            ('C-band', 4.0e9, 8.0e9),
            ('X-band', 8.0e9, 12.0e9),
            ('Ku-band', 12.0e9, 18.0e9),
            ('K-band', 18.0e9, 26.5e9),
            ('Ka-band', 26.5e9, 40.0e9),
            ('Q-band', 40.0e9, 50.0e9),
        ],
        # No vetted persistent-RFI frequency list for VLA is included
        # here (unlike GMRT's, which comes straight from CAPTURE's own
        # code). RFI at VLA is also much more site/band/time dependent
        # (GPS/Iridium in L-band, etc.). Add ranges here as you build
        # up your own list, left empty rather than guessed.
        'rfi_bands': [],
        # Perley & Butler (2013) standard flux density calibrators.
        # Included as both the common "3C" name and the J2000
        # coordinate-based name VLA field tables often use instead.
        'std_flux_cals': {
            '3C48':  ['3C48', 'J0137+3309', '0137+331'],
            '3C138': ['3C138', 'J0521+1638', '0521+166'],
            '3C147': ['3C147', 'J0542+4951', '0542+498'],
            '3C286': ['3C286', 'J1331+3030', '1331+305'],
        },
        # VLA/NRAO CASA pipeline quack conventions vary by project/band;
        # there isn't a single standard value to check against here.
        'default_quack_s': None,
        # Known physical dish diameter, used only in FITS mode as a
        # fallback (see the GMRT profile's comment above). VLA dishes are 25m.
        'dish_diameter_m': 25.0,
    },
    # To add another telescope (e.g. MeerKAT, ATCA), copy one of the
    # blocks above and fill in what you actually know. Leave
    # rfi_bands=[] and default_quack_s=None if unverified rather than
    # guessing.
}


def flatten_cals(entry):
    """std_flux_cals can be a flat list of names, or a dict of
    canonical_name -> [aliases]. Either way, return the flat set of
    strings that should be matched as a flux/bandpass calibrator."""
    if isinstance(entry, dict):
        out = set()
        for aliases in entry.values():
            out.update(aliases)
        return out
    return set(entry)


# Standard Stokes/correlation type codes, CASA MS convention (this is what
# the MS's POLARIZATION/CORR_TYPE column uses).
STOKES_CODES = {
    1: 'I', 2: 'Q', 3: 'U', 4: 'V',
    5: 'RR', 6: 'RL', 7: 'LR', 8: 'LL',
    9: 'XX', 10: 'XY', 11: 'YX', 12: 'YY',
}

# Stokes/correlation codes as they appear in a raw FITS header's STOKES
# axis (CRVAL/CDELT/CRPIX), AIPS Memo 117 convention. This is a DIFFERENT
# numbering from STOKES_CODES above (note the negative values), used only
# by the FITS-mode functions below, never mix the two up.
FITS_STOKES_CODES = {
    1: 'I', 2: 'Q', 3: 'U', 4: 'V',
    -1: 'RR', -2: 'LL', -3: 'RL', -4: 'LR',
    -5: 'XX', -6: 'YY', -7: 'XY', -8: 'YX',
}


def corr_type_strings(msname):
    """Actual correlation product names (e.g. ['RR','LL']), read directly
    from the POLARIZATION subtable rather than just a NUM_CORR count."""
    tb_t = table()
    tb_t.open(msname + '/POLARIZATION')
    corr_type = tb_t.getcol('CORR_TYPE')
    tb_t.close()
    codes = corr_type[:, 0] if corr_type.ndim > 1 else corr_type
    return [STOKES_CODES.get(int(c), 'code%d' % c) for c in codes]


def _max_baseline_from_positions(positions):
    """Longest pairwise separation among a set of (n, 3) antenna
    positions, in whatever units the positions are in. O(n^2) over
    antenna pairs, trivial for any real array size. Shared by both the
    MS-mode and FITS-mode max-baseline functions below."""
    n = positions.shape[0]
    maxb = 0.0
    for i in range(n):
        d = np.linalg.norm(positions[i + 1:] - positions[i], axis=1)
        if d.size:
            maxb = max(maxb, float(d.max()))
    return maxb


def max_baseline_m(msname):
    """Longest antenna-antenna separation in the array, in meters, from
    the ANTENNA subtable's ECEF positions."""
    tb_t = table()
    tb_t.open(msname + '/ANTENNA')
    positions = tb_t.getcol('POSITION').T  # (nant, 3)
    tb_t.close()
    return _max_baseline_from_positions(positions)


def dish_diameter_m(msname):
    """Representative single-dish diameter, read from the ANTENNA subtable
    (DISH_DIAMETER column) rather than hardcoded per telescope, so this
    works for GMRT (45m) and VLA (25m) from the same code path."""
    tb_t = table()
    tb_t.open(msname + '/ANTENNA')
    diam = tb_t.getcol('DISH_DIAMETER')
    tb_t.close()
    diam = diam[diam > 0]
    return float(np.median(diam)) if diam.size else None


def next_pow2(n):
    return int(2 ** np.ceil(np.log2(max(n, 1))))


# ---------------------------------------------------------------------
# FITS-mode helpers
#
# These read a raw UVFITS "random groups" file directly (the AIPS
# convention GMRT/uGMRT correlator output uses), without ever building a
# CASA MS. They cover only what's genuinely available before MS import:
# antenna positions/names (AN table), frequency setup (header + FQ
# table), polarization (header STOKES axis), source list/positions (SU
# table, or the header's OBJECT/OBSRA/OBSDEC for single-source files),
# and the observation time span (DATE/INTTIM random-group parameters,
# read without touching the visibility data array itself, so this stays
# fast regardless of file size).
#
# What's deliberately NOT attempted here: scan boundaries (CASA's own
# importer infers these from time gaps per source, a nontrivial piece of
# logic not worth silently re-implementing and risking a mismatch),
# per-antenna/per-channel flagging status (needs the same import step,
# and a freshly-correlated FITS file has essentially nothing flagged yet
# anyway), and anything derived from either of those.
# ---------------------------------------------------------------------

def _require_astropy():
    try:
        from astropy.io import fits as pyfits
        return pyfits
    except ImportError:
        raise ImportError(
            "FITS mode needs the 'astropy' package to read UVFITS binary tables "
            "(CASA's own msmetadata/table tools only work on Measurement Sets, "
            "not raw FITS). Install it with: pip install astropy")


def _fits_axis_for_ctype(header, ctype_prefix):
    """Find which NAXISn/CRVALn/etc. axis number has a given CTYPE
    (e.g. 'FREQ', 'STOKES'), rather than assuming a fixed axis order,
    since that order isn't guaranteed to be the same across files."""
    naxis = header.get('NAXIS', 0)
    for i in range(1, naxis + 1):
        ctype = str(header.get('CTYPE%d' % i, '')).strip().upper()
        if ctype.startswith(ctype_prefix.upper()):
            return i
    return None


def fits_telescope_name(fitsfile):
    """TELESCOP header keyword, the FITS-mode equivalent of the MS's
    OBSERVATION table telescope name."""
    pyfits = _require_astropy()
    with pyfits.open(fitsfile) as hdul:
        return str(hdul[0].header.get('TELESCOP', '')).strip()


def fits_antenna_info(fitsfile):
    """Antenna names and positions from the FITS AN ('AIPS AN') table.
    Positions are only ever used as PAIRWISE DIFFERENCES below (for max
    baseline), so the ARRAYX/ARRAYY/ARRAYZ array-center offset in the AN
    table's own header, which would convert these relative STABXYZ
    values into true absolute geocentric coordinates, cancels out
    exactly and doesn't need to be added back in."""
    pyfits = _require_astropy()
    with pyfits.open(fitsfile) as hdul:
        an = hdul['AIPS AN']
        names = [str(n).strip() for n in an.data['ANNAME']]
        positions = np.array(an.data['STABXYZ'], dtype=float)
    return names, positions


def fits_frequency_info(fitsfile):
    """Per-channel frequencies (Hz) for each IF/subband, built from the
    primary header's FREQ axis (CRVAL/CDELT/CRPIX/NAXIS) plus the FQ
    ('AIPS FQ') table's per-IF frequency offset, if present. Returns a
    list of (if_index, freq_array_hz). Falls back to a single IF at the
    header's own reference frequency if there's no FQ table or it's an
    unexpected shape, rather than guessing a multi-IF layout."""
    pyfits = _require_astropy()
    with pyfits.open(fitsfile) as hdul:
        hdr = hdul[0].header
        freq_axis = _fits_axis_for_ctype(hdr, 'FREQ')
        if freq_axis is None:
            raise ValueError("No FREQ axis (CTYPEn) found in the primary header.")
        crval = float(hdr['CRVAL%d' % freq_axis])
        cdelt = float(hdr['CDELT%d' % freq_axis])
        crpix = float(hdr['CRPIX%d' % freq_axis])
        nchan = int(hdr['NAXIS%d' % freq_axis])
        chan_idx = np.arange(1, nchan + 1)
        base_freqs = crval + (chan_idx - crpix) * cdelt  # Hz, one IF's channels

        if_offsets = np.array([0.0])
        try:
            fq = hdul['AIPS FQ']
            iffreq = np.array(fq.data['IF FREQ'], dtype=float)
            if_offsets = np.atleast_1d(iffreq[0]) if iffreq.ndim > 1 else np.atleast_1d(iffreq)
        except Exception:
            pass  # no FQ table, or an unexpected layout: treat as single-IF at CRVAL

        return [(i, base_freqs + off) for i, off in enumerate(if_offsets)]


def fits_polarization_products(fitsfile):
    """Polarization products (e.g. ['RR','LL']) from the primary
    header's STOKES axis, mapped through FITS_STOKES_CODES (the FITS
    header convention, not the MS one, see the comment on that dict)."""
    pyfits = _require_astropy()
    with pyfits.open(fitsfile) as hdul:
        hdr = hdul[0].header
        stokes_axis = _fits_axis_for_ctype(hdr, 'STOKES')
        if stokes_axis is None:
            raise ValueError("No STOKES axis (CTYPEn) found in the primary header.")
        crval = float(hdr['CRVAL%d' % stokes_axis])
        cdelt = float(hdr['CDELT%d' % stokes_axis])
        crpix = float(hdr['CRPIX%d' % stokes_axis])
        npol = int(hdr['NAXIS%d' % stokes_axis])
        codes = [int(round(crval + (i + 1 - crpix) * cdelt)) for i in range(npol)]
    return [FITS_STOKES_CODES.get(c, 'code%d' % c) for c in codes]


def fits_source_list(fitsfile):
    """List of (name, ra_deg, dec_deg). Multi-source files carry this in
    the SU ('AIPS SU') table; single-source files (no SU table) carry it
    in the primary header's OBJECT/OBSRA/OBSDEC keywords instead."""
    pyfits = _require_astropy()
    with pyfits.open(fitsfile) as hdul:
        try:
            su = hdul['AIPS SU']
            names = [str(n).strip() for n in su.data['SOURCE']]
            ra = np.array(su.data['RAEPO'], dtype=float)
            dec = np.array(su.data['DECEPO'], dtype=float)
            return list(zip(names, ra, dec))
        except KeyError:
            hdr = hdul[0].header
            name = str(hdr.get('OBJECT', 'UNKNOWN')).strip()
            ra = hdr.get('OBSRA')
            dec = hdr.get('OBSDEC')
            if ra is None or dec is None:
                raise ValueError("No AIPS SU table and no OBSRA/OBSDEC in the primary header.")
            return [(name, float(ra), float(dec))]


def fits_time_range_mjd_sec(fitsfile):
    """Observation start/end as MJD-seconds (matching what the
    elevation/sun functions below expect), plus the INTTIM random
    parameter if present. Reads only the random-group PARAMETER columns
    (DATE, INTTIM), not the visibility data array itself, so this stays
    fast regardless of file size. The split-precision DATE/DATE
    convention (two random parameters both named DATE, for JD integer +
    fractional parts) is handled automatically: astropy's GroupData.par()
    sums same-named parameters together, verified against a test file."""
    pyfits = _require_astropy()
    with pyfits.open(fitsfile) as hdul:
        data = hdul[0].data
        jd = np.array(data.par('DATE'), dtype=float)
        inttim = None
        try:
            inttim = np.array(data.par('INTTIM'), dtype=float)
        except Exception:
            pass
    mjd_sec = (jd - 2400000.5) * 86400.0
    return float(mjd_sec.min()), float(mjd_sec.max()), inttim


def get_integration_time(msname):
    """Median sampling interval in seconds. Reads only a small chunk of
    the MAIN table (not the whole column) so this stays fast on large MS."""
    tb_t = table()
    tb_t.open(msname)
    n = min(2000, tb_t.nrows())
    intervals = tb_t.getcol('INTERVAL', 0, n)
    tb_t.close()
    return float(np.median(intervals))


def deg_to_hms(ra_deg):
    h = (ra_deg % 360.0) / 15.0
    hh = int(h)
    m = (h - hh) * 60
    mm = int(m)
    s = (m - mm) * 60
    return "%02dh%02dm%05.2fs" % (hh, mm, s)


def deg_to_dms(dec_deg):
    sign = '-' if dec_deg < 0 else '+'
    d = abs(dec_deg)
    dd = int(d)
    m = (d - dd) * 60
    mm = int(m)
    s = (m - mm) * 60
    return "%s%02dd%02dm%05.2fs" % (sign, dd, mm, s)


def angsep_deg(ra1, dec1, ra2, dec2):
    """Great-circle separation between two sky positions, in degrees."""
    ra1r, dec1r, ra2r, dec2r = (np.radians(x) for x in (ra1, dec1, ra2, dec2))
    cosang = (np.sin(dec1r) * np.sin(dec2r)
              + np.cos(dec1r) * np.cos(dec2r) * np.cos(ra1r - ra2r))
    return np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))


def get_observatory_position(me_t, candidate_names):
    """Try each candidate name against CASA's own Observatories table
    (me.observatory()) and return the first that resolves, or None.
    No coordinates are hardcoded here, CASA already ships this database."""
    for name in candidate_names:
        if not name:
            continue
        try:
            return me_t.observatory(name)
        except Exception:
            continue
    return None


def elevation_deg(me_t, qa_t, obspos, ra_deg, dec_deg, mjd_sec):
    me_t.doframe(obspos)
    me_t.doframe(me_t.epoch('utc', qa_t.quantity(mjd_sec, 's')))
    d = me_t.direction('J2000', qa_t.quantity(ra_deg, 'deg'), qa_t.quantity(dec_deg, 'deg'))
    azel = me_t.measure(d, 'AZEL')
    return float(np.degrees(azel['m1']['value']))


def sun_elevation_deg(me_t, qa_t, obspos, mjd_sec):
    me_t.doframe(obspos)
    me_t.doframe(me_t.epoch('utc', qa_t.quantity(mjd_sec, 's')))
    azel = me_t.measure(me_t.direction('SUN'), 'AZEL')
    return float(np.degrees(azel['m1']['value']))


def sun_condition(elev_deg):
    if elev_deg > 0:
        return "daytime"
    elif elev_deg > -12:
        return "twilight (dawn/dusk terminator, historically the roughest ionospheric period)"
    else:
        return "night"


def parse_args():
    input_path = None
    telescope = None
    perchannel = None  # None = use CHECK_PERCHANNEL_FLAGS default, True/False = override
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        a = args[i]
        if a == '--telescope' and i + 1 < len(args):
            telescope = args[i + 1].lower()
            i += 2
            continue
        if a.startswith('--telescope='):
            telescope = a.split('=', 1)[1].lower()
            i += 1
            continue
        if a == '--nochanflags':
            perchannel = False
            i += 1
            continue
        if a == '--chanflags':
            perchannel = True
            i += 1
            continue
        if a.endswith('.ms') or a.endswith('.MS') or a.endswith('.fits') or a.endswith('.FITS'):
            input_path = a
        i += 1
    return input_path, telescope, perchannel


def find_input(input_arg):
    """Resolve the input path: an explicit argument, the MSNAME variable
    at the top of the script, or (if exactly one exists) whatever single
    .ms/.MS/.fits/.FITS is sitting in the current directory."""
    if input_arg:
        return input_arg
    if MSNAME:
        return MSNAME
    candidates = [d for d in os.listdir('.')
                  if d.endswith('.ms') or d.endswith('.MS')
                  or d.endswith('.fits') or d.endswith('.FITS')]
    if len(candidates) == 1:
        return candidates[0]
    elif len(candidates) > 1:
        raise ValueError("Multiple candidate inputs found (%s). Pass one explicitly."
                          % ', '.join(candidates))
    else:
        raise ValueError("No .ms/.MS or .fits/.FITS file found. Set MSNAME at the top "
                          "of the script or pass a path as an argument.")


def input_kind(path):
    if path.endswith('.ms') or path.endswith('.MS'):
        return 'ms'
    if path.endswith('.fits') or path.endswith('.FITS'):
        return 'fits'
    raise ValueError("Don't recognize the input type for '%s' (expected .ms/.MS or "
                      ".fits/.FITS)." % path)


def detect_telescope(obs_name):
    if not obs_name:
        return None
    upper = obs_name.upper()
    for key, prof in TELESCOPE_PROFILES.items():
        if any(m.upper() in upper for m in prof['match_names']):
            return key
    return None


def resolve_telescope(telescope_arg, obs_name, out):
    detected = detect_telescope(obs_name)
    chosen = telescope_arg or TELESCOPE or detected

    if chosen is None:
        raise ValueError(
            "Could not determine which telescope this is (OBSERVATION table says '%s'). "
            "Pass --telescope gmrt or --telescope vla explicitly. "
            "Known profiles: %s" % (obs_name, ', '.join(TELESCOPE_PROFILES)))

    if chosen not in TELESCOPE_PROFILES:
        raise ValueError("Unknown telescope '%s'. Known profiles: %s"
                          % (chosen, ', '.join(TELESCOPE_PROFILES)))

    if telescope_arg and detected and telescope_arg != detected:
        out("  ! You specified --telescope %s but the MS OBSERVATION table looks like %s (%s)."
            % (telescope_arg, detected, obs_name))
        out("    Proceeding with your choice (%s), but double-check this is really the right MS." % telescope_arg)
    elif telescope_arg:
        out("  Telescope: %s (user-specified)" % TELESCOPE_PROFILES[chosen]['label'])
    else:
        out("  Telescope: %s (auto-detected from OBSERVATION table: '%s')" % (TELESCOPE_PROFILES[chosen]['label'], obs_name))

    return chosen


def classify_fields(fields, phase_cals, flux_cals):
    ampcals, pcals, targets = [], [], []
    for f in fields:
        if f in flux_cals:
            ampcals.append(f)
        elif f in phase_cals:
            pcals.append(f)
        else:
            targets.append(f)
    return ampcals, pcals, targets


def load_phase_cals(path):
    if not os.path.isfile(path):
        return [], "WARNING: %s not found in this directory. Phase calibrators may be misclassified as targets." % path
    try:
        cals = np.loadtxt(path, dtype=str)
        return list(np.atleast_1d(cals)), None
    except Exception as e:
        return [], "WARNING: could not read %s (%s)" % (path, e)


def guess_band(center_hz, band_ranges):
    for name, lo, hi in band_ranges:
        if lo <= center_hz <= hi:
            return name
    return "unrecognized (%.1f MHz center, not in this profile's band_ranges)" % (center_hz / 1e6)


def safe(fn, default=None):
    try:
        return fn()
    except Exception as e:
        return "unavailable (%s)" % e


def main_ms(msname, telescope_arg, perchannel_arg):
    check_perchannel = CHECK_PERCHANNEL_FLAGS if perchannel_arg is None else perchannel_arg
    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out(" read_listobs.py SUMMARY: %s" % msname)
    out("=" * 72)

    msmd_t = msmetadata()
    msmd_t.open(msname)

    # ---------------- basic info + telescope profile ----------------
    out("\nBASIC INFO")
    n_scans = safe(lambda: len(msmd_t.scannumbers()))
    n_fields = safe(lambda: msmd_t.nfields())
    obsname = safe(lambda: msmd_t.observatorynames()[0], "")
    if not isinstance(obsname, str):
        obsname = ""
    out("  Observatory (from MS)  : %s" % (obsname or "(unknown)"))
    out("  Total scans            : %s" % n_scans)
    out("  Total fields           : %s" % n_fields)

    telescope_key = resolve_telescope(telescope_arg, obsname, out)
    profile = TELESCOPE_PROFILES[telescope_key]

    # ---------------- observation time range ----------------
    qa_t = quanta()
    me_t = measures()
    obs_t0 = obs_t1 = None
    try:
        tb_obs = table()
        tb_obs.open(msname + '/OBSERVATION')
        trange = tb_obs.getcol('TIME_RANGE')
        tb_obs.close()
        obs_t0, obs_t1 = float(trange.min()), float(trange.max())
        d0 = qa_t.time(qa_t.quantity(obs_t0, 's'), form='ymd')[0]
        d1 = qa_t.time(qa_t.quantity(obs_t1, 's'), form='ymd')[0]
        out("  Observation window     : %s to %s (%.1f min)" % (d0, d1, (obs_t1 - obs_t0) / 60.0))
    except Exception as e:
        out("  Observation window     : unavailable (%s)" % e)

    # observatory position, for elevation/day-night calculations below.
    # Looked up from CASA's own Observatories table (no coordinates
    # hardcoded here), tried under the MS's own name first, then the
    # profile's known CASA-registered name as a fallback.
    obspos = get_observatory_position(me_t, [obsname, profile.get('measures_site_name')])

    # ---------------- observing conditions (time of day) ----------------
    # Day/night/twilight at the telescope during the track matters a lot
    # for how badly the ionosphere (mainly) will have corrupted phases,
    # which is one of the drivers of how much selfcal will need to do.
    out("\nOBSERVING CONDITIONS")
    sun_conditions = []  # (label, elevation_deg, condition_str), reused in the risk-factor summary
    if obspos is not None and obs_t0 is not None:
        try:
            for label, t in (("start", obs_t0), ("mid", (obs_t0 + obs_t1) / 2.0), ("end", obs_t1)):
                sel = sun_elevation_deg(me_t, qa_t, obspos, t)
                cond = sun_condition(sel)
                sun_conditions.append((label, sel, cond))
                out("  Sun elevation at %-5s: %6.1f deg (%s)" % (label, sel, cond))
        except Exception as e:
            out("  Sun elevation unavailable (%s)" % e)
    else:
        out("  Could not resolve an observatory position for '%s' in CASA's Observatories"
            % (obsname or profile['label']))
        out("  table, so day/night/elevation checks are skipped.")

    # ---------------- antennas ----------------
    out("\nANTENNAS")
    ant_names = safe(lambda: list(msmd_t.antennanames()), [])
    if isinstance(ant_names, list):
        out("  Count: %d" % len(ant_names))
        out("  Names: %s" % ', '.join(ant_names))
        odd = [a for a in ant_names if not profile['antenna_pattern'].match(a)]
        if odd:
            out("  ! Name(s) that don't match the expected %s naming pattern (informational, verify manually): %s"
                % (profile['label'], ', '.join(odd)))
    else:
        out("  %s" % ant_names)

    # ---------------- frequency setup ----------------
    out("\nFREQUENCY SETUP")
    nspw = safe(lambda: msmd_t.nspw(), 0)
    out("  Number of spectral windows: %s" % nspw)
    all_freqs = {}  # spw index -> frequency array, keyed (not positional) so a
                    # failed spw never misaligns lookups done later by spw index
    cell_lo = cell_hi = None  # set below if a suggested cell size could be computed
    if isinstance(nspw, int) and nspw > 0:
        for spw in range(nspw):
            freqs = safe(lambda spw=spw: msmd_t.chanfreqs(spw))
            nchan = safe(lambda spw=spw: msmd_t.nchan(spw))
            bw = safe(lambda spw=spw: msmd_t.bandwidths(spw))
            if isinstance(freqs, np.ndarray):
                all_freqs[spw] = freqs
                out("  spw %d: %d channels, bandwidth %.3f MHz, range %.3f-%.3f MHz"
                    % (spw, nchan, bw / 1e6, freqs.min() / 1e6, freqs.max() / 1e6))
            else:
                out("  spw %d: %s" % (spw, freqs))
    center_freq_hz = None
    band_name = None
    if all_freqs:
        allf = np.concatenate(list(all_freqs.values()))
        center_freq_hz = (allf.min() + allf.max()) / 2.0
        band_name = guess_band(center_freq_hz, profile['band_ranges'])
        out("  Overall range: %.3f-%.3f MHz, center %.3f MHz" % (allf.min() / 1e6, allf.max() / 1e6, center_freq_hz / 1e6))
        out("  Likely band (%s): %s" % (profile['label'], band_name))

        # Suggested imaging cell size: synthesized beam ~ wavelength / max
        # baseline (the standard back-of-envelope interferometer resolution
        # estimate), then Nyquist-sample that beam. 3-5 pixels/beam is a
        # range, not a single rule: toward 5 for precise point-source
        # photometry or an elongated/non-Gaussian beam (poor uv-coverage,
        # e.g. a short low-elevation track); toward 3 when field of view
        # and compute cost are the binding constraint (imsize scales as
        # FOV/cell, and FOV gets large at low frequency). beam/4 is given
        # as a single practical default if you just want one number.
        cell_lo = cell_hi = cell_default = None
        try:
            bmax = max_baseline_m(msname)
            wavelength_m = 299792458.0 / center_freq_hz
            beam_arcsec = 206265.0 * wavelength_m / bmax
            cell_lo, cell_hi = beam_arcsec / 5.0, beam_arcsec / 3.0
            cell_default = beam_arcsec / 4.0
            out("  Max baseline: %.0f m -> synthesized beam ~%.2f arcsec (natural weighting, "
                "no tapering)" % (bmax, beam_arcsec))
            out("  Suggested cell size: ~%.2f-%.2f arcsec (3-5 px/beam), %.2f arcsec if you just want one number (4 px/beam)"
                % (cell_lo, cell_hi, cell_default))
        except Exception as e:
            out("  Suggested cell size unavailable (%s)" % e)

        # Field of view: single-dish primary beam ~ wavelength / dish
        # diameter (simple lambda/D convention, same one used in
        # vla_imaging_parameter_calc.ipynb; some conventions add a 1.22
        # factor for the first Airy null, that would widen this ~20%).
        # Dish diameter is read from the MS itself (ANTENNA/DISH_DIAMETER),
        # not hardcoded, so this is the same for GMRT (45m) and VLA (25m).
        if cell_default is not None:
            try:
                ddiam = dish_diameter_m(msname)
                if ddiam:
                    fov_arcsec = 206265.0 * wavelength_m / ddiam
                    imsize_needed = fov_arcsec / cell_default
                    imsize_pow2 = next_pow2(imsize_needed)
                    out("  Dish diameter: %.1f m -> field of view ~%.1f arcsec (%.2f arcmin)"
                        % (ddiam, fov_arcsec, fov_arcsec / 60.0))
                    out("  Suggested image size: %d px (next power of 2 above the %.0f px the FOV/cell needs)"
                        % (imsize_pow2, imsize_needed))
                else:
                    out("  Field of view / image size: no valid DISH_DIAMETER found in ANTENNA table.")
            except Exception as e:
                out("  Field of view / image size unavailable (%s)" % e)

        if profile['rfi_bands']:
            hits = [(lo, hi) for lo, hi in profile['rfi_bands'] if np.any((allf >= lo) & (allf <= hi))]
            if hits:
                out("  ! Known RFI-prone frequency range(s) fall inside your band:")
                for lo, hi in hits:
                    out("      %.1f-%.1f MHz" % (lo / 1e6, hi / 1e6))
            else:
                out("  No overlap with this profile's known persistent-RFI list.")
        else:
            out("  (No curated RFI-frequency list for %s in this script yet, not checked.)" % profile['label'])

    # polarization products
    pol_products = []
    try:
        pol_products = corr_type_strings(msname)
        out("  Polarization products: %s" % ', '.join(pol_products))
    except Exception as e:
        out("  Polarization products: unavailable (%s)" % e)

    # integration time
    try:
        dt = get_integration_time(msname)
        out("  Integration time (median sample interval): %.2f s" % dt)
    except Exception as e:
        out("  Integration time: unavailable (%s)" % e)

    # ---------------- field classification ----------------
    out("\nFIELD CLASSIFICATION")
    phase_cals, phase_cal_warn = load_phase_cals(PHASE_CAL_FILE)
    if phase_cal_warn:
        out("  " + phase_cal_warn)
    flux_cals = flatten_cals(profile['std_flux_cals'])
    fields = safe(lambda: list(msmd_t.fieldnames()), [])
    ampcals, pcals, targets = classify_fields(fields, phase_cals, flux_cals)
    out("  Flux/bandpass calibrators : %s" % (', '.join(ampcals) if ampcals else '(none found)'))
    out("  Phase calibrators         : %s" % (', '.join(pcals) if pcals else '(none found)'))
    out("  Targets                   : %s" % (', '.join(targets) if targets else '(none found)'))

    callike_targets = [t for t in targets if CALLIKE_NAME.match(t)]
    if callike_targets:
        out("  ! These 'targets' have coordinate-style names typical of calibrators,")
        out("    but aren't in %s or this profile's standard list. Verify they aren't" % PHASE_CAL_FILE)
        out("    meant to be a phase calibrator (a missing catalog entry will make a")
        out("    pipeline treat them as the science target):")
        out("      %s" % ', '.join(callike_targets))

    if ampcals and not pcals:
        out("  ! No separate phase calibrator found. Pipelines that fall back to using")
        out("    the amplitude/bandpass calibrator as the phase calibrator (e.g. CAPTURE)")
        out("    will pick the one with the most scans. Confirm that's what you want.")

    # ---------------- field positions ----------------
    out("\nFIELD POSITIONS")
    field_radec = {}
    for fid, fname in enumerate(fields):
        try:
            d = msmd_t.phasecenter(fid)
            ra_deg = np.degrees(d['m0']['value']) % 360.0
            dec_deg = np.degrees(d['m1']['value'])
            field_radec[fname] = (ra_deg, dec_deg)
            out("  %-20s RA %s   Dec %s" % (fname, deg_to_hms(ra_deg), deg_to_dms(dec_deg)))
        except Exception as e:
            out("  %-20s position unavailable (%s)" % (fname, e))

    max_cal_sep_deg = None
    if ampcals or pcals:
        out("\n  Calibrator-target separations:")
        for cal in ampcals + pcals:
            if cal not in field_radec:
                continue
            for tgt in targets:
                if tgt not in field_radec:
                    continue
                sep = angsep_deg(field_radec[cal][0], field_radec[cal][1],
                                  field_radec[tgt][0], field_radec[tgt][1])
                out("    %s <-> %s : %.1f deg" % (cal, tgt, sep))
                if max_cal_sep_deg is None or sep > max_cal_sep_deg:
                    max_cal_sep_deg = sep
        out("    (No universal cutoff here, but the larger this separation is, the less")
        out("     that calibrator's gain/phase solutions represent the sky above your")
        out("     target, especially if it's serving as your only calibrator.)")

    # ---------------- scan timing per field ----------------
    out("\nSCAN TIMING PER FIELD")
    out("  %-20s %8s %14s %12s" % ("Field", "Scans", "OnSourceTime(s)", "AvgScan(s)"))
    quack_s = profile['default_quack_s']
    field_tspan = {}    # field -> (min_time, max_time) across all its scans, for elevation range below
    field_maxgap = {}   # field -> longest gap between consecutive scans on that field
    field_avg_scan = {} # field -> average scan duration, reused in the risk-factor summary
    for f in fields:
        try:
            scans = msmd_t.scansforfield(f)
            durations = []
            scan_bounds = []
            for sc in scans:
                times = msmd_t.timesforscan(sc)
                if len(times) > 1:
                    durations.append(times.max() - times.min())
                    scan_bounds.append((times.min(), times.max()))
            total = sum(durations)
            avg = total / len(durations) if durations else 0.0
            field_avg_scan[f] = avg
            out("  %-20s %8d %14.1f %12.1f" % (f, len(scans), total, avg))
            if scan_bounds:
                field_tspan[f] = (min(b[0] for b in scan_bounds), max(b[1] for b in scan_bounds))
                scan_bounds.sort()
                gaps = [scan_bounds[i + 1][0] - scan_bounds[i][1] for i in range(len(scan_bounds) - 1)]
                if gaps:
                    field_maxgap[f] = max(gaps)
            if quack_s is not None:
                short = [d for d in durations if d < 4 * quack_s]
                if short:
                    out("    ! %d scan(s) on this field are short (<%.0fs) relative to a %.0fs quack"
                        % (len(short), 4 * quack_s, quack_s))
                    out("      interval on each end. Edge-flagging may remove most of the useful data.")
        except Exception as e:
            out("  %-20s error reading scan timing (%s)" % (f, e))
    if quack_s is None:
        out("  (No standard quack interval defined for %s in this script; short-scan check skipped.)" % profile['label'])

    # ---------------- elevation during each field's track ----------------
    # Low elevation means more troposphere/ionosphere path length and, for
    # GMRT's Central Square, more shadowing between antennas.
    target_min_elev = None
    if obspos is not None and field_tspan:
        out("\nELEVATION DURING TRACK")
        for f, (tmin, tmax) in field_tspan.items():
            if f not in field_radec:
                continue
            try:
                ra_deg, dec_deg = field_radec[f]
                e0 = elevation_deg(me_t, qa_t, obspos, ra_deg, dec_deg, tmin)
                e1 = elevation_deg(me_t, qa_t, obspos, ra_deg, dec_deg, tmax)
                lo, hi = min(e0, e1), max(e0, e1)
                # airmass ~ 1/sin(elevation): a plane-parallel-atmosphere approximation,
                # good enough here to turn "elevation dropped" into an actual multiple
                # of the sky path length you get looking straight up (zenith, airmass 1).
                airmass_lo = 1.0 / np.sin(np.radians(lo)) if lo > 0 else float('inf')
                out("  %-20s %5.1f to %5.1f deg  (airmass at lowest point: %.2fx zenith)"
                    % (f, lo, hi, airmass_lo))
                if f in targets and (target_min_elev is None or lo < target_min_elev):
                    target_min_elev = lo
            except Exception as e:
                out("  %-20s elevation unavailable (%s)" % (f, e))

    msmd_t.close()

    # ---------------- current flag status ----------------
    # This is the one section that reads actual visibility data (the FLAG
    # column) instead of lightweight metadata tables, so it's the slowest
    # part of the script by far, and scales with scans x baselines x
    # channels x polarizations. spwchan=True (the per-channel/RFI view)
    # is what makes it expensive: skip it with --nochanflags, or set
    # CHECK_PERCHANNEL_FLAGS = False above, for a much faster run when
    # you don't need that breakdown.
    out("\nCURRENT FLAGGING STATUS (before any pipeline flagging)")
    if not check_perchannel:
        out("  (per-channel breakdown skipped: --nochanflags / CHECK_PERCHANNEL_FLAGS=False)")
    summ = None
    frac = None
    try:
        summ = flagdata(vis=msname, mode='summary', spwchan=check_perchannel)
        frac = summ.get('flagged', 0) / summ.get('total', 1) * 100.0
        out("  Total data flagged: %.2f%%" % frac)
    except Exception as e:
        out("  Could not compute flag summary (%s)" % e)

    # per-antenna: an antenna already carrying a much higher flagged
    # fraction than the rest, before you've even started calibrating,
    # is the standard first clue that it was dead/bad during the track.
    if summ:
        try:
            ant_summ = summ.get('antenna', {})
            ant_fracs = [(a, c['flagged'] / c['total'] * 100.0)
                         for a, c in ant_summ.items() if c.get('total', 0) > 0]
            if ant_fracs:
                med = float(np.median([f for _, f in ant_fracs]))
                out("\n  Per-antenna flagged %% (median %.2f%%):" % med)
                for a, f in sorted(ant_fracs, key=lambda x: -x[1]):
                    tag = "  ! elevated, check this antenna" if f > max(2 * med, med + 10.0) else ""
                    out("    %-6s %6.2f%%%s" % (a, f, tag))
            else:
                out("  Per-antenna flag breakdown: no antennas with data.")
        except Exception as e:
            out("  Per-antenna flag breakdown unavailable (%s)" % e)

        # per-channel: a channel (or narrow group of channels) with a
        # flagged fraction that spikes well above its neighbors is the
        # standard signature of narrowband RFI. summ won't have a
        # 'spw:channel' key at all when check_perchannel was False
        # (spwchan=False), so .get() below degrades gracefully to the
        # "not available" branch without needing a separate guard here.
        try:
            chan_summ = summ.get('spw:channel', {})
            per_spw = {}
            for key, c in chan_summ.items():
                spw_str, chan_str = key.split(':')
                if c.get('total', 0) > 0:
                    per_spw.setdefault(int(spw_str), {})[int(chan_str)] = c['flagged'] / c['total'] * 100.0
            if per_spw:
                out("\n  Per-channel flagged %% (every channel above its spw's own median is listed, so")
                out("  nothing unusual is hidden behind a cutoff; '!' marks the ones that also clear the")
                out("  RFI_SPIKE_SIGMA/RFI_SPIKE_FLOOR_PCT heuristic cutoff defined near the top of this")
                out("  script, that heuristic is not a validated RFI detector, just a way to flag the")
                out("  channels that stand out the most within each spw):")
                for spw_i in sorted(per_spw):
                    chan_map = per_spw[spw_i]
                    chans_sorted = sorted(chan_map)
                    fracs = np.array([chan_map[c] for c in chans_sorted])
                    med = float(np.median(fracs))
                    spread = float(np.std(fracs))
                    thresh = med + max(RFI_SPIKE_SIGMA * spread, RFI_SPIKE_FLOOR_PCT)
                    above_med = sorted((c for c in chans_sorted if chan_map[c] > med),
                                        key=lambda c: -chan_map[c])
                    n_over_thresh = sum(1 for c in above_med if chan_map[c] > thresh)
                    out("    spw %d: median %.2f%%, heuristic cutoff %.1f%% (%d channel(s) clear it, "
                        "%d channel(s) above median overall)" % (spw_i, med, thresh, n_over_thresh, len(above_med)))
                    freqs = all_freqs.get(spw_i)
                    if not above_med:
                        out("      All channels at or below the median, nothing elevated to show.")
                    for c in above_med[:20]:
                        freq_str = "%.3f MHz" % (freqs[c] / 1e6) if freqs is not None and c < len(freqs) else "freq n/a"
                        marker = "!" if chan_map[c] > thresh else " "
                        out("      %s chan %4d (%s): %.1f%% flagged" % (marker, c, freq_str, chan_map[c]))
                    if len(above_med) > 20:
                        out("      ... and %d more channel(s) above the median (showing the top 20 by flagged %%)"
                            % (len(above_med) - 20))
            else:
                out("  Per-channel flag breakdown: not available from this flagdata summary.")
        except Exception as e:
            out("  Per-channel flag breakdown unavailable (%s)" % e)

    # ---------------- self-calibration risk factors ----------------
    # None of this can tell you whether selfcal will actually work, that
    # depends on the target being bright enough to solve stable gains
    # against, which only shows up once you look at real visibility SNR
    # or a first image. What metadata *can* do is flag the standard risk
    # factors that make selfcal more likely to be needed in the first
    # place. Treat this as a checklist of "reasons to expect trouble",
    # not a verdict.
    out("\nSELF-CALIBRATION RISK FACTORS (heuristic, not a verdict)")
    risk_factors = []

    if band_name:
        low_freq = ('Legacy' in band_name) or band_name in ('uGMRT Band 2', 'uGMRT Band 3', 'uGMRT Band 4',
                                                              '4-band', 'P-band')
        if low_freq:
            risk_factors.append("Low-frequency band (%s): ionospheric phase noise scales as 1/freq^2, "
                                 "usually the single biggest driver of selfcal need at these frequencies "
                                 "on %s." % (band_name, profile['label']))

    for label, sel, cond in sun_conditions:
        if 'twilight' in cond:
            risk_factors.append("Sun near horizon at track %s (%.1f deg): terminator crossings bring "
                                 "rapid ionospheric TEC gradients, often the roughest part of a day." % (label, sel))
        elif cond == 'daytime':
            risk_factors.append("Track %s is in local daytime (sun %.1f deg): daytime ionosphere is "
                                 "more disturbed than night." % (label, sel))

    if ampcals and not pcals:
        risk_factors.append("No separate phase calibrator: nothing is cycling in to track time-variable "
                             "phase during the target scans, selfcal is doing that job instead.")

    # Calibrator-target separation limit, from the NRAO/VLBA calibration
    # manual's phase-referencing guidance: ~4 deg below 5 GHz, ~5.7 deg
    # (0.1 rad) above 5 GHz. GMRT/uGMRT is always below 5 GHz.
    # (science.nrao.edu/facilities/vlba/docs/manuals/obsvlba/calibration)
    sep_limit_deg = 5.7 if (center_freq_hz and center_freq_hz >= 5e9) else 4.0
    if max_cal_sep_deg is not None and max_cal_sep_deg > sep_limit_deg:
        risk_factors.append("Calibrator-target separation is %.1f deg, beyond the ~%.1f deg the NRAO/VLBA "
                             "calibration guide recommends at this frequency for phase referencing: phase "
                             "solutions from a calibrator this far away represent a different patch of "
                             "sky/atmosphere." % (max_cal_sep_deg, sep_limit_deg))

    if target_min_elev is not None and target_min_elev < 30.0:
        airmass = 1.0 / np.sin(np.radians(target_min_elev)) if target_min_elev > 0 else float('inf')
        risk_factors.append("Target elevation drops to %.1f deg during the track (~%.2fx zenith airmass): "
                             "more troposphere/ionosphere path length, and more antenna shadowing risk."
                             % (target_min_elev, airmass))

    # Coherence-time-based cycle-time limit, same NRAO/VLBA guide: expect
    # ~120s coherence time at 300-700 MHz (ionosphere-dominated), ~300s
    # across 1-8 GHz. A calibrator cycle (cal-target-cal) shouldn't run
    # longer than this; an uninterrupted scan far past it has had that
    # long to drift with nothing correcting it in between.
    if center_freq_hz and center_freq_hz < 700e6:
        coherence_s = 120.0
    else:
        coherence_s = 300.0
    for f in targets:
        avg = field_avg_scan.get(f, 0)
        if avg > coherence_s:
            risk_factors.append("%s has long uninterrupted scans (avg %.0fs), %.1fx the ~%.0fs coherence-time-"
                                 "based cycle time the NRAO/VLBA guide suggests at this frequency: that's a lot "
                                 "of unattended integration time for phase to drift before any correction."
                                 % (f, avg, avg / coherence_s, coherence_s))

    if risk_factors:
        for r in risk_factors:
            out("  ! " + r)
    else:
        out("  No red flags found in the metadata checked here.")
    out("  This is a heuristic from metadata alone, not a measurement. The actual test is")
    out("  looking at the real data: phase vs time on the calibrator/target after the first")
    out("  calibration pass, or dynamic range and artifacts in a first dirty/CLEANed image.")

    # ---------------- checklist ----------------
    # Only for things this script can't compute a definite answer to
    # itself (things above already give a concrete number or verdict,
    # repeating them here as a vague yes/no checkbox would just be
    # confusing, e.g. see the separation and scan-length risk factors
    # above, which already say exactly what the situation is).
    out("\nCHECKLIST")
    out("  [ ] Calibrator names match IAU standard form or an entry in %s" % PHASE_CAL_FILE)
    out("  [ ] %s is present and up to date in this directory" % PHASE_CAL_FILE)
    out("  [ ] ref_ant in your pipeline config was actually present in THIS MS (check the Names list under")
    out("      ANTENNAS above; configs get reused across epochs and a real antenna can be down/unused")
    out("      for maintenance in any given observation, this isn't just a general telescope fact)")
    if cell_lo is not None:
        out("  [ ] when imaging the data, use a cell size in the ~%.2f-%.2f arcsec range suggested above"
            " (or %.2f as a single value)" % (cell_lo, cell_hi, cell_default))
    else:
        out("  [ ] when imaging the data, pick a cell size matching the band identified above"
            " (suggestion unavailable this run)")
    if callike_targets:
        out("  [ ] the coordinate-named 'target' field(s) flagged above have been double-checked")
    out("=" * 72)

    with open('read_listobs_report.txt', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print("\nReport saved to read_listobs_report.txt")


def main_fits(fitspath, telescope_arg):
    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 72)
    out(" read_listobs.py SUMMARY (FITS mode): %s" % fitspath)
    out("=" * 72)
    out("\n! FITS MODE: built directly from the FITS file, no MS was created and nothing")
    out("  was written to disk except this report. A raw FITS file doesn't have \"scans\" or")
    out("  \"flags\" yet (those only get built during MS import, e.g. CASA's importgmrt), so")
    out("  the following MS-mode sections are NOT available here and are skipped: SCAN TIMING")
    out("  PER FIELD, ELEVATION DURING TRACK, CURRENT FLAGGING STATUS, and the elevation-/")
    out("  scan-length-dependent self-cal risk factors. Everything else below (antennas,")
    out("  frequency/band setup, suggested cell size/FOV/image size, field classification and")
    out("  positions, polarization products, the risk factors that don't need scan timing) is")
    out("  read directly from the FITS file's own header and binary tables (AN/FQ/SU).")

    # ---------------- basic info + telescope profile ----------------
    out("\nBASIC INFO")
    try:
        obsname = fits_telescope_name(fitspath)
    except Exception as e:
        obsname = ""
        out("  Observatory (from FITS) : unavailable (%s)" % e)
    else:
        out("  Observatory (from FITS) : %s" % (obsname or "(unknown)"))

    try:
        sources = fits_source_list(fitspath)
    except Exception as e:
        sources = []
        out("  Total fields             : unavailable (%s)" % e)
    else:
        out("  Total fields             : %d" % len(sources))
    out("  Total scans              : not available in FITS mode (needs MS import)")

    telescope_key = resolve_telescope(telescope_arg, obsname, out)
    profile = TELESCOPE_PROFILES[telescope_key]

    # ---------------- observation time range / observing conditions ----------------
    qa_t = quanta()
    me_t = measures()
    out("\nOBSERVATION TIME RANGE / OBSERVING CONDITIONS")
    obs_t0 = obs_t1 = None
    inttim_arr = None
    try:
        obs_t0, obs_t1, inttim_arr = fits_time_range_mjd_sec(fitspath)
        d0 = qa_t.time(qa_t.quantity(obs_t0, 's'), form='ymd')[0]
        d1 = qa_t.time(qa_t.quantity(obs_t1, 's'), form='ymd')[0]
        out("  Observation window      : %s to %s (%.1f min)" % (d0, d1, (obs_t1 - obs_t0) / 60.0))
        out("  (Full span covered by the file's DATE parameters, not scan-by-scan timing,")
        out("  that needs MS import.)")
    except Exception as e:
        out("  Observation window: unavailable (%s)" % e)

    obspos = get_observatory_position(me_t, [obsname, profile.get('measures_site_name')])

    sun_conditions = []
    if obspos is not None and obs_t0 is not None:
        try:
            for label, t in (("start", obs_t0), ("mid", (obs_t0 + obs_t1) / 2.0), ("end", obs_t1)):
                sel = sun_elevation_deg(me_t, qa_t, obspos, t)
                cond = sun_condition(sel)
                sun_conditions.append((label, sel, cond))
                out("  Sun elevation at %-5s: %6.1f deg (%s)" % (label, sel, cond))
        except Exception as e:
            out("  Sun elevation unavailable (%s)" % e)
    else:
        out("  Could not resolve an observatory position for '%s' in CASA's Observatories"
            % (obsname or profile['label']))
        out("  table, so day/night/elevation checks are skipped.")

    # ---------------- antennas ----------------
    out("\nANTENNAS")
    ant_names = []
    ant_positions = None
    try:
        ant_names, ant_positions = fits_antenna_info(fitspath)
        out("  Count: %d" % len(ant_names))
        out("  Names: %s" % ', '.join(ant_names))
        odd = [a for a in ant_names if not profile['antenna_pattern'].match(a)]
        if odd:
            out("  ! Name(s) that don't match the expected %s naming pattern (informational, verify manually): %s"
                % (profile['label'], ', '.join(odd)))
    except Exception as e:
        out("  unavailable (%s)" % e)

    # ---------------- frequency setup ----------------
    out("\nFREQUENCY SETUP")
    center_freq_hz = None
    band_name = None
    cell_lo = cell_hi = cell_default = None
    try:
        subbands = fits_frequency_info(fitspath)
        out("  Number of IFs/subbands: %d" % len(subbands))
        all_subband_freqs = []
        for i, freqs in subbands:
            all_subband_freqs.append(freqs)
            out("  IF %d: %d channels, range %.3f-%.3f MHz" % (i, len(freqs), freqs.min() / 1e6, freqs.max() / 1e6))
        allf = np.concatenate(all_subband_freqs)
        center_freq_hz = (allf.min() + allf.max()) / 2.0
        band_name = guess_band(center_freq_hz, profile['band_ranges'])
        out("  Overall range: %.3f-%.3f MHz, center %.3f MHz" % (allf.min() / 1e6, allf.max() / 1e6, center_freq_hz / 1e6))
        out("  Likely band (%s): %s" % (profile['label'], band_name))

        if ant_positions is not None and len(ant_positions) >= 2:
            try:
                bmax = _max_baseline_from_positions(ant_positions)
                wavelength_m = 299792458.0 / center_freq_hz
                beam_arcsec = 206265.0 * wavelength_m / bmax
                cell_lo, cell_hi = beam_arcsec / 5.0, beam_arcsec / 3.0
                cell_default = beam_arcsec / 4.0
                out("  Max baseline: %.0f m -> synthesized beam ~%.2f arcsec (natural weighting, "
                    "no tapering)" % (bmax, beam_arcsec))
                out("  Suggested cell size: ~%.2f-%.2f arcsec (3-5 px/beam), %.2f arcsec if you just want one number (4 px/beam)"
                    % (cell_lo, cell_hi, cell_default))

                ddiam = profile.get('dish_diameter_m')
                if ddiam and cell_default:
                    fov_arcsec = 206265.0 * wavelength_m / ddiam
                    imsize_needed = fov_arcsec / cell_default
                    imsize_pow2 = next_pow2(imsize_needed)
                    out("  Dish diameter: %.1f m (known %s constant, not read from the FITS file) -> field of "
                        "view ~%.1f arcsec (%.2f arcmin)" % (ddiam, profile['label'], fov_arcsec, fov_arcsec / 60.0))
                    out("  Suggested image size: %d px (next power of 2 above the %.0f px the FOV/cell needs)"
                        % (imsize_pow2, imsize_needed))
            except Exception as e:
                out("  Suggested cell size / FOV / image size unavailable (%s)" % e)

        if profile['rfi_bands']:
            hits = [(lo, hi) for lo, hi in profile['rfi_bands'] if np.any((allf >= lo) & (allf <= hi))]
            if hits:
                out("  ! Known RFI-prone frequency range(s) fall inside your band:")
                for lo, hi in hits:
                    out("      %.1f-%.1f MHz" % (lo / 1e6, hi / 1e6))
            else:
                out("  No overlap with this profile's known persistent-RFI list.")
        else:
            out("  (No curated RFI-frequency list for %s in this script yet, not checked.)" % profile['label'])
    except Exception as e:
        out("  unavailable (%s)" % e)

    # polarization products
    pol_products = []
    try:
        pol_products = fits_polarization_products(fitspath)
        out("  Polarization products: %s" % ', '.join(pol_products))
    except Exception as e:
        out("  Polarization products: unavailable (%s)" % e)

    # integration time
    try:
        if inttim_arr is not None and len(inttim_arr):
            out("  Integration time (median INTTIM parameter): %.2f s" % float(np.median(inttim_arr)))
        else:
            out("  Integration time: no INTTIM parameter in this file, unavailable without MS import.")
    except Exception as e:
        out("  Integration time: unavailable (%s)" % e)

    # ---------------- field classification ----------------
    out("\nFIELD CLASSIFICATION")
    phase_cals, phase_cal_warn = load_phase_cals(PHASE_CAL_FILE)
    if phase_cal_warn:
        out("  " + phase_cal_warn)
    flux_cals = flatten_cals(profile['std_flux_cals'])
    fields = [s[0] for s in sources]
    field_radec = {s[0]: (s[1], s[2]) for s in sources}
    ampcals, pcals, targets = classify_fields(fields, phase_cals, flux_cals)
    out("  Flux/bandpass calibrators : %s" % (', '.join(ampcals) if ampcals else '(none found)'))
    out("  Phase calibrators         : %s" % (', '.join(pcals) if pcals else '(none found)'))
    out("  Targets                   : %s" % (', '.join(targets) if targets else '(none found)'))

    callike_targets = [t for t in targets if CALLIKE_NAME.match(t)]
    if callike_targets:
        out("  ! These 'targets' have coordinate-style names typical of calibrators,")
        out("    but aren't in %s or this profile's standard list. Verify they aren't" % PHASE_CAL_FILE)
        out("    meant to be a phase calibrator (a missing catalog entry will make a")
        out("    pipeline treat them as the science target):")
        out("      %s" % ', '.join(callike_targets))

    if ampcals and not pcals:
        out("  ! No separate phase calibrator found. Pipelines that fall back to using")
        out("    the amplitude/bandpass calibrator as the phase calibrator (e.g. CAPTURE)")
        out("    will pick the one with the most scans. Confirm that's what you want.")

    # ---------------- field positions ----------------
    out("\nFIELD POSITIONS")
    for fname, (ra_deg, dec_deg) in field_radec.items():
        out("  %-20s RA %s   Dec %s" % (fname, deg_to_hms(ra_deg), deg_to_dms(dec_deg)))

    max_cal_sep_deg = None
    if ampcals or pcals:
        out("\n  Calibrator-target separations:")
        for cal in ampcals + pcals:
            if cal not in field_radec:
                continue
            for tgt in targets:
                if tgt not in field_radec:
                    continue
                sep = angsep_deg(field_radec[cal][0], field_radec[cal][1],
                                  field_radec[tgt][0], field_radec[tgt][1])
                out("    %s <-> %s : %.1f deg" % (cal, tgt, sep))
                if max_cal_sep_deg is None or sep > max_cal_sep_deg:
                    max_cal_sep_deg = sep
        out("    (No universal cutoff here, but the larger this separation is, the less")
        out("     that calibrator's gain/phase solutions represent the sky above your")
        out("     target, especially if it's serving as your only calibrator.)")

    out("\nSCAN TIMING PER FIELD: not available in FITS mode (needs MS import).")
    out("\nELEVATION DURING TRACK: not available in FITS mode (needs scan timing from an MS).")
    out("\nCURRENT FLAGGING STATUS: not available in FITS mode (needs MS import; a freshly")
    out("  correlated file also has essentially nothing flagged yet anyway).")

    # ---------------- self-calibration risk factors (partial in FITS mode) ----------------
    out("\nSELF-CALIBRATION RISK FACTORS (heuristic, not a verdict, partial in FITS mode)")
    risk_factors = []

    if band_name:
        low_freq = ('Legacy' in band_name) or band_name in ('uGMRT Band 2', 'uGMRT Band 3', 'uGMRT Band 4',
                                                              '4-band', 'P-band')
        if low_freq:
            risk_factors.append("Low-frequency band (%s): ionospheric phase noise scales as 1/freq^2, "
                                 "usually the single biggest driver of selfcal need at these frequencies "
                                 "on %s." % (band_name, profile['label']))

    for label, sel, cond in sun_conditions:
        if 'twilight' in cond:
            risk_factors.append("Sun near horizon at track %s (%.1f deg): terminator crossings bring "
                                 "rapid ionospheric TEC gradients, often the roughest part of a day." % (label, sel))
        elif cond == 'daytime':
            risk_factors.append("Track %s is in local daytime (sun %.1f deg): daytime ionosphere is "
                                 "more disturbed than night." % (label, sel))

    if ampcals and not pcals:
        risk_factors.append("No separate phase calibrator: nothing is cycling in to track time-variable "
                             "phase during the target scans, selfcal is doing that job instead.")

    sep_limit_deg = 5.7 if (center_freq_hz and center_freq_hz >= 5e9) else 4.0
    if max_cal_sep_deg is not None and max_cal_sep_deg > sep_limit_deg:
        risk_factors.append("Calibrator-target separation is %.1f deg, beyond the ~%.1f deg the NRAO/VLBA "
                             "calibration guide recommends at this frequency for phase referencing: phase "
                             "solutions from a calibrator this far away represent a different patch of "
                             "sky/atmosphere." % (max_cal_sep_deg, sep_limit_deg))

    out("  (Elevation-during-track and long-uninterrupted-scan risk factors need scan timing")
    out("  from an MS, skipped here.)")

    if risk_factors:
        for r in risk_factors:
            out("  ! " + r)
    else:
        out("  No red flags found in what FITS mode could check.")
    out("  This is a heuristic from metadata alone, not a measurement. The actual test is")
    out("  looking at the real data: phase vs time on the calibrator/target after the first")
    out("  calibration pass, or dynamic range and artifacts in a first dirty/CLEANed image.")

    # ---------------- checklist ----------------
    out("\nCHECKLIST")
    out("  [ ] Calibrator names match IAU standard form or an entry in %s" % PHASE_CAL_FILE)
    out("  [ ] %s is present and up to date in this directory" % PHASE_CAL_FILE)
    out("  [ ] ref_ant in your pipeline config will actually be present once this becomes an MS")
    out("      (check the Names list under ANTENNAS above; MS import keeps the same antenna list)")
    if cell_lo is not None:
        out("  [ ] when imaging the data, use a cell size in the ~%.2f-%.2f arcsec range suggested above"
            " (or %.2f as a single value)" % (cell_lo, cell_hi, cell_default))
    else:
        out("  [ ] when imaging the data, pick a cell size matching the band identified above"
            " (suggestion unavailable this run)")
    if callike_targets:
        out("  [ ] the coordinate-named 'target' field(s) flagged above have been double-checked")
    out("=" * 72)

    with open('read_listobs_report.txt', 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print("\nReport saved to read_listobs_report.txt")


def main():
    input_arg, telescope_arg, perchannel_arg = parse_args()
    input_path = find_input(input_arg)
    kind = input_kind(input_path)
    if kind == 'fits':
        main_fits(input_path, telescope_arg)
    else:
        main_ms(input_path, telescope_arg, perchannel_arg)


if __name__ == '__main__':
    main()