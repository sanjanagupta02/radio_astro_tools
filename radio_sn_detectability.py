#!/usr/bin/env python
"""
radio_sn_detectability.py

Proposal-planning sensitivity calculator for GMRT/uGMRT and VLA continuum
observations of radio supernovae (or any point source with a known/assumed
flux at some reference frequency).

This is a PLANNING tool, separate from read_listobs.py (the MS preflight
script): it does not read any measurement set. It takes your own assumed
source flux and your planned observing setup, and tells you the predicted
noise level and signal-to-noise, so you can judge feasibility before you
apply for time or configure an observation.


FORMULA
This script instead uses the standard NRAO array-sensitivity equation 
(NRAO VLA Observational Status Summary, Table 3.2.1, eq. 1):

        sigma = SEFD_dish / (eta_c * sqrt(n_pol * N*(N-1) * t_int * bw))

    with a genuine PER-DISH SEFD and the actual antenna count N.

WHAT THIS DOES NOT INCLUDE (check these yourself for a real proposal)
    - Confusion noise. At low frequencies and/or long integrations, the
      confusion limit can dominate over thermal noise well before the
      thermal floor gets anywhere close to it. No confusion-limit table
      is included here (not verified for this script), check your band's
      literature value if you're planning a deep/low-frequency integration.
    - Weather, elevation-dependent gain/opacity loss, and the ~1.2x
      robust-vs-natural weighting penalty (NRAO guidance) that shows up
      once you actually image rather than just compute a thermal floor.
    - Real RFI/flagging loss: --loss-fraction defaults to 0.0 (no assumed
      loss) rather than a guessed number, set it yourself if you want a
      margin.

USAGE
    python radio_sn_detectability.py --telescope gmrt --band "uGMRT Band 5" \\
        --flux-ref-mjy <YOUR_FLUX_MJY> --freq-ref-mhz <YOUR_FREQ_MHZ> \\
        --tint-hr <YOUR_HOURS> [--alpha -0.7] [--bandwidth-mhz 200] \\
        [--n-ant 26] [--n-pol 2] [--snr-threshold 5.0] [--loss-fraction 0.0]

    python radio_sn_detectability.py --telescope vla --list-bands

    python radio_sn_detectability.py --telescope gmrt --compare-bands \\
        --flux-ref-mjy <YOUR_FLUX_MJY> --freq-ref-mhz <YOUR_FREQ_MHZ> \\
        --tint-hr <YOUR_HOURS>

ADDING A NEW TELESCOPE
    Add one entry to TELESCOPE_PROFILES with the same keys as the
    existing ones (same pattern as read_listobs.py's TELESCOPE_PROFILES).
    Leave 'deepest_rms_ujy' as None per band if you don't have a verified
    value, don't guess it.
"""

import argparse
import math
import sys

# ---------------------------------------------------------------------
# TELESCOPE PROFILES
#
# sefd_jy is the PER-DISH (single antenna) system-equivalent flux
# density, NOT reduced for array size, that reduction happens in
# predicted_rms_jy() using the actual antenna count you supply.
#
# deepest_rms_ujy (per band, where known) is the deepest RMS ever
# actually reported for that band in continuum mode, a sanity ceiling:
# your predicted floor should sit above this, not below it.
# ---------------------------------------------------------------------
TELESCOPE_PROFILES = {
    'gmrt': {
        'label': 'GMRT / uGMRT',
        'correlator_efficiency': 0.92,  # FX correlator (GWB), standard assumption, not GMRT-specific-cited
        'n_ant_full': 30,        # total antennas: 14 central square + 3 arms of ~6 (well known)
        'n_ant_guaranteed': 26,  # NCRA GTAC Cycle 47 status doc: "guarantee a minimum of 26 working antennas"
        'n_ant_default': 26,     # conservative default: the documented guaranteed minimum
        'source': 'NCRA "The GMRT: System Parameters and Current Status", GTAC Cycle 47 status doc, '
                   '15 June 2024, Table 1 (gmrt.ncra.tifr.res.in/doc/gtac_47_status_doc.pdf)',
        'bands': {
            # sefd_jy = Tsys/gain, using the midpoint of each band's quoted Tsys range in Table 1.
            'uGMRT Band 2': {'freq_range_mhz': (120, 250),   'sefd_jy': 1515.0, 'deepest_rms_ujy': 500.0},
            'uGMRT Band 3': {'freq_range_mhz': (250, 500),   'sefd_jy': 348.7,  'deepest_rms_ujy': 10.0},
            'uGMRT Band 4': {'freq_range_mhz': (550, 850),   'sefd_jy': 285.7,  'deepest_rms_ujy': 6.0},
            'uGMRT Band 5': {'freq_range_mhz': (1000, 1460), 'sefd_jy': 300.0,  'deepest_rms_ujy': 2.5},
        },
    },
    'vla': {
        'label': 'VLA / EVLA',
        'correlator_efficiency': 0.93,  # NRAO VLA OSS, 8-bit samplers
        'n_ant_full': 27,        # total antennas (well known)
        'n_ant_guaranteed': None,  # no NRAO-documented "guaranteed minimum" verified for this script;
                                   # in practice expect 1 antenna occasionally down for maintenance
        'n_ant_default': 27,
        'source': 'NRAO VLA Observational Status Summary, Table 3.2.1 '
                   '(science.nrao.edu/facilities/vla/docs/manuals/oss/performance/sensitivity)',
        'bands': {
            'P':  {'freq_range_mhz': (230, 470),      'sefd_jy': 2790.0, 'deepest_rms_ujy': None},
            'L':  {'freq_range_mhz': (1000, 2000),    'sefd_jy': 420.0,  'deepest_rms_ujy': None},
            'S':  {'freq_range_mhz': (2000, 4000),    'sefd_jy': 370.0,  'deepest_rms_ujy': None},
            'C':  {'freq_range_mhz': (4000, 8000),    'sefd_jy': 310.0,  'deepest_rms_ujy': None},
            'X':  {'freq_range_mhz': (8000, 12000),   'sefd_jy': 250.0,  'deepest_rms_ujy': None},
            'Ku': {'freq_range_mhz': (12000, 18000),  'sefd_jy': 320.0,  'deepest_rms_ujy': None},
            'K':  {'freq_range_mhz': (18000, 26500),  'sefd_jy': 500.0,  'deepest_rms_ujy': None},
            'Ka': {'freq_range_mhz': (26500, 40000),  'sefd_jy': 600.0,  'deepest_rms_ujy': None},
            'Q':  {'freq_range_mhz': (40000, 50000),  'sefd_jy': 1300.0, 'deepest_rms_ujy': None},
        },
    },
    # To add another telescope (e.g. MeerKAT, ATCA), copy one of the
    # blocks above and fill in what you actually know. Leave
    # 'deepest_rms_ujy': None if unverified rather than guessing.
}


def extrapolate_flux_mjy(flux_ref_mjy, freq_ref_mhz, freq_target_mhz, alpha):
    """Power-law flux extrapolation: S(nu) = S_ref * (nu/nu_ref)^alpha.
    Warns (doesn't block) when extrapolating more than a factor of 3 in
    frequency, a single spectral index is a weaker assumption that far
    from your reference point (spectral breaks, turnover, etc. can hide
    there)."""
    if freq_ref_mhz <= 0 or freq_target_mhz <= 0:
        raise ValueError("Frequencies must be positive.")
    ratio = freq_target_mhz / freq_ref_mhz
    if not (1.0 / 3.0 <= ratio <= 3.0):
        print("  WARNING: extrapolating by a factor of %.2fx in frequency (%.1f -> %.1f MHz). "
              "A single power-law spectral index is a weaker assumption this far from your "
              "reference frequency, watch for a real spectral break or turnover." % (ratio, freq_ref_mhz, freq_target_mhz))
    return flux_ref_mjy * (ratio ** alpha)


def predicted_rms_jy(sefd_dish_jy, eta_corr, n_ant, n_pol, bandwidth_hz, t_int_s):
    """NRAO array-sensitivity equation (NRAO VLA OSS, Table 3.2.1, eq. 1):
    sigma = SEFD_dish / (eta_c * sqrt(n_pol * N*(N-1) * t_int * bandwidth)).
    sefd_dish_jy is PER ANTENNA; N*(N-1) is what accounts for the array,
    the term GMRT_detectability.ipynb's original formula omitted. Returns
    None if inputs can't support a meaningful estimate."""
    if sefd_dish_jy is None or n_ant is None or n_ant < 2 or bandwidth_hz <= 0 or t_int_s <= 0:
        return None
    denom = eta_corr * math.sqrt(n_pol * n_ant * (n_ant - 1) * t_int_s * bandwidth_hz)
    if denom <= 0:
        return None
    return sefd_dish_jy / denom


def verdict(snr, snr_threshold):
    if snr >= 3 * snr_threshold:
        return "SAFELY DETECTABLE"
    elif snr >= snr_threshold:
        return "DETECTABLE (S/N above threshold, but not by a wide margin)"
    elif snr >= 0.5 * snr_threshold:
        return "MARGINAL (below your S/N threshold as set up; more time or a different band would help)"
    else:
        return "NOT DETECTABLE at this integration time"


def run_one_band(telescope_key, band_name, flux_ref_mjy, freq_ref_mhz, tint_hr, alpha,
                  target_freq_mhz=None, bandwidth_mhz=None, n_ant=None, n_pol=2,
                  snr_threshold=5.0, loss_fraction=0.0, verbose=True):
    """Compute predicted RMS, extrapolated flux, and S/N for one
    telescope/band combination. Returns a dict of the results; prints a
    human-readable report if verbose=True."""
    profile = TELESCOPE_PROFILES[telescope_key]
    band = profile['bands'][band_name]
    lo_mhz, hi_mhz = band['freq_range_mhz']
    f_target = target_freq_mhz if target_freq_mhz is not None else (lo_mhz + hi_mhz) / 2.0
    bw_hz = (bandwidth_mhz if bandwidth_mhz is not None else (hi_mhz - lo_mhz)) * 1e6
    n_ant_use = n_ant if n_ant is not None else profile['n_ant_default']
    eta_corr = profile['correlator_efficiency']
    sefd_dish = band['sefd_jy']

    t_int_s = tint_hr * 3600.0 * (1.0 - loss_fraction)

    flux_target_mjy = extrapolate_flux_mjy(flux_ref_mjy, freq_ref_mhz, f_target, alpha)
    rms_jy = predicted_rms_jy(sefd_dish, eta_corr, n_ant_use, n_pol, bw_hz, t_int_s)

    result = {
        'telescope': telescope_key, 'band': band_name, 'target_freq_mhz': f_target,
        'flux_target_mjy': flux_target_mjy, 'rms_jy': rms_jy, 'snr': None, 'verdict': None,
    }
    if rms_jy:
        snr = (flux_target_mjy / 1000.0) / rms_jy
        result['snr'] = snr
        result['verdict'] = verdict(snr, snr_threshold)

    if verbose:
        print("\n%s / %s (%.0f-%.0f MHz, predicting at %.1f MHz)" % (profile['label'], band_name, lo_mhz, hi_mhz, f_target))
        print("  Extrapolated flux at %.1f MHz: %.3f mJy (from %.3f mJy at %.1f MHz, alpha=%.2f)"
              % (f_target, flux_target_mjy, flux_ref_mjy, freq_ref_mhz, alpha))
        print("  Per-dish SEFD: %.0f Jy | correlator efficiency: %.2f | antennas: %d | pol products: %d"
              % (sefd_dish, eta_corr, n_ant_use, n_pol))
        print("  Integration time: %.2f hr usable (%.0f%% assumed loss) | bandwidth: %.1f MHz"
              % (t_int_s / 3600.0, loss_fraction * 100.0, bw_hz / 1e6))
        if f_target < 2000.0:
            print("  Note: no confusion-noise limit is included here. At this frequency, confusion can "
                  "dominate over thermal noise well before the number below, check your band's literature "
                  "confusion limit for a long integration.")
        if rms_jy:
            print("  Predicted RMS (theoretical, natural weighting): %.2f uJy/beam" % (rms_jy * 1e6))
            print("  Predicted S/N: %.1f  ->  %s (threshold S/N=%.1f)" % (snr, result['verdict'], snr_threshold))
            deepest = band.get('deepest_rms_ujy')
            if deepest is not None:
                if rms_jy * 1e6 < deepest:
                    print("  ! This predicted floor is BELOW the deepest RMS (%.1f uJy/beam) ever reported "
                          "for this band (%s). Treat it with real skepticism, actual noise is very unlikely "
                          "to beat the best ever achieved." % (deepest, profile['source']))
                else:
                    print("  For reference, the deepest RMS ever reported for this band is ~%.1f uJy/beam "
                          "(%s)." % (deepest, profile['source']))
            print("  Real image RMS will typically be worse than this theoretical floor: ~1.2x for robust")
            print("  vs. natural weighting (NRAO guidance), plus whatever confusion, weather, and elevation-")
            print("  dependent gain loss apply on top, none of which are included above.")
        else:
            print("  Predicted RMS unavailable (need >=2 antennas and positive bandwidth/time).")
    return result


def compare_bands(telescope_key, flux_ref_mjy, freq_ref_mhz, tint_hr, alpha, **kwargs):
    profile = TELESCOPE_PROFILES[telescope_key]
    print("=" * 72)
    print(" BAND COMPARISON: %s, %.3f mJy @ %.1f MHz (alpha=%.2f), %.2f hr"
          % (profile['label'], flux_ref_mjy, freq_ref_mhz, alpha, tint_hr))
    print("=" * 72)
    rows = []
    for band_name in profile['bands']:
        r = run_one_band(telescope_key, band_name, flux_ref_mjy, freq_ref_mhz, tint_hr, alpha,
                          verbose=True, **kwargs)
        rows.append(r)

    print("\nSUMMARY")
    print("  %-16s %10s %10s %8s  %s" % ("Band", "Freq(MHz)", "RMS(uJy)", "S/N", "Verdict"))
    for r in sorted(rows, key=lambda x: (x['snr'] is None, -(x['snr'] or 0))):
        rms_str = "%.2f" % (r['rms_jy'] * 1e6) if r['rms_jy'] else "n/a"
        snr_str = "%.1f" % r['snr'] if r['snr'] is not None else "n/a"
        print("  %-16s %10.1f %10s %8s  %s" % (r['band'], r['target_freq_mhz'], rms_str, snr_str, r['verdict'] or ''))

    best = max((r for r in rows if r['snr'] is not None), key=lambda x: x['snr'], default=None)
    if best:
        print("\n  Best predicted S/N: %s (%.1f)" % (best['band'], best['snr']))
    print("=" * 72)
    return rows


def build_parser():
    p = argparse.ArgumentParser(
        description="Proposal-planning sensitivity/detectability calculator for GMRT/uGMRT and VLA. "
                     "Reference flux, reference frequency, and integration time must be supplied "
                     "explicitly on the command line, there is no built-in example value for any of them.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--telescope', required=True, choices=sorted(TELESCOPE_PROFILES.keys()),
                    help="gmrt or vla")
    p.add_argument('--band', required=False,
                    help="Band name (e.g. 'uGMRT Band 5', 'L'). Required unless --list-bands or --compare-bands.")
    p.add_argument('--flux-ref-mjy', type=float, required=True,
                    help="Your own assumed/measured reference flux density, in mJy. REQUIRED, no default: "
                         "a light curve point, a detection at another frequency, or a model prediction "
                         "that you supply, not an example.")
    p.add_argument('--freq-ref-mhz', type=float, required=True,
                    help="The frequency (MHz) that --flux-ref-mjy applies at. REQUIRED, no default.")
    p.add_argument('--tint-hr', type=float, required=True,
                    help="Planned on-source integration time, in hours. REQUIRED, no default.")
    p.add_argument('--alpha', type=float, default=-0.7,
                    help="Spectral index, S ~ nu^alpha. Default -0.7 (a standard optically-thin synchrotron "
                         "planning assumption); override with your own source's measured or expected index "
                         "if you have one.")
    p.add_argument('--target-freq-mhz', type=float, default=None,
                    help="Frequency to predict the detection at. Default: the chosen band's center frequency.")
    p.add_argument('--bandwidth-mhz', type=float, default=None,
                    help="Usable processed bandwidth, in MHz. Default: the full receiver bandwidth for the "
                         "chosen band (an optimistic upper bound; real observers often process less than the "
                         "full band, set this to match your actual backend setup).")
    p.add_argument('--n-ant', type=int, default=None,
                    help="Antenna count to assume. Default: the telescope's documented guaranteed minimum "
                         "(GMRT: 26) or nominal full array (VLA: 27; see TELESCOPE_PROFILES comments for "
                         "caveats on that number).")
    p.add_argument('--n-pol', type=int, default=2,
                    help="Number of polarization products summed. Default 2 (both hands combined).")
    p.add_argument('--snr-threshold', type=float, default=5.0,
                    help="S/N threshold for a 'detection' verdict. Default 5.0 (standard radio-astronomy "
                         "convention for a secure detection).")
    p.add_argument('--loss-fraction', type=float, default=0.0,
                    help="Fraction of integration time to assume is lost to RFI/flagging, as a margin you "
                         "choose yourself. Default 0.0 (no assumed loss): guessing a number for you here "
                         "would be exactly the kind of invented assumption this tool is trying to avoid. "
                         "Set this explicitly (e.g. 0.2 for 20%%) if you want to build in a margin.")
    p.add_argument('--list-bands', action='store_true',
                    help="List the available bands (and their per-dish SEFD) for --telescope, then exit.")
    p.add_argument('--compare-bands', action='store_true',
                    help="Run the calculation across every band for --telescope and print a comparison "
                         "table, instead of a single --band.")
    return p


def main():
    args = build_parser().parse_args()
    profile = TELESCOPE_PROFILES[args.telescope]

    if args.list_bands:
        print("Bands available for %s (source: %s):" % (profile['label'], profile['source']))
        for name, b in profile['bands'].items():
            lo, hi = b['freq_range_mhz']
            print("  %-16s %8.0f-%-8.0f MHz   SEFD %8.1f Jy" % (name, lo, hi, b['sefd_jy']))
        return

    kwargs = dict(target_freq_mhz=args.target_freq_mhz, bandwidth_mhz=args.bandwidth_mhz,
                  n_ant=args.n_ant, n_pol=args.n_pol, snr_threshold=args.snr_threshold,
                  loss_fraction=args.loss_fraction)

    if args.compare_bands:
        compare_bands(args.telescope, args.flux_ref_mjy, args.freq_ref_mhz, args.tint_hr, args.alpha, **kwargs)
        return

    if not args.band:
        sys.exit("error: --band is required unless you pass --list-bands or --compare-bands. "
                  "Run with --list-bands to see the options for --telescope %s." % args.telescope)
    if args.band not in profile['bands']:
        sys.exit("error: unknown band '%s' for %s. Run with --list-bands to see the options."
                  % (args.band, profile['label']))

    run_one_band(args.telescope, args.band, args.flux_ref_mjy, args.freq_ref_mhz, args.tint_hr,
                 args.alpha, verbose=True, **kwargs)


if __name__ == '__main__':
    main()