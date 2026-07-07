# radio_astro_tools

Small, standalone tools for GMRT/uGMRT and VLA radio observations. Nothing here writes to a measurement set or requires a full pipeline install, each script is a single file you can drop into a working directory and run.

## Tools

- [`read_listobs.py`](#read_listobspy) — read-only pre-flight check on an existing measurement set (antennas, frequency setup, flagging, RFI, scan timing, self-cal risk factors)
- [`radio_sn_detectability.py`](#radio_sn_detectabilitypy) — proposal-planning sensitivity/detectability calculator, for before you have data
- [`sort_gmrt_files.py`](#sort_gmrt_filespy) — sorts raw GMRT/uGMRT downloads (lta-derived logs + FITS) into per-object folders, with a summary CSV

---

## `read_listobs.py`

### What it does

Runs against a measurement set you already have and prints a report covering:

- basic observation info and telescope identification (GMRT/uGMRT or VLA, auto-detected from the MS or set explicitly)
- observing conditions (sun elevation, day/night/twilight during the track)
- antennas present, flagged against the expected naming convention for the telescope
- frequency setup: spectral windows, likely receiver band, known RFI overlap, and suggested imaging cell size / field of view / image size
- polarization products and integration time
- field classification (flux/bandpass calibrator, phase calibrator, target), using a shared VLA calibrator-list file
- field positions and calibrator-target angular separation
- scan timing per field, including gaps and short-scan warnings
- elevation (and airmass) during each field's track
- current flagging status, with per-antenna and per-channel breakdowns to help spot a dead antenna or narrowband RFI
- self-calibration risk factors (low-frequency band, daytime/twilight track, no separate phase calibrator, large calibrator separation, low elevation, long uninterrupted scans), each backed by a real threshold (NRAO/VLBA calibration guide), not an arbitrary cutoff
- a closing checklist of things to verify by hand

It never modifies the MS. All CASA table/tool opens are read-only, and it never calls any flagging or calibration task that writes.

### How it operates

Everything except the flagging-status section is read from lightweight CASA metadata tables (`msmetadata`, small subtable reads via `table`), so it's fast even on a large MS. The flagging-status section calls `flagdata(mode='summary')`, which does read the actual FLAG column and is the slow part of the script, see the `--nochanflags` flag below.

A telescope "profile" (`TELESCOPE_PROFILES` in the script) supplies the instrument-specific knowledge: antenna naming pattern, receiver band edges, known RFI ranges, and standard calibrator names. GMRT/uGMRT and VLA/EVLA are included; adding another telescope means adding one more entry to that dict, nothing else in the script needs to change.

### Usage

```bash
casa -c read_listobs.py yourfile.ms
casa -c read_listobs.py --telescope gmrt yourfile.ms
casa -c read_listobs.py --telescope vla yourfile.ms
```

Or from the CASA prompt:

```python
MSNAME = 'yourfile.ms'
execfile('read_listobs.py')
```

If no filename is given and `MSNAME` is left as `None` at the top of the script, it looks for a single `*.ms`/`*.MS` directory in the current folder.

Requires a `vla-cals.list` file (the VLA calibrator manual list) in the working directory to correctly identify phase calibrators, used by GMRT/uGMRT observers for the same purpose.

Output is printed to the terminal and also written to `read_listobs_report.txt` in the current directory.

### Optional flags

| Flag | Effect |
|---|---|
| `--telescope gmrt` / `--telescope vla` | Override telescope auto-detection |
| `--nochanflags` | Skip the per-channel flagging breakdown (much faster; skips the one part of the script that reads real visibility data instead of metadata) |
| `--chanflags` | Force the per-channel breakdown on, overriding `CHECK_PERCHANNEL_FLAGS = False` in the script |

`CHECK_PERCHANNEL_FLAGS` (top of script) sets the default for the per-channel breakdown if you don't pass either flag.

---

## `radio_sn_detectability.py`

### What it does

A proposal-planning calculator: it doesn't read any measurement set. You give it a source's flux at some reference frequency, and it tells you the predicted thermal noise and signal-to-noise for a planned GMRT or VLA observation, so you can judge feasibility before applying for time.

There is no built-in example flux value. `--flux-ref-mjy`, `--freq-ref-mhz`, and `--tint-hr` are required arguments with no defaults, on purpose, so nothing here can be mistaken for a real result if you forget to set them.

### How it operates

Uses the standard NRAO array-sensitivity equation (NRAO VLA Observational Status Summary, Table 3.2.1, eq. 1):

```
sigma = SEFD_dish / (eta_c * sqrt(n_pol * N*(N-1) * t_int * bandwidth))
```

with a genuine per-dish SEFD and the actual antenna count `N`. Per-dish SEFDs for GMRT (Bands 2-5) come from NCRA's GTAC Cycle 47 status doc; for VLA (P through Q) from the NRAO VLA OSS. Flux at your target frequency is extrapolated from your reference flux with a power-law spectral index (`S ~ nu^alpha`), with a warning if you're extrapolating more than a factor of 3 in frequency.

Not included: confusion noise, weather, elevation-dependent gain loss, and the robust-vs-natural weighting penalty you'd see once you actually image. `--loss-fraction` (RFI/flagging margin) defaults to 0.0 rather than an assumed value, set it yourself if you want a margin.

### Usage

```bash
python radio_sn_detectability.py \
    --telescope gmrt --band "uGMRT Band 5" \
    --flux-ref-mjy <YOUR_FLUX_MJY> --freq-ref-mhz <YOUR_FREQ_MHZ> \
    --tint-hr <YOUR_HOURS>

python radio_sn_detectability.py --telescope vla --list-bands

python radio_sn_detectability.py --telescope gmrt --compare-bands \
    --flux-ref-mjy <YOUR_FLUX_MJY> --freq-ref-mhz <YOUR_FREQ_MHZ> \
    --tint-hr <YOUR_HOURS>
```

### Optional flags

| Flag | Default | Effect |
|---|---|---|
| `--alpha` | `-0.7` | Spectral index for flux extrapolation |
| `--target-freq-mhz` | band center | Frequency to predict the detection at |
| `--bandwidth-mhz` | full receiver bandwidth for the band | Usable processed bandwidth |
| `--n-ant` | GMRT: 26 (documented guaranteed minimum), VLA: 27 (nominal full array) | Antenna count assumed |
| `--n-pol` | `2` | Polarization products summed |
| `--snr-threshold` | `5.0` | S/N threshold for a "detection" verdict |
| `--loss-fraction` | `0.0` | Fraction of integration time assumed lost to RFI/flagging |
| `--list-bands` | — | List available bands and per-dish SEFD for `--telescope`, then exit |
| `--compare-bands` | — | Run the calculation across every band for `--telescope` and print a comparison table, instead of a single `--band` |

---

## `sort_gmrt_files.py`

### What it does

Sorts raw GMRT/uGMRT download files (`.obslog`, `.gvfits.log`, `.listscan.log`, `.listscan.plan`, and the converted `.fits`) into one folder per astronomical object: `<OBJECT_NAME>/gmrt_data/`. Every session also gets a row in a summary CSV: date of observation, band/RF frequency, calibrators seen, approximate time on target, project code/name, and which files were moved.

It doesn't touch CASA or a measurement set. It groups files by the proposal/date/band tag embedded in GMRT's raw filenames, works out which source in each session is a calibrator vs. the target, and files things accordingly.

### How it operates

Target identification reads the small text logs first: `.listscan.log` (or `.gvfits.log` as a fallback) already lists every source observed in a session, scan by scan, along with the true UT observation date. Anything matching a standard GMRT/VLA flux calibrator name or the coordinate-style naming convention used for phase calibrators is treated as a calibrator (extendable via `--extra-calibrators`); whatever's left is the target. Because these logs are a few KB, sorting normally happens without ever opening the FITS file (which can run tens of GB) or needing `astropy` installed. If no text logs are present for a session, it falls back to opening the FITS file directly and reading the AIPS SU source table, or the `OBJECT` header keyword, which does require `astropy`.

If a session's target can't be identified (only calibrators visible, or nothing to read), it gets its own placeholder folder (`object_1`, `object_2`, ...) instead of being lumped into one shared "unknown" bucket, since two unidentified sessions aren't necessarily the same real object.

Object folders are matched case-insensitively, so `2020bvc` and `2020BVC` from different sessions land in the same folder instead of creating a duplicate. Re-running the script against the same destination appends new rows to the CSV rather than overwriting it, and skips (with a warning) any file that would otherwise clobber an identically-named file already sitting in the destination.

### Usage

```bash
python sort_gmrt_files.py                                   # sorts the current folder in place
python sort_gmrt_files.py --src /path/to/downloads --dst /path/to/sorted
python sort_gmrt_files.py --dry-run                          # preview without touching any files
python sort_gmrt_files.py --copy                             # copy instead of move
```

### Optional flags

| Flag | Default | Effect |
|---|---|---|
| `--src` | current folder | Folder containing the raw downloaded files |
| `--dst` | same as `--src` | Destination base folder for the sorted object folders |
| `--copy` | off (moves) | Copy files instead of moving them |
| `--dry-run` | off | Show what would happen without touching any files |
| `--extra-calibrators` | — | Text file, one calibrator name per line, to extend the built-in calibrator list |
| `--csv` | `gmrt_sort_summary.csv` | Filename for the summary CSV (written into `--dst`) |
| `--unknown-prefix` | `object` | Prefix used to name sessions where no target could be identified |
