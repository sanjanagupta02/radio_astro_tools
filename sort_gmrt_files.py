#!/usr/bin/env python3
"""
sort_gmrt_files.py

Sorts raw uGMRT/GMRT download files (lta-derived .obslog, .gvfits.log,
.listscan.log, .listscan.plan, and the converted .fits file) into one
folder per astronomical object.

Object identification reads the small text logs first:
  - .listscan.log gives the observation date and a scan-by-scan list of
    source names (calibrators and target together), which is all that's
    needed to work out which source is the target.
  - .gvfits.log gives the same source list plus the RF frequency, used
    as a fallback if .listscan.log isn't present.
  - .obslog gives the project code/name and integration time, used to
    add a few extra columns to the summary CSV.
If none of those text logs are present, the script falls back to opening
the .fits file itself with astropy and reading its AIPS SU (source)
table, or its OBJECT header keyword. Because the text logs are only a
few KB (vs. the FITS file, which can be tens of GB), sorting normally
doesn't require astropy or the FITS file to be present at all.

Layout created:
    <dest>/<OBJECT_NAME>/gmrt_data/<all files for that session>

Also writes a summary CSV (default: gmrt_sort_summary.csv) with one row
per observing session (a session = one lta file and its derived
products): date of observation, band/RF frequency, detected object,
calibrators seen, approximate time on target, project code/name, and
which files were moved.

Usage:
    python sort_gmrt_files.py
        (no flags: sorts the files in the current folder in place, i.e.
        the <OBJECT_NAME>/gmrt_data/ folders are created right there)
    python sort_gmrt_files.py --src /path/to/downloads --dst /path/to/sorted
    python sort_gmrt_files.py --dry-run
    python sort_gmrt_files.py --copy

If --src is omitted, the current folder is used. If --dst is omitted, the
object folders are created directly inside --src itself, no extra
wrapping folder is added. Pass --dst explicitly if you want the sorted
output to go somewhere else.

astropy is only needed as a fallback, for sessions where none of the
text logs are present and the FITS file has to be opened directly. If
you need it:
    python3 -m venv gmrt_env
    source gmrt_env/bin/activate
    pip install astropy
"""

import argparse
import csv
import logging
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

try:
    from astropy.io import fits
except ImportError:
    fits = None

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger("sort_gmrt")

# ----------------------------------------------------------------------
# Known calibrator names. Anything matching these is NOT treated as the
# science target. Extend this list (or pass --extra-calibrators a text
# file, one name per line) if your projects use other secondary
# calibrators not covered here.
# ----------------------------------------------------------------------
STANDARD_FLUX_CALS = {
    "3C48", "3C147", "3C286", "3C468.1", "3C295",
    "0542+498", "1331+305", "0137+331",
}

# Phase/secondary calibrators at GMRT/VLA are named after their J2000-ish
# coordinates, e.g. "0834+555", "2246-121", "2052-474". This regex catches
# that naming convention generically.
CAL_COORD_PATTERN = re.compile(r"^\d{4}[+-]\d{2,4}$")

# Session key: proposal_code + sub-code + date + band, e.g.
#   47_021_25oct2024_b4
SESSION_RE = re.compile(
    r"(\d+_\d+_\d{2}[a-z]{3}\d{4}_b\d)", re.IGNORECASE
)

MONTHS = {
    "jan": "01", "feb": "02", "mar": "03", "apr": "04",
    "may": "05", "jun": "06", "jul": "07", "aug": "08",
    "sep": "09", "oct": "10", "nov": "11", "dec": "12",
}
DATE_IN_KEY_RE = re.compile(r"(\d{2})([a-z]{3})(\d{4})", re.IGNORECASE)


def load_extra_calibrators(path):
    if not path:
        return set()
    p = Path(path)
    if not p.exists():
        log.warning("Extra calibrators file %s not found, ignoring.", path)
        return set()
    names = set()
    for line in p.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            names.add(line.upper())
    return names


def is_calibrator(name, extra_cals):
    n = name.strip().upper()
    if n in STANDARD_FLUX_CALS:
        return True
    if n in extra_cals:
        return True
    if CAL_COORD_PATTERN.match(n):
        return True
    return False


def classify_sources(sources_ordered, extra_cals):
    targets, calibrators = [], []
    for name in sources_ordered:
        if is_calibrator(name, extra_cals):
            calibrators.append(name)
        else:
            targets.append(name)
    return targets, calibrators


# ----------------------------------------------------------------------
# Text log parsers (preferred: fast, no astropy/FITS access needed)
# ----------------------------------------------------------------------

SCAN_LINE_RE = re.compile(
    r"^Scan\s+\d+\s+(\S+)\s+(\d{2}:\d{2}:\d{2})\s+to\s+(\d{2}:\d{2}:\d{2})\s+(\d+)\s+recs",
    re.MULTILINE,
)


def parse_listscan_log(path):
    """Parse a .listscan.log file. Returns date_obs (YYYY-MM-DD or None)
    and a list of scan dicts: {source, start, end, recs}, in scan order."""
    text = path.read_text(errors="ignore")
    date_obs = None
    m = re.search(r"^DATE_OBS\s+(\S+)", text, re.MULTILINE)
    if m:
        date_obs = m.group(1).split("T")[0]
    scans = []
    for m in SCAN_LINE_RE.finditer(text):
        scans.append({
            "source": m.group(1), "start": m.group(2),
            "end": m.group(3), "recs": int(m.group(4)),
        })
    return date_obs, scans


GVFITS_SOURCE_RE = re.compile(r"^Source\s+\d+\s*:\s*(\S+)", re.MULTILINE)
GVFITS_FREQ_RE = re.compile(r"freq\s+([\d.]+)\s+[\d.]+\s+Hz")


def parse_gvfits_log(path):
    """Parse a .gvfits.log file. Returns (rf_freq_mhz, source_list)."""
    text = path.read_text(errors="ignore")
    rf_freq_mhz = None
    fm = GVFITS_FREQ_RE.search(text)
    if fm:
        try:
            rf_freq_mhz = round(float(fm.group(1)) / 1e6, 2)
        except ValueError:
            pass
    sources = GVFITS_SOURCE_RE.findall(text)
    return rf_freq_mhz, sources


def parse_obslog(path):
    """Parse a GTAC .obslog file for a few useful extra fields."""
    text = path.read_text(errors="ignore")
    info = {}
    m = re.search(r"Project Name\s*:\s*(.+)", text)
    if m:
        info["project_name"] = m.group(1).strip()
    m = re.search(r"Project Code\s*:\s*(\S+)", text)
    if m:
        info["project_code"] = m.group(1).strip()
    m = re.search(r"Start and End Time of Obs\.?\s*\(IST\)\s*:\s*(.+)", text)
    if m:
        info["obs_time_ist"] = m.group(1).strip()
    m = re.search(r"Integration Time\s*:\s*([\d.]+)\s*Sec", text)
    if m:
        try:
            info["integration_time_sec"] = float(m.group(1))
        except ValueError:
            pass
    return info


# ----------------------------------------------------------------------
# FITS fallback (only used when no text logs are present)
# ----------------------------------------------------------------------

def find_su_table(hdulist):
    """Return the AIPS SU (source) table HDU if present, else None."""
    for hdu in hdulist:
        extname = str(hdu.header.get("EXTNAME", "")).upper()
        if "SU" in extname and hasattr(hdu, "columns"):
            colnames = [c.upper() for c in hdu.columns.names]
            if any("SOURCE" in c for c in colnames):
                return hdu
    return None


def read_sources_from_fits(fits_path, extra_cals):
    """Fallback source/date reader used only when no .listscan.log or
    .gvfits.log is available. Returns (sources_ordered, date_obs, method)."""
    if fits is None:
        log.error("astropy is not installed, and no .listscan.log/.gvfits.log "
                   "was found for this session, so the target can't be identified. "
                   "Set up a venv and run: pip install astropy")
        return [], None, "unavailable (astropy missing, no logs found)"

    sources, date_obs, method = [], None, "unavailable"
    try:
        with fits.open(fits_path, ignore_missing_simple=True) as hdul:
            primary_header = hdul[0].header
            do = primary_header.get("DATE-OBS") or primary_header.get("DATOBS")
            if do:
                date_obs = str(do)[:10]

            su_hdu = find_su_table(hdul)
            if su_hdu is not None:
                col = next((c for c in su_hdu.columns.names if "SOURCE" in c.upper()), None)
                names = [str(x).strip() for x in su_hdu.data[col]]
                sources = sorted(set(n for n in names if n))
                method = "AIPS SU table (FITS fallback)"
            else:
                obj = primary_header.get("OBJECT")
                if obj:
                    sources = [str(obj).strip()]
                    method = "OBJECT header (FITS fallback)"
    except Exception as exc:
        log.error("Could not read %s: %s", fits_path, exc)

    return sources, date_obs, method


def date_from_session_key(session_key):
    m = DATE_IN_KEY_RE.search(session_key)
    if not m:
        return None
    dd, mon, yyyy = m.groups()
    mon_num = MONTHS.get(mon.lower())
    if not mon_num:
        return None
    return f"{yyyy}-{mon_num}-{dd}"


def band_from_session_key(session_key):
    m = re.search(r"_b(\d)$", session_key, re.IGNORECASE)
    return f"Band {m.group(1)}" if m else None


def group_files(src_dir):
    """Group all files in src_dir by session key. Files that don't match
    the expected naming pattern are returned separately under 'unmatched'."""
    groups = defaultdict(list)
    unmatched = []
    for f in sorted(Path(src_dir).iterdir()):
        if not f.is_file():
            continue
        m = SESSION_RE.search(f.name.lower())
        if m:
            groups[m.group(1).lower()].append(f)
        else:
            unmatched.append(f)
    return groups, unmatched


def find_by_suffix(files, suffix):
    for f in files:
        if f.name.lower().endswith(suffix):
            return f
    return None


def sanitize_folder_name(name):
    return re.sub(r"[^A-Za-z0-9_.+-]", "_", name.strip())


def existing_object_folders(dst_dir):
    """Map lowercased folder name -> actual on-disk folder name, for every
    object folder already sitting in dst_dir from a previous run. Used so
    a target found again later (possibly in different letter case, since
    the source name comes from whichever file/log happened to record it)
    lands in the SAME folder instead of creating a case-variant duplicate."""
    mapping = {}
    if dst_dir.exists():
        for d in dst_dir.iterdir():
            if d.is_dir():
                mapping[d.name.lower()] = d.name
    return mapping


def next_placeholder_name(dst_dir, prefix, used_this_run):
    """Pick the next free "<prefix>_N" name for a session whose target
    couldn't be identified automatically. Each unidentified session gets
    its own number (object_1, object_2, ...) rather than being dumped
    into one shared folder, since two unidentified sessions are not
    necessarily the same real object. Checks both names already handed
    out earlier in this run and folders left over from a previous run,
    so re-running the script won't collide with or overwrite old output."""
    i = 1
    while True:
        candidate = f"{prefix}_{i}"
        if candidate not in used_this_run and not (dst_dir / candidate).exists():
            used_this_run.add(candidate)
            return candidate
        i += 1


def identify_session(files, extra_cals):
    """Work out date, sources, and detection method for one session,
    preferring the small text logs over the (possibly absent/huge) FITS
    file. Returns a dict with everything the CSV/sorting logic needs."""
    listscan_log = find_by_suffix(files, ".listscan.log")
    gvfits_log = find_by_suffix(files, ".gvfits.log")
    obslog = find_by_suffix(files, ".obslog")
    fits_file = find_by_suffix(files, ".fits")

    sources_ordered = []
    scans_by_source = {}
    date_obs = None
    date_source = None
    method = "no source information found"
    rf_freq_mhz = None

    if listscan_log is not None:
        date_obs, scans = parse_listscan_log(listscan_log)
        if date_obs:
            date_source = "listscan.log"
        for sc in scans:
            if sc["source"] not in sources_ordered:
                sources_ordered.append(sc["source"])
            scans_by_source.setdefault(sc["source"], []).append(sc)
        if sources_ordered:
            method = "listscan.log"

    if gvfits_log is not None:
        freq, gv_sources = parse_gvfits_log(gvfits_log)
        rf_freq_mhz = freq
        if not sources_ordered and gv_sources:
            sources_ordered = list(dict.fromkeys(gv_sources))
            method = "gvfits.log"

    if not sources_ordered and fits_file is not None:
        fits_sources, fits_date, fits_method = read_sources_from_fits(fits_file, extra_cals)
        sources_ordered = fits_sources
        method = fits_method
        if fits_sources:
            method = fits_method
        if not date_obs and fits_date:
            date_obs = fits_date
            date_source = "FITS header"

    obslog_info = parse_obslog(obslog) if obslog is not None else {}

    targets, calibrators = classify_sources(sources_ordered, extra_cals)

    return {
        "sources": sources_ordered,
        "targets": targets,
        "calibrators": calibrators,
        "date_obs": date_obs,
        "date_source": date_source,
        "method": method,
        "rf_freq_mhz": rf_freq_mhz,
        "scans_by_source": scans_by_source,
        "obslog_info": obslog_info,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=".", help="Folder containing the raw downloaded files (default: current folder)")
    ap.add_argument("--dst", default=None, help="Destination base folder for the sorted object folders (default: sort in place inside --src)")
    ap.add_argument("--copy", action="store_true", help="Copy files instead of moving them")
    ap.add_argument("--dry-run", action="store_true", help="Show what would happen, without touching any files")
    ap.add_argument("--extra-calibrators", default=None, help="Text file, one calibrator name per line, to extend the built-in list")
    ap.add_argument("--csv", default="gmrt_sort_summary.csv", help="Filename for the summary CSV (written into --dst)")
    ap.add_argument("--unknown-prefix", default="object", help="Prefix used to name sessions where no target could be identified (object_1, object_2, ...)")
    args = ap.parse_args()

    src_dir = Path(args.src)
    if not src_dir.is_dir():
        log.error("Source folder %s does not exist.", src_dir)
        sys.exit(1)

    dst_dir = Path(args.dst) if args.dst else src_dir
    if not args.dst:
        log.info("No --dst given, sorting in place inside: %s", dst_dir.resolve())

    extra_cals = load_extra_calibrators(args.extra_calibrators)

    groups, unmatched = group_files(src_dir)
    if unmatched:
        log.warning("%d file(s) did not match the expected GMRT naming pattern and were left in place:", len(unmatched))
        for f in unmatched:
            log.warning("  %s", f.name)

    csv_rows = []
    action_word = "Would copy" if args.dry_run and args.copy else \
                  "Would move" if args.dry_run else \
                  "Copying" if args.copy else "Moving"
    used_placeholders = set()
    known_folders = existing_object_folders(dst_dir)

    for session_key, files in sorted(groups.items()):
        info = identify_session(files, extra_cals)

        date_obs = info["date_obs"] or date_from_session_key(session_key)
        date_source = info["date_source"] or "filename"
        band = band_from_session_key(session_key)

        was_identified = True
        if len(info["targets"]) == 1:
            object_name = info["targets"][0]
            note = ""
        elif len(info["targets"]) > 1:
            object_name = info["targets"][0]
            note = f"Multiple non-calibrator sources found ({', '.join(info['targets'])}); filed under the first one, please check."
        else:
            was_identified = False
            object_name = next_placeholder_name(dst_dir, args.unknown_prefix, used_placeholders)
            note = "No target identified automatically; check the raw files by hand and rename this folder."

        n_scans_on_target = ""
        approx_time_on_target_min = ""
        if was_identified:
            target_scans = info["scans_by_source"].get(object_name, [])
            if target_scans:
                n_scans_on_target = len(target_scans)
                total_recs = sum(sc["recs"] for sc in target_scans)
                itime = info["obslog_info"].get("integration_time_sec")
                if itime:
                    approx_time_on_target_min = round(total_recs * itime / 60, 1)

        folder_name = sanitize_folder_name(object_name)
        existing = known_folders.get(folder_name.lower())
        if existing and existing != folder_name:
            log.info("'%s' matches existing folder '%s' (different case), filing under the existing name.", folder_name, existing)
            folder_name = existing
        else:
            known_folders[folder_name.lower()] = folder_name
        dest_folder = dst_dir / folder_name / "gmrt_data"

        log.info("[%s] session=%s -> object=%s (%s), date=%s, band=%s, %d file(s)",
                  action_word, session_key, object_name, info["method"], date_obs, band, len(files))

        skipped_collisions = []
        if not args.dry_run:
            dest_folder.mkdir(parents=True, exist_ok=True)
            for f in files:
                dest_path = dest_folder / f.name
                if dest_path.exists():
                    log.warning("%s already exists in %s, leaving the copy already there in place (possible duplicate transfer).", f.name, dest_folder)
                    skipped_collisions.append(f.name)
                    continue
                if args.copy:
                    shutil.copy2(f, dest_path)
                else:
                    shutil.move(str(f), str(dest_path))
        if skipped_collisions:
            note = (note + " " if note else "") + f"{len(skipped_collisions)} file(s) already existed in the destination and were not overwritten: {', '.join(skipped_collisions)}."

        csv_rows.append({
            "object": object_name,
            "session_key": session_key,
            "date_obs": date_obs or "",
            "date_source": date_source,
            "band": band or "",
            "rf_freq_mhz": info["rf_freq_mhz"] if info["rf_freq_mhz"] is not None else "",
            "detection_method": info["method"],
            "calibrators_seen": "; ".join(info["calibrators"]),
            "all_sources_seen": "; ".join(info["sources"]),
            "n_scans_on_target": n_scans_on_target,
            "approx_time_on_target_min": approx_time_on_target_min,
            "project_code": info["obslog_info"].get("project_code", ""),
            "project_name": info["obslog_info"].get("project_name", ""),
            "obs_time_ist": info["obslog_info"].get("obs_time_ist", ""),
            "n_files": len(files),
            "files": "; ".join(f.name for f in files),
            "notes": note,
        })

    if not args.dry_run:
        dst_dir.mkdir(parents=True, exist_ok=True)
    csv_path = dst_dir / args.csv
    fieldnames = [
        "object", "session_key", "date_obs", "date_source", "band", "rf_freq_mhz",
        "detection_method", "calibrators_seen", "all_sources_seen",
        "n_scans_on_target", "approx_time_on_target_min",
        "project_code", "project_name", "obs_time_ist",
        "n_files", "files", "notes",
    ]
    if not args.dry_run:
        appending = csv_path.exists() and csv_path.stat().st_size > 0
        with open(csv_path, "a" if appending else "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            if not appending:
                writer.writeheader()
            for row in csv_rows:
                writer.writerow(row)
        log.info("%s %d row(s) %s %s", "Appended" if appending else "Wrote", len(csv_rows),
                  "to" if appending else "to new file", csv_path)
    else:
        already = " (an existing CSV there would be appended to, not overwritten)" if csv_path.exists() else ""
        log.info("Dry run: summary CSV would be written to %s (%d rows)%s", csv_path, len(csv_rows), already)

    n_flagged = sum(1 for r in csv_rows if r["notes"])
    if n_flagged:
        log.warning("%d session(s) flagged for manual review, see the 'notes' column in the CSV.", n_flagged)


if __name__ == "__main__":
    main()