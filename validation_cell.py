# ======================================================================
# UTMF v5.1 TCI/MCI NULL-TEST & PLOT GENERATOR  v3.0
#
# - Reads FULL_DETAILS_utmf_v*.json (UTMF v5.1 output)
# - Recomputes TCI/MCI null-tests *without* running the UTMF pipeline
# - Reconstructs structure figures + null-test figures from JSON only
# - Builds dendrogram + heatmap over ALL datasets (mean h(q))
#
# This cell is designed as a standalone reconstruction module:
# given a single FULL_DETAILS_utmf_v5.1_*.json file, all figures and
# scalar diagnostics used in the UTMF v5.1 paper can be regenerated.
# ======================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, glob, json, re
from datetime import datetime, timezone
from tqdm import tqdm
from scipy.stats import pearsonr
import scipy.signal
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from matplotlib.colors import ListedColormap

# ----------------------------------------------------------------------
# Global style settings (colors for real/synthetic and heatmaps)
# ----------------------------------------------------------------------
COLOR_REAL   = "#111111"   # dark grey / almost black
COLOR_SYNTH  = "#d62728"   # red
COLOR_OTHER1 = "#1f77b4"   # blue (spare color if needed)
HEATMAP_CMAP = "inferno"

# ======================================================================
# CLEAN CONSOLE OUTPUT HELPERS (PAPER-GRADE)
# ======================================================================

def banner(title: str):
    bar = "━" * 80
    print(f"\n{bar}\n{title.center(80)}\n{bar}")

def section(title: str):
    bar = "─" * 80
    print(f"\n{bar}\n▶ {title}\n{bar}")

def info(msg: str):
    print(f"  • {msg}")

def warn(msg: str):
    print(f"  [WARNING] {msg}")

def ok(msg: str):
    print(f"  [OK] {msg}")

def kv(label: str, value):
    print(f"  {label:<28}: {value}")

def list_preview(title: str, items, max_items: int = 10):
    print(f"  {title}: {len(items)} items")
    for x in items[:max_items]:
        print(f"     - {x}")
    if len(items) > max_items:
        print(f"     ... ({len(items)-max_items} more)")

# ----------------------------------------------------------------------
# Output directory for all figures / CSV used in the UTMF v5.1 paper
# ----------------------------------------------------------------------
OUTDIR = "/content/drive/MyDrive/Datasets_UTMF/UTMF_paper2_plots/"
os.makedirs(OUTDIR, exist_ok=True)

# ======================================================================
# LOAD LATEST FULL_DETAILS UTMF JSON (v5.1 OUTPUT)
# ======================================================================

json_paths = glob.glob(
    "/content/drive/MyDrive/Datasets_UTMF/UTMF_outputs/FULL_DETAILS_utmf_v*.json"
)
if not json_paths:
    raise FileNotFoundError("No FULL_DETAILS_utmf_v*.json found.")

# Most recent JSON by modification time
json_path = sorted(json_paths, key=os.path.getmtime)[-1]

with open(json_path, "r") as f:
    meta = json.load(f)

utmf_json = meta  # alias for clarity

# ======================================================================
# MFDFA IMPORTER
#
# JSON structure:
#   mfdfa_code = {
#       "version": "v5.2",
#       "sha256": "...",
#       "block": "<def polyfit_linear ... def jedi_mfdfa ...>"
#   }
#
# We restore the pure-Python (no @jit) implementation here so that
# jedi_mfdfa can be used for null-tests without re-running UTMF.
# ======================================================================

try:
    from numba import jit
except ImportError:
    # Dummy decorator if numba is not available
    def jit(*args, **kwargs):
        def wrap(f): 
            return f
        return wrap

if "mfdfa_code" not in meta or "block" not in meta["mfdfa_code"]:
    raise RuntimeError("JSON missing mfdfa_code['block'] — cannot restore jedi_mfdfa.")

mfdfa_block = meta["mfdfa_code"]["block"]

# Execute the stored source into an isolated namespace
_mfdfa_ns = {"np": np}
exec(mfdfa_block, _mfdfa_ns)

polyfit_linear = _mfdfa_ns["polyfit_linear"]
polyval_linear = _mfdfa_ns["polyval_linear"]
jedi_mfdfa     = _mfdfa_ns["jedi_mfdfa"]

# ======================================================================
# EXTRACT RUN TIMESTAMP (ROBUST: JSON FIELD → METADATA → FILENAME)
# ======================================================================

def extract_run_timestamp(json_path, meta_dict):
    """
    Robustly extract the original UTMF run_timestamp.

    Priority:
      1) meta['run_timestamp']
      2) meta['metadata']['run_timestamp']
      3) parse from filename: YYYYMMDD_HHMMSS
      4) 'Unknown run_timestamp'
    """
    # 1) direct field
    if isinstance(meta_dict, dict) and "run_timestamp" in meta_dict:
        rt = meta_dict["run_timestamp"]
        if isinstance(rt, str) and rt.strip():
            return rt

    # 2) nested metadata
    if isinstance(meta_dict, dict) and "metadata" in meta_dict:
        md = meta_dict["metadata"]
        if isinstance(md, dict) and "run_timestamp" in md:
            rt = md["run_timestamp"]
            if isinstance(rt, str) and rt.strip():
                return rt

    # 3) parse from filename
    fname = os.path.basename(json_path)
    m = re.search(r"(\d{8}_\d{6})", fname)
    if m:
        ts = m.group(1)
        try:
            dt = datetime.strptime(ts, "%Y%m%d_%H%M%S")
            return dt.replace(tzinfo=timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        except Exception:
            pass

    # 4) fallback
    return "Unknown run_timestamp"

RUN_TIMESTAMP = extract_run_timestamp(json_path, meta)
NOW_TS = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
PLOT_TIMESTAMP = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

banner("UTMF v5.1 TCI/MCI NULL-TEST & PLOT GENERATOR  v3.0")
section("Loading FULL_DETAILS JSON Metadata")

kv("Selected JSON file", json_path)
kv("UTMF run_timestamp", RUN_TIMESTAMP)
kv("Plot generation time", PLOT_TIMESTAMP)
print("  Loaded jedi_mfdfa from JSON (pure Python, no JIT).")

# ----------------------------------------------------------------------
# Helper: add a footer with JSON run timestamp to each figure
# ----------------------------------------------------------------------
def add_footer_timestamp(fig):
    txt = (
        f"UTMF v5.1 TCI/MCI null-test based on FULL_DETAILS JSON "
        f"(run_timestamp: {RUN_TIMESTAMP})"
    )
    fig.text(0.5, 0.01, txt, ha="center", va="bottom", fontsize=8.5, alpha=0.7)

# Unified show+save helper
def save_fig(fig, name, show=True):
    if show:
        plt.show()
    path = os.path.join(OUTDIR, name)
    fig.savefig(path, dpi=500, bbox_inches="tight")
    print(f"[FIG SAVED] {path}")
    plt.close(fig)

# ======================================================================
# DOMAIN MAPPING (FULL VERSION) — SHARED WITH FAT
# ======================================================================

ELEMENTS = [
    "AC", "AG", "AL", "AR", "B", "BA", "BE", "BI", "BR", "C",
    "CA", "CD", "CE", "CL", "CO", "CR", "CS", "CU", "DY", "EU",
    "F", "FE", "H", "HE", "HF", "HG", "I", "IN", "IR", "K",
    "KR", "LA", "LI", "MG", "MN", "MO", "N", "NA", "NB", "ND",
    "NE", "NI", "O", "OS", "P", "PR", "PT", "RB", "RE", "RH",
    "S", "SC", "SI", "SN", "SR", "TA", "TI", "TM", "V", "W",
    "XE", "Y", "ZR"
]

def get_domain(name: str) -> str:
    """
    Map a dataset name to a coarse physical domain label.
    This is used for coloring/clustering in the plots.
    """
    name_u = name.upper()

    # Explicit domains first
    if "QRNG" in name_u:
        return "NIST-QRNG"
    if "LIGO" in name_u:
        return "LIGO-L1"
    if "PLANCK" in name_u or "CMB" in name_u:
        return "Planck"
    if "NANOGRAV" in name_u:
        return "NANOGrav"
    if "GAIA" in name_u:
        return "Gaia"
    if "DESI" in name_u:
        return "DESI"
    if "CERN" in name_u or "2LEPTON" in name_u:
        return "CERN"

    # All atomic spectra lines from the NIST database
    if name_u.startswith("NIST_3_"):
        return "NIST-Elements"

    # Legacy element-fallback (harmless to keep)
    for el in ELEMENTS:
        patterns = [
            f"_{el}_", f"-{el}-", f"_{el}-", f"-{el}_",
            f"_{el}",  f"-{el}", f"{el}_", f"{el}-",
            f" {el} ", f" {el}_", f"_{el} ", el
        ]
        if any(p in name_u for p in patterns):
            return "NIST-Elements"

    return "Other"

domain_colors = {
    "LIGO-L1":        "#1f77b4",
    "Planck":         "#9467bd",
    "DESI":           "#ff7f0e",
    "CERN":           "#8c564b",
    "NIST-Elements":  "#d62728",
    "NIST-QRNG":      "#bcbd22",
    "NANOGrav":       "#2ca02c",
    "Gaia":           "#17becf",
    "Other":          "#7f7f7f",
}

def color_for_domain(dom: str) -> str:
    return domain_colors.get(dom, "#7f7f7f")

# ======================================================================
# q-GRID HELPERS (FOR h(q) CURVES FROM JSON)
# ======================================================================

def infer_q_grid_from_config(utmf_json, len_q: int | None = None):
    """
    Try to recover the q-grid (the x-axis for h(q)) from the JSON config.

    Search order:
      1) config_snapshot keys like 'hq_q', 'q_list', 'mfdfa_q'
      2) fallback: linspace(-8, 8, len_q) if length is known
    """
    cfg = utmf_json.get("config_snapshot", {}) or {}

    for key in ["hq_q", "q_list", "mfdfa_q", "mfdfa_q_list"]:
        if key in cfg:
            q = np.asarray(cfg[key], dtype=float)
            if q.size > 0:
                return q

    if len_q is not None and len_q > 0:
        return np.linspace(-8, 8, len_q)

    return None

def extract_hq_curves_from_json(utmf_json):
    """
    Extract per-dataset mean h(q) curves from the 'subsets_processed'
    block of the FULL_DETAILS JSON.

    Returns
    -------
    curves : dict
        Mapping name -> {"q": q_grid, "h": mean_hq}
    """
    subsets = utmf_json.get("subsets_processed", {}) or {}
    curves = {}

    if not subsets:
        print("[extract_hq_curves_from_json] No 'subsets_processed' found.")
        return curves

    # Determine q-length from the first dataset
    first_ds = next(iter(subsets.values()))
    first_hq = first_ds.get("hq_values")
    len_q = len(first_hq[0]) if first_hq and len(first_hq) > 0 else None

    q_grid = infer_q_grid_from_config(utmf_json, len_q=len_q)
    if q_grid is None:
        print("[extract_hq_curves_from_json] Could not infer q-grid.")
        return curves

    for ds_name, info in subsets.items():
        hq_values = info.get("hq_values")
        if not hq_values:
            continue

        hq_arr = np.asarray(hq_values, dtype=float)  # shape: (n_subsets, n_q)
        if hq_arr.ndim != 2 or hq_arr.shape[1] != q_grid.size:
            print(f"[extract_hq_curves_from_json] Skip {ds_name}: invalid hq-shape {hq_arr.shape}")
            continue

        mean_curve = np.nanmean(hq_arr, axis=0)

        # Apply q=0 repair for the mean curve
        mean_curve = fix_hq_q0(mean_curve, q_grid)

        curves[ds_name] = {
            "q": q_grid,
            "h": mean_curve,
        }

    return curves

# ======================================================================
# 1. DATA LOADING FOR TCI (tci_timeseries_snapshot → tci_datasets)
# ======================================================================

FIXED_LEN = 6144  # TCI analysis uses a fixed window length
print(f"\n-- Fixed TCI length used: {FIXED_LEN} samples --")
tci_datasets = {}

tci_snapshot = meta.get("tci_timeseries_snapshot", {})
if not isinstance(tci_snapshot, dict):
    warn("meta['tci_timeseries_snapshot'] is not a dict; TCI reconstruction may fail.")

def extract_series(obj):
    """
    Robust extractor for any JSON-encoded time-series structure.

    It tries several common patterns:
      - raw lists / arrays
      - dicts with keys like 'data', 'series', 'raw', ...
    """
    if isinstance(obj, (list, tuple, np.ndarray)):
        try:
            return np.asarray(obj, dtype=float)
        except Exception:
            return None

    if isinstance(obj, dict):
        # Preferred keys
        for key in ["series", "raw", "cleaned", "denoised", "data", "time_series"]:
            if key in obj:
                try:
                    return np.asarray(obj[key], dtype=float)
                except Exception:
                    pass

        # Fallback: first array-like value
        for _, v in obj.items():
            if isinstance(v, (list, tuple, np.ndarray)):
                try:
                    return np.asarray(v, dtype=float)
                except Exception:
                    pass

    return None

for name, obj in tci_snapshot.items():
    arr = extract_series(obj)
    if arr is None:
        print(f"[TCI WARNING] Could not extract usable series for dataset: {name}")
        continue

    if arr.size < FIXED_LEN:
        print(f"[TCI INFO] Skipping {name} (series too short: {arr.size} < {FIXED_LEN})")
        continue

    tci_datasets[name] = arr

section("TCI Data Extraction")
kv("TCI fixed length", FIXED_LEN)
kv("Usable TCI series", len(tci_datasets))
list_preview("Datasets included for TCI", list(tci_datasets.keys()))

# ======================================================================
# 2. DATA LOADING FOR MCI (subsets_processed → results_all-like)
# ======================================================================

results_all = {}

subsets = meta.get("subsets_processed", {})
if not isinstance(subsets, dict):
    warn("meta['subsets_processed'] is not a dict; MCI reconstruction may fail.")

for ds_name, ds_info in subsets.items():
    if not isinstance(ds_info, dict):
        continue

    if "hq_values" in ds_info:
        hq_list = ds_info["hq_values"]
    elif "h_q_values" in ds_info:
        hq_list = ds_info["h_q_values"]
    else:
        continue

    if not isinstance(hq_list, (list, tuple)) or len(hq_list) == 0:
        continue

    hq_list_np = [np.asarray(h, dtype=float) for h in hq_list]
    results_all[ds_name] = {"hq_values": hq_list_np}

# ======================================================================
# 3. TCI NULL-TEST – REAL VS SYNTHETIC h(q) COHERENCE
# ======================================================================

real_series = []
real_names = []

for name, series in tci_datasets.items():
    if series is None or len(series) < FIXED_LEN:
        continue

    s = np.asarray(series[:FIXED_LEN], dtype=np.float64)
    lname = name.lower()

    # Dataset-dependent preprocessing (must match UTMF conventions)
    if "ligo" in lname:
        s = s * 1e16
        s = scipy.signal.savgol_filter(s, window_length=7, polyorder=1, mode="nearest")
    elif "cmb" in lname or "planck" in lname:
        s = s / np.std(s)
    elif "qrng" in lname:
        s = (s - np.mean(s)) / (np.std(s) + 1e-10)
    else:
        s = s / np.std(s)

    real_series.append(s)
    real_names.append(name)

section("TCI Null-Test: Real-Series Preparation")
kv("Real TCI series (post-filter)", len(real_series))
list_preview("Real TCI datasets", real_names)

# ----------------------------------------------------------------------
# Synthetic TCI series (null models)
# ----------------------------------------------------------------------
synth_series = []
kinds = ["white", "ar1", "pink", "fgn", "phase_surrogate"]

for kind_idx, kind in enumerate(kinds):
    for rep in range(100):  # 5 × 100 = 500 synthetic series
        seed = rep + kind_idx * 100
        rng = np.random.default_rng(seed)

        if kind == "white":
            x = rng.normal(0, 1, FIXED_LEN)

        elif kind == "ar1":
            x = np.cumsum(rng.normal(0, 0.3, FIXED_LEN)) + rng.normal(0, 1, FIXED_LEN)

        elif kind == "pink":
            f = np.fft.rfftfreq(FIXED_LEN)
            f[0] = 1
            x = np.fft.irfft(
                np.fft.rfft(rng.normal(0, 1, FIXED_LEN)) / np.sqrt(f + 1e-12),
                n=FIXED_LEN,
            )

        elif kind == "fgn":
            f = np.fft.rfftfreq(FIXED_LEN)
            f[0] = 1
            x = np.fft.irfft(
                np.fft.rfft(rng.normal(0, 1, FIXED_LEN)) * (f + 1e-12) ** (-0.4),
                n=FIXED_LEN,
            )

        elif kind == "phase_surrogate":
            fft = np.fft.rfft(rng.normal(0, 1, FIXED_LEN))
            mag = np.abs(fft)
            phase = rng.uniform(0, 2 * np.pi, len(fft))
            x = np.fft.irfft(mag * np.exp(1j * phase), n=FIXED_LEN)

        else:
            x = rng.normal(0, 1, FIXED_LEN)

        x = x / np.std(x)
        synth_series.append(x)

info(f"Synthetic TCI series generated: {len(synth_series)} total\n")

# ======================================================================
# --- h(q) computation: JSON-driven, stable, with q=0 repair ---
# ======================================================================

q_values = np.asarray(meta["config_snapshot"]["mfdfa"]["q_values"], dtype=float)
# Mask for plotting (optionally exclude q = 0 from axes if desired)
Q_MASK = q_values != 0
Q_PLOT = q_values[Q_MASK]

def fix_hq_q0(hq, q_values):
    """
    Smooth the artificial kink around q = 0:

      - locate q ≈ 0 using np.isclose (instead of exact equality)
      - set h(0) to the average of the nearest valid neighbours
        (or to the nearest side if only one is available)
    """
    q_values = np.asarray(q_values, dtype=float)
    hq = np.asarray(hq, dtype=float).copy()

    idx0 = np.where(np.isclose(q_values, 0.0))[0]
    if idx0.size == 0:
        return hq

    k = idx0[0]

    left = hq[:k][np.isfinite(hq[:k])]
    right = hq[k+1:][np.isfinite(hq[k+1:])]

    if left.size and right.size:
        hq[k] = 0.5 * (left[-1] + right[0])
    elif left.size:
        hq[k] = left[-1]
    elif right.size:
        hq[k] = right[0]
    else:
        # If everything is NaN, leave as-is
        pass

    return hq

def get_hq(s):
    """
    Compute h(q) for a single time series `s` using jedi_mfdfa:

      - logspace scales from 8 to len(s)/8
      - q-grid taken from the JSON config
      - q=0 repair applied
      - interpolation over remaining NaNs where possible
    """
    scales = np.logspace(
        np.log10(8),
        np.log10(len(s) // 8),
        15
    ).astype(int)
    scales = scales[scales > 4]

    _, hq, _, _, _ = jedi_mfdfa(s, scales, q_values)
    hq = np.asarray(hq, dtype=float)

    # q=0 smoothing
    hq = fix_hq_q0(hq, q_values)

    # Fill remaining NaNs by interpolation over q
    if np.any(np.isnan(hq)):
        valid = np.isfinite(hq)
        if np.sum(valid) >= 3:
            hq = np.interp(q_values, q_values[valid], hq[valid])

    return hq

# ======================================================================
# Compute h(q) for real & synthetic TCI series
# ======================================================================

real_hq = []
for i, s in enumerate(tqdm(real_series, desc="Real TCI h(q)")):
    hq = get_hq(s)
    if hq is None:
        print(f"[TCI WARNING] h(q) for real[{i}] is None — skipped.")
        continue
    hq = np.asarray(hq, dtype=float)

    if np.all(np.isnan(hq)):
        print(f"[TCI WARNING] h(q) for real[{i}] is all-NaN — skipped.")
        continue
    if hq.ndim != 1:
        print(f"[TCI WARNING] h(q) for real[{i}] malformed (shape {hq.shape}) — skipped.")
        continue
    if len(hq) != len(q_values):
        print(
            f"[TCI WARNING] h(q) for real[{i}] wrong length ({len(hq)}), "
            f"expected {len(q_values)} — skipped."
        )
        continue

    real_hq.append(hq)

synth_hq = []
for i, s in enumerate(tqdm(synth_series, desc="Synthetic TCI h(q)")):
    hq = get_hq(s)
    if hq is None:
        print(f"[TCI WARNING] h(q) for synth[{i}] is None — skipped.")
        continue
    hq = np.asarray(hq, dtype=float)

    if np.all(np.isnan(hq)):
        print(f"[TCI WARNING] h(q) for synth[{i}] is all-NaN — skipped.")
        continue
    if hq.ndim != 1:
        print(f"[TCI WARNING] h(q) for synth[{i}] malformed (shape {hq.shape}) — skipped.")
        continue
    if len(hq) != len(q_values):
        print(
            f"[TCI WARNING] h(q) for synth[{i}] wrong length ({len(hq)}), "
            f"expected {len(q_values)} — skipped."
        )
        continue

    synth_hq.append(hq)

def get_all_pairwise_corrs(hq_list):
    """
    Compute all pairwise |corr| across a list of h(q)-vectors.
    Returns a flat list of absolute Pearson correlations.
    """
    corrs = []
    for i in range(len(hq_list)):
        for j in range(i + 1, len(hq_list)):
            c, _ = pearsonr(hq_list[i], hq_list[j])
            if np.isfinite(c):
                corrs.append(abs(c))
    return corrs if corrs else [np.nan]

real_pairwise_tci  = get_all_pairwise_corrs(real_hq)
synth_pairwise_tci = get_all_pairwise_corrs(synth_hq)

# ======================================================================
# 4. MCI NULL-TEST – DATASET-LEVEL COHERENCE OVER h(q) MEANS
# ======================================================================

section("MCI Data Extraction")

real_mci = []
real_names_mci = []

for name, res in results_all.items():
    if not isinstance(res, dict) or "hq_values" not in res:
        continue

    hq_list = res["hq_values"]
    if not isinstance(hq_list, (list, tuple)) or len(hq_list) == 0:
        continue

    hq_arr = np.asarray(hq_list, dtype=float)
    if hq_arr.ndim != 2:
        continue

    if np.all(np.isnan(hq_arr)):
        print(f"[MCI WARNING] Dataset '{name}' has only NaN h(q) entries — skipped.")
        continue

    col_has_value = np.any(np.isfinite(hq_arr), axis=0)
    if not np.any(col_has_value):
        print(f"[MCI WARNING] Dataset '{name}' has no finite h(q) columns — skipped.")
        continue

    hq_reduced = hq_arr[:, col_has_value]
    if hq_reduced.size == 0:
        print(f"[MCI WARNING] Dataset '{name}' has no usable h(q) columns — skipped.")
        continue

    # Mean over q's that are valid for this dataset
    mean_valid = np.nanmean(hq_reduced, axis=0)
    if np.all(np.isnan(mean_valid)):
        print(f"[MCI WARNING] Dataset '{name}' has NaN mean — skipped.")
        continue

    # Corresponding q-grid for these columns
    q_valid = q_values[col_has_value]

    # Interpolate to full q-grid so that all datasets share the same length
    try:
        mean_full = np.interp(q_values, q_valid, mean_valid)
    except Exception:
        print(f"[MCI WARNING] Interpolation failed for '{name}' — skipped.")
        continue

    # Smooth around q = 0
    mean_full = fix_hq_q0(mean_full, q_values)

    real_mci.append(mean_full)
    real_names_mci.append(name)

kv("Valid MCI datasets", len(real_mci))
list_preview("Datasets included for MCI", real_names_mci, max_items=10)
print()

# ----------------------------------------------------------------------
# Synthetic MCI curves – generated from a template curve
# ----------------------------------------------------------------------

def generate_synthetic(template, kind, seed=42):
    """
    Generate synthetic 1D series with the same length as `template`
    for MCI null-tests. Different `kind` values correspond to:
      - 'white'  : white Gaussian noise
      - 'ar1'    : AR(1)-like process
      - 'pink'   : 1/f noise
      - 'fgn'    : fractional Gaussian noise
      - 'phase_surrogate' : phase-scrambled surrogate of template
    """
    rng = np.random.default_rng(seed)
    n = len(template)

    if kind == "white":
        return rng.normal(0, 1, n)

    if kind == "ar1":
        out = np.zeros(n)
        out[0] = template[0]
        for i in range(1, n):
            out[i] = 0.8 * out[i - 1] + rng.normal(0, 1)
        return out

    if kind == "pink":
        w = rng.normal(0, 1, n)
        fft = np.fft.rfft(w)
        k = np.fft.rfftfreq(n)
        k[0] = 1
        fft /= np.sqrt(k)
        return np.fft.irfft(fft, n=n)

    if kind == "fgn":
        f = np.fft.rfftfreq(n)
        f[0] = 1
        psd = np.abs(f) ** (-0.8)
        w = rng.normal(0, 1, len(f)) + 1j * rng.normal(0, 1, len(f))
        return np.fft.irfft(w * np.sqrt(psd), n=n)

    if kind == "phase_surrogate":
        fft = np.fft.rfft(template)
        mag = np.abs(fft)
        phase = rng.uniform(0, 2 * np.pi, len(fft))
        return np.fft.irfft(mag * np.exp(1j * phase), n=n)

    # Default: white noise
    return rng.normal(0, 1, n)

if real_mci:
    template = real_mci[0]
else:
    template = np.zeros(1000)

q_values_cfg = np.asarray(meta["config_snapshot"]["mfdfa"]["q_values"], dtype=float)
scales_cfg = np.logspace(np.log10(8), np.log10(len(template) // 8), 15).astype(int)

synth_kinds = ["white", "ar1", "pink", "fgn", "phase_surrogate"]
synth_mci = []
synth_names = []

for kind in synth_kinds:
    for rep in range(100):
        series = generate_synthetic(template, kind, seed=rep)
        _, hq, _, _, _ = jedi_mfdfa(series, scales_cfg, q_values_cfg, detrend_order=0)
        hq = np.where(np.isfinite(hq), hq, 0.5)
        synth_mci.append(hq)
        synth_names.append(f"SYN_{kind}_{rep+1}")

section("MCI Synthetic Generation")
kv("Synthetic MCI curves", len(synth_mci))
ok("MCI synthetic generation complete")

# ======================================================================
# MCI SCALAR: COHERENCE INDEX PER DATASET (REAL + SYNTHETIC)
# ======================================================================

def dataset_mci_from_curves(hq_list):
    """
    Given a list of mean h(q) curves (one per dataset), compute:

        MCI_i = mean_j |corr(hq_i, hq_j)|,  j != i

    for each dataset i, plus the full correlation matrix.

    - Filters invalid/NaN curves.
    - Suppresses numerical warnings inside corrcoef.
    - Avoids 'mean of empty slice' by explicit checks.
    """
    curves = []
    for h in hq_list:
        arr = np.asarray(h, dtype=float)
        if arr.ndim != 1:
            continue
        if not np.any(np.isfinite(arr)):
            continue
        curves.append(arr)

    n = len(curves)
    if n < 2:
        return [np.nan] * n, np.full((n, n), np.nan)

    H = np.vstack(curves)  # shape: (N_datasets, N_q)

    # Corrcoef with suppressed warnings
    with np.errstate(divide='ignore', invalid='ignore'):
        C = np.corrcoef(H)

    C = np.abs(C)
    np.fill_diagonal(C, np.nan)

    mci_vals = np.full(n, np.nan)
    for i in range(n):
        row = C[i, :]
        finite = np.isfinite(row)
        if np.sum(finite) == 0:
            continue
        mci_vals[i] = np.nanmean(row[finite])

    return mci_vals.tolist(), C

# Compute real + synthetic MCI scalars
real_mci_scalar,  real_mci_corrmat  = dataset_mci_from_curves(real_mci)
synth_mci_scalar, synth_mci_corrmat = dataset_mci_from_curves(synth_mci)

# ----------------------------------------------------------------------
# Clean NaNs before exporting / plotting
# ----------------------------------------------------------------------
real_mci_scalar  = np.asarray(real_mci_scalar, dtype=float)
synth_mci_scalar = np.asarray(synth_mci_scalar, dtype=float)

mask_real  = np.isfinite(real_mci_scalar)
mask_synth = np.isfinite(synth_mci_scalar)

if not np.all(mask_real):
    warn(f"MCI: dropped {np.sum(~mask_real)} real datasets with NaN MCI.")
if not np.all(mask_synth):
    warn(f"MCI: dropped {np.sum(~mask_synth)} synthetic datasets with NaN MCI.")

real_mci_scalar  = real_mci_scalar[mask_real].tolist()
synth_mci_scalar = synth_mci_scalar[mask_synth].tolist()

real_names_mci = [name for name, keep in zip(real_names_mci, mask_real) if keep]
synth_names    = [name for name, keep in zip(synth_names, mask_synth)   if keep]

# ======================================================================
# 5. SUMMARY (MCI + TCI) AND CSV EXPORT
# ======================================================================

section("MCI Null-Test Results")

mci_real_mean = float(np.nanmean(real_mci_scalar)) if real_mci_scalar else np.nan
mci_real_std  = float(np.nanstd(real_mci_scalar))  if real_mci_scalar else np.nan
mci_syn_mean  = float(np.nanmean(synth_mci_scalar)) if synth_mci_scalar else np.nan
mci_syn_std   = float(np.nanstd(synth_mci_scalar))  if synth_mci_scalar else np.nan

kv("Real MCI mean",      f"{mci_real_mean:.4f} ± {mci_real_std:.4f}")
kv("Synthetic MCI mean", f"{mci_syn_mean:.4f} ± {mci_syn_std:.4f}")
ok("MCI null-test complete\n")

# ----------------------------------------------------------------------

section("TCI Null-Test Results")

tci_real_mean = float(np.nanmean(real_pairwise_tci))
tci_real_std  = float(np.nanstd(real_pairwise_tci))
tci_syn_mean  = float(np.nanmean(synth_pairwise_tci))
tci_syn_std   = float(np.nanstd(synth_pairwise_tci))

kv("Real TCI mean",      f"{tci_real_mean:.4f} ± {tci_real_std:.4f}")
kv("Synthetic TCI mean", f"{tci_syn_mean:.4f} ± {tci_syn_std:.4f}")
ok("TCI null-test complete\n")

# ----------------------------------------------------------------------

section("Pairwise TCI Distributions")
kv("Real correlations",      len(real_pairwise_tci))
kv("Synthetic correlations", len(synth_pairwise_tci))

df_tci = pd.DataFrame({
    "Index":   "TCI",
    "Type":    (["Real"] * len(real_pairwise_tci)) +
               (["Synthetic"] * len(synth_pairwise_tci)),
    "Value":   np.concatenate([real_pairwise_tci, synth_pairwise_tci]),
    "Dataset": (["REAL_PAIR"] * len(real_pairwise_tci)) +
               (["SYNTH_PAIR"] * len(synth_pairwise_tci)),
})

df_mci = pd.DataFrame({
    "Index":   "MCI",
    "Type":    (["Real"] * len(real_mci_scalar)) +
               (["Synthetic"] * len(synth_mci_scalar)),
    "Value":   real_mci_scalar + synth_mci_scalar,
    "Dataset": real_names_mci + synth_names,
})

df_all = pd.concat([df_tci, df_mci], ignore_index=True)

csv_path = os.path.join(
    OUTDIR,
    f"UTMF_TCI_MCI_combined_nulltest_{NOW_TS}_UTC.csv"
)

df_all.to_csv(csv_path, index=False, float_format="%.6f")

# --------------------------------------------------------------------
# Append META row to CSV: includes fixed TCI time-series length
# --------------------------------------------------------------------
meta_row = pd.DataFrame([{
    "Index": "META",
    "Type": "TCI_length",
    "Value": FIXED_LEN,
    "Dataset": "GLOBAL"
}])

df_all_with_meta = pd.concat([df_all, meta_row], ignore_index=True)

csv_path = os.path.join(
    OUTDIR,
    f"UTMF_TCI_MCI_combined_nulltest_{NOW_TS}_UTC.csv"
)

df_all_with_meta.to_csv(csv_path, index=False, float_format="%.6f")

section("Saving Outputs")
ok(f"CSV written to: {csv_path}")


section("Saving Outputs")
ok(f"CSV written to: {csv_path}")

# ======================================================================
# 5b. PER-DATASET TCI BASED ON TCI h(q)-CURVES
# ======================================================================

real_tci_per_dataset, real_tci_corrmat = dataset_mci_from_curves(real_hq)
real_tci_per_dataset = np.asarray(real_tci_per_dataset, dtype=float)

# Map MCI values by dataset name
mci_by_name = {name: val for name, val in zip(real_names_mci, real_mci_scalar)}

# Intersect datasets for which both TCI and MCI exist
tci_vals = []
mci_vals = []
common_names = []

for name, tci_val in zip(real_names, real_tci_per_dataset):
    if name in mci_by_name:
        tci_vals.append(tci_val)
        mci_vals.append(mci_by_name[name])
        common_names.append(name)

tci_vals = np.asarray(tci_vals, dtype=float)
mci_vals = np.asarray(mci_vals, dtype=float)

info(f"TCI–MCI dataset intersection: {len(common_names)} datasets.\n")

# ======================================================================
# 6. MAIN TCI + MCI VIOLIN PLOT (FIGURE 1)
# ======================================================================

fig, axes = plt.subplots(1, 2, figsize=(15, 6))

# ----------- TCI violin ------------------------------------------------
sns.violinplot(
    ax=axes[0],
    data=df_tci,
    x="Type",
    y="Value",
    hue="Type",
    dodge=False,
    palette=["#1f77b4", "#d62728"],
    inner="quartile",
)

# Downsample synthetic points for stripplot clarity
df_tci_strip = df_tci.copy()
real_mask  = df_tci_strip["Type"] == "Real"
synth_mask = ~real_mask

df_tci_synth_sample = df_tci_strip[synth_mask].sample(
    n=1500, random_state=0
)

df_tci_strip = pd.concat(
    [df_tci_strip[real_mask], df_tci_synth_sample],
    ignore_index=True
)

sns.stripplot(
    ax=axes[0],
    data=df_tci_strip,
    x="Type",
    y="Value",
    dodge=False,
    jitter=True,
    size=3,
    color="black",
    alpha=0.7,
)

axes[0].set_title("TCI Null-Test", fontsize=12)
axes[0].set_ylabel("Temporal Coherence Index")
axes[0].legend([], [], frameon=False)

axes[0].text(
    0.98, 0.98,
    f"Real = {tci_real_mean:.3f} ± {tci_real_std:.3f}\n"
    f"Synth = {tci_syn_mean:.3f} ± {tci_syn_std:.3f}",
    transform=axes[0].transAxes,
    ha="right", va="top",
    fontsize=9,
    bbox=dict(facecolor="white", alpha=0.95,
              edgecolor="black", linewidth=0.5)
)

# ----------- MCI violin ------------------------------------------------
sns.violinplot(
    ax=axes[1],
    data=df_mci,
    x="Type",
    y="Value",
    hue="Type",
    dodge=False,
    palette=["#1f77b4", "#d62728"],
    inner="quartile",
)

sns.stripplot(
    ax=axes[1],
    data=df_mci,
    x="Type",
    y="Value",
    dodge=False,
    jitter=True,
    size=4,
    color="black",
    alpha=0.8,
)

axes[1].set_title("MCI Null-Test", fontsize=12)
axes[1].set_ylabel("Measurement Coherence Index")
axes[1].legend([], [], frameon=False)

axes[1].text(
    0.98, 0.98,
    f"Real = {mci_real_mean:.3f} ± {mci_real_std:.3f}\n"
    f"Synth = {mci_syn_mean:.3f} ± {mci_syn_std:.3f}",
    transform=axes[1].transAxes,
    ha="right", va="top",
    fontsize=9,
    bbox=dict(facecolor="white", alpha=0.95,
              edgecolor="black", linewidth=0.5),
)

fig.suptitle(
    f"UTMF v5.1 Null Tests — {PLOT_TIMESTAMP}",
    fontsize=13
)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
add_footer_timestamp(fig)
save_fig(fig, f"UTMF_TCI_MCI_nulltest_{NOW_TS}_UTC.png")

# ======================================================================
# 7. ADDITIONAL NULL-TEST FIGURES: HISTOGRAMS, CDFs, KDE, SCATTER
# ======================================================================

def empirical_cdf(x):
    x = np.sort(x)
    y = np.linspace(0, 1, len(x))
    return x, y

# -------- FIG A — Histogram: TCI ---------------------------------------
fig = plt.figure(figsize=(10, 6))

plt.hist(real_pairwise_tci,  bins=60, density=True, alpha=0.7,
         label="Real TCI",      color="black")
plt.hist(synth_pairwise_tci, bins=60, density=True, alpha=0.5,
         label="Synthetic TCI", color="steelblue")

plt.xlabel("Pairwise TCI  |corr(h(q)_i, h(q)_j)|")
plt.ylabel("Density")
plt.title("Distribution of Pairwise TCI Values — Real vs Synthetic")
plt.legend()
plt.grid(alpha=0.3)

add_footer_timestamp(fig)
save_fig(fig, f"fig_tci_hist_{NOW_TS}.png")

# -------- FIG B — Histogram: MCI ---------------------------------------
fig = plt.figure(figsize=(10, 6))

bins = np.linspace(0.0, 1.0, 40)

plt.hist(synth_mci_scalar, bins=bins, density=False, alpha=0.5,
         label="Synthetic MCI", color="tab:red")

plt.hist(real_mci_scalar,  bins=bins, density=False, alpha=0.9,
         label="Real MCI",      color="black")

plt.xlabel("Measurement Coherence Index")
plt.ylabel("Count")
plt.title("MCI Distribution — Real vs Synthetic")
plt.xlim(0.0, 1.0)
plt.legend()
plt.grid(alpha=0.3)

add_footer_timestamp(fig)
save_fig(fig, f"fig_mci_hist_{NOW_TS}.png")

# -------- FIG C — CDF: TCI --------------------------------------------
fig = plt.figure(figsize=(10, 6))

xr, yr = empirical_cdf(real_pairwise_tci)
xs, ys = empirical_cdf(synth_pairwise_tci)

plt.plot(xr, yr, label="Real TCI",      color="black",      linewidth=2)
plt.plot(xs, ys, label="Synthetic TCI", color="steelblue",  linewidth=2, alpha=0.8)

plt.xlabel("Pairwise TCI  |corr(h(q)_i, h(q)_j)|")
plt.ylabel("CDF")
plt.title("TCI CDF Comparison — Real vs Synthetic")
plt.legend()
plt.grid(alpha=0.3)

add_footer_timestamp(fig)
save_fig(fig, f"fig_tci_cdf_{NOW_TS}.png")

# -------- FIG D — CDF: MCI --------------------------------------------
fig = plt.figure(figsize=(10, 6))

xr, yr = empirical_cdf(real_mci_scalar)
xs, ys = empirical_cdf(synth_mci_scalar)

plt.plot(xr, yr, label="Real MCI",      color="black",    linewidth=2)
plt.plot(xs, ys, label="Synthetic MCI", color="darkred",  linewidth=2, alpha=0.8)

plt.xlabel("Measurement Coherence Index")
plt.ylabel("CDF")
plt.title("MCI CDF Comparison — Real vs Synthetic")
plt.legend()
plt.grid(alpha=0.3)

add_footer_timestamp(fig)
save_fig(fig, f"fig_mci_cdf_{NOW_TS}.png")

# -------- FIG E — Scatter: Per-Dataset TCI vs MCI ----------------------
fig = plt.figure(figsize=(10, 6))

colors = [color_for_domain(get_domain(name)) for name in common_names]

plt.scatter(
    tci_vals, mci_vals,
    c=colors,
    s=60,
    alpha=0.9,
    edgecolors="white",
    linewidths=0.5,
)

plt.xlabel("Per-dataset TCI (mean |corr(h(q)_i, h(q)_j)|)")
plt.ylabel("Per-dataset MCI (mean |corr(mean h(q)_i, mean h(q)_j)|)")
plt.title("TCI–MCI Relationship Across Real Datasets")
plt.grid(alpha=0.3)

# Domain legend
domain_handles = {}
for name in common_names:
    dom = get_domain(name)
    if dom not in domain_handles:
        domain_handles[dom] = plt.Line2D(
            [0], [0],
            marker="o", linestyle="",
            markersize=7,
            markerfacecolor=color_for_domain(dom),
            markeredgecolor="white",
            markeredgewidth=0.5,
            label=dom,
        )

plt.legend(
    handles=list(domain_handles.values()),
    title="Domain",
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
)

add_footer_timestamp(fig)
save_fig(fig, f"fig_tci_vs_mci_per_dataset_{NOW_TS}.png")

# -------- FIG F — KDE: TCI --------------------------------------------
fig = plt.figure(figsize=(10, 6))

sns.kdeplot(real_pairwise_tci,  label="Real",      color="black",     linewidth=2)
sns.kdeplot(synth_pairwise_tci, label="Synthetic", color="darkgreen", linewidth=2)

plt.xlabel("Pairwise |corr(h(q))|")
plt.ylabel("Density")
plt.title("TCI Kernel Density — Real vs Synthetic")
plt.grid(alpha=0.3)
plt.legend()

add_footer_timestamp(fig)
save_fig(fig, f"fig_tci_kernel_{NOW_TS}.png")

# ============================================================
# 8. REAL-DATA STRUCTURE PLOTS (REAL DATA ONLY)
# ============================================================

# R1 — Build H_real matrix from mean h(q) curves (real_mci)
if real_mci:
    H_real = np.vstack(real_mci)        # shape: (N_real, N_q)
    N_real = H_real.shape[0]
else:
    H_real = np.empty((0, len(q_values)))
    N_real = 0

# R2 — Euclidean distance distribution between real datasets
dists = []
if N_real >= 2:
    for i in range(N_real):
        for j in range(i + 1, N_real):
            d = np.linalg.norm(H_real[i] - H_real[j])
            dists.append(d)

    fig = plt.figure(figsize=(10, 6))
    sns.kdeplot(dists, color="black", linewidth=2)
    plt.hist(dists, bins=40, density=True, alpha=0.4, color="gray")
    plt.xlabel("Euclidean distance between mean h(q) curves")
    plt.ylabel("Density")
    plt.title("Real dataset distance structure in h(q)-space")
    plt.grid(alpha=0.3)
    add_footer_timestamp(fig)
    save_fig(fig, f"fig_real_distance_kde_{NOW_TS}.png")
else:
    warn(
        "Too few real MCI curves to construct a distance distribution — "
        "skipping fig_real_distance_kde."
    )

# R3 — Domain-level MCI scalar means (bar plot)
domain_mci = {}
for name, val in zip(real_names_mci, real_mci_scalar):
    dom = get_domain(name)
    domain_mci.setdefault(dom, []).append(val)

dom_labels = sorted(domain_mci.keys())
dom_means  = [np.mean(domain_mci[d]) for d in dom_labels]
dom_stds   = [np.std(domain_mci[d])  for d in dom_labels]
dom_colors = [color_for_domain(d) for d in dom_labels]

fig = plt.figure(figsize=(10, 6))
x = np.arange(len(dom_labels))
plt.bar(x, dom_means, yerr=dom_stds, capsize=5, color=dom_colors, alpha=0.9)
plt.xticks(x, dom_labels, rotation=45, ha="right")
plt.ylabel("Mean MCI (per domain)")
plt.title("Domain-level mean MCI across real datasets")
plt.grid(axis="y", alpha=0.3)
add_footer_timestamp(fig)
save_fig(fig, f"fig_real_domain_means_{NOW_TS}.png")

# ============================================================
# 9. ENSEMBLE MEAN h(q) — REAL VS SYNTHETIC (WITH q = 0 FIX)
# ============================================================

fig = plt.figure(figsize=(10, 6))

if real_hq:
    real_mean = np.nanmean(np.vstack(real_hq), axis=0)
else:
    real_mean = np.zeros_like(q_values)

if synth_hq:
    synth_mean = np.nanmean(np.vstack(synth_hq), axis=0)
else:
    synth_mean = np.zeros_like(q_values)

real_mean  = fix_hq_q0(real_mean,  q_values)
synth_mean = fix_hq_q0(synth_mean, q_values)

plt.plot(Q_PLOT, real_mean[Q_MASK],  label="Real mean h(q)",      color="black",  linewidth=3)
plt.plot(Q_PLOT, synth_mean[Q_MASK], label="Synthetic mean h(q)", color="purple", linewidth=3, alpha=0.7)
plt.xlabel("q")
plt.ylabel("h(q)")
plt.title("Mean h(q) — real vs synthetic ensembles")
plt.grid(alpha=0.3)
plt.legend()
add_footer_timestamp(fig)
save_fig(fig, f"fig_mean_hq_overlay_{NOW_TS}.png")

# ============================================================
# 10. ALL-DATA MEAN h(q) CURVES FOR DENDROGRAM / HEATMAP
# ============================================================

all_hq = {}
for name, res in results_all.items():
    if isinstance(res, tuple):
        res = res[0]
    if "hq_values" not in res:
        continue

    hq_list = res["hq_values"]
    if len(hq_list) == 0:
        continue

    hq_arr = np.asarray(hq_list, dtype=float)
    if hq_arr.ndim != 2:
        continue

    if np.all(np.isnan(hq_arr)):
        continue

    col_has_value = np.any(np.isfinite(hq_arr), axis=0)
    if not np.any(col_has_value):
        continue

    hq_reduced = hq_arr[:, col_has_value]
    if hq_reduced.size == 0:
        continue

    # Mean over q_valid
    mean_valid = np.nanmean(hq_reduced, axis=0)
    q_valid    = q_values[col_has_value]

    try:
        mean_full = np.interp(q_values, q_valid, mean_valid)
    except Exception:
        continue

    mean_full = fix_hq_q0(mean_full, q_values)
    if np.all(np.isnan(mean_full)):
        continue

    all_hq[name] = mean_full

names_all = list(all_hq.keys())
H_all = np.vstack([all_hq[n] for n in names_all])

# D1 — Dendrogram (Ward linkage on mean h(q))
D = pdist(H_all[:, Q_MASK], metric="euclidean")
Z = linkage(D, method="ward")

fig = plt.figure(figsize=(15, 7))
dendrogram(Z, labels=names_all, leaf_rotation=90)
plt.title("Hierarchical clustering of all datasets\n(Ward linkage on mean h(q))")
plt.ylabel("Cluster distance")
plt.tight_layout()
add_footer_timestamp(fig)
save_fig(fig, f"fig_real_dendrogram_{NOW_TS}.png")

from matplotlib.colors import PowerNorm

with np.errstate(invalid="ignore", divide="ignore"):
    # Correlation matrix over mean h(q), using only q ≠ 0
    corr_mat = np.corrcoef(H_all[:, Q_MASK])

# Define a distance: 0 = identical, larger = less coherent
dist_mat = 1.0 - np.abs(corr_mat)

# Color scale range: use percentile clipping to reduce outlier impact
vmin = np.nanmin(dist_mat)
vmax = np.nanpercentile(dist_mat, 95)

fig = plt.figure(figsize=(12, 10))

# PowerNorm enhances contrast at small distances
norm = PowerNorm(gamma=0.5, vmin=vmin, vmax=vmax)

im = plt.imshow(dist_mat, cmap="magma", norm=norm)
plt.colorbar(im, label="1 - |corr(mean h(q))|")

plt.xticks(range(len(names_all)), names_all, rotation=90)
plt.yticks(range(len(names_all)), names_all)
plt.title("Pairwise distance 1 - |corr(mean h(q))| across all datasets")

plt.tight_layout()
add_footer_timestamp(fig)
save_fig(fig, f"fig_real_tci_heatmap_{NOW_TS}.png")

# ============================================================
# 10b. PCA EMBEDDING OF MEAN h(q) CURVES (REAL DATASETS)
# ============================================================

try:
    from sklearn.decomposition import PCA
except ImportError:
    warn("scikit-learn is not installed — PCA plot will be skipped.")
else:
    if H_all.shape[0] < 2:
        warn("Too few datasets for PCA — fig_real_pca_hq will be skipped.")
    else:
        # PCA on q ≠ 0 only (consistent with the heatmap)
        H_pca = H_all[:, Q_MASK]

        pca = PCA(n_components=2)
        X2 = pca.fit_transform(H_pca)

        fig = plt.figure(figsize=(10, 7))

        # Color points by domain
        point_colors = [color_for_domain(get_domain(nm)) for nm in names_all]

        plt.scatter(
            X2[:, 0], X2[:, 1],
            c=point_colors,
            s=70,
            alpha=0.9,
            edgecolors="white",
            linewidths=0.6,
        )

        plt.xlabel(f"PC1  ({pca.explained_variance_ratio_[0]*100:.1f}% var)")
        plt.ylabel(f"PC2  ({pca.explained_variance_ratio_[1]*100:.1f}% var)")
        plt.title("PCA of mean h(q) curves (q ≠ 0)")

        # Compact domain legend
        domain_handles = {}
        for nm in names_all:
            dom = get_domain(nm)
            if dom not in domain_handles:
                domain_handles[dom] = plt.Line2D(
                    [0], [0],
                    marker="o", linestyle="",
                    markersize=7,
                    markerfacecolor=color_for_domain(dom),
                    markeredgecolor="white",
                    markeredgewidth=0.5,
                    label=dom,
                )

        plt.legend(
            handles=list(domain_handles.values()),
            title="Domain",
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0,
        )

        plt.grid(alpha=0.3)
        add_footer_timestamp(fig)
        plt.tight_layout()
        save_fig(fig, f"fig_real_pca_hq_{NOW_TS}.png")

def plot_hq_overlay_raw_from_json(utmf_json):
    """
    Raw mean h(q) overlays per dataset, colored by domain.

    This variant:
    - uses the subset-level h(q) stored in `subsets_processed`,
    - averages h(q) per dataset,
    - applies the q = 0 repair,
    - plots all datasets with domain-colored curves,
    - and builds a domain-level legend.
    """

    curves = extract_hq_curves_from_json(utmf_json)
    if not curves:
        warn(
            "No h(q) curves found in JSON (subsets_processed is empty or unusable). "
            "Raw overlay will not be generated."
        )
        return None, None

    fig, ax = plt.subplots(figsize=(10, 6))

    domain_handles = {}
    all_h_values = []

    for ds_name, ch in curves.items():
        q = ch["q"]
        h = ch["h"]

        # Skip datasets with no finite values
        if not np.any(np.isfinite(h)):
            continue

        dom   = get_domain(ds_name)
        color = color_for_domain(dom)

        # One representative handle per domain for the legend
        if dom not in domain_handles:
            (handle,) = ax.plot(
                q, h,
                color=color,
                alpha=0.30,
                linewidth=1.3,
                label=dom
            )
            domain_handles[dom] = handle
        else:
            ax.plot(
                q, h,
                color=color,
                alpha=0.30,
                linewidth=1.0
            )

        all_h_values.append(h)

    if not all_h_values:
        warn("All h(q) curves were fully NaN — nothing to plot in raw overlay.")
        return fig, ax

    all_h_values = np.concatenate(all_h_values)

    # Use only finite values for y-axis limits
    finite_mask = np.isfinite(all_h_values)
    if np.any(finite_mask):
        vmin = all_h_values[finite_mask].min()
        vmax = all_h_values[finite_mask].max()
        yrange = vmax - vmin
        margin = 0.05 * yrange if yrange > 0 else 0.1
        ax.set_ylim(vmin - margin, vmax + margin)
    else:
        warn(
            "All h(q) values are NaN/Inf when determining y-limits — "
            "falling back to Matplotlib defaults."
        )
        # let Matplotlib choose default y-limits

    ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.5, linestyle="--")

    ax.set_xlabel("q")
    ax.set_ylabel("h(q)")
    ax.set_title("Mean h(q) per dataset (raw overlay)")
    ax.grid(alpha=0.3)

    # Domain-based legend from the collected handles
    handles = list(domain_handles.values())
    if handles:
        ax.legend(handles=handles, loc="best", title="Domain")

    fig.tight_layout()
    return fig, ax

fig, ax = plot_hq_overlay_raw_from_json(utmf_json)
if fig is not None:
    add_footer_timestamp(fig)
    save_fig(fig, f"fig_real_hq_overlay_raw_{NOW_TS}.png")
else:
    print("No figure generated for raw h(q) overlay.")

print("\n[OK] All supplementary figures generated.\n")
banner("UTMF v5.1 TCI/MCI NULL-TEST AND STRUCTURE PLOTS COMPLETED SUCCESSFULLY")
