# ================================================================
# UTMF v5.1 – Unified Temporal–Measurement Framework
# Configuration & Data Loading Module
#
# Author: Jedi Markus Strive
# Framework Version: 5.1
# Last Updated: December 11, 2025
#
# This module provides:
#   • Dataset loading and preprocessing for all supported domains
#   • Initial multifractal analysis via MFDFA (jedi_mfdfa, order = 0)
#   • Optional computation of TCI / MCI / TMCI indices
#   • Centralised metadata logging for reproducibility
#
# UTMF v5.1 introduces:
#   • A unified, numerically stable TCI/MCI/TMCI pipeline
#   • Fully reproducible JSON snapshots (FULL_DETAILS)
#   • Embedded MFDFA implementation with SHA-256 integrity check
#   • Improved, memory-safe NANOGrav processing with detailed diagnostics
#
# Workflow (Phase 1):
#   1. Mount storage and initialise paths
#   2. Load configuration; enable/disable datasets via 'utmf_use'
#   3. Load → preprocess/denoise → MFDFA → (optional) compute indices
#   4. Persist outputs:
#        - Lightweight CSV with subset summaries and global indices
#        - Full JSON artifact (if enabled), containing:
#              * configuration snapshot
#              * TCI time-series snapshots
#              * MFDFA subset results for all datasets
#              * TCI/MCI/TMCI + bootstrap + cross-validation
#              * embedded jedi_mfdfa source
#
# This file is executed once per UTMF run; all downstream analysis
# (Phase 2 notebooks) can be performed using ONLY the JSON snapshot.
# ================================================================

# Required packages (install once per session)
!pip install healpy joblib h5py uproot pywavelets pint-pulsar numpy pandas scipy matplotlib tqdm dtaidistance tensorly dask tensorflow torch plotly gudhi ripser scikit-learn seaborn pytz networkx

# Import libraries
import networkx as nx
import numpy as np
import pandas as pd
import h5py
import copy
from astropy.io import fits
import healpy as hp
import scipy.signal
import dask.array as da
import gc
from tqdm import tqdm
from numba import jit
import uproot
import inspect, hashlib
import pywt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold  # Voor cross-val
from scipy.stats import pearsonr, norm, entropy
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
from google.colab import drive
from joblib import Parallel, delayed
import os
import glob
import warnings
import json  # For serializing nested metadata to CSV
from datetime import datetime
import pytz
import requests
import time
from bs4 import BeautifulSoup
import warnings
warnings.filterwarnings("ignore")
try:
    from pint.models import get_model, get_model_and_toas
    from pint.toa import get_TOAs
except ImportError:
    print("Error: pint-pulsar not installed. Install with `pip install pint-pulsar`.")

# Mount Google Drive
drive.mount('/content/drive', force_remount=True)

# Output directory for CSVs
OUTPUT_DIR = '/content/drive/MyDrive/Datasets_UTMF/UTMF_outputs/'
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Output dir: {OUTPUT_DIR}")

# Element names (for NIST_3 dataset)
ELEMENT_NAMES = {
    'Ac': 'Actinium', 'Ag': 'Silver', 'Al': 'Aluminium', 'Ar': 'Argon', 'B': 'Boron',
    'Ba': 'Barium', 'Be': 'Beryllium', 'Bi': 'Bismuth', 'Br': 'Bromine', 'C': 'Carbon',
    'Ca': 'Calcium', 'Cd': 'Cadmium', 'Ce': 'Cerium', 'Cl': 'Chlorine', 'Co': 'Cobalt',
    'Cr': 'Chromium', 'Cs': 'Cesium', 'Cu': 'Copper', 'Dy': 'Dysprosium', 'Eu': 'Europium',
    'F': 'Fluorine', 'Fe': 'Iron', 'H': 'Hydrogen', 'He': 'Helium', 'Hf': 'Hafnium',
    'Hg': 'Mercury', 'I': 'Iodine', 'In': 'Indium', 'Ir': 'Iridium', 'K': 'Potassium',
    'Kr': 'Krypton', 'La': 'Lanthanum', 'Li': 'Lithium', 'Mg': 'Magnesium', 'Mn': 'Manganese',
    'Mo': 'Molybdenum', 'N': 'Nitrogen', 'Na': 'Sodium', 'Nb': 'Niobium', 'Nd': 'Neodymium',
    'Ne': 'Neon', 'Ni': 'Nickel', 'O': 'Oxygen', 'Os': 'Osmium', 'P': 'Phosphorus',
    'Pr': 'Praseodymium', 'Pt': 'Platinum', 'Rb': 'Rubidium', 'Re': 'Rhenium', 'Rh': 'Rhodium',
    'S': 'Sulfur', 'Sc': 'Scandium', 'Si': 'Silicon', 'Sn': 'Tin', 'Sr': 'Strontium',
    'Ta': 'Tantalum', 'Ti': 'Titanium', 'Tm': 'Thulium', 'V': 'Vanadium', 'W': 'Tungsten',
    'Xe': 'Xenon', 'Y': 'Yttrium', 'Zr': 'Zirconium'
}

# Global metadata logger (for save_run_metadata) - Enhanced for full subset logging
metadata_log = {
    "run_timestamp": None,
    "config_snapshot": None,
    "datasets_loaded": {},      # existing
    "subsets_processed": {},    # existing
    "cca_results": None,
    "tci_pairs": [],
    "mci_measurements": [],
    "tmci_folds": [],
    "errors": [],
    "tci": None,
    "mci": None,
    "tmci": None,
    "tmci_ci": None,
    "tmci_corr": None,
    "tmci_std_cv": None,
    "tci_meta": {},             # TCI dataset diagnostics
    "mci_meta": {},             # MCI pair metadata
    "tmci_bootstrap": {},       # tci_samples, mci_samples, tmci_samples
    "load_meta": {},            # dataset loading diagnostics
    "subset_warnings": {}       # per-subset warnings/errors
}

MAX_SUBSETS_IN_CSV = 50000        # How many rows in CSV.
SAVE_FULL_DETAILS_JSON = True    # Turn to True to save the JSON

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S UTC")
print(f"Run started: {timestamp}")

# Configuration for datasets
CONFIG = {
    'ligo_files': [
        '/content/drive/MyDrive/Datasets_UTMF/L-L1_GWOSC_O4a_16KHZ_R1-1368350720-4096.hdf5',
        '/content/drive/MyDrive/Datasets_UTMF/L-L1_GWOSC_04a_16KHZ_R1-1384779776-4096.hdf5',
        '/content/drive/MyDrive/Datasets_UTMF/L-L1_GWOSC_O4a_16KHZ_R1-1370202112-4096.hdf5',
        '/content/drive/MyDrive/Datasets_UTMF/L-L1_GWOSC_O4a_16KHZ_R1-1389420544-4096.hdf5'
    ],
    'ligo_names': ['LIGO-L1_1368350720', 'LIGO-L1_1384779776(GW231123)', 'LIGO-L1_1370202112', 'LIGO-L1_1389420544'],
    'ligo': [
        {'sample_rate': 16384,     # LIGO-L1_1368350720
         'total_duration': 4096,
         'subset_duration': 4,
         'n_subsets': 100,
         'freq_range': [1, 30],
         'expected_D_f': 1.22,
         'sigma_D_f': 0.05,
         'min_std': 1e-5,
         'scales': np.logspace(np.log10(8), np.log10(4 * 16384 / 16), 20, dtype=np.int32),
         'utmf_use': True
         },
        {'sample_rate': 16384,     # LIGO-L1_1384779776(GW231123)
         'total_duration': 4096,
         'subset_duration': 4,
         'n_subsets': 100,
         'freq_range': [1, 30],
         'expected_D_f': 1.22,
         'sigma_D_f': 0.05,
         'min_std': 1e-5,
         'scales': np.logspace(np.log10(8), np.log10(4 * 16384 / 16), 20, dtype=np.int32),
         'utmf_use': False
         },
        {'sample_rate': 16384,     # LIGO-L1_1370202112
         'total_duration': 4096,
         'subset_duration': 4,
         'n_subsets': 100,
         'freq_range': [1, 30],
         'expected_D_f': 1.22,
         'sigma_D_f': 0.05,
         'min_std': 1e-5,
         'scales': np.logspace(np.log10(8), np.log10(4 * 16384 / 16), 20, dtype=np.int32),
         'utmf_use':False
         },
        {'sample_rate': 16384,    # LIGO-L1_1389420544
         'total_duration': 4096,
         'subset_duration': 4,
         'n_subsets': 100,
         'freq_range': [1, 30],
         'expected_D_f': 1.22,
         'sigma_D_f': 0.05,
         'min_std': 1e-5,
         'scales': np.logspace(np.log10(8), np.log10(4 * 16384 / 16), 20, dtype=np.int32),
         'utmf_use': False
         }
    ],
    'cmb_files': [
        '/content/drive/MyDrive/Datasets_UTMF/COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits',
        '/content/drive/MyDrive/Datasets_UTMF/COM_CMB_IQU-smica_2048_R3.00_full.fits',
        '/content/drive/MyDrive/Datasets_UTMF/LFI_SkyMap_070_1024_R3.00_survey-1.fits'
    ],
    'cmb_names': ['Planck_CMB_I-Stokes_nosz', 'Planck_CMB_I-Stokes', 'Planck_LFI_70GHz'],
    'cmb': [
        {'nside': 2048,          # Planck_CMB_I-Stokes_nosz
         'subset_size': 100000,
         'n_subsets': 100,
         'expected_D_f': 1.19,
         'sigma_D_f': 0.04,
         'galactic_mask': False,
         'disc_radius': np.radians(20),
         'min_std': 1e-5,
         'fields': [0],
         'scales': np.array([2, 4, 8, 16, 32, 64, 128], dtype=np.int32),
         'utmf_use': True
         },
        {'nside': 2048,          # Planck_CMB_I-Stokes
         'subset_size': 100000,
         'n_subsets': 100,
         'expected_D_f': 1.19,
         'sigma_D_f': 0.04,
         'galactic_mask': True,
         'disc_radius': np.radians(20),
         'min_std': 1e-5,
         'fields': [0],
         'scales': np.array([2, 4, 8, 16, 32, 64, 128], dtype=np.int32),
         'utmf_use': False
         },
        {'nside': 1024,           # Planck_LFI_70GHz
         'subset_size': 45000,
         'n_subsets': 250,
         'expected_D_f': 1.19,
         'sigma_D_f': 0.04,
         'galactic_mask': True,
         'disc_radius': np.radians(30),
         'min_std': 1e-5,
         'fields': [0],
         'scales': np.array([2, 4, 8, 16, 32, 64, 128], dtype=np.int32),
         'utmf_use': False
         }
    ],
    'desi_file': '/content/drive/MyDrive/Datasets_UTMF/LRG_full.dat.fits',
    'desi_name': 'DESI_LRG',
    'desi': {'subset_size': 3700,
             'n_subsets': 100,
             'expected_D_f': 1.19,
             'sigma_D_f': 0.04,
             'min_std': 1e-5,
             'columns': ['FLUX_Z', 'FLUX_G', 'FLUX_R'],
             'scales': np.array([2, 4, 8, 16, 32, 64, 128], dtype=np.int32),
             'utmf_use': True
             },
    'cern_file': '/content/drive/MyDrive/Datasets_UTMF/data_B.exactly2lep.root',
    'cern_name': 'CERN_2Lepton',
    'cern': {'subset_size': 5000,
             'n_subsets': 100,
             'expected_D_f': 1.19,
             'sigma_D_f': 0.04,
             'min_std': 1e-5, 'tree': 'mini', 'columns': ['lep_pt', 'lep_eta', 'lep_phi'],
             'scales': np.array([2, 4, 8, 16, 32, 64, 128], dtype=np.int32),
             'utmf_use': True
             },
    'nist_name': 'NIST_elements',  # Contains 63 elements
    'nist': {
        'elements': [],  # Initialize empty for loop
        'subset_size': lambda n: 50 if n >= 100 else max(10, n // 5),                                 # ← dynamic
        'scales': lambda n: np.logspace(np.log10(4), np.log10(min(n // 4, 500)), 15, dtype=np.int32), # ← dynamic
        'n_subsets': 100,
        'expected_D_f': 1.16,
        'sigma_D_f': 0.04,
        'min_std': 1e-5,
        'columns': ['intens', 'wn(cm-1)', 'Aki(10^8 s^-1)', 'obs_wl_vac(nm)'],
        'elements_list': [  # All elements in NIST_elements file (63)
            ['Ac'], ['Ag'], ['Al'], ['Ar'], ['B'], ['Ba'], ['Be'], ['Bi'], ['Br'], ['C'],
            ['Ca'], ['Cd'], ['Ce'], ['Cl'], ['Co'], ['Cr'], ['Cs'], ['Cu'], ['Dy'], ['Eu'],
            ['F'], ['Fe'], ['H'], ['He'], ['Hf'], ['Hg'], ['I'], ['In'], ['Ir'], ['K'],
            ['Kr'], ['La'], ['Li'], ['Mg'], ['Mn'], ['Mo'], ['N'], ['Na'], ['Nb'], ['Nd'],
            ['Ne'], ['Ni'], ['O'], ['Os'], ['P'], ['Pr'], ['Pt'], ['Rb'], ['Re'], ['Rh'],
            ['S'], ['Sc'], ['Si'], ['Sn'], ['Sr'], ['Ta'], ['Ti'], ['Tm'], ['V'], ['W'],
            ['Xe'], ['Y'], ['Zr']
        ],
        'elements_list_utmf': [  # Elements for UTMF analysis
           
            ['Ca'], ['Cd'], ['Ce']           
        ],
        'utmf_use': True  # For elements in elements_list_utmf
    },
    'nanograv': {  # Contains 83 pulsars
        'base_dir': '/content/drive/MyDrive/Datasets_UTMF/NANOGrav15yr_PulsarTiming_v2.1.0',
        'file_templates': {
            'tim': 'narrowband/tim/{pulsar_name}_PINT_{date}.nb.tim',
            'par': 'narrowband/par/{pulsar_name}_PINT_{date}.nb.par',
            'res_full': 'residuals/{pulsar_name}_NG15yr_nb.full.res',
            'res_avg': 'residuals/{pulsar_name}_NG15yr_nb.avg.res',
            'dmx': 'narrowband/dmx/{pulsar_name}_dmxparse.nb.out',
            'noise': 'narrowband/noise/{pulsar_name}.nb.pars.txt',
            'template': 'narrowband/template/{pulsar_name}.sum.sm',
            'clock': 'clock/time_vla.dat'
        },
        'default_date': '20220302',
        'date_exceptions': {'J1713+0747': '20220309'},
        'pulsar_list': [  # Pulsars in NANOGrav file
            'J0030+0451', 'J0340+4130', 'J0406+3039', 'J0437-4715', 'J0509+0856',
            'J0557+1551', 'J0605+3757', 'J0610-2100', 'J0613-0200', 'J0614-3329',
            'J0636+5128', 'J0645+5158', 'J0709+0458', 'J0740+6620', 'J0931-1902',
            'J1012+5307', 'J1012-4235', 'J1022+1001', 'J1024-0719', 'J1125+7819',
            'J1312+0051', 'J1453+1902', 'J1455-3330', 'J1600-3053', 'J1614-2230',
            'J1630+3734', 'J1640+2224', 'J1643-1224', 'J1705-1903', 'J1713+0747',
            'J1719-1438', 'J1730-2304', 'J1738+0333', 'J1741+1351', 'J1744-1134',
            'J1745+1017', 'J1747-4036', 'J1802-2124', 'J1811-2405', 'J1832-0836',
            'J1843-1113', 'J1853+1303', 'J1903+0327', 'J1909-3744', 'J1910+1256',
            'J1911+1347', 'J1918-0642', 'J1923+2515', 'J1944+0907', 'J1946+3417',
            'J2010-1323', 'J2017+0603', 'J2022+2534', 'J2033+1734', 'J2043+1711',
            'J2124-3358', 'J2214+3000', 'J2229+2643', 'J2234+0611', 'J2234+0944',
            'J2302+4442', 'J2317+1439', 'J2322+2057', 'B1855+09', 'B1937+21',
            'B1953+29', 'J0751+1807', 'J0023+0923', 'J1751-2857', 'J0125-2327',
            'J0732+2314', 'J1221-0633', 'J1400-1431', 'J1630+3550', 'J2039-3616',
            'J0218+4232', 'J0337+1715', 'J0621+2514', 'J0721-2038', 'J1803+1358',
            'B1257+12', 'J1327+3423', 'J0154+1833'
        ],
        'pulsar_list_utmf': [  # Pulsars for UTMF analys
            'J1843-1113', 'J1853+1303', 'J1903+0327'
        ],
         'subset_size': 250,
         'n_subsets': 100,
         'expected_D_f': 1.19,
         'sigma_D_f': 0.04,
         'min_std': 1e-5,
         'columns': ['residuals', 'dmx', 'red_noise'],
         'scales': np.array([2, 4, 8, 16, 32, 64, 128, 256], dtype=np.int32),
         'utmf_use': True  # For pulsars in pulsar_list_utmf
    },
    'qrng': {'subset_size': 2560,
             'n_subsets': 100,
             'expected_D_f': 1.21,
             'sigma_D_f': 0.04,
             'min_std': 1e-5,
             'scales': np.array([2, 4, 8, 16, 32, 64, 128, 256, 512, 1024], dtype=np.int32),
             'utmf_use': True
             },
    'gaia': {'file': '/content/drive/MyDrive/Datasets_UTMF/gaia_dr3.tsv',
             'name': 'Gaia_DR3',
             'subset_size': 7500,
             'n_subsets': 100,
             'expected_D_f': 1.19,
             'sigma_D_f': 0.04,
             'min_std': 1e-5,
             'columns': ['RA_ICRS', 'DE_ICRS', 'pmRA', 'pmDE', 'Gmag'],
             'sep': '\t',
             'scales': np.array([1, 2, 4, 6, 8, 16, 32, 64], dtype=np.int32),
             'utmf_use': True
             },
    'mfdfa': {
              'q_values': np.arange(-8, 8.2, 0.2),
              'detrend_order': 0
              },
    # UTMF thresholds
    'utmf': {
             'd_f_min_threshold': 1.0,
             'd_f_max_threshold': 2.0
             },
    # Nieuw: Flags voor computations
    'compute_indices': True,  # Turn TCI/MCI/TMCI on/off
    'cross_val': True,       # Simple k=3 fold for TMCI std (voor p<0.05)
    # Metadata
    'metadata': {'save_flag': True, 'log_level': 'info'}
}

# Numba-compatible linear detrending
@jit(nopython=True)
def polyfit_linear(x, y, lambda_reg=1e-5):
    # Linear polyfit with regularization for stability
    n = len(x)
    sum_x = np.sum(x)
    sum_y = np.sum(y)
    sum_xy = np.sum(x * y)
    sum_x2 = np.sum(x * x)
    denom = n * sum_x2 - sum_x**2 + lambda_reg
    if abs(denom) < 1e-10:
        return np.array([0.0, sum_y / n])
    m = (n * sum_xy - sum_x * sum_y) / denom
    b = (sum_y * sum_x2 - sum_x * sum_xy) / denom
    return np.array([m, b])

@jit(nopython=True)
def polyval_linear(coeffs, x):
    return coeffs[0] * x + coeffs[1]

# Unified MFDFA function
@jit(nopython=True)
def jedi_mfdfa(data, scales, q_values, detrend_order=0):
    n = len(data)
    fluct = np.zeros((len(q_values), len(scales)))
    rms_values = []
    slopes = np.zeros(len(q_values))

    for i in range(len(scales)):
        s = scales[i]
        segments = n // s
        if segments < 2:
            fluct[:, i] = np.nan
            continue
        rms = np.zeros(segments)
        valid_segments = 0
        for v in range(segments):
            segment = data[v*s:(v+1)*s]
            if len(segment) != s or np.std(segment) < 1e-10:
                continue
            x = np.arange(s, dtype=np.float64)
            if detrend_order > 0:
                try:
                    coeffs = polyfit_linear(x, segment)
                    trend = polyval_linear(coeffs, x)
                    detrended = segment - trend
                except:
                    detrended = segment - np.sum(segment) / s
            else:
                detrended = segment - np.sum(segment) / s
            sum_squares = 0.0
            for j in range(s):
                sum_squares += detrended[j]**2
            rms_val = np.sqrt(sum_squares / s + 1e-12)
            if rms_val > 1e-10:
                rms[valid_segments] = rms_val
                valid_segments += 1
        if valid_segments < 2:
            fluct[:, i] = np.nan
            continue
        rms = rms[:valid_segments]
        rms_values.append(rms)
        for j in range(len(q_values)):
            q = q_values[j]
            if q == 0:
                sum_log = 0.0
                count = 0
                for k in range(valid_segments):
                    if rms[k] > 1e-10:
                        sum_log += np.log(rms[k]**2 + 1e-12)
                        count += 1
                fluct[j, i] = np.exp(0.5 * (sum_log / count)) if count > 0 else np.nan
            else:
                sum_power = 0.0
                count = 0
                for k in range(valid_segments):
                    if rms[k] > 1e-10:
                        sum_power += (rms[k] + 1e-12)**q
                        count += 1
                fluct[j, i] = (sum_power / count)**(1/q) if count > 0 else np.nan
                if not np.isfinite(fluct[j, i]) or fluct[j, i] <= 0:
                    fluct[j, i] = np.nan
    valid_scales = np.sum(np.isfinite(fluct), axis=0)
    if np.max(valid_scales) < 4:
        return np.nan, np.full(len(q_values), np.nan), rms_values, fluct, slopes
    for j in range(len(q_values)):
        valid = np.isfinite(fluct[j, :]) & (fluct[j, :] > 0)
        if np.sum(valid) < 4:
            slopes[j] = np.nan
            continue
        coeffs = np.zeros(2)
        X = np.log(scales[valid])
        Y = np.log(fluct[j, valid] + 1e-12)
        n_valid = len(X)
        sum_x = np.sum(X)
        sum_y = np.sum(Y)
        sum_xy = np.sum(X * Y)
        sum_x2 = np.sum(X * X)
        denom = n_valid * sum_x2 - sum_x**2 + 1e-5
        if abs(denom) > 1e-10:
            coeffs[0] = (n_valid * sum_xy - sum_x * sum_y) / denom
            coeffs[1] = (sum_y * sum_x2 - sum_x * sum_xy) / denom
        slopes[j] = coeffs[0] if coeffs[0] > 0 else np.nan
    hq = slopes
    valid_hq = np.isfinite(hq)
    if np.sum(valid_hq) >= 2:
        tau = hq * q_values - 1
        alpha = np.diff(tau[valid_hq]) / np.diff(q_values[valid_hq])
        f_alpha = q_values[valid_hq][1:] * alpha - tau[valid_hq][1:]
        D_f = np.nanmean(alpha) if np.isfinite(alpha).any() else np.nan
    else:
        D_f = np.nan
    return D_f, hq, rms_values, fluct, slopes

# Denoising function
def denoise_data(data, data_type, ligo_idx=None, cmb_idx=None):
    try:
        std_before = np.std(data)
       # print(f"Standard deviation before denoising ({data_type}): {std_before:.2e}") # Remove first # to see denoising
        if data_type == 'ligo':
            denoised = data * 1e16
            denoised = scipy.signal.savgol_filter(denoised, window_length=7, polyorder=1, mode='nearest')
        elif data_type == 'nist':
            data = np.log1p(np.abs(data))
            coeffs = pywt.wavedec(data, 'db4', level=4 if data_type == 'nist' else 5)
            threshold = np.std(coeffs[-1]) * np.sqrt(2 * np.log(len(data)))
            coeffs = [pywt.threshold(c, threshold, mode='soft') for c in coeffs]
            denoised = pywt.waverec(coeffs, 'db4')
            denoised = denoised[:len(data)]
            denoised = scipy.signal.savgol_filter(denoised, window_length=5, polyorder=1, mode='nearest')
        elif data_type == 'nanograv':
            coeffs = pywt.wavedec(data, 'db4', level=5)
            threshold = np.std(coeffs[-1]) * np.sqrt(2 * np.log(len(data)))
            coeffs = [pywt.threshold(c, threshold, mode='soft') for c in coeffs]
            denoised = pywt.waverec(coeffs, 'db4')
            denoised = denoised[:len(data)]
            denoised = scipy.signal.savgol_filter(denoised, window_length=11, polyorder=2, mode='nearest')
        elif data_type == 'qrng':
            # Normalize randomness (no trend, but std=1 for coherence)
            denoised = (data - np.mean(data)) / (np.std(data) + 1e-10)
        else:  # cmb, desi, cern, gaia
            denoised = data / np.std(data)
        std_after = np.std(denoised)
        min_std = CONFIG[data_type]['min_std'] if data_type not in ['ligo', 'cmb'] else \
                  CONFIG['ligo'][ligo_idx]['min_std'] if data_type == 'ligo' else \
                  CONFIG['cmb'][cmb_idx]['min_std']
       # print(f"Standard deviation after denoising ({data_type}): {std_after:.2e}") # Remove first # to see denoising
        if std_after < min_std:
            print(f"Warning: Low variability after denoising ({data_type}, std={std_after:.2e})")
        return denoised.astype(np.float64)
    except Exception as e:
        print(f"Error in denoising {data_type}: {e}")
        return data.astype(np.float64)

        # Extraction tool for TCI (time-series) and MCI (measurement vectors)
def extract_tci_mci_data(raw_data, data_type, ligo_idx=None, cmb_idx=None):
    """
    Extracts dataset-specific inputs for:
        - TCI: a single denoised time series (typically the longest valid sequence)
        - MCI: one or more measurement channels (DMX, red-noise, PSD, or other columns)

    This function supports ALL UTMF dataset types, including NANOGrav.

    Parameters
    ----------
    raw_data : array-like, dict, or structured object
        Raw input data for the dataset.
        For NANOGrav this is a dictionary containing 'residuals', 'dmx', etc.
    data_type : str
        One of: 'nanograv', 'ligo', 'cmb', 'desi', 'cern', 'nist', 'qrng', 'gaia', ...
    ligo_idx, cmb_idx : int or None
        Optional dataset indices for selecting LIGO/CMB parameters.

    Returns
    -------
    tci_extracted : ndarray or None
        A single denoised time series (fixed-length if needed), or None if unavailable.
    mci_extracted : ndarray or None
        A measurement matrix aligned with the TCI length, or None if unavailable.
    """
    MAX_LEN = 1_000_000  # Upper bound for downsampling

    # ==================================================================
    # 1. NANOGrav special handling (raw_data is a dictionary)
    # ==================================================================
    if data_type == 'nanograv':

        # Residuals are mandatory for TCI
        if not isinstance(raw_data, dict) or 'residuals' not in raw_data:
            print("   → NANOGrav: no residuals found → TCI/MCI = None")
            return None, None

        # --- TCI extraction ---
        tci_raw = np.asarray(raw_data['residuals'], dtype=np.float64)

        if len(tci_raw) < 1024 or not np.isfinite(tci_raw).any():
            print("   → NANOGrav: residuals too short or invalid → skipping TCI")
            tci_extracted = None
        else:
            tci_extracted = denoise_data(tci_raw, data_type)
            print(f"   → TCI uses residuals (length = {len(tci_extracted):,})")

        # --- MCI extraction ---
        # Prefer DMX + red-noise channels when available
        mci_parts = []
        if 'dmx' in raw_data and len(raw_data['dmx']) > 10:
            mci_parts.append(denoise_data(np.asarray(raw_data['dmx'], dtype=np.float64), data_type))
        if 'red_noise' in raw_data and len(raw_data['red_noise']) > 10:
            mci_parts.append(denoise_data(np.asarray(raw_data['red_noise'], dtype=np.float64), data_type))

        if mci_parts:
            # Align all channels to the TCI length (or MAX_LEN if TCI missing)
            target_len = len(tci_extracted) if tci_extracted is not None else MAX_LEN
            padded = []
            for m in mci_parts:
                if len(m) < target_len:
                    m = np.pad(m, (0, target_len - len(m)), mode='constant')
                else:
                    m = m[:target_len]
                padded.append(m)
            mci_extracted = (
                np.column_stack(padded) if len(padded) > 1
                else padded[0].reshape(-1, 1)
            )
        else:
            # Fallback: PSD of the residuals
            data_for_psd = tci_extracted if tci_extracted is not None else tci_raw
            f, psd = scipy.signal.welch(
                data_for_psd,
                fs=1.0,
                nperseg=min(1024, len(data_for_psd) // 4)
            )
            psd_norm = psd / (np.max(psd) + 1e-10)

            # Tile PSD to match the TCI length
            repeat = (len(data_for_psd) // len(psd_norm)) + 1
            mci_extracted = np.tile(psd_norm, (repeat, 1))[:len(data_for_psd), :8]

        # Downsample if overly long
        if tci_extracted is not None and len(tci_extracted) > MAX_LEN:
            idx = np.linspace(0, len(tci_extracted) - 1, MAX_LEN, dtype=int)
            tci_extracted = tci_extracted[idx]

        gc.collect()
        return tci_extracted, mci_extracted

    # ==================================================================
    # 2. General case for all other dataset types
    # ==================================================================

    # Standardise input to list of 1D signals
    if isinstance(raw_data, (list, tuple)):
        signals = [np.asarray(s, dtype=np.float64) for s in raw_data]
    elif hasattr(raw_data, 'ndim') and raw_data.ndim > 1:
        signals = [raw_data[:, i] for i in range(raw_data.shape[1])]
    else:
        signals = [np.asarray(raw_data, dtype=np.float64)]

    # Filter out invalid or too-short signals
    valid_signals = [s for s in signals if len(s) >= 1024 and np.isfinite(s).any()]
    if not valid_signals:
        return None, None

    # --- Select TCI as the longest valid signal ---
    tci_candidates = [(len(s), i, s) for i, s in enumerate(valid_signals)]
    tci_candidates.sort(reverse=True)
    best_len, best_idx, tci_raw = tci_candidates[0]
    print(f"   → TCI uses column {best_idx} (length = {best_len:,})")

    tci_extracted = denoise_data(tci_raw, data_type, ligo_idx, cmb_idx)

    # Downsample if excessively long
    if len(tci_extracted) > MAX_LEN:
        idx = np.linspace(0, len(tci_extracted) - 1, MAX_LEN, dtype=int)
        tci_extracted = tci_extracted[idx]

    # --- MCI extraction from remaining columns ---
    remaining = [s for i, s in enumerate(valid_signals) if i != best_idx]

    if data_type in ['ligo', 'cmb', 'qrng'] or not remaining:
        # Use PSD-based MCI
        fs = 1.0
        if data_type == 'ligo':
            fs = CONFIG['ligo'][ligo_idx]['sample_rate']

        f, psd = scipy.signal.welch(
            tci_extracted,
            fs=fs,
            nperseg=min(1024, len(tci_extracted) // 4)
        )
        psd_norm = psd / (np.max(psd) + 1e-10)

        repeat = (len(tci_extracted) // len(psd_norm)) + 1
        mci_extracted = np.tile(psd_norm, (repeat, 1))[:len(tci_extracted), :8]

    else:
        # Use additional channels as measurement features
        mci_list = [denoise_data(s, data_type, ligo_idx, cmb_idx) for s in remaining]
        target_len = len(tci_extracted)

        padded = [
            np.pad(m, (0, target_len - len(m)), mode='constant')
            if len(m) < target_len else m[:target_len]
            for m in mci_list
        ]

        mci_extracted = (
            np.column_stack(padded) if len(padded) > 1
            else padded[0].reshape(-1, 1)
        )

    gc.collect()
    return tci_extracted, mci_extracted

# Dynamic date selection for NANOGrav
def get_pulsar_date(pulsar_name, base_dir, template):
    files = glob.glob(os.path.join(base_dir, template.format(pulsar_name=pulsar_name, date='*')))
    if files:
        try:
            return files[0].split('_PINT_')[1].split('.nb.tim')[0]
        except IndexError:
            print(f"Warning: Cannot extract date for {pulsar_name}, using default")
    return CONFIG['nanograv']['default_date']

# Data loading (full rebuild: snippet-based gaia, 5 cols, aligned dropna, scalar clip)
def load_data(file_path, data_type, ligo_idx=None, cmb_idx=None, pulsar_name=None):
    try:
        if data_type == 'ligo':
            print(f"Loading LIGO file: {file_path}")
            with h5py.File(file_path, 'r') as f:
                if 'strain' not in f or 'Strain' not in f['strain']:
                    raise KeyError("Key 'strain/Strain' not found")
                strain = f['strain']['Strain'][:]
                sos = scipy.signal.butter(2, CONFIG['ligo'][ligo_idx]['freq_range'], btype='band',
                                         fs=CONFIG['ligo'][ligo_idx]['sample_rate'], output='sos')
                strain_filtered = scipy.signal.sosfilt(sos, strain)
                data = da.from_array(strain_filtered, chunks='auto').compute().astype(np.float64)
                for_mfdfa = data
                tci_extracted, mci_extracted = extract_tci_mci_data(data, data_type, ligo_idx=ligo_idx)
                return for_mfdfa, tci_extracted, mci_extracted
        elif data_type == 'cmb':
            print(f"Loading CMB file: {file_path}")
            cmb_maps = hp.read_map(file_path, field=CONFIG['cmb'][cmb_idx]['fields'], verbose=False)
            if isinstance(cmb_maps, np.ndarray):
                cmb_maps = [cmb_maps]
            if CONFIG['cmb'][cmb_idx]['galactic_mask']:
                npix = hp.nside2npix(CONFIG['cmb'][cmb_idx]['nside'])
                mask = np.ones(npix, dtype=bool)
                galactic_pixels = hp.query_strip(CONFIG['cmb'][cmb_idx]['nside'], np.radians(60), np.radians(120))
                mask[galactic_pixels] = False
                cmb_maps = [cmb_map[mask] for cmb_map in cmb_maps]
            data = da.from_array(cmb_maps[0], chunks='auto').compute().astype(np.float64)
            for_mfdfa = data
            tci_extracted, mci_extracted = extract_tci_mci_data(data, data_type, cmb_idx=cmb_idx)
            return for_mfdfa, tci_extracted, mci_extracted
        elif data_type == 'desi':
            print(f"Loading DESI file: {file_path}")
            with fits.open(file_path) as hdul:
                signals = []
                for column in CONFIG['desi']['columns']:
                    data_col = hdul[1].data[column].astype(np.float64)
                    data_col = np.nan_to_num(data_col, nan=np.nanmedian(data_col), posinf=np.nanmedian(data_col), neginf=np.nanmedian(data_col))
                    data_col = data_col[np.isfinite(data_col)]
                    if len(data_col) < 10:
                        print(f"Too few data for DESI column {column}: {len(data_col)}")
                        continue
                    signals.append(data_col)
                if not signals:
                    print("No valid data for DESI")
                    return None, None, None
                for_mfdfa = signals[0]
                min_length = min(len(s) for s in signals)
                signals = [s[:min_length] for s in signals]
                tci_extracted, mci_extracted = extract_tci_mci_data(signals, data_type)
                return for_mfdfa, tci_extracted, mci_extracted
        elif data_type == 'gaia':
            print(f"Loading Gaia DR3 TSV: {file_path}")
            try:
                skip_rows = 257  # Skip to data (post-header/units/dashes, from snippet)
                data = pd.read_csv(
                    file_path,
                    sep='\t',
                    skiprows=skip_rows,
                    nrows=CONFIG['gaia']['subset_size'] * 2,  # *2 for sample, then drop
                    header=None,  # No header
                    engine='python',
                    on_bad_lines='skip',
                    quoting=3,
                    comment=None
                )
                print(f"Loaded Gaia DR3: {len(data)} rows, {len(data.columns)} cols")
                col_map = {
                    'RA_ICRS': 1,
                    'DE_ICRS': 2,
                    'pmRA': 12,  # From snippet
                    'pmDE': 14,
                    'Gmag': 54
                }
                available_cols = []
                signals = []
                for col_name, idx in col_map.items():
                    if idx < len(data.columns):
                        col_data = pd.to_numeric(data.iloc[:, idx], errors='coerce')
                        signals.append(col_data.dropna().values.astype(np.float64))
                        available_cols.append(col_name)
                        print(f"Loaded {col_name}: {len(signals[-1])} values, mean={np.mean(signals[-1]):.3f}")
                    else:
                        print(f"Warning: Col {col_name} (idx {idx}) beyond {len(data.columns)} cols; skip")
                if not signals:
                    raise ValueError("No signals loaded – check snippet indices")
                print(f"Available cols: {available_cols}")
                min_length = min(len(s) for s in signals)
                signals = [s[:min_length] for s in signals]
                data_array = np.column_stack(signals)
                for i, col in enumerate(available_cols):
                    col_data = data_array[:, i]
                    q_low = np.quantile(col_data, 0.01)
                    q_high = np.quantile(col_data, 0.99)
                    data_array[:, i] = np.clip(col_data, q_low, q_high)
                    print(f"Clipped {col}: [{q_low:.3f}, {q_high:.3f}]")
                if 'RA_ICRS' in available_cols and 'DE_ICRS' in available_cols:
                    ra_idx = available_cols.index('RA_ICRS')
                    de_idx = available_cols.index('DE_ICRS')
                    data_array[:, ra_idx] = np.deg2rad(data_array[:, ra_idx])
                    data_array[:, de_idx] = np.deg2rad(data_array[:, de_idx])
                    print("RA/DE converted to radians")
                if min_length < 10:
                    print("Too few valid data for Gaia")
                    return None, None, None
                for_mfdfa = np.mean(data_array, axis=1)
                tci_extracted, mci_extracted = extract_tci_mci_data(data_array, data_type)

                print(f"Gaia loaded: {len(signals)} cols, length {min_length}")
                return for_mfdfa, tci_extracted, mci_extracted

            except Exception as e:
                print(f"Error loading Gaia: {e}")
                import traceback
                traceback.print_exc()
                return None, None, None
        elif data_type == 'cern':
            print(f"Loading CERN file: {file_path}")
            with uproot.open(file_path) as f:
                signals = []
                for column in CONFIG['cern']['columns']:
                    data = f[CONFIG['cern']['tree']][column].array(library='np')
                    if isinstance(data, np.ndarray):
                        if column == 'lep_pt':
                            data = np.concatenate([np.array(x).flatten() for x in data if len(x) > 0])
                        else:
                            data = np.array([np.mean(x) if len(x) > 0 else np.nan for x in data])
                    data = np.nan_to_num(data, nan=np.nanmedian(data), posinf=np.nanmedian(data), neginf=np.nanmedian(data))
                    data = data[np.isfinite(data)]
                    if len(data) < 10:
                        print(f"Not enough data for CERN column {column}: {len(data)}")
                        continue
                    signals.append(data)
                if not signals:
                    print("No valid data for CERN")
                    return None, None, None
                for_mfdfa = signals[0]
                min_length = min(len(s) for s in signals)
                signals = [s[:min_length] for s in signals]
                tci_extracted, mci_extracted = extract_tci_mci_data(signals, data_type)
                return for_mfdfa, tci_extracted, mci_extracted
        elif data_type == 'nist':
            print(f"Loading NIST file: {file_path}")
            dtypes = {
                'sp_num': str, 'Aki(10^8 s^-1)': str, 'wn(cm-1)': str, 'intens': str,
                'tp_ref': str, 'line_ref': str
            }
            data = pd.read_csv(file_path, dtype=dtypes, na_values=['', 'nan', 'NaN', '"'])
            elements = CONFIG['nist']['elements']
            elements_str = '_'.join(elements)  # e.g., 'Ac'
            full_element = ELEMENT_NAMES.get(elements[0], elements[0])  # e.g., 'Actinium'
            print(f"[NIST] Loading for element: {elements_str} ({full_element})")
            data = data[data['element'].isin(elements)]
            signals = []
            for column in CONFIG['nist']['columns']:
                data[column] = data[column].str.strip('="').replace('', np.nan)
                data[column] = pd.to_numeric(data[column], errors='coerce')
                q_low, q_high = data[column].quantile([0.025, 0.975])
                signal = data[(data[column] >= q_low) & (data[column] <= q_high)][column].dropna().values
                if len(signal) < 10:
                    print(f"Too few data for {full_element} ({elements_str}), column {column}: {len(signal)}")
                    continue
                signals.append(signal)
            if not signals:
                print(f"No valid data for {full_element} ({elements_str})")
                return None, None, None, None  # 4-tuple for consistency
            for_mfdfa = signals[0]
            min_length = min(len(s) for s in signals)
            signals = [s[:min_length] for s in signals]
            tci_extracted, mci_extracted = extract_tci_mci_data(signals, data_type)
            return for_mfdfa, tci_extracted, mci_extracted, full_element
        elif data_type == 'nanograv':
            print(f"Loading NANOGrav data for pulsar: {pulsar_name}")
            base_dir = CONFIG['nanograv']['base_dir']
            date = get_pulsar_date(pulsar_name, base_dir, CONFIG['nanograv']['file_templates']['tim'])
            files = {
                'tim': os.path.join(base_dir, CONFIG['nanograv']['file_templates']['tim'].format(pulsar_name=pulsar_name, date=date)),
                'par': os.path.join(base_dir, CONFIG['nanograv']['file_templates']['par'].format(pulsar_name=pulsar_name, date=date)),
                'res_full': os.path.join(base_dir, CONFIG['nanograv']['file_templates']['res_full'].format(pulsar_name=pulsar_name)),
                'res_avg': os.path.join(base_dir, CONFIG['nanograv']['file_templates']['res_avg'].format(pulsar_name=pulsar_name)),
                'dmx': os.path.join(base_dir, CONFIG['nanograv']['file_templates']['dmx'].format(pulsar_name=pulsar_name)),
                'noise': os.path.join(base_dir, CONFIG['nanograv']['file_templates']['noise'].format(pulsar_name=pulsar_name)),
                'clock': os.path.join(base_dir, CONFIG['nanograv']['file_templates']['clock'])
            }
            residuals = None
            for res_file, res_type in [(files['res_full'], 'full'), (files['res_avg'], 'avg')]:
                if os.path.exists(res_file):
                    try:
                        df = pd.read_csv(res_file, delim_whitespace=True, comment='#', header=None)
                        res_col = 1
                        if df.shape[1] <= res_col:
                            print(f"Error: No residual column ({res_col}) in {res_file}")
                            continue
                        residuals = pd.to_numeric(df.iloc[:, res_col], errors='coerce').dropna().values
                        print(f"Loaded {len(residuals)} residuals from {res_file}")
                        break
                    except Exception as e:
                        print(f"Error loading {res_type} residuals: {e}")
            if residuals is None:
                if not os.path.exists(files['tim']) or not os.path.exists(files['par']):
                    print(f"Error: .tim or .par file not found for {pulsar_name}")
                    return None, None, None
                try:
                    model, toas = get_model_and_toas(files['par'], files['tim'], planets=True)
                    print(f"Model components for {pulsar_name}: {list(model.components.keys())}")
                    model.validate()
                    if os.path.exists(files['clock']):
                        with open(files['clock'], 'r') as f:
                            clock_data = pd.read_csv(f, delim_whitespace=True, comment='#', header=None)
                            clock_corrections = clock_data.iloc[:, 1].values
                            toas.adjust_TOAs(clock_corrections * 1e6)
                    residuals = model.residuals(toas, use_weighted_mean=False).to('s').value
                    residuals = residuals[np.isfinite(residuals)]
                    print(f"Loaded residuals: {len(residuals)} points")
                except Exception as e:
                    print(f"Error loading with pint: {e}")
                    return None, None, None
            if len(residuals) < CONFIG['nanograv']['subset_size']:
                print(f"Error: Insufficient residuals: {len(residuals)} < {CONFIG['nanograv']['subset_size']}")
                return None, None, None
            dmx_data = None
            if os.path.exists(files['dmx']):
                try:
                    df_dmx = pd.read_csv(files['dmx'], delim_whitespace=True, comment='#', header=None)
                    dmx_data = pd.to_numeric(df_dmx.iloc[:, 1], errors='coerce').dropna().values
                    print(f"Loaded DMX data: {len(dmx_data)} points")
                except Exception as e:
                    print(f"Error loading DMX data: {e}")
            red_noise = None
            if os.path.exists(files['noise']):
                try:
                    with open(files['noise'], 'r') as f:
                        noise_lines = f.readlines()
                    for line in noise_lines:
                        if 'EFAC' in line or 'EQUAD' in line or 'ECORR' in line:
                            red_noise = float(line.split()[1]) if red_noise is None else red_noise
                    print(f"Loaded red noise parameter: {red_noise}")
                except Exception as e:
                    print(f"Error loading noise parameters: {e}")
            data_dict = {'residuals': residuals}
            if dmx_data is not None:
                data_dict['dmx'] = dmx_data
            if red_noise is not None:
                data_dict['red_noise'] = np.full_like(residuals, red_noise)
            min_length = min(len(data) for data in data_dict.values())
            data_array = np.column_stack([data[:min_length] for data in data_dict.values()])
            for_mfdfa = residuals
            tci_extracted, mci_extracted = extract_tci_mci_data(data_dict, data_type)
            return for_mfdfa, tci_extracted, mci_extracted
        elif data_type == 'qrng':
            print("Loading LFDR QRNG via API...")
            try:
                data_list = []
                for i in range(50):  # 50 calls for ~100k bits
                    for retry in range(3):
                        url = 'https://lfdr.de/qrng_api/qrng?length=256&format=HEX'
                        response = requests.get(url, timeout=5)
                        if response.status_code == 200:
                            try:
                                json_data = response.json()
                                hex_string = json_data['qrn']
                                bits = []
                                for char in hex_string:
                                    byte = int(char, 16)
                                    bits.extend([int(b) for b in f"{byte:04b}"])
                                data_list.extend(bits)
                                print(f"QRNG call {i+1}: {len(bits)} bits geladen")
                                break
                            except (json.JSONDecodeError, KeyError) as je:
                                print(f"JSON/Key error, retry {retry+1}: {je}. Response: {response.text[:100]}...")
                                if retry == 2:
                                    raise ValueError("API parse failed")
                        else:
                            print(f"HTTP {response.status_code}, retry {retry+1}. Response: {response.text[:100]}...")
                            if retry < 2:
                                time.sleep(1)
                            else:
                                raise ValueError("API failed after retries")
                if len(data_list) < 100:
                    raise ValueError(f"Not enough QRNG-data: {len(data_list)} bits")
                data = np.array(data_list).astype(np.float64)
                data += np.random.normal(0, 1e-6, len(data))
                print(f"Geladen {len(data)} quantum random bits")
                for_mfdfa = data
                tci_extracted, mci_extracted = extract_tci_mci_data(data, data_type)
                return for_mfdfa, tci_extracted, mci_extracted
            except Exception as e:
                print(f"Error loading QRNG (fallback to simulation): {e}")
                data = np.random.randint(0, 2, 102400).astype(np.float64)
                data += np.random.normal(0, 1e-6, len(data))
                print(f"Fallback: Simulated {len(data)} random samples")
                for_mfdfa = data
                tci_extracted, mci_extracted = extract_tci_mci_data(data, data_type)
                return for_mfdfa, tci_extracted, mci_extracted
        else:
            raise ValueError("Invalid data_type")
    except Exception as e:
        print(f"Error loading {data_type}: {e}")
        return None, None, None

# Z-test validation for MFDFA
def validate_df(D_f_values, expected_D_f, sigma_D_f):
    valid_D_f = [x for x in D_f_values if np.isfinite(x)]
    if not valid_D_f:
        print("No valid D_f values for validation")
        return np.nan
    mean_D_f = np.nanmean(valid_D_f)
    std_D_f = np.nanstd(valid_D_f)
    n = len(valid_D_f)
    if std_D_f < 1e-10 or n < 2:
        print("Insufficient variability or data for Z-test")
        return np.nan
    z_score = (mean_D_f - expected_D_f) / (std_D_f / np.sqrt(n))
    p_value = 2 * (1 - norm.cdf(abs(z_score)))
    return p_value

# Check reliability of h_q
def is_valid_hq(hq, q_values, min_valid_ratio=0.6, monotonicity_tolerance=0.05, require_monotonicity=True):
    if isinstance(hq, np.floating):
        return False
    valid_hq = np.isfinite(hq)
    valid_ratio = np.sum(valid_hq) / len(hq)
    if valid_ratio < min_valid_ratio:
        return False
    if not require_monotonicity:
        return True
    valid_indices = np.where(valid_hq)[0]
    if len(valid_indices) < 2:
        return False
    hq_valid = hq[valid_indices]
    q_valid = q_values[valid_indices]
    diff_hq = np.diff(hq_valid) / np.diff(q_valid)
    monotonic = np.all(diff_hq <= monotonicity_tolerance)
    return monotonic

# Subset processing - Enhanced logging for every subset
def process_subset(subset_idx, data, data_type, dataset_name, scales, ligo_idx=None, cmb_idx=None):
    """
    Processes a single subset: Extract → Denoise → MFDFA.
    Returns: D_f, hq, fluct, slopes, subset_data_denoised (no logging here - done in process_dataset)
    """
    try:
        if data_type == 'ligo':
            n_samples = int(CONFIG['ligo'][ligo_idx]['subset_duration'] * CONFIG['ligo'][ligo_idx]['sample_rate'])
            max_start = len(data) - n_samples
            subset_start = np.random.randint(0, max_start + 1)
            subset_data = data[subset_start:subset_start + n_samples]
            if len(subset_data) != n_samples:
                subset_data = np.pad(subset_data, (0, n_samples - len(subset_data)), mode='constant')
        elif data_type == 'cmb':
            nside = CONFIG['cmb'][cmb_idx]['nside']
            subset_size = CONFIG['cmb'][cmb_idx]['subset_size']
            npix = len(data)
            valid_indices = np.arange(npix)
            center_pix = np.random.choice(valid_indices)
            try:
                subset_indices = hp.query_disc(nside, hp.pix2vec(nside, center_pix, nest=False), radius=CONFIG['cmb'][cmb_idx]['disc_radius'])
                valid_mask = subset_indices < npix
                subset_indices = subset_indices[valid_mask]
                if len(subset_indices) < subset_size:
                    subset_indices = np.random.choice(valid_indices, size=subset_size, replace=False)
                elif len(subset_indices) > subset_size:
                    subset_indices = np.random.choice(subset_indices, size=subset_size, replace=False)
                subset_data = data[subset_indices]
                if np.std(subset_data) < CONFIG['cmb'][cmb_idx]['min_std']:
                    subset_indices = np.random.choice(valid_indices, size=subset_size, replace=False)
                    subset_data = data[subset_indices]
            except Exception:
                subset_indices = np.random.choice(valid_indices, size=subset_size, replace=False)
                subset_data = data[subset_indices]
        elif data_type in ['desi', 'cern', 'nist', 'nanograv', 'qrng']:
            subset_size = CONFIG[data_type]['subset_size'] if 'subset_size' in CONFIG[data_type] else len(data) // 4
            if data_type == 'nist':
                subset_size = CONFIG[data_type]['subset_size'](len(data))
            max_start = len(data) - subset_size
            if max_start < 0:
                print(f"Error: Data too short for {dataset_name} (length={len(data)}, required={subset_size})")
                q_len = len(CONFIG['mfdfa']['q_values'])
                return np.nan, np.full(q_len, np.nan), np.full((q_len, len(scales)), np.nan), np.full(q_len, np.nan), None
            subset_start = np.random.randint(0, max_start + 1)
            subset_data = data[subset_start:subset_start + subset_size]
            if len(subset_data) != subset_size:
                subset_data = np.pad(subset_data, (0, subset_size - len(subset_data)), mode='constant')
        elif data_type == 'gaia':
            subset_size = CONFIG['gaia']['subset_size']
            max_start = len(data) - subset_size
            if max_start < 0:
                q_len = len(CONFIG['mfdfa']['q_values'])
                return np.nan, np.full(q_len, np.nan), np.full((q_len, len(scales)), np.nan), np.full(q_len, np.nan), None
            subset_start = np.random.randint(0, max_start + 1)
            subset_data = data[subset_start:subset_start + subset_size]
            if len(subset_data) != subset_size:
                subset_data = np.pad(subset_data, (0, subset_size - len(subset_data)), mode='constant', constant_values=0)
        else:
            raise ValueError("Invalid data_type")

        # Denoise
        subset_data_denoised = denoise_data(subset_data, data_type, ligo_idx, cmb_idx)

        # MFDFA input (for gaia: pairwise dists)
        if data_type == 'gaia' and subset_data_denoised.ndim > 1:
            ra = subset_data_denoised[:, 0]  # RA
            de = subset_data_denoised[:, 1]  # DE
            pmra = subset_data_denoised[:, 2]  # pmRA
            pmde = subset_data_denoised[:, 3]  # pmDE
            sample_idx = np.random.choice(len(ra), 50, replace=False)
            ra_s, de_s, pmra_s, pmde_s = ra[sample_idx], de[sample_idx], pmra[sample_idx], pmde[sample_idx]
            dists = []
            for i in range(len(ra_s)):
                for j in range(i+1, len(ra_s)):
                    d_pos = np.sqrt((ra_s[i] - ra_s[j])**2 + (de_s[i] - de_s[j])**2)
                    d_vel = np.sqrt((pmra_s[i] - pmra_s[j])**2 + (pmde_s[i] - pmde_s[j])**2)
                    d = np.sqrt(d_pos**2 + d_vel**2)  # 3D (pos + vel)
                    dists.append(d)
            dists = np.array(dists)
            if len(dists) < 100:
                dists = np.pad(dists, (0, 100 - len(dists)), mode='constant')
            mfdfa_input = np.log(dists + 1)
        else:
            mfdfa_input = subset_data_denoised

        # MFDFA
        D_f, hq, rms_values, fluct, slopes = jedi_mfdfa(mfdfa_input, scales, CONFIG['mfdfa']['q_values'], CONFIG['mfdfa']['detrend_order'])

        return D_f, hq, fluct, slopes, subset_data_denoised
    except Exception as e:
        metadata_log['errors'].append({'subset': subset_idx, 'dataset': dataset_name, 'error': str(e)})  # Dit is main-safe, want append na return
        # PATCH 3
        if dataset_name not in metadata_log["subset_warnings"]:
            metadata_log["subset_warnings"][dataset_name] = []
        metadata_log["subset_warnings"][dataset_name].append({
            "subset_idx": subset_idx,
            "error": str(e)
        })

        print(f"Error in subset {subset_idx+1} of {dataset_name}: {e}")
        q_len = len(CONFIG['mfdfa']['q_values'])
        return np.nan, np.full(q_len, np.nan), np.full((q_len, len(scales)), np.nan), np.full(q_len, np.nan), None
    finally:
        gc.collect()

def process_dataset(data, data_type, dataset_name, scales, expected_D_f, sigma_D_f,
                    ligo_idx=None, cmb_idx=None):
    """
    Run MFDFA subset processing for a single dataset.

    This function:
        • Dynamically determines subset size and scales for NIST (v5.1 logic)
        • Applies fixed configuration-based settings for all other data types
        • Executes subset-level MFDFA in parallel
        • Computes robust mean D_f, std, and mean h(q)
        • Performs a Z-test relative to expected_D_f
        • Logs all results into metadata_log['subsets_processed']

    Returns
    -------
    results_dict : dict or None
        Summary of MFDFA metrics for the dataset, or None if invalid.
    subset_data_list : list
        List of raw subset results returned by process_subset().
    """

    # ================================================================
    # 0. Input Validation
    # ================================================================
    if data is None:
        print(f"No valid data for {dataset_name}. Skipping.")
        metadata_log['errors'].append({'dataset': dataset_name, 'error': 'No data'})
        return None, []

    # ================================================================
    # 1. Dynamic handling for NIST datasets (small-sequence stabilisation)
    # ================================================================
    if data_type == 'nist':
        n_total = len(data)

        # Adaptive subset size (~20% of series, lower bound = 40)
        subset_size = max(40, n_total // 5)

        # Adaptive scales (fallback if config provided function/lambda)
        if isinstance(scales, (list, np.ndarray)):
            pass
        else:
            scales = np.logspace(
                np.log10(4),
                np.log10(max(8, n_total // 4)),
                15
            ).astype(int)

        # Maximum number of effectively non-overlapping subsets
        max_nonoverlap = max(1, n_total // subset_size)

        # Allow ~3× coverage for robustness, with a hard cap
        max_reasonable = max_nonoverlap * 3

        if n_total < 500:
            # Very small NIST sequences: keep subset count modest
            n_subsets = min(30, max(10, max_reasonable))
        else:
            # Larger NIST sequences: use config upper bound
            n_subsets = min(CONFIG['nist']['n_subsets'], max(10, max_reasonable))

        print(
            f"   → NIST dynamic configuration: {n_total:,} samples → "
            f"subset_size={subset_size}, scales={len(scales)}, "
            f"n_subsets={n_subsets} (max_nonoverlap={max_nonoverlap})"
        )

    # ================================================================
    # 2. Fixed configuration logic for all other dataset types
    # ================================================================
    else:
        if data_type == 'ligo':
            n_subsets = CONFIG['ligo'][ligo_idx]['n_subsets']
        elif data_type == 'cmb':
            n_subsets = CONFIG['cmb'][cmb_idx]['n_subsets']
        elif data_type == 'gaia':
            n_subsets = CONFIG['gaia']['n_subsets']
        else:
            n_subsets = CONFIG[data_type]['n_subsets']

    # ================================================================
    # 3. Parallel subset-level processing (MFDFA)
    # ================================================================
    results = Parallel(n_jobs=-1)(
        delayed(process_subset)(
            subset_idx, data, data_type, dataset_name, scales,
            ligo_idx=ligo_idx, cmb_idx=cmb_idx
        )
        for subset_idx in tqdm(range(n_subsets), desc=f"Subsets for {dataset_name}")
    )

    # Unpack results
    D_f_values       = [r[0] for r in results]
    hq_values        = [r[1] for r in results if r[1] is not None]
    fluct_values     = [r[2] for r in results]
    slopes_values    = [r[3] for r in results]
    subset_data_list = [r[4] for r in results]

    # ================================================================
    # 4. Summary and Statistical Evaluation
    # ================================================================
    print(f"\nSummary for {dataset_name}:")

    if D_f_values and np.any(np.isfinite(D_f_values)):
        D_f_arr   = np.array(D_f_values, dtype=float)
        D_f_valid = D_f_arr[np.isfinite(D_f_arr)]

        if len(D_f_valid) == 0:
            mean_D_f = np.nan
            std_D_f  = np.nan
        else:
            # Robust core selection using median ± 3*MAD
            med = np.nanmedian(D_f_valid)
            mad = np.nanmedian(np.abs(D_f_valid - med)) + 1e-12
            mask = np.abs(D_f_valid - med) <= 3 * mad
            core = D_f_valid[mask]

            # If filtering removes too many subsets, fallback to all valid values
            if len(core) < max(3, len(D_f_valid) // 3):
                core = D_f_valid

            mean_D_f = float(np.nanmean(core))
            std_D_f  = float(np.nanstd(core))

        # Mean h(q) across valid subsets
        valid_hq_means = [np.nanmean(hq) for hq in hq_values if np.any(np.isfinite(hq))]
        mean_hq = np.nanmean(valid_hq_means) if valid_hq_means else np.nan

        # Z-test against expected D_f
        p_value = validate_df(D_f_valid, expected_D_f, sigma_D_f)

        # Console summary (kept concise)
        print(f"Mean D_f: {mean_D_f:.3f} ± {std_D_f:.3f}")
        print(f"Expected D_f: {expected_D_f} ± {sigma_D_f}")
        print(f"|Δ D_f| = {abs(mean_D_f - expected_D_f):.3f}")
        print(f"Mean h(q): {mean_hq:.3f}")
        print(f"Valid subsets: {np.sum(np.isfinite(D_f_arr))}")
        print(f"Z-test p-value: {p_value:.3f}")

        # ============================================================
        # 5. Metadata Logging (JSON + CSV compatibility)
        # ============================================================
        metadata_log['subsets_processed'][dataset_name] = {
            'n_subsets': n_subsets,
            'D_f_values': D_f_values,
            'hq_values': [
                hq.tolist() if np.any(np.isfinite(hq)) else [] 
                for hq in hq_values
            ],
            'fluct_values': [
                f.tolist() if np.any(np.isfinite(f)) else [] 
                for f in fluct_values
            ],
            'slopes_values': [
                s.tolist() if np.any(np.isfinite(s)) else [] 
                for s in slopes_values
            ],
            'mean_D_f': float(mean_D_f),
            'std_D_f': float(std_D_f),
            'mean_hq': float(mean_hq),
            'p_value': float(p_value),   # Needed by NANOGrav pipeline
            'n_valid_subsets': int(np.sum(np.isfinite(D_f_arr))),
        }

        return {
            'D_f': D_f_values,
            'hq_values': hq_values,
            'fluct': fluct_values,
            'slopes': slopes_values,
            'mean_D_f': float(mean_D_f),
            'std_D_f': float(std_D_f),
            'mean_hq': float(mean_hq),
            'p_value': float(p_value),
            'best_subsets': subset_data_list
        }, subset_data_list

    # ================================================================
    # 6. Failure Case: No valid D_f
    # ================================================================
    print("No valid D_f values computed.")
    metadata_log['subsets_processed'][dataset_name] = {
        'n_subsets': n_subsets,
        'mean_D_f': np.nan,
        'std_D_f': np.nan,
        'mean_hq': np.nan,
        'p_value': np.nan
    }
    metadata_log['errors'].append({'dataset': dataset_name, 'error': 'No valid D_f'})
    return None, []

def convert_numpy(obj):
    """
    Recursive convert NumPy-types to Python-natives for JSON-serialisatie.
    """
    if isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: convert_numpy(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy(item) for item in obj]
    else:
        return obj

def save_run_metadata(save_flag):
    if not save_flag:
        return

    # 1. Light CSV
    df = pd.DataFrame(columns=['type', 'pair', 'dataset', 'subset_idx', 'D_f', 'hq_mean',
                               'corr', 'wavelet_corr', 'weight', 'tci_global', 'mci_global',
                               'fold_idx', 'tmci_value', 'tmci_global', 'tmci_std_cv'])

    # TCI
    for pair, d in metadata_log.get('tci_pairs', {}).items():
        df = pd.concat([df, pd.DataFrame([{
            'type': 'TCI', 'pair': pair, 'corr': d.get('corr'),
            'wavelet_corr': d.get('wavelet_corr'), 'weight': d.get('weight')
        }])], ignore_index=True)

    # MCI
    for item in metadata_log.get('mci_measurements', []):
        df = pd.concat([df, pd.DataFrame([{
            'type': 'MCI', 'pair_idx': item.get('pair_idx'),
            'corr': item.get('corr'), 'psd_corr': item.get('psd_corr')
        }])], ignore_index=True)

    # TMCI folds
    for i, val in enumerate(metadata_log.get('tmci_folds', [])):
        df = pd.concat([df, pd.DataFrame([{
            'type': 'TMCI_fold', 'fold_idx': i, 'tmci_value': val
        }])], ignore_index=True)

    # Max number of subsets
    subset_counter = 0
    for dataset, proc in metadata_log.get('subsets_processed', {}).items():
        for i in range(len(proc.get('D_f_values', []))):
            if subset_counter >= 5000:
                break
            hq_mean = np.nanmean(proc['hq_values'][i]) if i < len(proc['hq_values']) else np.nan
            df = pd.concat([df, pd.DataFrame([{
                'type': 'MFDFA_subset', 'dataset': dataset, 'subset_idx': i,
                'D_f': proc['D_f_values'][i], 'hq_mean': hq_mean
            }])], ignore_index=True)
            subset_counter += 1
        if subset_counter >= 5000:
            break

    # Global metrics
    df = pd.concat([df, pd.DataFrame([{
        'type': 'GLOBAL',
        'tci_global': metadata_log.get('tci'),
        'mci_global': metadata_log.get('mci'),
        'tmci_global': metadata_log.get('tmci'),
        'tmci_std_cv': metadata_log.get('tmci_std_cv')
    }])], ignore_index=True)

    timestamp = metadata_log['run_timestamp']
    csv_path = f"{OUTPUT_DIR}full_melted_utmf_v5.1_{timestamp}.csv"
    df.to_csv(csv_path, index=False)
    print(f"→ Light CSV saved ({len(df)} rows): {csv_path}")

# ================================================================
# Temporal Coherence Index (TCI)
# ================================================================
def calculate_tci_multivariate(tci_datasets):
    """
    Compute the Temporal Coherence Index (TCI) across all datasets.

    TCI measures the mean absolute Pearson correlation between the
    multifractal h(q) spectra of multiple time-series. Each series is:

        • trimmed to ≤ 6144 samples
        • analysed using the embedded jedi_mfdfa()
        • interpolated if h(q) contains NaNs
        • compared against every other dataset

    Parameters
    ----------
    tci_datasets : dict[str, ndarray]
        Mapping from dataset names to extracted TCI time series.

    Returns
    -------
    tci : float
        Temporal coherence index (mean |corr| across all pairs).
    tci_meta : dict
        Placeholder for compatibility (v5.1 logs directly to metadata_log).
    """
    spectra = []
    valid_names = []

    for name, series in tci_datasets.items():
        if series is None or len(series) < 1000:
            continue

        # Reduce multi-column inputs
        if isinstance(series, np.ndarray) and series.ndim > 1:
            series = np.mean(series, axis=1)

        # Limit maximum analysis length
        series = series[:6144]

        try:
            scales = np.logspace(
                np.log10(8),
                np.log10(len(series) // 8),
                15
            ).astype(int)
            scales = scales[scales > 4]
            q_values = CONFIG['mfdfa']['q_values']

            D_f, hq, _, _, _ = jedi_mfdfa(series, scales, q_values)

            # Interpolate missing values in h(q)
            if np.any(np.isnan(hq)):
                valid = np.isfinite(hq)
                if np.sum(valid) < 2:
                    continue
                hq = np.interp(q_values, q_values[valid], hq[valid])

            spectra.append(hq)
            valid_names.append(name)
            print(f"[TCI] {name}: h(q) computed successfully")
        except Exception:
            continue

    # Require >1 dataset to compute TCI
    if len(spectra) < 2:
        print(f"[TCI] Insufficient spectra ({len(spectra)}) → TCI = 0.0")
        return 0.0, {}

    # Pairwise |corr| over all h(q)
    all_corrs = []
    for i in range(len(spectra)):
        for j in range(i + 1, len(spectra)):
            corr, _ = pearsonr(spectra[i], spectra[j])
            all_corrs.append(abs(corr))

    tci = np.mean(all_corrs)
    print(f"[TCI] Fractal TCI = {tci:.3f} ({len(all_corrs)} pairs from {len(valid_names)} datasets)")

    # Metadata logging (v5.1 JSON snapshot)
    for name, spec in zip(valid_names, spectra):
        metadata_log["tci_meta"][name] = {
            "hq_length": len(spec),
            "valid": True
        }

    return tci, {}

#================================================================
# Measurement Coherence Index (MCI)
#================================================================
def calculate_mci_multivariate(results_all):
    """
    Compute the Measurement Coherence Index (MCI) across datasets.

    MCI is the mean absolute correlation between the *mean* h(q) spectra
    derived from the MFDFA subset statistics of each dataset.

    Parameters
    ----------
    results_all : dict[str, dict]
        Contains MFDFA results for all datasets:
        each entry must include 'hq_values' (list of h(q) arrays).

    Returns
    -------
    mci : float
        Measurement coherence index.
    """
    spectra = []
    valid_datasets = []

    for dataset_name, res in results_all.items():
        # Allow tuple outputs of process_dataset
        if isinstance(res, tuple) and len(res) > 0:
            res = res[0]

        if not isinstance(res, dict):
            continue
        if 'hq_values' not in res:
            continue

        hq_list = res['hq_values']
        if len(hq_list) == 0:
            continue

        try:
            hq_array = np.stack(hq_list)          # shape: (n_subsets, n_q)
            mean_hq = np.nanmean(hq_array, axis=0)

            if len(mean_hq) < 5 or np.all(~np.isfinite(mean_hq)):
                continue

            spectra.append(mean_hq)
            valid_datasets.append(dataset_name)
        except Exception:
            continue

    print(f"[MCI] Datasets with valid h(q): {len(valid_datasets)}")

    if len(spectra) < 2:
        print("[MCI] Fewer than 2 valid spectra → MCI = NaN")
        return float("nan")

    all_corrs = []
    for i in range(len(spectra)):
        for j in range(i + 1, len(spectra)):
            a = np.asarray(spectra[i])
            b = np.asarray(spectra[j])

            mask = np.isfinite(a) & np.isfinite(b)
            if np.sum(mask) < 2:
                continue

            corr, _ = pearsonr(a[mask], b[mask])
            if np.isfinite(corr):
                all_corrs.append(abs(corr))

    if not all_corrs:
        print("[MCI] No finite pairwise correlations → MCI = NaN")
        mci = float("nan")
    else:
        mci = float(np.mean(all_corrs))

    print(f"[MCI] Measurement MCI = {mci:.3f} ({len(all_corrs)} pairs from {len(valid_datasets)} datasets)")

    metadata_log["mci_meta"] = {
        "datasets": valid_datasets,
        "n_pairs": len(all_corrs)
    }

    return mci

#================================================================
# Temporal–Measurement Coherence Index (TMCI)
#================================================================
def calculate_tmci_empirical(tci, mci, n_boot=500):
    """
    Compute the empirical TMCI metric:
        • bootstrap resampling (Gaussian jitter)
        • correlation between TCI and MCI
        • adaptive weighting based on redundancy

    The weight scheme:
        w = (1 − |corr|) / 2
    ensures that:
        • when TCI and MCI are independent → equal weight
        • when they are highly correlated → redundancy reduces influence

    Parameters
    ----------
    tci : float
        Temporal Coherence Index.
    mci : float
        Measurement Coherence Index.
    n_boot : int
        Number of bootstrap samples.

    Returns
    -------
    tmci_mean : float
    tmci_ci : tuple(float, float)
        95% confidence interval of TMCI samples.
    corr_tm : float
        Correlation between bootstrapped TCI and MCI.
    tci_samples, mci_samples, tmci_samples : ndarray
    """
    # Edge case: undefined TMCI
    if np.isnan(tci) or np.isnan(mci):
        return (
            np.nan,
            (np.nan, np.nan),
            np.nan,
            np.array([]),
            np.array([]),
            np.array([])
        )

    # Bootstrap jitter
    tci_samples = tci + np.random.normal(0, 0.02, n_boot)
    mci_samples = mci + np.random.normal(0, 0.02, n_boot)

    # Correlation between the bootstrapped sets
    corr_tm, _ = pearsonr(tci_samples, mci_samples)

    # Adaptive weighting
    w = (1 - abs(corr_tm)) / 2

    tmci_samples = (
        w * tci_samples +
        w * mci_samples +
        (1 - 2 * w) * ((tci_samples + mci_samples) / 2)
    )

    tmci_mean = np.mean(tmci_samples)
    tmci_ci   = (
        np.percentile(tmci_samples, 2.5),
        np.percentile(tmci_samples, 97.5)
    )

    return tmci_mean, tmci_ci, corr_tm, tci_samples, mci_samples, tmci_samples

def load_utmf_data():
    """
    Load and preprocess all UTMF datasets across domains (LIGO, CMB, DESI, CERN, NIST, NANOGrav, QRNG, Gaia).

    For each dataset that has `utmf_use = True`:
        • load_data() extracts data for MFDFA, TCI, and MCI
        • process_dataset() runs MFDFA subset analysis
        • results are written to:
              - results_all      (MFDFA summary)
              - tci_datasets     (TCI time series)
              - mci_datasets     (measurement matrices)
              - metadata_log     (diagnostics, shapes, subset info)

    Returns
    -------
    results_all : dict
        Mapping: dataset_name → MFDFA summary dictionary.
    tci_datasets : dict
        Mapping: dataset_name → extracted TCI time series.
    mci_datasets : dict
        Mapping: dataset_name → extracted measurement vectors.
    """
    results_all = {}
    tci_datasets = {}
    mci_datasets = {}

    # Store a shallow configuration snapshot
    metadata_log['config_snapshot'] = CONFIG.copy()

    # ================================================================
    #                       L I G O   (Strain)
    # ================================================================
    for idx, (file_path, dataset_name) in enumerate(zip(CONFIG['ligo_files'],
                                                        CONFIG['ligo_names'])):
        if not CONFIG['ligo'][idx]['utmf_use']:
            continue

        print(f"\n[UTMF] Loading LIGO dataset: {dataset_name}")
        for_mfdfa, tci_extracted, mci_extracted = load_data(
            file_path, 'ligo', ligo_idx=idx
        )

        if for_mfdfa is None:
            print(f"   → No valid data for {dataset_name}; skipping.")
            metadata_log['errors'].append(
                {'dataset': dataset_name, 'error': 'No data'}
            )
            continue

        # Run MFDFA subset pipeline
        results, _ = process_dataset(
            for_mfdfa,
            'ligo',
            dataset_name,
            CONFIG['ligo'][idx]['scales'],
            CONFIG['ligo'][idx]['expected_D_f'],
            CONFIG['ligo'][idx]['sigma_D_f'],
            ligo_idx=idx
        )

        if results:
            results_all[dataset_name] = results
            tci_datasets[dataset_name] = tci_extracted
            mci_datasets[dataset_name] = mci_extracted

            if tci_extracted is not None and len(tci_extracted) >= 12288:
                print(f"   → TCI {dataset_name}: using longest valid column "
                      f"({len(tci_extracted):,} samples, LIGO strain)")
            else:
                print(f"   → TCI {dataset_name}: no suitable column → excluded from TCI")

            metadata_log['datasets_loaded'][dataset_name] = {
                'n_subsets': CONFIG['ligo'][idx]['n_subsets'],
                'status': 'OK'
            }
            metadata_log['load_meta'][dataset_name] = {
                'n_rows': len(for_mfdfa),
                'tci_length': len(tci_extracted) if tci_extracted is not None else None,
                'mci_shape': (
                    list(mci_extracted.shape)
                    if isinstance(mci_extracted, np.ndarray) else None
                )
            }

    gc.collect()

    # ================================================================
    #                    C M B   /   P L A N C K
    # ================================================================
    for idx, (file_path, dataset_name) in enumerate(zip(CONFIG['cmb_files'],
                                                        CONFIG['cmb_names'])):
        if not CONFIG['cmb'][idx]['utmf_use']:
            continue

        print(f"\n[UTMF] Loading CMB dataset: {dataset_name}")
        for_mfdfa, tci_extracted, mci_extracted = load_data(
            file_path, 'cmb', cmb_idx=idx
        )

        if for_mfdfa is None:
            print(f"   → No valid data for {dataset_name}; skipping.")
            metadata_log['errors'].append(
                {'dataset': dataset_name, 'error': 'No data'}
            )
            continue

        results, _ = process_dataset(
            for_mfdfa,
            'cmb',
            dataset_name,
            CONFIG['cmb'][idx]['scales'],
            CONFIG['cmb'][idx]['expected_D_f'],
            CONFIG['cmb'][idx]['sigma_D_f'],
            cmb_idx=idx
        )

        if results:
            results_all[dataset_name] = results
            tci_datasets[dataset_name] = tci_extracted
            mci_datasets[dataset_name] = mci_extracted

            if tci_extracted is not None and len(tci_extracted) >= 12288:
                print(f"   → TCI {dataset_name}: using longest valid column "
                      f"({len(tci_extracted):,} samples, CMB map)")
            else:
                print(f"   → TCI {dataset_name}: no suitable column → excluded from TCI")

            metadata_log['datasets_loaded'][dataset_name] = {
                'n_subsets': CONFIG['cmb'][idx]['n_subsets'],
                'status': 'OK'
            }

            metadata_log['load_meta'][dataset_name] = {
                'n_rows': len(for_mfdfa),
                'tci_length': len(tci_extracted) if tci_extracted is not None else None,
                'mci_shape': (
                    list(mci_extracted.shape)
                    if isinstance(mci_extracted, np.ndarray) else None
                )
            }

    gc.collect()

    # ================================================================
    #                            D E S I
    # ================================================================
    if CONFIG['desi'].get('utmf_use', False):
        dataset_name = CONFIG['desi_name']
        print(f"\n[UTMF] Loading DESI dataset: {dataset_name}")

        for_mfdfa, tci_extracted, mci_extracted = load_data(
            CONFIG['desi_file'], 'desi'
        )

        if for_mfdfa is not None:
            results, _ = process_dataset(
                for_mfdfa,
                'desi',
                dataset_name,
                CONFIG['desi']['scales'],
                CONFIG['desi']['expected_D_f'],
                CONFIG['desi']['sigma_D_f']
            )

            if results:
                results_all[dataset_name] = results
                tci_datasets[dataset_name] = tci_extracted
                mci_datasets[dataset_name] = mci_extracted

                if tci_extracted is not None and len(tci_extracted) >= 12288:
                    print(f"   → TCI {dataset_name}: using longest valid column "
                          f"({len(tci_extracted):,} samples)")
                else:
                    print(f"   → TCI {dataset_name}: no suitable column → excluded from TCI")

                metadata_log['datasets_loaded'][dataset_name] = {
                    'n_subsets': CONFIG['desi']['n_subsets'],
                    'status': 'OK'
                }
                metadata_log['load_meta'][dataset_name] = {
                    'n_rows': len(for_mfdfa),
                    'tci_length': len(tci_extracted) if tci_extracted is not None else None,
                    'mci_shape': (
                        list(mci_extracted.shape)
                        if isinstance(mci_extracted, np.ndarray) else None
                    )
                }

    gc.collect()

    # ================================================================
    #                             C E R N
    # ================================================================
    if CONFIG['cern'].get('utmf_use', False):
        dataset_name = CONFIG['cern_name']
        print(f"\n[UTMF] Loading CERN dataset: {dataset_name}")

        for_mfdfa, tci_extracted, mci_extracted = load_data(
            CONFIG['cern_file'], 'cern'
        )

        if for_mfdfa is not None:
            results, _ = process_dataset(
                for_mfdfa,
                'cern',
                dataset_name,
                CONFIG['cern']['scales'],
                CONFIG['cern']['expected_D_f'],
                CONFIG['cern']['sigma_D_f']
            )

            if results:
                results_all[dataset_name] = results
                tci_datasets[dataset_name] = tci_extracted
                mci_datasets[dataset_name] = mci_extracted

                if tci_extracted is not None and len(tci_extracted) >= 12288:
                    print(f"   → TCI {dataset_name}: using longest valid column "
                          f"({len(tci_extracted):,} samples)")
                else:
                    print(f"   → TCI {dataset_name}: no suitable column → excluded from TCI")

                metadata_log['datasets_loaded'][dataset_name] = {
                    'n_subsets': CONFIG['cern']['n_subsets'],
                    'status': 'OK'
                }
                metadata_log['load_meta'][dataset_name] = {
                    'n_rows': len(for_mfdfa),
                    'tci_length': len(tci_extracted) if tci_extracted is not None else None,
                    'mci_shape': (
                        list(mci_extracted.shape)
                        if isinstance(mci_extracted, np.ndarray) else None
                    )
                }

    gc.collect()

    # ================================================================
    #                    N I S T   E L E M E N T S
    #      (Robust parsing, outlier-trimming, heavy GC per element)
    # ================================================================
    if CONFIG['nist'].get('utmf_use', False):
        print(f"\n[UTMF] Processing {len(CONFIG['nist']['elements_list_utmf'])} NIST elements "
              f"(robust parsing + aggressive memory cleanup).")

        NIST_ELEMENTS_DIR = "/content/drive/MyDrive/Datasets_UTMF/NIST_elements/"
        n_elements = len(CONFIG['nist']['elements_list_utmf'])

        for elem_idx, elem_list in enumerate(CONFIG['nist']['elements_list_utmf']):
            elements_str = "_".join(elem_list)
            full_element_name = ELEMENT_NAMES.get(elem_list[0], elem_list[0])
            dataset_name = f"NIST_3_{elements_str}_{full_element_name}"
            csv_path = os.path.join(NIST_ELEMENTS_DIR, f"NIST_{elem_list[0]}.csv")

            # Predeclare local vars for cleanup
            data = None
            signals = []
            for_mfdfa = None
            tci_extracted = None
            mci_extracted = None
            results = None

            print(f"\n→ [{elem_idx+1}/{n_elements}] Processing: {elements_str.ljust(4)} ({full_element_name})")

            if not os.path.exists(csv_path):
                print(f"   → File not found: {csv_path}")
                gc.collect()
                continue

            try:
                # Load CSV with robust string → numeric conversion
                data = pd.read_csv(
                    csv_path,
                    dtype=str,
                    na_values=['', 'nan', 'NaN', '"', ' '],
                    engine='python',
                    on_bad_lines='skip',
                    comment='#'
                )

                if data.empty or len(data) < 10:
                    print(f"   → Empty or too small file ({len(data)} rows) → skipping")
                    continue

                data = data.dropna(axis=1, how='all')

                # Extract spectral columns
                signals = []
                for column in CONFIG['nist']['columns']:
                    if column not in data.columns:
                        continue
                    col = data[column].astype(str).str.strip('="\' ')
                    col = pd.to_numeric(col, errors='coerce')

                    if col.notna().sum() < 10:
                        continue

                    # Robust trimming (2.5%–97.5%)
                    q_low, q_high = col.quantile([0.025, 0.975])
                    signal = col[(col >= q_low) & (col <= q_high)].dropna().values

                    if len(signal) >= 10:
                        signals.append(signal)

                if len(signals) == 0:
                    print(f"   → No valid numeric signals found for {elements_str}")
                    continue

                # Equalize lengths for MFDFA
                for_mfdfa = signals[0]
                min_len = min(len(s) for s in signals)
                signals = [s[:min_len] for s in signals]

                # TCI/MCI extraction
                tci_extracted, mci_extracted = extract_tci_mci_data(signals, "nist")

                if for_mfdfa is None or len(for_mfdfa) < 50:
                    print(f"   → Insufficient data after preprocessing "
                          f"({len(for_mfdfa) if for_mfdfa is not None else 0} points)")
                    continue

                n_total = len(for_mfdfa)
                subset_size = max(50, n_total // 5)
                scales = CONFIG['nist']['scales'](n_total)
                n_subsets = min(100, max(10, n_total // subset_size))

                print(f"   → {n_total:,} rows → {n_subsets} subsets of ~{subset_size} each (dynamic)")

                results, _ = process_dataset(
                    data=for_mfdfa,
                    data_type="nist",
                    dataset_name=dataset_name,
                    scales=scales,
                    expected_D_f=CONFIG['nist']['expected_D_f'],
                    sigma_D_f=CONFIG['nist']['sigma_D_f']
                )

                if results is not None:
                    results_all[dataset_name] = results

                    if tci_extracted is not None and len(tci_extracted) >= 500:
                        tci_datasets[dataset_name] = tci_extracted
                    if mci_extracted is not None:
                        mci_datasets[dataset_name] = mci_extracted

                    metadata_log["datasets_loaded"][dataset_name] = {
                        "n_subsets": n_subsets,
                        "lines": len(for_mfdfa),
                        "file": os.path.basename(csv_path),
                        "status": "OK"
                    }
                    metadata_log["load_meta"][dataset_name] = {
                        "n_rows": len(for_mfdfa),
                        "tci_length": len(tci_extracted) if tci_extracted is not None else None,
                        "mci_shape": (
                            list(mci_extracted.shape)
                            if isinstance(mci_extracted, np.ndarray) else None
                        )
                    }

                    print(f"   → SUCCESS: D_f = {results['mean_D_f']:.3f} ± {results['std_D_f']:.3f}")
                else:
                    print(f"   → No valid MFDFA results for {elements_str}")

            except Exception as e:
                print(f"   → Error processing {elements_str}: {e}")

            finally:
                # Heavy cleanup per element
                for var in ["data", "signals", "for_mfdfa", "tci_extracted", "mci_extracted", "results"]:
                    if var in locals():
                        try:
                            del locals()[var]
                        except Exception:
                            pass
                gc.collect()
                gc.collect()
                print(f"   → Memory cleared for {elements_str}")

        print(f"\n[UTMF] Finished all {n_elements} NIST elements")

    # ================================================================
    #                N A N O G R A V   P U L S A R S
    #      (Stable, memory-safe loop with progress reporting)
    # ================================================================
    if CONFIG['nanograv'].get('utmf_use', False):
        pulsars = CONFIG['nanograv']['pulsar_list_utmf']
        n_pulsars = len(pulsars)

        print(f"\n[UTMF] Processing {n_pulsars} NANOGrav pulsars.")

        for p_idx, pulsar_name in enumerate(pulsars):
            dataset_name = f"NANOGrav_{pulsar_name}"
            print(f"\n→ [{p_idx+1}/{n_pulsars}] Loading & analyzing pulsar: "
                  f"{pulsar_name} → {dataset_name}")

            for_mfdfa = None
            tci_extracted = None
            mci_extracted = None
            results = None

            try:
                # Load pulsar time-series (residuals + additional channels)
                for_mfdfa, tci_extracted, mci_extracted = load_data(
                    file_path=None,
                    data_type="nanograv",
                    pulsar_name=pulsar_name
                )

                if for_mfdfa is None or len(for_mfdfa) < 1024:
                    print(f"   → Skipped {pulsar_name}: insufficient or corrupt data")
                    metadata_log["errors"].append(
                        {"dataset": dataset_name, "error": "insufficient data"}
                    )
                    continue

                print(f"   → Residual length = {len(for_mfdfa):,} points → "
                      f"{CONFIG['nanograv']['n_subsets']} subsets of "
                      f"{CONFIG['nanograv']['subset_size']}")

                # MFDFA computation
                results, _ = process_dataset(
                    data=for_mfdfa,
                    data_type="nanograv",
                    dataset_name=dataset_name,
                    scales=CONFIG['nanograv']['scales'],
                    expected_D_f=CONFIG['nanograv']['expected_D_f'],
                    sigma_D_f=CONFIG['nanograv']['sigma_D_f']
                )

                if results is None:
                    print(f"   → No valid MFDFA results for {pulsar_name}")
                    continue

                results_all[dataset_name] = results
                mci_datasets[dataset_name] = mci_extracted

                if tci_extracted is not None and len(tci_extracted) >= 1024:
                    tci_datasets[dataset_name] = tci_extracted
                    print(f"   → TCI added (length = {len(tci_extracted):,})")
                else:
                    print(f"   → No valid TCI for {pulsar_name}")

                metadata_log["datasets_loaded"][dataset_name] = {
                    "n_subsets": CONFIG["nanograv"]["n_subsets"],
                    "residuals_length": len(for_mfdfa),
                    "status": "OK"
                }
                metadata_log["load_meta"][dataset_name] = {
                    "n_rows": len(for_mfdfa),
                    "tci_length": len(tci_extracted) if tci_extracted is not None else None,
                    "mci_shape": (
                        list(mci_extracted.shape)
                        if isinstance(mci_extracted, np.ndarray) else None
                    )
                }

            except Exception as e:
                print(f"   → Error in pulsar {pulsar_name}: {e}")
                metadata_log["errors"].append(
                    {"dataset": dataset_name, "error": str(e)}
                )

            finally:
                # Memory cleanup
                for var in ["for_mfdfa", "tci_extracted", "mci_extracted", "results"]:
                    if var in locals():
                        try:
                            del locals()[var]
                        except Exception:
                            pass
                gc.collect()
                gc.collect()
                print(f"   → Memory cleared for {pulsar_name}")

        print(f"\n[UTMF] Finished all {n_pulsars} NANOGrav pulsars")

    # ================================================================
    #                        Q R N G  (Quantum RNG)
    # ================================================================
    if CONFIG['qrng'].get('utmf_use', False):
        dataset_name = "NIST_QRNG"
        print(f"\n[UTMF] Loading QRNG dataset: {dataset_name}")

        for_mfdfa, tci_extracted, mci_extracted = load_data(None, "qrng")

        if for_mfdfa is not None:
            results, _ = process_dataset(
                for_mfdfa,
                "qrng",
                dataset_name,
                CONFIG["qrng"]["scales"],
                CONFIG["qrng"]["expected_D_f"],
                CONFIG["qrng"]["sigma_D_f"]
            )

            if results:
                results_all[dataset_name] = results
                tci_datasets[dataset_name] = tci_extracted
                mci_datasets[dataset_name] = mci_extracted

                if tci_extracted is not None and len(tci_extracted) >= 12288:
                    print(f"   → TCI {dataset_name}: using longest valid column "
                          f"({len(tci_extracted):,} samples, quantum bitstream)")
                else:
                    print(f"   → TCI {dataset_name}: no suitable column → excluded from TCI")

                metadata_log["datasets_loaded"][dataset_name] = {
                    "n_subsets": CONFIG["qrng"]["n_subsets"],
                    "status": "OK"
                }
                metadata_log["load_meta"][dataset_name] = {
                    "n_rows": len(for_mfdfa),
                    "tci_length": len(tci_extracted) if tci_extracted is not None else None,
                    "mci_shape": (
                        list(mci_extracted.shape)
                        if isinstance(mci_extracted, np.ndarray) else None
                    )
                }

    gc.collect()

    # ================================================================
    #                     G A I A   D R 3  (Astrometry)
    # ================================================================
    if CONFIG['gaia'].get('utmf_use', False):
        dataset_name = CONFIG["gaia"]["name"]
        print(f"\n[UTMF] Loading Gaia dataset: {dataset_name}")

        for_mfdfa, tci_extracted, mci_extracted = load_data(
            CONFIG["gaia"]["file"], "gaia"
        )

        if for_mfdfa is not None:
            results, _ = process_dataset(
                for_mfdfa,
                "gaia",
                dataset_name,
                CONFIG["gaia"]["scales"],
                CONFIG["gaia"]["expected_D_f"],
                CONFIG["gaia"]["sigma_D_f"]
            )

            if results:
                results_all[dataset_name] = results
                tci_datasets[dataset_name] = tci_extracted
                mci_datasets[dataset_name] = mci_extracted

                if tci_extracted is not None and len(tci_extracted) >= 12288:
                    print(f"   → TCI {dataset_name}: using longest valid column "
                          f"({len(tci_extracted):,} samples)")
                else:
                    print(f"   → TCI {dataset_name}: no suitable column → excluded from TCI")

                metadata_log["datasets_loaded"][dataset_name] = {
                    "n_subsets": CONFIG["gaia"]["n_subsets"],
                    "status": "OK"
                }
                metadata_log["load_meta"][dataset_name] = {
                    "n_rows": len(for_mfdfa),
                    "tci_length": len(tci_extracted) if tci_extracted is not None else None,
                    "mci_shape": (
                        list(mci_extracted.shape)
                        if isinstance(mci_extracted, np.ndarray) else None
                    )
                }

    gc.collect()

    # ================================================================
    #                           F I N I S H E D
    # ================================================================
    print(f"[UTMF] Loaded {len(results_all)} datasets.")
    return results_all, tci_datasets, mci_datasets

def build_dataset_health(metadata_log):
    """
    Construct a reconstruction-ready dataset health summary for each dataset,
    based on 'subsets_processed', 'subset_warnings', and 'errors'.

    This enables full run reconstruction from the JSON snapshot, without
    needing access to the original raw data.
    """
    out = {}

    subsets = metadata_log.get("subsets_processed", {})
    warnings = metadata_log.get("subset_warnings", {})
    errors = metadata_log.get("errors", [])

    # Build quick error lookup per dataset
    error_map = {}
    for e in errors:
        ds = e.get("dataset")
        if ds is None:
            continue
        error_map.setdefault(ds, []).append(e.get("error"))

    for dataset_name, proc in subsets.items():

        Df_vals = np.array(proc.get("D_f_values", []), dtype=float)
        hq_vals = proc.get("hq_values", [])

        # h(q) spectrum summary
        if hq_vals and len(hq_vals) > 0:
            try:
                hq_stack = np.vstack([
                    np.asarray(h) for h in hq_vals if len(h) > 0
                ])
                hq_mean = np.nanmean(hq_stack, axis=0)
                hq_std  = np.nanstd(hq_stack, axis=0)
                hq_min  = np.nanmin(hq_stack, axis=0)
                hq_max  = np.nanmax(hq_stack, axis=0)
                # Multifractal spectrum width Δh
                hq_width = float(np.nanmax(hq_mean) - np.nanmin(hq_mean))
            except Exception:
                hq_mean = []
                hq_std  = []
                hq_min  = []
                hq_max  = []
                hq_width = np.nan
        else:
            hq_mean = []
            hq_std  = []
            hq_min  = []
            hq_max  = []
            hq_width = np.nan

        # D_f statistics
        if len(Df_vals) > 0:
            Df_mean   = float(np.nanmean(Df_vals))
            Df_std    = float(np.nanstd(Df_vals))
            nan_ratio = float(np.mean(~np.isfinite(Df_vals)))
        else:
            Df_mean   = np.nan
            Df_std    = np.nan
            nan_ratio = 1.0

        # Subset counts
        n_total = int(proc.get("n_subsets", len(Df_vals)))
        n_valid = int(proc.get("n_valid_subsets", np.sum(np.isfinite(Df_vals))))

        out[dataset_name] = {
            "summary": {
                "n_subsets_total": n_total,
                "n_subsets_valid": n_valid,
                "nan_ratio_Df": nan_ratio,
                "Df_mean": Df_mean,
                "Df_std": Df_std,
                "hq_mean": convert_numpy(hq_mean),
                "hq_std": convert_numpy(hq_std),
                "hq_min": convert_numpy(hq_min),
                "hq_max": convert_numpy(hq_max),
                "hq_width": hq_width,
                "p_value": proc.get("p_value", np.nan),
            },

            "subset_size_used": proc.get("subset_size", None),
            "scales_used": proc.get("scales_used", None),

            "warnings": warnings.get(dataset_name, []),
            "errors": error_map.get(dataset_name, []),
        }

    return out

def build_tci_timeseries_snapshot(tci_datasets, max_len=32768):
    """
    Build a compact snapshot of all TCI time series for JSON export.

    For each dataset:
        - Store at most `max_len` samples starting at index 0
        - Record whether truncation occurred
        - Accept both 1D and 2D inputs (flattened to 1D)

    The resulting structure is suitable for direct inclusion in the
    FULL_DETAILS JSON; convert_numpy() can be applied later to ensure
    full JSON-serialisability.
    """
    snapshot = {}

    for name, arr in tci_datasets.items():
        if arr is None:
            continue

        ts = np.asarray(arr)

        # 1D/2D/whatever → 1D representation
        ts_1d = ts.ravel()
        full_len = int(len(ts_1d))

        if full_len == 0:
            continue

        if full_len > max_len:
            used = ts_1d[:max_len]
            truncated = True
        else:
            used = ts_1d
            truncated = False

        snapshot[name] = {
            "length_full": full_len,
            "length_used": int(len(used)),
            "truncated": bool(truncated),
            # convert_numpy can later turn this into a plain list
            "data": used.astype(float)
        }

    return snapshot

# ============================================================
# Unified TCI/MCI/TMCI computation (clean + structured + JSON-ready)
# ============================================================

def compute_all_indices(tci_datasets, results_all, config, metadata_log):
    """
    Full unified computation pipeline for coherence indices.

    This function computes:
        • TCI  – Temporal Coherence Index
        • MCI  – Measurement Coherence Index
        • TMCI – Combined temporal–measurement coherence index
        • Bootstrap samples for (TCI, MCI, TMCI)
        • Cross-validation folds for TMCI (k = 3, if enabled)
        • A JSON-friendly nested structure of all outputs

    Parameters
    ----------
    tci_datasets : dict
        Dataset → time series used for TCI.
    results_all : dict
        Dataset → MFDFA result dictionaries (with 'hq_values' for MCI).
    config : dict
        UTMF configuration dictionary (includes 'cross_val' flag).
    metadata_log : dict
        Global metadata log to be extended with index-related entries.

    Returns
    -------
    tmci_mean : float
        Mean TMCI value from the bootstrap ensemble.
    tmci_ci : tuple(float, float)
        95% confidence interval for TMCI.
    tmci_corr : float
        Correlation between bootstrapped TCI and MCI.
    tci_samples : ndarray
        Bootstrap samples for TCI.
    mci_samples : ndarray
        Bootstrap samples for MCI.
    tmci_samples : ndarray
        Bootstrap samples for TMCI.
    """
    print("\n[UTMF] Computing indices (TCI, MCI, TMCI)...")

    # -----------------------------------------------------------
    # 1. TCI
    # -----------------------------------------------------------
    tci, tci_pairs = calculate_tci_multivariate(tci_datasets)
    print(f"TCI = {tci:.3f}")

    metadata_log["tci_pairs"] = tci_pairs

    # -----------------------------------------------------------
    # 2. MCI
    # -----------------------------------------------------------
    mci = calculate_mci_multivariate(results_all)
    print(f"MCI = {mci:.3f}")

    # -----------------------------------------------------------
    # 3. TMCI + Bootstrap
    # -----------------------------------------------------------
    (
        tmci_mean, tmci_ci, tmci_corr,
        tci_samples, mci_samples, tmci_samples
    ) = calculate_tmci_empirical(tci, mci)

    print(f"TMCI = {tmci_mean:.3f}   (95% CI = [{tmci_ci[0]:.3f}, {tmci_ci[1]:.3f}])")
    print(f"TCI–MCI corr = {tmci_corr:.3f}")

    # -----------------------------------------------------------
    # 4. Cross-validation (k = 3)
    # -----------------------------------------------------------
    crossval_folds = []
    crossval_ci = []

    if config.get("cross_val", False):
        print("\n[UTMF] Running TMCI cross-validation (k=3)...")

        for fold in range(3):
            fold_results = {}

            # Thinned h(q) fields per dataset
            for name, res in results_all.items():
                if not isinstance(res, dict) or "hq_values" not in res:
                    continue

                hq = res["hq_values"]
                if len(hq) < 3:
                    continue

                # Randomly select ~2/3 of the subsets for this fold
                idx = np.random.choice(len(hq), size=int(len(hq) * 0.67), replace=False)
                folded_hq = np.stack([hq[i] for i in idx])
                folded_mean = np.nanmean(folded_hq, axis=0)

                fold_results[name] = {"hq_values": [folded_mean]}

            # Restrict TCI datasets to those present in this fold
            fold_tci_ds = {
                name: tci_datasets[name]
                for name in fold_results.keys()
                if name in tci_datasets
            }

            fold_tci = calculate_tci_multivariate(fold_tci_ds)[0] if fold_tci_ds else np.nan
            fold_mci = calculate_mci_multivariate(fold_results)

            (
                tmci_f_mean, tmci_f_ci, tmci_f_corr,
                _, _, _
            ) = calculate_tmci_empirical(fold_tci, fold_mci)

            crossval_folds.append(tmci_f_mean)
            crossval_ci.append(tmci_f_ci)

            print(
                f"  Fold {fold+1}: TMCI = {tmci_f_mean:.3f}, "
                f"CI = [{tmci_f_ci[0]:.3f}, {tmci_f_ci[1]:.3f}]"
            )

    crossval_std = float(np.nanstd(crossval_folds)) if crossval_folds else np.nan

    if crossval_folds:
        print(f"TMCI cross-val std = {crossval_std:.3f}")
    else:
        print("No valid cross-val folds.")

    # -----------------------------------------------------------
    # 5. STRUCTURED JSON OUTPUT
    # -----------------------------------------------------------
    indices = {
        "tci": float(tci),
        "mci": float(mci),
        "tmci": {
            "mean": float(tmci_mean),
            "ci": [float(tmci_ci[0]), float(tmci_ci[1])],
            "corr": float(tmci_corr),
            "bootstrap": {
                "tci": convert_numpy(tci_samples),
                "mci": convert_numpy(mci_samples),
                "tmci": convert_numpy(tmci_samples)
            },
            "crossval": {
                "folds": convert_numpy(crossval_folds),
                "fold_ci": convert_numpy(crossval_ci),
                "std": float(crossval_std)
            }
        }
    }

    # -----------------------------------------------------------
    # 6. Insert into metadata_log (block + backward compatibility)
    # -----------------------------------------------------------
    metadata_log["indices"] = indices

    # Backwards-compatible top-level fields (for CSV and legacy code)
    metadata_log["tci"] = indices["tci"]
    metadata_log["mci"] = indices["mci"]
    metadata_log["tmci"] = indices["tmci"]["mean"]
    metadata_log["tmci_std_cv"] = indices["tmci"]["crossval"]["std"]
    metadata_log["tmci_folds"] = indices["tmci"]["crossval"]["folds"]

    return indices

# ============================================================
# FINAL EXECUTION ORDER FOR UTMF v5
# ============================================================

# 1. Load all datasets
results_all, tci_datasets, mci_datasets = load_utmf_data()

# 2. Compute full TCI / MCI / TMCI (with bootstrap and cross-validation)
indices = compute_all_indices(tci_datasets, results_all, CONFIG, metadata_log)

# 2b. Add compact TCI time-series snapshots (≤ 32k samples) to metadata_log
metadata_log["tci_timeseries_snapshot"] = build_tci_timeseries_snapshot(
    tci_datasets,
    max_len=32758
)

# 3. Dataset health summary for reconstruction (no raw data required)
metadata_log["dataset_health"] = build_dataset_health(metadata_log)

print("\n[UTMF] Final indices (compact):")
print(f"  TCI  = {indices['tci']:.3f}")
print(f"  MCI  = {indices['mci']:.3f}")
print(
    f"  TMCI = {indices['tmci']['mean']:.3f} "
    f"(95% CI = {indices['tmci']['ci'][0]:.3f} – {indices['tmci']['ci'][1]:.3f})"
)
print(f"  Corr(TCI, MCI) = {indices['tmci']['corr']:.3f}")
print(f"  Cross-validation std = {indices['tmci']['crossval']['std']:.3f}")

# 4. Write metadata to disk (CSV + JSON snapshots)
save_run_metadata(CONFIG['metadata']['save_flag'])

# ============================================================
# Full JSON export (optional)
# ============================================================
if SAVE_FULL_DETAILS_JSON:
    json_path = f"{OUTPUT_DIR}FULL_DETAILS_utmf_v5.1_{timestamp}.json"

    safe_log = copy.deepcopy(metadata_log)

    # Remove lambda objects
    safe_config = {k: v for k, v in CONFIG.items() if not callable(v)}
    if 'scales' in safe_config.get('nist', {}):
        safe_config['nist']['scales'] = "lambda removed for JSON"
    if 'subset_size' in safe_config.get('nist', {}):
        safe_config['nist']['subset_size'] = "lambda removed for JSON"

    safe_log['config_snapshot'] = safe_config
    safe_log = convert_numpy(safe_log)

    # ============================================================
    # Embed full MFDFA implementation (jit-free) into JSON
    # ============================================================
    import inspect
    import hashlib

    def strip_jit_decorators(source):
        """Remove lines starting with @jit(...) for portability."""
        clean = []
        for line in source.split("\n"):
            if line.strip().startswith("@jit"):
                continue
            clean.append(line)
        return "\n".join(clean)

    def get_clean_function_source(func):
        """Return source code with @jit decorators stripped."""
        raw = inspect.getsource(func)
        return strip_jit_decorators(raw)

    # Collect MF-DFA helpers + core routine
    src_polyfit  = get_clean_function_source(polyfit_linear)
    src_polyval  = get_clean_function_source(polyval_linear)
    src_mfdfa    = get_clean_function_source(jedi_mfdfa)

    mfdfa_block = (
        src_polyfit.strip()  + "\n\n" +
        src_polyval.strip()  + "\n\n" +
        src_mfdfa.strip()    + "\n"
    )

    # Hash for integrity checking
    mf_sha = hashlib.sha256(mfdfa_block.encode("utf-8")).hexdigest()

    safe_log["mfdfa_code"] = {
        "version": "v5.2",
        "sha256": mf_sha,
        "block": mfdfa_block
    }

    print(f"[UTMF] Embedded jedi_mfdfa block into JSON (SHA256={mf_sha[:16]}...)")

    with open(json_path, "w") as f:
        json.dump(safe_log, f, indent=2, default=str)

    size_mb = os.path.getsize(json_path) / 1e6
    print(f"→ Full JSON saved ({size_mb:.1f} MB): {json_path}")
else:
    print("→ Full JSON skipped (SAVE_FULL_DETAILS_JSON = False)")

# ============================================================
# Export globals for interactive Cell 2 analysis
# ============================================================
globals().update({
    'utmf_results_all_full': results_all,
    'utmf_tci_datasets_full': tci_datasets,
    'utmf_mci_datasets_full': mci_datasets,
    'utmf_config_snapshot': CONFIG.copy(),
    'utmf_metadata_log': metadata_log,
    'utmf_tci_pairs_full': metadata_log["tci_pairs"],
    'utmf_mci_measurements_full': metadata_log.get("mci_measurements", []),
    'utmf_mci_full': indices["mci"],
    'utmf_tci_full': indices["tci"],
    'utmf_tmci_full': indices["tmci"]["mean"],
    'utmf_tmci_std_cv_full': indices["tmci"]["crossval"]["std"],
    'utmf_tmci_folds_full': indices["tmci"]["crossval"]["folds"],
})

print("Globals updated — ready for interactive analysis (Cell 2).")
print(f"TMCI = {indices['tmci']['mean']:.3f}, MCI = {indices['mci']:.3f}")
