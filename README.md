# README.md

````markdown
# UTMF v5.1 — Universal Temporal–Measurement Framework  
## Cross-Domain Multifractal Coherence Test (TCI/MCI) — Fully Reproducible Release

This repository contains the complete, self-contained reproducible package used in:

> **M. Eversdijk (2025).**  
> *Multifractal Coherence in Physical Measurements:  
> A Fully Reproducible Cross-Domain Analysis Using UTMF v5.1.*  
> DOI: https://doi.org/10.5281/zenodo.17927721

The aim of this release is to provide a **transparent, inspectable computational archive** of the UTMF v5.1 multifractal coherence analysis across gravitational-wave interferometry, CMB maps, large-scale structure surveys, atomic spectra, collider data, pulsar timing residuals, quantum random number generators, and astrometric catalogues.

---

## 1. What is UTMF v5.1?

UTMF v5.1 is a **domain-agnostic framework** for extracting, comparing, and validating multifractal structure in arbitrary measurement systems. It focuses on two coherence measures:

### 1.1 Temporal Coherence Index (TCI)

TCI measures the **stability of multifractal exponents over time** within a single dataset.

- Computed from fixed-length time series (here: 6144 samples).  
- Uses pairwise absolute Pearson correlations between reconstructed \(h(q)\) curves:  
  \[
  \mathrm{TCI}_{ij} = \bigl|\,\rho(h_i(q), h_j(q)) \,\bigr|.
  \]

### 1.2 Measurement Coherence Index (MCI)

MCI measures **cross-dataset similarity** of mean multifractal spectra:

- For each dataset \(k\), UTMF computes subset-level \(h_{k,n}(q)\).  
- The dataset-level mean spectrum \(\bar{h}_k(q)\) is formed by averaging over subsets.  
- MCI is the mean absolute correlation between these mean spectra across datasets.

Both indices are derived from the **same MFDFA implementation** and the same \(q\)-grid and scale configuration embedded directly in the JSON archive.

---
UTMF-v5.1-Coherence/
│
├── README.md
│   ← Overview of the project, results, and reproducibility
│
├── LICENSE
│   ← MIT License
│
├── utmf_v5.1.py
│   ← Complete UTMF v5.1 pipeline as a single executable script
│     (exactly as used in the paper; provided as-is)
│
├── FULL_DETAILS_utmf_v5.1_*.json
│   ← Self-contained computational archive:
│     • embedded MFDFA implementation
│     • configuration snapshot
│     • stored time series
│     • all h(q) spectra
│     • all metadata required for full reconstruction
│
├── results/
│   ├── README.md
│   │   ← Explanation of all figures and CSV files used in the paper
│   │
│   ├── UTMF_TCI_MCI_combined_nulltest_*.csv
│   │   ← Primary null-test output (real vs synthetic distributions)
│   │
│   ├── UTMF_VALIDATION_metrics_*.csv
│   │   ← Scalar statistical validation metrics (13-test suite)
│   │
│   ├── UTMF_VALIDATION_summary_*.png
│   │   ← Multi-panel validation figure (paper-grade)
│   │
│   └── figures/
│       ├── fig_tci_hist_*.png
│       ├── fig_mci_hist_*.png
│       ├── fig_tci_cdf_*.png
│       ├── fig_mci_cdf_*.png
│       ├── fig_real_distance_kde_*.png
│       ├── fig_real_dendrogram_*.png
│       ├── fig_real_tci_heatmap_*.png
│       ├── fig_real_pca_hq_*.png
│       └── fig_mean_hq_overlay_*.png
│
└── DATASET_LINKS.md
    ← Direct download links and references to all datasets used


> **Note:** Only `UTMF_FULL_DETAILS.json` is required to reproduce all figures and statistics from the paper.
> The original `utmf_v5.1.py` is included as a **single cell pipeline** for transparency; users are free to refactor it into modules.

---

## 3. The `UTMF_FULL_DETAILS.json` Cell

The file `UTMF_FULL_DETAILS.json` (~85 MB) is a **single, self-contained computational cell** produced by UTMF v5.1. It contains:

* All time series used for TCI (truncated to a maximum of 32,768 samples per dataset).
* All subset-level MFDFA outputs used for MCI.
* The complete `jedi_mfdfa` implementation (detrending order 0) and linear helper routines.
* The scale grid and (q)-grid used for all analyses.
* All per-dataset multifractal spectra (h(q)) and derived indices.
* Configuration snapshot, random seeds, and internal diagnostics.

Because the JSON embeds the **exact MFDFA algorithm and configuration**, no separate UTMF installation is required to reproduce the paper.  All key results are reconstructed directly from this file.

---

## 4. Reproducibility Workflow

### 4.1 Environment

The analysis was originally developed and executed in a Python environment with:

* Python ≥ 3.10
* NumPy, SciPy, pandas
* scikit-learn
* Matplotlib, seaborn

A minimal environment file (e.g. `environment.yml` or `requirements.txt`) can be added if desired, but the notebooks and scripts are standard scientific Python.

---

### 4.2 Step 1 — Reconstruct multifractal spectra from JSON

All multifractal spectra used in the paper are either stored directly in, or recomputed from, the JSON archive. A minimal example:

```python
import json
import numpy as np

with open("UTMF_FULL_DETAILS.json", "r") as f:
    meta = json.load(f)

# Example: access the config snapshot and h(q) values
config = meta["config_snapshot"]
subsets = meta["subsets_processed"]

# h(q) curves for a single dataset
dset_name = next(iter(subsets.keys()))
info = subsets[dset_name]
hq_values = np.array(info["hq_values"], dtype=float)  # shape: (n_subsets, n_q)
```

In the notebooks, this reconstruction is wrapped in higher-level utilities that:

* repair the (q=0) point if necessary,
* interpolate across any remaining NaNs,
* enforce a common (q)-grid across datasets.

---

### 4.3 Step 2 — Run the TCI/MCI null-test

Use:

* `notebooks/utmf_nulltest.ipynb` (interactive), or
* `scripts/run_nulltest.py` (non-interactive CLI)

This step:

1. Reads `UTMF_FULL_DETAILS.json`.
2. Extracts fixed-length time series for TCI (length = **6144 samples**).
3. Recomputes (h(q)) for each real series using the embedded `jedi_mfdfa` block.
4. Generates synthetic ensembles (white, AR(1), pink, fGn-like, phase-surrogate).
5. Computes TCI (pairwise (|\mathrm{corr}(h_i, h_j)|)) and MCI (mean (|\mathrm{corr}(\bar{h}_i, \bar{h}_j)|)).
6. Saves a combined CSV:

```text
UTMF_TCI_MCI_combined_nulltest_YYYYMMDD_HHMMSS_UTC.csv
```

The CSV includes real and synthetic samples for both indices, along with a metadata row:

```text
Index,Type,Value,Dataset
...
META,TCI_length,6144,GLOBAL
```

This `TCI_length` entry is used by the validation suite to keep track of the fixed TCI time-series length.

---

### 4.4 Step 3 — Run the 13-test validation suite

Use:

* `notebooks/utmf_validation.ipynb`, or
* `scripts/run_validation.py`

These read the most recent combined null-test CSV and compute a **canonical 13-test battery** comparing real vs synthetic distributions for:

* TCI (pairwise distribution)
* MCI (per-dataset scalars)
* Combined (TCI + MCI concatenated)

The suite includes:

* Cohen’s d (effect size)
* ROC AUC
* Cross-validated logistic regression accuracy
* Permutation p-value (mean difference)
* Approximate Bayes Factor (log10 BF10)
* Silhouette score
* Energy distance
* Kolmogorov–Smirnov p-value
* Mann–Whitney U p-value
* Cramér–von Mises p-value

Outputs:

* `UTMF_VALIDATION_metrics_YYYYMMDD_HHMMSS_UTC.csv`
  (table of scalar metrics: Index, Metric, Value)

* `UTMF_VALIDATION_summary_YYYYMMDD_HHMMSS_UTC.png`
  (multi-panel figure summarising distributions, CDFs, effect sizes, and p-values)

---

## 5. Dataset Download Information

The original UTMF v5.1 run used publicly available datasets from a range of experiments and surveys. For licensing and size reasons, the raw data are **not** bundled in this repository.

## **All datasets used for UTMF v5.1 run listed below.** (Total 11.2GB) 
#### - Direct downloadlinks.
#### - For Colab: Create: /MyDrive/Datasets_UTMF/UTMF_outputs/
#### - Place the datasets in folder: /Datasets_UTMF/
#### - Mount Drive
#### - Run UTMF v5.1
#### - Results are returned in folder: /UTMF_outputs/
-----
- **[LIGO – GWOSC](https://gwosc.org/archive/links/O4a_16KHZ_R1/L1/1368195220/1389456018/simple/)**  
  HDF5 strain files (e.g., `L-L1_GWOSC_O4a_16KHZ_R1-*.hdf5`).
                                                           
  **Datasets used in UTMF v5.1 configuration:**
- `L-L1_GWOSC_O4a_16KHZ_R1-1384779776-4096.hdf5` [Download](https://gwosc.org/archive/data/O4a_16KHZ_R1/1384120320/L-L1_GWOSC_O4a_16KHZ_R1-1384779776-4096.hdf5) (486MB)
- `L-L1_GWOSC_O4a_16KHZ_R1-1368350720-4096.hdf5` [Download](https://gwosc.org/archive/data/O4a_16KHZ_R1/1367343104/L-L1_GWOSC_O4a_16KHZ_R1-1368350720-4096.hdf5) (486MB)
- `L-L1_GWOSC_O4a_16KHZ_R1-1370202112-4096.hdf5` [Download](https://gwosc.org/archive/data/O4a_16KHZ_R1/1369440256/L-L1_GWOSC_O4a_16KHZ_R1-1370202112-4096.hdf5) (486MB)
- `L-L1_GWOSC_O4a_16KHZ_R1-1389420544-4096.hdf5` [Download](https://gwosc.org/archive/data/O4a_16KHZ_R1/1389363200/L-L1_GWOSC_O4a_16KHZ_R1-1389420544-4096.hdf5) (486MB)
---
- **[Planck – ESA Archive](https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/)**  
  FITS CMB maps (e.g., SMICA IQU maps such as `COM_CMB_IQU-smica_2048_R3.00_full.fits`).

  **Datasets used in UTMF v5.1 configuration:**
- `COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits` [Download](https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/cmb/COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits) (384MB)
- `COM_CMB_IQU-smica_2048_R3.00_full.fits`      [Download](https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/cmb/COM_CMB_IQU-smica_2048_R3.00_full.fits) (1.88GB)
- `LFI_SkyMap_070_1024_R3.00_survey-1.fits`     [Download](https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/LFI_SkyMap_070_1024_R3.00_full.fits) (480MB)
---
- **[DESI – Data Release Portal](https://data.desi.lbl.gov/doc/releases/dr1/)**  
  LRG FITS catalogs (e.g., `LRG_full.dat.fits`).

  **Dataset used in UTMF v5.1 configuration:**
- `LRG_full.dat.fits` [Download](https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.2/LRG_full.dat.fits) (2.77GB)
---  
- **[CERN Open Data](https://opendata.cern.ch/record/15007)**  
  ROOT event files (e.g., `data_B.exactly2lep.root`).

  **Dataset used in UTMF v5.1 configuration:**
- `data_B.exactly2lep.root` [Download:](https://opendata.cern.ch/record/15007/files/data_B.exactly2lep.root) (451MB)
                                                                                              
    ➕ This repository includes a helpfile: `data_B.exactly2lep.h5` [Download:](https://github.com/Jedi-Markus-Strive/UTMF-CRISP/raw/refs/heads/main/Datasets/data_B.exactly2lep.h5) (315KB) (helpfile for .root, store it at the same location as .root-file.)
---
- **[NIST Atomic Spectra Database](https://physics.nist.gov/PhysRefData/ASD/lines_form.html)**                          
  CSV spectra 
  
  **Dataset used in UTMF v5.1 configuration:**                                                                        
- `NIST_elements` [Download](https://github.com/Jedi-Markus-Strive/UTMF-CRISP/raw/refs/heads/main/Datasets/NIST_3.zip)** (3.3MB) (Complete dataset as used in UTMF v5.0, unzip and use the CSV
(17.4MB) for UTMF analysis.)
---
- **[NANOGrav Data Releases](https://zenodo.org/records/16051178)**  
- Pulsar timing residuals (e.g., `NG15yr narrowband` files).

  **Dataset used in UTMF v5.1 configuration:**
- `NANOGrav15yr_PulsarTiming_v2.1.0` [Download:](https://zenodo.org/records/16051178/files/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz?download=1) (639MB) ([Unzip](https://github.com/Jedi-Markus-Strive/UTMF-CRISP/raw/refs/heads/main/prepare_NANOGrav_15yr_data.ipynb) the file, use for UTMF analysis.)
---
- **[Gaia Archive (DR3)](https://vizier.cds.unistra.fr/viz-bin/VizieR-4)**  
  Source catalogs in TSV format (e.g., `gaia_dr3.tsv`).                                                                 
  **Dataset used in UTMF v5.1:**                                                                                      
      **Select:**                                                                                                       
        1- 'gaiadr3'                                                                                                   
        2. Table: `I/355/gaiadr3`  
        3. Rows: `1-999999`  
        4. Format: Tab-Separated Values  
        5. Columns: All                                                                                
        6. Rename file to 'gaia_dr3' or update path in config.                                                          
---
- **[ANU Quantum Random Numbers (QRNG)](https://qrng.anu.edu.au/)**  
  API-based quantum random sequences (no download required, incorporated in UTMF v5.1 configuration).
---

                                                                                            
Users who wish to **re-run the full UTMF pipeline** (not required for reproducing this paper’s figures) can download the raw data following that manifest and use `utmf_v5.1.py` as a reference pipeline.

---

## 6. Citation

If you use this repository, the JSON archive, or the coherence indices in your own work, please cite both the paper and the software release. Once the Zenodo DOI is available, a citation entry will look like:

```bibtex
@misc{eversdijk2025_utmf51,
  author       = {Eversdijk, M.},
  title        = {UTMF v5.1: Cross-Domain Multifractal Coherence Test},
  year         = {2025},
  doi          = {10.xxxx/zenodo.xxxxx},
  publisher    = {Zenodo},
}
```

and the accompanying paper:

```bibtex
@article{eversdijk2025_crossdomain,
  author  = {Eversdijk, M.},
  title   = {Cross-Domain Multifractal Coherence in Physical Measurements: \\ A Fully Reproducible Test Using UTMF v5.1},
  journal = {TBA},
  year    = {2025},
  volume  = {TBA},
  pages   = {TBA},
  doi     = {TBA}
}
```

---

## 7. Philosophy of this Release

This repository is intended as a **reference implementation and a starting point**, not as a large software ecosystem.

* The central object is the `UTMF_FULL_DETAILS.json` cell, which acts as a **portable measurement archive**.
* The Python notebooks and scripts are deliberately straightforward and transparent.
* Researchers are encouraged to:

  * swap in their own datasets,
  * modify or extend the null models,
  * change the fixed TCI length,
  * or plug in alternative multifractal estimators,
    while keeping the coherence indices conceptually comparable.

The core idea is simple: **treat multifractal structure as a measurable, cross-domain coherence signal**, and make the entire pipeline **open, inspectable, and exactly reproducible**.

---

## 8. License

This project is released under the MIT License (see `LICENSE`).

````

---

# LICENSE

```text
MIT License

Copyright (c) 2025 M. Eversdijk

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
````

---

# utmf_v5.1.py

```python
"""UTMF v5.1 — unified pipeline cell

This file contains the original UTMF v5.1 analysis cell as a single Python
module. In the paper, this cell is executed once to generate the
`UTMF_FULL_DETAILS.json` archive, which then serves as the sole input for all
reconstruction, null tests, and statistical validation.

For transparency, the entire cell is kept here in its original, monolithic
form so that users can:

- inspect every processing step used to generate the JSON;
- re-run the full UTMF v5.1 pipeline if they have access to the same datasets;
- or refactor the code into modular components if desired.

