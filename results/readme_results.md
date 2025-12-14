# UTMF v5.1 — Results & Validation Outputs

This folder contains **all numerical results, figures, and validation artifacts** used in the paper:

> **"Cross-Domain Multifractal Coherence in Physical Measurements: A Fully Reproducible Test Using UTMF v5.1"**

The contents of this directory are generated **entirely from a single UTMF v5.1 run** and can be reproduced exactly from the accompanying `FULL_DETAILS` JSON file.

---

## 📦 What This Folder Contains

### 1. Null-Test CSV (Core Statistical Artifact)

**`UTMF_TCI_MCI_combined_nulltest_*.csv`**

This CSV is the *authoritative input* for all statistical validation and plots. It contains:

- Pairwise **Temporal Coherence Index (TCI)** values
- Dataset-level **Measurement Coherence Index (MCI)** values
- Corresponding **synthetic null-model values**
- A global metadata row specifying the fixed TCI time-series length

**Column structure:**
- `Index`   → `TCI`, `MCI`, or `META`
- `Type`    → `Real`, `Synthetic`, or `TCI_length`
- `Value`   → scalar coherence value or metadata
- `Dataset` → dataset identifier or `GLOBAL`

This CSV alone is sufficient to reproduce:
- All hypothesis tests
- All validation figures
- All summary statistics reported in the paper

---

### 2. Validation Metrics CSV

**`UTMF_VALIDATION_metrics_*.csv`**

Contains the scalar outputs of the **canonical 13-test validation suite**, for:
- TCI
- MCI
- Combined TCI+MCI distribution

Each row has the form:
- `Index`  → `TCI`, `MCI`, or `Combined`
- `Metric` → name of statistical diagnostic
- `Value`  → numerical result

Metrics include:
- Cohen’s *d*
- ROC AUC
- Cross-validated accuracy
- Permutation *p*-value
- Bayes factor (log-scale)
- Silhouette score
- Energy distance
- KS / Mann–Whitney / Cramér–von Mises *p*-values

---

### 3. Validation Figures (Paper-Grade)

All figures are generated automatically from the null-test CSV.

#### Core Validation Summary
- **`UTMF_VALIDATION_summary_*.png`**  
  Multi-panel figure combining:
  - TCI distribution
  - MCI distribution
  - Combined distribution
  - CDFs
  - Effect sizes
  - Significance levels

#### Individual Diagnostic Figures
- `fig_tci_hist_*.png` — TCI histogram
- `fig_tci_cdf_*.png` — TCI empirical CDF
- `fig_tci_kernel_*.png` — TCI KDE
- `fig_mci_hist_*.png` — MCI histogram
- `fig_mci_cdf_*.png` — MCI empirical CDF
- `fig_tci_vs_mci_per_dataset_*.png` — Dataset-level TCI vs MCI

---

### 4. Structural Analysis Figures (Real Data)

These figures analyse **only the real datasets**, independent of null models:

- `fig_real_distance_kde.png` — Distance distribution in *h(q)* space
- `fig_real_dendrogram.png` — Hierarchical clustering of datasets
- `fig_real_tci_heatmap.png` — Pairwise coherence heatmap
- `fig_real_pca_hq.png` — PCA of mean *h(q)* curves
- `fig_real_hq_overlay_raw.png` — Overlay of mean *h(q)* curves
- `fig_real_domain_means.png` — Domain-averaged MCI values

These figures correspond to the structural analyses described in the *Results from Real Measurements* section of the paper.

---

## 🔁 Reproducibility Guarantee

All files in this directory can be regenerated **exactly** using:

- The provided `FULL_DETAILS_utmf_v5.1_*.json`
- The UTMF v5.1 null-test cell
- The UTMF v5.1 validation cell

No additional preprocessing, parameter tuning, or external configuration is required.

The JSON file embeds:
- The full `jedi_mfdfa` implementation (order 0)
- All scale and *q*-grid definitions
- Archived time series used for TCI
- All subset-level *h(q)* spectra used for MCI

---

## 📖 Intended Use

This folder is designed to support:
- Independent verification of the paper’s claims
- Re-analysis using alternative statistical tests
- Comparison with new datasets or null models
- Methodological reuse of the TCI/MCI framework

Researchers are encouraged to treat the CSV files as the **primary scientific record**, with figures serving as visual summaries.

---

## 📌 Notes

- File timestamps encode the UTC generation time.
- Minor filename differences across runs reflect updated null-model realizations.
- All results reported in the paper correspond to the **latest CSV in this folder** at the time of DOI registration.

---

**UTMF v5.1 — Unified Temporal–Measurement Framework**  
Fully reproducible, domain-agnostic multifractal coherence analysis.

