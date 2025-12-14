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
- `NIST_elements` [Download](https://github.com/Jedi-Markus-Strive/UTMF-v5.1-Coherence/blob/master/downloads/NIST_elements.zip)** (8.5MB) (Complete dataset as used in UTMF v5.1, unzip and use the CSV's 57.5MB) for UTMF analysis.)
---
- **[NANOGrav Data Releases](https://zenodo.org/records/16051178)**  
- Pulsar timing residuals (e.g., `NG15yr narrowband` files).

  **Dataset used in UTMF v5.1 configuration:**
- `NANOGrav15yr_PulsarTiming_v2.1.0` [Download:](https://zenodo.org/records/16051178/files/NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz?download=1) (639MB) ([Unzip](https://github.com/Jedi-Markus-Strive/UTMF-v5.1-Coherence/blob/master/downloads/prepare_NANOGrav_15yr_data.ipynb) the file, use for UTMF analysis.)
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
