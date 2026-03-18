# Analysis code — Multi-modal myocardial tissue characterization in cats

Code repository accompanying:

> Spanjersberg et al. "Multi-modal myocardial tissue characterization in cats" *(eBioMedicine, in revision)*

---

## Contents

### QuPath scripts (QuPath version 0.5.1)

All `.groovy` scripts run inside [QuPath](https://qupath.github.io). Each script requires the relevant pixel classifier(s) or model file(s) to be present — see **Models** below.

| Script | Description |
|--------|-------------|
| `Adipocyte_detection.groovy` | Quantifies adipose tissue surface area and individual adipocyte size in HE-stained sections. Applies two pixel classifiers: one for tissue/adipocyte area quantification at low resolution and one for adipocyte membrane detection at high resolution. Includes automated filtering based on circularity and size thresholds. |
| `Fibrosis_detection.groovy` | Quantifies interstitial fibrosis in Masson's trichrome-stained sections. Applies an automated 700 µm inward erosion step to exclude epicardial and endocardial connective tissue before quantification. **Requires manual review** of excluded regions before final measurements are applied. |
| `Custom_Cellpose_nucleus_detection_HE_.groovy` | Detects myocardial nuclei in HE-stained sections using a custom-trained Cellpose model. Extracts the hematoxylin channel via color deconvolution prior to detection. Returns detected nuclei as annotations with shape and intensity measurements. Requires the [QuPath Cellpose extension](https://github.com/BIOP/qupath-extension-cellpose). Adapted from a template by Olivier Burri (BIOP, EPFL). |
| `Cellpose_training.groovy` | Helper script used to prepare training data for the custom Cellpose nucleus detection model. |

### R script

| Script | Description |
|--------|-------------|
| `Statistical_analysis.R` | Statistical analysis and figure generation for all main and supplemental figures. Includes morphometric comparisons across clinical outcome groups (ATE, HF, control), fibrosis quantification, adipocyte size analysis, and nuclear morphometry. Written in R; uses ggplot2 and ComplexHeatmap. |

---

## Models

Pixel classifier files (`.json`) and the custom Cellpose model are archived on Zenodo:

> **[INSERT ZENODO DOI]**

Download and place files as follows:

- QuPath pixel classifiers → `classifiers/pixel_classifiers/` inside your QuPath project folder
- Cellpose model → `models/` and update the `pathModel` path in `Custom_Cellpose_nucleus_detection_HE_.groovy`

---

## Requirements

- [QuPath](https://qupath.github.io) v0.5.1
- [QuPath Cellpose extension](https://github.com/BIOP/qupath-extension-cellpose) (for nucleus detection only)
- [Cellpose](https://github.com/MouseLand/cellpose) (installed and accessible from QuPath)
- R ≥ 4.0, with packages: `ggplot2`, `ComplexHeatmap`, `dplyr`

---

## Contact

For questions about the analysis code, contact the corresponding author via the journal article.
