# BWTO2D - Bayesian Wavelet Topology Optimization 2D

> **A MATLAB GUI for Unsupervised Extraction of Geological Lineaments**

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Input/Output Formats](#inputoutput-formats)
- [Citation](#citation)
- [Contact](#contact)
- [License](#license)

---

## Overview

**BWTO2D** is a MATLAB-based graphical user interface (GUI) for automated, unsupervised extraction and characterization of curvilinear geological lineaments from potential field data. It integrates 2D Continuous Wavelet Transform (CWT) and a novel Bayesian Topology Optimization (BTO) engine for semi-automatic or fully automated detection.

---

## Features

- **Advanced Wavelet Analysis:**  
  Uses custom mother wavelets (e.g., Derivatives of Poisson) for multiscale, multi-orientation analysis.
- **Intelligent Feature Selection:**  
  Algorithm combines local feature saliency and global weighted dissimilarity.
- **Bayesian Topology Optimization (BTO):**  
  Automatically finds optimal hyperparameters for lineament extraction via a Graph Representativeness Metric (GRM).
- **Comprehensive Visualization & Export:**  
  Stepwise visualization; export to CSV, GeoTIFF, Shapefile.
- **Typical Applications:**  
  Mineral exploration, geophysical data analysis, and geological mapping.

---

## System Requirements

- **OS:** Windows  
- **RAM:** 8 GB minimum (more recommended for large data)  
- **Storage:** SSD (NVMe recommended)  
- **MATLAB:** R2024b  
- **MATLAB Toolboxes:**  
  - Image Processing Toolbox  
  - Mapping Toolbox  
  - Statistics and Machine Learning Toolbox  
  - Signal Processing Toolbox  
  - Wavelet Toolbox

---

## Installation

1. Clone or download this repository.
2. Place all `.m` files in your MATLAB current folder.
3. In MATLAB Command Window, run:
    ```matlab
    BWTO2D
    ```

---

## Usage

The GUI is organized into four tabs:  
**Data Sets, Spectral, Lineaments, Export**

### Data Sets Tab

- Define study area and load input data.
- Set geographic bounds:  
  - `.txt` file (min/max lon/lat)  
  - or automatic from `.csv`/`.xyz`
- Choose gridding method and resolution.
- Load point data (`.csv`/`.xyz`).

### Spectral Tab

- Extract features via CWT (choose wavelet type, parameters).
- Feature selection to reduce redundancy.

### Lineaments Tab

- Detect and map lineaments:
  - **Semi-automatic:** Manually tune parameters.
  - **Fully automatic:** Run BTO for hyperparameter optimization.
- Load ground-truth for validation (optional).

### Export Tab

- Save results as CSV, GeoTIFF, or Shapefile.
- Export all figures as `.fig`, `.png`, or `.jpg`.

---

## Input/Output Formats

### Input

- **Point Data:** `.csv` or Geosoft `.xyz` (Longitude, Latitude, Value)
- **Coordinate Definition:**  
  - `.txt`: 4 lines (min lon, max lon, min lat, max lat)  
  - or `.csv`/`.xyz`: auto-boundary detection
- **Lineament Validation:**  
  - Image files (`.jpg`, `.png`, `.bmp`) + matching `.txt` for georeferencing  
  - Point data files (`.csv`, `.xyz`)

### Output

- **Raster Data:** GeoTIFF (`.tif`)
- **Vector Data:** Esri Shapefile (`.shp`)
- **Table:** CSV (`.csv`)
- **Figures:** `.fig`, `.png`, `.jpg`

---

## Citation

If you use this software, please cite:

> Abbassi, B., & Cheng, L.Z. (2025). Bayesian Wavelet Topology Optimization for Curvilinear Pattern Recognition in Potential Field Data. *Computers & Geosciences*.

---

## Contact

**Bahman Abbassi**  
Université du Québec en Abitibi-Témiscamingue  
📧 [bahman.abbassi@uqat.ca](mailto:bahman.abbassi@uqat.ca)

---

## License

This project is licensed under the [GNU GPL v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).

> This program is distributed in the hope that it will be useful, but **WITHOUT ANY WARRANTY**; without even the implied warranty of **MERCHANTABILITY** or **FITNESS FOR A PARTICULAR PURPOSE**.

---

