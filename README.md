# BWTO2D - Bayesian Wavelet Topology Optimization 2D: A MATLAB GUI

## What is this repository for?

BWTO2D is a MATLAB-based graphical user interface (GUI) that provides an unsupervised framework for the automated extraction and characterization of curvilinear geological lineaments from potential field data. The application integrates a sophisticated workflow that combines the 2D Continuous Wavelet Transform (CWT) with a novel Bayesian Topology Optimization (BTO) engine, allowing for both semi-automatic and fully automated, unsupervised detection of lineaments.

---

## Key features include:

- **Advanced Wavelet Analysis:** Utilizes specialized mother wavelets, including Derivatives of Poisson, to analyze data at multiple scales and orientations.

- **Intelligent Feature Selection:** Employs an algorithm based on localized feature saliency and global weighted dissimilarity to select the most informative CWT-derived features.

- **Bayesian Topology Optimization (BTO):** Automatically finds the optimal set of hyperparameters for lineament extraction in an unsupervised manner by maximizing a Graph Representativeness Metric (GRM).

- **Comprehensive Visualization and Export:** Allows for viewing results at every stage and exporting them to common GIS formats including CSV, GeoTIFF, and Shapefiles.

The framework is designed for applications in mineral exploration, geophysical data analysis, and geological mapping.

---

## Please cite the software as:

> Abbassi, B., & Cheng, L.Z. (2025). Bayesian Wavelet Topology Optimization for Curvilinear Pattern Recognition in Potential Field Data. Computers & Geosciences.

---

## How do I get set up?

This program is designed to run on any Windows-based personal computer.

- **Memory (RAM):** A minimum of 8 GB of RAM is required. More RAM is recommended as it allows for processing larger images.
- **Storage:** A solid-state drive (SSD) with a non-volatile memory express (NVMe) interface is recommended due to the large matrix operations performed by the program.
- **MATLAB Version:** BWTO2D requires MATLAB R2024b.

**Setup:**  
To use the program, locate the M-Files in the current folder of MATLAB and then type `BWTO2D` in the MATLAB Command Window.

### Required MATLAB Toolboxes:

- Image Processing Toolbox
- Mapping Toolbox
- Statistics and Machine Learning Toolbox
- Signal Processing Toolbox
- Wavelet Toolbox

---

## Usage

The BWTO2D interface is organized into four main tabs:

- **Data Sets**
- **Spectral**
- **Lineaments**
- **Export**

### Data Sets Tab

This is the starting point for defining the study area and loading data.

#### Coordinates Panel

- **Coord Method:** Define the geographic extent of your project.
    - **Rectangular coordinates:** Uses a `.txt` file with four values (Min Lon, Max Lon, Min Lat, Max Lat) to define the boundaries.
    - **Automatic polygonization:** Automatically calculates the bounding box from a `.csv` or `.xyz` file containing data points.
- **Polygonize:** Click this button after selecting a method to load the coordinates.
- **Spacing (Arc-Seconds):** Determines the resolution of the analysis grid.
- **Data Filter:** Applies a Gaussian smoothing filter to the data for visualization.

#### Data Sets Panel

- **Point Data:** Click to load scattered point data from a `.csv` or `.xyz` file.
- **Interpolation Method:** Choose between Regular Interp (Nearest Neighbor) and Scattered Interp (Natural Neighbor) for gridding the data.
- **Xn / Yn:** Displays the resulting grid dimensions after polygonization.

### Spectral Tab

This tab is for feature extraction using the Continuous Wavelet Transform (CWT).

#### Spectral Feature Extraction Inputs

- **Point Data:** Check this box to include the loaded dataset for CWT analysis.
- **Merge:** Combines the selected datasets into a single input matrix for the CWT.

#### Continuous Wavelet Transform (CWT) Panel

- **Mother Wavelet:** Select the wavelet for the analysis. Options include Derivatives of Gaussian and Derivatives of Poisson.  
  *Note: The Merge button must be clicked again after changing the mother wavelet.*
- **CWT Parameters:** Configure the transform with Number of Scales (Na), Scale Dilation Factor (δ), Wavelet Smoothness Filter Ratio (σ), Number of CWT Angles (Nθ), Derivatives Orders (n and m), and Poisson Upward Continuation (Z).
- **CWT Button:** Executes the CWT process.

#### Spectral Feature Selection Panel

This panel reduces the number of CWT feature maps to an informative, non-redundant subset.

- **Parameters:** Set the Feature selection dimension (d), Dissimilarity (γ), and Feature Saliency Percentage (λ).
- **Feature Selection Button:** Runs the selection algorithm.

### Lineaments Tab

This section uses the extracted spectral features to detect and map lineaments via semi-automatic or fully automatic workflows.

#### Digitized Lineaments Panel

- Load pre-existing, known lineaments from an image or point data file for calculating the F-beta score to validate automatic results.

#### Lineament Extraction (Semi-Automatic)

- **Merge Inputs:** Select the features to use (P for Point Data, CWT Selected for spectral features) and click Merge.
- **Set Parameters:** Manually control Resolution Factor (R), Image Filter Ratio (σI), Step Filtering Width (w), Variability of Step Filtering Width (υ), Number of Step Filtering Angles (Nϕ), and Lineament confidence threshold (ξ).
- **Lineament Extraction Button:** Runs the detection process with the specified manual parameters.

#### Lineament Extraction by Bayesian Topology Optimization (BTO)

This is a fully automatic workflow that optimizes all critical hyperparameters. Using the Spectral tab beforehand is not required.

- **Merge Inputs:** Select the P (Point Data) checkbox and click Merge.
- **Select Mother Wavelet:** Choose the wavelet for the internal CWT calculation.
- **BTO Max Iterations:** Set the number of iterations for the optimizer.
- **Lineament Extraction (BTO) Button:** Starts the automated optimization and extraction process. The optimal hyperparameters found will populate the fields marked with an asterisk (*).

#### Display Panel

Visualize all results, including target data, manual and automatic lineaments, rose diagrams, and BTO diagnostic plots like BHO Performance and Hyperparameter Importance.

### Export Tab

Save raw data, feature maps, and final lineament products.

- **Data Selection:** Use the radio buttons to choose the dataset to export (e.g., P Data, CWT (Manual), All Lines (Auto)).
- **Export to CSV:** Saves gridded data as a three-column (X, Y, Value) text file, or lineaments as vertex coordinates.
- **Export to GeoTIFF / Shapefile:**
    - If a gridded dataset is selected (e.g., CWT maps, Densities), it saves as a GeoTIFF (.tif) file with spatial referencing.
    - If a lineament dataset is selected (e.g., All Lines, Deep Lines), it saves the vector data as an Esri Shapefile (.shp).

---

## Input/output formats

### Input Formats

- **Point Data:** `.csv` or Geosoft `.xyz` files containing three columns: Longitude, Latitude, and Value.
- **Coordinate Definition:**
    - `.txt` file containing four lines for minimum longitude, maximum longitude, minimum latitude, and maximum latitude.
    - `.csv` or `.xyz` file with point data for automatic boundary detection.
- **Target Lineaments (for validation):**
    - Image files (`.jpg`, `.png`, `.bmp`) accompanied by a `.txt` file of the same name containing georeferencing coordinates.
    - Point data files (`.csv`, `.xyz`) representing known lineaments.

### Output Formats

- **GIS Formats:**
    - **GeoTIFF (.tif):** For all gridded raster data, such as CWT feature maps and lineament density maps. This format preserves georeferencing information.
    - **Esri Shapefile (.shp):** For all vector lineament data (All, Deep, and Shallow Lines). This is the industry-standard format for vector GIS data.
- **Text Format:**
    - **CSV (.csv):** For exporting gridded data (X, Y, Value) or the vertex coordinates of vector lineaments.
- **MATLAB Figures:**
    - All plots generated within the application can be saved from their respective figure windows in standard MATLAB formats (`.fig`) or other image formats (`.png`, `.jpg`).

---

## Who do I talk to?

Bahman Abbassi  
Université du Québec en Abitibi-Témiscamingue  
Email: bahman.abbassi@uqat.ca

---

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
