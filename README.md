# BWTO2D - Bayesian Wavelet Topology Optimization 2D: A MATLAB GUI

**Dr. Bahman Abbassi**
*Postdoctoral Fellow in Mining Engineering, Université du Québec en Abitibi-Témiscamingue (UQAT)*
*Email: [bahman.abbassi@uqat.ca](mailto:bahman.abbassi@uqat.ca)* [cite: 2, 463, 469]

---

### What is this repository for?

**BWTO2D** is a MATLAB-based graphical user interface (GUI) that provides an unsupervised framework for the automated extraction and characterization of curvilinear geological lineaments from potential field data[cite: 16, 282, 502]. The application integrates a sophisticated workflow that combines the 2D Continuous Wavelet Transform (CWT) with a novel Bayesian Topology Optimization (BTO) engine, allowing for both semi-automatic and fully automated, unsupervised detection of lineaments[cite: 17, 503, 504].

Key features include:
* **Advanced Wavelet Analysis**: Utilizes specialized mother wavelets, including Derivatives of Poisson, to analyze data at multiple scales and orientations[cite: 50, 507].
* **Intelligent Feature Selection**: Employs an algorithm based on localized feature saliency and global weighted dissimilarity to select the most informative CWT-derived features[cite: 18, 86, 508].
* **Bayesian Topology Optimization (BTO)**: Automatically finds the optimal set of hyperparameters for lineament extraction in an unsupervised manner by maximizing a Graph Representativeness Metric (GRM)[cite: 21, 172, 482].
* **Comprehensive Visualization and Export**: Allows for viewing results at every stage and exporting them to common GIS formats including CSV, GeoTIFF, and Shapefiles[cite: 512, 818].

The framework is designed for applications in mineral exploration, geophysical data analysis, and geological mapping [cite: 22, 513-517].

Please cite the software as:
> Abbassi, B., & Cheng, L.Z. (2025). Bayesian Wavelet Topology Optimization for Curvilinear Pattern Recognition in Potential Field Data. *Computers & Geosciences*. [cite: 1, 2]

---
### How do I get set up?

This program is designed to run on any Windows-based personal computer[cite: 297, 522].
* **Memory (RAM)**: A minimum of **8 GB of RAM** is required[cite: 297, 522]. More RAM is recommended as it allows for processing larger images[cite: 298, 523].
* **Storage**: A **solid-state drive (SSD)** with a non-volatile memory express (NVMe) interface is recommended due to the large matrix operations performed by the program[cite: 299, 525].
* **MATLAB Version**: BWTO2D requires **MATLAB R2024b**[cite: 301, 527].
* **Setup**: To use the program, locate the M-Files in the current folder of MATLAB and then type `BWTO2D` in the MATLAB Command Window[cite: 528].

**Required MATLAB Toolboxes**[cite: 544]:
* Image Processing Toolbox [cite: 545]
* Mapping Toolbox [cite: 546]
* Statistics and Machine Learning Toolbox [cite: 547]
* Signal Processing Toolbox [cite: 548]
* Wavelet Toolbox [cite: 549]

---
### Usage

The BWTO2D interface is organized into four main tabs: **Data Sets**, **Spectral**, **Lineaments**, and **Export**[cite: 559, 561].

#### **Data Sets Tab** [cite: 562]
This is the starting point for defining the study area and loading data[cite: 629].
* **Coordinates Panel**[cite: 633]:
    * **Coord Method**: Define the geographic extent of your project[cite: 634].
        * `Rectangular coordinates`: Uses a `.txt` file with four values (Min Lon, Max Lon, Min Lat, Max Lat) to define the boundaries [cite: 636-637].
        * `Automatic polygonization`: Automatically calculates the bounding box from a `.csv` or `.xyz` file containing data points [cite: 644-645].
    * **Polygonize**: Click this button after selecting a method to load the coordinates[cite: 640, 650].
    * **Spacing (Arc-Seconds)**: Determines the resolution of the analysis grid[cite: 653].
    * **Data Filter**: Applies a Gaussian smoothing filter to the data for visualization[cite: 672].
* **Data Sets Panel**[cite: 656]:
    * **Point Data**: Click to load scattered point data from a `.csv` or `.xyz` file[cite: 657, 659].
    * **Interpolation Method**: Choose between `Regular Interp` (Nearest Neighbor) and `Scattered Interp` (Natural Neighbor) for gridding the data [cite: 661-663].
    * **Xn / Yn**: Displays the resulting grid dimensions after polygonization[cite: 654].

#### **Spectral Tab** [cite: 563]
This tab is for feature extraction using the Continuous Wavelet Transform (CWT)[cite: 679].
* **Spectral Feature Extraction Inputs**[cite: 684]:
    * **Point Data**: Check this box to include the loaded dataset for CWT analysis[cite: 686].
    * **Merge**: Combines the selected datasets into a single input matrix for the CWT[cite: 689].
* **Continuous Wavelet Transform (CWT) Panel**[cite: 692]:
    * **Mother Wavelet**: Select the wavelet for the analysis. Options include `Derivatives of Gaussian` and `Derivatives of Poisson` [cite: 696-698]. Note: The `Merge` button must be clicked again after changing the mother wavelet[cite: 700].
    * **CWT Parameters**: Configure the transform with `Number of Scales (Na)` [cite: 702], `Scale Dilation Factor (δ)` [cite: 705], `Wavelet Smoothness Filter Ratio (σ)` [cite: 708], `Number of CWT Angles (Nθ)` [cite: 709], `Derivatives Orders (n and m)` [cite: 711], and `Poisson Upward Continuation (Z)`[cite: 181, 714].
    * **CWT Button**: Executes the CWT process[cite: 718].
* **Spectral Feature Selection Panel**[cite: 722]:
    * This panel reduces the number of CWT feature maps to an informative, non-redundant subset[cite: 723].
    * **Parameters**: Set the `Feature selection dimension (d)` [cite: 181, 727], `Dissimilarity (γ)` [cite: 181, 728], and `Feature Saliency Percentage (λ)`[cite: 181, 729].
    * **Feature Selection Button**: Runs the selection algorithm[cite: 734].

#### **Lineaments Tab** [cite: 565]
This section uses the extracted spectral features to detect and map lineaments via semi-automatic or fully automatic workflows[cite: 743].
* **Digitized Lineaments Panel**[cite: 747]:
    * Load pre-existing, known lineaments from an image or point data file for calculating the **F-beta score** to validate automatic results [cite: 748-749].
* **Lineament Extraction (Semi-Automatic)**[cite: 765]:
    * **Merge Inputs**: Select the features to use (`P` for Point Data, `CWT Selected` for spectral features) and click `Merge` [cite: 766-769].
    * **Set Parameters**: Manually control `Resolution Factor (R)` [cite: 771], `Image Filter Ratio (σI)` [cite: 773], `Step Filtering Width (w)` [cite: 774], `Variability of Step Filtering Width (υ)` [cite: 777], `Number of Step Filtering Angles (Nϕ)` [cite: 779], and `Lineament confidence threshold (ξ)`[cite: 780].
    * **Lineament Extraction Button**: Runs the detection process with the specified manual parameters[cite: 781].
* **Lineament Extraction by Bayesian Topology Optimization (BTO)**[cite: 784]:
    * This is a **fully automatic** workflow that optimizes all critical hyperparameters[cite: 785]. Using the Spectral tab beforehand is not required[cite: 786].
    * **Merge Inputs**: Select the `P` (Point Data) checkbox and click `Merge`[cite: 791].
    * **Select Mother Wavelet**: Choose the wavelet for the internal CWT calculation[cite: 792].
    * **BTO Max Iterations**: Set the number of iterations for the optimizer[cite: 795].
    * **Lineament Extraction (BTO) Button**: Starts the automated optimization and extraction process[cite: 802]. The optimal hyperparameters found will populate the fields marked with an asterisk (*)[cite: 804].
* **Display Panel**:
    * Visualize all results, including target data, manual and automatic lineaments, rose diagrams, and BTO diagnostic plots like `BHO Performance` and `Hyperparameter Importance` [cite: 807-815].

#### **Export Tab** [cite: 567]
Save raw data, feature maps, and final lineament products[cite: 818].
* **Data Selection**: Use the radio buttons to choose the dataset to export (e.g., `P Data`, `CWT (Manual)`, `All Lines (Auto)`) [cite: 821-822].
* **Export to CSV**: Saves gridded data as a three-column (X, Y, Value) text file, or lineaments as vertex coordinates [cite: 830-833].
* **Export to GeoTIFF / Shapefile**:
    * If a **gridded dataset** is selected (e.g., CWT maps, Densities), it saves as a **GeoTIFF (.tif)** file with spatial referencing[cite: 846].
    * If a **lineament dataset** is selected (e.g., All Lines, Deep Lines), it saves the vector data as an **Esri Shapefile (.shp)**[cite: 848].

---
### Input/output formats

#### **Input Formats**
* **Point Data**: `.csv` or Geosoft `.xyz` files containing three columns: **Longitude**, **Latitude**, and **Value**[cite: 657].
* **Coordinate Definition**:
    * A `.txt` file containing four lines for minimum longitude, maximum longitude, minimum latitude, and maximum latitude[cite: 637].
    * A `.csv` or `.xyz` file with point data for automatic boundary detection[cite: 645].
* **Target Lineaments (for validation)**:
    * Image files (`.jpg`, `.png`, `.bmp`) accompanied by a `.txt` file of the same name containing georeferencing coordinates [cite: 753-754].
    * Point data files (`.csv`, `.xyz`) representing known lineaments[cite: 755].

#### **Output Formats**
* **GIS Formats**:
    * **GeoTIFF (.tif)**: For all gridded raster data, such as CWT feature maps and lineament density maps. This format preserves georeferencing information[cite: 839, 846].
    * **Esri Shapefile (.shp)**: For all vector lineament data (All, Deep, and Shallow Lines). This is the industry-standard format for vector GIS data [cite: 848-849].
* **Text Format**:
    * **CSV (.csv)**: For exporting gridded data (X, Y, Value) or the vertex coordinates of vector lineaments [cite: 830-833].
* **MATLAB Figures**:
    * All plots generated within the application can be saved from their respective figure windows in standard MATLAB formats (`.fig`) or other image formats (`.png`, `.jpg`).

---
### Who do I talk to?

* **Bahman Abbassi**, Université du Québec en Abitibi-Témiscamingue, [bahman.abbassi@uqat.ca](mailto:bahman.abbassi@uqat.ca) [cite: 463, 466, 469]

---
### License

This program is free software: you can redistribute it and/or modify it under the terms of the **GNU General Public License** as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but **WITHOUT ANY WARRANTY**; without even the implied warranty of **MERCHANTABILITY** or **FITNESS FOR A PARTICULAR PURPOSE**. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
