# BWTO2D - Bayesian Wavelet Topology Optimization 2D: A MATLAB GUI

**Dr. Bahman Abbassi**
*Postdoctoral Fellow in Mining Engineering, Université du Québec en Abitibi-Témiscamingue (UQAT)*
*Email: [bahman.abbassi@uqat.ca](mailto:bahman.abbassi@uqat.ca)* [cite: 2, 463, 469]

---

### What is this repository for?

[cite_start]**BWTO2D** is a MATLAB-based graphical user interface (GUI) that provides an unsupervised framework for the automated extraction and characterization of curvilinear geological lineaments from potential field data[cite: 16, 282, 502]. [cite_start]The application integrates a sophisticated workflow that combines the 2D Continuous Wavelet Transform (CWT) with a novel Bayesian Topology Optimization (BTO) engine, allowing for both semi-automatic and fully automated, unsupervised detection of lineaments[cite: 17, 503, 504].

Key features include:
* [cite_start]**Advanced Wavelet Analysis**: Utilizes specialized mother wavelets, including Derivatives of Poisson, to analyze data at multiple scales and orientations[cite: 50, 507].
* [cite_start]**Intelligent Feature Selection**: Employs an algorithm based on localized feature saliency and global weighted dissimilarity to select the most informative CWT-derived features[cite: 18, 86, 508].
* [cite_start]**Bayesian Topology Optimization (BTO)**: Automatically finds the optimal set of hyperparameters for lineament extraction in an unsupervised manner by maximizing a Graph Representativeness Metric (GRM)[cite: 21, 172, 482].
* [cite_start]**Comprehensive Visualization and Export**: Allows for viewing results at every stage and exporting them to common GIS formats including CSV, GeoTIFF, and Shapefiles[cite: 512, 818].

[cite_start]The framework is designed for applications in mineral exploration, geophysical data analysis, and geological mapping [cite: 22, 513-517].

Please cite the software as:
> Abbassi, B., & Cheng, L.Z. (2025). Bayesian Wavelet Topology Optimization for Curvilinear Pattern Recognition in Potential Field Data. [cite_start]*Computers & Geosciences*. [cite: 1, 2]

---
### How do I get set up?

[cite_start]This program is designed to run on any Windows-based personal computer[cite: 297, 522].
* [cite_start]**Memory (RAM)**: A minimum of **8 GB of RAM** is required[cite: 297, 522]. [cite_start]More RAM is recommended as it allows for processing larger images[cite: 298, 523].
* [cite_start]**Storage**: A **solid-state drive (SSD)** with a non-volatile memory express (NVMe) interface is recommended due to the large matrix operations performed by the program[cite: 299, 525].
* [cite_start]**MATLAB Version**: BWTO2D requires **MATLAB R2024b**[cite: 301, 527].
* [cite_start]**Setup**: To use the program, locate the M-Files in the current folder of MATLAB and then type `BWTO2D` in the MATLAB Command Window[cite: 528].

[cite_start]**Required MATLAB Toolboxes**[cite: 544]:
* [cite_start]Image Processing Toolbox [cite: 545]
* [cite_start]Mapping Toolbox [cite: 546]
* [cite_start]Statistics and Machine Learning Toolbox [cite: 547]
* [cite_start]Signal Processing Toolbox [cite: 548]
* [cite_start]Wavelet Toolbox [cite: 549]

---
### Usage

[cite_start]The BWTO2D interface is organized into four main tabs: **Data Sets**, **Spectral**, **Lineaments**, and **Export**[cite: 559, 561].

#### [cite_start]**Data Sets Tab** [cite: 562]
[cite_start]This is the starting point for defining the study area and loading data[cite: 629].
* [cite_start]**Coordinates Panel**[cite: 633]:
    * [cite_start]**Coord Method**: Define the geographic extent of your project[cite: 634].
        * [cite_start]`Rectangular coordinates`: Uses a `.txt` file with four values (Min Lon, Max Lon, Min Lat, Max Lat) to define the boundaries [cite: 636-637].
        * [cite_start]`Automatic polygonization`: Automatically calculates the bounding box from a `.csv` or `.xyz` file containing data points [cite: 644-645].
    * [cite_start]**Polygonize**: Click this button after selecting a method to load the coordinates[cite: 640, 650].
    * [cite_start]**Spacing (Arc-Seconds)**: Determines the resolution of the analysis grid[cite: 653].
    * [cite_start]**Data Filter**: Applies a Gaussian smoothing filter to the data for visualization[cite: 672].
* [cite_start]**Data Sets Panel**[cite: 656]:
    * [cite_start]**Point Data**: Click to load scattered point data from a `.csv` or `.xyz` file[cite: 657, 659].
    * [cite_start]**Interpolation Method**: Choose between `Regular Interp` (Nearest Neighbor) and `Scattered Interp` (Natural Neighbor) for gridding the data [cite: 661-663].
    * [cite_start]**Xn / Yn**: Displays the resulting grid dimensions after polygonization[cite: 654].

#### [cite_start]**Spectral Tab** [cite: 563]
[cite_start]This tab is for feature extraction using the Continuous Wavelet Transform (CWT)[cite: 679].
* [cite_start]**Spectral Feature Extraction Inputs**[cite: 684]:
    * [cite_start]**Point Data**: Check this box to include the loaded dataset for CWT analysis[cite: 686].
    * [cite_start]**Merge**: Combines the selected datasets into a single input matrix for the CWT[cite: 689].
* [cite_start]**Continuous Wavelet Transform (CWT) Panel**[cite: 692]:
    * **Mother Wavelet**: Select the wavelet for the analysis. [cite_start]Options include `Derivatives of Gaussian` and `Derivatives of Poisson` [cite: 696-698]. [cite_start]Note: The `Merge` button must be clicked again after changing the mother wavelet[cite: 700].
    * [cite_start]**CWT Parameters**: Configure the transform with `Number of Scales (Na)` [cite: 702][cite_start], `Scale Dilation Factor (δ)` [cite: 705][cite_start], `Wavelet Smoothness Filter Ratio (σ)` [cite: 708][cite_start], `Number of CWT Angles (Nθ)` [cite: 709][cite_start], `Derivatives Orders (n and m)` [cite: 711][cite_start], and `Poisson Upward Continuation (Z)`[cite: 181, 714].
    * [cite_start]**CWT Button**: Executes the CWT process[cite: 718].
* [cite_start]**Spectral Feature Selection Panel**[cite: 722]:
    * [cite_start]This panel reduces the number of CWT feature maps to an informative, non-redundant subset[cite: 723].
    * [cite_start]**Parameters**: Set the `Feature selection dimension (d)` [cite: 181, 727][cite_start], `Dissimilarity (γ)` [cite: 181, 728][cite_start], and `Feature Saliency Percentage (λ)`[cite: 181, 729].
    * [cite_start]**Feature Selection Button**: Runs the selection algorithm[cite: 734].

#### [cite_start]**Lineaments Tab** [cite: 565]
[cite_start]This section uses the extracted spectral features to detect and map lineaments via semi-automatic or fully automatic workflows[cite: 743].
* [cite_start]**Digitized Lineaments Panel**[cite: 747]:
    * [cite_start]Load pre-existing, known lineaments from an image or point data file for calculating the **F-beta score** to validate automatic results [cite: 748-749].
* [cite_start]**Lineament Extraction (Semi-Automatic)**[cite: 765]:
    * [cite_start]**Merge Inputs**: Select the features to use (`P` for Point Data, `CWT Selected` for spectral features) and click `Merge` [cite: 766-769].
    * [cite_start]**Set Parameters**: Manually control `Resolution Factor (R)` [cite: 771][cite_start], `Image Filter Ratio (σI)` [cite: 773][cite_start], `Step Filtering Width (w)` [cite: 774][cite_start], `Variability of Step Filtering Width (υ)` [cite: 777][cite_start], `Number of Step Filtering Angles (Nϕ)` [cite: 779][cite_start], and `Lineament confidence threshold (ξ)`[cite: 780].
    * [cite_start]**Lineament Extraction Button**: Runs the detection process with the specified manual parameters[cite: 781].
* [cite_start]**Lineament Extraction by Bayesian Topology Optimization (BTO)**[cite: 784]:
    * [cite_start]This is a **fully automatic** workflow that optimizes all critical hyperparameters[cite: 785]. [cite_start]Using the Spectral tab beforehand is not required[cite: 786].
    * [cite_start]**Merge Inputs**: Select the `P` (Point Data) checkbox and click `Merge`[cite: 791].
    * [cite_start]**Select Mother Wavelet**: Choose the wavelet for the internal CWT calculation[cite: 792].
    * [cite_start]**BTO Max Iterations**: Set the number of iterations for the optimizer[cite: 795].
    * [cite_start]**Lineament Extraction (BTO) Button**: Starts the automated optimization and extraction process[cite: 802]. [cite_start]The optimal hyperparameters found will populate the fields marked with an asterisk (*)[cite: 804].
* **Display Panel**:
    * [cite_start]Visualize all results, including target data, manual and automatic lineaments, rose diagrams, and BTO diagnostic plots like `BHO Performance` and `Hyperparameter Importance` [cite: 807-815].

#### [cite_start]**Export Tab** [cite: 567]
[cite_start]Save raw data, feature maps, and final lineament products[cite: 818].
* [cite_start]**Data Selection**: Use the radio buttons to choose the dataset to export (e.g., `P Data`, `CWT (Manual)`, `All Lines (Auto)`) [cite: 821-822].
* [cite_start]**Export to CSV**: Saves gridded data as a three-column (X, Y, Value) text file, or lineaments as vertex coordinates [cite: 830-833].
* **Export to GeoTIFF / Shapefile**:
    * [cite_start]If a **gridded dataset** is selected (e.g., CWT maps, Densities), it saves as a **GeoTIFF (.tif)** file with spatial referencing[cite: 846].
    * [cite_start]If a **lineament dataset** is selected (e.g., All Lines, Deep Lines), it saves the vector data as an **Esri Shapefile (.shp)**[cite: 848].

---
### Input/output formats

#### **Input Formats**
* [cite_start]**Point Data**: `.csv` or Geosoft `.xyz` files containing three columns: **Longitude**, **Latitude**, and **Value**[cite: 657].
* **Coordinate Definition**:
    * [cite_start]A `.txt` file containing four lines for minimum longitude, maximum longitude, minimum latitude, and maximum latitude[cite: 637].
    * [cite_start]A `.csv` or `.xyz` file with point data for automatic boundary detection[cite: 645].
* **Target Lineaments (for validation)**:
    * [cite_start]Image files (`.jpg`, `.png`, `.bmp`) accompanied by a `.txt` file of the same name containing georeferencing coordinates [cite: 753-754].
    * [cite_start]Point data files (`.csv`, `.xyz`) representing known lineaments[cite: 755].

#### **Output Formats**
* **GIS Formats**:
    * **GeoTIFF (.tif)**: For all gridded raster data, such as CWT feature maps and lineament density maps. [cite_start]This format preserves georeferencing information[cite: 839, 846].
    * **Esri Shapefile (.shp)**: For all vector lineament data (All, Deep, and Shallow Lines). [cite_start]This is the industry-standard format for vector GIS data [cite: 848-849].
* **Text Format**:
    * [cite_start]**CSV (.csv)**: For exporting gridded data (X, Y, Value) or the vertex coordinates of vector lineaments [cite: 830-833].
* **MATLAB Figures**:
    * All plots generated within the application can be saved from their respective figure windows in standard MATLAB formats (`.fig`) or other image formats (`.png`, `.jpg`).

---
### Who do I talk to?

* [cite_start]**Bahman Abbassi**, Université du Québec en Abitibi-Témiscamingue, [bahman.abbassi@uqat.ca](mailto:bahman.abbassi@uqat.ca) [cite: 463, 466, 469]

---
### License

This program is free software: you can redistribute it and/or modify it under the terms of the **GNU General Public License** as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but **WITHOUT ANY WARRANTY**; without even the implied warranty of **MERCHANTABILITY** or **FITNESS FOR A PARTICULAR PURPOSE**. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
