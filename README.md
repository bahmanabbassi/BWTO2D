# BWTO2D: Bayesian Wavelet Topology Optimization for Curvilinear Pattern Recognition in Magnetic Data

## Publication

> Abbassi, B., Cheng, L.-Z., Legault, M. (2026). Bayesian Wavelet Topology Optimization for Curvilinear Pattern Recognition in Magnetic Data. *Computers & Geosciences*. <!-- doi: TBD -->

**Authors:** Bahman Abbassi, Li-Zhen Cheng, Marc Legault

**Affiliation:** Université du Québec en Abitibi-Témiscamingue (UQAT), QC J9X 5E4, Canada

**Contact:** [bahman.abbassi@uqat.ca](mailto:bahman.abbassi@uqat.ca)

**Version:** 1.2 | **Language:** MATLAB | **Size:** ~590 KB (source code)

**Repository:** [https://github.com/bahmanabbassi/BWTO2D](https://github.com/bahmanabbassi/BWTO2D)

## Abstract

We present an unsupervised framework for automated extraction and characterization of curvilinear geological patterns from magnetic data. The method integrates a Poisson kernel-based Continuous Wavelet Transform (kCWT), feature selection, and Bayesian Topology Optimization (BTO). The Poisson kCWT generates multidirectional coefficient maps and a feature selection algorithm identifies a compact, informative, and non-redundant subset by balancing feature saliency and dissimilarity. An ablation study demonstrates that this targeted selection outperforms Principal Component Analysis (PCA) by preserving pure directional boundaries essential for topological coherence. The selected features are processed with directional step filtering and adaptive thresholding to delineate curvilinear patterns. These patterns are transformed into a graph representation and evaluated using a novel Graph Representativeness Metric (GRM), an unsupervised objective function that quantifies structural coherence without requiring ground truth. BTO systematically optimizes critical hyperparameters to maximize GRM, eliminating manual tuning. Validation against expert-digitized curvilinear patterns demonstrates a strong correlation between GRM and supervised metrics, including Fᵦ score and fracture intensity ratio (R_P21), confirming GRM's effectiveness as a proxy for detection quality. The framework is evaluated using synthetic magnetic data derived from the Widgiemooltha Dome geological model (Australia) for algorithmic assessment, sensitivity analysis, and noise robustness testing, demonstrating graceful performance degradation under noise, while maintaining strong GRM–Fᵦ correlations. Application to aeromagnetic datasets from the Blake River Group (BRG) area in Abitibi Greenstone Belt (Quebec, Canada) confirms that the method supports fast, reproducible unsupervised structural mapping in complex geological settings.

**Keywords:** curvilinear pattern recognition · magnetic data · continuous wavelet transform · Bayesian optimization · graph topology

---

## Table of Contents

1. [Description](#description)
2. [Requirements](#requirements)
3. [Repository Structure](#repository-structure)
4. [Quick Start](#quick-start)
5. [Detailed Workflow](#detailed-workflow)
   - [Step 1 — Define Coordinates](#step-1--define-coordinates)
   - [Step 2 — Load Point Data](#step-2--load-point-data)
   - [Step 3 — Prepare kCWT Inputs (Auto)](#step-3--prepare-kcwt-inputs-auto)
   - [Step 4 — Load Ground Truth (Optional)](#step-4--load-ground-truth-optional)
   - [Step 5 — Configure BTO Hyperparameters](#step-5--configure-bto-hyperparameters)
   - [Step 6 — Run BTO](#step-6--run-bto)
   - [Step 7 — Inspect Results](#step-7--inspect-results)
   - [Step 8 — Export](#step-8--export)
6. [BTO Objective Function Pipeline](#bto-objective-function-pipeline)
7. [Parameter Reference](#parameter-reference)
   - [Fixed Parameters](#fixed-parameters)
   - [Optimizable Hyperparameters](#optimizable-hyperparameters)
   - [Bayesian Optimization Settings](#bayesian-optimization-settings)
   - [Validation Parameters](#validation-parameters)
8. [Display Options](#display-options)
9. [Export Formats](#export-formats)
10. [File Descriptions](#file-descriptions)
11. [Data Directories](#data-directories)
12. [Troubleshooting](#troubleshooting)
13. [Citation](#citation)
14. [License](#license)

---

## Description

BWTO2D is a MATLAB App Designer application that implements the unsupervised framework described in the paper above. It couples a four-stage image-processing pipeline — **Poisson kernel-based Continuous Wavelet Transform (kCWT)**, **kCWT coefficient selection and fusion**, **directional step filtering with adaptive thresholding**, and **graph-based topological representation** — with MATLAB's Bayesian optimization (`bayesopt`) to automatically tune 10 hyperparameters, maximizing a novel **Graph Representativeness Metric (GRM)**.

When ground-truth curvilinear patterns are available, the tool also tracks the **Fᵦ score**, **fracture intensity ratio (R_P21)**, **Hellinger similarity (HS)**, and **Jensen–Shannon divergence (JSD)** across iterations, enabling post-hoc evaluation of the unsupervised GRM proxy against supervised metrics.

The accompanying data directories contain results for the two case studies presented in the paper: synthetic magnetic data from the **Widgiemooltha Dome** geological model (including noise robustness and sensitivity analyses) and aeromagnetic data from the **Blake River Group (BRG)** in the Abitibi Greenstone Belt, Quebec, Canada.

### Key capabilities

- Geographic coordinate handling with automatic polygonization (shapefile or GeoTIFF)
- Poisson kCWT with Poisson Derivatives Mother Wavelet (PDMW) for multidirectional coefficient maps (multiscale optional; the paper uses a single scale)
- Saliency- and dissimilarity-based kCWT coefficient selection and fusion (Algorithm 1 in the paper)
- Directional step filtering and adaptive thresholding for curvilinear pattern delineation
- Graph-based topological representation with GRM scoring (Eqs. 16–19 in the paper)
- Bayesian Topology Optimization with configurable acquisition function, exploration ratio, and stopping criteria (Eq. 20)
- Ablation study comparing feature selection against PCA and an unreduced kCWT stack (`FeatureReductionComparison.m`, Figures 18 and 19)
- Live iteration visualization (curvilinear pattern map + GRM vs Fᵦ scatter)
- Best-so-far progression gallery
- Export to CSV, GeoTIFF, and ESRI Shapefile

---

## Requirements

| Requirement | Details |
|---|---|
| **MATLAB** | R2026a or later
| **Toolboxes** | Statistics and Machine Learning Toolbox (for `bayesopt`), Image Processing Toolbox, Mapping Toolbox (for GeoTIFF/shapefile I/O) |
| **OS** | Windows (developed and tested on Windows; see manuscript) |
| **RAM** | 8 GB minimum; 16+ GB recommended for large grids or many kCWT angles |
| **Program size** | ~590 KB (source code); ~5.6 MB including the bundled datasets |

---

## Repository Structure

The repository root holds four archives plus the licence and this file:

```
BWTO2D/
├── Matlab Codes.zip                # All source code (App, pipeline, metrics, utilities)
├── Feature Reduction Comparisons.zip # Ablation study: FS vs PCA vs unreduced stack (Figures 18-19)
├── Data 1 Synthetic.zip            # Widgiemooltha Dome synthetic benchmark data and results
├── Data 2 Blake River.zip          # Blake River Group (BRG) case study data and results
├── LICENSE                         # GNU General Public License v3.0
└── README.md                       # This file
```

Unzip all four archives in place before running anything. `Matlab Codes.zip` expands to:

```
Matlab Codes/
├── BWTO2D.mlapp                    # Main GUI application (App Designer)
├── Lineaments_Auto7.m              # BTO objective function (4-stage pipeline)
├── BHO_SensitivityAnalysis.m       # OAT sensitivity analysis driver (Figure 10)
├── add_noise_to_xyz.m              # Gaussian noise contamination for robustness tests
│
├── cwt2d.m                         # Core 2D kCWT engine (YAWTB-derived)
├── cwt2D_DerGus.m                  # kCWT — derivative of Gaussian (manual path)
├── cwt2D_DerPoisson.m              # kCWT — Poisson PDMW (manual path)
├── cwt2D_DerPoisson_BTO.m          # kCWT — Poisson PDMW (BTO path, Section 1.1)
├── dergauss2d.m                    # 2D derivative-of-Gaussian wavelet kernel
├── derpoisson2d.m                  # 2D Poisson Derivatives Mother Wavelet (PDMW)
│
├── CWTFeatureSelection.m           # kCWT coefficient selection and fusion (Algorithm 1)
├── step_Filtering.m                # Multi-angle directional step filter
├── getFaultDetection.m             # Hysteresis thresholding + morphological cleanup
├── FilterB.m                       # Scale-dependent Gaussian smoothing
├── nDstrb2D.m                      # 2D statistical distribution normalizer
│
├── calculateFbetaScore.m           # Bidirectional Fᵦ with tolerance buffer
├── calculateStructuralMetrics.m    # R_P21, orientation, HS, JSD metrics
├── calculateAngularSimilarity.m    # HS and JSD from orientation histograms
├── estimateTolerancePreBTO.m       # Automatic tolerance estimation before BTO
├── plotGRMvsGroundTruth.m          # Post-BTO correlation plots
├── createRoseDiagrams.m            # Orientation extraction from raster
├── createRoseDiagramsFromVectors.m # Orientation extraction from vector shapes
│
├── useBHO_BestOutputs.m            # Reconstruct globals from best iteration
├── convertRasterToShape.m          # Raster lineaments → shape structure
├── plotNewColorMap0_.m             # Colored lineament overlay renderer
├── toggleGTVisibilityCallback.m    # UI callback — toggle ground truth display
├── toggleGTVisibilityForTab.m      # UI callback — per-tab ground truth toggle
├── toggleVectorLineamentsVisibility.m # UI callback — vector patterns toggle
├── toggleImageVisibility.m         # UI callback — density/raster toggle
│
├── Surf1.m                         # Feature map gallery viewer
├── Surf1_Enhanced.m                # Enhanced point-data viewer
├── Custom_csvwrite.m               # CSV exporter with numbered suffix
├── geotiff.m                       # GeoTIFF exporter
│
├── getopts.m                       # YAWTB option parser
├── getyawtbprefs.m                 # YAWTB preference reader
├── yashow.m / yashow_cwt2d.m / …   # YAWTB visualization utilities
├── yapuls.m / yapuls2.m             # YAWTB pulse grid generators
└── vect.m / sphgrid.m / …          # YAWTB numeric helpers
```

`Feature Reduction Comparisons.zip` contains `FeatureReductionComparison.m`, the self-contained
ablation script that reproduces Figures 18 and 19, together with its saved outputs.

---

## Quick Start

Unzip `Matlab Codes.zip` (and the two data archives, if you want to reproduce the case studies), then:

```matlab
% 1. Add the unzipped code folder (and subfolders) to MATLAB's path
addpath(genpath('path/to/Matlab Codes'));

% 2. Launch the GUI
BWTO2D
```

The application opens a tabbed window with three tabs: **Data Sets**, **Lineaments**, and **Export**.

---

## Detailed Workflow

### Step 1 — Define Coordinates

On the **Data Sets** tab, under the **Coordinates** panel:

1. **Choose a coordinate method** from the dropdown:
   - **Rectangular coordinates** — Load a `.txt` file containing `[MinLon, MaxLon, MinLat, MaxLat]`.
   - **Automatic polygonization** — Load a `.shp` polygon or a georeferenced `.tif` / `.tiff` file. The bounding box and study-area polygon are extracted automatically.

2. **Set Spacing (Arc-Seconds)** — Controls the raster resolution. Smaller values produce finer grids but increase memory and computation time. Default: `1`.

3. **Set Data Filter** — Optional Gaussian pre-filter sigma for the coordinate grid.

4. Click **Polygonize** to build the internal coordinate meshgrids (`XI`, `YI`) and the polygon mask (`in`).

### Step 2 — Load Point Data

Still on the **Data Sets** tab, under the **Data Sets** panel:

1. Select an **interpolation method**:
   - `Regular Interp` — uses `griddata` with nearest-neighbor.
   - `Scattered Interp` — uses `scatteredInterpolant` with natural-neighbor.

2. Click **Point Data** and select a `.csv` or `.xyz` file with columns `[X, Y, Value]`.

3. The data is interpolated onto the coordinate grid, normalized (`nDstrb2D`), filtered (`FilterB`), and stored in `All_Data_Matrix`.

4. Use **Plot Data Sets** (under the **Display** panel) to visualize the loaded data.

### Step 3 — Prepare kCWT Inputs (Auto)

Switch to the **Lineaments** tab. In the BTO parameter panel:

1. Click **Merge** to prepare `CWT_Input_Features_Matrix_Auto` by selecting which loaded data layers to include as input for the Poisson kCWT.

### Step 4 — Load Ground Truth (Optional)

Under the target-lineaments panel on the **Lineaments** tab:

1. Set **Spacing** and grid dimensions for the target layer.
2. Click **Lineaments Shape File** to load a shapefile of known/mapped curvilinear patterns (e.g., expert-digitized lineaments).
3. The shapefile is rasterized onto the coordinate grid and stored in `All_Trg_Matrix`.

This enables Fᵦ score tracking during BTO and post-optimization validation.

### Step 5 — Configure BTO Hyperparameters

On the **Lineaments** tab, the left panel contains two groups:

#### Fixed Parameters (set once, not optimized)

| GUI Field | Paper symbol | Description | Value used in the paper |
|---|---|---|---|
| Number of Scales | *N*ₐ | kCWT scale count | 1 |
| Scale Dilation Factor | δ | Step between consecutive scales | 1 |
| Line Resolution | *R* | Resize factor for feature maps before extraction | 2 |
| CWT Reduction To | *d* | Target dimensionality after kCWT coefficient selection | 3 |
| Image Filter | σ_I | Gaussian sigma for pre-filtering (`FilterB`) | 0 |

`Nₐ = 1` and `δ = 1` restrict the kCWT to a single scale so that discrete scale jumps do not mix sources from different depths; depth control is left to the Poisson upward continuation `Z`. Set `Nₐ > 1` for coarser multiscale reconnaissance.

#### Hyperparameter Search Ranges (optimized by BTO)

Each field accepts a `[min, max]` range vector:

| GUI Field | Variable | Paper symbol | Type | Range in Table 1 | Description |
|---|---|---|---|---|---|
| CWT Angles | `NAngles` | *N*_θ | integer | {12 : 72} | Number of kCWT angular orientations |
| Order | `orderx` | *m*, *n* | categorical | {1,0} or {0,1} | Poisson derivative orders |
| Poisson Z | `Z` | *Z* | real | {1 : 5} | Upward continuation of the Poisson kernel |
| kCWT Filtering Ratio | `WSFR` | σ | real | {0 : 1} | kCWT support scaling factor |
| Center / Range (w) | `w` | *w̄* | real | {0.75·*w*_c : 2.5·*w*_c} | Average step-filter width |
| Step Filter Width Variability | `VSFW` | υ | real | {0 : 1} | Variance-based scaling of the filter window |
| Step Filter Angles | `NSFA` | *N*_φ | integer | {12 : 72} | Number of step-filter directional angles |
| Hysteresis Threshold | `SLC` | ξ | integer | {35 : 85} | Global confidence threshold (%) |
| Dissimilarity | `FDissimilarity` | γ | real | {0 : 1} | Feature selection dissimilarity weight |
| Saliency | `FSaliency` | λ | real | {0 : 100} | Feature selection saliency percentile |

The step-filter centre is formed internally as `w_c = 0.1 · min(X_Pix, Y_Pix) · R`, so the bounds of `w` scale with the grid. The adaptive width actually applied to a feature map is `w = w̄ + (υ·w̄)·z_var` (Eq. 14), which can exceed `2.5·w_c` where the local variance is low.

#### Bayesian Optimization Settings

| GUI Field | Description | Default | Recommended (Section 2.2.2) |
|---|---|---|---|
| BTO Iterations | Maximum objective evaluations | 50 | 50 |
| Acquisition Function | Surrogate-based selection strategy | EI per Second Plus | EI per Second Plus |
| Exploration Ratio | Explore (1.0) vs exploit (0.0) balance | 0.5 | 0.8 |
| Seed Points | Initial random evaluations | 4 | 3–6 |
| GP Active Set Size | Gaussian process approximation size | 300 | 150–250 |
| Max Time (min) | Wall-clock budget; 0 = unlimited | 0 | 0 |
| Deterministic Objective | If checked, suppresses noise modeling | off | off |

The recommended column reports the settings that maximized mean GRM and Fᵦ across the 231 sensitivity runs described in Section 2.2.2 of the paper (`BHO_SensitivityAnalysis.m`).

#### Validation Parameters

| GUI Field | Description | Default |
|---|---|---|
| Tolerance | Spatial buffer (coordinate units) for Fᵦ matching; leave **empty** for automatic estimation | (auto) |
| F-beta (β) | Weighting of recall vs precision; β < 1 emphasizes precision | 0.2 |

The **Tolerance** field applies to Fᵦ matching only. The R_P21 matching distance is fixed at 70 pixels inside `calculateStructuralMetrics.m` and is not exposed in the GUI.

### Step 6 — Run BTO

Click **Lineaments (Auto)** to launch the optimization loop. The pipeline:

1. Clears previous state.
2. Builds the `optimizableVariable` search space from GUI ranges.
3. Wraps `Lineaments_Auto7` as an anonymous objective for `bayesopt`.
4. Runs `bayesopt` with the configured settings.
5. Writes optimized values back to the GUI and saves `BTO_Progress.xlsx`.

A live figure updates each iteration showing the current curvilinear pattern map (left) and a GRM vs Fᵦ scatter plot (right).

### Step 7 — Inspect Results

Use the **Display** dropdown on the Lineaments tab to view:

| Display Option | Description |
|---|---|
| Target Lines | Ground-truth curvilinear patterns (vector overlay) |
| CWT (Auto) | kCWT coefficient maps from the best iteration |
| Selected (Auto) | Feature-selected kCWT maps |
| All Lineaments (Auto) | Per-feature curvilinear pattern overlays (tabbed gallery) |
| Step Filtering (Auto) | Average step-filter width visualization |
| Lineament Densities (Auto) | Convolution-based density heat map |
| Lineaments (Auto) | Final stacked binary curvilinear pattern map |
| Deep Lineaments (Auto) | Subset from features with above-average filter width |
| Shallow Lineaments (Auto) | Subset from features with below-average filter width |
| Fbeta Score (Auto) | Overlay with Fᵦ, precision, recall, structural metrics |
| GRM vs Ground Truth | Correlation scatter (GRM vs Fᵦ, R_P21, HS, JSD) |
| BTO Progress (Auto) | Objective trace, parameter evolution, convergence plots |
| BTO Best So Far Lineaments (Auto) | Iteration-by-iteration gallery of best curvilinear pattern maps |

### Step 8 — Export

Switch to the **Export** tab. Select a data product via the radio buttons:

| Radio Button | Outputs |
|---|---|
| Point Data | Original interpolated data |
| CWT (Auto) | CWT coefficient grids |
| Selected (Auto) | Feature-selected grids |
| Densities (Auto) | Lineament density map |
| All Lines (Auto) | All curvilinear patterns (vector) |
| Deep Lines (Auto) | Deep curvilinear patterns (vector) |
| Shallow Lines (Auto) | Shallow curvilinear patterns (vector) |

Then click:
- **CSV** — Exports as `[X, Y, Value]` columns (raster) or `[ID, X, Y]` (vector).
- **GRD** — Exports as GeoTIFF (`.tif`) with geographic referencing.
- For vector outputs, a **Shapefile (.shp)** dialog appears automatically.

---

## BTO Objective Function Pipeline

`Lineaments_Auto7.m` implements the four-stage pipeline that `bayesopt` evaluates at each iteration. This corresponds to the methodology described in Sections 1.1–1.4 of the paper.

```
Input Image(s) + Hyperparameters h
        │
        ▼
┌──────────────────────────────────────┐
│  Stage 1: Poisson kCWT               │
│  Frequency-domain implementation     │
│  of Poisson kernel-based CWT         │
│  → multidirectional coefficient      │
│    maps                              │
│  (Section 1.1)                       │
└──────────────┬───────────────────────┘
               │
               ▼
┌──────────────────────────────────────┐
│  Stage 2: kCWT Coefficient Selection │
│  and Fusion                          │
│  Saliency + dissimilarity-based      │
│  selection of compact, non-redundant │
│  feature subset                      │
│  (Section 1.2 / Algorithm 1)         │
└──────────────┬───────────────────────┘
               │
               ▼
┌──────────────────────────────────────┐
│  Stage 3: Curvilinear Pattern        │
│  Extraction                          │
│  Per selected feature:               │
│    FilterB → step_Filtering →        │
│    getFaultDetection →               │
│    component strength filter         │
│  Fuse maps: intersection during      │
│  BO scoring, union for display       │
└──────────────┬───────────────────────┘
               │
               ▼
┌──────────────────────────────────────┐
│  Stage 4: Topological Representation │
│  and Graph Representativeness        │
│  Graph construction → GRM scoring    │
│  (Section 1.3, Eq. 18)               │
│  BTO minimizes F(h) = −GRM           │
│  (Section 1.4, Eq. 20)               │
└──────────────────────────────────────┘
```

### Graph Representativeness Metric (GRM)

GRM scores an extracted network without ground truth (Eq. 18 in the paper):

\[ \text{GRM} = \sum_{P=1}^{N} w_P \cdot \left(1 - \min_{\substack{Q=1 \\ Q \neq P}}^{N} \text{sim}(S_P, S_Q)\right) \]

where **w_P** is the length of segment *P* in pixels and **sim(S_P, S_Q) = 1 / (1 + d_PQ)** is the spatial similarity between segments *P* and *Q*, with *d_PQ* their separation (Eq. 19). Maximizing GRM favors maps that retain a large total length of continuous traces, since each segment contributes in proportion to the pixels it occupies; broken or noisy maps collapse into many short fragments and score low. A multi-term similarity adding orientation and overlap components is identified as future work in the paper and is **not** used here.

### BTO Objective

BTO explores the hyperparameter space **χ** for the configuration **h\*** that maximizes GRM, equivalently minimizing (Eq. 20 in the paper):

\[ \mathbf{h}^* = \arg\min_{\mathbf{h} \in \chi} F(\mathbf{h}), \quad \text{where } F(\mathbf{h}) = -\text{GRM} \]

Since `bayesopt` minimizes, a more negative objective corresponds to higher GRM (better structural representativeness of extracted curvilinear patterns). The raw, unnormalized score is used: every iteration within a run evaluates the same spatial grid, so the pixel-count score is sufficient for ranking hyperparameter configurations.

If ground truth is loaded, Fᵦ and structural metrics are computed at each iteration for tracking — they do **not** influence the objective.

---

## Parameter Reference

### Fixed Parameters

| Parameter | Variable | Paper symbol | Description |
|---|---|---|---|
| Number of Scales | `N_Scales` | *N*ₐ | How many kCWT scales to compute. Fixed at 1 in the paper so that discrete scale jumps do not mix sources from different depths. |
| Scale Dilation Factor | `SDF` | δ | Step between scales: `Scales = 1 : SDF : N_Scales×SDF`. |
| Line Resolution | `Line_Res` | *R* | Image resize factor applied before curvilinear pattern extraction. `1` = original resolution; the paper uses `2`. All pixel-based tolerances and filter widths are expressed on this resized grid. |
| CWT Reduction | `DR` | *d* | Number of kCWT feature maps retained after coefficient selection and fusion. |
| Image Filter Ratio | `sigma2` | σ_I | Gaussian blur sigma applied to each feature before step filtering. |

### Optimizable Hyperparameters

| Variable | Type | Range used in the paper | Role in Pipeline |
|---|---|---|---|
| `WSFR` | real | [0, 1] | kCWT support scaling factor — controls spatial extent of the mother wavelet |
| `NAngles` | integer | [12, 72] | Number of kCWT angular orientations — higher values detect more orientations but increase computation |
| `orderx` | categorical | {0, 1} | Derivative order — selects between horizontal (1,0) and vertical (0,1) emphasis |
| `Z` | real | [1, 5] | Poisson upward continuation — controls the depth sensitivity of the PDMW |
| `w` | real | [0.75·w_c, 2.5·w_c] | Average step-filter window size — larger windows detect broader/deeper features |
| `VSFW` | real | [0, 1] | Variance-based window scaling — adapts filter width per feature based on local signal variance |
| `NSFA` | integer | [12, 72] | Number of step-filter angles — directional coverage of the enhancement filter |
| `SLC` | integer | [35, 85] | Global confidence threshold — percentile cutoff for retaining connected components |
| `FDissimilarity` | real | [0, 1] | Dissimilarity weight in feature selection — higher values prefer less-correlated features |
| `FSaliency` | real | [0, 100] | Saliency percentile — edge-strength threshold during feature selection |

### Bayesian Optimization Settings

| Setting | MATLAB Parameter | Notes |
|---|---|---|
| Max Evaluations | `MaxObjectiveEvaluations` | Total iterations including seed points |
| Acquisition Function | `AcquisitionFunctionName` | Options: `expected-improvement-per-second-plus`, `expected-improvement`, `expected-improvement-plus`, `expected-improvement-per-second`, `lower-confidence-bound`, `probability-of-improvement` |
| Exploration Ratio | `ExplorationRatio` | 0 = pure exploitation, 1 = pure exploration |
| Seed Points | `NumSeedPoints` | Random initial evaluations before GP-guided search |
| GP Active Set | `GPActiveSetSize` | Subset size for scalable GP approximation |
| Max Time | `MaxTime` | Wall-clock limit in seconds (GUI input is in minutes) |
| Deterministic | `IsObjectiveDeterministic` | Set `true` if objective has no stochastic component |

### Validation Parameters

| Parameter | Description |
|---|---|
| **Tolerance** | Spatial buffer distance (in coordinate units, e.g., arc-seconds) for matching detected curvilinear patterns to ground truth. If left empty, `estimateTolerancePreBTO` runs a pilot detection to find a tolerance yielding Fᵦ in [0.4, 0.6]. |
| **F-beta (β)** | Controls the precision–recall trade-off. `β < 1` emphasizes precision (fewer false positives); `β > 1` emphasizes recall (fewer missed features). Default `0.2` strongly favors precision. |

---

## Display Options

All display options are accessible from the **Lineaments** tab dropdown + **Plot** button:

| Option | What It Shows |
|---|---|
| **Target Lines** | Ground-truth curvilinear patterns as vector lines over the study area with rose diagram |
| **CWT (Auto)** | Tabbed gallery of kCWT coefficient maps from the best BTO iteration |
| **Selected (Auto)** | Feature-selected kCWT maps after saliency/dissimilarity reduction |
| **All Lineaments (Auto)** | Per-feature colored curvilinear pattern overlays in a tabbed figure, each with density map, raster, and ground-truth toggles |
| **Step Filtering (Auto)** | Average adaptive step-filter width used across features |
| **Lineament Densities (Auto)** | Convolution-based density heat map of all stacked curvilinear patterns |
| **Lineaments (Auto)** | Final stacked binary curvilinear pattern map (all features combined) |
| **Deep Lineaments (Auto)** | Patterns from features whose step-filter width exceeds the average (deeper structures) |
| **Shallow Lineaments (Auto)** | Patterns from features with below-average filter width (shallower structures) |
| **Fbeta Score (Auto)** | Detected vs ground-truth overlay with Fᵦ, precision, recall, TP, FP, FN, and structural metrics |
| **GRM vs Ground Truth** | Four scatter panels: GRM vs Fᵦ, R_P21, HS, JSD, with Pearson correlation coefficients (*r*) |
| **BTO Progress (Auto)** | Objective trace, minimum trace, per-variable evolution, and convergence diagnostics |
| **BTO Best So Far Lineaments (Auto)** | Iteration gallery showing the best curvilinear pattern map at each improvement step |

---

## Export Formats

| Format | Method | Content |
|---|---|---|
| **CSV** | `Custom_csvwrite` | Three-column `[X, Y, Value]` for raster products; `[ID, X, Y]` for vector lineaments |
| **GeoTIFF** | `geotiff` | Georeferenced `.tif` with coordinate bounds from the study area |
| **Shapefile** | `shapewrite` (MATLAB) | ESRI `.shp` + `.shx` + `.dbf` for vector curvilinear pattern geometries |
| **Excel** | Auto-generated | `BTO_Progress.xlsx` — full parameter trace, objective values, and convergence data per iteration |

---

## File Descriptions

### Core Application

| File | Description |
|---|---|
| `BWTO2D.mlapp` | Main GUI application (App Designer). Contains all UI layout, callbacks, and orchestration logic. |
| `Lineaments_Auto7.m` | BTO objective function. Implements the four-stage pipeline (Poisson kCWT → kCWT coefficient selection and fusion → curvilinear pattern extraction → GRM scoring). Called by `bayesopt` at each iteration. |

### Reproducibility Scripts

| File | Description |
|---|---|
| `BHO_SensitivityAnalysis.m` | Runs the one-at-a-time sensitivity analysis over exploration ratio, seed points, GP active set size, and acquisition function (231 runs, Figure 10 and Section 2.2.2). |
| `FeatureReductionComparison.m` | Ablation study of the dimensionality-reduction stage, distributed in `Feature Reduction Comparisons.zip`. Compares the proposed feature selection against PCA and the unreduced kCWT stack under identical hyperparameters, and regenerates Figures 18 and 19. Calls the pipeline functions unchanged. |
| `add_noise_to_xyz.m` | Adds zero-mean Gaussian noise scaled to a percentage of the data range, used to build the 10% and 25% noise datasets of Section 2.2.3. |

### Poisson kCWT Engine

| File | Description |
|---|---|
| `cwt2d.m` | Core 2D CWT computation engine — frequency-domain implementation (derived from YAWTB). |
| `cwt2D_DerGus.m` | kCWT with derivative-of-Gaussian kernel — used by the manual CWT path. Includes wavelet visualization. |
| `cwt2D_DerPoisson.m` | kCWT with Poisson Derivatives Mother Wavelet (PDMW) — used by the manual CWT path. Includes wavelet visualization. |
| `cwt2D_DerPoisson_BTO.m` | kCWT with PDMW — streamlined for BTO (no visualization overhead). This is the primary kernel used in the paper (Section 1.1). |
| `dergauss2d.m` | Generates the 2D derivative-of-Gaussian wavelet kernel in the frequency domain. |
| `derpoisson2d.m` | Generates the 2D Poisson Derivatives Mother Wavelet (PDMW) kernel in the frequency domain. |

### Image Processing

| File | Description |
|---|---|
| `CWTFeatureSelection.m` | kCWT Coefficient Selection and Fusion — iterative selection balancing feature saliency and dissimilarity (Section 1.2 / Algorithm 1 in the paper). |
| `step_Filtering.m` | Directional step filtering — multi-angle, multi-bandwidth zero-mean step filter for curvilinear pattern enhancement (Panagiotakis and Kokinou, 2015). |
| `getFaultDetection.m` | Adaptive thresholding + morphological cleanup for binary curvilinear pattern detection. |
| `FilterB.m` | Gaussian image filter; applies `imgaussfilt` if `sigma > 0`, otherwise passes through. |
| `nDstrb2D.m` | Normalizes a 2D matrix's statistical distribution (handles NaN). |
| `plotNewColorMap0_.m` | Renders colored curvilinear pattern overlays for visualization (used in BTO path). |

### Validation and Metrics

| File | Description |
|---|---|
| `calculateFbetaScore.m` | Computes Fᵦ score with configurable tolerance buffer, bidirectional matching, and length weighting. |
| `calculateStructuralMetrics.m` | Computes the fracture intensity ratio R_P21, orientation distributions, Hellinger similarity (HS), and Jensen–Shannon divergence (JSD). R_P21 is the fraction of ground-truth skeleton length falling within a 70-pixel dilation of the detected skeleton (Eq. 23), a recall-type measure, not the ratio of the two P₂₁ intensities. |
| `calculateAngularSimilarity.m` | Angular histogram comparison between detected and target orientations (optional). |
| `estimateTolerancePreBTO.m` | Runs a quick pilot detection to auto-estimate a tolerance that yields Fᵦ ∈ [0.4, 0.6]. |
| `plotGRMvsGroundTruth.m` | Generates correlation scatter plots (GRM vs each ground-truth metric). |
| `createRoseDiagrams.m` | Extracts curvilinear pattern orientations from a binary raster image. |
| `createRoseDiagramsFromVectors.m` | Extracts curvilinear pattern orientations from a vector shape structure. |

### Result Management and Display

| File | Description |
|---|---|
| `useBHO_BestOutputs.m` | Reconstructs display globals (`plotH1_rgb`, density maps, deep/shallow splits) from stored best BTO iteration outputs. |
| `convertRasterToShape.m` | Converts a binary raster curvilinear pattern image into a MATLAB shape structure. |
| `Surf1.m` | Tabbed gallery viewer for a stack of 2D feature maps. |
| `Surf1_Enhanced.m` | Enhanced multi-panel viewer for loaded point data. |
| `Custom_csvwrite.m` | Exports a matrix to CSV with a user-specified numbered suffix. |
| `geotiff.m` | Writes a matrix as a georeferenced GeoTIFF file. |

### UI Toggle Callbacks

| File | Description |
|---|---|
| `toggleGTVisibilityCallback.m` | Toggles ground-truth overlay on main figures. |
| `toggleGTVisibilityForTab.m` | Toggles ground-truth overlay per tab in the best-so-far gallery. |
| `toggleVectorLineamentsVisibility.m` | Toggles vector lineament visibility per tab. |
| `toggleImageVisibility.m` | Toggles density map or raster lineament visibility per tab. |

### YAWTB Library (Bundled)

Files from the [Yet Another Wavelet Toolbox](http://www.fyma.ucl.ac.be/projects/yawtb/) required by `cwt2d.m`:

`getopts.m`, `getyawtbprefs.m`, `yashow.m`, `yashow_cwt2d.m`, `yashow_matrix.m`, `yashow_spheric.m`, `yashow_sphvf.m`, `yashow_timeseq.m`, `yashow_volume.m`, `yapuls.m`, `yapuls2.m`, `yapbar.m`, `yamax.m`, `yahelp.m`, `yadiro.m`, `yastrfind.m`, `yawopts.m`, `yathresh.m`, `vect.m`, `sphgrid.m`, `list_elem.m`

---

## Data Directories

| Directory | Contents |
|---|---|
| `Data 1 Synthetic/Without Noise/` | Synthetic benchmark results — Widgiemooltha Dome model, no noise (Section 2.2.1) |
| `Data 1 Synthetic/Guassian Noise/Exp_5p/` | Noise robustness test — 5% Gaussian noise (additional run, not reported in the paper) |
| `Data 1 Synthetic/Guassian Noise/Exp_10p/` | Noise robustness test — 10% Gaussian noise (Section 2.2.3) |
| `Data 1 Synthetic/Guassian Noise/Exp_25p/` | Noise robustness test — 25% Gaussian noise (Section 2.2.3) |
| `Data 1 Synthetic/Sensitivty/` | BTO sensitivity analysis — exploration ratio, seed points, GP active set, acquisition function (Section 2.2.2) |
| `Data 2 Blake River/Run 1–4/` | Real-world case study — Blake River Group (BRG), Abitibi Greenstone Belt, Quebec (Section 3) |

These directories are distributed as the archives `Data 1 Synthetic.zip` and `Data 2 Blake River.zip` in the repository root; unzip them in place to obtain the paths above.

The synthetic data is derived from the **Widgiemooltha Dome** geological model (Yilgarn Craton, Western Australia) using the Noddy 3D Modelling System. The real-world data consists of aeromagnetic datasets from the volcanogenic massive sulphide-endowed **Blake River Group** area.

Each run directory typically contains: `Data_1.csv` (input), `Selected_Auto_*.csv` (selected features), `All/Deep/Shallow Lineaments (Auto).csv` (results), `Lineaments Densities Auto_*.csv` (density map), and a `Log*.txt` file.

---

## Troubleshooting

### "Undefined function cwt2D_DerGus_BTO"

The Gaussian-kernel BTO path calls `cwt2D_DerGus_BTO`, which must be on your MATLAB path. The paper uses the **Poisson kCWT** exclusively (`derpoisson2d`), so this function is not needed for the default pipeline. To use Gaussian kCWT in BTO mode, create `cwt2D_DerGus_BTO.m` analogous to `cwt2D_DerPoisson_BTO.m` (stripped of visualization code from `cwt2D_DerGus.m`).

### Out of memory during kCWT

- Reduce **Number of Scales** or **Number of CWT Angles**.
- Increase **Spacing (Arc-Seconds)** to coarsen the grid.
- Reduce **CWT Reduction To** to select fewer features.

### BTO runs slowly

- Reduce **BTO Iterations**.
- Set a **Max Time** limit (in minutes).
- Use **EI per Second Plus** (default) as the acquisition function — it accounts for evaluation time.
- Reduce grid resolution via larger arc-second spacing.

### Fᵦ is always 0 or NaN

- Verify that ground truth is loaded (click **Lineaments Shape File** and check the MATLAB console).
- Check that the ground-truth spacing matches the study-area grid.
- Try increasing the **Tolerance** value or leaving it empty for automatic estimation.

### Figures do not appear

- Close any stale figures: `close all` in the MATLAB Command Window.
- Ensure the Lineaments tab display dropdown has a valid selection before clicking **Plot**.

---

## Citation

If you use BWTO2D in your research, please cite:

> Abbassi, B., Cheng, L.-Z., Legault, M. (2026). Bayesian Wavelet Topology Optimization for Curvilinear Pattern Recognition in Magnetic Data. *Computers & Geosciences*. <!-- doi: https://doi.org/TBD -->

```bibtex
@article{Abbassi2026BWTO2D,
  title     = {Bayesian Wavelet Topology Optimization for Curvilinear Pattern Recognition in Magnetic Data},
  author    = {Abbassi, Bahman and Cheng, Li-Zhen and Legault, Marc},
  journal   = {Computers \& Geosciences},
  year      = {2026},
  publisher = {Elsevier}
}
```

Source code: [https://github.com/bahmanabbassi/BWTO2D](https://github.com/bahmanabbassi/BWTO2D)

### Related work

> Abbassi, B., Cheng, L.-Z. (2025). Curvilinear lineament extraction: Bayesian optimization of principal component wavelet analysis and hysteresis thresholding. *Computers & Geosciences* 194, 105768. https://doi.org/10.1016/j.cageo.2024.105768

---

## License

Copyright (c) Bahman Abbassi, Li-Zhen Cheng, Marc Legault — Université du Québec en Abitibi-Témiscamingue (UQAT).

This software is distributed under the **GNU General Public License v3.0** (see the `LICENSE` file for the full terms).

Lead Developer: Bahman Abbassi


### Dependencies and Acknowledgments

This project uses components of the **Yet Another Wavelet Toolbox (YAWTB)**, copyright (C) 2001–2002 by the YAWTB team. The original license headers have been preserved in the relevant `.m` files, in accordance with the YAWTB licensing conditions. [YAWTB GitHub](https://github.com/jacquesdurden/yawtb).

BWTO2D also uses:

* A **geological fault detection** code (hysteresis thresholding / step filtering) developed by Costas Panagiotakis: *Costas Panagiotakis (2026). Detection of Geological Faults, MATLAB Central File Exchange. Retrieved March 16, 2026.* [Link](https://www.mathworks.com/matlabcentral/fileexchange/64693-detection-of-geological-faults)

*Note: In addition to citing BWTO2D itself, we recommend citing these packages when the corresponding modules are explicitly used.*
