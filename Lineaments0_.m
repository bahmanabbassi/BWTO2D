function [plotH_, Number_OF_Faults, BW1] = Lineaments0_(Grid, SLC)
%Lineaments0_ Detects, filters, and classifies linear features (faults) in a 2D data grid.
%
% Description:
%   This function executes a workflow to identify lineaments from a 2D grid.
%   The process includes:
%   1. Resizing the input grid to a standard resolution.
%   2. Applying a step-filter to enhance linear features across multiple orientations.
%   3. Performing an initial fault detection based on the filtered output.
%   4. Filtering the detected lineaments based on their "strength" to remove noise
%      and weak features. The 'SLC' parameter controls this threshold.
%   5. Classifying the remaining lineaments to identify a subset of "strong faults".
%   6. Generating a color-mapped RGB image to visualize the results.
%
%   The function relies on several global variables for processing parameters like
%   resolution, filter settings, and plot coordinates.
%
% Developer:
%   Bahman Abbassi (bahman.abbassi@uqat.ca)
%   University of Quebec (UQAT)
%   Rouyn-Noranda, QC, Canada
%   *Annotation Date: August 1, 2025*
%
% Inputs:
%   Grid - A 2D matrix representing the data (e.g., geophysical, satellite).
%   SLC  - A scalar percentile value (0-100) that controls the sensitivity of
%          lineament detection and classification. Higher values keep more features.
%
% Outputs:
%   plotH_           - An RGB image visualizing the detected lineaments.
%   Number_OF_Faults - The total count of features classified as "strong faults".
%   BW1              - The final binary image (mask) of detected lineaments after filtering.
%
% Example:
%   % --- Setup for Example ---
%   % This function uses global variables, so we must define them first.
%   global XI YI Line_Res grm2 GSF_Angles;
%   [XI, YI] = meshgrid(1:200, 1:200);
%   Line_Res = [128 128]; % Target resolution
%   grm2 = 8;             % Step filter width
%   GSF_Angles = 12;      % Number of filter angles
%
%   % Create sample data and call the function
%   sample_grid = peaks(200);
%   sensitivity_level = 75; % Keep top 25% of components
%   [Plot_Image, Strong_Fault_Count, Fault_Mask] = Lineaments0_(sample_grid, sensitivity_level);
%
%   % Display the results
%   figure;
%   imshow(Plot_Image);
%   title(['Detected Lineaments (', num2str(Strong_Fault_Count), ' Strong)']);


    % --- Declare Global Variables ---
    % These parameters are configured outside this function and control its behavior.
    global YI
    global XI
    global sigma2 % Note: sigma2 is declared but not used in this function.
    global X_Pix  % Note: X_Pix is declared but not used in this function.
    global Y_Pix  % Note: Y_Pix is declared but not used in this function.
    global Line_Res
    global grm2
    global GSF_Angles
    global XI_line
    global YI_line

    % --- Step 1: Data and Coordinate Preparation ---
    % Resize the input grid and coordinate matrices to a consistent resolution.
    Grid = imresize(Grid, Line_Res);
    XI_line = imresize(XI, Line_Res);
    YI_line = imresize(YI, Line_Res);

    % Convert grid to double precision for numerical stability in filtering.
    data1 = im2double(Grid);

    % --- Step 2: Step Filtering to Enhance Lineaments ---
    % Initialize parameters for the step-filtering process.
    I = data1;                    % Use the prepared grid as input.
    grammes = grm2;               % Set filter width from global parameter.
    Step_Filter_nAng = GSF_Angles;% Set number of filter angles from global.
    Step_Filter_step = pi / Step_Filter_nAng; % Calculate angle increment.
    Step_Filter_DTheta = 0:Step_Filter_step:pi - Step_Filter_step; % Create angle array.
    Step_Filter_Dbw = [2 4 6];    % Define a set of fixed widths for the filter.

    % Apply the step-filter to enhance linear features in the image.
    [Y] = step_Filtering(I, grammes, Step_Filter_DTheta, Step_Filter_Dbw, 2);

    % --- Step 3: Initial Fault Detection ---
    % Generate an initial binary map of potential faults from the filtered data.
    [BW1, ~] = getFaultDetection(Y, data1, SLC);

    % --- Step 4: Component-Level Strength Filtering ---
    % This section removes weak or noisy lineaments based on their strength.
    CC = bwconncomp(BW1); % Find all connected components (individual lineaments).
    % Extract the underlying pixel values from the step-filtered image (Y) for each component.
    stats = regionprops(CC, Y, 'PixelValues');
    % Calculate the "strength" of each component by summing its pixel values.
    strengths = cellfun(@sum, {stats.PixelValues});

    % Keep only the top percentage of strongest lineaments, controlled by SLC.
    pct = 100 - SLC; % Calculate the percentile threshold.
    th  = prctile(strengths, pct); % Find the strength value at that percentile.
    keep = find(strengths >= th); % Find the indices of components above the threshold.
    BW1  = ismember(labelmatrix(CC), keep); % Create a new binary mask with only the strong components.

    % --- Step 5: Feature Extraction from Filtered Lineaments ---
    C = bwconncomp(BW1); % Re-identify connected components from the cleaned binary image.
    L = C.NumObjects;    % Get the number of remaining lineaments.

    % Compute feature values for each remaining component.
    RS_Y = regionprops(C, Y, 'PixelValues'); % Get pixel values from Y again for the cleaned set.
    feat = zeros(L, 2); % Initialize a feature matrix.
    for i = 1:L
        vec = RS_Y(i).PixelValues;
        feat(i, 1) = sum(vec); % Feature 1 is the sum of pixel values (strength).
    end
    feat(:, 2) = feat(:, 1); % Feature 2 is a copy (can be replaced with another metric).

    % --- Step 6: Classify "Strong Faults" ---
    % Determine a threshold to classify the most significant faults.
    sv = sort(feat(:, 1), 'descend'); % Sort lineaments by strength.
    % Calculate the number of faults to classify as "strong" based on a small percentage of SLC.
    Number_OF_Faults = round(0.01* SLC * length(sv));
    if Number_OF_Faults < 1 && ~isempty(sv) % Ensure at least one fault is selected if possible
        Number_OF_Faults = 1;
    end
    Thresh = sv(Number_OF_Faults); % Set the threshold to the strength of the Nth fault.

    % Find the labels of the regions classified as strong faults.
    Labels = find(feat(:, 1) >= Thresh);

    % --- Step 7: Create Visualization ---
    % Generate a color-mapped RGB image of the detected lineaments,
    % highlighting the strong faults.
    Map = labelmatrix(C); % Create a label matrix from the final components.
    [plotH_] = plotNewColorMap0_(Y, zeros(size(data1)), Map, Labels, ...
                             'Detected Lineaments', YI_line(:), XI_line(:));
end