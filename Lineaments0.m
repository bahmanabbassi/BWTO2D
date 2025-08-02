function plotH = Lineaments0(Grid, SLC)
%Lineaments0 Detects and visualizes linear features (lineaments) in a 2D data grid.
%
% Description:
%   This function implements a workflow to identify and visualize lineaments from a
%   2D data grid. The process involves resizing the data to a standard resolution,
%   applying a directional step-filter to enhance linear features, and then
%   detecting an initial set of lineaments. These initial detections are filtered
%   based on their "strength" (derived from the step-filter response) to keep
%   only the most prominent components, controlled by the 'SLC' parameter. The
%   final output is an RGB image that visualizes these detected lineaments.
%
%   The function's behavior is controlled by several global variables that must be
%   set prior to calling it.
%
% Developer:
%   Bahman Abbassi (bahman.abbassi@uqat.ca)
%   University of Quebec (UQAT)
%   Rouyn-Noranda, QC, Canada
%   *Annotation Date: August 1, 2025*
%
% Inputs:
%   Grid (matrix) - A 2D matrix representing the data for lineament detection.
%   SLC (scalar)  - A percentile value (0-100) controlling the sensitivity. Higher
%                   values result in more lineaments being kept after strength filtering.
%
% Output:
%   plotH (image) - An RGB image visualizing the detected lineaments.
%
% Example:
%   % --- Setup for Example ---
%   % This function relies on global variables. Define them first.
%   global XI YI Line_Res grm2 GSF_Angles;
%   [XI, YI] = meshgrid(1:250, 1:250);
%   Line_Res = [150 150]; % Target resolution for processing
%   grm2 = 10;            % Step filter width
%   GSF_Angles = 16;      % Number of filter angles
%
%   % Create sample data and call the function 🛰️
%   sample_grid = peaks(250);
%   sensitivity_level = 80; % Keep the top 20% strongest components
%   Plot_Image = Lineaments0(sample_grid, sensitivity_level);
%
%   % Display the result
%   figure;
%   imshow(Plot_Image);
%   title('Final Detected Lineaments Visualization');


    % --- Declare Global Variables ---
    % These parameters are configured externally and control the function's operation.
    global YI
    global XI
    global X_Pix
    global Y_Pix
    global GSF_Angles
    global grm2
    global Line_Res
    global XI_line_size
    global XI_line
    global YI_line

    % --- Step 1: Data and Coordinate Preparation ---
    % Resize the input grid and coordinate matrices to a consistent resolution.
    Grid = imresize(Grid, Line_Res);
    XI_line = imresize(XI, Line_Res);
    YI_line = imresize(YI, Line_Res);
    XI_line_size = size(XI_line);

    % Convert grid to double precision for numerical processing.
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
    % Generate an initial binary map of potential faults from the filtered result.
    [BW, ~] = getFaultDetection(Y, data1, SLC);

    % --- Step 4: Component-Level Strength Filtering ---
    % This section removes weak or noisy lineaments.
    CC = bwconncomp(BW); % Find all connected components (individual lineaments).
    % Extract pixel values from the step-filtered image 'Y' for each component.
    stats = regionprops(CC, Y, 'PixelValues');
    % Calculate the "strength" of each component by summing its pixel values.
    strengths = cellfun(@sum, {stats.PixelValues});

    % Keep only the top percentage of the strongest lineaments, controlled by SLC.
    pct = 100 - SLC; % Calculate the percentile threshold.
    th  = prctile(strengths, pct); % Find the strength value at that percentile.
    keep = find(strengths >= th); % Find indices of components above the threshold.
    BW  = ismember(labelmatrix(CC), keep); % Create a new binary mask with only the strong components.

    % --- Step 5: Feature Extraction and Classification ---
    C = bwconncomp(BW); % Re-identify connected components from the cleaned binary image.
    L = C.NumObjects; % Get the number of remaining lineaments.
    Map = labelmatrix(C); % Create a label matrix from the final components.

    % Compute feature values for each remaining component.
    RS_Y = regionprops(C, Y, 'PixelValues');
    feat = zeros(L, 2); % Initialize a feature matrix.
    for i = 1:L
        vec = RS_Y(i).PixelValues;
        feat(i, 1) = sum(vec); % Feature 1 is the component strength.
    end
    feat(:, 2) = feat(:, 1); % Feature 2 is a copy (can be replaced with another metric).

    % Determine a threshold to classify the "strongest" faults for visualization.
    sv = sort(feat(:, 1), 'descend'); % Sort lineaments by strength.
    % Calculate the number of faults to classify as "strong".
    NUMBER_OF_Faults = round(0.01 * SLC * length(sv));
    if NUMBER_OF_Faults < 1 && ~isempty(sv)
        NUMBER_OF_Faults = 1; % Ensure at least one is selected.
    end
    Thresh = sv(NUMBER_OF_Faults); % Set the threshold to the strength of the Nth fault.

    % Find the labels of regions classified as strong faults.
    Labels = find(feat(:, 1) >= Thresh);

    % --- Step 6: Create Final Visualization ---
    % Generate a color-mapped RGB image. The step-filtered image 'Y' is used as
    % the background to provide context for the detected lineaments.
    [plotH] = plotNewColorMap0(Y, Y, Map, Labels, 'Detected Lineaments', YI_line(:), XI_line(:));
end