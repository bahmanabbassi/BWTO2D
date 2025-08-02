function [L_val, N_CWT] = Lineaments_Auto7(CWT_Input_Features_Matrix, All_Data_Matrix, WSFR, N_Angles, orderx, Z, w, VSFW, NSFA, SLC, FDissimilarity, FSaliency, N_Scales, Line_Res, DR, SDF, sigma2)
warning('off','all');
%--------------------------------------------------------------------------
% FUNCTION: Lineaments_Auto7
%
% PURPOSE:
%   Automates the detection of geological lineaments from geophysical or
%   remote sensing data using a multi-stage process. The function is
%   designed to be an objective function for Bayesian optimization, where
%   it evaluates a set of hyperparameters and returns a score (`L_val`)
%   representing the quality of the detected lineaments.
%
%   The process involves four main stages:
%   1.  **CWT Feature Generation**: Applies a 2D Continuous Wavelet Transform
%       to generate feature maps at multiple scales and angles.
%   2.  **Feature Selection**: Reduces the dimensionality of the CWT features
%       using a custom selection method based on saliency and dissimilarity.
%   3.  **Lineament Extraction**: Identifies lineament segments from the
%       selected feature maps using filtering and morphological operations.
%   4.  **Graph Representativeness**: Constructs a graph from the extracted
%       lineaments and computes a "representativeness" score, which serves
%       as the final objective value for optimization.
%
% AUTHOR:
%   Bahman Abbassi
%   University of Quebec (UQAT)
%   Email: bahman.abbassi@uqat.ca
%   GitHub: https://github.com/bahmanabbassi
%   MathWorks: https://www.mathworks.com/matlabcentral/profile/authors/31336389
%
% INPUTS:
%   CWT_Input_Features_Matrix (double): A 2D or 3D matrix of input data for CWT.
%   All_Data_Matrix (double): A 3D matrix of additional data to be combined
%                             with CWT features before lineament extraction.
%   WSFR (double):      Wavelet support scaling factor for CWT.
%   N_Angles (integer): Number of CWT angles to compute.
%   orderx (integer):   Derivative order in the x-direction for the CWT kernel.
%   Z (double):         Parameter for the Poisson wavelet.
%   w (double):         Base window size for the step-filter.
%   VSFW (double):      Variance-based scaling factor for the step-filter window.
%   NSFA (integer):     Number of angles for the step-filter.
%   SLC (double):       Strength level cutoff for lineament component filtering (%).
%   FDissimilarity (double): Dissimilarity weight for CWT feature selection.
%   FSaliency (double): Saliency percentile for CWT feature selection.
%   N_Scales (integer): Number of CWT scales to compute.
%   Line_Res (double):  Resize factor for the feature maps before extraction.
%   DR (integer):       Dimensionality reduction target (number of features).
%   SDF (double):       Scale step for CWT.
%   sigma2 (double):    Sigma for the Gaussian filter (`FilterB`).
%
% OUTPUTS:
%   L_val (double):     The final objective function value (Graph Representativeness).
%                       Lower is better for Bayesian optimization.
%   N_CWT (integer):    Total number of CWT features generated before selection.
%
% EXAMPLE:
%   % This function is typically called by a Bayesian optimizer. A direct
%   % call would require setting up numerous global variables and inputs.
%   %
%   % --- Setup (Illustrative) ---
%   % global XI YI X_Pix Y_Pix in pop_choice3 YI_Real_Ratio;
%   % load('input_data.mat'); % Load CWT_Input, All_Data, globals, etc.
%   %
%   % % Define a sample hyperparameter set
%   % params.WSFR = 1; params.N_Angles = 8; params.orderx = 1;
%   % params.Z = 2; params.w = 5; params.VSFW = 0.5;
%   % params.NSFA = 8; params.SLC = 80; params.FDissimilarity = 0.5;
%   % params.FSaliency = 75; params.N_Scales = 10; params.Line_Res = 1;
%   % params.DR = 10; params.SDF = 1.5; params.sigma2 = 1;
%   %
%   % [objectiveValue, numFeatures] = Lineaments_Auto7(...
%   %     CWT_Input_Features_Matrix, All_Data_Matrix, params.WSFR, ...
%   %     params.N_Angles, params.orderx, params.Z, params.w, params.VSFW, ...
%   %     params.NSFA, params.SLC, params.FDissimilarity, params.FSaliency, ...
%   %     params.N_Scales, params.Line_Res, params.DR, params.SDF, params.sigma2);
%   %
%   % fprintf('Objective Value (L_val): %.4f\n', objectiveValue);
%--------------------------------------------------------------------------



% Declare global variables used within this function or by called functions
global XI YI X_Pix Y_Pix in; % Essential for geometric context and masking
global pop_choice3;         % Determines the CWT kernel
global YI_Real_Ratio;       % Used for cropping operations
global XI_line_Auto YI_line_Auto; % For resized coordinate grids


if orderx == 1
    ordery = 0;
elseif orderx == 0
    ordery = 1;
end

% --- Stage 1: CWT Feature Generation ---
if orderx == ordery && orderx ~= 0 % Assuming non-zero equal orders might use 90 deg range
    betad = 90;
else
    betad = 180; % Standard angular range for CWT
end

Scales = 1:SDF:(N_Scales*SDF);
if N_Angles > 0
    angle_step = betad / N_Angles;
    Angless = 0 : angle_step : (betad - angle_step); % Ensure angles are within [0, betad)
else
    Angless = 0; % Default to a single angle if N_Angles is zero or invalid
end
Angles = Angless * pi / 180;

if ndims(CWT_Input_Features_Matrix) == 2
    N_Fs = 1;
else
    szCWT = size(CWT_Input_Features_Matrix);
    N_Fs = szCWT(3);
end
N_CWT = N_Scales * N_Angles * N_Fs;

h_cwt = waitbar(0.5,'Kernel CWT...', 'WindowStyle', 'modal', 'Name', 'CWT Computation');
frames_cwt = java.awt.Frame.getFrames();
if ~isempty(frames_cwt)
    frames_cwt(end).setAlwaysOnTop(1);
end

% Ensure pop_choice3 is set (e.g., from app state or default)
if isempty(pop_choice3)
    pop_choice3 = 'gaus'; % Default to 'gaus' if not set
    disp("Warning: global 'pop_choice3' for CWT type was empty, defaulted to 'gaus'.");
end

% Call the appropriate CWT function (assuming cwt2D_DerGus_BTO can handle derpoisson2d logic)
% Ensure cwt2D_DerGus_BTO is on the path
if (strcmp(pop_choice3,'gaus'))
    Fe_All = cwt2D_DerGus_BTO(CWT_Input_Features_Matrix, Scales, Angles, orderx, ordery, WSFR);

    elseif (strcmp(pop_choice3,'derpoisson2d'))
    Fe_All = cwt2D_DerPoisson_BTO(CWT_Input_Features_Matrix,Scales,Angles,orderx,ordery,WSFR,Z);   
    
else
    close(h_cwt);
    error(['Unsupported pop_choice3: ', pop_choice3]);
end
if ishandle(h_cwt), close(h_cwt); end
pause(0.01);

h_save_cwt = waitbar(0.5,'Save Kernel CWT Results...', 'WindowStyle', 'modal', 'Name', 'Saving CWT');
frames_save_cwt = java.awt.Frame.getFrames();
if ~isempty(frames_save_cwt)
    frames_save_cwt(end).setAlwaysOnTop(1);
end
[XI_siz, YI_siz, nn, mm, kk] = size(Fe_All);
% Result3D = reshape(Fe_All, XI_siz, YI_siz, nn*mm*kk);
Result3D = reshape(Fe_All, XI_siz, YI_siz, nn*mm, kk);

close(h_save_cwt);
pause(0.01);


Result3D_Size = size(Result3D, 3);
% DR is passed as a character string from optimizableVariable, convert if not already numeric
if ischar(DR)
    DR_numeric = str2double(DR);
else
    DR_numeric = DR;
end

if Result3D_Size >= DR_numeric
    DR2 = DR_numeric;
else
    DR2 = Result3D_Size;
end

if DR2 <= 0 % Ensure DR2 is at least 1 if there are features
    DR2 = 1;
end



            % Get the number of principal components to retain
            % D_Reduct0 = app.CWTReductionto.Value;
            % RudimentaryEdgePercentile = app.FeatureSaliencyEditField.Value;
            % MapCorrelationWeight = app.DissimilarityEditField.Value;

            input_ndims = ndims(Result3D);
            image_count = 1; % Default for a single image or 2D input

            if input_ndims == 2
                 % Handle single 2D input - unlikely for CWT which expects at least 3D (image x scales)
                 % This part might need re-evaluation based on how Result3D is populated for 2D.
                 % Assuming Result3D is at least 3D for meaningful CWT.
                 Spectral_PCA_Features_Matrix(:, :, 1) = Result3D;
                 % warning('Input is 2D. CWT typically operates on 3D data (image x scales). Please check input format.');
                 % Further processing for 2D would be different, skipping for now based on the core request.
                 % close(h); % Close waitbar if we exit early
                 % return; % Exit function
            elseif input_ndims == 3 % Single 3D image (m x n x scales/features)
                image_count = 1;
                % Process the single image
                % SPCA_Features = CWTFeatureSelection(Result3D, D_Reduct0, RudimentaryEdgePercentile, MapCorrelationWeight);
                SPCA_Features = CWTFeatureSelection(Result3D, DR2, FSaliency, FDissimilarity);

                % Initialize the output matrix
                [rows, cols, ~] = size(SPCA_Features);
                Spectral_PCA_Features_Matrix = zeros(rows, cols, image_count);

                % Apply masking and combine selected features
                Grid = zeros(rows, cols, size(SPCA_Features, 3)); % Temporary storage for interpolated maps
                SPCA_Features(SPCA_Features == 0) = 1e-12; % Handle zero values

                for k = 1:size(SPCA_Features, 3)
                    SPCA_Featuresk = in .* SPCA_Features(:, :, k);
                    SPCA_Featuresk(SPCA_Featuresk == 0) = NaN; % Set masked values to NaN
                    % Store the interpolated result
                    Grid(:, :, k) = nDstrb2D(SPCA_Featuresk); % Assuming nDstrb2D handles NaN and interpolates
                end

                % Sum the interpolated maps
                multiplied_map = Grid(:,:,1);
                for i = 2:size(Grid, 3)
                    multiplied_map = multiplied_map + Grid(:,:,i);
                end
                % multiplied_map = nDstrb2D(double(multiplied_map)); % This final interpolation might be redundant if nDstrb2D is applied per map

                % Store the combined map for this image
                Spectral_PCA_Features_Matrix(:, :, 1) = multiplied_map;

            elseif input_ndims == 4 % Multiple 3D images (m x n x scales/features x num_images)
                image_count = size(Result3D, 4);
                [rows, cols, ~] = size(Result3D(:,:,:,1)); % Get dimensions of a single image volume

                % Initialize the output matrix to store combined maps for each image
                Spectral_PCA_Features_Matrix = zeros(rows, cols, image_count);

                for img_idx = 1:image_count
                    % waitbar(img_idx / image_count, h, sprintf('Processing Image %d of %d...', img_idx, image_count));

                    % Extract the current image volume
                    current_image_volume = Result3D(:, :, :, img_idx);

                    % Apply CWT Feature Selection to the current image volume
                    % Assuming CWTFeatureSelection can handle a 3D input (m x n x scales/features)
                    SPCA_Features_current_image = CWTFeatureSelection(current_image_volume, DR2, FSaliency, FDissimilarity);

                    % --- Processing for the current image (similar to the ndims==3 case) ---
                    SPCA_Features_current_image(SPCA_Features_current_image == 0) = 1e-12; % Handle zero values

                    % Initialize temporary storage for interpolated maps for this image
                    Grid_current_image = zeros(rows, cols, size(SPCA_Features_current_image, 3));

                    for k = 1:size(SPCA_Features_current_image, 3)
                        SPCA_Featuresk = in .* SPCA_Features_current_image(:, :, k);
                        SPCA_Featuresk(SPCA_Featuresk == 0) = NaN; % Set masked values to NaN
                        % Store the interpolated result
                        Grid_current_image(:, :, k) = nDstrb2D(SPCA_Featuresk); % Assuming nDstrb2D handles NaN and interpolates
                    end

                    % Sum the interpolated maps for the current image
                    multiplied_map_current_image = Grid_current_image(:,:,1);
                    for i = 2:size(Grid_current_image, 3)
                        multiplied_map_current_image = multiplied_map_current_image + Grid_current_image(:,:,i);
                    end
                    % multiplied_map_current_image = nDstrb2D(double(multiplied_map_current_image)); % Final interpolation?

                    % Store the combined map for the current image in the output matrix
                    Spectral_PCA_Features_Matrix(:, :, img_idx) = multiplied_map_current_image;
                    % --- End of processing for the current image ---
                end
            end

            % SPCA_Features_Size is now the number of *features* for a single image,
            % not the number of images. You might want a different variable
            % to report the number of images processed.
            % app.CWTReductionto.Value = SPCA_Features_Size; % This line should probably be removed or adjusted


% --- Stage 3: Lineament Extraction per Feature ---
h_le = waitbar(0.66, 'Lineament Extraction...', 'WindowStyle', 'modal', 'Name', 'Lineament Extraction');
frames_le = java.awt.Frame.getFrames();
if ~isempty(frames_le)
    frames_le(end).setAlwaysOnTop(1);
end

originalArrayNoNaNPages = cat(3, Spectral_PCA_Features_Matrix, All_Data_Matrix);
clear Spectral_PCA_Features_Matrix; % Free memory

% Remove pages that are all NaN (if any, though masking might prevent this)
if ndims(originalArrayNoNaNPages) == 3
    nanPageIndices = all(all(isnan(originalArrayNoNaNPages), 1), 2);
    originalArrayNoNaNPages(:,:,squeeze(nanPageIndices)) = [];
end

if isempty(originalArrayNoNaNPages) || all(isnan(originalArrayNoNaNPages(:)))
    disp('Warning: All feature maps are NaN or empty after masking/selection.');
    L_val = Inf;
    if ishandle(h_le), close(h_le); end
    return;
end

if ndims(originalArrayNoNaNPages) == 2
    Number_of_Features_used_for_Lineaments_Detection = 1;
    % Ensure it's 3D for consistency in the loop
    temp_arr = originalArrayNoNaNPages;
    clear originalArrayNoNaNPages;
    originalArrayNoNaNPages(:,:,1) = temp_arr;
else
    Number_of_Features_used_for_Lineaments_Detection = size(originalArrayNoNaNPages, 3);
end

% Ensure X_Pix and Y_Pix are valid
if isempty(X_Pix) || isempty(Y_Pix) || X_Pix == 0 || Y_Pix == 0
    disp("Warning: X_Pix or Y_Pix is not valid. Using image dimensions.");
    X_Pix = YI_siz; % Assuming YI_siz corresponds to X dimension of image based on typical meshgrid
    Y_Pix = XI_siz; % Assuming XI_siz corresponds to Y dimension
end

% Ensure Line_Res is valid
if Line_Res <=0
    disp("Warning: Line_Res is invalid. Defaulting to 256.");
    Line_Res = 1;
end

resize_factor = Line_Res;

% Initialize XI_line_Auto and YI_line_Auto based on the *original* grid dimensions (XI, YI)
% These should be available as globals from the main app scope
if isempty(XI) || isempty(YI)
    disp("Error: Global XI or YI is empty. Cannot proceed with resizing coordinates.");
    L_val = Inf;
    if ishandle(h_le), close(h_le); end
    return;
end
XI_line_Auto = imresize(XI, resize_factor);
YI_line_Auto = imresize(YI, resize_factor);

% Calculate variance for dynamic grm2 adjustment
vari1 = zeros(Number_of_Features_used_for_Lineaments_Detection,1);
for ttt6_var = 1:Number_of_Features_used_for_Lineaments_Detection
    Grid0_var = originalArrayNoNaNPages(:,:,ttt6_var);
    % Ensure FilterB is on the path
    Grid2_var = FilterB(Grid0_var, sigma2);
    Grid2_var(isnan(Grid2_var)) = 0;
    if var(double(Grid2_var(:))) == 0
        vari1(ttt6_var) = 1; % Avoid division by zero if variance is zero
    else
        vari1(ttt6_var) = 1/var(double(Grid2_var(:)));
    end
end
if length(vari1) > 1
    vari1_std = zscore(vari1);
else
    vari1_std = zeros(size(vari1)); % Zscore of a single element is NaN, handle this.
end

plotH0 = zeros(size(XI_line_Auto,1), size(XI_line_Auto,2), Number_of_Features_used_for_Lineaments_Detection);

for ttt6_le = 1:Number_of_Features_used_for_Lineaments_Detection
    Grid0_le = originalArrayNoNaNPages(:,:,ttt6_le);
    Grid2_le = FilterB(Grid0_le, sigma2);
    Grid2_le(isnan(Grid2_le)) = 0;

    current_grm = w; % w is the parameter from BHO
    grm2 = round(current_grm + (VSFW * current_grm) * vari1_std(ttt6_le));
    if grm2 <= 1
        grm2 = 1;
    end

    Grid_resized = imresize(Grid2_le, resize_factor);
    data1 = im2double(Grid_resized);

    current_GSF_Angles = NSFA; % NSFA is the parameter from BHO
    Step_Filter_nAng = current_GSF_Angles;
    if Step_Filter_nAng > 0
        Step_Filter_step = pi / Step_Filter_nAng;
        Step_Filter_DTheta = 0:Step_Filter_step:pi-(Step_Filter_step);
    else
        Step_Filter_DTheta = 0; % Default if NSFA is 0
    end
    Step_Filter_Dbw = [2 4 6]; % This is fixed in your original code

    % Ensure step_Filtering is on the path
    Y_filtered = step_Filtering(data1, grm2, Step_Filter_DTheta, Step_Filter_Dbw, 2);

    % Ensure getFaultDetection is on the path
    [BW_detected, ~] = getFaultDetection(Y_filtered, data1, SLC);

    % Component-level filter based on strength
    CC_detected = bwconncomp(BW_detected);
    if CC_detected.NumObjects > 0
        stats_detected = regionprops(CC_detected, Y_filtered, 'PixelValues');
        strengths_detected = cellfun(@(x) sum(x(:)), {stats_detected.PixelValues});

        pct_cutoff = 100 - SLC;
        if isempty(strengths_detected)
            th_strength = 0;
        else
            th_strength = prctile(strengths_detected, pct_cutoff);
        end
        keep_indices = find(strengths_detected >= th_strength);
        BW_filtered_strength = ismember(labelmatrix(CC_detected), keep_indices);
    else
        BW_filtered_strength = false(size(BW_detected));
    end

    % Visualization part (plotNewColorMap0_) and subsequent processing for plotH00
    % This part seems to generate a binary image from a colored plot, which is indirect.
    % A more direct approach from BW_filtered_strength might be better if plotH_ is just for visualization.
    % Assuming plotNewColorMap0_ generates an RGB image and you need to convert it.
    if CC_detected.NumObjects > 0 && ~isempty(keep_indices)
        Map_for_plot = bwlabel(BW_filtered_strength, 8);
        Labels_for_plot = unique(Map_for_plot(Map_for_plot > 0))';

        % Ensure plotNewColorMap0_ is on the path
        if isempty(Labels_for_plot) % If no labels satisfy condition
            plotH_colored = zeros(size(data1,1), size(data1,2), 3); % black image
        else
            plotH_colored = plotNewColorMap0_(Y_filtered, zeros(size(data1)), Map_for_plot, Labels_for_plot, ...
                'Detected Lineaments', YI_line_Auto(:), XI_line_Auto(:));
        end
    else
        plotH_colored = zeros(size(data1,1), size(data1,2), 3); % Render as black if no components
    end

    plotH00_gray = (plotH_colored(:,:,1) + plotH_colored(:,:,2) + plotH_colored(:,:,3))/3;
    plotH00_bin = imbinarize(plotH00_gray);

    CC_final_filter = bwconncomp(plotH00_bin);
    if CC_final_filter.NumObjects > 0
        componentStats_final = regionprops(CC_final_filter, 'Area');
        % Using a fixed area threshold of 50 as in your code
        validComponents_final = find([componentStats_final.Area] >= 50);
        plotH00_final = ismember(labelmatrix(CC_final_filter), validComponents_final);
    else
        plotH00_final = false(size(plotH00_bin));
    end

    plotH0(:,:,ttt6_le) = imcomplement(plotH00_final); % Store complemented binary image
end

% Stack lineaments from all features
if Number_of_Features_used_for_Lineaments_Detection > 0
    % Initialize stacked image based on the size of one processed feature map
    plotH1_stacked = false(size(plotH0,1), size(plotH0,2));
    for n_stack = 1:Number_of_Features_used_for_Lineaments_Detection
        % Ensure plotH0 layers are logical for imfuse with 'blend' like behavior on binary
        % Or use direct logical ORing for binary images
        plotH1_stacked = plotH1_stacked | plotH0(:,:,n_stack);
    end
else
    plotH1_stacked = false(size(XI_line_Auto,1), size(XI_line_Auto,2)); % empty if no features
end

% Ensure YI_Real_Ratio is sensible (e.g. aspect ratio of original image pixels if not square)
if isempty(YI_Real_Ratio) || YI_Real_Ratio <= 0
    YI_Real_Ratio = 1.0; % Default if not set or invalid
end

plotH1_rgb_cropped = im2double(plotH1_stacked); % Convert logical to double for cropping
cropPixels1 = round(0.025 * YI_Real_Ratio * size(plotH1_stacked,1));
cropPixels2 = round(0.025 * size(plotH1_stacked,2));

% Ensure crop pixels are not too large
cropPixels1 = min(cropPixels1, floor(size(plotH1_rgb_cropped,1)/2)-1);
cropPixels2 = min(cropPixels2, floor(size(plotH1_rgb_cropped,2)/2)-1);
if cropPixels1 < 0, cropPixels1 = 0; end
if cropPixels2 < 0, cropPixels2 = 0; end

if cropPixels1 > 0
    plotH1_rgb_cropped(1:cropPixels1, :) = 1;
    plotH1_rgb_cropped(end-cropPixels1+1:end, :) = 1;
end
if cropPixels2 > 0
    plotH1_rgb_cropped(:, 1:cropPixels2) = 1;
    plotH1_rgb_cropped(:, end-cropPixels2+1:end, :) = 1;
end

if ishandle(h_le), close(h_le); end
pause(0.01);

% --- Stage 4: Compute Graph Representativeness (Objective Value) ---
h_graph = waitbar(0.5, 'Building Graph & Computing Representativeness...', 'WindowStyle', 'modal', 'Name', 'Graph Measure');
frames_graph = java.awt.Frame.getFrames();
if ~isempty(frames_graph)
    frames_graph(end).setAlwaysOnTop(1);
end

lineamentBinary = imcomplement(plotH1_rgb_cropped); % Assuming 0 is lineament, 1 is background from previous step
lineamentBinary = imbinarize(lineamentBinary); % Ensure it's binary {0,1}

% Density Map (used for segment weighting if getSegmentWeight requires it)
% This density map is based on the final cropped and stacked lineaments
grayImage_density = lineamentBinary; % Already binary
windowSize_density = 3; % As in your code
filterKernel_density = ones(windowSize_density, windowSize_density);
densityMap_final = conv2(double(grayImage_density), filterKernel_density, 'same');
densityMap_final = densityMap_final / (windowSize_density^2);
% densityMap_final = imcomplement(densityMap_final); % Your original code complements. Adjust if needed.

CC_segments = bwconncomp(lineamentBinary);
L_segments_pixels = CC_segments.PixelIdxList;

Segments = struct('startCoord', {}, 'endCoord', {}, 'pixelList', {});
if CC_segments.NumObjects > 0
    for sIdx = 1:length(L_segments_pixels)
        pixList = L_segments_pixels{sIdx};
        if numel(pixList) < 2 % Skip single-pixel segments or very short ones
            continue;
        end
        [rs, cs] = ind2sub(size(lineamentBinary), pixList);

        % A simple way to define start/end for non-skeletonized segments:
        % Find extreme points. For more robust, skeletonize first.
        % This is a simplification; proper skeletonization and endpoint/branchpoint detection is better.
        [~, min_idx] = min(rs); % Top-most as start_y
        [~, max_idx] = max(rs); % Bottom-most as end_y

        seg_struct.startCoord = [rs(min_idx(1)), cs(min_idx(1))]; % Take the first if multiple min/max
        seg_struct.endCoord   = [rs(max_idx(1)), cs(max_idx(1))];
        seg_struct.pixelList  = pixList;
        Segments = [Segments; seg_struct];
    end
end

ME_value = 0;
N_graph_segments = numel(Segments);

if N_graph_segments == 0
    ME_value = 0; % Or a penalty if no segments are found
else
    % --- Build Graph (simplified for ME_value calculation if full graph structure not strictly needed here) ---
    % The graph structure (G) was for visualization or other metrics in your example.
    % For ME_value, you primarily need segments and their pairwise similarities.

    % 1) Precompute segment weights (W_vals)
    % This requires a function `getSegmentWeight_auto`
    % Placeholder: using segment length as weight
    W_vals = zeros(N_graph_segments,1);
    for i = 1:N_graph_segments
        % W_vals(i) = getSegmentWeight_auto(Segments(i), densityMap_final); % Your custom function
        W_vals(i) = numel(Segments(i).pixelList); % Example: Length as weight
    end

    % 2) Precompute similarity matrix (simMatrix)
    % This requires a function `computeSimilarity_auto`
    simMatrix = zeros(N_graph_segments, N_graph_segments);
    if N_graph_segments > 1
        for i = 1:N_graph_segments
            for j = i+1:N_graph_segments % Only upper triangle
                % ST_val = computeSimilarity_auto(Segments(i), Segments(j), densityMap_final); % Your custom function
                % Placeholder: inverse of distance between centroids
                centroid_i = mean(ind2sub(size(lineamentBinary), Segments(i).pixelList),1);
                centroid_j = mean(ind2sub(size(lineamentBinary), Segments(j).pixelList),1);
                dist_centroids = norm(centroid_i - centroid_j);
                ST_val = 1 / (1 + dist_centroids); % Example similarity
                simMatrix(i,j) = ST_val;
                simMatrix(j,i) = ST_val; % Symmetric
            end
        end
    end

    % 3) Sum up for ME_value
    for segIndex = 1:N_graph_segments
        W_val_current = W_vals(segIndex);

        if N_graph_segments == 1
            hatS_T = 0; % No other segments to compare to
        else
            rowSimilarities = simMatrix(segIndex,:);
            rowSimilarities(segIndex) = Inf; % Exclude self-similarity by setting to Inf for min
            hatS_T = min(rowSimilarities);
            if isinf(hatS_T) % If only one segment, min(Inf) is Inf
                hatS_T = 0;
            end
        end

        ME_value = ME_value + W_val_current * (1 - hatS_T);
    end
end

if ishandle(h_graph), close(h_graph); end
pause(0.01);

zeta = 1; % Weighting factor for graph representativeness
L_val = -zeta * ME_value; % Bayesopt minimizes, so negate if ME_value is a "goodness" measure

% Handle cases where L_val might be NaN or Inf
if isnan(L_val) || isinf(L_val)
    L_val = Inf; % Penalize bad parameter sets heavily
end

% fprintf('Objective L_val: %f, N_CWT: %d\n', L_val, N_CWT);
end