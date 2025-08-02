function plotH1_rgb1 = Lineaments_All(Grid,numberoffeatures,sigma2, SLC)
%
% =========================================================================
% FUNCTION: Lineaments_All
% =========================================================================
%
% PURPOSE:
% This function processes a 3D grid of geophysical or remote sensing data
% to extract, classify, and visualize linear features (lineaments). It
% identifies lineaments in each 2D slice (feature), classifies them as
% 'deep' or 'shallow' based on data texture, and generates combined maps,
% orientation data (for rose diagrams), and exportable shapefiles.
%
% -------------------------------------------------------------------------
%
% AUTHOR:
% Bahman Abbassi, University of Quebec (UQAT)
% Email: bahman.abbassi@uqat.ca
% Github: https://github.com/bahmanabbassi
% Mathworks: https://www.mathworks.com/matlabcentral/profile/authors/31336389
%
% -------------------------------------------------------------------------
%
% INPUTS:
%   Grid             (3D matrix): A 3D array where each slice
%                                 (Grid(:,:,n)) represents a 2D data feature.
%   numberoffeatures (integer):   The number of features (slices) to process
%                                 from the Grid.
%   sigma2           (double):    The sigma value for the Gaussian filter,
%                                 controlling the degree of initial smoothing.
%   SLC              (double):    A parameter passed to the underlying
%                                 lineament detection functions (`Lineaments0`,
%                                 `Lineaments0_`).
%
% -------------------------------------------------------------------------
%
% OUTPUTS:
%   plotH1_rgb1 (3D matrix): An RGB image representing the combined map of
%                            all detected lineaments from all features.
%
% GLOBAL VARIABLES CREATED:
%   plotH1_rgb1        - Final combined image of all lineaments.
%   plotH1_rgb1_Deep   - Combined image of only 'deep' lineaments.
%   plotH1_rgb1_Shallow- Combined image of only 'shallow' lineaments.
%   Shape              - Shapefile structure for all lineaments.
%   Shape_Deep         - Shapefile structure for 'deep' lineaments.
%   Shape_Shallow      - Shapefile structure for 'shallow' lineaments.
%   Orient0            - A vector containing the orientations of all detected
%                        lineaments.
%   densityMap         - A 2D map showing the spatial density of lineaments.
%
% -------------------------------------------------------------------------
%
% EXAMPLE USAGE:
%   % Assume 'myGeophysicalData' is a 200x200x10 matrix
%   num_features = 10;
%   sigma_val = 1.5;
%   slc_param = 0.5;
%   all_lineaments_map = Lineaments_All(myGeophysicalData, num_features, sigma_val, slc_param);
%   figure; imshow(all_lineaments_map); title('All Detected Lineaments');
%
% =========================================================================

% --- Step 1: Initialization and Environment Setup ---

% Clear previous variables to prevent conflicts from prior runs.
clear plotH
clear global plotH
clear Shape
clear global Shape
clear Shape_Deep
clear global Shape_Deep
clear Shape_Shallow
clear global Shape_Shallow

% Declare variables as global to allow them to be accessed by other functions or workspaces.
global plotH
global Shape
global Shape_Deep
global Shape_Shallow

% Suppress all warnings for a cleaner command window output.
warning('off','all')

% Remove any slices from the input grid that consist entirely of NaN values.
originalArrayNoNaNPages = Grid;
nanPageIndices = all(all(isnan(originalArrayNoNaNPages), 1), 2);
originalArrayNoNaNPages(:,:,nanPageIndices) = [];

% Calculate potential subplot dimensions for visualizing individual features.
plt1 = ceil(numberoffeatures^0.5);
plt2 = plt1;

% Declare global variables to be used throughout the function.
global YI_Real_Ratio
global grm
global grm2
global VSFW

% Clear and declare a global vector to store filter widths.
clear grm2_vector
clear global grm2_vector
global grm2_vector

% Set curvilinearity control parameter.
curv = VSFW;

% --- Step 2: Pre-analysis and Variance Calculation ---
% This loop assesses the texture/variability of each feature before the main processing.

for tttt6 = 1:numberoffeatures
    % Extract the current 2D feature slice.
    Grid0 = (originalArrayNoNaNPages(:,:,tttt6));
    
    % Apply a Gaussian filter to smooth the feature data.
    Grid2 = FilterB(Grid0,sigma2);
    
    % Replace any remaining NaN values with 0.
    Grid2(isnan(Grid2)) = 0;
    
    % Calculate the inverse variance. High variance (heterogeneous) features get a low value.
    vari1(tttt6) = 1 / var(double(Grid2(:)));
end
% Standardize the variance scores to have a mean of 0 and a standard deviation of 1.
vari1_std = zscore(vari1);

% Initialize a counter for the total number of faults (lineaments).
total_faults = 0;

% Initialize a cell array to store orientation data from each feature.
allOrientations = {};

% --- Step 3: Main Lineament Extraction Loop ---
% This loop processes each feature to detect and analyze lineaments.

for tttt6 = 1:numberoffeatures
    % Extract the current 2D feature slice from the original grid.
    Grid0 = (Grid(:,:,tttt6));
    
    % Apply Gaussian filtering.
    Grid2 = FilterB(Grid0,sigma2);
    Grid2(isnan(Grid2)) = 0; % Replace NaNs.
    
    % Dynamically adjust the step filtering width ('grm2') based on the feature's variance.
    % Features with high variance (more complex) will have a smaller 'grm2'.
    grm2 = round(grm + (curv*grm)*vari1_std(tttt6));
    if grm2 <= 1
        grm2 = 1; % Ensure the width is at least 1.
    end
    grm2_vector(tttt6) = grm2; % Store the calculated width for this feature.
    
    % Perform initial lineament extraction (legacy or primary method).
    plotH(:,:,:,tttt6) = Lineaments0(Grid2, SLC);
    
    % Perform refined lineament extraction to get a binary map and fault count.
    [plotH_, Number_OF_Faults, BW1] = Lineaments0_(Grid2, SLC);
    
    % Generate orientation data from the binary lineament map for a rose diagram.
    [orientations360] = createRoseDiagrams(BW1);
    allOrientations{tttt6} = orientations360; % Store the orientations.
    
    % Add the number of faults from this feature to the total count.
    total_faults = total_faults + Number_OF_Faults;
    disp(['Number of faults in feature ', num2str(tttt6), ': ', num2str(Number_OF_Faults)]);
    
    % --- Post-process the extracted lineament image for clarity ---
    % Convert the multi-channel image to a single grayscale intensity image.
    plotH00 = (plotH_(:,:,1)+plotH_(:,:,2)+plotH_(:,:,3))/3;
    
    % Binarize the image to create a clean black and white map.
    plotH00 = imbinarize(plotH00);
    
    % Remove small, noisy components (artifacts) based on an area threshold of 10 pixels.
    CC = bwconncomp(plotH00);
    componentStats = regionprops(CC, 'Area');
    validComponents = find([componentStats.Area] >= 10);
    plotH00 = ismember(labelmatrix(CC), validComponents);
    
    % Invert the binary image so lineaments are black (0) and background is white (1).
    plotH00 = imcomplement(plotH00);
    
    % Store the final cleaned image for this feature.
    plotH0(:,:,tttt6) = plotH00;
end

% --- Step 4: Finalize Orientations and Display Summary ---

% Combine orientation data from all features into a single global vector.
global Orient0
Orient0 = vertcat(allOrientations{:});

% Display completion message and summary statistics.
disp('.')
disp('.')
disp('.')
disp('Manual Lineament Extraction Completed');
disp(['Total number of faults for all features: ', num2str(total_faults)]);

% --- Step 5: Create a Combined Map of ALL Lineaments ---

% Initialize a blank image (canvas) with the same dimensions as the processed features.
plotH1 = zeros(size(plotH_));

% Fuse (stack) the lineament maps from all features onto the canvas.
for n = 1:numberoffeatures
    plotH1 = imfuse(plotH1,plotH0(:,:,n),'blend');
    plotH1 = imbinarize(plotH1); % Binarize the result after each fusion.
end

% Calculate pixel dimensions to crop the border, removing potential edge artifacts.
global YI_Real_Ratio
cropPixels1 = round(0.025 * YI_Real_Ratio * size(plotH1,1), 0);
cropPixels2 = round(0.025 * size(plotH1,2), 0);

% Store the final cropped image in a global variable. This is the main function output.
clear global plotH1_rgb1
global plotH1_rgb1
plotH1_rgb1 = im2double(plotH1);
plotH1_rgb1(1:cropPixels1, :, :) = 1;             % Set top border to white.
plotH1_rgb1(end-cropPixels1+1:end, :, :) = 1;  % Set bottom border to white.
plotH1_rgb1(:, 1:cropPixels2, :) = 1;             % Set left border to white.
plotH1_rgb1(:, end-cropPixels2+1:end, :) = 1;  % Set right border to white.


% --- Step 6: Create a Combined Map of DEEP Lineaments ---
% 'Deep' lineaments are assumed to come from features where grm2_vector > average.

% Calculate the average step filtering width, which acts as a threshold.
avg_grm2 = mean(grm2_vector);

% Initialize a blank canvas for the deep lineaments map.
plotH1 = zeros(size(plotH_));

% Fuse only the lineaments from features classified as 'deep'.
for n = 1 : numberoffeatures
    if grm2_vector(n) > avg_grm2
        plotH1 = imfuse(plotH1, plotH0(:,:,n), 'blend');
        plotH1 = imbinarize(plotH1);
    end
end

% Crop the borders and store the result in a global variable.
clear global plotH1_rgb1_Deep
global plotH1_rgb1_Deep
plotH1_rgb1_Deep = im2double(plotH1);
plotH1_rgb1_Deep(1:cropPixels1, :, :) = 1;
plotH1_rgb1_Deep(end-cropPixels1+1:end, :, :) = 1;
plotH1_rgb1_Deep(:, 1:cropPixels2, :) = 1;
plotH1_rgb1_Deep(:, end-cropPixels2+1:end, :) = 1;

% --- Step 7: Create a Combined Map of SHALLOW Lineaments ---
% 'Shallow' lineaments are from features where grm2_vector <= average.

avg_grm2 = mean(grm2_vector);

% Initialize a blank canvas for the shallow lineaments map.
plotH1 = zeros(size(plotH_));

% Fuse only the lineaments from features classified as 'shallow'.
for n = 1 : numberoffeatures
    if grm2_vector(n) <= avg_grm2
        plotH1 = imfuse(plotH1, plotH0(:,:,n), 'blend');
        plotH1 = imbinarize(plotH1);
    end
end

% Crop the borders and store the result in a global variable.
clear global plotH1_rgb1_Shallow
global plotH1_rgb1_Shallow
plotH1_rgb1_Shallow = im2double(plotH1);
plotH1_rgb1_Shallow(1:cropPixels1, :, :) = 1;
plotH1_rgb1_Shallow(end-cropPixels1+1:end, :, :) = 1;
plotH1_rgb1_Shallow(:, 1:cropPixels2, :) = 1;
plotH1_rgb1_Shallow(:, end-cropPixels2+1:end, :) = 1;

% --- Step 8: Create Shapefiles from Lineaments for GIS ---

% --- Step 8a: Generate Shapefile for ALL Lineaments ---
Shape = struct('Geometry', {}, 'X', {}, 'Y', {}, 'Name', {}); % Initialize empty structure.
for n = 1:numberoffeatures
    clear BW
    BW = plotH0(:,:,n); % Get the binary lineament map.
    
    % Access global coordinate vectors.
    global XI_line
    global YI_line
    
    % --- Vectorize the raster lineaments ---
    % Invert image so lines are objects (1) on a background (0).
    BW = imcomplement(BW);
    % Skeletonize lines to ensure single-pixel width for accurate tracing.
    BW = bwmorph(BW, 'thin', Inf);
    
    % Trace the boundaries of each line object to get pixel coordinates.
    [B,L] = bwboundaries(BW, 'noholes');
    
    % Get the real-world coordinate vectors corresponding to the image dimensions.
    xVec = linspace(min(XI_line(:)), max(XI_line(:)), size(BW,2));
    yVec = linspace(min(YI_line(:)), max(YI_line(:)), size(BW,1));
    
    % Create a temporary structure for the current feature's shapes.
    currentShape = struct('Geometry', {}, 'X', {}, 'Y', {}, 'Name', {});
    for k = 1:length(B)
        boundary = B{k};
        rowIdx = boundary(:,1);  % Pixel Y
        colIdx = boundary(:,2);  % Pixel X
        
        % Convert pixel indices to real-world geographic/projected coordinates.
        xCoords = xVec(colIdx);
        yCoords = yVec(rowIdx);
        
        % Populate the shape structure fields.
        currentShape(k).Geometry = 'Line';
        currentShape(k).X = xCoords';
        currentShape(k).Y = yCoords';
        currentShape(k).Name = ['Fault_', num2str(k)];
    end
    currentShape = currentShape(:);  % Ensure it is a column vector.
    Shape = [Shape; currentShape]; % Append to the main shape structure.
end

% --- Step 8b: Generate Shapefile for DEEP Lineaments ---
Shape_Deep = struct('Geometry', {}, 'X', {}, 'Y', {}, 'Name', {}); % Initialize.
for n = 1:numberoffeatures
    % Process only if the feature is classified as 'deep'.
    if grm2_vector(n) > avg_grm2
        clear BW
        BW = plotH0(:,:,n); % Get the binary lineament map.
        
        % Vectorize the raster lineaments (same process as in 8a).
        BW = imcomplement(BW);
        BW = bwmorph(BW, 'thin', Inf);
        [B,L] = bwboundaries(BW, 'noholes');
        
        xVec = linspace(min(XI_line(:)), max(XI_line(:)), size(BW,2));
        yVec = linspace(min(YI_line(:)), max(YI_line(:)), size(BW,1));
        
        currentShape_Deep = struct('Geometry', {}, 'X', {}, 'Y', {}, 'Name', {});
        for k = 1:length(B)
            boundary = B{k};
            rowIdx = boundary(:,1);
            colIdx = boundary(:,2);
            
            xCoords = xVec(colIdx);
            yCoords = yVec(rowIdx);
            
            currentShape_Deep(k).Geometry = 'Line';
            currentShape_Deep(k).X = xCoords';
            currentShape_Deep(k).Y = yCoords';
            currentShape_Deep(k).Name = ['Fault_', num2str(k)];
        end
        currentShape_Deep = currentShape_Deep(:);
        Shape_Deep = [Shape_Deep; currentShape_Deep]; % Append to the 'Deep' shape structure.
    end
end

% --- Step 8c: Generate Shapefile for SHALLOW Lineaments ---
Shape_Shallow = struct('Geometry', {}, 'X', {}, 'Y', {}, 'Name', {}); % Initialize.
for n = 1:numberoffeatures
    % Process only if the feature is classified as 'shallow'.
    if grm2_vector(n) < avg_grm2
        clear BW
        BW = plotH0(:,:,n); % Get the binary lineament map.
        
        % Vectorize the raster lineaments (same process as in 8a).
        BW = imcomplement(BW);
        BW = bwmorph(BW, 'thin', Inf);
        [B,L] = bwboundaries(BW, 'noholes');
        
        xVec = linspace(min(XI_line(:)), max(XI_line(:)), size(BW,2));
        yVec = linspace(min(YI_line(:)), max(YI_line(:)), size(BW,1));
        
        currentShape_Shallow = struct('Geometry', {}, 'X', {}, 'Y', {}, 'Name', {});
        for k = 1:length(B)
            boundary = B{k};
            rowIdx = boundary(:,1);
            colIdx = boundary(:,2);
            
            xCoords = xVec(colIdx);
            yCoords = yVec(rowIdx);
            
            currentShape_Shallow(k).Geometry = 'Line';
            currentShape_Shallow(k).X = xCoords';
            currentShape_Shallow(k).Y = yCoords';
            currentShape_Shallow(k).Name = ['Fault_', num2str(k)];
        end
        currentShape_Shallow = currentShape_Shallow(:);
        Shape_Shallow = [Shape_Shallow; currentShape_Shallow]; % Append to the 'Shallow' structure.
    end
end

% --- Step 9: Create Lineament Density Map ---

% Convert the final combined RGB image to a grayscale image.
grayImage = rgb2gray(double(plotH1_rgb1));

% Binarize the image (lineaments should be 1, background 0).
binaryImage = imbinarize(grayImage);

% Define the size of the moving window for density calculation.
windowSize = 3;
% Create the kernel for convolution (a simple averaging filter).
filterKernel = ones(windowSize, windowSize);

% Use 2D convolution to count lineament pixels within the moving window at each point.
global densityMap
densityMap = conv2(double(binaryImage), filterKernel, 'same');

% Normalize the map by the window area to get a density value between 0 and 1.
densityMap = densityMap / (windowSize^2);

% Invert the map for visualization (high density appears dark, low density appears light).
densityMap = imcomplement(densityMap);

end