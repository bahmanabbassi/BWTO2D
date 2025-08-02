function [BW_object, plotSkel] = getFaultDetection(z, z_orig, SLC)
% getFaultDetection - Detects faults, cleans the resulting binary map, and generates visualizations.
% This function uses hysteresis thresholding via the 'binarize' sub-function and then refines
% the output by removing noise and reconnecting broken lineaments.
%
% Inputs:
%   z        (2D matrix): Processed data (e.g., filtered, gradient) used for detection.
%   z_orig   (2D matrix): Original, unprocessed data used for the background of the visualization.
%   SLC      (double):    A sensitivity parameter passed to the binarize function to set the initial threshold.
%
% Outputs:
%   BW_object (binary matrix): The final, cleaned binary map where '1' indicates a detected fault.
%   plotSkel  (RGB matrix):    An RGB image showing the detected faults
%   overlaid on the original data.  7.	
% Modified from:
% Panagiotakis, C., 2025. Filters for Curvilinear Enhancement. MATLAB Central File Exchange. 
% https://www.mathworks.com/matlabcentral/fileexchange/46541-filters-for-curvilinear-enhancement (accessed August 1, 2025).



% Suppress irrelevant warnings for a cleaner output.
warning('off','all');

% --- Step 1: Backup / Restore ---
% This step is a no-op in its current form and does not alter the data.
rotI = z; z = rotI;
clear rotI;

% --- Step 2: Hysteresis Thresholding & Morphological Cleanup
% Perform the core fault detection using hysteresis thresholding.
[BW] = binarize(z, SLC);
BW_object = BW;       % Store the raw binary map from the binarize function.

% --- Clean up the binary fault map ---
% Remove small, isolated pixel groups (noise) with an area less than 10 pixels.
BW_object = bwareaopen(BW_object, 10);
% Connect adjacent lineaments that have single-pixel breaks.
BW_object = bwmorph(BW_object,'bridge');
% Remove short, dead-end branches ("spurs") from the lineaments up to 50 pixels long.
BW_object = bwmorph(BW_object,'spur',50);
% --- End Cleanup ---

% --- Step 3: Prepare Background for Visualization
% Normalize the original data to create a grayscale background image for the overlay.
normal = max(z(:)); % Normalization factor from the processed data for fault coloring.
par    = max(z_orig(:)) - min(z_orig(:)); % Range of original data for background scaling.
plotH(:,:,1) = (z_orig - min(z_orig(:))) / par; % Red channel
plotH(:,:,2) = plotH(:,:,1); % Green channel
plotH(:,:,3) = plotH(:,:,1); % Blue channel
clear par;

% --- Step 4: Create Overlay Visualization ---
% Overlay the cleaned fault pixels onto the grayscale background image.
% The color of the fault varies from blue (low intensity) to red (high intensity).
for i=1:size(BW_object,1)
  for j=1:size(BW_object,2)
    % Check if the current pixel is part of a detected fault.
    if BW_object(i,j)
      plotH(i,j,1) = abs(z(i,j))/normal;      % Red component increases with fault intensity.
      plotH(i,j,2) = 0;                       % Green component is set to zero for faults.
      plotH(i,j,3) = 1 - abs(z(i,j))/normal;  % Blue component decreases with fault intensity.
    end
  end
end
% Assign the final overlay image to the output variable.
plotSkel = plotH;

% --- Step 5: Create Alternative Visualization (White Background) ---
% NOTE: This visualization is created but not returned by the function.
% It is useful for displaying faults without the underlying data.

% Create a blank white canvas.
plotH(:, :, 1) = ones(size(z_orig, 1), size(z_orig, 2)); % Red channel to 1.
plotH(:, :, 2) = ones(size(z_orig, 1), size(z_orig, 2)); % Green channel to 1.
plotH(:, :, 3) = ones(size(z_orig, 1), size(z_orig, 2)); % Blue channel to 1.

% Overlay faults onto the white background using the same coloring scheme.
for i = 1:size(BW, 1)
    for j = 1:size(BW, 2)
        if BW(i, j) > 0 % If the pixel belongs to a detected fault.
            plotH(i, j, 1) = abs(z(i, j)) / normal;         % Update red channel.
            plotH(i, j, 2) = 0;                             % Set green channel to 0.
            plotH(i, j, 3) = 1 - abs(z(i, j)) / normal;     % Update blue channel.
        end
    end
end
clear normal BW;
end





% This is the core helper function that implements the hysteresis thresholding algorithm. 
% It classifies pixels into strong and weak candidates based on initial thresholds and then 
% connects weak candidates that are adjacent to strong ones, resulting in a coherent binary fault map.

function [B1] = binarize(z, SLC)
% binarize - Converts input data into a binary fault map using hysteresis thresholding.
% This method identifies prominent features (strong faults) and then grows them into
% connected, weaker regions.
%
% Inputs:
%   z  (2D matrix): Processed input data where high values indicate potential faults.
%   SLC (double):   A sensitivity factor (e.g., 1.0) to control the initial threshold. Higher
%                   values make the initial selection less sensitive (requires stronger features).
%
% Outputs:
%   B1 (binary matrix): A raw binary map where '1' indicates a detected fault region.

% --- Step 1: Initialize Parameters 💡 ---
% Calculate a border/padding offset to avoid edge effects during window operations.
global YI_Real_Ratio
cropPixels1 = round(0.03 * YI_Real_Ratio * size(z,1), 0);
cropPixels2 = round(0.03 * size(z,2), 0);
apo = round((cropPixels1 + cropPixels2)/2,0);

% Initialize the output binary map with zeros.
B = 0 * z;

% --- Step 2: Determine Initial Threshold for Strong Features ---
% Flatten the data matrix into a single vector.
allD = reshape(z, 1, numel(z));
% Sort all pixel values in descending order.
allD = sort(allD, 'descend');
% Set a high threshold 'Tr' to select the most prominent features.
% This threshold is based on a top percentile of values, scaled by SLC.
Tr = allD(round(0.01* SLC * numel(z)));
clear allD;

% --- Step 3: Perform Initial Pixel Classification ---
% Classify pixels as 'strong' (1) or 'weak' (-1) based on the initial threshold 'Tr'.
for i = apo:size(z, 1) - apo
    for j = apo:size(z, 2) - apo
        % Extract a 3x3 window and sort its values.
        temp = reshape(z(i-1:i+1, j-1:j+1), 1, 9);
        st = sort(temp, 'descend');
        % Check if the center pixel exceeds the threshold.
        if z(i, j) > Tr
            % If it's a local maximum (top 2), classify as a strong candidate.
            if (z(i, j) >= st(2))
                B(i, j) = 1;
            % If it's prominent but not a peak (3rd or 4th), classify as a weak candidate.
            elseif (z(i, j) == st(4) || z(i, j) == st(3))
                B(i, j) = -1;
            end
        end
    end
end

% --- Step 4: Calculate High and Low Hysteresis Thresholds ---
% The actual high and low thresholds are determined from the mean values of the
% initially classified strong and weak pixels, respectively.
[Sx, Sy] = find(B == 1); % Find all strong candidate pixels.
y = zeros(1, length(Sx));
for i = 1:length(Sx)
    y(i) = z(Sx(i), Sy(i)); % Collect their values from the original data.
end
T_high = mean(y); % The high threshold is the mean of strong candidates.
clear Sx Sy y;

[Sx, Sy] = find(B == -1); % Find all weak candidate pixels.
y = zeros(1, length(Sx));
for i = 1:length(Sx)
    y(i) = z(Sx(i), Sy(i)); % Collect their values.
end
T_low = mean(y); % The low threshold is the mean of weak candidates.
clear Sx Sy y;

% --- Step 5: Reassign Binary Labels Based on Hysteresis Thresholds ---
% Re-classify all pixels to set up the final propagation step.
B = -2 * ones(size(z)); % Initialize map: -2 for non-fault, 1 for strong, 0 for undecided.
for i = apo:size(z, 1) - apo
    for j = apo:size(z, 2) - apo
        temp = reshape(z(i-1:i+1, j-1:j+1), 1, 9);
        st = sort(temp, 'descend');
        if z(i, j) <= T_low
            B(i, j) = -2; % Definitely not a fault.
        elseif z(i, j) >= T_high && (z(i, j) >= st(2))
            B(i, j) = 1; % Definitely a strong fault.
        elseif (z(i, j) < st(4))
            B(i, j) = -1; % Weak region.
        else
            B(i, j) = 0; % Undefined region between thresholds.
        end
    end
end

% --- Step 6: Propagate Strong Faults to Neighbors ---
% Iteratively connect undecided pixels (0) to adjacent strong fault pixels (1).
while 1
    [Sx, Sy] = find(B == 0); % Find all undecided pixels.
    change = 0; % Counter to track if any changes were made in an iteration.
    for i = 1:length(Sx)
        x = Sx(i); y = Sy(i);
        % Check if any neighbor in the 3x3 window is a strong fault.
        [u, v] = find(B(x-1:x+1, y-1:y+1) == 1);
        if ~isempty(u)
            B(x, y) = 1; % Promote the undecided pixel to a strong fault.
            change = change + 1;
        end
    end
    % If no pixels were promoted in a full pass, the propagation is complete.
    if change == 0
        break;
    end
    clear Sx Sy u v;
end

% --- Step 7: Finalize Binary Fault Map ---
% Set all remaining non-fault and weak pixels (values < 0) to 0.
B(B < 0) = 0;
% Assign the result to the output variable.
B1 = B;
clear B apo cropPixels1 cropPixels2 T_high T_low;
end