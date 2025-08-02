function [Y] = step_Filtering(I, grammes, DTheta, DPaxos, typeOfFilter)
%
% =========================================================================
% FUNCTION: step_Filtering
% =========================================================================
%
% PURPOSE: 🧭
% Applies a series of directional filters to an input image to enhance linear
% features (lineaments). It iterates through multiple orientations and scales,
% keeping the maximum filter response at each pixel.
%
% -------------------------------------------------------------------------
%
% AUTHOR:
% Bahman Abbassi, University of Quebec (UQAT)
% Email: bahman.abbassi@uqat.ca
%
% -------------------------------------------------------------------------
%
% INPUTS:
%   I            (2D matrix): The input image to be filtered.
%   grammes      (integer):   The size (width) of the filter kernel in pixels.
%   DTheta       (vector):    A vector of angles (in radians) for directional filtering.
%   DPaxos       (vector):    A vector of filter bandwidths (main wave widths).
%   typeOfFilter (integer):   Selects the filter type. 1 for a basic step filter,
%                             any other value for a smooth step filter.
%
% OUTPUTS:
%   Y (2D matrix): The filtered image where lineaments are enhanced.
%
% =========================================================================


% --- Step 1: Initialization ---
gabors = length(DTheta);    % Get the number of different angles to process.
gabors2 = length(DPaxos);   % Get the number of different filter widths to process.

% Pre-allocate the output matrix 'Y' with zeros. This is more efficient
% than growing the matrix inside the loop and avoids conditional checks.
Y = zeros(size(I), 'like', I);

% --- Step 2: Main Filtering Loop ---
% Loop over each specified filter width.
for j = 1:gabors2
    bw = DPaxos(j);     % Set the current filter bandwidth.
    
    % Loop over each specified angle.
    for i = 1:gabors
        theta = DTheta(i);  % Set the current filter orientation.
        
        % --- Generate the directional filter kernel ---
        if typeOfFilter == 1
            % Create a basic filter with sharp transitions.
            gb = mask_fn(grammes, bw, theta, 0);
        else
            % Create a filter with smooth (parabolic) transitions.
            gb = mask_fn_smooth(grammes, bw, theta, 0);
        end
        
        % Normalize the filter's energy to ensure responses from different
        % orientations are comparable and avoid amplification artifacts.
        gb = gb / sqrt(sum(sum(gb.^2)));
        
        % Apply the filter to the image using 2D convolution.
        % 'same' ensures the output image has the same size as the input.
        G = conv2(I, gb, 'same');
        
        % Update the output 'Y' by keeping the maximum response found so far
        % for each pixel. This is the key step for enhancing multi-directional lineaments.
        Y = max(Y, G);
        
        % Clear temporary variables in each iteration to manage memory.
        clear gb G
    end
end

% --- Step 3: Remove Border Artifacts ✂️ ---
% Convolution can create artifacts at the image edges. This step sets the
% border pixels to zero to produce a cleaner final output.
global YI_Real_Ratio
cropPixels1 = round(0.03 * YI_Real_Ratio * size(Y, 1), 0);
cropPixels2 = round(0.03 * size(Y, 2), 0);
Y(1:cropPixels1, :) = 0;                        % Zero out top border.
Y(:, 1:cropPixels2) = 0;                        % Zero out left border.
Y(size(I, 1)-[0:cropPixels1-1], :) = 0;        % Zero out bottom border.
Y(:, size(I, 2)-[0:cropPixels2-1]) = 0;        % Zero out right border.

% Clear remaining temporary variables.
clear gabors gabors2 cropPixels1 cropPixels2
end



function gb = mask_fn(grammes, bw, theta, toPlot)
% mask_fn - Creates a basic directional step filter kernel.
% This filter has a sharp +1 central band and -1 side bands.
%
% INPUTS:
%   grammes (integer): The size of the filter kernel.
%   bw      (double):  The bandwidth (width of the central positive band).
%   theta   (double):  The orientation of the filter in radians.
%   toPlot  (integer): Flag to plot the filter for debugging (1=yes, 0=no).
%
% OUTPUTS:
%   gb (2D matrix): The generated filter kernel.

% --- Step 1: Define Filter Grid ---
sz = grammes; % Set the kernel size.
if mod(sz, 2) == 0
    sz = sz + 1; % Ensure the kernel size is odd for a distinct center pixel.
end

% Create a coordinate grid for the filter.
[x, y] = meshgrid(-fix(sz/2):1:fix(sz/2), fix(sz/2):-1:fix(-sz/2));
% Rotate the coordinate grid by the specified angle 'theta'.
x_theta = x * cos(theta) + y * sin(theta);
y_theta = -x * sin(theta) + y * cos(theta);

% --- Step 2: Define the Filter Profile ---
gb = zeros(size(x));  % Pre-allocate the kernel matrix.
for i = 1:length(x)
    for j = 1:length(y)
        % Get the distance from the center along the rotated x-axis.
        d = abs(x_theta(i, j));
        if d < bw
            gb(i, j) = 1; % Assign +1 to the central positive band.
        elseif d <= 3 * bw
            gb(i, j) = -1; % Assign -1 to the outer negative bands.
        end
    end
end

% --- Step 3: Balance the Filter Weights ---
% This ensures the filter has a zero-mean response to uniform regions.
[n1, n2] = find(gb == -1); % Find indices of negative weights.
[p1, p2] = find(gb == 1);  % Find indices of positive weights.
% Calculate the normalization factor for negative weights to balance the sum.
nv_n1 = -length(p1) / length(n1);
% Apply the normalization.
for i = 1:length(n1)
    gb(n1(i), n2(i)) = nv_n1;
end

% --- Step 4: Normalize Filter Energy ---
gb = gb / sqrt(sum(sum(gb.^2)));

% Optional: Plot the filter for visual inspection.
if toPlot == 1
    mean(mean(gb));
end

clear x y x_theta y_theta n1 n2 p1 p2
end



function gb = mask_fn_smooth(grammes, bw, theta, toPlot)
% mask_fn_smooth - Creates a smooth directional step filter kernel.
% This filter has a parabolic profile for smoother transitions.
%
% INPUTS:
%   grammes (integer): The size of the filter kernel.
%   bw      (double):  The bandwidth (width of the central positive region).
%   theta   (double):  The orientation of the filter in radians.
%   toPlot  (integer): Flag to plot the filter for debugging (1=yes, 0=no).
%
% OUTPUTS:
%   gb (2D matrix): The generated smooth filter kernel.

% --- Step 1: Define Filter Grid ---
sz = grammes; % Set the kernel size.
if mod(sz, 2) == 0
    sz = sz + 1; % Ensure the kernel size is odd.
end

% Create and rotate the coordinate grid.
[x, y] = meshgrid(-fix(sz/2):1:fix(sz/2), fix(sz/2):-1:fix(-sz/2));
x_theta = x * cos(theta) + y * sin(theta);
y_theta = -x * sin(theta) + y * cos(theta);

% --- Step 2: Define the Smooth Filter Profile ---
gb = zeros(size(x));  % Pre-allocate the kernel matrix.
for i = 1:length(x)
    for j = 1:length(y)
        d = abs(x_theta(i, j)); % Distance from the center.
        if d < bw
            % Positive central region with a downward-opening parabola.
            gb(i, j) = 1 - ((d / bw)^2);
        elseif d <= 3 * bw
            % Negative outer regions with an upward-opening parabola.
            gb(i, j) = 0.5 * (-1 + (((d / bw) - 2)^2));
        end
    end
end

% --- Step 3: Balance the Filter Weights ---
% This ensures the sum of all weights is zero (zero-mean filter).
[u1, v1] = find(gb > 0); % Find indices of positive weights.
[u2, v2] = find(gb < 0); % Find indices of negative weights.

% Calculate the mean of the positive weights.
m1 = 0; N1 = length(u1);
for i = 1:length(u1)
    m1 = m1 + gb(u1(i), v1(i));
end
m1 = m1 / N1;

% Calculate the mean of the negative weights.
m2 = 0; N2 = length(u2);
for i = 1:length(u2)
    m2 = m2 + gb(u2(i), v2(i));
end
m2 = m2 / N2;

% Calculate and apply a scaling factor 'c' to the negative weights to balance the filter.
c = -m1 * N1 / (m2 * N2);
for i = 1:length(u2)
    gb(u2(i), v2(i)) = c * gb(u2(i), v2(i));
end

% --- Step 4: Normalize Filter Energy ---
gb = gb / sqrt(sum(sum(gb.^2)));

% Optional: Plot the filter for visual inspection.
if toPlot == 1
    mean(mean(gb));
end

clear x y x_theta y_theta u1 v1 u2 v2 m1 m2 N1 N2
end