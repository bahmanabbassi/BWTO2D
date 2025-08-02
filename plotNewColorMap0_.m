function [plotH_] = plotNewColorMap0_(RDepth, Depth, Map, Labels, tname, Y_cord, X_cord)
%
% =========================================================================
% FUNCTION: plotNewColorMap0_
% =========================================================================
%
% PURPOSE: 🗺️
% This function creates a visualization by overlaying specific, labeled regions
% onto a grayscale background image. The background is generated from the
% original 'Depth' data, while the color of the highlighted regions is
% determined by the corresponding values in the 'RDepth' (processed) data.
%
% -------------------------------------------------------------------------
%
% AUTHOR:
% Bahman Abbassi, University of Quebec (UQAT)
% Email: bahman.abbassi@uqat.ca
%
% Modified from:
% Panagiotakis, C., 2025. Filters for Curvilinear Enhancement. MATLAB Central File Exchange. 
% https://www.mathworks.com/matlabcentral/fileexchange/46541-filters-for-curvilinear-enhancement (accessed August 1, 2025).
%
% -------------------------------------------------------------------------
%
% INPUTS:
%   RDepth (2D matrix): Processed data whose values will determine the color
%                       of the highlighted regions (e.g., a blue-to-red gradient).
%   Depth  (2D matrix): The original data used to create the grayscale background.
%   Map    (2D matrix): A label matrix where each integer corresponds to a detected region.
%   Labels (vector):    A list of the specific integer labels from 'Map' that should be
%                       highlighted in the final image.
%   tname  (string):    (Not used) Intended for a plot title.
%   Y_cord (vector):    (Not used) Intended for y-axis coordinates.
%   X_cord (vector):    (Not used) Intended for x-axis coordinates.
%
% -------------------------------------------------------------------------
%
% OUTPUTS:
%   plotH_ (3D matrix): An RGB image representing the final visualization with
%                       highlighted regions overlaid on the grayscale background.
%
% -------------------------------------------------------------------------
%
% EXAMPLE USAGE:
%   % Assume R, D, M, and L are pre-defined matrices and vectors.
%   % R = processed data, D = original data, M = label map, L = [2, 5]
%   highlighted_image = plotNewColorMap0_(R, D, M, L, 'Title', [], []);
%   figure;
%   imshow(highlighted_image);
%   title('Highlighted Regions 2 and 5');
%
% =========================================================================

% --- Step 1: Create the Grayscale Background Image ---
% The original 'Depth' data is normalized to the range [0, 1] to create
% a standard grayscale image that will serve as the base layer.

% Calculate the dynamic range of the original depth data.
par = max(max(Depth)) - min(min(Depth));
% Normalize the data and assign it to all three (R, G, B) channels.
plotH_(:, :, 1) = (Depth - min(min(Depth))) / par; % Red channel
plotH_(:, :, 2) = (Depth - min(min(Depth))) / par; % Green channel
plotH_(:, :, 3) = (Depth - min(min(Depth))) / par; % Blue channel
clear par;

% --- Step 2: Get Normalization Factor for the Overlay Color ---
% The maximum value of the processed data ('RDepth') is used to scale the
% color of the highlighted pixels, ensuring the color gradient is consistent.
normal = max(max(RDepth));

% --- Step 3: Overlay Highlighted Regions 🎨 ---
% This loop iterates through every pixel of the label 'Map'. If a pixel
% belongs to one of the regions specified in 'Labels', its color is replaced.

for i = 1:size(Map, 1) % Iterate through each row.
    for j = 1:size(Map, 2) % Iterate through each column.

        % Check if the pixel's label exists in the list of labels to highlight.
        if Map(i, j) > 0 && ~isempty(find(Map(i, j) == Labels, 1))

            % If it's a target region, color the pixel using a blue-to-red gradient
            % based on the value in the 'RDepth' matrix.
            plotH_(i, j, 1) = abs(RDepth(i, j)) / normal;      % Red component increases with RDepth value.
            plotH_(i, j, 2) = 0;                              % Green component is set to 0.
            plotH_(i, j, 3) = 1 - abs(RDepth(i, j)) / normal;  % Blue component decreases with RDepth value.
        end
    end
end

% The final output 'plotH_' is the complete RGB image.
clear normal;
end