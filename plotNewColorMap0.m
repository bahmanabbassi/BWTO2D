function [plotH] = plotNewColorMap0(RDepth, Depth, Map, Labels, tname, Y_cord, X_cord)
%
% =========================================================================
% FUNCTION: plotNewColorMap0
% =========================================================================
%
% PURPOSE: 🗺️
% This function generates a visualization by overlaying specified regions onto
% a grayscale depth image. The background is created from the original 'Depth'
% data, while the highlighted regions are color-coded based on values from
% the processed 'RDepth' data.
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
% This is a core visualization routine used in various geophysical and
% remote sensing analysis workflows developed in my lab.
%
% -------------------------------------------------------------------------
%
% INPUTS:
%   RDepth (2D matrix): Processed data whose values determine the color for
%                       the highlighted regions (e.g., fault intensity).
%   Depth  (2D matrix): The original data used for the grayscale background.
%   Map    (2D matrix): A label matrix where each integer identifies a distinct region.
%   Labels (vector):    A list of the integer labels from 'Map' that are to be
%                       overlaid and highlighted.
%   tname  (string):    (Unused) Intended for a plot title.
%   Y_cord (vector):    (Unused) Intended for y-axis coordinates.
%   X_cord (vector):    (Unused) Intended for x-axis coordinates.
%
% -------------------------------------------------------------------------
%
% OUTPUTS:
%   plotH (3D matrix): An RGB image with the specified regions highlighted
%                      on top of the grayscale background.
%
% -------------------------------------------------------------------------
%
% EXAMPLE USAGE:
%   % Assume R, D, and M are 2D matrices, and L is a vector of labels.
%   % L = [3, 7, 12]; % Labels of regions to highlight.
%   overlay_image = plotNewColorMap0(R, D, M, L, 'My Title', [], []);
%   figure;
%   imshow(overlay_image);
%   title('Highlighted Geological Features');
%
% =========================================================================

% --- Step 1: Create the Grayscale Background Image ---
% The original 'Depth' data is normalized to a [0, 1] range to serve as a
% standard grayscale base layer for the visualization.

% Calculate the dynamic range of the original depth data.
par = max(max((Depth))) - min(min(Depth));
% Normalize and assign the result to all three (R, G, B) channels.
plotH(:, :, 1) = (Depth - min(min(Depth))) / par; % Red channel
plotH(:, :, 2) = (Depth - min(min(Depth))) / par; % Green channel
plotH(:, :, 3) = (Depth - min(min(Depth))) / par; % Blue channel
clear par;

% --- Step 2: Determine Normalization Factor for Overlay Color ---
% The maximum value from the processed 'RDepth' data is used to scale the
% overlay colors, ensuring a consistent gradient.
normal = max(max(RDepth));

% --- Step 3: Overlay Highlighted Regions 🎨 ---
% This double loop iterates through each pixel of the label 'Map'. If a pixel's
% label matches one in the 'Labels' list, its color in the output image is updated.

for i = 1:size(Map, 1) % Iterate over each row.
    for j = 1:size(Map, 2) % Iterate over each column.

        % Check if the current pixel's label is in the list of regions to highlight.
        if Map(i, j) > 0 && ~isempty(find(Map(i, j) == Labels, 1))

            % Calculate the color value 'c' based on the processed data.
            c = abs(RDepth(i,j))/normal;

            % Update the RGB channels to color the pixel using a blue-to-red gradient.
            plotH(i, j, 1) = c;      % Red component increases with 'c'.
            plotH(i, j, 2) = 0;      % Green component is set to zero for a pure gradient.
            plotH(i, j, 3) = 1 - c;  % Blue component decreases as 'c' increases.
        end
    end
end

% The final 'plotH' image is now ready for display.
clear normal;
end