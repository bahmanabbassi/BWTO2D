function [Grid2, plt1, plt2] = Surf1(Grid,numberoffeatures,sigma2)
%Surf1 Processes and plots a series of 2D grids as surface plots.
%
% Description:
%   This function visualizes a set of 2D feature maps from a 3D grid.
%   It arranges them in a subplot grid, applies an optional Gaussian filter
%   to each map, and then renders each as a top-down surface plot.
%   The function relies on several global variables for coordinates (XI, YI),
%   region-of-interest boundaries (xv, yv), and aspect ratio correction
%   (YI_Real_Ratio).
%
% Developed by:
%   Bahman Abbassi (bahman.abbassi@uqat.ca)
%   University of Quebec (UQAT)
%
% Inputs:
%   Grid             - A 3D matrix (height x width x numberoffeatures) where
%                      each slice is a 2D feature map.
%   numberoffeatures - The number of feature maps from the Grid to plot.
%   sigma2           - The standard deviation (sigma) for the Gaussian filter.
%                      If sigma2 is 0, no filtering is applied.
%
% Outputs:
%   Grid2 - The 3D matrix of feature maps after optional filtering.
%   plt1  - The calculated number of rows for the subplot grid.
%   plt2  - The calculated number of columns for the subplot grid.
%
% Example:
%   % --- Setup for Example ---
%   % This function uses global variables, so we must define them first.
%   global XI YI YI_Real_Ratio xv yv;
%   [XI, YI] = meshgrid(1:100, 1:80); % Define coordinate grids
%   YI_Real_Ratio = 80/100;           % Define aspect ratio
%   xv = [10 90 90 10 10];            % Define polygon boundary X-coords
%   yv = [10 10 70 70 10];            % Define polygon boundary Y-coords
%
%   % Create a sample 3D grid with 4 feature maps
%   GridData = cat(3, peaks(100), membrane(1,49), rand(100), -peaks(100));
%   GridData = GridData(1:80, :, :); % Match grid dimensions
%
%   % Call the function to create the plots
%   figure('Position', [100, 100, 800, 800]);
%   [FilteredGrid, subplot_rows, subplot_cols] = Surf1(GridData, 4, 1.5);


    % --- Step 1: Initialize Plot Layout and Global Variables ---

    % Calculate the number of rows and columns needed for the subplot grid
    % to arrange the plots in a nearly square layout.
    plt1 = ceil(numberoffeatures^0.5);
    plt2 = plt1;

    % Declare the global variables that this function will use for plotting.
    global XI
    global YI
    global YI_Real_Ratio
    global xv
    global yv

    % Define a logical mask for a polygonal region of interest.
    % NOTE: The 'in' mask is calculated here but is not used in the current
    % version of the code. It was likely intended to mask the data.
    in = inpolygon(XI,YI,xv,yv);

    % --- Step 2: Loop Through Features and Plot ---
    for tttt6 = 1:numberoffeatures
        % Extract the current 2D feature map (slice) from the 3D grid.
        Grid0 = (Grid(:,:,tttt6));

        % Apply a conditional Gaussian filter for smoothing using the helper function.
        Grid2(:,:,tttt6) = FilterB(Grid0,sigma2);

        % Create a subplot at the correct position in the grid.
        subplot(plt1,plt2,tttt6), surf(XI,YI,Grid2(:,:,tttt6))

        % --- Step 3: Format the Current Subplot ---
        title(['# ',num2str(tttt6),''])    % Add a simple title.
        box on                             % Draw a box around the plot.
        set(gca,'TickDir','out','linewidth',1,'Layer', 'top'); % Format axes.
        colormap jet;                      % Apply the 'jet' colormap.
        axis equal                         % Ensure x and y axes have equal spacing.
        shading interp                     % Interpolate colors for a smooth surface.
        view(0,90)                         % Set the view to top-down.
        grid off                           % Turn off the grid lines.

        % Adjust the data aspect ratio to correct for non-square real-world
        % coordinates, preventing spatial distortion.
        daspect([YI_Real_Ratio 1 1])
    end
end