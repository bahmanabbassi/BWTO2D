function F_Grid = FilterB(Grid,sigma2)
%FilterB Conditionally applies a Gaussian filter to a 2D grid.
%
% Description:
%   This function applies a Gaussian filter with a specified standard
%   deviation (sigma) to an input grid. If the provided sigma value is
%   zero or negative, the function returns the original grid unmodified.
%   It serves as a convenient wrapper to bypass filtering when not needed.
%
% Developed by:
%   Bahman Abbassi (bahman.abbassi@uqat.ca)
%   University of Quebec (UQAT)
%
% Inputs:
%   Grid   - A 2D numerical matrix (e.g., an image or data grid).
%   sigma2 - The standard deviation (sigma) for the Gaussian filter.
%
% Output:
%   F_Grid - The filtered grid, or the original grid if sigma2 <= 0.
%
% Example:
%   % 1. Create a sample grid and add some noise
%   my_grid = peaks(100);
%   noisy_grid = my_grid + 0.5 * randn(100);
%
%   % 2. Apply the filter with a positive sigma
%   filtered_grid = FilterB(noisy_grid, 2.0); % Filtering will be applied
%
%   % 3. "Apply" the filter with a zero sigma
%   unfiltered_grid = FilterB(noisy_grid, 0); % Filtering will be skipped
%
%   % 4. Visualize the results
%   figure;
%   subplot(1,3,1); imagesc(noisy_grid); title('Noisy Grid'); axis image;
%   subplot(1,3,2); imagesc(filtered_grid); title('Filtered (Sigma=2)'); axis image;
%   subplot(1,3,3); imagesc(unfiltered_grid); title('Filtered (Sigma=0)'); axis image;


    % Check if the provided sigma value is positive.
    % A Gaussian filter requires a sigma greater than zero.
    if sigma2 > 0
        % Apply the 2D Gaussian smoothing filter to the grid.
        F_Grid  = imgaussfilt(Grid,sigma2);
    else
        % If sigma is zero or negative, skip filtering and return the original grid.
        F_Grid = Grid;
    end
end