function ImD4 = nDstrb2D(ImD0)
%nDstrb2D Adjusts the distribution of a 2D matrix with NaN handling.
%
% Description:
%   This function iteratively transforms a 2D input matrix (ImD0), which may
%   contain NaN values, to produce an output matrix (ImD4) where the standard
%   deviation of the non-NaN elements falls within a target range [0.25, 0.30].
%   The process correctly ignores NaN values during all statistical calculations
%   and transformations.
%
% Developed by:
%   Bahman Abbassi (bahman.abbassi@uqat.ca)
%   University of Quebec (UQAT)
%
% Input:
%   ImD0 - A 2D numerical matrix, which can include NaN values.
%
% Output:
%   ImD4 - The transformed 2D matrix with the same dimensions and NaN locations
%          as the input.
%
% Example:
%   % 1. Create a sample 2D matrix (e.g., using peaks)
%   sample_data = peaks(100);
%
%   % 2. Introduce some NaN values to simulate missing data
%   sample_data(20:40, 50:70) = NaN;
%
%   % 3. Apply the distribution adjustment
%   transformed_data = nDstrb2D(sample_data);
%
%   % 4. Visualize the results
%   figure;
%   subplot(1, 2, 1);
%   imagesc(sample_data);
%   title('Original Data');
%   axis image;
%   colorbar;
%
%   subplot(1, 2, 2);
%   imagesc(transformed_data);
%   title('Transformed Data');
%   axis image;
%   colorbar;


% --- Step 1: Input Validation and Initialization ---

    % Ensure the input is a 2D matrix, as the logic is tailored for it.
    if ndims(ImD0) ~= 2
        error('Input must be a 2D matrix.');
    end

    % Check for an all-NaN matrix, which cannot be processed.
    if all(isnan(ImD0(:)))
        error('The input image is all NaN; cannot proceed.');
    end

    % Create a logical mask to identify the locations of valid (non-NaN) numbers.
    % This mask is essential for all subsequent NaN-aware calculations.
    validMask = ~isnan(ImD0);

    % Extract only the valid values into a vector for initial analysis.
    validValues = ImD0(validMask);

    % Calculate initial statistics using only the non-NaN values.
    ImD0_STD = std(validValues);
    ImD0_Mean = mean(validValues);

    % --- Step 2: Define Control Parameters ---
    max_iterations = 50;  % Set the maximum number of attempts to reach the target.
    LT_STD = 0.25;        % Lower threshold for the target standard deviation.
    HT_STD = 0.30;        % Higher threshold for the target standard deviation.
    realmin_val = realmin;% Cache the value of realmin to avoid repeated calls.

    % --- Step 3: Prepare Data for Iteration ---

    % Center the data by subtracting the mean of valid values.
    ImD0_Centered = ImD0 - ImD0_Mean;
    % Temporarily set original NaN locations to zero in the centered data
    % to prevent NaN propagation in arithmetic operations.
    ImD0_Centered(~validMask) = 0;

    % Initialize the output matrix. If the loop finishes without success,
    % this version will be returned.
    ImD4 = ImD0;

    % Precompute scaling weights for each potential iteration to optimize the loop.
    i = 1:max_iterations;
    weights = 1 ./ (i * ImD0_STD);

    % --- Step 4: Iteratively Transform the Data ---
    for iter = 1:max_iterations
        % 4a. Normalize the centered data using the precomputed weight for this iteration.
        current_weight = weights(iter);
        ImD1 = ImD0_Centered * current_weight;

        % 4b. Shift the data to ensure all values are positive before the log transform.
        %     Calculations are performed only on the valid data region.
        min_val = min(ImD1(validMask));
        std_val = std(ImD1(validMask));
        ImD2 = ImD1 - min_val + (iter * std_val);

        % 4c. Apply a logarithmic transformation. `max(ImD2, 0)` prevents taking the
        %     log of negative numbers, and `realmin` prevents `log(0)`.
        ImD3 = log10(max(ImD2, 0) + realmin_val);

        % 4d. Apply histogram equalization to the valid data region to flatten its distribution.
        %     `real()` is a safeguard against any potential complex numbers.
        %     The result is placed back into the valid region of the output matrix.
        ImD4(validMask) = histeq(real(ImD3(validMask))) + realmin_val;

        % 4e. Calculate the standard deviation of the current result (valid data only).
        ImD4_STD = std(ImD4(validMask));

        % 4f. Check if the result's standard deviation is within the target range.
        if (ImD4_STD > LT_STD) && (ImD4_STD < HT_STD)
            break; % If successful, exit the loop.
        end
    end

    % --- Step 5: Finalize Output ---

    % Restore NaN values to the output matrix in their original locations.
    % This ensures the output has the same data mask as the input.
    ImD4(~validMask) = NaN;
end