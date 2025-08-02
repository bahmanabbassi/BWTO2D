%--------------------------------------------------------------------------
% nDstrb1D: A function for 1D data normalization and distribution adjustment
%--------------------------------------------------------------------------
%
% Developed by: Bahman Abbassi
% University:   University of Quebec (UQAT)
% Email:        bahman.abbassi@uqat.ca
% Github:       https://github.com/bahmanabbassi
% MathWorks:    https://www.mathworks.com/matlabcentral/profile/authors/31336389
%
%--------------------------------------------------------------------------
%
% Description:
%   This function iteratively transforms a 1D input data array (ImD0) to 
%   achieve an output array (ImD4) with a standard deviation within a 
%   predefined range [0.25, 0.30]. The process involves normalization,
%   logarithmic transformation, and histogram equalization.
%
%--------------------------------------------------------------------------
%
% Input:
%   ImD0 - A 1D numerical array (e.g., a vector or a flattened matrix).
%
% Output:
%   ImD4 - The transformed 1D numerical array.
%
%--------------------------------------------------------------------------
%
% Example:
%   % Generate some sample data
%   original_data = randn(1000, 1) * 5 + 10;
%
%   % Apply the transformation
%   transformed_data = nDstrb1D(original_data);
%
%   % Display the standard deviation of the output
%   disp(['Standard Deviation of Output: ', num2str(std(transformed_data(:)))]);
%
%--------------------------------------------------------------------------

function ImD4 = nDstrb1D(ImD0)
    %-- Step 1: Calculate initial statistics of the input data
    ImD0_STD = std(ImD0(:));   % Standard deviation of the input data
    ImD0_Mean = mean(ImD0(:)); % Mean of the input data

    %-- Step 2: Define control parameters for the iterative process
    max_iterations = 10; % Set the maximum number of loop iterations
    LT_STD = 0.25;       % Lower threshold for the target standard deviation
    HT_STD = 0.30;       % Higher threshold for the target standard deviation

%-- Step 3: Start the iterative transformation loop
for i = 1:1:max_iterations
    
    % Step 3a: Normalize the data based on its mean and standard deviation
    % The divisor 'i' adjusts the scaling factor in each iteration
    ImD1 = (ImD0 - ImD0_Mean)./(i*ImD0_STD);
    
    % Step 3b: Shift the data to ensure all values are positive
    ImD2 = -min(ImD1(:)) + ImD1 + (i*std(ImD1(:)));
    
    % Step 3c: Apply a logarithmic transformation to compress the data range
    ImD3 = log10(ImD2);
    
    % Step 3d: Apply histogram equalization to flatten the data distribution
    % 'realmin' is added to avoid potential zero values after equalization
    ImD4 = histeq(ImD3)+realmin;
    
    % Step 3e: Calculate the standard deviation of the newly transformed data
    ImD4_STD = std(ImD4(:));
    
    % Step 3f: Check if the standard deviation is within the target range
    if (ImD4_STD > LT_STD) && (ImD4_STD < HT_STD)
        break % Exit the loop if the condition is met
    end
end