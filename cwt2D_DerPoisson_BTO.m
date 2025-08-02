function Fe_All = cwt2D_DerPoisson_BTO(CWT_Input_Features_Matrix, Scales, Angles, orderx, ordery, Sigma3, Z)
%--------------------------------------------------------------------------
% FUNCTION: cwt2D_DerPoisson_BTO
%
% PURPOSE:
%   Performs a 2D Continuous Wavelet Transform (CWT) using a 2D derivative
%   of a Poisson wavelet ('derpoisson2d'). This function is analogous to
%   `cwt2D_DerGus_BTO` but uses a different mother wavelet, which requires an
%   additional parameter, Z. It processes input data across multiple scales
%   and angles to generate feature maps for further analysis.
%
% AUTHOR:
%   Bahman Abbassi
%   University of Quebec (UQAT)
%   Email: bahman.abbassi@uqat.ca
%   GitHub: https://github.com/bahmanabbassi
%   MathWorks: https://www.mathworks.com/matlabcentral/profile/authors/31336389
%
% INPUTS:
%   CWT_Input_Features_Matrix (double): A 2D (H x W) or 3D (H x W x F)
%                                       matrix of input data.
%   Scales (double array):              A 1D array of scale values for the CWT.
%   Angles (double array):              A 1D array of angles (in radians) for the CWT.
%   orderx (integer):                   The derivative order in the x-direction.
%   ordery (integer):                   The derivative order in the y-direction.
%   Sigma3 (double):                    Scaling factor for an optional Gaussian
%                                       filter. Set to 0 to disable.
%   Z (double):                         A specific parameter required by the
%                                       Poisson wavelet.
%
% OUTPUTS:
%   Fe_All (double):                    A 5D matrix (H x W x N_Scales x N_Angles x N_Features)
%                                       containing the absolute values of the CWT
%                                       coefficients.
%--------------------------------------------------------------------------

%% --- 1. Pre-processing and Initialization ---

% Replace any NaN values in the input matrix with 0.
CWT_Input_Features_Matrix(isnan(CWT_Input_Features_Matrix)) = 0;

% Get the number of scales and angles from the input arrays.
N_Scales = numel(Scales);
N_Angles = numel(Angles);

% Clear temporary variables to ensure a clean state before the loops.
clear Feature cwtmor cwtout cwtout_abs WL WL2 WLFs Fe_All WLFs_allscales WLFs_allangles

% Determine the number of input feature maps.
if ndims(CWT_Input_Features_Matrix) == 2
    N_Fs = 1; % A single 2D feature map.
elseif ndims(CWT_Input_Features_Matrix) == 3
    N_Fs = size(CWT_Input_Features_Matrix, 3); % A stack of feature maps.
end

%% --- 2. CWT Computation Loop ---
% Iterate through each feature map in the input matrix.
for kkk = 1:N_Fs
    % Select the current feature map.
    Feature = CWT_Input_Features_Matrix(:, :, kkk);
    
    % Pre-compute the 2D Fast Fourier Transform (FFT) for speed.
    tX = fft2(double(Feature));
    
    % Perform the 2D CWT using the 'derpoisson2d' wavelet.
    % Note the additional input parameter 'Z' required by this wavelet.
    cwtmor = cwt2d(tX, 'derpoisson2d', Scales, Angles, orderx, ordery, Z);
    
    % Extract the complex CWT coefficient data.
    cwtout = cwtmor.data; % Dimensions: H x W x N_Scales x N_Angles
    
    % Compute the absolute value (magnitude) of the coefficients.
    cwtout_abs = abs(cwtout);
    
    %% --- 3. Post-processing and Assembly ---
    % Loop through each angle and scale to apply optional filtering.
    for jj = 1:N_Angles
        for jjj = 1:N_Scales
            % Extract the CWT coefficients for the current scale and angle.
            WL = cwtout_abs(:, :, jjj, jj);
            
            % If Sigma3 > 0, apply a Gaussian filter. The filter's sigma
            % is scaled by the current CWT scale (jjj) for adaptive smoothing.
            if Sigma3 > 0
                WL2 = nDstrb2D(FilterB(WL, Sigma3 * jjj));
            else
                % Otherwise, just ensure there are no NaNs.
                WL2 = nDstrb2D(WL);
            end
            
            % Store the processed map for the current scale.
            WLFs_allscales(:, :, jjj) = WL2;
        end
        % Assemble all processed scales for the current angle.
        WLFs_allangles(:, :, :, jj) = WLFs_allscales;
    end
    
    % Assemble the results for the current feature into the final 5D output matrix.
    Fe_All(:, :, :, :, kkk) = WLFs_allangles;
end

end