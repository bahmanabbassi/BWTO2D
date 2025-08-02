function Fe_All = cwt2D_DerGus(CWT_Input_Features_Matrix,Scales,Angles,orderx,ordery,Sigma3,valu_auto)
%cwt2D_DerGus Performs 2D CWT with a derivative-of-Gaussian wavelet and post-processing.
%
% Description:
%   This function applies a 2D Continuous Wavelet Transform (CWT) to one or
%   more input feature maps using a 2D derivative-of-Gaussian as the mother
%   wavelet. It processes the CWT coefficients for each specified scale and
%   angle by optionally applying a scale-dependent Gaussian filter and then
%   normalizing the statistical distribution of the coefficient magnitudes.
%
%   Finally, it generates a graphical user interface (GUI) with tabbed plots
%   to visualize the mother wavelet in the position domain, 2D frequency
%   domain, and 3D frequency domain.
%
% Developer:
%   Bahman Abbassi (bahman.abbassi@uqat.ca)
%   University of Quebec (UQAT)
%
% Inputs:
%   CWT_Input_Features_Matrix (double) - A 2D (H x W) or 3D (H x W x N) matrix
%                                        containing the input feature map(s).
%   Scales (vector)                    - A 1D vector of scales for the CWT.
%   Angles (vector)                    - A 1D vector of angles (in radians) for the CWT.
%   orderx (scalar)                    - The derivative order in the x-direction for the wavelet.
%   ordery (scalar)                    - The derivative order in the y-direction for the wavelet.
%   Sigma3 (scalar)                    - Sigma for the post-CWT Gaussian filter. If > 0, a
%                                        scale-dependent filter (sigma = Sigma3 * current_scale)
%                                        is applied. If <= 0, filtering is skipped.
%   valu_auto (any)                    - This input is not used within the function.
%
% Output:
%   Fe_All (double) - A 5D matrix (H x W x num_scales x num_angles x num_features)
%                     containing the final processed CWT coefficients.
%
% Example:
%   % 1. Create a sample 128x128 feature map
%   feature_map = peaks(128);
%
%   % 2. Define CWT parameters
%   scales_vec = 2.^(2:5);      % 4 scales from 4 to 32
%   angles_vec = 0:pi/4:3*pi/4; % 4 angles from 0 to 135 degrees
%   x_derivative_order = 2;
%   y_derivative_order = 1;
%   post_filter_sigma = 0.5;
%
%   % 3. Run the CWT and post-processing function
%   Processed_Coeffs = cwt2D_DerGus(feature_map, scales_vec, angles_vec, ...
%                                   x_derivative_order, y_derivative_order, ...
%                                   post_filter_sigma, []);
%
%   % The function will also generate a figure with plots of the mother wavelet.

% --- Step 1: Initialization and Pre-computation ---

% Replace any NaN values in the input with 0 to prevent computation errors.
CWT_Input_Features_Matrix(isnan(CWT_Input_Features_Matrix)) = 0;

% Note: Global variables XI and YI are declared but not used in this function.
global XI
global YI

% Get the number of scales and angles from the input vectors.
N_Scales = numel(Scales);
N_Angles = numel(Angles);

% Determine the number of input features (N_Fs).
if ndims(CWT_Input_Features_Matrix) == 2
    N_Fs = 1; % A single 2D feature map.
elseif ndims(CWT_Input_Features_Matrix) == 3
    N_Fs = size(CWT_Input_Features_Matrix, 3); % Multiple feature maps in a 3D stack.
end

% Clear temporary variables to ensure a clean state (optional).
clear Feature cwtmor cwtout cwtout_abs WL WL2 WLFs Fe_All WLFs_allscales WLFs_allangles

% --- Step 2: Main Loop for CWT and Post-Processing ---
% This loop iterates through each input feature map.
for kkk = 1:N_Fs
    % Select the current feature map.
    Feature = CWT_Input_Features_Matrix(:,:,kkk);

    % --- 2a. Perform the 2D Continuous Wavelet Transform ---
    tX = fft2(double(Feature)); % Precompute the 2D FFT of the feature.
    % Compute the CWT using the 'dergauss2d' wavelet definition.
    cwtmor = cwt2d(tX,'dergauss2d',Scales,Angles,orderx,ordery);
    cwtout = cwtmor.data; % Extract the raw CWT coefficient data.
    cwtout_abs = abs(cwtout); % Calculate the magnitude of the coefficients.

    % --- 2b. Post-process CWT Coefficients ---
    % Loop through each angle of the CWT result.
    for jj = 1:N_Angles
        % Loop through each scale for the current angle.
        for jjj = 1:N_Scales
            % Extract the 2D coefficient magnitude map for the current scale and angle.
            WL = cwtout_abs(:,:,jjj,jj);
            % Apply conditional, scale-dependent filtering and distribution normalization.
            if Sigma3 > 0
                % Apply a Gaussian filter with sigma increasing with scale, then normalize.
                WL2 = nDstrb2D(FilterB(WL,Sigma3*jjj));
            else
                % Normalize without filtering.
                WL2 = nDstrb2D(WL);
            end
            % Store the processed map for the current scale.
            WLFs_allscales(:,:,jjj) = WL2;
        end
        % Assemble the processed maps for all scales into an angle-specific 4D matrix.
        WLFs_allangles(:,:,:,jj) = WLFs_allscales;
    end
    % Assemble the results for the current feature into the final 5D output matrix.
    Fe_All(:,:,:,:,kkk) = WLFs_allangles;
end

% --- Step 3: Visualize the Mother Wavelet ---
disp('Mother Wavelet: 2D Derivatives of Gaussian Kernel');

% Close any pre-existing figures with specific numbers to avoid clutter.
figHandles = findobj('Type', 'figure', 'Number', [2, 3, 4]);
if ~isempty(figHandles)
    close(figHandles);
end

% Create the main GUI window for wavelet visualization.
f = figure('Name', 'Derivative of Gaussian Wavelet Viewer', ...
           'Position', [300, 200, 500, 500]);
f.WindowState = 'maximized'; % Open the figure in a maximized state.

% Create a tab group within the figure.
tg = uitabgroup(f);

% === Tab 1: Wavelet in the Position (Spatial) Domain ===
tab1 = uitab(tg, 'Title', 'Position Domain');
ax1 = axes('Parent', tab1);
yashow_cwt2d(cwtmor, 'filter', 'pos'); % Use YAWTB plotting tool.
xlim([-200, 200]);
ylim([-200, 200]);
set(gca, 'FontSize', 18, 'FontName', 'Times New Roman');
grid on;
colorbar;
zoom(20); % Zoom in for a better view.

% === Tab 2: Wavelet in the Frequency Domain (2D Top-Down View) ===
tab2 = uitab(tg, 'Title', 'Frequency Domain');
ax2 = axes('Parent', tab2);
yashow_cwt2d(cwtmor, 'filter');
xlim([-3, 3]);
ylim([-3, 3]);
axis equal;
set(gca, 'FontSize', 30, 'FontName', 'Times New Roman');
grid on;
colorbar;
zoom(2.1);

% === Tab 3: Wavelet in the Frequency Domain (3D Surface View) ===
tab3 = uitab(tg, 'Title', 'Frequency Domain (3D)');
ax3 = axes('Parent', tab3);
yashow_cwt2d(cwtmor, 'filter', 'surf'); % Plot as a surface.
xlim([-3, 3]);
ylim([-3, 3]);
axis equal;
set(gca, 'FontSize', 30, 'FontName', 'Times New Roman');
grid on;
colorbar;
view(45, 10)  % Set a custom 3D viewing angle (azimuth, elevation).
zoom(1);

end