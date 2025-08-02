function [out] = derpoisson2d(kx,ky,orderx,ordery,Z_param)
%derpoisson2d Computes the 2D derivative of a Poisson wavelet in the frequency domain.
%
% Description:
%   This function computes the Fourier transform of a 2D wavelet based on the
%   derivatives of a Poisson function. The formula for the wavelet in the
%   frequency domain is:
%
%   PSIHAT(kx,ky) = (i*kx)^orderx * (i*ky)^ordery * ||k|| * exp(-Z * ||k||)
%
%   where ||k|| = sqrt(kx^2 + ky^2) is the radial frequency.
%
%   This function is intended to be used as a custom wavelet definition within a
%   2D Continuous Wavelet Transform (CWT) framework, such as one based on the
%   YAWTB (Yet Another Wavelet Toolbox).
%
% Developer:
%   Bahman Abbassi (bahman.abbassi@uqat.ca)
%   University of Quebec (UQAT)
%
% Inputs:
%   kx, ky    (matrix) - Real-valued matrices representing the frequency plane,
%                        typically created using `meshgrid`.
%   orderx    (scalar) - The derivative order for the x-direction.
%   ordery    (scalar) - The derivative order for the y-direction.
%   Z_param   (scalar) - The decay parameter `Z` for the exponential term.
%
% Output:
%   out       (matrix) - The complex-valued wavelet in the frequency domain.
%
% Example:
%   % 1. Define a frequency grid
%   step = 2*pi/256;
%   [kx, ky] = meshgrid(-pi:step:(pi-step));
%
%   % 2. Define wavelet parameters
%   x_order = 1;
%   y_order = 1;
%   z_val = 1.0;
%
%   % 3. Compute the Poisson derivative wavelet
%   wav = derpoisson2d(kx, ky, x_order, y_order, z_val);
%
%   % 4. Visualize the real part of the wavelet
%   figure;
%   imagesc(real(wav));
%   title('Real Part of 2D Poisson Derivative Wavelet');
%   axis image;
%   colorbar;


    % --- Parameter Definition for Toolbox Integration ---

    % Define the wavelet's parameters and their default values. This structure
    % is often used by higher-level CWT functions to manage inputs.
    wavparval = {'orderx',1,'ordery',1,'Z',1.0};

    % --- Handle Special Input Case: Parameter Query ---
    % If the function is called with a single, empty input, it returns the
    % parameter definitions. This allows a toolbox to query the function's parameters.
    if ( (nargin == 1) && isempty(kx) )
      out = wavparval;
      return
    end

    % --- Set Default Value for Z_param ---
    % This block robustly sets the value for Z_param. If it's not provided in the
    % function call, it retrieves the default value from the `wavparval` cell array.
    if nargin < 5
        idx = find(strcmp(wavparval, 'Z'));
        if ~isempty(idx) && length(wavparval) >= idx+1
            Z_param_to_use = wavparval{idx+1}; % Use default from wavparval.
        else
            Z_param_to_use = 1.0; % Fallback default.
        end
    else
        Z_param_to_use = Z_param; % Use the provided value.
    end

    % --- Compute the Wavelet ---
    % Calculate the radial frequency (magnitude of the wave vector).
    R_val = sqrt(kx.^2 + ky.^2);

    % Implement the core Poisson derivative wavelet formula in the frequency domain.
    out = (1i * kx) .^ orderx .* (1i * ky) .^ ordery .* R_val .* exp(-Z_param_to_use * R_val);

end