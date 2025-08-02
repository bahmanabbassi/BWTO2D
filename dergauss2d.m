function [out] = dergauss2d(kx,ky,orderx,ordery)
%dergauss2d Computes the 2D derivative of a Gaussian wavelet in the frequency domain.
%
% Description:
%   This function computes the Fourier transform of a 2D wavelet derived
%   from the partial derivatives of a Gaussian function. The formula in the
%   frequency domain is:
%
%   PSIHAT(kx,ky) = (i*kx)^orderx * (i*ky)^ordery * exp(-(kx^2 + ky^2)/2)
%
%   This function is often used as a wavelet definition within a 2D Continuous
%   Wavelet Transform (CWT) framework, such as the `cwt2d` routine in the
%   original YAWTB toolbox.
%
% Developer Note:
%   This file has been annotated by Bahman Abbassi (bahman.abbassi@uqat.ca)
%   for clarity. The core algorithm and structure are from the YAWTB Team.
%
% Original Toolbox:
%   Yet Another Wavelet Toolbox (YAWTB)
%   http://www.fyma.ucl.ac.be/projects/yawtb
%
% Inputs:
%   kx, ky    (matrix) - Real-valued matrices representing the frequency plane,
%                        typically created using `meshgrid`.
%   orderx    (scalar) - The derivative order for the x-direction.
%   ordery    (scalar) - The derivative order for the y-direction.
%
% Output:
%   out       (matrix) - The complex-valued wavelet in the frequency domain.
%
% Example:
%   % 1. Define a frequency grid
%   step = 2*pi/128;
%   [kx, ky] = meshgrid(-pi:step:(pi-step));
%
%   % 2. Compute the wavelet for derivative orders (x=6, y=1)
%   wav = dergauss2d(kx, ky, 6, 1);
%
%   % 3. Visualize the real part of the wavelet
%   figure;
%   imagesc(real(wav));
%   title('Real Part of 2D Derivative of Gaussian Wavelet');
%   axis image;
%   colorbar;

% --- Original YAWTB Documentation (preserved) ---
% \manchap
%
% Compute the 2D multiple derivative of Gaussian
%
% \mansecSyntax
% [out] = polgauss2d(kx,ky,orderx,ordery)
%
% \mansecDescription
%
% This function computes the 2D multiple derivative of Gaussian.
% That is, the wavelet given by
% \begin{verbatim}
%   PSIHAT (kx,ky) = (i*kx).^orderx .* (i*ky)^ordery .* ...
%                              exp( - (kx.^2 + ky.^2) / 2 )
% \end{verbatim}
% where PSIHAT is the Fourier transform of PSI;
%
% This wavelet depends of two parameters: orderx and ordery.
% This function is used by the cwt2d routine which compute
% continuous wavelet transform in 2D.
%
% \mansubsecInputData
% \begin{description}
%
% \item[kx,ky] [REAL MATRICES]: The frequency plane. Use meshgrid
% to create it.
%
% \item[orderx, ordery] [REAL SCALARS]: The wavelet
% parameters.
%
% \end{description}
%
% \mansubsecOutputData
% \begin{description}
% \item[out] [REAL MATRIX]: The wavelet in frequency plane.
% \end{description}
%
% \mansecExample
%
% \begin{code}
% >> step = 2*pi/128;
% >> [kx,ky] = meshgrid( -pi : step : (pi-step) );
% >> wav = dergauss2d(kx,ky,6,1);
% >> imagesc(wav);
% \end{code}
%
% \mansecReference
%
% \mansecSeeAlso
%
% cwt2d meshgrid
%
% \mansecLicense
%
% This file is part of YAW Toolbox (Yet Another Wavelet Toolbox)
% You can get it at
% \url{"http://www.fyma.ucl.ac.be/projects/yawtb"}{"yawtb homepage"}
%
% $Header: /home/cvs/yawtb/continuous/2d/wave_defs/dergauss2d.m,v 1.4 2001-10-21 21:04:15 coron Exp $
%
% Copyright (C) 2001, the YAWTB Team (see the file AUTHORS distributed with
% this library) (See the notice at the end of the file.)

%% Define default parameter names and values for toolbox integration.
wavparval = {'orderx',1,'ordery',1};

%% Handle special input case for querying parameters.
% If the function is called with a single empty input, it returns the
% default parameter definitions. This is a common design pattern in toolboxes.
if ( (nargin == 1) && isempty(kx) )
  out = wavparval;
  return
end

%% Compute the wavelet in the frequency domain.
% This line implements the core formula:
% (i*kx)^orderx * (i*ky)^ordery * exp(-(kx^2 + ky^2)/2)
out = (1i*kx) .^ orderx .* (1i*ky) .^ ordery ...
      .* exp( - (kx.^2 + ky.^2) / 2 );

% %% The following is an alternative implementation for a Poisson wavelet,
% %% which is currently commented out and not active.
% out = (1i * kx) .^ orderx .* (1i * ky) .^ ordery .* sqrt(kx.^2 + ky.^2) .* exp(-sqrt(kx.^2 + ky.^2));

% --- Original GPL License Information (preserved) ---
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA