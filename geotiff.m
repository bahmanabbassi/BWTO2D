function geotiff(matrix, xmin, xmax, ymin, ymax, defaultFilename)

%
% =========================================================================
% FUNCTION: geotiff
% =========================================================================
%
% PURPOSE: 💾
% Exports a MATLAB matrix to a georeferenced TIFF (.tif) file. The function
% prompts the user for a save location and formats the output to be
% broadly compatible with GIS software like QGIS, ArcGIS, and Geosoft.
%
% -------------------------------------------------------------------------
%
% AUTHOR:
% Bahman Abbassi, University of Quebec (UQAT)
% Email: bahman.abbassi@uqat.ca
%
% -------------------------------------------------------------------------
%
% INPUTS:
%   matrix          (2D matrix): The numerical data grid to be exported.
%   xmin, xmax      (double):    The longitude (X-coordinate) limits of the grid.
%   ymin, ymax      (double):    The latitude (Y-coordinate) limits of the grid.
%   defaultFilename (string):    A suggested filename to present in the save dialog.
%
% -------------------------------------------------------------------------
%
% OUTPUTS:
%   This function does not return any variables to the MATLAB workspace.
%   It writes one file to the disk:
%     - A .tif file containing the georeferenced raster data.
%
%   NOTE: The console message mentions .xml and .gi files for Geosoft
%   compatibility, but this version of the code does not generate them.
%
% -------------------------------------------------------------------------
%
% EXAMPLE USAGE:
%   % Create some sample data and define its geographic extent
%   myData = peaks(100);
%   x_min = -79.5; x_max = -78.5; % Longitude
%   y_min = 48.0;  y_max = 49.0;  % Latitude
%   geotiff(myData, x_min, x_max, y_min, y_max, 'Sample_Peaks_Data.tif');
%
% =========================================================================



    % Exports a matrix to Geosoft-compatible GeoTIFF with .xml and .gi metadata.

    % Ensure the default filename has a .tif extension
    [filepath, name, ext] = fileparts(defaultFilename);
    if isempty(ext) || ~strcmpi(ext, '.tif')
        defaultFilename = fullfile(filepath, [name, '.tif']);
    end

    % Prompt the user to select the output file location
    [filename, filepath] = uiputfile('*.tif', 'Save GeoTIFF File As', defaultFilename);
    if isequal(filename, 0) || isequal(filepath, 0)
        disp('File save canceled by user.');
        return;
    end

    % Full base path for output files
    basePath = fullfile(filepath, filename);
    % [~, name, ~] = fileparts(basePath);

    % Convert matrix to single precision
    if ~isa(matrix, 'single')
        matrix = single(matrix); % Convert to float32
    end

    % Replace NaN and Inf with -99999
    noDataValue = -99999;
    matrix(isnan(matrix)) = noDataValue;
    matrix(isinf(matrix)) = noDataValue;

    % Calculate cell size and validate extent
    nx = size(matrix, 2);
    ny = size(matrix, 1);
    cellSizeX = (xmax - xmin) / (nx - 1);
    cellSizeY = (ymax - ymin) / (ny - 1);

    % Create spatial referencing object
    R = georasterref('RasterSize', size(matrix), ...
                     'Latlim', [ymin, ymax], ...
                     'Lonlim', [xmin, xmax], ...
                     'ColumnsStartFrom', 'north');

    % Write the GeoTIFF file
    try
        geotiffwrite(basePath, matrix, R, 'CoordRefSysCode', 'EPSG:4326'); % WGS 84
        disp(['GeoTIFF saved to: ', basePath]);
    catch ME
        disp(['Error saving GeoTIFF: ', ME.message]);
        return;
    end

    disp(['Metadata files saved to: ', basePath, '.xml and .gi']);
end