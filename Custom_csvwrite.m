function Custom_csvwrite(defaultFilename, matrix, ExpN)
    % Custom_csvwrite saves a matrix to a CSV file with a specified numerical suffix
    %
    % Inputs:
    %   defaultFilename - Default name for the CSV file (e.g., 'Selected.csv')
    %   matrix - The matrix to be saved
    %   ExpN - Numerical suffix to append to the file name
    %
    % Example usage:
    %   Custom_csvwrite('Selected.csv', rand(5,5), 3);
    % AUTHOR:
    % Bahman Abbassi, University of Quebec (UQAT)
    % Email: bahman.abbassi@uqat.ca

    % Check if the input matrix is valid
    if ~ismatrix(matrix)
        error('Input must be a matrix.');
    end

    % Check if ExpN is a valid positive integer
    if ~isnumeric(ExpN) || ExpN < 0 || floor(ExpN) ~= ExpN
        error('ExpN must be a non-negative integer.');
    end

    % Open a file save dialog with the default filename
    [filename, filepath] = uiputfile('*.csv', 'Save Matrix as CSV File', defaultFilename);

    % If the user cancels the dialog, exit the function
    if isequal(filename, 0) || isequal(filepath, 0)
        disp('File save canceled by user.');
        return;
    end

    % Extract the file name and extension
    [~, baseName, ext] = fileparts(filename);

    % Ensure the extension is .csv
    if isempty(ext)
        ext = '.csv';
    elseif ~strcmpi(ext, '.csv')
        error('File extension must be .csv');
    end

    % Construct the file name with the given suffix
    finalFileName = fullfile(filepath, [baseName, '_', num2str(ExpN), ext]);

    % Save the matrix as a CSV file
    writematrix(matrix, finalFileName);

    % Confirm the save
    disp(['Matrix saved to: ', finalFileName]);
end
