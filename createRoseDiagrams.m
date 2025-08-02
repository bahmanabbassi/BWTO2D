function [orientations360, numOrientations] = createRoseDiagrams(BW)
%createRoseDiagrams Calculates the orientation of line segments in a binary image.
%
% Description:
%   This function takes a binary image, assumed to contain skeletonized
%   lineaments (e.g., faults, fractures), and calculates the orientation of each
%   line segment. It fits a straight line to each segment, computes its angle,
%   and prepares the data for plotting a standard rose diagram by mirroring the
%   orientations across 180 degrees.
%
% Developed by:
%   Bahman Abbassi (bahman.abbassi@uqat.ca)
%   University of Quebec (UQAT)
%
% Input:
%   BW - A 2D binary matrix where non-zero pixels represent lineaments.
%
% Outputs:
%   orientations360 - A column vector of orientations in degrees [0-360].
%                     Each line segment is represented twice (angle and angle+180)
%                     to create a symmetric rose diagram.
%   numOrientations - The total number of unique line segments found.
%
% Example:
%   % 1. Create a simple binary image with a diagonal line
%   my_image = false(100, 100);
%   for i = 10:90, my_image(i, i) = true; end
%
%   % 2. Calculate the orientations
%   [angles, count] = createRoseDiagrams(my_image);
%
%   % 3. Display the results
%   fprintf('Found %d line segment(s).\n', count);
%   % >> Found 1 line segment(s).
%   % The output 'angles' will contain [45; 225].
%
%   % 4. Plot the rose diagram (requires Mapping Toolbox or compatible function)
%   % figure;
%   % rose(deg2rad(angles));
%   % title('Orientation Rose Diagram');


    % --- Step 1: Pre-process the Image ---

    % The input BW is assumed to be a skeletonized binary image.
    skeletonImage = BW;

    % Define a pixel width for cropping the border. This helps remove
    % lineaments that touch the edge of the image, which may be incomplete.
    cropPixels = 25;

    % Set the border pixels to 1 (white).
    % NOTE: This assumes lineaments are 0 (black) and the background is 1 (white).
    % If lineaments are 1, this operation may connect them at the border.
    skeletonImage(1:cropPixels, :) = 1;
    skeletonImage(end-cropPixels+1:end, :) = 1;
    skeletonImage(:, 1:cropPixels) = 1;
    skeletonImage(:, end-cropPixels+1:end) = 1;

    % --- Step 2: Identify and Analyze Line Segments ---

    % Label each distinct (connected) line segment in the image.
    % 'num' will be the total count of separate segments.
    [labeledImage, num] = bwlabel(skeletonImage);

    % Initialize an empty array to store the orientation of each segment.
    orientations = [];
    numSegments = num; % Store the total number of segments found.

    % Loop through each labeled segment to calculate its orientation.
    for k = 1:num
        % Find the (row, col) pixel coordinates for the current segment.
        [rows, cols] = find(labeledImage == k);

        % Proceed only if the segment is a line (more than one pixel).
        if length(rows) > 1
            % Fit a 1st-degree polynomial (a straight line) to the pixels.
            coeffs = polyfit(cols, rows, 1);

            % Calculate the angle from the line's slope (coeffs(1)).
            angle = atan(coeffs(1)) * (180 / pi); % Convert from radians to degrees.

            % Adjust the angle to the standard orientation range of [0, 180).
            if angle < 0
                angle = angle + 180;
            end

            % Add the calculated orientation to the list.
            orientations = [orientations; angle];
        end
    end

    % --- Step 3: Prepare Data for Rose Diagram ---

    % A rose diagram is symmetric. To create this, duplicate all orientations
    % and add 180 degrees to the copies. This represents each line by two
    % opposite vectors (e.g., a line at 20° is also shown at 200°).
    orientations360 = [orientations; orientations + 180];

    % Get the final count of unique orientations (line segments).
    numOrientations = length(orientations);
end