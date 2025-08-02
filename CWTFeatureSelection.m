function [Selected_Features_3D, selected_indices_final, debug_info] = CWTFeatureSelection(CWT_Coefficients_3D, TargetNumFeatures, RudimentaryEdgePercentile, MapCorrelationWeight)
%CWTFeatureSelection Selects a distinct subset of CWT coefficient maps.
%
% Description:
%   This function selects a specified number of Continuous Wavelet Transform (CWT)
%   maps from a larger 3D set. The selection algorithm iteratively builds a
%   subset of features by prioritizing maps that are both high-variance and
%   dissimilar to maps already chosen. Dissimilarity is a weighted combination
%   of low spatial correlation and low overlap of their most prominent edge features.
%   This ensures a rich, non-redundant set of features.
%
% Developer:
%   Bahman Abbassi (bahman.abbassi@uqat.ca)
%
% Inputs:
%   CWT_Coefficients_3D (double)       - 3D matrix of CWT coefficients (height x width x num_maps).
%   TargetNumFeatures (integer)        - The desired number of maps to select.
%   RudimentaryEdgePercentile (double) - Percentile (0-100) used to define "edges" by
%                                        thresholding the positive CWT coefficients. E.g., 95 means
%                                        pixels with coefficients in the top 5% are considered edges.
%   MapCorrelationWeight (double)      - A weight (0 to 1) that balances the two dissimilarity criteria.
%                                        * Weight for map correlation dissimilarity = MapCorrelationWeight
%                                        * Weight for edge overlap dissimilarity = 1 - MapCorrelationWeight
%
% Outputs:
%   Selected_Features_3D (double)    - 3D matrix containing only the selected CWT maps.
%   selected_indices_final (integer) - A row vector of the indices of the selected maps.
%   debug_info (struct)              - A struct array providing detailed metrics for each selection step,
%                                      useful for tuning and analysis.
%
% Example:
%   % 1. Create a sample 3D CWT coefficient matrix (e.g., 100x100 with 50 maps)
%   num_maps = 50;
%   base_data = peaks(100);
%   CWT_Coefficients_3D = zeros(100, 100, num_maps);
%   for i = 1:num_maps
%       % Create diverse maps by rotating and adding noise
%       CWT_Coefficients_3D(:,:,i) = imrotate(base_data, i*7, 'bilinear', 'crop') + randn(100)*0.1;
%   end
%
%   % 2. Define selection parameters
%   TargetNumFeatures = 10;
%   RudimentaryEdgePercentile = 98; % Use top 2% of positive coeffs as edges
%   MapCorrelationWeight = 0.5;      % Give equal weight to edge and correlation dissimilarity
%
%   % 3. Run the feature selection
%   [Selected_Maps, Selected_Indices, Debug_Report] = CWTFeatureSelection(...
%       CWT_Coefficients_3D, TargetNumFeatures, RudimentaryEdgePercentile, MapCorrelationWeight);
%
%   % 4. Display results
%   fprintf('Selected %d features.\n', size(Selected_Maps, 3));
%   fprintf('Selected Indices: %s\n', num2str(Selected_Indices));
%   disp('Debug Information for the first selected feature:');
%   disp(Debug_Report(1));


% --- Argument and State Initialization ---

    % Get dimensions of the input CWT data
    [height, width, num_total_maps] = size(CWT_Coefficients_3D);
    
    % --- Input Validation ---
    if nargin < 4
        error('All four input arguments: CWT_Coefficients_3D, TargetNumFeatures, RudimentaryEdgePercentile, and MapCorrelationWeight must be provided.');
    end
    if TargetNumFeatures <= 0, error('TargetNumFeatures must be a positive integer.'); end
    if RudimentaryEdgePercentile < 0 || RudimentaryEdgePercentile > 100
        error('RudimentaryEdgePercentile must be between 0 and 100.');
    end
    if MapCorrelationWeight < 0 || MapCorrelationWeight > 1
        error('MapCorrelationWeight must be between 0 and 1.');
    end
    % Handle cases where requested features exceed available maps
    if TargetNumFeatures > num_total_maps
        warning('TargetNumFeatures (%d) is > available maps (%d). Returning all maps.', TargetNumFeatures, num_total_maps);
        Selected_Features_3D = CWT_Coefficients_3D; selected_indices_final = 1:num_total_maps; debug_info = []; return;
    end
    % Handle empty input
    if num_total_maps == 0, Selected_Features_3D = []; selected_indices_final = []; debug_info = []; warning('Input CWT_Coefficients_3D is empty.'); return; end

    % Define the weight for edge dissimilarity based on the correlation weight
    EDGE_DISSIMILARITY_WEIGHT = 1 - MapCorrelationWeight;
    
    % --- Pre-computation for Efficiency ---
    
    % Pre-calculate variance and create flattened (vector) versions of all maps
    % This avoids redundant calculations inside the main selection loop.
    original_map_variances = zeros(1, num_total_maps);
    flattened_cwt_maps = zeros(height * width, num_total_maps);
    for i = 1:num_total_maps
        current_map_data = CWT_Coefficients_3D(:, :, i);
        original_map_variances(i) = var(current_map_data(:), 1); % '1' for normalization by N
        flattened_cwt_maps(:, i) = current_map_data(:);
    end

    % Handle edge case where all maps have no variance (e.g., are constant)
    if all(original_map_variances < eps)
        warning('All maps have near-zero variance. Selecting first %d maps.', TargetNumFeatures);
        Selected_Features_3D = CWT_Coefficients_3D(:,:,1:TargetNumFeatures); selected_indices_final = 1:TargetNumFeatures; debug_info = []; return;
    end
    
    % --- Initialize Selection Loop Variables ---

    % Array to store the indices of the maps selected at each step
    selected_indices_final = zeros(1, TargetNumFeatures);
    % Caches to store data for already-selected maps to speed up comparisons
    selected_rudimentary_edges_cache = false(height, width, TargetNumFeatures);
    selected_flattened_cwt_cache = zeros(height*width, TargetNumFeatures);

    % Parameters for the internal edge extraction helper function
    edgeParams.method = 'percentile_positive';
    edgeParams.percentile_value = RudimentaryEdgePercentile;
    
    % Pre-allocate a struct array to hold debugging information for each selection
    debug_info = repmat(struct('selected_idx',0,'variance',0,'max_edge_overlap',0,'max_cwt_corr',0,'score',0), TargetNumFeatures, 1);
    
    % Create a temporary copy of variances to be modified during selection
    temp_map_variances = original_map_variances;

    % --- Iteratively Select TargetNumFeatures Maps ---
    for k = 1:TargetNumFeatures
        best_candidate_idx = -1;
        best_candidate_score = -inf;

        % --- Step 1: Select the first map ---
        % The first map is chosen based purely on having the highest variance.
        if k == 1
            [max_var_val, best_candidate_idx] = max(temp_map_variances);
             % Check if any valid map was found
             if isinf(max_var_val) && max_var_val < 0
                warning('No valid map with positive variance found for first selection. Aborting.');
                Selected_Features_3D = []; selected_indices_final = []; debug_info = []; return;
            end
            % Mark this variance as -inf so it won't be picked again
            temp_map_variances(best_candidate_idx) = -inf;
            
            % Log debug information for the first selection
            debug_info(k).selected_idx = best_candidate_idx;
            debug_info(k).variance = max_var_val;
            debug_info(k).max_edge_overlap = 0; % No previous maps to compare against
            debug_info(k).max_cwt_corr = 0;     % No previous maps to compare against
            debug_info(k).score = max_var_val;  % The score is just its variance
        
        % --- Steps 2 to N: Select subsequent maps ---
        else
            % Identify unselected maps to find the max variance for normalization
            unselected_indices_mask = true(1, num_total_maps);
            unselected_indices_mask(selected_indices_final(1:(k-1))) = false;
            
            % Find the maximum variance among the remaining, valid candidates
            available_variances = original_map_variances(unselected_indices_mask & (original_map_variances >= 0) );
            if isempty(available_variances) || all(available_variances < eps)
                 max_overall_variance_for_norm = 1; % Fallback normalization factor
            else
                 max_overall_variance_for_norm = max(available_variances);
                 if max_overall_variance_for_norm < eps, max_overall_variance_for_norm = 1; end
            end

            % Iterate through all maps to find the best *next* candidate
            for candidate_idx = 1:num_total_maps
                % Skip maps that have already been selected
                if ismember(candidate_idx, selected_indices_final(1:(k-1))), continue; end
                
                var_candidate = original_map_variances(candidate_idx);
                % Skip maps with no significant variance
                if var_candidate < eps, continue; end
                
                % Get candidate's flattened data and extract its edge map
                candidate_cwt_map_flat = flattened_cwt_maps(:, candidate_idx);
                candidate_edges = ExtractRudimentaryEdgesFixedInternal(CWT_Coefficients_3D(:, :, candidate_idx), edgeParams);
                
                % --- Calculate Dissimilarity against all previously selected maps ---
                num_already_selected = k - 1;
                edge_overlaps = zeros(1, num_already_selected);
                cwt_corrs = zeros(1, num_already_selected);
                
                for j = 1:num_already_selected
                    % A) Calculate Edge Overlap (Jaccard Index)
                    selected_edge_map_j = selected_rudimentary_edges_cache(:, :, j);
                    if ~any(candidate_edges(:)) || ~any(selected_edge_map_j(:))
                        edge_overlaps(j) = 0; % No edges means no overlap
                    else
                        intersection_val = sum(candidate_edges & selected_edge_map_j, 'all');
                        union_val = sum(candidate_edges | selected_edge_map_j, 'all');
                        if union_val == 0, edge_overlaps(j) = 0; else, edge_overlaps(j) = intersection_val / union_val; end
                    end
                    
                    % B) Calculate CWT Map Correlation (Absolute Pearson Correlation)
                    selected_cwt_flat_j = selected_flattened_cwt_cache(:, j);
                    % Handle zero-variance cases for correlation
                    if std(candidate_cwt_map_flat,"omitnan") < eps || std(selected_cwt_flat_j,"omitnan") < eps
                        % If both are constant and equal, correlation is 1; otherwise 0.
                        if std(candidate_cwt_map_flat,"omitnan") < eps && std(selected_cwt_flat_j,"omitnan") < eps && abs(mean(candidate_cwt_map_flat,"omitnan") - mean(selected_cwt_flat_j,"omitnan")) < eps
                            cwt_corrs(j) = 1; else, cwt_corrs(j) = 0; end
                    else
                        % Standard correlation calculation for maps with variance
                        common_valid_indices = ~isnan(candidate_cwt_map_flat) & ~isnan(selected_cwt_flat_j);
                        if sum(common_valid_indices) < 2 
                             cwt_corrs(j) = 0; % Not enough data to correlate
                        else
                            map_A_for_corr = candidate_cwt_map_flat(common_valid_indices);
                            map_B_for_corr = selected_cwt_flat_j(common_valid_indices);
                            % Check again for constant series after removing NaNs
                            if std(map_A_for_corr) < eps || std(map_B_for_corr) < eps 
                                cwt_corrs(j) = 0;
                            else
                                corr_matrix = corrcoef(map_A_for_corr, map_B_for_corr);
                                cwt_corrs(j) = abs(corr_matrix(1,2));
                            end
                        end
                    end
                end
                
                % Find the *maximum* similarity (worst-case overlap) to any already selected map
                max_edge_overlap_with_selected = max([0, edge_overlaps]);
                max_cwt_corr_with_selected = max([0, cwt_corrs]);
                
                % Convert similarity to dissimilarity (higher is better)
                edge_dissimilarity = 1 - max_edge_overlap_with_selected;
                map_corr_dissimilarity = 1 - max_cwt_corr_with_selected;
                
                % --- Calculate Final Score for the Candidate ---
                % Weighted average of the two dissimilarity types
                combined_dissimilarity = EDGE_DISSIMILARITY_WEIGHT * edge_dissimilarity + MapCorrelationWeight * map_corr_dissimilarity;
                
                % Normalize the candidate's variance against the max available variance
                norm_var_candidate = var_candidate / max_overall_variance_for_norm;
                
                % Final score combines variance and dissimilarity. We want to maximize this.
                current_score = norm_var_candidate * combined_dissimilarity;

                % If this candidate has the best score so far, update our choice
                if current_score > best_candidate_score
                    best_candidate_score = current_score;
                    best_candidate_idx = candidate_idx;
                    % Temporarily store debug info for the best candidate
                    debug_info(k).variance = var_candidate;
                    debug_info(k).max_edge_overlap = max_edge_overlap_with_selected;
                    debug_info(k).max_cwt_corr = max_cwt_corr_with_selected;
                    debug_info(k).score = current_score;
                end
            end % End of loop through candidate maps
            
            % --- Fallback Selection ---
            % If no suitable candidate was found (e.g., all remaining maps are identical),
            % select the one with the highest remaining variance as a fallback.
            if best_candidate_idx == -1
                warning('Could not find a distinct candidate for map %d. Filling with highest remaining unselected variance.', k);
                temp_variances_fallback = original_map_variances;
                temp_variances_fallback(selected_indices_final(1:(k-1))) = -inf;
                [fallback_var_val, best_candidate_idx] = max(temp_variances_fallback);
                
                % If no maps are left, exit gracefully
                if isinf(fallback_var_val) && fallback_var_val < 0
                    warning('No more valid maps to select for map %d. Returning %d maps.', k, k-1);
                    Selected_Features_3D = CWT_Coefficients_3D(:, :, selected_indices_final(1:(k-1)));
                    selected_indices_final = selected_indices_final(1:(k-1));
                    debug_info = debug_info(1:(k-1));
                    return; % Exit function early
                end
                % Update debug info for the fallback choice
                debug_info(k).variance = fallback_var_val;
                debug_info(k).max_edge_overlap = NaN;
                debug_info(k).max_cwt_corr = NaN;
                debug_info(k).score = fallback_var_val;
            end
            
            % Mark the chosen map so it is not considered in future variance checks
            temp_map_variances(best_candidate_idx) = -inf;
            % Finalize debug info for the selected map
            debug_info(k).selected_idx = best_candidate_idx;
        end % End of if/else for k=1 vs k>1
        
        % --- Update State with the Newly Selected Map ---
        
        % Add the winning index to our final list
        selected_indices_final(k) = best_candidate_idx;
        
        % Retrieve the selected map's data
        selected_map_data = CWT_Coefficients_3D(:, :, best_candidate_idx);
        
        % Store its edges and flattened data in the cache for the next iteration
        selected_rudimentary_edges_cache(:, :, k) = ExtractRudimentaryEdgesFixedInternal(selected_map_data, edgeParams);
        selected_flattened_cwt_cache(:, k) = flattened_cwt_maps(:, best_candidate_idx);
        
    end % End of main selection loop (for k)
    
    % --- Final Output Generation ---
    % Use the final list of indices to extract the selected maps from the original data
    Selected_Features_3D = CWT_Coefficients_3D(:, :, selected_indices_final);
end

%--------------------------------------------------------------------------
% Helper Function
%--------------------------------------------------------------------------
function edge_map = ExtractRudimentaryEdgesFixedInternal(cwt_map, params)
% Extracts a binary edge map by thresholding a CWT map.
% It identifies edges as the pixels with the highest positive coefficient values.
    [h, w] = size(cwt_map);
    edge_map = false(h, w);
    
    % Default parameters if not provided
    if ~isfield(params, 'method') || ~strcmpi(params.method, 'percentile_positive'), params.method = 'percentile_positive'; end
    if ~isfield(params, 'percentile_value'), params.percentile_value = 95; end
    
    % Return an empty edge map if the input map is blank or all NaN
    if sum(abs(cwt_map(:)),"omitnan") < eps || all(isnan(cwt_map(:))), return; end
    
    % Isolate positive coefficients, as they represent the "ridges" or features
    positive_coeffs = cwt_map(cwt_map > 0 & ~isnan(cwt_map));
    if isempty(positive_coeffs), return; end % No positive features found
    
    % Calculate the threshold using the specified percentile
    threshold = prctile(positive_coeffs, params.percentile_value);
    
    % --- Robustness checks for the calculated threshold ---
    if isempty(threshold) || ~isscalar(threshold) || isnan(threshold)
        if ~isempty(positive_coeffs), threshold = median(positive_coeffs);
            if isempty(threshold) || threshold <=0 || isnan(threshold), return; end
        else, return; end
    end
    % Handle cases where the percentile results in a near-zero threshold
    if threshold <= 1e-9 && any(positive_coeffs > 1e-9)
        min_positive_val = min(positive_coeffs(positive_coeffs > 1e-9));
        if ~isempty(min_positive_val), threshold = min_positive_val; else, return; end
    end
    
    % Create the final binary edge map
    edge_map = (cwt_map >= threshold) & (cwt_map > 0) & ~isnan(cwt_map);
end