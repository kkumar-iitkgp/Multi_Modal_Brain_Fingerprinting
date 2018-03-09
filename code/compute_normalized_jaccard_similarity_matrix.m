function NxN_data_matrix = compute_normalized_jaccard_similarity_matrix(cell_feature_count_data)

%%
% Summary:
%         1. MATLAB function to compute Normalized Jaccard Similarity
%         matrix
%         2. Note: make the diagonal values to be explicitly 1
%
%%
% Function Parameters:
%         Input:
%               1. cell_feature_count_data: cell array containing the
%               feature match count data for one or more modalities
%         Output:
%               1. NxN_data_matrix: subject proximity graph ( Nsub x Nsub )
%
%%
% Reference: 
%           Multi-modal brain fingerprinting: a manifold approximation based framework
% Authors: 
%          Kuldeep Kumar (kkumar@livia.etsmtl.ca), 
%          Laurent Chauvin
%          Matthew Toews (Matthew.Toews@etsmtl.ca) 
%          Olivier Colliot and 
%          Christian Desrosiers (christian.desrosiers@etsmtl.ca)
%     
% LIVIA, ETS Montreal, Canada
% January 2018
%
%%


% Add feature match count for multiple modalities
temp_mat = cell_feature_count_data{1,1};
for loop_i=2:length(cell_feature_count_data)
    temp_mat = temp_mat + cell_feature_count_data{loop_i,1};
end

% compute Normalized Jaccard Similarity
temp_feat_count_array = temp_mat(1,:);
temp_norm_mat = temp_mat./(repmat(temp_feat_count_array,size(temp_mat,1),1)+ repmat([ 1; temp_feat_count_array(:);],1,size(temp_mat,2)) - temp_mat );

% Normalized Jaccard Similarity Matrix
NxN_data_matrix = temp_norm_mat(2:end,:);

NxN_data_matrix(1:(size(NxN_data_matrix,1)+1):end) = 1;

% make it symmetric (explicitly)
NxN_data_matrix = (NxN_data_matrix +NxN_data_matrix')/2 ;

end