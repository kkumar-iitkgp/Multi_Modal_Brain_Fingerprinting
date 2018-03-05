function [nearest_neighbor_matrix] = compute_nearest_neighbor_matrix(NxN_data_matrix,knn_value,flag_distance_or_similarity)

%%
% Summary:
%         1. MATLAB function to compute nearest neighbor matrix based on
%         data matrix
%         2. An index odered matrix is generated 
%         3. Data matrix can have distance values ( e.g. pairwise Euclidean distance)
%           or similarity values (e.g. pairwise correlation or normalzied Jaccard similarity) 
%         4. Note: self similarity will be excluded
%
%%
% Function Parameters:
%         Input:
%               1. NxN_data_matrix: subject proximity graph ( Nsub x Nsub)
%               2. knn_value: number of nearest neighbors required (default
%               = Nsub - 1)
%               3. flag_distance_or_similarity: 0-distance; 1-similarity
%         Output:
%               1. nearest_neighbor_matrix: nearest neighbor index ordered
%               matrix ( Nsub x knn_value )
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


% check knn_value to be less than available data
temp_max_knn_value = size(NxN_data_matrix,2);
if(knn_value > temp_max_knn_value)
    knn_value = temp_max_knn_value -1;
end

if(flag_distance_or_similarity)
        [~,temp_sort_index_set] = sort(NxN_data_matrix,2,'descend');
else
        [~,temp_sort_index_set] = sort(NxN_data_matrix,2,'ascend');     
end

% check if we need to remove self
% if the first element correspond to diagonal indices
if(isempty(find(temp_sort_index_set(:,1)-(1:temp_max_knn_value)', 1)))
    start_index = 2;
     end_index = knn_value+1 ;
else
     start_index = 1;
     end_index = knn_value ;
end
     
nearest_neighbor_matrix = temp_sort_index_set(:,start_index:end_index) ;

end