function [mean_avg_precision,mean_recall_at_10,array_avg_precision,mat_recall_at_k,mat_precision_at_k] = compute_rank_reterival_measures(sub_id_set,pair_id_set,nearest_neighbor_matrix)

%%
% Summary:
%         1. MATLAB function to compute rank reterival measures:
%         mean_avg_precision (MAP); mean_recall_at_10;
%         array_avg_precision ; mat_recall_at_k; and mat_precision_at_k
%         2. Note: It's a single function to compute all rank reterival measures for a
%         given twin/sibling type and a given nearest_neighbor_matrix
%         3. Avg Prec = Sum_(k=1) ^ knn_value (P(k)*rel(k)) / Number of siblings
%           Where: 
%                   P(k) is precision at k
%                   rel(k) is an indicator function equaling 1 if the subject at rank k is a twin/sibling, zero otherwise
%           Essentialy sum of precision at locations where we find the twin/sibling pair match
%          4. Since we consider pairs of twin/sibling, Number of siblings
%          will be 1 by default.
%
%          For definitions See:
%          *https://en.wikipedia.org/wiki/Information_retrieval#Mean_average_precision
%          *https://en.wikipedia.org/wiki/Information_retrieval#Average_precision
%          *https://en.wikipedia.org/wiki/Information_retrieval#Precision_at_K
%
%%
% Function Parameters:
%         Input:
%               1. sub_id_set: array containing the twin/sib id set
%               2. pair_id_set: array containing the twin/sib pair id set
%               3. nearest_neighbor_matrix_mod1: nearest neighbor index ordered
%               matrix ( Nsub x knn_value ) for modality 1
%         Output:
%               1. mean_avg_precision: mean average precision value (mean
%               over average precision array); (size: 1 x 1 )
%               2. mean_recall_at_10: mean recall at k (=10) (size: 1 x 1)
%               3. array_avg_precision: array containing average precision
%                  value for each twin/sinling pair 
%                 (size: length(sub_id_set) x 1)
%               4. mat_recall_at_k: matrix containing recall at k values
%               for each twin/sinling pair 
%                 (size: length(sub_id_set) x number nearest neighbors)
%               5. mat_precision_at_k: matrix containing precision at k values
%               for each twin/sinling pair 
%                 (size: length(sub_id_set) x number nearest neighbors)
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


temp_num_sub = length(sub_id_set);   
max_nn_value = size(nearest_neighbor_matrix,2);  % maximum number of nearest neighbors

% get the subset nearest neighbor matrix
temp_knn_mat = nearest_neighbor_matrix(sub_id_set(:),:);

% initialize recall@k matrix, precision@k matrix, and avg precision array
mat_precision_at_k = zeros(temp_num_sub,max_nn_value);
mat_recall_at_k = zeros(temp_num_sub,max_nn_value);
array_avg_precision = zeros(temp_num_sub,1);
    
for loop_i=1:temp_num_sub

    % compute precision@k and recall@k
    for loop_k =1:max_nn_value
        mat_precision_at_k(loop_i,loop_k) = length(find(temp_knn_mat(loop_i,1:loop_k)==repmat(pair_id_set(loop_i),1,loop_k)))/loop_k ;
        mat_recall_at_k(loop_i,loop_k) = length(find(temp_knn_mat(loop_i,1:loop_k)==repmat(pair_id_set(loop_i),1,loop_k)))/(size(pair_id_set(loop_i),1)) ;   
    end

    % compute avg precision over all nearest neighbors
    temp_ind = find(temp_knn_mat(loop_i,1:max_nn_value)==pair_id_set(loop_i));
    if(~isempty(temp_ind))
        array_avg_precision(loop_i) = mat_precision_at_k(loop_i,temp_ind); % since number of sib/twin is 1 it will be a single value only 
    end

end
    
% compute mean average precision (MAP) 
mean_avg_precision = mean(array_avg_precision);

% compute mean recall@10
temp_mean_recall_array = mat_recall_at_k(:,10); % change it manually for the desired experiment
mean_recall_at_10 = mean(temp_mean_recall_array);

end