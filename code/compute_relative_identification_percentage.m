function [array_relative_identification_percentage] = compute_relative_identification_percentage(sub_id_set,pair_id_set,nearest_neighbor_matrix_mod1,nearest_neighbor_matrix_mod2)

%%
% Summary:
%         1. MATLAB function to compute RELATIVE identification percentage
%         2. Identification is considered a success if the twin/sib pair is
%         identified within k nearest neighbors (for example k = 10 ); this
%         parameter is by default the number of columns of
%         nearest_neighbor_matrices
%         3. Relative percentages indicate insatnces identified by 
%                i: both modalities, 
%               ii: first modality only, 
%              iii: second modality only,
%           and iv. none of the modalities
%         4. Note: it's relative percentage, 
%         as such depends on the pair of modalities being compared and 
%         on the number of neighbors considered for identification task 
%         5. Note: Identification is considered over all subjects as
%         nearest neighbor matrix lists all Nsub - 1 neighbors first
%
%%
% Function Parameters:
%         Input:
%               1. sub_id_set: array containing the twin/sib id set
%               2. pair_id_set: array containing the twin/sib pair id set
%               3. nearest_neighbor_matrix_mod1: nearest neighbor index ordered
%               matrix ( Nsub x knn_value ) for modality 1
%               4. nearest_neighbor_matrix_mod2: nearest neighbor index ordered
%               matrix ( Nsub x knn_value ) for modality 2
%         Output:
%               1. array_relative_identification_percentage: array
%               containing the relative identification percentages (size: 4
%               x 1): with percentages for iderntification by i. both ii.
%               mod1 only iii. modality 2 only iv. none
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


array_relative_identification_percentage = zeros(4,1);  % value 1: both, 2: first only, 3: second only, 4: none
temp_num_subs = length(sub_id_set);
num_idfn_tasks = temp_num_subs ;


idfn_task_recall_array_1 = zeros(num_idfn_tasks,1);
idfn_task_recall_array_2 = zeros(num_idfn_tasks,1);

% here we assume over all knn_value
% but it can be modified to consider the order also
for i=1:num_idfn_tasks
    if(~isempty(find(nearest_neighbor_matrix_mod1(sub_id_set(i),:) == pair_id_set(i),1)))
        idfn_task_recall_array_1(i) = 1;
    end
    if(~isempty(find(nearest_neighbor_matrix_mod2(sub_id_set(i),:) == pair_id_set(i),1)))
        idfn_task_recall_array_2(i) = 1;
    end
end

a1 = find(idfn_task_recall_array_1 ==1);
b1 = find(idfn_task_recall_array_2 ==1);

% both
array_relative_identification_percentage(1) = length(find(idfn_task_recall_array_1 & idfn_task_recall_array_2))/num_idfn_tasks; 
% 1st only
array_relative_identification_percentage(2) = length(setdiff(a1,b1))/num_idfn_tasks;
% second only
array_relative_identification_percentage(3) = length(setdiff(b1,a1))/num_idfn_tasks;
% none
array_relative_identification_percentage(4) = 1 - (length(find(idfn_task_recall_array_1 | idfn_task_recall_array_2))/num_idfn_tasks);

% Convert to percentage
array_relative_identification_percentage = array_relative_identification_percentage*100 ;

end