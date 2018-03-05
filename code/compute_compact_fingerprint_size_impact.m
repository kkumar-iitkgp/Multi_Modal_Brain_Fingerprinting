function cell_compact_fingerprint = compute_compact_fingerprint_size_impact(NxN_data_matrix,array_num_spect_component)

%%
% Summary:
%         1. MATLAB function to compute compact fingerprint for a range of 
%            num of spectral components (size imapct analysis)
%            for a given subject proximity graph (Nsub x Nsub)
%         2. Embedding is performed using Laplacian Eigenmaps approach
%         3. Note: Each compact fingerprint will be converted to z-score (0 mean
%         and unit variance)
%         4. This script is useful for testing the impact of number of
%         eigen vectors used for compact fingerprint
%
%%
% Function Parameters:
%         Input:
%               1. NxN_data_matrix: subject proximity graph ( Nsub x Nsub)
%               2. array_num_spect_component: array for the set of values of 
%                  the number of components in the compact fingerprint fingerprint
%         Output:
%               1. cell_compact_fingerprint: cell array containing
%               matrices of compact fingerprints (Nsub x num_spect_component)
%               correpoding to the array_num_spect_component
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

W_mat = NxN_data_matrix ;
row_sum = sum(W_mat,2);

D_inv_sqrt_mat = diag(1./sqrt(row_sum)) ;     
W_mat_norm = D_inv_sqrt_mat*W_mat*D_inv_sqrt_mat ; %  The effect of this normalization is to divide W_mat each element of by the square root of the row and column sum at that location.

% ensure symmetry
W_mat_norm = (W_mat_norm + W_mat_norm')/2 ;

% compute svd
[U,S,~] = svd(W_mat_norm); %#ok<ASGLU>


temp_d = diag(D_inv_sqrt_mat);
temp_Sqrt_Djj_mat = repmat(temp_d,1,size(W_mat,2));


cell_compact_fingerprint =  cell(length(array_num_spect_component),1);

for loop_i = 1:length(array_num_spect_component)
    % each component is normalized by the inverse of sqrt of row sum
    Temp_embed_mat = temp_Sqrt_Djj_mat(:,2:(array_num_spect_component(loop_i)+1)).*U(:,2:(array_num_spect_component(loop_i)+1));
    Z_test = (zscore(Temp_embed_mat'))';
    
    % compact fingerprints (Nsub x num_spect_component)
    cell_compact_fingerprint{loop_i,1} = Z_test ;
end
    
end