function [array_pairwise_distance] = compute_pairwise_fingerprint_distance_sibling_pairs(sub_id_set,pair_id_set,mat_compact_fingerprint)

%%
% Summary:
%         1. MATLAB function to compute Eudlidean distance between 
%            pairs of compact fingerprints for a given twin/sibling type
%
%%
% Function Parameters:
%         Input:
%               1. sub_id_set: array containing the twin/sib id set
%               2. pair_id_set: array containing the twin/sib pair id set
%               3. mat_compact_fingerprint: matrix containing compact
%               fingerprints for each subject (Nsub x num_spect_component)
%         Output:
%               1. array_pairwise_distance: array containing Euclidean
%               distance between compact fingerprints of pairs of
%               twins/siblings  (size: length(sub_id_set) x 1 )
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


array_pairwise_distance = sqrt( sum((mat_compact_fingerprint(sub_id_set(:),:) - mat_compact_fingerprint(pair_id_set(:),:)).^2,2));

end