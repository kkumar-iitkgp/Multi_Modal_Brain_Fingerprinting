function [pval_ttest2,pval_ranksum] = compute_pval_ttest2_ranksum(data_array_1,data_array_2)

%%
% Summary:
%         1. MATLAB function to compute p-value between two distributions
%         2. Compute ttest2 and ranksum
%         3. Match sample size (default) using random sampling
%         4. Note: Sample size match can also be performed using a simple bipartite
%         matching based on age or other parameter
%         (requires age information and can be done as preprocessing)
%
%%
% Function Parameters:
%         Input:
%               1. data_array_1: data from distribution 1
%               2. data_array_2: data from distribution 2
%         Output:
%               1. pval_ttest2: p-value for ttest2 (unequal variance)
%               2. pval_ranksum: p-value for ranksum
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

rand_seed = 1;

if(length(data_array_1) == length(data_array_2))  % sample size matched
    temp_array_1 = data_array_1;
    temp_array_2 = data_array_2;
elseif(length(data_array_1) > length(data_array_2))
    rng(rand_seed);
    temp_rand_perm_set = randperm(length(data_array_1),length(data_array_2));
    
    temp_array_1 = data_array_1(temp_rand_perm_set(:));
    temp_array_2 = data_array_2;
else
    rng(rand_seed);
    temp_rand_perm_set = randperm(length(data_array_2),length(data_array_1));
    
    temp_array_1 = data_array_1;
    temp_array_2 = data_array_2(temp_rand_perm_set(:));
end


[pval_ranksum,~]=ranksum(temp_array_1,temp_array_2);
[~,pval_ttest2] = ttest2(temp_array_1,temp_array_2,'Vartype','unequal');

end