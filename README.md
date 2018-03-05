## Multi_Modal_Brain_Fingerprinting

The code in this repository reproduces most of the analysis reported in "Multi-modal brain fingerprinting: a manifold approximation based framework".

__________________________________________________________________
#Following analysis can be performed:
                
        1. Comapct fingerprint generation using Laplacian Eigenmaps based embedding.
        
        2. Analysis of impact of number of eigen vectors used for compact fingerprint generation.
        
        3. Rank reterival analysis for twin/sibling types using randomized subject IDs and sample subject proximity graphs. Following measures can be generated: mean average precision (MAP), mean recall@10.
        
        4. Illustration for T1w, FA, T1w+FA, T1w+rfMRI, FA+rfMRI, and T1w+T2w+FA+rfMRI based fingerprints.
        
        5. Relative informativeness of modality: a modality vs modality comparison for the task of twin/sibling identification.

        NOTE: a sample set of report files and twin/sibling pair ids
        (randomized) have been included in data folder for illustration
        purposes. 
__________________________________________________________________
#Pre-processing Steps Summary:

        (* REQUIREMENTS for complete analysis, starting from NiFTI volume files*)
        
        1. Obtain Bag of feature (BoF) files for each image for a given modality.
        
        2. Obtain report.txt file containing approximate feature correspondence counts, for each modality (see below). This will be converted to Normalized Jaccard Similarity measure and used as subject proximity graph (manifold approximation).
        
        3. In case of rfMRI or FreeSurfer based measures or using full image as fingerprint: compute pairwise (images x images) similarity for example, using pearson correlation. This images x images matrix will be used as subject proximity graph (manifold approximation).
        
        4. For modality combinations two cases arise: a) combining feature correspondence counts: add the two/more modality report.txt files, and then compute Normalized Jaccard Similarity. b) For combining Normalized Jaccard Similarity with correlation based similarity matrix, use linear combination (it can be further optimized, future work). 
        
        5. twin/sibling pair ids are numeric indices between 1 to number of
        subjects/images based on the subject ordering/reading list, for
        example: MZ sub id set (1, 2, 5, ..) and MZ pair id set (4, 8, 10,..) etc... Cross check the ids and the file read list for report.txt.
        
        
__________________________________________________________________
%
#Brief Overview of Bag of Feature (BoF) pre-processing:

        1. We assume that Bag-of-features (BoF) have been generated and 
           approximate matching has been performed using 
           the code available at http://www.matthewtoews.com/ 
           THE PRE-PROCESSING involves running:
           
           i. featExtract.exe for each image (.nii files): output will be
           .key files containing 3D SIFT features
           
           ii. featMatchMultiple -t 2 -r- -n 20 *.key (for all .key files
           of a given modality): output will be multiple files including
           report.txt
           
         2. report.txt generated for each modality will act as a input
         for the following analysis
         
         3. We have included sample report files (Random permutation
         applied on our set)
         
         4. Note: cross check the file order read in featMatchMultiple.exe 
          and the sibling/twin pair IDs

__________________________________________________________________
%
 #ADDITIONAL processing/files:
 
         1. Pairwise feature match visualization requires 3D SIFT feature files.
         (not included here)
         
         2. Significance values (-log10 pvalues) in paper are based on
         sample size matches using bipartite matching (age-based);         
         However, sample illustration uses random sample size matching.
         
         3. Similarly hemisphere asymmetry analysis and 
         pairwise feature correspondence analysis requires restricted data
         and thus have been excluded here.  


__________________________________________________________________
%
#Reference: 

          Multi-modal brain fingerprinting: a manifold approximation based framework
          
#Authors: 

         Kuldeep Kumar (kkumar@livia.etsmtl.ca), 
         
         Laurent Chauvin,
         
         Matthew Toews (Matthew.Toews@etsmtl.ca),
         
         Olivier Colliot, and 
         
         Christian Desrosiers (christian.desrosiers@etsmtl.ca)
    
#LIVIA, ETS Montreal, Canada

#January 2018
