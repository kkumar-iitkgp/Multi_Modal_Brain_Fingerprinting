# Multi_Modal_Brain_Fingerprinting

The code in this repository reproduces most of the analysis reported in "Multi-modal brain fingerprinting: a manifold approximation based framework".

Preprocessing required:
        * REQUIREMENTS: 
        1. report.txt files for each modality (see below)
        2. twin/sibling pair ids (numeric id's between 1 to number of
        subjects/images based on the subject ordering/reading list) for
        example: MZ sub id set and MZ pair id set etc.. 
        3. Note: The randomized sample includes pairs of MZ; DZ; 
          Full-Siblings (FS); Maternal Half Siblings (MHS) and 
          Paternal Half Siblings (PHS) 

        NOTE: a sample set of report files and twin/sibling pair ids
        (randomized) have been included in data folder for illustration
        purposes

%
Brief Overview of pre-processing:
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

%
 ADDITIONAL processing/files:
         1. feature match visualization requires 3D SIFT feature files
         (not included here)
         2. Significance values (-log10 pvalues) in paper are based on
         sample size matches using bipartite matching (age-based);
         However, sample illustration uses random sample size matching
         3. Similarly hemisphere asymmetry analysis and 
         pairwise feature correspondence analysis requires restricted data
         and thus have been excluded here      


%
Reference: 
          Multi-modal brain fingerprinting: a manifold approximation based framework
Authors: 
         Kuldeep Kumar (kkumar@livia.etsmtl.ca), 
         Laurent Chauvin
         Matthew Toews (Matthew.Toews@etsmtl.ca) 
         Olivier Colliot and 
         Christian Desrosiers (christian.desrosiers@etsmtl.ca)
    
LIVIA, ETS Montreal, Canada
January 2018
