# Sex classification using long-term temporal dependence of rs-fMRI time series

This repository contains all of the MATLAB and Jupyter Notebook based implementation of the work described in: 
Dhamala E., Jamison, K., Sabuncu, M.R., Kuceyeski, A. (2019). Sex classification using long-term temporal dependence of resting state fMRI time series.   
  
Keywords: fMRI, resting-state fMRI, hurst exponent, temporal dependence, temporal dynamics, sex classification, sex differences  
  
This implementation is catered to the HCP-S1200 dataset where we focus on:  
(1) computation of hurst exponent (measure of long-term temporal dependence) from resting-state fMRI data using multiple parcellations (atlases) (MATLAB)  
(2) computation of subject-specific volume using mulitple parcellations (MATALB)  
(2) statistical analysis of sex differences at an ROI-level using different parcellations (MATLAB)  
(3) machine learning sex classification using hurst exponent and volume data using different parcellations (Jupyter Notebook - Python 3)  
  
Description of each MATLAB file:   
he_all_atlases - computes and saves the atlas-specific HE for each subject - uses the subject-specific atlas  
get_volume - computes and saves the atlas-specific regional volumes for each subject  
hurst_analysis - function that generates mean HE across all segments within a scan, and all scans, and analyses group-level sex differences in mean HE; calls the following functions: he_mean_across_scans, he_mean_across_segments, sex_analysis, p_correction  
he_mean_across_scans - computes the atlas-specific mean HE for each subject across the 4 scans  
he_mean_across_segments - compute the the atlas-specific mean HE for each subject across the segments within a scan  
sex_analysis - computes the group-level sex differences in mean HE  
p_correction - performs FDR correction using the Bonferonni-Hochberg method on the p-values generated by sex_analysis; correction is performed across all atlases  
  
Description of each Jupyter Notebook:   
he_sexclf_traintestsplit - generates n sets of indices for training/testing splits so they can be kept consistent over all models (i.e. across different atlases, for the HE-based vs. volume-based predictions)  
he_sexclf_permutations - performs n permutations of sex classification using a SVM with nested cross validation, uses the train/test indices generated by he_sexclf_traintestsplit - saves the optimised hyperparameter, prediction probabilities for the test subjects, accuracy, AUC, and feature importance for each permutation  





