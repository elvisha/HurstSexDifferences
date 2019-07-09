function hurst_analysis(start_seg,end_seg)

%choose start and end segments to use when computing the mean of the hurst
%across segments within each scan
start_seg=start_seg;
end_seg=end_seg;

%compute hurst mean across segments from rfMRI_REST1_LR for all atlases
load('rfMRI_REST1_LR_hurst.mat')
he_aal_1LR = he_mean_across_segments(he_aal, start_seg, end_seg);
he_ez_1LR = he_mean_across_segments(he_ez, start_seg, end_seg);
he_ho_1LR = he_mean_across_segments(he_ho, start_seg, end_seg);
he_cc200_1LR = he_mean_across_segments(he_cc200, start_seg, end_seg);
he_cc400_1LR = he_mean_across_segments(he_cc400, start_seg, end_seg);
he_fs86_1LR = he_mean_across_segments(he_fs86, start_seg, end_seg);
he_tt_1LR = he_mean_across_segments(he_tt, start_seg, end_seg);
clear he_aal he_cc200 he_cc400 he_ez he_ho he_tt fmriname subj_str i;

%compute hurst mean across segments from rfMRI_REST1_RL for all atlases
load('rfMRI_REST1_RL_hurst.mat')
he_aal_1RL = he_mean_across_segments(he_aal, start_seg, end_seg);
he_ez_1RL = he_mean_across_segments(he_ez, start_seg, end_seg);
he_ho_1RL = he_mean_across_segments(he_ho, start_seg, end_seg);
he_cc200_1RL = he_mean_across_segments(he_cc200, start_seg, end_seg);
he_cc400_1RL = he_mean_across_segments(he_cc400, start_seg, end_seg);
he_fs86_1RL = he_mean_across_segments(he_fs86, start_seg, end_seg);
he_tt_1RL = he_mean_across_segments(he_tt, start_seg, end_seg);
clear he_aal he_cc200 he_cc400 he_ez he_ho he_tt fmriname subj_str i;

%compute hurst mean across segments from rfMRI_REST2_LR for all atlases
load('rfMRI_REST2_LR_hurst.mat')
he_aal_2LR = he_mean_across_segments(he_aal, start_seg, end_seg);
he_ez_2LR = he_mean_across_segments(he_ez, start_seg, end_seg);
he_ho_2LR = he_mean_across_segments(he_ho, start_seg, end_seg);
he_cc200_2LR = he_mean_across_segments(he_cc200, start_seg, end_seg);
he_cc400_2LR = he_mean_across_segments(he_cc400, start_seg, end_seg);
he_fs86_2LR = he_mean_across_segments(he_fs86, start_seg, end_seg);
he_tt_2LR = he_mean_across_segments(he_tt, start_seg, end_seg);
clear he_aal he_cc200 he_cc400 he_ez he_ho he_tt fmriname subj_str i;

%compute hurst mean across segments from rfMRI_REST2_RL for all atlases
load('rfMRI_REST2_RL_hurst.mat')
he_aal_2RL = he_mean_across_segments(he_aal, start_seg, end_seg);
he_ez_2RL = he_mean_across_segments(he_ez, start_seg, end_seg);
he_ho_2RL = he_mean_across_segments(he_ho, start_seg, end_seg);
he_cc200_2RL = he_mean_across_segments(he_cc200, start_seg, end_seg);
he_cc400_2RL = he_mean_across_segments(he_cc400, start_seg, end_seg);
he_fs86_2RL = he_mean_across_segments(he_fs86, start_seg, end_seg);
he_tt_2RL = he_mean_across_segments(he_tt, start_seg, end_seg);
clear he_aal he_cc200 he_cc400 he_ez he_ho he_tt fmriname subj_str i;

%compute hurst mean across all scans for all atlases
he_aal=he_mean_across_scans(he_aal_1LR,he_aal_1RL,he_aal_2LR,he_aal_2RL);
he_cc200=he_mean_across_scans(he_cc200_1LR,he_cc200_1RL,he_cc200_2LR,he_cc200_2RL);
he_cc400=he_mean_across_scans(he_cc400_1LR,he_cc400_1RL,he_cc400_2LR,he_cc400_2RL);
he_ho=he_mean_across_scans(he_ho_1LR,he_ho_1RL,he_ho_2LR,he_ho_2RL);
he_tt=he_mean_across_scans(he_tt_1LR,he_tt_1RL,he_tt_2LR,he_tt_2RL);
he_ez=he_mean_across_scans(he_ez_1LR,he_ez_1RL,he_ez_2LR,he_ez_2RL);
he_fs86=he_mean_across_scans(he_fs86_1LR,he_fs86_1RL,he_fs86_2LR,he_fs86_2RL);

%remove ROIs that are missing/not matching with label files
he_tt(:,67)=[];
he_tt(:,66)=[];
he_tt(:,63)=[];
he_ho(:,83)=[];

%load in subject sex data
subj_sex=load('subj_sex.txt');

%compute sex difference stats for each atlas
%generates t-stat and p-values (not corrected)
[aal_t,aal_p]=sex_analysis(he_aal,subj_sex);
[cc200_t,cc200_p]=sex_analysis(he_cc200,subj_sex);
[cc400_t,cc400_p]=sex_analysis(he_cc400,subj_sex);
[ez_t,ez_p]=sex_analysis(he_ez,subj_sex);
[tt_t,tt_p]=sex_analysis(he_tt,subj_sex);
[ho_t,ho_p]=sex_analysis(he_ho,subj_sex);
[fs86_t,fs86_p]=sex_analysis(he_fs86,subj_sex);

%correct p-values across all atlases
[p_corr, aal_p_corr, cc200_p_corr, cc400_p_corr, ez_p_corr, fs86_p_corr, ho_p_corr, tt_p_corr]=p_correction(aal_p, cc200_p, cc400_p, ez_p, fs86_p, ho_p, tt_p);

%save files
save(sprintf('hurst_analysis_%d_%d_bhfdr.mat', start_seg, end_seg));
end

