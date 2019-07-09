function [t_stat,p] = sex_analysis(data,subj_sex)
%to study whether sex differences exist in data

%data is nxp where n is number of subjects and p is number of variables
%(i.e. ROIs, time points, time segments, etc)
%subj_sex is nx1 double where 1==male and 0==female

%t_stat is t_stats generated for the comparison of male vs female in data

%p is the uncorrected p values generated for the comparison
%p values need to be corrected for comparisons across all atlases (as
%opposed to at a single atlas level) - this is done using the p_correction
%function

male=subj_sex==1;
male_data = cell(size(male));
male_data = data(male,:);

female=subj_sex==0;
female_data= cell(size(female));
female_data = data(female,:);

[h,p,ci,stats]=ttest2(male_data, female_data, 'vartype', 'unequal');
t_stat=stats.tstat;
%p_corr=mafdr(p,'BHFDR',true);

end

