function [p_corr, p_corr1, p_corr2, p_corr3, p_corr4, p_corr5, p_corr6, p_corr7] = p_correction(p_1,p_2,p_3,p_4,p_5,p_6,p_7)
%corrected p values across all atlases used

%generates corrected p values across all atlases used using benjamini
%hochberg procedure

%p_1, p_2, ... are p values generated from each of the atlases; all are 1xn
%where n is the number of ROIs in that given atlas

%p_corr is array of all corrected values (across all atlases); is 1xp where
%p is the sum of the number of ROIs across all atlases

%p_corr1, p_corr2, ... are corrected p values for each of the atlases; all
%are 1xn where n is the number of ROIs in that given atlas

p_all = [p_1, p_2, p_3, p_4, p_5, p_6, p_7];

%do this both ways with BHFDR true and false
%BHFDR true uses Benjamini-Hochberg (1995) method - more stringent; use
%this to report all results
%BHFDR false (default) uses Storey (2002) method - less stringent
p_corr = mafdr(p_all,'BHFDR',true);

temp = p_corr;

p_corr1 = temp(1:length(p_1));
temp(1:length(p_1)) = [];
p_corr2 = temp(1:length(p_2));
temp(1:length(p_2)) = [];
p_corr3 = temp(1:length(p_3));
temp(1:length(p_3)) = [];
p_corr4 = temp(1:length(p_4));
temp(1:length(p_4)) = [];
p_corr5 = temp(1:length(p_5));
temp(1:length(p_5)) = [];
p_corr6 = temp(1:length(p_6));
temp(1:length(p_6)) = [];
p_corr7 = temp(1:length(p_7));
temp(1:length(p_7)) = [];

end

