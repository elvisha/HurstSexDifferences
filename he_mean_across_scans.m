function [he_mean] = he_mean_across_scans(he_scan1, he_scan2, he_scan3, he_scan4)
%compute mean of HE across scans

%takes HE value from 4 scans and computes overall mean
%each of the he_scan files are nxp where n is number of subjects and p is
%number of ROIs

%he_mean is nxp where n is number of subjects and p is number of ROIs
for i = 1:length(he_scan1)
    temp(1,:) = he_scan1(i,:);
    temp(2,:) = he_scan2(i,:);
    temp(3,:) = he_scan3(i,:);
    temp(4,:) = he_scan4(i,:);
    he_mean(i,:) = nanmean(temp);
    clear temp;
end
