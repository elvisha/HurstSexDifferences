function [he_mean] = he_mean_across_segments(he_data, j,k)
%calculate HE mean across segments within same scan

%he_data is 1xn cell where n is number of subjects
%each cell is pxq where p is number of ROIs and q is number of segments

%i and j specify the start and end segments you want to use to compute the
%mean, i.e. to get the mean for segments 4 through 8, i=4, j=8

%he_mean is nxp where n is number of subjects and p is number of ROIs
for i = 1:length(he_data)
    subj_temp = cell2mat(he_data(i));
    he_mean(i,:) = nanmean(subj_temp(:,j:k),2);
end

