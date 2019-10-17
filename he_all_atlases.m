function he_all_atlases(subject,fmriname)

%read in time series data for each of the HCP subjects 
%subject is a string containing the subject ID 
%fmriname is a string containing the fmri name (i.e. rfMRI_REST1_LR)
V=read_avw_img(sprintf('%s_hp2000_clean',fmriname));
V1=reshape(V,[],size(V,4));

%read in ROI specific data for different atlases
Vparc_wmparc=read_avw_img('wmparc.2');
Vparc1(:,1)=reshape(Vparc_wmparc,[],1);
Vparc_aal=read_avw_img('aal_roi_atlas_2mm_hcp');
Vparc1(:,2)=reshape(Vparc_aal,[],1);
Vparc_cc200=read_avw_img('cc200_roi_atlas_2mm_hcp');
Vparc1(:,3)=reshape(Vparc_cc200,[],1);
Vparc_cc400=read_avw_img('cc400_roi_atlas_2mm_hcp');
Vparc1(:,4)=reshape(Vparc_cc400,[],1);
Vparc_ez=read_avw_img('ez_roi_atlas_2mm_hcp');
Vparc1(:,5)=reshape(Vparc_ez,[],1);
Vparc_ho=read_avw_img('ho_roi_atlas_2mm_hcp');
Vparc1(:,6)=reshape(Vparc_ho,[],1);
Vparc_tt=read_avw_img('tt_roi_atlas_2mm_hcp');
Vparc1(:,7)=reshape(Vparc_tt,[],1);

%remove all voxels that are not present in any of the atlas masks
zero_rows = find(all(Vparc1==0,2));
V_nonzero=V1;
V_nonzero(zero_rows,:)=[];
Vparc1(zero_rows,:)=[];

%get new copy of non-zero voxels to feed into dfa fast
V_nonzero1=V_nonzero';

%get rid of all zero voxels still present in the time series data 
%need to do this specifically for the dfa_fast to work properly on octave
%if you have zero voxels present, it otherwise results in the entire
%H_voxel matrix being entirely nan
zero_cols = find(all(V_nonzero1==0));
V_nonzero1(:,zero_cols)=[];
Vparc2=Vparc1;
Vparc2(zero_cols,:)=[];

%need to preallocate size to make it efficient
H_voxel=nan(size(V_nonzero1,2),8);

%calculate HE on a voxelwise basis
H_voxel(:,1)=dfa_fast(V_nonzero1,11,160,[10,15,25,30]);
H_voxel(:,2)=dfa_fast(V_nonzero1,159,308,[10,15,25,30]);
H_voxel(:,3)=dfa_fast(V_nonzero1,307,456,[10,15,25,30]);
H_voxel(:,4)=dfa_fast(V_nonzero1,455,604,[10,15,25,30]);
H_voxel(:,5)=dfa_fast(V_nonzero1,603,752,[10,15,25,30]);
H_voxel(:,6)=dfa_fast(V_nonzero1,751,900,[10,15,25,30]);
H_voxel(:,7)=dfa_fast(V_nonzero1,899,1048,[10,15,25,30]);
H_voxel(:,8)=dfa_fast(V_nonzero1,1047,1196,[10,15,25,30]);

%save('-mat7-binary',sprintf('%s_%s_voxelwise_he.mat',subject,fmriname),'H_voxel');
save(sprintf('%s_%s_voxelwise_he.mat',subject,fmriname),'H_voxel');

clear Vparc_wmparc Vparc_aal Vparc_cc200 Vparc_cc400 Vparc_ez Vparc_ho Vparc_tt;
clear V_nonzero1;
clear zero_cols zero_rows;

%get all voxels in ROIs from wmparc.2 atlas
%roi=unique(Vparc1(:,1));
%roi=roi(roi~=0);
gmlabels=[8;10;11;12;13;17;18;26;28;47;49;50;51;52;53;54;58;60;1001;
    1002;1003;1005;1006;1007;1008;1009;1010;1011;1012;1013;1014;1015;
    1016;1017;1018;1019;1020;1021;1022;1023;1024;1025;1026;1027;1028;
    1029;1030;1031;1032;1033;1034;1035;2001;2002;2003;2005;2006;2007;
    2008;2009;2010;2011;2012;2013;2014;2015;2016;2017;2018;2019;2020;
    2021;2022;2023;2024;2025;2026;2027;2028;2029;2030;2031;2032;2033;
    2034;2035];
roi=gmlabels;

for i=1:length(roi);
    %get time series for roi
    Vparc_roi=Vparc1(:,1)==roi(i);
    Vparc_roi2=Vparc2(:,1)==roi(i);
    V_roi=(V_nonzero(Vparc_roi,:));
    H_roi{i}=(H_voxel(Vparc_roi2,:));
    
    %get roi size and number of zero columns
    roi_size(i,1) = size(V_roi,1);
    zero_cols = find(all(V_roi==0,2));
    
    %remove all zero columns
    V_roi_nonzero=V_roi;
    V_roi_nonzero(zero_cols,:)=[];
    
    %compute roi size after removing zero columns and # of roi columns
    %roi_size columns 2 + column 3 should equal roi_size column 1
    roi_size(i,2) = size(V_roi_nonzero,1);
    roi_size(i,3) = length(zero_cols);
    roi_ts(i,:)= mean(V_roi_nonzero);
    
    %save indices of zero columns
    zero_col_idx{i}=zero_cols;
    
    %calculate voxelwise HE
    
    H_roi_mean(i,:) = mean(cell2mat(H_roi(i)));

    clear V_roi Vparc_roi V_roi_nonzero;
    clear zero_cols;
    clear gm_labels;
end

%save('-mat7-binary',sprintf('%s_%s_fs86_he.mat',subject,fmriname),'H_roi','H_roi_mean', 'roi_size', 'zero_col_idx');
%save('-mat7-binary',sprintf('%s_%s_fs86_ts.mat',subject,fmriname),'roi_ts', 'roi_size', 'zero_col_idx');
save(sprintf('%s_%s_fs86_he.mat',subject,fmriname),'H_roi','H_roi_mean', 'roi_size', 'zero_col_idx');
save(sprintf('%s_%s_fs86_ts.mat',subject,fmriname),'roi_ts', 'roi_size', 'zero_col_idx');
clear roi roi_ts roi_size zero_col_idx;
clear H_roi H_roi_mean;
clear i; 


%get ROIs from aal atlas
roi=unique(Vparc1(:,2));
roi=roi(roi~=0);

for i=1:length(roi);
    %get time series for roi
    Vparc_roi=Vparc1(:,2)==roi(i);
    Vparc_roi2=Vparc2(:,2)==roi(i);
    V_roi=(V_nonzero(Vparc_roi,:));
    H_roi{i}=(H_voxel(Vparc_roi2,:));
    
    %get roi size and number of zero columns
    roi_size(i,1) = size(V_roi,1);
    zero_cols = find(all(V_roi==0,2));
    
    %remove al zero columns
    V_roi_nonzero=V_roi;
    V_roi_nonzero(zero_cols,:)=[];
    
    %compute roi size after removing zero columns and # of roi columns
    %roi_size columns 2 + column 3 should equal roi_size column 1
    roi_size(i,2) = size(V_roi_nonzero,1);
    roi_size(i,3) = length(zero_cols);
    roi_ts(i,:)= mean(V_roi_nonzero);
    
    %save indices of zero columns
    zero_col_idx{i}=zero_cols;
    
    %calculate roi HE mean
    H_roi_mean(i,:) = mean(cell2mat(H_roi(i)));

    clear V_roi Vparc_roi V_roi_nonzero;
    clear zero_cols;
end

%save('-mat7-binary',sprintf('%s_%s_aal_he.mat',subject,fmriname),'H_roi','H_roi_mean', 'roi_size', 'zero_col_idx');
%save('-mat7-binary',sprintf('%s_%s_aal_ts.mat',subject,fmriname),'roi_ts', 'roi_size', 'zero_col_idx');
save(sprintf('%s_%s_aal_he.mat',subject,fmriname),'H_roi','H_roi_mean', 'roi_size', 'zero_col_idx');
save(sprintf('%s_%s_aal_ts.mat',subject,fmriname),'roi_ts', 'roi_size', 'zero_col_idx');
clear roi roi_ts roi_size zero_col_idx;
clear H_roi H_roi_mean;
clear i; 

%get ROIs from cc200 atlas
roi=unique(Vparc1(:,3));
roi=roi(roi~=0);

for i=1:length(roi);
    %get time series for roi
    Vparc_roi=Vparc1(:,3)==roi(i);
    Vparc_roi2=Vparc2(:,3)==roi(i);
    V_roi=(V_nonzero(Vparc_roi,:));
    H_roi{i}=(H_voxel(Vparc_roi2,:));
    
    %get roi size and number of zero columns
    roi_size(i,1) = size(V_roi,1);
    zero_cols = find(all(V_roi==0,2));
    
    %remove al zero columns
    V_roi_nonzero=V_roi;
    V_roi_nonzero(zero_cols,:)=[];
    
    %compute roi size after removing zero columns and # of roi columns
    %roi_size columns 2 + column 3 should equal roi_size column 1
    roi_size(i,2) = size(V_roi_nonzero,1);
    roi_size(i,3) = length(zero_cols);
    roi_ts(i,:)= mean(V_roi_nonzero);
    
    %save indices of zero columns
    zero_col_idx{i}=zero_cols;
    
    %calculate roi HE mean
    H_roi_mean(i,:) = mean(cell2mat(H_roi(i)));

    clear V_roi Vparc_roi V_roi_nonzero;
    clear zero_cols;
end

%save('-mat7-binary',sprintf('%s_%s_cc200_he.mat',subject,fmriname),'H_roi','H_roi_mean', 'roi_size', 'zero_col_idx');
%save('-mat7-binary',sprintf('%s_%s_cc200_ts.mat',subject,fmriname),'roi_ts', 'roi_size', 'zero_col_idx');
save(sprintf('%s_%s_cc200_he.mat',subject,fmriname),'H_roi','H_roi_mean', 'roi_size', 'zero_col_idx');
save(sprintf('%s_%s_cc200_ts.mat',subject,fmriname),'roi_ts', 'roi_size', 'zero_col_idx');

clear roi roi_ts roi_size zero_col_idx;
clear H_roi H_roi_mean;
clear i; 

%get ROIs fromcc400 atlas
roi=unique(Vparc1(:,4));
roi=roi(roi~=0);

for i=1:length(roi);
    %get time series for roi
    Vparc_roi=Vparc1(:,4)==roi(i);
    Vparc_roi2=Vparc2(:,4)==roi(i);
    V_roi=(V_nonzero(Vparc_roi,:));
    H_roi{i}=(H_voxel(Vparc_roi2,:));
    
    %get roi size and number of zero columns
    roi_size(i,1) = size(V_roi,1);
    zero_cols = find(all(V_roi==0,2));
    
    %remove al zero columns
    V_roi_nonzero=V_roi;
    V_roi_nonzero(zero_cols,:)=[];
    
    %compute roi size after removing zero columns and # of roi columns
    %roi_size columns 2 + column 3 should equal roi_size column 1
    roi_size(i,2) = size(V_roi_nonzero,1);
    roi_size(i,3) = length(zero_cols);
    roi_ts(i,:)= mean(V_roi_nonzero);
    
    %save indices of zero columns
    zero_col_idx{i}=zero_cols;
    
    %calculate roi HE mean
    H_roi_mean(i,:) = mean(cell2mat(H_roi(i)));

    clear V_roi Vparc_roi V_roi_nonzero;
    clear zero_cols;
end

%save('-mat7-binary',sprintf('%s_%s_cc400_he.mat',subject,fmriname),'H_roi','H_roi_mean', 'roi_size', 'zero_col_idx');
%save('-mat7-binary',sprintf('%s_%s_cc400_ts.mat',subject,fmriname),'roi_ts', 'roi_size', 'zero_col_idx');
save(sprintf('%s_%s_cc400_he.mat',subject,fmriname),'H_roi','H_roi_mean', 'roi_size', 'zero_col_idx');
save(sprintf('%s_%s_cc400_ts.mat',subject,fmriname),'roi_ts', 'roi_size', 'zero_col_idx');

clear roi roi_ts roi_size zero_col_idx;
clear H_roi H_roi_mean;
clear i; 

%get ROIs from ez atlas
roi=unique(Vparc1(:,5));
roi=roi(roi~=0);

for i=1:length(roi);
    %get time series for roi
    Vparc_roi=Vparc1(:,5)==roi(i);
    Vparc_roi2=Vparc2(:,5)==roi(i);
    V_roi=(V_nonzero(Vparc_roi,:));
    H_roi{i}=(H_voxel(Vparc_roi2,:));
    
    %get roi size and number of zero columns
    roi_size(i,1) = size(V_roi,1);
    zero_cols = find(all(V_roi==0,2));
    
    %remove al zero columns
    V_roi_nonzero=V_roi;
    V_roi_nonzero(zero_cols,:)=[];
    
    %compute roi size after removing zero columns and # of roi columns
    %roi_size columns 2 + column 3 should equal roi_size column 1
    roi_size(i,2) = size(V_roi_nonzero,1);
    roi_size(i,3) = length(zero_cols);
    roi_ts(i,:)= mean(V_roi_nonzero);
    
    %save indices of zero columns
    zero_col_idx{i}=zero_cols;
    
    %calculate roi HE mean
    H_roi_mean(i,:) = mean(cell2mat(H_roi(i)));

    clear V_roi Vparc_roi V_roi_nonzero;
    clear zero_cols;
end

%save('-mat7-binary',sprintf('%s_%s_ez_he.mat',subject,fmriname),'H_roi','H_roi_mean', 'roi_size', 'zero_col_idx');
%save('-mat7-binary',sprintf('%s_%s_ez_ts.mat',subject,fmriname),'roi_ts', 'roi_size', 'zero_col_idx');
save(sprintf('%s_%s_ez_he.mat',subject,fmriname),'H_roi','H_roi_mean', 'roi_size', 'zero_col_idx');
save(sprintf('%s_%s_ez_ts.mat',subject,fmriname),'roi_ts', 'roi_size', 'zero_col_idx');

clear roi roi_ts roi_size zero_col_idx;
clear H_roi H_roi_mean;
clear i; 

%get ROIs from ho atlas
roi=unique(Vparc1(:,6));
roi=roi(roi~=0);

for i=1:length(roi);
    %get time series for roi
    Vparc_roi=Vparc1(:,6)==roi(i);
    Vparc_roi2=Vparc2(:,6)==roi(i);
    V_roi=(V_nonzero(Vparc_roi,:));
    H_roi{i}=(H_voxel(Vparc_roi2,:));
    
   %get roi size and number of zero columns
    roi_size(i,1) = size(V_roi,1);
    zero_cols = find(all(V_roi==0,2));
    
    %remove al zero columns
    V_roi_nonzero=V_roi;
    V_roi_nonzero(zero_cols,:)=[];
    
    %compute roi size after removing zero columns and # of roi columns
    %roi_size columns 2 + column 3 should equal roi_size column 1
    roi_size(i,2) = size(V_roi_nonzero,1);
    roi_size(i,3) = length(zero_cols);
    roi_ts(i,:)= mean(V_roi_nonzero);
    
    %save indices of zero columns
    zero_col_idx{i}=zero_cols;
    
    %calculate roi HE mean
    H_roi_mean(i,:) = mean(cell2mat(H_roi(i)));

    clear V_roi Vparc_roi V_roi_nonzero;
    clear zero_cols;
end
%save('-mat7-binary',sprintf('%s_%s_ho_he.mat',subject,fmriname),'H_roi','H_roi_mean', 'roi_size', 'zero_col_idx');
%save('-mat7-binary',sprintf('%s_%s_ho_ts.mat',subject,fmriname),'roi_ts', 'roi_size', 'zero_col_idx');
save(sprintf('%s_%s_ho_he.mat',subject,fmriname),'H_roi','H_roi_mean', 'roi_size', 'zero_col_idx');
save(sprintf('%s_%s_ho_ts.mat',subject,fmriname),'roi_ts', 'roi_size', 'zero_col_idx');

clear roi roi_ts roi_size zero_col_idx;
clear H_roi H_roi_mean;
clear i; 

%get ROIs from tt atlas
roi=unique(Vparc1(:,7));
roi=roi(roi~=0);

for i=1:length(roi);
    %get time series for roi
    Vparc_roi=Vparc1(:,7)==roi(i);
    Vparc_roi2=Vparc2(:,7)==roi(i);
    V_roi=(V_nonzero(Vparc_roi,:));
    H_roi{i}=(H_voxel(Vparc_roi2,:));
    
    %get roi size and number of zero columns
    roi_size(i,1) = size(V_roi,1);
    zero_cols = find(all(V_roi==0,2));
    
    %remove al zero columns
    V_roi_nonzero=V_roi;
    V_roi_nonzero(zero_cols,:)=[];
    
    %compute roi size after removing zero columns and # of roi columns
    %roi_size columns 2 + column 3 should equal roi_size column 1
    roi_size(i,2) = size(V_roi_nonzero,1);
    roi_size(i,3) = length(zero_cols);
    roi_ts(i,:)= mean(V_roi_nonzero);
    
    %save indices of zero columns
    zero_col_idx{i}=zero_cols;
    
    %calculate roi HE mean
    H_roi_mean(i,:) = mean(cell2mat(H_roi(i)));

    clear V_roi Vparc_roi V_roi_nonzero;
    clear zero_cols;
end

%save('-mat7-binary',sprintf('%s_%s_tt_he.mat',subject,fmriname),'H_roi','H_roi_mean', 'roi_size', 'zero_col_idx');
%save('-mat7-binary',sprintf('%s_%s_tt_ts.mat',subject,fmriname),'roi_ts', 'roi_size', 'zero_col_idx');
save(sprintf('%s_%s_tt_he.mat',subject,fmriname),'H_roi','H_roi_mean', 'roi_size', 'zero_col_idx');
save(sprintf('%s_%s_tt_ts.mat',subject,fmriname),'roi_ts', 'roi_size', 'zero_col_idx');

clear roi roi_ts roi_size zero_col_idx;
clear H_roi H_roi_mean;
clear i; 
