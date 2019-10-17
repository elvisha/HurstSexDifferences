function get_volume(subject)
%setenv('FSLDIR','/usr/local/fsl')

V_aal=read_avw(sprintf('%s_aal',subject));
V1_aal=reshape(V_aal,[],size(V_aal,4));

V_cc200=read_avw(sprintf('%s_cc200',subject));
V1_cc200=reshape(V_cc200,[],size(V_cc200,4));

V_cc400=read_avw(sprintf('%s_cc400',subject));
V1_cc400=reshape(V_cc400,[],size(V_cc400,4));

V_ez=read_avw(sprintf('%s_ez',subject));
V1_ez=reshape(V_ez,[],size(V_ez,4));

V_ho=read_avw(sprintf('%s_ho',subject));
V1_ho=reshape(V_ho,[],size(V_ho,4));

V_tt=read_avw(sprintf('%s_tt',subject));
V1_tt=reshape(V_tt,[],size(V_tt,4));

% Vparc_aal=read_avw('aal_roi_atlas_2mm_hcp');
% Vparc1(:,1)=reshape(Vparc_aal,[],1);
% Vparc_cc200=read_avw('cc200_roi_atlas_2mm_hcp');
% Vparc1(:,2)=reshape(Vparc_cc200,[],1);
% Vparc_cc400=read_avw('cc400_roi_atlas_2mm_hcp');
% Vparc1(:,3)=reshape(Vparc_cc400,[],1);
% Vparc_ez=read_avw('ez_roi_atlas_2mm_hcp');
% Vparc1(:,4)=reshape(Vparc_ez,[],1);
% Vparc_ho=read_avw('ho_roi_atlas_2mm_hcp');
% Vparc1(:,5)=reshape(Vparc_ho,[],1);
% Vparc_tt=read_avw('tt_roi_atlas_2mm_hcp');
% Vparc1(:,6)=reshape(Vparc_tt,[],1);

% get roi size for aal
roi=unique(V1_aal);
roi=roi(roi~=0);
for i=1:length(roi)
    Vparc_roi=V1_aal==roi(i);
    V_roi=(V1_aal(Vparc_roi,:));
    roi_size_aal(i,1) = size(V_roi,1);
end
clear roi V_roi Vparc_roi;

% get roi size for cc200
roi=unique(V1_cc200);
roi=roi(roi~=0);
for i=1:length(roi)
    Vparc_roi=V1_cc200==roi(i);
    V_roi=(V1_cc200(Vparc_roi,:));
    roi_size_cc200(i,1) = size(V_roi,1);
end
clear roi V_roi Vparc_roi;

% get roi size for cc400
roi=unique(V1_cc400);
roi=roi(roi~=0);
for i=1:length(roi)
    Vparc_roi=V1_cc400==roi(i);
    V_roi=(V1_cc400(Vparc_roi,:));
    roi_size_cc400(i,1) = size(V_roi,1);
end
clear roi V_roi Vparc_roi;

% get roi size for ez
roi=unique(V1_ez);
roi=roi(roi~=0);
for i=1:length(roi)
    Vparc_roi=V1_ez==roi(i);
    V_roi=(V1_ez(Vparc_roi,:));
    roi_size_ez(i,1) = size(V_roi,1);
end
clear roi V_roi Vparc_roi;

% get roi size for ho
roi=unique(V1_ho);
roi=roi(roi~=0);
for i=1:length(roi)
    Vparc_roi=V1_ho==roi(i);
    V_roi=(V1_ho(Vparc_roi,:));
    roi_size_ho(i,1) = size(V_roi,1);
end
clear roi V_roi Vparc_roi;

% get roi size for tt
roi=unique(V1_tt);
roi=roi(roi~=0);
for i=1:length(roi)
    Vparc_roi=V1_tt==roi(i);
    V_roi=(V1_tt(Vparc_roi,:));
    roi_size_tt(i,1) = size(V_roi,1);
end
clear roi V_roi Vparc_roi;

save(sprintf('%s_atlas_vol.mat',subject),'roi_size_aal','roi_size_cc200','roi_size_cc400','roi_size_ez','roi_size_ho','roi_size_tt');


