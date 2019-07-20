function generate_brain_imgs(atlas, atlas_roi, atlas_mask, atlas_values, filename)
%generate nifti images of the brain with values projected onto the ROIs
%i.e. t stat or feature importances

%atlas is name of specific atlas you want to plot onto
%atlas_roi is double containing specific rois you're looking to plot onto - 
%is 1xn where n is thenumber of ROIs
%atlas_mask is logical mask that you want to apply - is 1xn
%if all values are to be plotted, atlas_mask can be entirely 1
%atlas_values are specific values you want to plot - is 1xn

%filename is specific name you want to call the nifti

%function does not produce any outputs - values are plotted onto desired
%atlas and nifti file is directly saved by calling save_avw//5re67u

%for some reason read_avw works here but read_avw_img does not produce
%images as intended when saving  later on
Vparc=read_avw(sprintf('%s',atlas));
Vparc1=reshape(Vparc,[],1);

for i=1:length(Vparc1)
    if ismember(Vparc1(i),atlas_roi)
        idx=find(atlas_roi==Vparc1(i));
        if atlas_mask(idx)==1
            Vparc2(i)=atlas_values(idx);
            
        else
            Vparc2(i)=NaN;
        end
    else
        Vparc2(i)=NaN;
    end
end

Vparc3=reshape(Vparc2,91,109,91);

save_avw(Vparc3,sprintf('%s',filename),'f',[2.0 2.0 2.0 1])
end

