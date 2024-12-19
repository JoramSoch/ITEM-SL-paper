function normalize_subj_rad_img(subj_id, meth_str, SL_rad, ana_name, img_name, vox_size)
% _
% Create and perform normalization for given analysis
%     subj_id  - subject ID (e.g. 'AAA-03')
%     meth_str - method string (i.e. 'ITEM', 'LS-A', 'LS-S' or 'GLM-single')
%     SL_rad   - searchlight radius (in mm)
%     ana_name - analysis name (e.g. 'sects-all')
%     img_name - image filename (e.g. 'cvCC')
%     vox_size - normalized voxel size (1 x 3 vector)
% 
% written by Joram Soch <joram.soch@bccn-berlin.de>, 07/12/2021, 17:08;
% finalized: 15/11/2023, 16:19; updated: 19/12/2024, 09:07


% specify sectors
sect_inds = [1:48];
num_sect  = numel(sect_inds);
MS_name   = 'item';
GLM_names ={'full'};
meth_str(strfind(meth_str,'-')) = '_';

% load directories and files
load project_directories.mat
GLM_dir   = strcat(bids_dir,'sub-',subj_id,'/mods/','glms-',MS_name,'/','glm-',GLM_names{1},'/');
anat_img  = strcat(bids_dir,'sub-',subj_id,'/anat/','y_rof_sub-',subj_id,'_T1w.nii');
func_imgs = cell(num_sect,1);
for j = 1:num_sect
    img_num = sect_inds(j)+1;
    img_str = MF_int2str0(img_num,4);
    if strcmp(meth_str,'ITEM')
        func_imgs{j} = strcat(GLM_dir,'ITEM_dec_recon/ITEM_',...
                              ana_name,'_SL-',num2str(SL_rad),'mm/',...
                              img_name,'_',img_str,'.nii');
    elseif strncmp(meth_str,'LS',2)
        func_imgs{j} = strcat(GLM_dir,'ITEM_dec_recon_',meth_str,'/ITEM_',...
                              ana_name,'_SL-',num2str(SL_rad),'mm/',...
                              img_name,'_',img_str,'.nii');
    elseif strcmp(meth_str,'GLM_single')
        func_imgs{j} = strcat(GLM_dir,meth_str,'_SVR','/SVR_',...
                              ana_name,'_SL-',num2str(SL_rad),'mm/',...
                              img_name,'_',img_str,'.nii');
    end;
end;

% configure segmentation
matlabbatch{1}.spm.spatial.preproc.channel.vols     = {anat_img};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write    = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm    = {'C:\spm\spm12_r7771\tpm\TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus  = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm    = {'C:\spm\spm12_r7771\tpm\TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus  = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm    = {'C:\spm\spm12_r7771\tpm\TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus  = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm    = {'C:\spm\spm12_r7771\tpm\TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus  = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm    = {'C:\spm\spm12_r7771\tpm\TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus  = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm    = {'C:\spm\spm12_r7771\tpm\TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus  = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf         = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup     = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg         = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg      = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm        = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp        = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write       = [0 1];
matlabbatch{1}.spm.spatial.preproc.warp.vox         = NaN;
matlabbatch{1}.spm.spatial.preproc.warp.bb          = [NaN NaN NaN
                                                       NaN NaN NaN];

% configure normalization
matlabbatch{2}.spm.spatial.normalise.write.subj.def(1)     = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{2}.spm.spatial.normalise.write.subj.resample   = func_imgs;
matlabbatch{2}.spm.spatial.normalise.write.woptions.bb     = [-78 -112 -70
                                                               78   76  85];
matlabbatch{2}.spm.spatial.normalise.write.woptions.vox    = vox_size;
matlabbatch{2}.spm.spatial.normalise.write.woptions.interp = 4;
matlabbatch{2}.spm.spatial.normalise.write.woptions.prefix = 'w';

% save and execute batch
filename = strcat(bids_dir,'sub-',subj_id,'/','normalize','_',subj_id,'_',num2str(SL_rad),'mm','_',img_name,'.mat');
save(filename, 'matlabbatch');
spm_jobman('run', filename);