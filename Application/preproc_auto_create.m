function preproc_auto_create(work_dir, subj_id, ro_para)
% _
% This function auto-creates a preprocessing batch for a selected subject.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 14/11/2017, 15:15


% load batch file
load preproc_template.mat

% set session folders
sess_ids = {'CircRun', 'CircRun', 'CircRun', 'CircRun', 'CircRun', 'CircRun', 'CircRun', 'CircRun', 'CircCon'};
num_sess = numel(sess_ids);
num_vols = zeros(size(sess_ids));

% specify reorient images
images{1} = strcat(work_dir,'/','sub-',subj_id,'/anat/','sub-',subj_id,'_T1w.nii');
num_imgs  = 1;
for j = 1:num_sess
    % get 4D NIfTI filename
    if j < num_sess
        func_imgs = strcat(work_dir,'/','sub-',subj_id,'/func/','sub-',subj_id,'_task-',sess_ids{j},'_run-0',num2str(j),'_bold.nii');
    else
        func_imgs = strcat(work_dir,'/','sub-',subj_id,'/func/','sub-',subj_id,'_task-',sess_ids{j},'_bold.nii');
    end;
    % get number of volumes
    epi_hdr = spm_vol(func_imgs);
    num_vol = numel(epi_hdr);
    % edit volumes in batch
    for k = 1:num_vol
        images{num_imgs+k,1} = strcat(func_imgs,',',num2str(k));
    end;
    num_imgs = num_imgs + num_vol;
    num_vols(j) = num_vol;
end;
matlabbatch{1}.spm.util.reorient.srcfiles = images;
matlabbatch{1}.spm.util.reorient.transform.transprm = ro_para;
ro_pref = matlabbatch{1}.spm.util.reorient.prefix;
clear func_imgs epi_hdr num_vol images

% specify slice timing
sessions = cell(num_sess,1);
for j = 1:num_sess
    % get 4D NIfTI filename
    if j < num_sess
        func_imgs = strcat(work_dir,'/','sub-',subj_id,'/func/',ro_pref,'sub-',subj_id,'_task-',sess_ids{j},'_run-0',num2str(j),'_bold.nii');
    else
        func_imgs = strcat(work_dir,'/','sub-',subj_id,'/func/',ro_pref,'sub-',subj_id,'_task-',sess_ids{j},'_bold.nii');
    end;
    % edit volumes in batch
    scans = cell(num_vols(j),1);
    for k = 1:num_vols(j)
        scans{k} = strcat(func_imgs,',',num2str(k));
    end;
    sessions{j} = scans;
end;
matlabbatch{2}.spm.temporal.st.scans = sessions;
clear func_imgs sessions scans

% specify spatial realignment
for j = 1:num_sess
    % edit dependency in batch
    matlabbatch{3}.spm.spatial.realign.estwrite.data{j}(1) = cfg_dep(sprintf('Slice Timing: Slice Timing Corr. Images (Sess %d)',j), substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{j}, '.','files'));
end;

% save batch file
filename = strcat(work_dir,'/','sub-',subj_id,'/','preproc_',subj_id,'.mat');
save(filename, 'matlabbatch');