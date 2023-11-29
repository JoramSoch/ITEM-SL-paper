function modelling_auto_create(work_dir, subj_id, MS_name, GLM_names)
% _
% This function auto-creates a first-level batch for a selected subject.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 10/11/2017, 16:20 (1st version)
%         17/11/2017, 15:05 (MS "test")
%         30/11/2017, 20:10 (MS "sectors")
%         15/03/2018, 10:15 (MS "trials")
%         30/11/2018, 11:45 (MS "item")
%         15/11/2023, 15:38 (for upload)


% load batch file
load modelling_template.mat

% prepare modelling
task_str = 'CircRun';           % name of the task
num_sess = 8;                   % number of sessions

% specify data
for j = 1:num_sess
    
    % get 4D NIfTI filename
    func_imgs = strcat(work_dir,'/','sub-',subj_id,'/func/','rarof_','sub-',subj_id,'_task-',task_str,'_run-0',num2str(j),'_bold.nii');
    matlabbatch{1}.spm.stats.fmri_spec.sess(j).scans{1} = func_imgs;
    
    % get RP TXT filename
    real_para = strcat(work_dir,'/','sub-',subj_id,'/func/','rp_arof_','sub-',subj_id,'_task-',task_str,'_run-0',num2str(j),'_bold.txt');
    matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi_reg{1} = real_para;
    
end;
clear func_imgs real_para
matlabbatch_temp = matlabbatch;

% specify designs
for i = 1:numel(GLM_names)
    
    % edit stats directory
    matlabbatch = matlabbatch_temp;
    GLM_dir = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','glm-',GLM_names{i},'/');
    matlabbatch{1}.spm.stats.fmri_spec.dir{1} = GLM_dir;
    
    % edit onset files
    for j = 1:num_sess
        
        %%% TEST MODEL SPACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(MS_name,'test')
            if strcmp(GLM_names{i},'base-JS')
                mat_mat = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{i},'_run-0',num2str(j),'_matrix.mat');
                matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi_reg{1} = mat_mat;
            else
                ons_mat = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{i},'_run-0',num2str(j),'_onsets.mat');
                RPs_txt = strcat(work_dir,'/','sub-',subj_id,'/func/','rp_arof_','sub-',subj_id,'_task-',task_str,'_run-0',num2str(j),'_bold.txt');
                matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi{1} = ons_mat;
                matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi_reg{1} = RPs_txt;
            end;
        end;
        
        %%% MODEL SPACE "SECTORS" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(MS_name,'sectors')
            if strncmp(GLM_names{i},'sect-',5)
                mat_mat = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{i},'_run-0',num2str(j),'_matrix.mat');
                matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi_reg{1} = mat_mat;
            end;
        end;
        
        %%% MODEL SPACE "TRIALS" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(MS_name,'trials')
            mat_mat = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{i},'_run-0',num2str(j),'_matrix.mat');
            matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi = {''};
            matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi_reg{1} = mat_mat;
        end;
        
        %%% MODEL SPACE "ITEM" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(MS_name,'item')
            ons_mat = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{i},'_run-0',num2str(j),'_onsets.mat');
            RPs_txt = strcat(work_dir,'/','sub-',subj_id,'/func/','rp_arof_','sub-',subj_id,'_task-',task_str,'_run-0',num2str(j),'_bold.txt');
            matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi{1} = ons_mat;
            matlabbatch{1}.spm.stats.fmri_spec.sess(j).multi_reg{1} = RPs_txt;
        end;
        
    end;
    
    % remove model estimation
    if strcmp(MS_name,'sectors')
        if ~strcmp(GLM_names{i},'sect-all') && ~strcmp(GLM_names{i},'sect-none')
            matlabbatch = matlabbatch(1);
        end;
    end;
    
    % save batch file
    filename = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{i},'_design.mat');
    save(filename,'matlabbatch');
    
end;
clear GLM_dir mat_mat ons_mat RPs_txt