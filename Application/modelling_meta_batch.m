function modelling_meta_batch(subj_file, MS_file, proc_mode)
% _
% This function executes first-level batches for selected subjects.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 14/11/2017, 16:05; 15/11/2023, 15:36


tic
fprintf('\n');

% load directories and subjects
load project_directories.mat
load(subj_file);
load(MS_file);
num_subj = numel(subj_ids);
num_mods = numel(GLM_names);

% create batches, if mode is "prepare"
if strcmp(proc_mode,'prepare')
    
    % create model directories
    for i = 1:num_subj
        fprintf('-> Creating model directories for %s (%d out of %d) ... ', subj_ids{i}, i, num_subj);
        % analysis folder
        mods_dir = strcat(bids_dir,'/','sub-',subj_ids{i},'/mods/');
        if ~exist(mods_dir,'dir'), mkdir(mods_dir); end;
        % model space
        MS_dir = strcat(bids_dir,'/','sub-',subj_ids{i},'/mods/','glms-',MS_name,'/');
        if ~exist(MS_dir,'dir'), mkdir(MS_dir); end;
        % single models
        for j = 1:num_mods
            GLM_dir = strcat(bids_dir,'/','sub-',subj_ids{i},'/mods/','glms-',MS_name,'/','glm-',GLM_names{j},'/');
            if ~exist(GLM_dir,'dir'), mkdir(GLM_dir); end;
        end;
        fprintf('successful!\n');
    end;
    fprintf('\n');
    
    % create onsets and durations
    for i = 1:num_subj
        fprintf('-> Creating onsets and durations for %s (%d out of %d) ... ', subj_ids{i}, i, num_subj);
        create_onset_files(bids_dir, subj_ids{i}, MS_name);
        fprintf('successful!\n');
    end;
    fprintf('\n');
    
    % create modelling batches
    for i = 1:num_subj
        fprintf('-> Creating modelling batches for %s (%d out of %d) ... ', subj_ids{i}, i, num_subj);
        modelling_auto_create(bids_dir, subj_ids{i}, MS_name, GLM_names);
        fprintf('successful!\n');
    end;
    fprintf('\n');
    
end;

% execute batches, if mode is "perform"
if strcmp(proc_mode,'perform')
    
    % execute modelling batches
    num_errs = 0;
    for i = 1:num_subj
        for j = 1:num_mods
            fprintf('-> Executing modelling batch "%s" for "%s" (%d out of %d in %d out of %d) ... \n', GLM_names{j}, subj_ids{i}, j, num_mods, i, num_subj);
            filename = strcat(bids_dir,'/','sub-',subj_ids{i},'/mods/','glms-',MS_name,'/','sub-',subj_ids{i},'_glm-',GLM_names{j},'_design.mat');
            try
               spm_jobman('run', filename);
            catch
               num_errs = num_errs + 1;
               err_subj{num_errs,1} = subj_nos{i};
               err_subj{num_errs,2} = MS_mods{j};
            end
        end;
    end;
    fprintf('\n');
    
    % display error subjects
    if num_errs > 0
        fprintf('-> The following subjects did not run properly:\n');
        disp(err_subj);
        fprintf('\n');
    end;
    
    % change to tools directory
    cd(tool_dir);
    
end;

toc
fprintf('\n');