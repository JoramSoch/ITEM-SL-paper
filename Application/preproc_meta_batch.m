function preproc_meta_batch(subj_file, proc_mode)
% _
% This function executes preprocessing batches for selected subjects.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 14/11/2017, 15:10


tic
fprintf('\n');

% load directories and subjects
load project_directories.mat
load(subj_file);
num_subj = numel(subj_ids);

% load reorientation parameters
[num, txt, raw] = xlsread('subjects_ro.xls');
ro_paras = NaN(num_subj,12);
for i = 1:num_subj
    for j = 2:size(raw,1)
        if strcmp(raw{j,2},subj_ids{i})
            ro_paras(i,:) = cell2mat(raw(j,3:end));
        end;
    end;
end;    

% create batches, if mode is "prepare"
if strcmp(proc_mode,'prepare')
        
    % create preproc batches
    for i = 1:num_subj
        fprintf('-> Creating preprocessing batch for %s (%d out of %d) ... ', subj_ids{i}, i, num_subj);
        preproc_auto_create(bids_dir, subj_ids{i}, ro_paras(i,:));
        fprintf('successful!\n');
    end;
    fprintf('\n');
    
end;

% execute batches, if mode is "perform"
if strcmp(proc_mode,'perform')
    
    % execute preproc batches
    num_errs = 0;
    for i = 1:num_subj
        fprintf('-> Executing preprocessing batch for %s (%d out of %d) ... \n', subj_ids{i}, i, num_subj);
        filename = strcat(bids_dir,'/','sub-',subj_ids{i},'/','preproc_',subj_ids{i},'.mat');
        try 
           spm_jobman('run', filename);
        catch
           num_errs = num_errs + 1;
           err_subj{num_errs} = subj_ids{i};
        end
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