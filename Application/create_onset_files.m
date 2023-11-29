function create_onset_files(work_dir, subj_id, MS_name)
% _
% This function creates onset files for first-level GLMs for fMRI.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 10/11/2017, 15:30 (1st version)
%         17/11/2017, 16:55 (MS "test")
%         30/11/2017, 20:05 (MS "sectors")
%         15/03/2018, 10:10 (MS "trials")
%         30/11/2018, 11:40 (MS "item")
%         15/11/2023, 15:37 (for upload)


%%% LOAD ONSETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% paradigm dimensions
S = 8;                          % number of sessions
T = 100;                        % trials per session
C = 48;                         % categories per trial
L = 4;                          % levels per category

% scanning parameters
settings.n  = 220;              % number of scans [1]
settings.TR = 1.5;              % inter-scan interval [s]
settings.dt = 0.01;             % microtime resolution [ms]

% prepare loading
task_str = 'CircRun';           % name of the task
design   = struct([]);          % design information

% load onsets and durations
for j = 1:S
    filename = strcat(work_dir,'/','sub-',subj_id,'/func/','sub-',subj_id,'_task-',task_str,'_run-0',num2str(j),'_events.mat');
    load(filename);
    orth = num2cell(false(size(names)));
    design(j).names     = names;
    design(j).onsets    = onsets;
    design(j).durations = durations;
    design(j).pmod      = pmod;
    design(j).orth      = orth;
end;

% load realignment parameters
for j = 1:S
    filename = strcat(work_dir,'/','sub-',subj_id,'/func/','rp_arof_','sub-',subj_id,'_task-',task_str,'_run-0',num2str(j),'_bold.txt');
    RP = load(filename);
    design(j).RP = RP;
end;


%%% TEST MODEL SPACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(MS_name,'test')

% set models
GLM_names = {'base', 'base-orth', 'base-mc', 'base-JS', 'fix-stim', 'fix-resp'}';

% save onsets
for j = 1:S

    % base model: assemble onsets
    names     = design(j).names(1);
    onsets    = design(j).onsets(1);
    durations = design(j).durations(1);
    pmod      = design(j).pmod(1);
    orth      = design(j).orth(1);
    RP        = design(j).RP;
    
    % base model: save onsets
    filename = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{1},'_run-0',num2str(j),'_onsets.mat');
    save(filename, 'names', 'onsets', 'durations', 'pmod', 'orth');
    
    % base with orthogonalization
    orth_orig = orth;
    orth = num2cell(true(size(names)));
    filename = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{2},'_run-0',num2str(j),'_onsets.mat');
    save(filename, 'names', 'onsets', 'durations', 'pmod', 'orth');
    
    % base with mean-centering
    orth = orth_orig;
    pmod_orig = pmod;
    for i = 1:C
        pmod.param{i} = pmod.param{i} - mean(pmod.param{i});
    end;
    filename = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{3},'_run-0',num2str(j),'_onsets.mat');
    save(filename, 'names', 'onsets', 'durations', 'pmod', 'orth');
    
    % base with own design matrix
    pmod = pmod_orig;
    [R, L] = create_design_matrix(names, onsets, durations, pmod, [], RP, settings);
    filename = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{4},'_run-0',num2str(j),'_matrix.mat');
    save(filename, 'R');
    
    % fixation stimulus model: assemble onsets
    names     = design(j).names(2:3);
    onsets    = design(j).onsets(2:3);
    durations = design(j).durations(2:3);
    pmod      = design(j).pmod(2:3);
    orth      = design(j).orth(2:3);
    
    % fixation stimulus model: save onsets
    filename = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{5},'_run-0',num2str(j),'_onsets.mat');
    save(filename, 'names', 'onsets', 'durations', 'pmod', 'orth');
    
    % fixation response model: assemble onsets
    names     = design(j).names(4:5);
    onsets    = design(j).onsets(4:5);
    durations = design(j).durations(4:5);
    pmod      = design(j).pmod(4:5);
    orth      = design(j).orth(4:5);
    
    % fixation response model: save onsets
    orth = num2cell(true(size(names)));
    filename = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{6},'_run-0',num2str(j),'_onsets.mat');
    save(filename, 'names', 'onsets', 'durations', 'pmod', 'orth');
    
end;
    
end;


%%% MODEL SPACE "SECTORS" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(MS_name,'sectors')

% set models
GLM_names = cell(C+2,1);
for i = 1:C
    GLM_names{i} = strcat('sect-',MF_int2str0(i,2));
end;
GLM_names{C+1} = 'sect-all';
GLM_names{C+2} = 'sect-none';

% save onsets
for j = 1:S
    
    % assemble names, onsets, durations
    conds     = [1,2,3];
    names     = design(j).names(conds);
    onsets    = design(j).onsets(conds);
    durations = design(j).durations(conds);
    pmod_temp = design(j).pmod(conds);
    orth      = design(j).orth(conds);
    RP        = design(j).RP;
    
    % all sectors model: save onsets
    pmod = pmod_temp;
    [R, L] = create_design_matrix(names, onsets, durations, pmod, [], RP, settings);
    filename = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{C+1},'_run-0',num2str(j),'_matrix.mat');
    save(filename, 'R');
    
    % no sectors model: save onsets
    pmod = pmod_temp;
    pmod(1).poly  = [];
    pmod(1).name  = [];
    pmod(1).param = []; 
    [R, L] = create_design_matrix(names, onsets, durations, pmod, [], RP, settings);
    filename = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{C+2},'_run-0',num2str(j),'_matrix.mat');
    save(filename, 'R');
    
    % sector models: save onsets
    for i = 1:C
        pmod = pmod_temp;
        pmod(1).poly  = pmod(1).poly(i);
        pmod(1).name  = pmod(1).name(i);
        pmod(1).param = pmod(1).param(i);
        [R, L] = create_design_matrix(names, onsets, durations, pmod, [], RP, settings);
        filename = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{i},'_run-0',num2str(j),'_matrix.mat');
        save(filename, 'R');
    end;
    
end;

end;


%%% MODEL SPACE "TRIALS" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(MS_name,'trials')

% set models
GLM_names = {'sects', 'trls'}';

% save onsets
for j = 1:S
    
    % assemble names, onsets, durations
    conds     = [1,2,3];
    names     = design(j).names(conds);
    onsets    = design(j).onsets(conds);
    durations = design(j).durations(conds);
    pmod_temp = design(j).pmod(conds);
    orth      = design(j).orth(conds);
    RP        = design(j).RP;
    
    % sector model: create design matrix
    pmod = pmod_temp;
    [R, L] = create_design_matrix(names, onsets, durations, pmod, [], RP, settings);
    filename = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{1},'_run-0',num2str(j),'_matrix.mat');
    save(filename, 'R');
    
    % expand names, onsets, durations
    onsets1    = onsets{1};
    durations1 = durations{1};
    names(T+[1,2])     = names([2,3]);
    onsets(T+[1,2])    = onsets([2,3]);
    durations(T+[1,2]) = durations([2,3]);
    for k = 1:T
        names{k}     = strcat('trl-',MF_int2str0(k,3));
        onsets{k}    = onsets1(k);
        durations{k} = durations1(k);
    end;
    
    % trials model: create design matrix
    pmod = [];
    [R, L] = create_design_matrix(names, onsets, durations, pmod, [], RP, settings);
    filename = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{2},'_run-0',num2str(j),'_matrix.mat');
    save(filename, 'R');
    
end;

end;


%%% MODEL SPACE "ITEM" %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(MS_name,'item')

% set models
GLM_names = {'full'};

% save onsets
for j = 1:S
    
    % full model: assemble onsets
    conds     = [1,2,3];
    names     = design(j).names(conds);
    onsets    = design(j).onsets(conds);
    durations = design(j).durations(conds);
    pmod      = design(j).pmod(conds);
    orth      = design(j).orth(conds);
    
    % full model: save onsets
    filename = strcat(work_dir,'/','sub-',subj_id,'/mods/','glms-',MS_name,'/','sub-',subj_id,'_glm-',GLM_names{1},'_run-0',num2str(j),'_onsets.mat');
    save(filename, 'names', 'onsets', 'durations', 'pmod', 'orth');
    
end;

end;