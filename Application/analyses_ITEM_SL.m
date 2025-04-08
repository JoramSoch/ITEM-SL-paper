% Searchlight-based ITEM analysis: pipeline script
% _
% This script performs the whole data analysis underlying the results
% reported in Soch (2025), revised for Imaging Neuroscience.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% Version History:
% - 14/11/2017: Step 2
% - 19/12/2018: Step 1,3,4
% - 24/11/2021: Step 5
% - 07/12/2021: Step 7
% - 09/12/2021: Step 8
% - 28/02/2023: Step 0
% - 01/03/2023: Step 6
% - 31/10/2023, 14:41: prepared for upload
% - 15/11/2023, 15:53: prepared for upload
% - 21/11/2023, 16:02: finalized for upload
% - 18/12/2024, 19:01: added LS-A, LS-S, GLMsingle
% - 19/12/2024, 09:08: added group analyses


clear
close all

%%% Step 0: download dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify directories
project_directories

% download dataset
download_dataset

% launch SPM
spm fmri
spm_jobman('initcfg');
load project_directories.mat


%%% Step 1: specify analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% specify subjects
sub2do = [1:4]; % AAA03, AAA02, AAX01, AAX02

% specify SL radii
SL_rad = 6;     % in millimeters

% specify methods
methods ={'ITEM', 'LS-A', 'LS-S', 'GLM-single'};

% select subjects
if numel(sub2do) == 1, subj_file = 'subjects_1st.mat'; end;
if numel(sub2do) == 3, subj_file = 'subjects_2-4.mat'; end;
if numel(sub2do) == 4, subj_file = 'subjects_all.mat'; end;

% select models
MS_file = 'model_spaces/glms-item.mat';

% load subjects and models
load(subj_file);
load(MS_file);
num_subj = numel(subj_ids);     % number of subjects
num_mods = numel(GLM_names);    % number of models

% set other parameters
num_sess = 8;                   % number of sessions
num_scan = 220;                 % number of scans per session
num_trls = 100;                 % number of trials per session
num_sect = 48;                  % number of sectors per trial


%%% Step 2: fMRI preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preprocess subjects
preproc_meta_batch(subj_file, 'prepare')
preproc_meta_batch(subj_file, 'perform')


%%% Step 3: estimate standard GLMs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% analyze subjects
modelling_meta_batch(subj_file, MS_file, 'prepare')
modelling_meta_batch(subj_file, MS_file, 'perform')


%%% Step 4: estimate trial-wise GLMs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display message
fprintf('\n\n-> Estimate trial-wise GLMs:\n');

% for all subjects
for i = 1:num_subj
    % load SPM.mat
    fprintf('   - Subject %s (%d out of %d):\n', subj_ids{i}, i, num_subj);
    SPM_mat = strcat(bids_dir,'/','sub-',subj_ids{i},'/mods/','glms-',MS_name,'/','glm-',GLM_names{1},'/','SPM.mat');
    load(SPM_mat);
    % for all methods
    for j = 1:numel(methods)
        if strcmp(methods{j},'ITEM')                    % ITEM
            ITEM_est_1st_lvl(SPM, 'DCT', [2 3]);
        elseif strcmp(methods{j},'LS-A')                % LS-A
            ITEM_est_1st_lvl_LS_A(SPM, 'DCT', [2 3]);   
        elseif strcmp(methods{j},'LS-S')                % LS-S
            ITEM_est_1st_lvl_LS_S(SPM, 'DCT', [2 3]);   
        elseif strcmp(methods{j},'GLM-single')          % GLMsingle
            GLMsingle_est_1st_lvl(SPM, [2 3]);          
        end;
    end;
end;


%%% Step 5: perform SL-based ITEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display message
fprintf('\n-> Perform searchlight-based ITEM:\n');
    
% for all subjects
for i = 1:num_subj
    % load SPM.mat
    fprintf('   - Subject %s (%d out of %d):\n', subj_ids{i}, i, num_subj);
    SPM_mat = strcat(bids_dir,'/','sub-',subj_ids{i},'/mods/','glms-',MS_name,'/','glm-',GLM_names{1},'/','SPM.mat');
    load(SPM_mat);
    % for all methods
    for j = 1:numel(methods)
        fprintf('     - %s, searchlight radius r = %d mm:\n', methods{j}, SL_rad);
        % specify analysis
        rad  = SL_rad;
        c    = [0, ones(1,num_sect)];
        con  = 'sects-all';
        meth = methods{j};
        % perform decoding
        if strcmp(methods{j},'ITEM')                    % ITEM
            ITEM_dec_recon_SL(SPM, rad, c, con);
        elseif strncmp(methods{j},'LS',2)               % LS-A or LS-S
            ITEM_dec_recon_SL_SVR(SPM, rad, c, con, meth);
        elseif strcmp(methods{j},'GLM-single')          % GLMsingle
            GLMsingle_dec_recon_SL_SVR(SPM, rad, c, con, meth)
        end;
    end;
end;


%%% Step 6: delete unnecessary files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display message
fprintf('\n-> Delete unnecessary files:\n');  
  
% for all subjects
for i = 1:num_subj
    % diplay message
    fprintf('   - Subject %s (%d out of %d):\n', subj_ids{i}, i, num_subj);
    % for all methods
    for j = 1:1 % ITEM only
        % specify radius/contrast
        fprintf('     - %s, searchlight radius r = %d mm:\n', methods{j}, SL_rad);
        rad = SL_rad;
        con = 'sects-all';
        % delete images
        fprintf('       - delete images ... ');
        GLM_dir = strcat(bids_dir,'/','sub-',subj_ids{i},'/mods/','glms-',MS_name,'/','glm-',GLM_names{1},'/');
        ana_dir = strcat(GLM_dir,'ITEM_dec_recon/','ITEM_',con,'_SL-',num2str(rad),'mm/');
        files   = [dir(strcat(ana_dir,'avgCC_*.nii')); dir(strcat(ana_dir,'oosCC_*.nii'))];
        for k = 1:numel(files)
            delete(strcat(files(k).folder,'/',files(k).name));
        end;
        fprintf('done.\n');
        % delete decodings
        fprintf('       - delete decodings ... ');
        ITEM_mat  = strcat(ana_dir,'/','ITEM.mat');
        load(ITEM_mat);
        ITEM      = rmfield(ITEM, 'VoosCC');
        ITEM      = rmfield(ITEM, 'VavgCC');
        ITEM.Sess = rmfield(ITEM.Sess,'Yp');
        save(ITEM_mat, 'ITEM');
        fprintf('done.\n');
    end;
end;


%%% Step 7: normalize ITEM-SL maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display message
fprintf('\n-> Normalize ITEM-SL maps:\n');  

% for all subjects
for i = 1:num_subj
    % diplay message
    fprintf('   - Subject %s (%d out of %d):\n', subj_ids{i}, i, num_subj);
    % for all methods
    for j = 1:numel(methods)
        % diplay message
        fprintf('     - %s, searchlight radius r = %d mm:\n', methods{j}, SL_rad);
        normalize_subj_rad_img(subj_ids{i}, methods{j}, SL_rad, 'sects-all', 'cvCC', [3 3 3]);
    end;
end;


%%% Step 8: perform group-level analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display message
fprintf('\n->Perform group-level analysis:\n');  

% for all methods
for j = 1:numel(methods)
    
    % diplay message
    meth_str = methods{j};
    fprintf('   - Analysis "%s" (%d out of %d):\n', meth_str, j, numel(methods));
    
    % run group-level analyses
    SL_rad   = 6;                   % searchlight radius: 6 mm
    ana_name = 'sects-all';         % analysis name: "all sectors"
    img_name = 'cvCC';              % cvCC = cross-validated correlation coefficient
    create_full_factorial(meth_str, SL_rad, ana_name, strcat('w',img_name));
    
    % save thresholded SPMs
    SPM_mat = strcat(stat_dir,'glms-',MS_name,'_','glm-',GLM_names{1},'_',...
                     meth_str,'_',ana_name,'_SL-',num2str(rad),'mm','_w',img_name,'/','SPM.mat');
    load(SPM_mat);
    con = 53+[1:26];
    for i = 1:numel(con)
        % diplay message
        con_name = SPM.xCon(con(i)).name;
        fprintf('     - Contrast "%s" (%d out of %d):\n', con_name, i, numel(con));
        % threshold SPM
        if contains(con_name,'<>')
            con_name = strcat(con_name(1:strfind(con_name,'<>')-1),'_neq_',con_name(strfind(con_name,'<>')+2:end));
        end;
        if contains(con_name,'>')
            con_name = strcat(con_name(1:strfind(con_name,'>')-1),'_gr_',con_name(strfind(con_name,'>')+1:end));
        end;
        filename = strcat('con_',MF_int2str0(con(i),4),'_',con_name,'_FWE_0.05_0.nii');
        spm_save_thr_SPM(SPM, con(i), true, 0.05, 0, filename);
        filename = strcat('con_',MF_int2str0(con(i),4),'_',con_name,'_unc_0.001_10.nii');
        spm_save_thr_SPM(SPM, con(i), false, 0.001, 10, filename);
    end;
    
end;