function fact_des_full_fact(subj_file, MS_file, img_paths, var_names, con_info, stat_suff, run_batch)
% _
% Create batch for factorial design full factorial analysis
%     subj_file - a MAT file containing subject IDs and subject info
%     MS_file   - a MAT file containg model space name and GLM names
%     img_paths - image path names specifying within-subject factors
%     var_names - variable indices specifying between-subject factors
%     con_info  - a struct variable specifying contrast vectors
%     stat_suff - a string appended as suffix to analysis directory
%     run_batch - a logical indicating whether the batch is executed
% 
% written by Joram Soch <joram.soch@bccn-berlin.de>, 07/12/2021, 10:37;
% adapted: 09/12/2021, 13:23; finalized: 21/11/2023, 15:39


% load directories and subjects
load project_directories.mat
load(subj_file);
load(MS_file);
subj_ids = subj_ids';
num_subj = numel(subj_ids);
num_mods = numel(GLM_names);

% get factors and cells
num_wsf  = sum(size(img_paths)>1);
num_bsf  = numel(var_names);
num_fact = num_wsf + num_bsf;
for i = 1:num_wsf
    facts(i).name = strcat('cond',num2str(i));
    if num_wsf == 1, facts(i).name = 'cond'; end;
    facts(i).lvls = size(img_paths,i);
    facts(i).dep  = 1; % dependence yes
    facts(i).var  = 0; % variance equal
end;
for i = (num_wsf+1):num_fact
    facts(i).name = var_names{i-num_wsf};
    facts(i).name(strfind(facts(i).name,'_')) = '-';
    facts(i).vals = cell2mat(subj_info(:,var_names(i)));
    facts(i).lvls = max(facts(i).vals);
    facts(i).dep  = 0; % dependence no
    facts(i).var  = 1; % variance unequal
end;
num_lvls = horzcat(facts.lvls);
num_cell = prod(num_lvls);
cells = zeros(num_cell,num_fact);
for i = 1:num_fact
    if i == 1
        cells(:,i) = kron([1:num_lvls(i)]',ones(num_cell/num_lvls(i),1));
    elseif i == num_fact
        cells(:,i) = kron(ones(num_cell/num_lvls(i),1),[1:num_lvls(i)]');
    else
        cells(:,i) = repmat(kron([1:num_lvls(i)]',ones(prod(num_lvls(i+1:end)),1)), [prod(num_lvls(1:i-1)) 1]);
    end;
end;

% specify directory
ana_dir = strcat(stat_dir,'/','glms-',MS_name,'_','glm-',GLM_names{1},'_',stat_suff,'/');
matlabbatch{1}.spm.stats.factorial_design.dir = {ana_dir};

% specify factors
for i = 1:num_fact
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(i).name     = facts(i).name;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(i).levels   = facts(i).lvls;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(i).dept     = facts(i).dep;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(i).variance = facts(i).var;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(i).gmsca    = 0; % no GM scaling
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(i).ancova   = 0; % no ANCOVA
end;

% specify cells
for j = 1:num_cell
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(j).levels = cells(j,:)';
    % check between-subject factors
    if num_bsf == 0
        subj_inds = true(num_subj,1);
    else
        subj_vars = horzcat(facts(num_wsf+1:end).vals);
        req_vars  = repmat(cells(j,num_wsf+1:end),[num_subj 1]);
        subj_inds = sum(subj_vars==req_vars,2)==num_bsf;
    end;
    subj_cell = subj_ids(subj_inds);
    % check within-subject factors
    if num_wsf == 0
        rpj = img_paths{1};
    elseif num_wsf == 1
        rpj = img_paths{cells(j,1)};
    elseif num_wsf == 2
        rpj = img_paths{cells(j,1),cells(j,2)};
    elseif num_wsf == 3
        rpj = img_paths{cells(j,1),cells(j,2),cells(j,3)};
    end;
    % collect subject scans
    subj_scans = cell(numel(subj_cell),1);
    for k = 1:numel(subj_cell)
        rpk = rpj;
        if contains(rpk,'*')
            rpk = strcat(rpk(1:strfind(rpk,'*')-1),subj_cell{k},rpk(strfind(rpk,'*')+1:end));
        end;
        if contains(rpk,'+')
            rpk = strcat(rpk(1:strfind(rpk,'+')-1),MS_name,rpk(strfind(rpk,'+')+1:end));
        end;
        if contains(rpk,'#')
            rpk = strcat(rpk(1:strfind(rpk,'#')-1),GLM_names{1},rpk(strfind(rpk,'#')+1:end));
        end;
        subj_scans{k} = strcat(bids_dir,rpk);
    end;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(j).scans = subj_scans;
end;
clear subj_cell subj_scans rpj rpk

% specify all the rest
matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts       = 1;
matlabbatch{1}.spm.stats.factorial_design.multi_cov              = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none     = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im             = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em             = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit         = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm        = 1;

% specify model estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat(1)        = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals  = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% specify contrast manager
matlabbatch{3}.spm.stats.con.spmmat(1)               = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.fcon.name    = 'EOI';
matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = eye(num_cell);
matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete                  = 0;
if ~isempty(con_info)
    for j = 1:numel(con_info)
        if strcmp(con_info(j).type,'F')
            matlabbatch{3}.spm.stats.con.consess{1+j}.fcon.name    = con_info(j).name;
            matlabbatch{3}.spm.stats.con.consess{1+j}.fcon.weights = con_info(j).vec;
            matlabbatch{3}.spm.stats.con.consess{1+j}.fcon.sessrep = 'none';
        end;
        if strcmp(con_info(j).type,'t')
            matlabbatch{3}.spm.stats.con.consess{1+j}.tcon.name    = con_info(j).name;
            matlabbatch{3}.spm.stats.con.consess{1+j}.tcon.weights = con_info(j).vec;
            matlabbatch{3}.spm.stats.con.consess{1+j}.tcon.sessrep = 'none';
        end;
    end;
end;

% save batch
if ~exist(ana_dir,'dir'), mkdir(ana_dir); end;
filename = strcat(ana_dir,'design.mat');
save(filename,'matlabbatch');

% display message
fprintf('\n');
fprintf('-> Thank you! The following files have been created:\n');
fprintf('   - SPM batch: %s.\n', strcat(ana_dir,'design.mat'));
fprintf('\n');

% run batch
if run_batch
    spm_jobman('run', filename);
end;