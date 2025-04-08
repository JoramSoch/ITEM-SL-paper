function GLMsingle_est_1st_lvl(SPM, nind)
% _
% Estimate Single-Trial Responses using GLMsingle
% FORMAT GLMsingle_est_1st_lvl(SPM, nind)
%     SPM  - a structure specifying an estimated GLM
%     nind - a vector indexing non-trial-wise nuisance conditions
% 
% FORMAT GLMsingle_est_1st_lvl(SPM, nind) estimates trial-wise responses
% using GLMsingle, leaving out nuisance conditions of no interest.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 25/09/2024, 08:10 (V0.4)
%  Last edit: 07/03/2025, 14:05 (V0.4)


%%% directories and files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get directories
orig_dir = pwd;
SPM_dir  = SPM.swd;
GLMs_dir = strcat(SPM_dir,'/GLM_single/');


%%% GLMsingle trial-wise estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% configure data
S     = numel(SPM.Sess);
Y     = MA_load_data(SPM);
data  = cell(1,S);
for h = 1:S
    data{h} = Y(SPM.Sess(h).row,:)';
    % disp(data{h});
end;
clear Y

% configure timing
dur = SPM.Sess(1).U(1).dur(1);
TR  = SPM.xY.RT;

% configure design
design= cell(1,S);
for h = 1:S
    n = numel(SPM.Sess(h).row);
    p = numel(SPM.Sess(h).U);
    X = zeros(n,p);
    for j = 1:p
        if ~ismember(j,nind)
            for k = 1:numel(SPM.Sess(h).U(j).ons)
                i = round(SPM.Sess(h).U(j).ons(k)/TR)+1;
                X(i,j) = 1;
            end;
        end;
    end;
    X = X(:,sum(X)>0);
    design{h} = X;
    % disp(size(X));
end;
clear X

% configure covariates
regs  = cell(1,S);
for h = 1:S
    Xr= [];
    Xh= SPM.xX.X(SPM.Sess(h).row, SPM.Sess(h).col);
    for j = nind
        Xr = [Xr, Xh(:,SPM.Sess(h).Fc(j).i)];
    end;
    Xr= [Xr, SPM.Sess(h).C.C];
    regs{h} = Xr;
    % disp(size(Xr));
end;
clear Xr

% configure options
options.wantlibrary       = 1;
options.wantglmdenoise    = 1;
options.wantfracridge     = 1;
options.wantfileoutputs   = [false, false, false, true];
options.wantmemoryoutputs = [false, false, false, false];
options.extraregressors   = regs;

% estimate trials
[results, design] = GLMestimatesingletrial(design, data, dur, TR, GLMs_dir, options);
cd(orig_dir);