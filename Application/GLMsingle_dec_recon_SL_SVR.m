function GLMsingle_dec_recon_SL_SVR(SPM, rad, c, con, meth)
% _
% Decoding from Trials for Reconstruction (searchlight-based SVR)
% FORMAT GLMsingle_dec_recon_SL_SVR(SPM, rad, c, con, meth)
%     SPM  - a structure specifying an estimated GLM
%     rad  - a scalar specifying the searchlight radius in mm (e.g. 6)
%     c    - a 1 x p vector with +1s indicating variables to reconstruct or
%            a q x p matrix with one +1 in each row indexing each variable
%     con  - a string without spaces describing the contrast (e.g. 'PE')
%     meth - a string specifying trial-wise betas to use (i.e. 'GLM-single')
% 
% FORMAT GLMsingle_dec_recon_SL_SVR(SPM, rad, c, con, meth) performs searchlight
% decoding using support vector regression for reconstruction of selected
% variables indicated by c using searchlights of size rad.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 25/09/2024, 16:08 (Vn/a)
%  Last edit: 25/09/2024, 16:08 (Vn/a)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get SPM.mat if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    SPM_mat = spm_select(1,'^SPM\.mat$','Select SPM.mat!');
    SPM_dir = fileparts(SPM_mat); load(SPM_mat);
    SPM.swd = SPM_dir;
    GLMsingle_dec_recon_SL_SVR(SPM);
    return
end;

% Set searchlight radius if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(rad)
    rad = 2*abs(SPM.xY.VY(1).mat(1,1));
end;

% Set decoding contrast if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(c)
    c = zeros(1,0);
    k = 0;
    while isempty(c) && k < numel(SPM.Sess(1).U(1))
    	k = k + 1;
        if ~strcmp(SPM.Sess(1).U(k).P(1).name,'none')
            c(1,1) = SPM.Sess(1).U(k).P(1).i(2);
        end;
    end;
    if isempty(c)
        c(1,1) = 1;
    end;
end;

% Set contrast name if necessary
%-------------------------------------------------------------------------%
if nargin < 4 || isempty(con)
    if size(c,1) > 1, d = sum(c,1); else, d = c; end;
    con = '';
    ind = find(d);
    for l = 1:numel(ind)
        con = strcat(con,int2str(ind(l)));
        if l < numel(ind), con = strcat(con,','); end;
    end;
    clear d ind
end;

% Set beta method if necessary
%-------------------------------------------------------------------------%
if nargin < 5 || isempty(meth)
    meth = 'GLM-single';
end;

% Change to SPM.swd if specified
%-------------------------------------------------------------------------%
orig_dir = pwd;
try
    cd(SPM.swd);
catch
    SPM.swd = pwd;
end

% Get number of sessions
%-------------------------------------------------------------------------%
s = numel(SPM.Sess);

% Load GLM.mat in sub-directory
%-------------------------------------------------------------------------%
suff = meth; suff(4) = '_';
GLMs_dir = strcat(SPM.swd,'/',suff);
load(strcat(SPM.swd,'/','ITEM_est_1st_lvl','/','GLM1.mat'));

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','GLMsingle_dec_recon_SL_SVR: load');
spm_progress_bar('Init',100,'Load trial-wise parameter estimates...','');

% Load mask image
%-------------------------------------------------------------------------%
[M, m_dim, m_ind] = MA_load_mask(SPM);
[m_img, m_xyz]    = spm_read_vols(SPM.VM);
clear m_img

% Load restricted mask
%-------------------------------------------------------------------------%
mr_str = strcat(GLMs_dir,'/','mask_SVR.nii');
if exist(mr_str,'file')
    mr_hdr = spm_vol(mr_str);
    mr_dim = mr_hdr.dim;
   [mr_img, mr_xyz] = spm_read_vols(mr_hdr);
    mr_img = reshape(mr_img,[1 prod(mr_dim)]);
    mr_ind = find(mr_img~=0);
    Mr     = mr_img;
else
    Mr     = M;
    mr_ind = m_ind;
    mr_xyz = m_xyz;
end;
clear mr_str mr_hdr mr_img

% Load trial-wise responses
%-------------------------------------------------------------------------%
GLMs_mat = strcat(GLMs_dir,'/','TYPED_FITHRF_GLMDENOISE_RR.mat'); load(GLMs_mat);
G = double(squeeze(modelmd)');

% Get number of voxels
%-------------------------------------------------------------------------%
V = numel(m_ind);
v = numel(mr_ind);
d = floor(v/100);

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% E S T I M A T I O N   ( 1 ) :   P R E - C A L C U L A T I O N S         %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','GLMsingle_dec_recon_SL_SVR: estimate (1)');
spm_progress_bar('Init', 100, 'Determine searchlight voxels...', '');

% Get voxels per searchlight
%-------------------------------------------------------------------------%
SLs  = cell(v,1);
VpSL = NaN(size(Mr));
XYZ  = m_xyz(:,m_ind);
for j = 1:v
    xyz_cent = mr_xyz(:,mr_ind(j));
    vox_ind  = find(sqrt(sum((XYZ - repmat(xyz_cent,[1 V])).^2)) <= rad);
    SLs{j}   = vox_ind;
    VpSL(mr_ind(j)) = numel(vox_ind);
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
clear xyz_cent vox_ind

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% E S T I M A T I O N   ( 2 ) :   I N V E R S E   M O D E L               %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','GLMsingle_dec_recon_SL_SVR: estimate (2)');

% Augment decoding contrast if necessary
%-------------------------------------------------------------------------%
if size(c,1) > 1, c = sum(c,1); end;
c = [c, zeros(1, GLM1.p(1)-numel(c))];
q = numel(find(c));

% Cycle through recording sessions
%-------------------------------------------------------------------------%
for h = 1:s
    
    % "targets" - the X matrix
    %---------------------------------------------------------------------%
    Th = GLM1.Sess(h).T(1:GLM1.t(h),1:GLM1.p(h));
    Xh = Th(:,c==1);
    ITEM.Sess(h).X = Xh;
    clear Xh Th
    
    % "features" - gamma estimates
    %---------------------------------------------------------------------%
    ih = [sum(GLM1.t(1:(h-1)))+[1:GLM1.t(h)]];
    Yh = G(ih,:);                                       % standardize:
    Yh = Yh - repmat(mean(Yh),[size(Yh,1) 1]);          % - subtract mean
    Yh = Yh./ repmat(std(Yh), [size(Yh,1) 1]);          % - divide by std
    ITEM.Sess(h).Y = Yh;
    clear Yh
    
end;


%=========================================================================%
% E S T I M A T I O N   ( 3 ) :   C R O S S - V A L I D A T I O N         %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','GLMsingle_dec_recon_SL_SVR: estimate (3)');

% Cycle through cross-validation folds
%-------------------------------------------------------------------------%
for g = 1:s
    
    % Init progress bar
    %---------------------------------------------------------------------%
    spm_progress_bar('Init', 100, sprintf('Searchlight-based SVR analysis of session %d',g), '');
    
    % List sessions for this fold
    %---------------------------------------------------------------------%
    fold = 1:s;
    fold = fold(fold~=g);
    
    % Establish (in-sample) training data set
    %---------------------------------------------------------------------%
    X_in = vertcat(ITEM.Sess(fold).X);
    Y_in = vertcat(ITEM.Sess(fold).Y);
    
    % Establish (out-of-sample) test data set
    %---------------------------------------------------------------------%
    X_out = ITEM.Sess(g).X;
    Y_out = ITEM.Sess(g).Y;
    n_out = size(Y_out,1);
    
    % Perform searchlight-based SVR analysis
    %---------------------------------------------------------------------%
    C   = 1;                    % SVM cost parameter
    opt = sprintf('-s 4 -t 0 -c %s -q', num2str(C));
    Xp  = zeros(n_out,q,v);
    for j = 1:v                 % for each searchlight
        vj= SLs{j};
        for k = 1:q             % for each target
            SVM_in    = svmtrain(X_in(:,k), Y_in(:,vj), opt);
            Xp(:,k,j) = svmpredict(X_out(:,k), Y_out(:,vj), SVM_in, '-q');
        end;
        if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
    end;
    ITEM.Sess(g).Xp = Xp;
    clear X_in Y_in X_out Y_out SVM_in Xp
    
end;

% Remove feature matrix from ITEM structure
%-------------------------------------------------------------------------%
ITEM.Sess = rmfield(ITEM.Sess,'Y');

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% E S T I M A T I O N   ( 4 ) :   D E C O D I N G   A C C U R A C Y       %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','GLMsingle_dec_recon_SL_SVR: estimate (4)');

% Calculate out-of-sample correlation coefficients
%-------------------------------------------------------------------------%
% oosCC = NaN(q,numel(M),s);
% for g = 1:s
%     X_true  = ITEM.Sess(g).X;
%     X_recon = ITEM.Sess(g).Xp;
%     for k = 1:q
%         spm_progress_bar('Init', 100, sprintf('Calculate correlation coefficient for session %d, variable %d',g,k), '');
%         i_eff = find(X_true(:,k)~=0)';
%         for j = 1:v
%             oosCC(k,mr_ind(j),g) = corr(X_true(i_eff,k),X_recon(i_eff,k,j));
%             if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
%         end;
%     end;
% end;
% avgCC = mean(oosCC,3);
% clear X_true X_recon i_eff

% Calculate cross-validated correlation coefficients
%-------------------------------------------------------------------------%
cvCC    = NaN(q,numel(M));
X_true  = vertcat(ITEM.Sess(1:s).X);
X_recon = vertcat(ITEM.Sess(1:s).Xp);
for k = 1:q
    spm_progress_bar('Init', 100, sprintf('Calculate correlation coefficient across all sessions, variable %d',k), '');
    i_eff = find(X_true(:,k)~=0)';
    for j = 1:v
        cvCC(k,mr_ind(j)) = corr(X_true(i_eff,k),X_recon(i_eff,k,j));
        if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
    end;
end;
clear X_true X_recon i_eff

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','GLMsingle_dec_recon_SL_SVR: save');

% Initialize image files
%-------------------------------------------------------------------------%
i = find(c);
H = MA_init_header(SPM, false);
ITEM.swd = strcat(SPM.swd,'/',suff,'_SVR','/','SVR_',con,'_SL-',num2str(rad),'mm','/');
if ~exist(ITEM.swd,'dir'), mkdir(ITEM.swd); end;
cd(ITEM.swd);

% Save out-of-sample correlations
%-------------------------------------------------------------------------%
% d = floor(log10(s))+1;
% for h = 1:s
%     for k = 1:q
%         H.fname   = strcat('oosCC_',MF_int2str0(i(k),4),'_S',MF_int2str0(h,d),'.nii');
%         H.descrip = sprintf('GLMsingle_dec_recon_SL_SVR: out-of-sample correlation coefficient; session %d, target %d', h, i(k));
%         spm_write_vol(H,reshape(oosCC(k,:,h),m_dim));
%         ITEM.VoosCC(h,k) = H;
%     end;
% end;

% Save averaged correlations
%-------------------------------------------------------------------------%
% for k = 1:q
%     H.fname   = strcat('avgCC_',MF_int2str0(i(k),4),'.nii');
%     H.descrip = sprintf('GLMsingle_dec_recon_SL_SVR: averaged correlation coefficient; %d sessions, target %d', s, i(k));
%     spm_write_vol(H,reshape(avgCC(k,:),m_dim));
%     ITEM.VavgCC(1,k) = H;
% end;

% Save cross-validated correlations
%-------------------------------------------------------------------------%
for k = 1:q
    H.fname   = strcat('cvCC_',MF_int2str0(i(k),4),'.nii');
    H.descrip = sprintf('GLMsingle_dec_recon_SL_SVR: cross-validated correlation coefficient; %d sessions, target %d', s, i(k));
    spm_write_vol(H,reshape(cvCC(k,:),m_dim));
    ITEM.VcvCC(1,k) = H;
end;

% Save voxels per searchlight
%-------------------------------------------------------------------------%
H.fname   = strcat('VpSL.nii');
H.descrip = sprintf('GLMsingle_dec_recon_SL_SVR: voxels per searchlight; %s mm radius', num2str(rad));
spm_write_vol(H,reshape(VpSL,m_dim));
ITEM.VVpSL = H;

% Complete ITEM structure
%-------------------------------------------------------------------------%
ITEM.rad       = rad;
ITEM.SLs       = SLs;
ITEM.Recon.c   = c;
ITEM.Recon.con = con;
ITEM.Sess      = rmfield(ITEM.Sess,'Xp');

% Save ITEM structure
%-------------------------------------------------------------------------%
save(strcat(ITEM.swd,'ITEM.mat'),'ITEM','-v7.3');

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);