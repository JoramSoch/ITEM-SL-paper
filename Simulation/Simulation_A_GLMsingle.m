% ITEM-SL: Simulation
% _
% This script performs a conceptual analogue of the simulation in Mumford
% et al. (2012), extended to multivariate signals and supplemented with an
% inverse transformed encoding model (ITEM) approach.
% 
% In the entire script, we use the following suffixes:
% - "*G": fractional ridge regression a.k.a.
%         Prince & Kay's "GLM single" (GLMs)
% 
% Author: Joram Soch
% E-Mail: joram.soch@bccn-berlin.de
% 
% Version History:
% - 29/11/18,22/07/19: univariate simulation
% - 24/04/2023, 22:35: multivariate simulation
% - 28/11/2024, 10:55: added GLMsingle
% - 07/03/2025, 07:39: modified GLMsingle
% - 07/03/2025, 08:18: removed LS-A, LS-S, ITEM
% - 07/03/2025, 10:41: standardized features
% - 08/04/2025, 16:18: finalized for upload


clear
close all

%%% Step 1: set simulation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n-> Step 1: ');

% set parameters
rng(1);
N   = 1000;                     % number of simulations
S   = 2;                        % number of sessions per test
t   = 100;                      % number of trials per session
v   = 33;                       % number of voxels per trial [1]
r   = 1/5;                      % proportion of voxels with information
mu  = 0;                        % mean of betas across voxels
s2v = 1;                        % between-voxel variance
s2t = 0.5^2;                    % between-trial variance
s2n = [0.8, 1.6, 3.2];          % between-scan variance
t0  = 10;                       % stimulus onset
dt  = 2;                        % stimulus duration
TR  = 2;                        % repetition time
ISI = [0,4; 2,6; 4,8];          % inter-stimulus-intervals
rho = 0.12;                     % temporal auto-correlation
ny  = 0.48;                     % spatial auto-correlation
% [1] https://oeis.org/A000605:  1, 7, 33, 123, 257

% specify settings
C     = 1;                      % cost parameter for support vector classification
chlvl = 0.5;                    % chance level for classification task
alpha = 0.05;                   % significance level for binomial test

% preallocate results
Sim   = struct([]);             % simulations structure
Res   = struct([]);             % sim results structure

fprintf('end.');


%%% Step 2: generate designs and data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n-> Step 2: ');

% for each noise level
for g = 1:numel(s2n)
    
    fprintf('%1.1f: ', s2n(g));
    
    % for each ISI range
    for h = 1:size(ISI,1)
        
        fprintf('[%d,%d], ', ISI(h,1), ISI(h,2));
        
        %%% Step 2a: set design fundamentals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % calculate number of scans needed for ISI
        n = ceil((t0 + t*dt + (t-1)*((ISI(h,1)+ISI(h,2))/2) + 5*t0)/TR);
        z = zeros(n,1);
        Sim(g,h).base.n = n;
        
        % create temporal auto-correlation accordingly
        V = toeplitz(rho.^[0:1:(n-1)]);
        W = sqrtm(inv(V));
        Sim(g,h).base.V = V;
        Sim(g,h).base.W = W;
        
        % multiply with variance to get covariance matrix
        s2V = s2n(g)^2 * V;
        Sim(g,h).base.s2V = s2V;
        
        % configure settings for HRF convolution
        settings.n  = n;
        settings.TR = TR;
        
        %%% Step 2b: create design matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % for each session
        for j = 1:S
        
            % sample trial types
            ttm = [[1*ones(t/2,1); 2*ones(t/2,1)], rand(t,1)];
            ttm = sortrows(ttm, 2);
            tt  = ttm(:,1);
            Sim(g,h).des(j).tt = tt;

            % sample trial onsets
            isi = MD_unirnd(ISI(h,1), ISI(h,2), t-1)';
            ons = [t0, t0 + cumsum(isi) + dt*[1:(t-1)]]';
            dur = dt*ones(t,1);
            Sim(g,h).des(j).isi = isi';
            Sim(g,h).des(j).ons = ons;
            Sim(g,h).des(j).dur = dur;

            % create design matrix X_S (Mumford: LS-A)
            clear names onsets durations
            for k = 1:t
                names{k}     = strcat('trl-',MF_int2str0(k,2));
                onsets{k}    = ons(k);
                durations{k} = dur(k);
            end;
            [X,L] = ITEM_get_des_mat(names, onsets, durations, [], [], [], settings);
            Sim(g,h).des(j).DM.X_t = X;

            % create design matrices X_T (Mumford: LS-S)
            X_t = X;
            for k = 1:t
                T = zeros(t,2);
                T(       k,1) = 1;
                T([1:t]~=k,2) = 1;
                X = X_t * T;
                Sim(g,h).des(j).DM.X_s{k} = X;
            end;

            % create design matrix T (Soch: LS-T)
            X_t = Sim(g,h).des(j).DM.X_t;
            U   = (X_t' * V^(-1) * X_t)^(-1);
            T   = zeros(t,2);
            for k = 1:2
                T(tt==k,k) = 1;
            end;
            X = X_t * T;
            Sim(g,h).des(j).DM.X = X;
            Sim(g,h).des(j).DM.T = T;
            Sim(g,h).des(j).DM.U = U;
        
        end;
        
        %%% Step 2c: sample data matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % for each session
        for j = 1:S
            
            % create identity matrices
            ov = ones(2,v);
            ot = ones(t/2,1);
            zn = zeros(n,v);
            Sv = toeplitz(ny.^[0:1:(v-1)]);
            Iv = eye(v);
            It = eye(t/2);
            In = eye(n);
            
            % sample voxel parameters
            if j == 1
                % sample voxel-wise beta parameters
                b = MD_matnrnd(mu*ov, s2v*eye(2), Iv, N);
                for i = 1:N
                    % eliminate differences from
                    % voxels without information
                    k = find(rand(1,v)>r);
                    b(2,k,i) = b(1,k,i);
                end;
            else
                b = Sim(g,h).data(1).b;
            end;
            Sim(g,h).data(j).b = b;
            
            % sample trial parameters
            tt = Sim(g,h).des(j).tt;
            B  = zeros(t,v,N);
            for i = 1:N
                b1 = MD_matnrnd(ot*b(1,:,i), s2t*It, Iv);
                b2 = MD_matnrnd(ot*b(2,:,i), s2t*It, Iv);
                B(tt==1,:,i) = b1;
                B(tt==2,:,i) = b2;
            end;
            Sim(g,h).data(j).B = B;
            clear b1 b2

            % sample fMRI data
            X = Sim(g,h).des(j).DM.X_t;
            E = MD_matnrnd(zn, s2V, Sv, N);
            Y = zeros(n,v,N);
            for i = 1:N
                Y(:,:,i) = X*B(:,:,i) + E(:,:,i);
            end;
            Sim(g,h).data(j).Y = Y;
            clear E Y

        end;
        
    end;
    
    fprintf('done; ');
    
end;

fprintf('end.');


%%% Step 3: estimate models and test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fprintf('\n\n-> Step 3: ');

% for each noise level
for g = 1:numel(s2n)
    
    fprintf('%1.1f: ', s2n(g));
    
    % for each ISI range
    for h = 1:size(ISI,1)
        
        fprintf('[%d,%d]: \n', ISI(h,1), ISI(h,2));
        
        %%% Step 3a: estimate model parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % configure timing for GLMsingle
        n   = Sim(g,h).base.n;
        p   = max(Sim(g,h).des(1).tt);
        dur = mean(Sim(g,h).des(1).dur);
        %TR = TR;
        
        % configure data for GLMsingle
        data  = cell(1,S);
        for j = 1:S
            Y = zeros(n,v*N);
            for i = 1:N
                ji= [(i-1)*v+[1:v]];
                Y(:,ji) = Sim(g,h).data(j).Y(:,:,i);
            end;
            data{j} = Y';
        end;
        
        % configure design for GLMsingle
        design= cell(1,S);
        for j = 1:S
            X = zeros(n,p);
            tt= Sim(g,h).des(j).tt;
            for k = 1:numel(tt)
                ik= round(Sim(g,h).des(j).ons(k)/TR)+1;
                X(ik,tt(k)) = 1;
            end;
            design{j} = X;
        end;

        % perform GLMsingle estimation
        orig_dir = pwd;
        GLMs_dir = {NaN, NaN};
        options.wantlibrary       = 0;
        options.wantglmdenoise    = 0;
        options.wantfracridge     = 1;
        options.wantfileoutputs   = [false, false, false, false];
        options.wantmemoryoutputs = [false, false, false, true];
        [results, design] = GLMestimatesingletrial(design, data, dur, TR, GLMs_dir, options);
        cd(orig_dir);

        % extract GLMsingle estimates
        bG = double(squeeze(results{4}.modelmd)');
        for j = 1:S
            for i = 1:N
                ij= [(j-1)*t+[1:t]];
                ji= [(i-1)*v+[1:v]];
                BG= bG(ij,ji);                          % standardize:
                BG= BG - repmat(mean(BG),[size(BG,1) 1]);%- subtract mean
                BG= BG./ repmat(std(BG), [size(BG,1) 1]);%- divide by std
                Sim(g,h).est(j).BG(:,:,i) = BG;
            end;
        end;
        clear BG
        
        %%% Step 3c: decode/classify trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % for each session
        for j = 1:S
            
            % get train/test labels
            x_train = Sim(g,h).des([1:S]~=j).tt;
            x_test  = Sim(g,h).des(j).tt;
            
            % for each simulation
            for i = 1:N
            
                % specify SVM options
                opt = sprintf('-s 0 -t 0 -c %s -q', num2str(C));
                
                % SVC using GLMsingle parameter estimates
                x_train = Sim(g,h).des([1:S]~=j).tt;
                x_test  = Sim(g,h).des(j).tt;
                Y_train = Sim(g,h).est([1:S]~=j).BG(:,:,i);
                Y_test  = Sim(g,h).est(j).BG(:,:,i);
                svm_tr  = svmtrain(x_train, Y_train, opt);
                x_pred  = svmpredict(x_test, Y_test, svm_tr, '-q');
                a = sum(x_test==x_pred)/t;
                Sim(g,h).pred.aG(j,i) = a;

            end;
            
        end;
        
    end;
    
    fprintf('done; ');
    
end;

fprintf('end.');
total_time = toc;


%%% Step 4: compute results summaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n-> Step 4: ');
Y_11  = Sim(1,1).data(1).Y(:,:,1);
BG_11 = Sim(1,1).est(1).BG(:,:,1);

% for each noise level
for g = 1:numel(s2n)
    
    % for each ISI range
    for h = 1:size(ISI,1)
        
        % true positive rates
        nT = S*t;
        cG = round(sum(Sim(g,h).pred.aG*t,1));
        [phat, pciG] = binofit(cG, nT, alpha);
        Res(g,h).TPR = [zeros(1,3), mean(pciG(:,1)>chlvl)];
        clear nT cG phat pciG
        
        % decoding accuracies
        Res(g,h).DA = [zeros(3,N); ...
                       mean(Sim(g,h).pred.aG)];
                   
        % remove fields
        Sim(g,h).data = rmfield(Sim(g,h).data,{'B','Y'});
        Sim(g,h).est  = rmfield(Sim(g,h).est,{'BG'});
        
    end;
    
end;

save('Simulation_A_GLMsingle.mat', 'Sim', 'Res', 'total_time');

fprintf('end.');


%%% Step 5: plot simulation results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n');


% plot design matrices
figure('Name', 'design matrices', 'Color', [1 1 1], 'Position', [50 50 900 900]);

subplot(5,3,[1 4 7]);
imagesc(Sim(1,1).des(1).DM.X); axis off;
title('X_{ } [n x p]', 'FontSize', 16);

subplot(5,3,[2 5 8]);
imagesc(Sim(1,1).des(1).DM.X_t); axis off;
title('X_t [n x t]', 'FontSize', 16);

subplot(5,3,[3 6 9]);
imagesc(Sim(1,1).des(1).DM.T); axis off;
title('T_{ } [t x p]', 'FontSize', 16);

subplot(5,3,[11 14]);
imagesc(Sim(1,1).des(1).DM.U); axis off; axis square;
title('U = (X_{t}^{T}V^{-1}X_{t})^{-1}', 'FontSize', 16);


% plot linear models
figure('Name', 'linear models', 'Color', [1 1 1], 'Position', [50 50 900 900]);

% standard model
subplot(3,4,1);
axis off
text(1/2, 1/2, 'standard model', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');

subplot(3,4,2);
imagesc(Y_11); axis off;
title('Y', 'FontSize', 16);

subplot(3,4,3);
imagesc(Sim(1,1).des(1).DM.X); axis off;
title('X', 'FontSize', 16);

subplot(3,4,4);
imagesc(Sim(1,1).base.V); axis off; axis square;
title('V', 'FontSize', 16);

% first-level model
subplot(3,4,5);
axis off
text(1/2, 1/2, '"first-level" model', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');

subplot(3,4,6);
imagesc(Y_11); axis off;
title('Y', 'FontSize', 16);

subplot(3,4,7);
imagesc(Sim(1,1).des(1).DM.X_t); axis off;
title('X_t', 'FontSize', 16);

subplot(3,4,8);
imagesc(Sim(1,1).base.V); axis off; axis square;
title('V', 'FontSize', 16);

% second-level model
subplot(3,4,9);
axis off
text(1/2, 1/2, '"second-level" model', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');

subplot(3,4,10);
imagesc(BG_11); axis off;
title('\gamma', 'FontSize', 16);

subplot(3,4,11);
imagesc(Sim(1,1).des(1).DM.T); axis off;
title('T', 'FontSize', 16);

subplot(3,4,12);
imagesc(Sim(1,1).des(1).DM.U); axis off; axis square;
title('U', 'FontSize', 16);


% plot test performances
figure('Name', 'test performances', 'Color', [1 1 1], 'Position', [50 50 1000 1000]);

for g = 1:numel(s2n)
    for h = 1:size(ISI,1)
        subplot(numel(s2n), size(ISI,1), (g-1)*size(ISI,1)+h);
        hold on;
        bar(1, Res(g,h).TPR(1), 'r');
        bar(2, Res(g,h).TPR(2), 'b');
        bar(3, Res(g,h).TPR(3), 'g');
        bar(4, Res(g,h).TPR(4), 'm');
        if r > 0
            axis([(1-0.5), (4+0.5), 0.5, 1.01]);
        else
            axis([(1-0.5), (4+0.5), -0.01, 0.5]);
        end;
        set(gca,'Box','On');
        set(gca,'XTick',[1:4],'XTickLabel',{'LS-A', 'LS-S', 'ITEM', 'GLMs'});
        if g == numel(s2n)
            xlabel(['t_{isi} \sim', sprintf(' U(%d,%d)', ISI(h,1), ISI(h,2))], 'FontSize', 16);
        end;
        if h == 1
            ylabel(['\sigma^2 =', sprintf(' %1.1f', s2n(g))], 'FontSize', 16);
        end;
    end;
end;


% plot decoding accuracies
figure('Name', 'decoding accuracies', 'Color', [1 1 1], 'Position', [50 50 1000 1000]);

for g = 1:numel(s2n)
    for h = 1:size(ISI,1)
        subplot(numel(s2n), size(ISI,1), (g-1)*size(ISI,1)+h);
        hBP = boxplot(Res(g,h).DA', 'Positions', [1:4], 'Width', 2/3, ...
                     'Colors', 'rbgm', 'Symbol', '+k', ...
                     'Labels', {'LS-A', 'LS-S', 'ITEM', 'GLMs'});
        set(hBP,'LineWidth',2);
        if r > 0
            axis([(1-0.5), (4+0.5), 0.5, 1]);
        else
            axis([(1-0.5), (4+0.5), 0.25, 0.75]);
        end;
        set(gca,'Box','On');
        if g == numel(s2n)
            xlabel(['t_{isi} \sim', sprintf(' U(%d,%d)', ISI(h,1), ISI(h,2))], 'FontSize', 16);
        end;
        if h == 1
            ylabel(['\sigma^2 =', sprintf(' %1.1f', s2n(g))], 'FontSize', 16);
        end;
    end;
end;