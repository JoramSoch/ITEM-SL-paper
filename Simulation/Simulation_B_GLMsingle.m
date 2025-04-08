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
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% Version History:
% - 29/11/18,22/07/19: univariate simulation
% - 24/04/2023, 22:35: multivariate simulation
% - 31/10/2023, 13:04: prepared for upload
% - 29/11/2023, 09:43: finalized for upload
% - 12/11/2024, 10:04: additional simulation
% - 28/11/2024, 01:59: added GLMsingle
% - 07/03/2025, 08:09: modified GLMsingle
% - 07/03/2025, 08:38: removed LS-A, LS-S, ITEM
% - 07/03/2025, 10:46: standardized features
% - 08/04/2025, 16:20: finalized for upload


clear
close all

%%% Step 1: set simulation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n-> Step 1: ');

% set parameters
rng(1);
N   = 100;                      % number of simulations
S   = 2;                        % number of sessions per test
t   =[100, 200, 400];           % number of trials per session
v   =[16, 32, 64];              % number of voxels per trial
r   =[0:0.05:1];                % proportion of voxels with information
mu  = 0;                        % mean of betas across voxels
s2v = 1;                        % between-voxel variance
s2t = 0.5^2;                    % between-trial variance
s2n = 1.6;                      % between-scan variance
t0  = 10;                       % stimulus onset
dt  = 2;                        % stimulus duration
TR  = 2;                        % repetition time
ISI =[0,4];                     % inter-stimulus-intervals
rho = 0.12;                     % temporal auto-correlation
ny  = 0.48;                     % spatial auto-correlation

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

% for each SL size
for g = 1:numel(v)
    
    fprintf('v = %d: t = ', v(g));
    
    % for each sample size
    for h = 1:numel(t)
        
        fprintf('%d, ', t(h));
        
        %%% Step 2a: set design fundamentals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % calculate number of scans needed for ISI
        n = ceil((t0 + t(h)*dt + (t(h)-1)*((ISI(1)+ISI(2))/2) + 5*t0)/TR);
        z = zeros(n,1);
        Sim(g,h).base.n = n;
        
        % create temporal auto-correlation accordingly
        V = toeplitz(rho.^[0:1:(n-1)]);
        W = sqrtm(inv(V));
        Sim(g,h).base.V = V;
        Sim(g,h).base.W = W;
        
        % multiply with variance to get covariance matrix
        s2V = s2n^2 * V;
        Sim(g,h).base.s2V = s2V;
        
        % configure settings for HRF convolution
        settings.n  = n;
        settings.TR = TR;
        
        %%% Step 2b: create design matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % for each session
        for j = 1:S
        
            % sample trial types
            ttm = [[1*ones(t(h)/2,1); 2*ones(t(h)/2,1)], rand(t(h),1)];
            ttm = sortrows(ttm, 2);
            tt  = ttm(:,1);
            Sim(g,h).des(j).tt = tt;

            % sample trial onsets
            isi = MD_unirnd(ISI(1), ISI(2), t(h)-1)';
            ons = [t0, t0 + cumsum(isi) + dt*[1:(t(h)-1)]]';
            dur = dt*ones(t(h),1);
            Sim(g,h).des(j).isi = isi';
            Sim(g,h).des(j).ons = ons;
            Sim(g,h).des(j).dur = dur;

            % create design matrix X_S (Mumford: LS-A)
            clear names onsets durations
            for k = 1:t(h)
                names{k}     = strcat('trl-',MF_int2str0(k,2));
                onsets{k}    = ons(k);
                durations{k} = dur(k);
            end;
            [X,L] = ITEM_get_des_mat(names, onsets, durations, [], [], [], settings);
            Sim(g,h).des(j).DM.X_t = X;

            % create design matrices X_T (Mumford: LS-S)
            X_t = X;
            for k = 1:t(h)
                T = zeros(t(h),2);
                T(          k,1) = 1;
                T([1:t(h)]~=k,2) = 1;
                X = X_t * T;
                Sim(g,h).des(j).DM.X_s{k} = X;
            end;

            % create design matrix T (Soch: LS-T)
            X_t = Sim(g,h).des(j).DM.X_t;
            U   = (X_t' * V^(-1) * X_t)^(-1);
            T   = zeros(t(h),2);
            for k = 1:2
                T(tt==k,k) = 1;
            end;
            X = X_t * T;
            Sim(g,h).des(j).DM.X = X;
            Sim(g,h).des(j).DM.T = T;
            Sim(g,h).des(j).DM.U = U;
        
        end;
        
        %%% Step 2c: sample data matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % for each voxel proportion
        for i = 1:numel(r)
            
            % for each session
            for j = 1:S

                % create identity matrices
                ov = ones(2,v(g));
                ot = ones(t(h)/2,1);
                zn = zeros(n,v(g));
                Sv = toeplitz(ny.^[0:1:(v(g)-1)]);
                Iv = eye(v(g));
                It = eye(t(h)/2);
                In = eye(n);

                % sample voxel parameters
                if j == 1
                    % sample voxel-wise beta parameters
                    b = MD_matnrnd(mu*ov, s2v*eye(2), Iv, N);
                    for s = 1:N
                        % eliminate differences from
                        % voxels without information
                        k = find(rand(1,v(g))>r(i));
                        b(2,k,s) = b(1,k,s);
                    end;
                else
                    b = Sim(g,h).data(1,i).b;
                end;
                Sim(g,h).data(j,i).b = b;

                % sample trial parameters
                tt = Sim(g,h).des(j).tt;
                B  = zeros(t(h),v(g),N);
                for s = 1:N
                    b1 = MD_matnrnd(ot*b(1,:,s), s2t*It, Iv);
                    b2 = MD_matnrnd(ot*b(2,:,s), s2t*It, Iv);
                    B(tt==1,:,s) = b1;
                    B(tt==2,:,s) = b2;
                end;
                Sim(g,h).data(j,i).B = B;
                clear b1 b2

                % sample fMRI data
                X = Sim(g,h).des(j).DM.X_t;
                E = MD_matnrnd(zn, s2V, Sv, N);
                Y = zeros(n,v(g),N);
                for s = 1:N
                    Y(:,:,s) = X*B(:,:,s) + E(:,:,s);
                end;
                Sim(g,h).data(j,i).Y = Y;
                clear E Y

            end;
            
        end;
        
    end;
    
    fprintf('done; ');
    
end;

fprintf('end.');


%%% Step 3: estimate models and test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fprintf('\n\n-> Step 3: ');

% for each SL size
for g = 1:numel(v)
    
    fprintf('v = %d: ', v(g));
    
    % for each sample size
    for h = 1:numel(t)
        
        fprintf('t = %d: \n ', t(h));
        
        %%% Step 3a: estimate model parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % for each voxel proportion
        for i = 1:numel(r)
            
            % configure timing for GLMsingle
            n   = Sim(g,h).base.n;
            p   = max(Sim(g,h).des(1).tt);
            dur = mean(Sim(g,h).des(1).dur);
            %TR = TR;
            
            % configure data for GLMsingle
            data  = cell(1,S);
            for j = 1:S
                Y = zeros(n,v(g)*N);
                for s = 1:N
                    js= [(s-1)*v(g)+[1:v(g)]];
                    Y(:,js) = Sim(g,h).data(j,i).Y(:,:,s);
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
                for s = 1:N
                    ij= [(j-1)*t(h)+[1:t(h)]];
                    js= [(s-1)*v(g)+[1:v(g)]];
                    BG= bG(ij,js);                                  % standardize:
                    BG= BG - repmat(mean(BG),[size(BG,1) 1]);       %- subtract mean
                    BG= BG./ repmat(std(BG), [size(BG,1) 1]);       %- divide by std
                    Sim(g,h).est(j,i).BG(:,:,s) = BG;
                end;
            end;
            clear BG
            
        end;
        
        %%% Step 3c: decode/classify trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % for each voxel proportion
        for i = 1:numel(r)
            
            % for each session
            for j = 1:S
                
                % get train/test labels
                x_train = Sim(g,h).des([1:S]~=j).tt;
                x_test  = Sim(g,h).des(j).tt;
                
                % for each simulation
                for s = 1:N
                    
                    % specify SVM options
                    opt = sprintf('-s 0 -t 0 -c %s -q', num2str(C));
                    
                    % SVC using GLMsingle parameter estimates
                    Y_train = Sim(g,h).est([1:S]~=j,i).BG(:,:,s);
                    Y_test  = Sim(g,h).est(j,i).BG(:,:,s);
                    svm_tr  = svmtrain(x_train, Y_train, opt);
                    x_pred  = svmpredict(x_test, Y_test, svm_tr, '-q');
                    a = sum(x_test==x_pred)/t(h);
                    Sim(g,h).pred.aG(s,i,j) = a;
                    
                end;

            end;
            
        end;
        
    end;
    
    fprintf('done; ');
    
end;

fprintf('end.');
total_time = toc;


%%% Step 4: compute results summaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n-> Step 4: ');
Y_11  = Sim(1,1).data(1,1).Y(:,:,1);
BG_11 = Sim(1,1).est(1,1).BG(:,:,1);

% for each SL size
for g = 1:numel(v)
    
    fprintf('v = %d: ', v(g));
    
    % for each sample size
    for h = 1:numel(t)
        
        fprintf('t = %d, ', t(h));
        
        % true positive rates
        for i = 1:numel(r)
            nT = S*t(h);
            cG = round(sum(Sim(g,h).pred.aG(:,i,:)*t(h),3));
            [phat, pciG]      = binofit(cG, nT, alpha);
            Res(g,h).TPR(:,i) = [zeros(1,3), mean(pciG(:,1)>chlvl)];
        end;
        clear nT cA cS cT cG phat pci*
        
        % decoding accuracies
        Res(g,h).DA(:,:,1) = zeros(N,numel(r));
        Res(g,h).DA(:,:,2) = zeros(N,numel(r));
        Res(g,h).DA(:,:,3) = zeros(N,numel(r));
        Res(g,h).DA(:,:,4) = mean(Sim(g,h).pred.aG,3);
        
        % remove fields
        Sim(g,h).data = rmfield(Sim(g,h).data,{'B','Y'});
        Sim(g,h).est  = rmfield(Sim(g,h).est,{'BG'});
        
    end;
    
end;

save('Simulation_B_GLMsingle.mat', 'Sim', 'Res', 'total_time');

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

for g = 1:numel(v)
    for h = 1:numel(t)
        subplot(numel(v), numel(t), (g-1)*numel(t)+h);
        hold on;
        plot(r, Res(g,h).TPR(1,:), '-r');
        plot(r, Res(g,h).TPR(2,:), '-b');
        plot(r, Res(g,h).TPR(4,:), '-m');
        plot(r, Res(g,h).TPR(3,:), '-g');
        axis([(0-0.01), (1+0.01), (0.5-0.01), (1+0.01)]);
        set(gca,'Box','On');
        if g == numel(v) && h == numel(t)
            legend({'LS-A', 'LS-S', 'ITEM', 'GLMsingle'}, 'Location', 'SouthEast');
        end;
        if g == numel(v) && h == 2
            xlabel('proportion of voxels with signal carrying information about experimental conditions', 'FontSize', 14);
        end;
        if h == 1
            yl  = ylabel('true positive rate', 'FontSize', 14);
            pos = get(yl,'Position');
            text(pos(1)-1/5, pos(2), sprintf('v = %d', v(g)), 'Rotation', 90, 'FontSize', 16, ...
                 'FontWeight', 'bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
            clear yl pos
        end;
        if g == 1
            title(sprintf('t = %d', t(h)), 'FontSize', 16);
        end;
    end;
end;


% plot decoding accuracies
figure('Name', 'decoding accuracies', 'Color', [1 1 1], 'Position', [50 50 1000 1000]);

for g = 1:numel(v)
    for h = 1:numel(t)
        subplot(numel(v), numel(t), (g-1)*numel(t)+h);
        hold on;
        plot(r, mean(Res(g,h).DA(:,:,1)), '-r', 'LineWidth', 2);
        plot(r, mean(Res(g,h).DA(:,:,2)), '-b', 'LineWidth', 2);
        plot(r, mean(Res(g,h).DA(:,:,3)), '-g', 'LineWidth', 2);
        plot(r, mean(Res(g,h).DA(:,:,4)), '-m', 'LineWidth', 2);
        plot(r, mean(Res(g,h).DA(:,:,1))+std(Res(g,h).DA(:,:,1)), '--r');
        plot(r, mean(Res(g,h).DA(:,:,1))-std(Res(g,h).DA(:,:,1)), '--r');
        plot(r, mean(Res(g,h).DA(:,:,2))+std(Res(g,h).DA(:,:,2)), '--b');
        plot(r, mean(Res(g,h).DA(:,:,2))-std(Res(g,h).DA(:,:,2)), '--b');
        plot(r, mean(Res(g,h).DA(:,:,3))+std(Res(g,h).DA(:,:,3)), '--g');
        plot(r, mean(Res(g,h).DA(:,:,3))-std(Res(g,h).DA(:,:,3)), '--g');
        plot(r, mean(Res(g,h).DA(:,:,4))+std(Res(g,h).DA(:,:,4)), '--m');
        plot(r, mean(Res(g,h).DA(:,:,4))-std(Res(g,h).DA(:,:,4)), '--m');
        axis([(0-0.01), (1+0.01), (0.5-0.01), (1+0.01)]);
        set(gca,'Box','On');
        if g == numel(v) && h == numel(t)
            legend({'LS-A', 'LS-S', 'GLMsingle', 'ITEM'}, 'Location', 'SouthEast');
        end;
        if g == numel(v) && h == 2
            xlabel('proportion of voxels with signal carrying information about experimental conditions', 'FontSize', 14);
        end;
        if h == 1
            yl  = ylabel('decoding accuracy', 'FontSize', 14);
            pos = get(yl,'Position');
            text(pos(1)-1/5, pos(2), sprintf('v = %d', v(g)), 'Rotation', 90, 'FontSize', 16, ...
                 'FontWeight', 'bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
            clear yl pos
        end;
        if g == 1
            title(sprintf('t = %d', t(h)), 'FontSize', 16);
        end;
    end;
end;