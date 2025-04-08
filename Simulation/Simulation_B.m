% ITEM-SL: Simulation B
% _
% This script performs a conceptual analogue of the simulation in Mumford
% et al. (2012), extended to multivariate signals and supplemented with an
% inverse transformed encoding model (ITEM) approach.
% 
% In the entire script, we use the following suffixes:
% - "*A": trial-wise design matrix X_S a.k.a.
%         Mumford's "least squares, all" (LS-A)
% - "*S": trial-based design matrix X_T a.k.a.
%         Mumford's "least squares, separate" (LS-S)
% - "*T": (inverse) transformed encoding model a.k.a.
%         Soch's "least squares, transformed" (LS-T)
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% Version History:
% - 29/11/18,22/07/19: univariate simulation
% - 24/04/2023, 22:35: multivariate simulation
% - 12/11/2024, 10:04: rewritten to Simulation B
% - 18/12/2024, 17:34: finalized for upload
% - 27/03/2023, 16:46: added time measurement


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
ReML  = true;                   % execution of restricted maximum likelihood
C     = 1;                      % cost parameter for support vector classification
chlvl = 0.5;                    % chance level for classification task
alpha = 0.05;                   % significance level for binomial test

% preallocate results
Sim    = struct([]);            % simulations structure
Res    = struct([]);            % sim results structure
TPR_bA = zeros(1,numel(r),numel(v),numel(t));
TPR_bS = zeros(1,numel(r),numel(v),numel(t));
TPR_bT = zeros(1,numel(r),numel(v),numel(t));

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
        
        if ReML
            fprintf('t = %d: \n ', t(h));
        else
            fprintf('t = %d, ', t(h));
        end;
        
        %%% Step 3a: estimate model parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % for each voxel proportion
        for i = 1:numel(r)
            
            % for each session
            for j = 1:S
        
                % for each simulation
                for s = 1:N
                
                    % get data and whitening matrix
                    Y = Sim(g,h).data(j,i).Y(:,:,s);
                    W = Sim(g,h).base.W;

                    % estimate using trial-wise design matrix
                    WY = W*Y;
                    WX = W*Sim(g,h).des(j).DM.X_t;
                    bA = (WX'*WX)^(-1) * WX'*WY;
                    Sim(g,h).est(j,i).BA(:,:,s) = bA;

                    % estimate using trial-based design matrices
                    bS = zeros(t(h),v(g));
                    for k = 1:t(h)
                        WX = W*Sim(g,h).des(j).DM.X_s{k};
                        bk = (WX'*WX)^(-1) * WX'*WY;
                        bS(k,:) = bk(1,:);
                    end;
                    Sim(g,h).est(j,i).BS(:,:,s) = bS;
                    
                end;
                
            end;
            
        end;
        
        %%% Step 3b: restricted maximum likelihood %%%%%%%%%%%%%%%%%%%%%%%%
        
        % for each voxel proportion
        for i = 1:numel(r)
            
            % for each session
            for j = 1:S

                % if ReML specified
                if ReML

                    % prepare ReML analysis
                    BA = reshape(Sim(g,h).est(j,i).BA,[t(h) v(g)*N]);
                    YY = (1/(v(g)*N)) * (BA*BA');
                    X  = Sim(g,h).des(j).DM.T;
                    Q{1} = eye(t(h));
                    Q{2} = Sim(g,h).des(j).DM.U;

                    % perform ReML analysis
                    [V, s2] = spm_reml(YY, X, Q);
                    Sim(g,h).est(j,i).Sg = V;
                    Sim(g,h).est(j,i).s2 = full(s2)';
                    % Sim(g,h).des(j).DM.Sg = V;
                    % Sim(g,h).des(j).DM.s2 = full(s2)';

                else

                    % no ReML analysis
                    Sim(g,h).est(j,i).Sg = Sim(g,h).des(j).DM.U;
                    Sim(g,h).est(j,i).s2 = 1;
                    % Sim(g,h).des(j).DM.Sg = Sim(g,h).des(j).DM.U;
                    % Sim(g,h).des(j).DM.s2 = 1;

                end;
                
            end;
            
        end;
        
        %%% Step 3c: decode/classify trials %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % for each voxel proportion
        for i = 1:numel(r)
            
            % for each session
            for j = 1:S
                
                % get design matrices
                x_train = Sim(g,h).des([1:S]~=j).tt;
                x_test  = Sim(g,h).des(j).tt;
                X_train = Sim(g,h).des([1:S]~=j).DM.T;
                X_test  = Sim(g,h).des(j).DM.T;
                
                % for each simulation
                for s = 1:N
                    
                    % specify SVM options
                    opt = sprintf('-s 0 -t 0 -c %s -q', num2str(C));

                    % SVC using trial-wise parameter estimates
                    Y_train = Sim(g,h).est([1:S]~=j,i).BA(:,:,s);
                    Y_test  = Sim(g,h).est(j,i).BA(:,:,s);
                    svm_tr  = svmtrain(x_train, Y_train, opt);
                    x_pred  = svmpredict(x_test, Y_test, svm_tr, '-q');
                    a = sum(x_test==x_pred)/t(h);
                    Sim(g,h).pred.aA(s,i,j) = a;

                    % SVC using trial-based parameter values
                    Y_train = Sim(g,h).est([1:S]~=j,i).BS(:,:,s);
                    Y_test  = Sim(g,h).est(j,i).BS(:,:,s);
                    svm_tr  = svmtrain(x_train, Y_train, opt);
                    x_pred  = svmpredict(x_test, Y_test, svm_tr, '-q');
                    a = sum(x_test==x_pred)/t(h);
                    Sim(g,h).pred.aS(s,i,j) = a;

                    % classification using SL-based ITEM
                    Y_train = [Sim(g,h).est([1:S]~=j,i).BA(:,:,s), ones(t(h),1)];
                    Y_test  = [Sim(g,h).est(j,i).BA(:,:,s), ones(t(h),1)];
                    V_train = Sim(g,h).est([1:S]~=j,i).Sg;
                    V_test  = Sim(g,h).est(j,i).Sg;
                    P_train = inv(V_train);
                    W_test  = sqrtm(inv(V_test));
                    b_train = (Y_train'*P_train*Y_train)^(-1) * Y_train'*P_train*X_train;
                    X_rec   = W_test * Y_test * b_train;
                    X_pred  = [(X_rec(:,1)>X_rec(:,2)), (X_rec(:,2)>X_rec(:,1))];
                    a = sum(X_test(:,1)==X_pred(:,1))/t(h);
                    Sim(g,h).pred.aT(s,i,j) = a;
                    
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
BA_11 = Sim(1,1).est(1,1).BA(:,:,1);

% for each SL size
for g = 1:numel(v)
    
    fprintf('v = %d: ', v(g));
    
    % for each sample size
    for h = 1:numel(t)
        
        fprintf('t = %d, ', t(h));
        
        % true positive rates
        nT = S*t(h);
        for i = 1:numel(r)
            cA = round(sum(Sim(g,h).pred.aA(:,i,:)*t(h),3));
            cS = round(sum(Sim(g,h).pred.aS(:,i,:)*t(h),3));
            cT = round(sum(Sim(g,h).pred.aT(:,i,:)*t(h),3));
            [phat, pciA] = binofit(cA, nT, alpha);
            [phat, pciS] = binofit(cS, nT, alpha);
            [phat, pciT] = binofit(cT, nT, alpha);
            Res(g,h).TPR(1,i) = mean(pciA(:,1)>chlvl);
            Res(g,h).TPR(2,i) = mean(pciS(:,1)>chlvl);
            Res(g,h).TPR(3,i) = mean(pciT(:,1)>chlvl);
        end;
        clear nT cA cS cT phat pci*
        
        % decoding accuracies
        Res(g,h).DA(:,:,1) = mean(Sim(g,h).pred.aA,3);
        Res(g,h).DA(:,:,2) = mean(Sim(g,h).pred.aS,3);
        Res(g,h).DA(:,:,3) = mean(Sim(g,h).pred.aT,3);
        
        % remove fields
        Sim(g,h).data = rmfield(Sim(g,h).data,{'B','Y'});
        Sim(g,h).est  = rmfield(Sim(g,h).est,{'BA','BS'});
        
    end;
    
end;

save('Simulation_B.mat', 'Sim', 'Res', 'total_time');

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
imagesc(BA_11); axis off;
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
        plot(r, Res(g,h).TPR(3,:), '-g');
        axis([(0-0.01), (1+0.01), (0.5-0.01), (1+0.01)]);
        set(gca,'Box','On');
        if g == numel(v) && h == numel(t)
            legend({'LS-A' 'LS-S' 'ITEM'}, 'Location', 'SouthEast');
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
        plot(r, mean(Res(g,h).DA(:,:,1))+std(Res(g,h).DA(:,:,1)), '--r');
        plot(r, mean(Res(g,h).DA(:,:,1))-std(Res(g,h).DA(:,:,1)), '--r');
        plot(r, mean(Res(g,h).DA(:,:,2))+std(Res(g,h).DA(:,:,2)), '--b');
        plot(r, mean(Res(g,h).DA(:,:,2))-std(Res(g,h).DA(:,:,2)), '--b');
        plot(r, mean(Res(g,h).DA(:,:,3))+std(Res(g,h).DA(:,:,3)), '--g');
        plot(r, mean(Res(g,h).DA(:,:,3))-std(Res(g,h).DA(:,:,3)), '--g');
        axis([(0-0.01), (1+0.01), (0.5-0.01), (1+0.01)]);
        set(gca,'Box','On');
        if g == numel(v) && h == numel(t)
            legend({'LS-A' 'LS-S' 'ITEM'}, 'Location', 'SouthEast');
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