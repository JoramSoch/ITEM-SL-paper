% ITEM-SL Paper, Figure 5

clear
close all

% load data
data     = load('Design_Matrices.mat');
num_mods = numel(data.X);
num_sect = numel(data.pmod.name);

% load models
X = cell(1,num_mods);
for j = 1:num_mods
    if j ~= 4
        X{j} = data.X{j};
    else
       [X{j}, L] = ITEM_get_des_mat(data.names, data.onsets, data.durations, data.pmod, data.orth, data.R, data.settings);
    end;
    X{j} = X{j}(:,1:(1+num_sect));
end;

% analyze models
r_ons  = zeros(num_mods,num_sect);
r_pmod = zeros(num_mods,num_sect);
r_mvar = zeros(num_mods,num_sect);
for j = 1:num_mods
    for k = 1:num_sect
        r_ons(j,k)  = corr(X{j}(:,1), X{j}(:,1+k));
        r_pmod(j,k) = mean(corr(X{j}(:,1+k), X{j}(:,[true, [1:num_sect]~=k])));
        r_mvar(j,k) = corr(X{j}(:,1+k), X{4}(:,1+k));
    end;
end;

% process design matrix
XD = X{4};
XD = XD./repmat(max(abs(XD)),[size(XD,1) 1]);
is = [3:3:num_sect];

% display correlations
figure('Name', 'parametric modulators', 'Color', [1 1 1], 'Position', [10 50 1600 800]);
titles = {'SPM: not orthogonalized', 'SPM: sequentially orthogonalized', ...
          'SPM: not orthogonalized, but mean-centered', 'ITEM\_get\_des\_mat: mean-centered'};
labels = {'correlation with onset regressor', 'correlation with other parametric modulators', 'correlation with underlying modulator variable'};

% for all models
for j = 1:num_mods
    subplot(2,2,j);
    h = bar([r_ons(j,is)', r_pmod(j,is)', r_mvar(j,is)'], 1.2, 'grouped');
    axis([(1-1), (numel(is)+1), (0-0.05), (1+0.15)]);
    if j == 2, legend(labels, 'Location', 'NorthEast'); end;
    set(gca,'FontSize',12);
    set(gca,'XTick',[1:numel(is)],'XTickLabel',is);
    xlabel('stimulus sector', 'FontSize', 16);
    ylabel('correlation coefficient', 'FontSize', 16);
    title(titles{j}, 'FontSize', 20);
end;

% display correlations
figure('Name', 'design matrix', 'Color', [1 1 1], 'Position', [10 50 600 800]);

imagesc(XD);
colormap gray;
set(gca,'FontSize',12);
set(gca,'XTick',[1, 1+is]);
xlabel('regressor index', 'FontSize', 16);
ylabel('fMRI scan', 'FontSize', 16);
title('SPM/ITEM: design matrix', 'FontSize', 20);