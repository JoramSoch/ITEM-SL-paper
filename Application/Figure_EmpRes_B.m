% ITEM-SL Paper, Figure 4B

clear
close all

% add StimDiscLib code
addpath('StimDiscLib/');

% L/R: left/right
figure('Name', 'ITEM: L/R', 'Color', [1 1 1], 'Position', [50 50 500 500]);
data = repmat([1*ones(1,6), 2*ones(1,6)],[1 4]);
cmap = [1, 0, 0; 0, 0, 1];
plot_dartboard;
plot_dataset_ITEM(data, cmap);

% T/B: top/bottom
figure('Name', 'ITEM: T/B', 'Color', [1 1 1], 'Position', [150 150 500 500]);
data = repmat([1*ones(1,3), 2*ones(1,6), 1*ones(1,3)],[1 4]);
cmap = [0, 1, 0; 1, 1, 0];
plot_dartboard;
plot_dataset_ITEM(data, cmap);

% E: eccentricity
figure('Name', 'ITEM: E', 'Color', [1 1 1], 'Position', [250 250 500 500]);
data = kron([1:4],ones(1,12));
cmap = [1, 1, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
plot_dartboard;
plot_dataset_ITEM(data, cmap);

% A: angle
figure('Name', 'ITEM: E', 'Color', [1 1 1], 'Position', [350 350 500 500]);
data = kron(ones(1,4),[1:+1:6, 6:-1:1]);
cmap = [1, 1, 0; 1, 0, 0; 1, 0, 1; 0, 1, 0; 0, 1, 1; 0, 0, 1];
plot_dartboard;
plot_dataset_ITEM(data, cmap);