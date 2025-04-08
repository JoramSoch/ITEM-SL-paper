function violin_plot(Y, y_bin, p_thr, x, cols, vw, lwd, ms)
% _
% Display violin (and sina) plots


% set missing parameters
col = [1,0,0; 0,1,0; 0,0,1; 0,1,1; 1,0,1; 1,1,0; 0,0,0];
if size(y_bin,1) == 1,           y_bin = y_bin';       end;
if nargin < 3 || isempty(p_thr), p_thr = 0.01;         end;
if nargin < 4 || isempty(x),     x     = [1:numel(Y)]; end;
if nargin < 5 || isempty(cols),  cols  = col;          end;
if nargin < 6 || isempty(vw),    vw    = 1;            end;
if nargin < 7 || isempty(lwd),   lwd   = 1;            end;
if nargin < 8 || isempty(ms),    ms    = 10;           end;
clear col

% extract dimensions
v = numel(Y);
% = numel(x);
n = zeros(1,v);
k = numel(y_bin);
for j = 1:v
    n(j) = numel(Y{j});
end;

% calculate density estimates
Y_ksd = zeros(k,v);
Y_min = zeros(1,v);
Y_max = zeros(1,v);
Y_mean= zeros(1,v);
for j = 1:v
    Y_mean(j)  = mean(Y{j});
    Y_ksd(:,j) = ksdensity(Y{j}, y_bin);
    Y_ksd(:,j) = Y_ksd(:,j)./max(Y_ksd(:,j));
    y_j        = y_bin(y_bin < Y_mean(j) & Y_ksd(:,j) < p_thr);
    if ~isempty(y_j), Y_min(j) = max(y_j);
    else,             Y_min(j) = min(y_bin); end;
    y_j        = y_bin(y_bin > Y_mean(j) & Y_ksd(:,j) < p_thr);
    if ~isempty(y_j), Y_max(j) = min(y_j);
    else,             Y_max(j) = max(y_bin); end;
end;

% plot density estimates
rng(1);
for j = 1:v
    x_j     = Y_ksd(y_bin>=Y_min(j) & y_bin<=Y_max(j),j)*(vw/2);
    y_j     = y_bin(y_bin>=Y_min(j) & y_bin<=Y_max(j));
    plot(x(j)-x_j, y_j, '-', 'Color', cols(j,:), 'LineWidth', lwd);
    plot(x(j)+x_j, y_j, '-', 'Color', cols(j,:), 'LineWidth', lwd);
    for l = 1:n(j)
        [m,ind] = min(abs(y_bin-Y{j}(l)));
        x_jl    = (rand-0.5)*Y_ksd(ind,j)*(vw*0.9);
        plot(x(j)+x_jl, Y{j}(l), '.', 'Color', cols(j,:), 'MarkerSize', ms);
    end;
    [m,ind] = min(abs(y_bin-median(Y{j})));
    x_j     = Y_ksd(ind,j)*(vw/2);
    plot([x(j)-x_j, x(j)+x_j], [median(Y{j}), median(Y{j})], '-k', 'LineWidth', lwd);
end;
clear x_j y_j x_jl