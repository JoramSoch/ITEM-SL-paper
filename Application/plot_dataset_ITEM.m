function plot_dataset(data, cmap, cr, ca, cb, ls, lw, dx);
% _
% Creates a colorplot from a one-by-sector data vector


% set default parameters
if nargin < 1 || isempty(data), data = zeros(1,48);            end;
if nargin < 2 || isempty(cmap), cmap = colormap('jet');        end;
if nargin < 3 || isempty(cr),   cr   = size(cmap,1);           end;
if nargin < 4 || isempty(ca),   ca   = [min(data), max(data)]; end;
if nargin < 5 || isempty(cb),   cb   = 'None';                 end;
if nargin < 6 || isempty(ls),   ls   = '-k';                   end;
if nargin < 7 || isempty(lw),   lw   = 2.5;                    end;
if nargin < 8 || isempty(dx),   dx   = 0.01;                   end;

% get radii/directions
[r, th, ns] = get_rad_dir;
clear r th

% calculate coordinates
[xb, yb] = calc_sect_area(dx);

% plot dataset
hold on;
for k = 1:ns
    ck = round(1+(data(k)-ca(1))/(ca(2)-ca(1))*(cr-1));
    if ck < 1,  ck = 1;  end;
    if ck > cr, ck = cr; end;
    fill(xb{k}, yb{k}, cmap(ck,:), 'EdgeColor', ls(end), 'LineWidth', lw);
end;
caxis(ca);
if ~strcmp(cb,'None')
    colorbar('Location', cb);
end;