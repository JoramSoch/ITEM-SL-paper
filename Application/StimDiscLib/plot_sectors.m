function plot_sectors(tc, fs, fw, cf);
% _
% Writes sector indices into the sector fields


% set default parameters
if nargin < 1 || isempty(tc), tc = 'k';    end;
if nargin < 2 || isempty(fs), fs = 18;     end;
if nargin < 3 || isempty(fw), fw = 'bold'; end;
if nargin < 4 || isempty(cf), cf = 5/4;    end;

% get radii/directions
[r, th, ns] = get_rad_dir;
clear r th

% calculate coordinates
[xs, ys] = calc_sect_cent;

% plot sectors
for k = 1:ns
    text(xs(k), ys(k), int2str(k), 'Color', tc, 'FontSize', fs, 'FontWeight', fw, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
end;