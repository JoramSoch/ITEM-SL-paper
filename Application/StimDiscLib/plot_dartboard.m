function plot_dartboard(ls, lw, dx);
% _
% Draws an empty dartboard into the current axis


% set default parameters
if nargin < 1 || isempty(ls), ls = '-k'; end;
if nargin < 2 || isempty(lw), lw = 2.5;  end;
if nargin < 3 || isempty(dx), dx = 0.01; end;

% get radii/directions
[r, th, ns] = get_rad_dir;
clear ns

% calculate coordinates
[xr, yr] = calc_rad(dx);
[xt, yt] = calc_dir;

% plot dartboard
hold on;
for i = 1:numel(r)
    ri = r(i)+1/8;
    xi = xr>=-ri & xr<=+ri;
    plot(xr(xi), yr(i,xi,1), ls, 'LineWidth', lw);
    plot(xr(xi), yr(i,xi,2), ls, 'LineWidth', lw);
end;
for j = 1:numel(th)
    plot(xt(j,:), yt(j,:), ls, 'LineWidth', lw);
end;
axis([min(xr)-dx, max(xr)+dx, min(xr)-dx, max(xr)+dx]);
axis equal off;