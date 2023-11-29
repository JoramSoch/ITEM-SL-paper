function [xs, ys] = calc_sect_cent(cf);
% _
% Returns (x,y)-coordinates of sector centers

% set correction factor
if nargin < 1 || isempty(cf), cf = 6/5; end;

% get radii/directions
[r, th, ns] = get_rad_dir;

% calculate centers
xs = zeros(1,ns);               % sectors
ys = zeros(1,ns);
for k = 1:ns
    rk = r(ceil(k/12));
    tk = th(mod(k,12) + (mod(k,12)==0)*12)*(pi/180);
    xs(k) = cos(tk) * rk;
    ys(k) = sin(tk) * rk;
    if rk == min(r)
        xs(k) = cf*xs(k);
        ys(k) = cf*ys(k);
    end;
end;
clear rk tk