function [xt, yt] = calc_dir;
% _
% Returns (x,y)-coordinates of directional boundaries

% get directions
[r, th, ns] = get_rad_dir;
clear r ns

% calculate boundaries
xt = zeros(numel(th),2);
yt = zeros(numel(th),2);
for j = 1:numel(th)
    tj = (th(j)+15)*(pi/180);
    xt(j,2) = cos(tj);
    yt(j,2) = sin(tj);
end;
clear tj