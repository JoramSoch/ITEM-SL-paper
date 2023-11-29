function [xr, yr] = calc_rad(dx)
% _
% Returns (x,y)-coordinates of radial boundaries

% get radii
[r, th, ns] = get_rad_dir;
clear th ns

% calculate boundaries
xr = [-1:dx:+1];
yr = zeros(numel(r),numel(xr),2);
for i = 1:numel(r)
    ri = r(i)+1/8;
    yr(i,:,1) = +sqrt(ri^2 - xr.^2);
    yr(i,:,2) = -sqrt(ri^2 - xr.^2);
end;
clear ri