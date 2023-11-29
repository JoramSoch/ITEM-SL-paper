function [xb, yb] = calc_sect_area(dx);
% _
% Returns (x,y)-coordinates of sector areas

% get radii/directions
[r, th, ns] = get_rad_dir;

% calculate areas
xc = zeros(ns,4);               % corners
yc = zeros(ns,4);
xb = cell(ns,1);                % boundaries
yb = cell(ns,1);
for k = 1:ns
    rk = r(ceil(k/12));
    tk = th(mod(k,12) + (mod(k,12)==0)*12)*(pi/180);
    xc(k,:) = [cos(tk-15*(pi/180)) * (rk-1/8), cos(tk+15*(pi/180)) * (rk-1/8), cos(tk+15*(pi/180)) * (rk+1/8), cos(tk-15*(pi/180)) * (rk+1/8)];
  % yc(k,:) = [sin(tk-15*(pi/180)) * (rk-1/8), sin(tk+15*(pi/180)) * (rk-1/8), sin(tk+15*(pi/180)) * (rk+1/8), sin(tk-15*(pi/180)) * (rk+1/8)];
    xk1 = [xc(k,1):(sign(xc(k,2)-xc(k,1))*dx):xc(k,2), xc(k,2)];
    yk1 = ((mod(k,12)<=3 | mod(k,12)>9) - (mod(k,12)>3 & mod(k,12)<=9)) * sqrt((rk-1/8)^2 - xk1.^2);
    xk2 = [xc(k,3):(sign(xc(k,4)-xc(k,3))*dx):xc(k,4), xc(k,4)];
    yk2 = ((mod(k,12)<=3 | mod(k,12)>9) - (mod(k,12)>3 & mod(k,12)<=9)) * sqrt((rk+1/8)^2 - xk2.^2);
    xb{k} = [xk1, xk2];
    yb{k} = [yk1, yk2];
end;
clear rk tk xk1 yk1 xk2 yk2