function [r, th, ns] = get_rad_dir;
% _
% Returns radii of 4 rings and directions of 12 segments

% radii and directions
r  = [1/8:1/4:7/8];             % radius, in units of 1
th = [105:30:345, 15:30:75];    % direction, in degrees
ns = numel(r)*numel(th);        % number of sectors