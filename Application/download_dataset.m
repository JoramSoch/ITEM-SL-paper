% OpenNeuro download script (ds002013)
% _
% This script takes the GitHub repository of an OpenNeuro dataset and
% downloads the entire dataset from OpenNeuro.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 09/07/2021, 10:00; 28/02/2023, 09:18; 15/11/2023, 15:30


% set directories
load project_directories.mat
data_dir = bids_dir;
temp_dir = 'C:\Joram\projects\BCCN\VisRec\Soch_et_al_2020\tools';
% This "template directory" is the directory which contains, as a
% sub-folder, a cloned version of the GitHub repository "ds002013"
% (URL: https://github.com/OpenNeuroDatasets/ds002013).

% specify dataset
dataset  = 'ds002013';
snapshot = '1.0.3';

% prepare source
url_base = strcat('https://openneuro.org/crn/datasets/',dataset,'/snapshots/',snapshot,'/files/');
fprintf('\n');

% LEVEL 0
if ~exist(strcat(data_dir),'dir')
    mkdir(strcat(data_dir));
end;
f0 = dir(strcat(temp_dir,'/',dataset,'/'));
for i0 = 1:numel(f0)
    if ~strncmp(f0(i0).name,'.',1)
        % file
        if ~f0(i0).isdir
            [fold, file, ext] = fileparts(strcat(temp_dir,'/',dataset,'/',f0(i0).name));
            URL = strcat(url_base,f0(i0).name);
            filename = strcat(data_dir,'/',f0(i0).name);
            fprintf('-> Download: "%s".\n', f0(i0).name);
            websave(filename, URL);
        % directory
        else
            % LEVEL 1
            if ~exist(strcat(data_dir,'/',f0(i0).name),'dir')
                mkdir(strcat(data_dir,'/',f0(i0).name));
            end;
            f1 = dir(strcat(temp_dir,'/',dataset,'/',f0(i0).name,'/'));
            for i1 = 1:numel(f1)
                if ~strncmp(f1(i1).name,'.',1)
                    % file
                    if ~f1(i1).isdir
                        [fold, file, ext] = fileparts(strcat(temp_dir,'/',dataset,'/',f0(i0).name,'/',f1(i1).name));
                        URL = strcat(url_base,f0(i0).name,':',f1(i1).name);
                        filename = strcat(data_dir,'/',f0(i0).name,'/',f1(i1).name);
                        fprintf('-> Download: "%s/%s".\n', f0(i0).name, f1(i1).name);
                        websave(filename, URL);
                    % directory
                    else
                        % LEVEL 2
                        if ~exist(strcat(data_dir,'/',f0(i0).name,'/',f1(i1).name),'dir')
                            mkdir(strcat(data_dir,'/',f0(i0).name,'/',f1(i1).name));
                        end;
                        f2 = dir(strcat(temp_dir,'/',dataset,'/',f0(i0).name,'/',f1(i1).name,'/'));
                        for i2 = 1:numel(f2)
                            if ~strncmp(f2(i2).name,'.',1)
                                % file
                                if ~f2(i2).isdir
                                    [fold, file, ext] = fileparts(strcat(temp_dir,'/',dataset,'/',f0(i0).name,'/',f1(i1).name,'/',f2(i2).name));
                                    URL = strcat(url_base,f0(i0).name,':',f1(i1).name,':',f2(i2).name);
                                    filename = strcat(data_dir,'/',f0(i0).name,'/',f1(i1).name,'/',f2(i2).name);
                                    fprintf('-> Download: "%s/%s/%s".\n', f0(i0).name, f1(i1).name, f2(i2).name);
                                    websave(filename, URL);
                                end;
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

% finalize download
fprintf('\n');