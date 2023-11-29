% Project directories script
% _
% This script generates project directories for the current study.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 10/11/2017, 12:50; 28/02/2023, 11:10; 31/10/2023, 13:44


% set directories
stud_dir = 'C:\Joram\projects\BCCN\VisRec\Soch_et_al_2020\';
bids_dir = strcat(stud_dir,'data/');
stat_dir = strcat(stud_dir,'stats/');
tool_dir = strcat(pwd,'/');

% save directories
save('project_directories.mat', 'stud_dir', 'bids_dir', 'stat_dir', 'tool_dir');