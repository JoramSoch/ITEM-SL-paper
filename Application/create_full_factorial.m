function create_full_factorial(SL_rad, ana_name, img_name)
% _
% Create full factorial batch for given searchlight radius
%     SL_rad   - searchlight radius, in mm
%     ana_name - name of the SL-based ITEM analysis
%     img_name - image filename ('avgCC', 'cvCC')
% 
% written by Joram Soch <joram.soch@bccn-berlin.de>, 07/12/2021, 11:21;
% adapted: 09/12/2021, 13:18; edited: 24/03/2023, 12:18;
% finalized: 21/11/2023, 15:41


% create sector matrix
sects = [ 1:12;
         13:24;
         25:36;
         37:48];
[n,m] = size(sects);
     
% specify basic information
subj_file = 'subjects_all.mat';
MS_file   = 'model_spaces/glms-item.mat';
stat_suff = strcat('ITEM_sects-all_SL-',num2str(SL_rad),'mm_',img_name);

% specify sector maps
img_paths = cell(size(sects));
for i = 1:n
    for j = 1:m
        img_num = sects(i,j)+1;
        img_str = MF_int2str0(img_num,4);
        img_paths{i,j} = strcat('sub-*/mods/','glms-+/','glm-#/',...
                                'ITEM_dec_recon/ITEM_',ana_name,'_SL-',num2str(SL_rad),'mm/',...
                                img_name,'_',img_str,'.nii');
    end;
end;

% specify contrast vectors
con_info(1).type  = 'F';
con_info(1).name  = 'L<>R';
con_info(1).vec   = repmat([+1*ones(1,6), -1*ones(1,6)],[1 4]);
con_info(2).type  = 't';
con_info(2).name  = 'L>R';
con_info(2).vec   = repmat([+1*ones(1,6), -1*ones(1,6)],[1 4]);
con_info(3).type  = 't';
con_info(3).name  = 'R>L';
con_info(3).vec   = repmat([-1*ones(1,6), +1*ones(1,6)],[1 4]);
con_info(4).type  = 'F';
con_info(4).name  = 'T<>B';
con_info(4).vec   = repmat([+1*ones(1,3), -1*ones(1,6), +1*ones(1,3)],[1 4]);
con_info(5).type  = 't';
con_info(5).name  = 'T>B';
con_info(5).vec   = repmat([+1*ones(1,3), -1*ones(1,6), +1*ones(1,3)],[1 4]);
con_info(6).type  = 't';
con_info(6).name  = 'B>T';
con_info(6).vec   = repmat([-1*ones(1,3), +1*ones(1,6), -1*ones(1,3)],[1 4]);
con_info(7).type  = 't';
con_info(7).name  = 'L';
con_info(7).vec   = repmat([+1*ones(1,6), zeros(1,6)],[1 4]);
con_info(8).type  = 't';
con_info(8).name  = 'R';
con_info(8).vec   = repmat([zeros(1,6), +1*ones(1,6)],[1 4]);
con_info(9).type  = 't';
con_info(9).name  = 'T';
con_info(9).vec   = repmat([+1*ones(1,3), zeros(1,6), +1*ones(1,3)],[1 4]);
con_info(10).type = 't';
con_info(10).name = 'B';
con_info(10).vec  = repmat([zeros(1,3), +1*ones(1,6), zeros(1,3)],[1 4]);
for i = 1:n
    con_info(10+i).type = 't';
    con_info(10+i).name = strcat('E',num2str(i));
    con_info(10+i).vec  = kron([zeros(1,i-1), 1, zeros(1,n-i)],ones(1,m));
end;
for j = 1:m
    con_info(10+n+j).type = 't';
    con_info(10+n+j).name = strcat('A',num2str(j));
    con_info(10+n+j).vec  = kron(ones(1,n),[zeros(1,j-1), 1, zeros(1,m-j)]);
end;

% create ITEM analyses
fact_des_full_fact(subj_file, MS_file, img_paths, {}, con_info, stat_suff, true);