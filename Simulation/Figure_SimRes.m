% ITEM-SL Paper, Figure 3

clear
close all

% specify simulation parameters
r   =  1/5;                     % proportion of voxels with information
s2n = [0.8, 1.6, 3.2];          % between-scan variance
ISI = [0,4; 2,6; 4,8];          % inter-stimulus-intervals
lab = {'LS-A', 'LS-S', 'ITEM', 'GLMsingle'};

% load decoding accuracies
S1  = load('Simulation_A.mat');
S2  = load('Simulation_A_GLMsingle.mat');
Res = S2.Res;
for g = 1:numel(s2n)
    for h = 1:size(ISI,1)
        Res(g,h).DA(1:3,:) = S1.Res(g,h).DA;
    end;
end;

% display decoding accuracies
figure('Name', 'decoding accuracies', 'Color', [1 1 1], 'Position', [50 50 1200 900]);
cols = [1,0,0; 0,0,1; 0,1,0; 1,0,1];

for g = 1:numel(s2n)
    for h = 1:size(ISI,1)
        subplot(numel(s2n), size(ISI,1), (g-1)*size(ISI,1)+h); hold on;
        Y = cell(1,numel(lab));
        for k = 1:numel(lab), Y{k} = Res(g,h).DA(k,:)'; end;
        violin_plot(Y, [0:0.01:1], 0.01, [1:numel(lab)], cols, 1, 1, 4);
        if r > 0
            axis([(1-0.5), (numel(lab)+0.5), 0.4, 1]);
        else
            axis([(1-0.5), (numel(lab)+0.5), 0.2, 0.8]);
        end;
        set(gca,'Box','On');
        set(gca,'XTick',[1:numel(lab)],'XTickLabel',lab);
        if g == numel(s2n)
            xlabel('analysis approach', 'FontSize', 14);
        end;
        if h == 1
            yl  = ylabel('decoding accuracy', 'FontSize', 14);
            pos = get(yl,'Position');
            text(pos(1)-2/3, pos(2), ['\sigma^2 =', sprintf(' %1.1f', s2n(g))], 'Rotation', 90, 'FontSize', 16, ...
                 'FontWeight', 'bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
            clear yl pos
        end;
        if g == 1
            title(['t_i \sim', sprintf(' U(%d,%d)', ISI(h,1), ISI(h,2))], 'FontSize', 16);
        end;
    end;
end;

% output decoding accuracies
fprintf('\n-> Decoding accuracies:\n');

for g = 1:numel(s2n)
    fprintf('   - sigma^2 = %1.1f:\n', s2n(g));
    for h = 1:size(ISI,1)
        fprintf('     - t_i ~ U(%d,%d): ', ISI(h,1), ISI(h,2));
        for k = 1:numel(lab)
            fprintf('%s: %0.3f; ', lab{k}, median(Res(g,h).DA(k,:)));
        end;
        fprintf('end.\n');
    end;
end;
fprintf('\n');