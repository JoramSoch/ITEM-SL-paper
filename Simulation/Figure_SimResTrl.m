% ITEM-SL Paper, Figure 4

clear
close all

% specify simulation parameters
t   = [100, 200, 400];          % number of trials per session
v   = [16, 32, 64];             % number of voxels per trial
r   = [0:0.05:1];               % proportion of voxels with information
lab = {'LS-A', 'LS-S', 'ITEM', 'GLMsingle'};
col =  'rbgm';

% load decoding accuracies
S1  = load('Simulation_B.mat');
S2  = load('Simulation_B_GLMsingle.mat');
Res = S2.Res;
for g = 1:numel(v)
    for h = 1:numel(t)
        Res(g,h).DA(:,:,1:3) = S1.Res(g,h).DA;
    end;
end;

% display decoding accuracies
figure('Name', 'decoding accuracies', 'Color', [1 1 1], 'Position', [50 50 1200 900]);

for g = 1:numel(v)
    for h = 1:numel(t)
        subplot(numel(v), numel(t), (g-1)*numel(t)+h);
        hold on;
        for k = 1:numel(lab)
            plot([-1, -1], [0, 1], strcat('-',col(k)), 'LineWidth', 2);
        end;
        for k = [1,2,4,3]
            plot(r, mean(Res(g,h).DA(:,:,k))+std(Res(g,h).DA(:,:,k)), strcat('--',col(k)));
            plot(r, mean(Res(g,h).DA(:,:,k))-std(Res(g,h).DA(:,:,k)), strcat('--',col(k)));
        end;
        for k = [1,2,4,3]
            plot(r, mean(Res(g,h).DA(:,:,k)), strcat('-',col(k)), 'LineWidth', 2);
        end;
        if v(g) == 32 && t(h) == 100
            plot([1/5,1/5], [-2,+2], '-k', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
        end;
        axis([(0-0.01), (1+0.01), (0.5-0.01), (1+0.01)]);
        set(gca,'Box','On');
        if g == numel(v) && h == numel(t)
            legend(lab, 'Location', 'SouthEast');
        end;
        if g == numel(v) && h == 2
            xlabel('proportion of voxels with signal carrying information about experimental conditions', 'FontSize', 14);
        end;
        if h == 1
            yl  = ylabel('decoding accuracy', 'FontSize', 14);
            pos = get(yl,'Position');
            text(pos(1)-1/6, pos(2), sprintf('v = %d', v(g)), 'Rotation', 90, 'FontSize', 16, ...
                 'FontWeight', 'bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
            clear yl pos
        end;
        if g == 1
            title(sprintf('t = %d', t(h)), 'FontSize', 16);
        end;
    end;
end;

% output decoding accuracies
fprintf('\n-> Decoding accuracies:\n');

for g = 1:numel(v)
	fprintf('   - v = %d:\n', v(g));
    for h = 1:numel(t)
		fprintf('     - t = %d:\n', t(h));
		for k = 1:numel(lab)
			fprintf('       - %s: ', lab{k});
            y = mean(Res(g,h).DA(:,:,k),1);
            for l = 1:numel(y)
                fprintf('%0.4f, ', y(l));
            end;
            fprintf('end.\n');
		end;
    end;
end;
fprintf('\n');