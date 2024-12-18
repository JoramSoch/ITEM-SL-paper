% ITEM-SL Paper, Figure 4

clear
close all

% load decoding accuracies
load Simulation_B.mat
t   = [100, 200, 400];          % number of trials per session
v   = [16, 32, 64];             % number of voxels per trial
r   = [0:0.05:1];               % proportion of voxels with information
lab = {'LS-A', 'LS-S', 'ITEM'};

% display decoding accuracies
figure('Name', 'decoding accuracies', 'Color', [1 1 1], 'Position', [50 50 860 860]);

for g = 1:numel(v)
    for h = 1:numel(t)
        subplot(numel(v), numel(t), (g-1)*numel(t)+h);
        hold on;
        plot(r, mean(Res(g,h).DA(:,:,1)), '-r', 'LineWidth', 2);
        plot(r, mean(Res(g,h).DA(:,:,2)), '-b', 'LineWidth', 2);
        plot(r, mean(Res(g,h).DA(:,:,3)), '-g', 'LineWidth', 2);
        plot(r, mean(Res(g,h).DA(:,:,1))+std(Res(g,h).DA(:,:,1)), '--r');
        plot(r, mean(Res(g,h).DA(:,:,1))-std(Res(g,h).DA(:,:,1)), '--r');
        plot(r, mean(Res(g,h).DA(:,:,2))+std(Res(g,h).DA(:,:,2)), '--b');
        plot(r, mean(Res(g,h).DA(:,:,2))-std(Res(g,h).DA(:,:,2)), '--b');
        plot(r, mean(Res(g,h).DA(:,:,3))+std(Res(g,h).DA(:,:,3)), '--g');
        plot(r, mean(Res(g,h).DA(:,:,3))-std(Res(g,h).DA(:,:,3)), '--g');
        if v(g) == 32 && t(h) == 100
            plot([1/5,1/5], [-2,+2], '-k', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2);
        end;
        axis([(0-0.01), (1+0.01), (0.5-0.01), (1+0.01)]);
        set(gca,'Box','On');
        if g == numel(v) && h == numel(t)
            legend({'LS-A', 'LS-S', 'ITEM'}, 'Location', 'SouthEast');
        end;
        if g == numel(v) && h == 2
            xlabel('proportion of voxels with signal carrying information about experimental conditions', 'FontSize', 14);
        end;
        if h == 1
            yl  = ylabel('decoding accuracy', 'FontSize', 14);
            pos = get(yl,'Position');
            text(pos(1)-1/5, pos(2), sprintf('v = %d', v(g)), 'Rotation', 90, 'FontSize', 16, ...
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
		for k = 1:3
			fprintf('       - %s: ', lab{k});
            y = Res(g,h).DA(:,:,k);
            for l = 1:numel(y)
                fprintf('%0.4f, ', y(l));
            end;
            fprintf('end.\n');
		end;
    end;
end;
fprintf('\n');