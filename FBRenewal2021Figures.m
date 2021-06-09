%%
load('F:/Post Saccadic Dynamics Modeling/Data/Simulated Activities/SacDB/UG - Noise & Grating Simulated Separately/figures - withInternalNoise/uni-no_bias-fa25_hit75 - durOffset=0/PerformanceData.mat');
iL = 4;
sen2 = shiftdim(Sensitivities(iL,1,:,:), 2)';
sen2SD = shiftdim(SensitivitiesSTD(iL,1,:,:), 2)';
sen10 = shiftdim(Sensitivities(iL,2,:,:), 2)';
sen10SD = shiftdim(SensitivitiesSTD(iL,2,:,:), 2)';

durs = durs - 57;

figure('color', 'w');
colors = {'r', 'g', 'b', [1 0.5 0.5], [0.5 1 0.5], [0.5 0.5 1], [0.5 0 0], [0 0.5 0], [0 0 0.5]};
subplot(2,2,1); hold on;
for(k = [1:8])
	fill([durs, fliplr(durs)], [sen2(:,k) - sen2SD(:,k); flipud(sen2(:,k) + sen2SD(:,k))], 'k', 'LineStyle', 'none', 'FaceColor', colors{k}, 'FaceAlpha', 0.5);
	plot(durs, sen2(:,k), 'color', colors{k}, 'LineWidth', 2, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('2 cpd');
YLim1 = ylim;

subplot(2,2,2); hold on;
for(k = [1:8])
	fill([durs, fliplr(durs)], [sen10(:,k) - sen10SD(:,k); flipud(sen10(:,k) + sen10SD(:,k))], 'k', 'LineStyle', 'none', 'FaceColor', colors{k}, 'FaceAlpha', 0.5);
	plot(durs, sen10(:,k), 'color', colors{k}, 'LineWidth', 2, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('10 cpd');
YLim2 = ylim;

subplot(2,2,3); hold on;
for(k = [1:8])
	fill([durs, fliplr(durs)], [sen2(:,k) - sen2SD(:,k); flipud(sen2(:,k) + sen2SD(:,k))] ./ sen2(end,k), 'k', 'LineStyle', 'none', 'FaceColor', colors{k}, 'FaceAlpha', 0.5);
	plot(durs, sen2(:,k) ./ sen2(end,k), 'color', colors{k}, 'LineWidth', 2, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('2 cpd (Normalized)');

subplot(2,2,4); hold on; h = [];
for(k = [1:8])
	fill([durs, fliplr(durs)], [sen10(:,k) - sen10SD(:,k); flipud(sen10(:,k) + sen10SD(:,k))] ./ sen10(end,k), 'k', 'LineStyle', 'none', 'FaceColor', colors{k}, 'FaceAlpha', 0.5);
	h(k) = plot(durs, sen10(:,k) ./ sen10(end,k), 'color', colors{k}, 'LineWidth', 2, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('10 cpd (Normalized)');
legend(h, 'location', 'SouthEast');

set(findobj(gcf,'type','axes'), 'xlim', [0 550], 'yscale', 'log', 'fontsize', 20, 'linewidth', 2, 'box', 'off');

%% un-normalized results
figure('color', 'w');  pause(0.1); jf = get(handle(gcf),'javaframe'); jf.setMaximized(1); pause(1);
colors = {'r', 'g', 'b', [1 0.5 0.5], [0.5 1 0.5], [0.5 0.5 1], [0.5 0 0], [0 0.5 0], [0 0 0.5]};
tIdx = durs >= 0;
iSF = 1;
h = [];
axes( 'Position', [0.17+0.33*(iSF-1), 0.35, 0.2620, 0.9/2-0.04], ...
	  'xlim', [-25 525], 'xtick', 0:150:500, 'ylim', [0.7 9], 'ytick', [1 2 4 8], 'YScale', 'log', 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
for(k = [1:8])
	fill([durs(tIdx), fliplr(durs(tIdx))], [sen2(tIdx,k) - sen2SD(tIdx,k); flipud(sen2(tIdx,k) + sen2SD(tIdx,k))], 'k', 'LineStyle', 'none', 'FaceColor', colors{k}, 'FaceAlpha', 0.5);
	h(k) = plot(durs(tIdx), sen2(tIdx,k), 'color', colors{k}, 'LineWidth', 2, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('2 cpd');
legend(h, 'Position', [0.7866 0.3565 0.1005 0.3777]);

iSF = 2;
axes( 'Position', [0.17+0.33*(iSF-1), 0.35, 0.2620, 0.9/2-0.04], ...
	  'xlim', [-25 525], 'xtick', 0:150:500, 'ylim', [0.7 9], 'ytick', [1 2 4 8], 'YScale', 'log', 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
for(k = [1:8])
	fill([durs(tIdx), fliplr(durs(tIdx))], [sen10(tIdx,k) - sen10SD(tIdx,k); flipud(sen10(tIdx,k) + sen10SD(tIdx,k))], 'k', 'LineStyle', 'none', 'FaceColor', colors{k}, 'FaceAlpha', 0.5);
	plot(durs(tIdx), sen10(tIdx,k), 'color', colors{k}, 'LineWidth', 2, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('10 cpd');


%% results for ecc=0,4,8 and t=50,150,500, un-normalized
figure('color', 'w');  pause(0.1); jf = get(handle(gcf),'javaframe'); jf.setMaximized(1); pause(1);
colors = {'r', 'g', 'b', [1 0.5 0.5], [0.5 1 0.5], [0.5 0.5 1], [0.5 0 0], [0 0.5 0], [0 0 0.5]};
tIdx = abs(durs-50) < 5 | abs(durs-150) < 5 | abs(durs-500) < 5;
iSF = 1;
h = [];
axes( 'Position', [0.17+0.33*(iSF-1), 0.35, 0.2620, 0.9/2-0.04], ...
	  'xlim', [-25 525], 'xtick', 0:150:500, 'ylim', [0.7 9], 'ytick', [1 2 4 8], 'YScale', 'log', 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
for(k = [1 3 5])
	h(k) = errorbar(durs(tIdx), sen2(tIdx,k), -sen2SD(tIdx,k), sen2SD(tIdx,k), 'color', colors{k}, 'LineWidth', 3, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('2 cpd');
legend(h([1 3 5]), 'location', 'SouthEast');

iSF = 2;
axes( 'Position', [0.17+0.33*(iSF-1), 0.35, 0.2620, 0.9/2-0.04], ...
	  'xlim', [-25 525], 'xtick', 0:150:500, 'ylim', [0.7 9], 'ytick', [1 2 4 8], 'YScale', 'log', 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
for(k = [1 3 5])
	errorbar(durs(tIdx), sen10(tIdx,k), -sen10SD(tIdx,k), sen10SD(tIdx,k), 'color', colors{k}, 'LineWidth', 3, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('10 cpd');

%% results for ecc=0,4,8 and t=50,150,500, un-normalized, used to compare with empirical results
figure('color', 'w');  pause(0.1); jf = get(handle(gcf),'javaframe'); jf.setMaximized(1); pause(1);
colors = {'r', [1 0.5 0.5], 'g', [0.5 1 0.5], 'b', [0.5 0.5 1], [0.5 0 0], [0 0.5 0], [0 0 0.5]};
tIdx = abs(durs-50) < 5 | abs(durs-500) < 5;
iSF = 1;
h = [];
axes( 'Position', [0.07+0.33*(iSF-1), 0.35, 0.1620, 0.9/2-0.04], ...
	  'xlim', [-50 600], 'xtick', [50 500], 'ylim', [0.5 14], 'ytick', [1 2 4 8], 'YScale', 'log', 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
for(k = [1 3 5])
	h(k) = errorbar(durs(tIdx), sen2(tIdx,k), -sen2SD(tIdx,k), sen2SD(tIdx,k), 'color', colors{k}, 'LineWidth', 3, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('2 cpd');
legend(h([1 3 5]), 'location', 'SouthEast');

iSF = 2;
axes( 'Position', [0.07+0.33*(iSF-1), 0.35, 0.1620, 0.9/2-0.04], ...
	  'xlim', [-50 600], 'xtick', [50 500], 'ylim', [0.5 14], 'ytick', [1 2 4 8], 'YScale', 'log', 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
for(k = [1 3 5])
	errorbar(durs(tIdx), sen10(tIdx,k), -sen10SD(tIdx,k), sen10SD(tIdx,k), 'color', colors{k}, 'LineWidth', 3, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('10 cpd');

iSF = 2;
tIdx = find(abs(durs-50) < 5 | abs(durs-150) < 5 | abs(durs-500) < 5);
for(k = [1 3 5])
	axes( 'Position', [0.74, 0.69-(k-1)*0.1, 0.15, 0.155], ...
	  'xlim', [0 550], 'xtick', [50 150 500], 'ylim', [0.5 1.1], 'ytick', [0.5 1], 'YScale', 'linear', 'xgrid', 'on', 'ygrid', 'on', 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
	errorbar(durs(tIdx), sen10(tIdx,k) ./ sen10(tIdx(end),k), -sen10SD(tIdx,k) ./ sen10(tIdx(end),k), sen10SD(tIdx,k) ./ sen10(tIdx(end),k), 'color', colors{k}, 'LineWidth', 3, 'DisplayName', sprintf('Ecc = %d', k*2-2));
	if(k ~= 5)
		set(gca, 'xticklabel', []);
	end
	if(k == 1)
		title('10 cpd');
	end
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');



%%
%%%%%% Decision proess %%%%%%
t = [0 550];
frBG = cat(1, encoder.ExampleCellsActivitiesOnCondition(encoder.layers(1).name, 0, 0, 0, 'saccadeOff', t, false, true),...
			  encoder.ExampleCellsActivitiesOnCondition(encoder.layers(2).name, 0, 0, 0, 'saccadeOff', t, false, true),...
			  encoder.ExampleCellsActivitiesOnCondition(encoder.layers(3).name, 0, 0, 0, 'saccadeOff', t, false, true),...
			  encoder.ExampleCellsActivitiesOnCondition(encoder.layers(4).name, 0, 0, 0, 'saccadeOff', t, false, true));
fr2 = cat(1, encoder.ExampleCellsActivitiesOnCondition(encoder.layers(1).name, 0, 2, 0.5, 'saccadeOff', t, false, true),...
			 encoder.ExampleCellsActivitiesOnCondition(encoder.layers(2).name, 0, 2, 0.5, 'saccadeOff', t, false, true),...
			 encoder.ExampleCellsActivitiesOnCondition(encoder.layers(3).name, 0, 2, 0.5, 'saccadeOff', t, false, true),...
			 encoder.ExampleCellsActivitiesOnCondition(encoder.layers(4).name, 0, 2, 0.5, 'saccadeOff', t, false, true));
%%
figure('color', 'w'); pause(0.1); jf = get(handle(gcf),'javaframe'); jf.setMaximized(1); pause(1);
colors = {'k', [0.0745    0.6235    1.0000], [1.0000    0.4118    0.1608]};
names = {'Absent', 'Present', 'Noise + 10cpd'};
FontSize = 30;

% average across cells for iTrial
idx = [27 50 84];
for(iTrial = idx)
	iPlot = find(idx == iTrial);
	hAxes(iPlot) = axes( 'Position', [0.07+0.96/3*(mod(iPlot-1,5)), 0.2+0.92/2, 0.96/3-0.03, 0.9/2-0.04], ...
			             'xlim', t, 'xtick', t(1):200:t(2), 'ylim', [0 1.5], 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
	data = cat(3, frBG(:,:,iTrial), fr2(:,:,iTrial));
    h = [];
	for(k = 1:2)
		m = mean(data(:,:,k), 1);
		sem = std(data(:,:,k), [], 1) ./ sqrt(size(data,1));
		h(k) = plot(t(1) : t(2), m, '-', 'color', colors{k}, 'lineWidth', 2, 'DisplayName', names{k});
		% fill([t(1):t(2), t(2):-1:t(1)], [m-sem, fliplr(m+sem)], 'k', 'FaceColor', colors{k}, 'LineStyle', 'none', 'FaceAlpha', 0.5 );
	end
	if(iTrial == idx(end))
		xlabel('Time from saccade off (ms)');
	else
		set(gca, 'XTickLabel', []);
	end
	ylabel('Average response (s^{-1})');
	title(sprintf('#Trial %d (r_%d)', iPlot, iPlot), 'FontSize', FontSize);
	if(iPlot == 1)
		legend(h, 'location', 'NorthEast', 'FontSize', FontSize, 'LineWidth', 2);
	end
	pos = get(gca, 'position');
	set(gca, 'position', [pos(1:2), pos(3)/1.16, pos(4)/2]);
end
%%
for(mu2 = 110:20:170)
    figure('color', 'w', 'position', [681   611   560   275], 'name', ['mu2 = ' num2str(mu2)]);
    colors = {'k', [0.0745    0.6235    1.0000], [1.0000    0.4118    0.1608]};
    mu1 = 90;
    % mu2 = 110;170;
    Sigma = (170 - 90) / (2 * norminv(0.90)); (mu2 - mu1) / (2 * norminv(0.90));
    th = (mu1 + mu2) / 2;
    x = mu1 - Sigma * 3.5 : 0.1 : mu2 + Sigma * 3.5;
    axes( 'Position', [0.1018    0.2218    0.8661    0.7527], ...
          'xlim', x([1 end]), 'xtick', [], 'YTick', [], 'fontsize', 26, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
    plot(x, normpdf((x-mu1)/Sigma), '-', 'color', colors{1}, 'lineWidth', 2, 'DisplayName', 'Absent');
    fill([x(x>th) fliplr(x(x>th))], [zeros(1,sum(x>th)), normpdf((fliplr(x(x>th))-mu1)/Sigma)], 'k', 'FaceColor', colors{1}, 'LineStyle', 'none', 'FaceAlpha', 0.5);
    plot(x, normpdf((x-mu2)/Sigma), '-', 'color', colors{2}, 'lineWidth', 2, 'DisplayName', 'Present');
    fill([x(x>th) fliplr(x(x>th))], [zeros(1,sum(x>th)), normpdf((fliplr(x(x>th))-mu2)/Sigma)], 'k', 'FaceColor', colors{2}, 'LineStyle', 'none', 'FaceAlpha', 0.5);
    ylim(ylim*1.2);
    plot([1 1]*th, ylim, 'k--', 'lineWidth', 2);
    xlabel('Accumulated response (\eta(t))');
    ylabel('Probability');
    xlim([-20 280]);
    saveas(gcf, sprintf('../../Manuscript/FB Renewal 2021/Accumulated Population Mean FR Distribution - mu2=%d.fig', mu2));
    saveas(gcf, sprintf('../../Manuscript/FB Renewal 2021/Accumulated Population Mean FR Distribution - mu2=%d.png', mu2));
end


%% spatial kernel: 2D version
figure('color', 'w'); hold on;
iL = 1;
cellIdx = encoder.layers(iL).idxExampleCells{5}(1);
rC = encoder.layers(iL).sRFParams(cellIdx).centerRadii * 1.5;
rS = encoder.layers(iL).sRFParams(cellIdx).surroundRadii;
sigmaC = rC/sqrt(2);
sigmaS = rS/sqrt(2);
peakC = encoder.layers(iL).sRFParams(cellIdx).centerPeakSensitivities;
peakS = encoder.layers(iL).sRFParams(cellIdx).surroundPeakSensitivities * 5;
x = linspace(-2.5*sigmaS, 2.5*sigmaS, 1000);
plot(x, peakC * exp(-x.^2/rC^2), 'r', 'lineWidth', 3);
plot(x, -peakS * exp(-x.^2/rS^2), 'b', 'lineWidth', 3);
plot(xlim, [0 0], 'k', 'lineWidth', 2);
plot([0 0], [-peakS-peakC*0.1, peakC*1.1], 'k--', 'lineWidth', 2);
plot([1 1] * rC, [0, peakC*exp(-1)], 'r--', 'lineWidth', 2);
plot([0 rC], [1 1]*peakC*exp(-1), 'r--', 'lineWidth', 2);
plot([1 1] * rS, [0, -peakS-peakC*0.1], 'b--', 'lineWidth', 2);
set(gca, 'ylim', [-peakS-peakC*0.2, peakC*1.2], 'fontsize', 20);