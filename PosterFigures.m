%%
sbj = 'A014';
encoder = Encoder([], true);
folder = 'Noise & Grating Simulated Separately';
withInternalNoise = true;
encoder.LoadExampleCellsActivities(fullfile('../../Data/Simulated Activities', sbj, folder, [folder '.mat']));
saveFolder = '../../Manuscript/poster figures';
trials = encoder.activityParams.trials;
egTrial = 3;

%%
%%%%%% stimulus + eye trace, retinal input, average eye traces, luminances %%%%%%
figure('color', 'w'); pause(0.1); jf = get(handle(gcf),'javaframe'); jf.setMaximized(1); pause(1);

% stimulus + eye trace
subplot(3, 3, [1 2 4 5]);
[noise, inputX, inputY] = encoder.LoadNoise( fullfile('../../data', trials(1).backgroundImage), trials(1).pixelAngle/60 );
[grating] = encoder.GenerateGrating( 0, 1.5, 2, 0, trials(1).pixelAngle/60 );
stimulus = (noise + grating) / 4 + 0.5;
imshow( stimulus, 'XData', inputX([1 end]), 'YData', inputY([1 end]) ); hold on;
plot(trials(egTrial).x.position(trials(egTrial).saccadeOn+(-100:500))/60,...
	 trials(egTrial).y.position(trials(egTrial).saccadeOn+(-100:500))/60, 'r', 'LineWidth', 4);
set(gca, 'visible', 'on', 'XLim', [-1.5 7], 'YLim', [-2.5 2.5], 'YTick', -2:2, 'FontSize', 16, 'LineWidth', 2, 'XColor', 'k', 'YColor', 'k');
xlabel('Horizontal position (\circ)');
ylabel('Vertical position (\circ)');

% eye traces
subplot(3, 3, 3); hold on; h = [];
offset = [-150 500];
eyeX = []; eyeY = [];
for( iTrial = size(trials,2) : -1 : 1 )
	x = trials(iTrial).x.position( offset(1)+trials(iTrial).saccadeOn : min(end, offset(2)+trials(iTrial).saccadeOn) );
	y = trials(iTrial).y.position( offset(1)+trials(iTrial).saccadeOn : min(end, offset(2)+trials(iTrial).saccadeOn) );
	eyeX( iTrial, 1:size(x,2) ) = x;
	eyeY( iTrial, 1:size(y,2) ) = y;
end
h(1) = plot( offset(1) : offset(2), eyeX(egTrial,:), 'b', 'lineWidth', 2, 'DisplayName', 'Horizontal' );
h(2) = plot( offset(1) : offset(2), eyeY(egTrial,:), 'r', 'lineWidth', 2, 'DisplayName', 'Vertical' );
% mX = mean(eyeX, 1);
% mY = mean(eyeY, 1);
% sdX = std(eyeX, [], 1);
% sdY = std(eyeY, [], 1);
% h(1) = plot( offset(1) : offset(2), mX, 'b', 'lineWidth', 2, 'DisplayName', 'Horizontal' );
% h(2) = plot( offset(1) : offset(2), mY, 'r', 'lineWidth', 2, 'DisplayName', 'Vertical' );
% fill( [offset(1) : offset(2), offset(2): -1 : offset(1)], [mX-sdX, mX(end:-1:1)+sdX(end:-1:1)], 'b', 'LineStyle', 'none', 'FaceAlpha', 0.5 );
% fill( [offset(1) : offset(2), offset(2): -1 : offset(1)], [mY-sdY, mY(end:-1:1)+sdX(end:-1:1)], 'r', 'LineStyle', 'none', 'FaceAlpha', 0.5 );
ylim(ylim);
h(3) = plot( [0 0], get(gca,'ylim'), 'k--', 'lineWidth', 2, 'DisplayName', 'Saccade on' );
xlabel('Time aligned to saccade on (ms)');
ylabel('Eye position (arcmin)');
legend(h, 'location', 'NorthEast');
set( gca, 'XLim', offset, 'lineWidth', 2, 'fontsize', 16, 'XColor', 'k', 'YColor', 'k' );

% luminance
subplot(3, 3, 6); hold on;
egTrial = 3;
lum = zeros(1, diff(offset)+1);
for(k = 1:size(lum,2))
	r = 0.1;	% 0.1 deg
	lum(k) = mean(stimulus((inputX - eyeX(egTrial,k)/60).^2 + (inputY' - eyeY(egTrial,k)/60).^2 <= r^2)) * 255 * 0.0633 + 0.1682;
end
plot( offset(1) : offset(2), lum, 'r', 'lineWidth', 2 );
ylim(ylim);
plot( [0 0], ylim, 'k--', 'lineWidth', 2 );
xlabel('Time aligned to saccade on (ms)');
ylabel('Luminance (cd/m^2');
set( gca, 'XLim', offset, 'lineWidth', 2, 'fontsize', 16, 'XColor', 'k', 'YColor', 'k' );

% retinal input
t = (0 : 7) * 70;
egTrial = 3;
for(k = 1:size(t,2))
	subplot(3, 8, 16+k); hold on;
	iX = find(inputX >= eyeX(egTrial, t(k)-offset(1)+1)/60 - 0.5 & inputX <= eyeX(egTrial, t(k)-offset(1)+1)/60 + 0.5);
	iY = find(inputY >= eyeY(egTrial, t(k)-offset(1)+1)/60 - 0.5 & inputY <= eyeY(egTrial, t(k)-offset(1)+1)/60 + 0.5);
    iX = iX(1 : min(size(iX,2), size(iY,2)));
    iY = iY(1 : min(size(iX,2), size(iY,2)));
	img = stimulus(iY,iX);
	y = 0.5 - size(img,1)/2 : size(img,1)/2 - 0.5;
	x = 0.5 - size(img,2)/2 : size(img,2)/2 - 0.5;
	img(y'.^2+x.^2 > size(img,1)^2/4) = 1;
	imshow(img, 'XData', x([1 end]), 'YData', y([1 end]));
	plot(cosd(0:359)*size(img,2)/2, sind(0:359)*size(img,2)/2, 'k', 'lineWidth', 2);
	drawnow;
	
	img = getframe(gca); img = img.cdata;
	mask = ones(size(img,1), size(img,2));
	y = 0.5 - size(img,1)/2 : size(img,1)/2 - 0.5;
	x = 0.5 - size(img,2)/2 : size(img,2)/2 - 0.5;
	mask(y'.^2+x.^2 > size(img,1)^2/4) = 0;
	imwrite(img, fullfile(saveFolder, sprintf('Retinal image - %02d.png', k)), 'alpha', mask);
end
saveas(gcf, fullfile(saveFolder, 'Retinal input.fig'));
saveas(gcf, fullfile(saveFolder, 'Retinal input.pdf'));
saveas(gcf, fullfile(saveFolder, 'Retinal input.png'));


%%
%%%%%% Radius VS Spacing %%%%%%
encoder.SpatialModel.DisplaySpacing2RadiusFitting();
set(findobj(gcf, 'type', 'axes'), 'FontSize', 26, 'lineWidth', 2);
subplot(2,2,1);
title('P Center');
ylabel('Radius (\circ)');
subplot(2,2,2);
title('M Center');
subplot(2,2,3);
title('P Surround');
xlabel('Spacing (\circ)');
ylabel('Radius (\circ)');
subplot(2,2,4);
title('M Surround');
xlabel('Spacing (\circ)');
for(k = 1:4)
	subplot(2,2,k);
	pos = get(gca, 'position');
	set(gca, 'position', [pos(1:2), pos(3)/1.8 pos(4)/1.2]);
end


%%
%%%%%% Example neuronal activity %%%%%%
t = [-150 550];
frBG = encoder.ExampleCellsActivitiesOnCondition(encoder.layers(iL).name, 0, 0, 0, 'saccadeOff', t, false, true);	% P On cells
fr2 = encoder.ExampleCellsActivitiesOnCondition(encoder.layers(iL).name, 0, 2, 0.5, 'saccadeOff', t, false, true);
fr10 = encoder.ExampleCellsActivitiesOnCondition(encoder.layers(iL).name, 0, 10, 0.5, 'saccadeOff', t, false, true);
noise = encoder.AddInternalNoise(zeros(size(fr2)), sbj);
%%
figure('color', 'w'); pause(0.1); jf = get(handle(gcf),'javaframe'); jf.setMaximized(1); pause(1);
colors = {'k', [0.0745    0.6235    1.0000], [1.0000    0.4118    0.1608]};
names = {'Natural noise', 'Noise + 2cpd', 'Noise + 10cpd'};
FontSize = 26;

% average across trials
subplot(2,2,1); hold on; h = [];
data = cat(1, mean(frBG,1), mean(fr2,1), mean(fr10,1));
for(k = 1:3)
	m = mean(data(k,:,:), 3);
	sem = std(data(k,:,:), [], 3) ./ sqrt(size(data,3));
	h(k) = plot(t(1) : t(2), m, '-', 'color', colors{k}, 'lineWidth', 1, 'DisplayName', names{k});
	fill([t(1):t(2), t(2):-1:t(1)], [m-sem, fliplr(m+sem)], 'k', 'FaceColor', colors{k}, 'LineStyle', 'none', 'FaceAlpha', 0.5 );
end
ylim(ylim);
plot([0 0], ylim, 'k--', 'lineWidth', 2);
plot([1 1]*mean([trials.saccadeOn] - [trials.saccadeOff])./trials(1).sRate*1000, ylim, '--', 'color', [0.5 0.5 0.5], 'lineWidth', 2);
xlabel('Time from sac. off (ms)');
ylabel('Spike rate (s^{-1})');
title('P On | Average Across Trials');
set(gca, 'XLim', t, 'FontSize', FontSize, 'lineWidth', 2, 'XColor', 'k', 'YColor', 'k');
legend(h, 'location', 'NorthEast', 'FontSize', 18, 'LineWidth', 1);
pos = get(gca, 'position');
set(gca, 'position', [pos(1:2), pos(3)/2, pos(4)/1.4]);

% average across cells for egTrial
subplot(2,2,2); hold on; h = [];
data = cat(3, frBG(:,:,egTrial), fr2(:,:,egTrial), fr10(:,:,egTrial));
for(k = 1:3)
	m = mean(data(:,:,k), 1);
	sem = std(data(:,:,k), [], 1) ./ sqrt(size(data,1));
	h(k) = plot(t(1) : t(2), m, '-', 'color', colors{k}, 'lineWidth', 1, 'DisplayName', names{k});
	fill([t(1):t(2), t(2):-1:t(1)], [m-sem, fliplr(m+sem)], 'k', 'FaceColor', colors{k}, 'LineStyle', 'none', 'FaceAlpha', 0.5 );
end
ylim(ylim);
plot([0 0], ylim, 'k--', 'lineWidth', 2);
plot([1 1]*(trials(egTrial).saccadeOn - trials(egTrial).saccadeOff)./trials(egTrial).sRate*1000, ylim, '--', 'color', [0.5 0.5 0.5], 'lineWidth', 2);
xlabel('Time from sac. off (ms)');
ylabel('Spike rate (s^{-1})');
title(sprintf('P On | Average Across Cells for iTrial=%d', egTrial));
set(gca, 'XLim', t, 'FontSize', FontSize, 'lineWidth', 2, 'XColor', 'k', 'YColor', 'k');
legend(h, 'location', 'NorthEast', 'FontSize', 18, 'LineWidth', 1);
pos = get(gca, 'position');
set(gca, 'position', [pos(1:2), pos(3)/2, pos(4)/1.4]);

%% internal noise
subplot(2,2,3); hold on; h = [];
m = mean(noise(:,:,egTrial), 1);
sem = std(noise(:,:,egTrial), [], 1) ./ sqrt(size(noise,1));
h = plot(t(1) : t(2), m, '-', 'color', 'b', 'lineWidth', 1, 'DisplayName', 'Internal noise');
fill([t(1):t(2), t(2):-1:t(1)], [m-sem, fliplr(m+sem)], 'k', 'FaceColor', 'b', 'LineStyle', 'none', 'FaceAlpha', 0.5 );

ylim([0.05 0.12]);
plot([0 0], ylim, 'k--', 'lineWidth', 2);
plot([1 1]*(trials(egTrial).saccadeOn - trials(egTrial).saccadeOff)./trials(egTrial).sRate*1000, ylim, '--', 'color', [0.5 0.5 0.5], 'lineWidth', 2);
xlabel('Time from sac. off (ms)');
ylabel('Spike rate (s^{-1})');
title(sprintf('Internal Noise | Average Across Cells for iTrial=%d', egTrial));
set(gca, 'XLim', t, 'FontSize', FontSize, 'lineWidth', 2, 'XColor', 'k', 'YColor', 'k');
legend(h, 'location', 'NorthEast', 'FontSize', 18, 'LineWidth', 1);
pos = get(gca, 'position');
set(gca, 'position', [pos(1:2), pos(3)/2, pos(4)/1.4]);

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
x = 0:0.1:270;
Sigma = 25;
mu1 = 90;
mu2 = 170;
th = 120;
axes( 'Position', [0.3635    0.0791    0.2984    0.3926], ...
      'xlim', x([1 end]), 'xtick', [], 'YTick', [], 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
plot(x, normpdf((x-mu1)/Sigma), '-', 'color', colors{1}, 'lineWidth', 2, 'DisplayName', names{1});
fill([x(x>th) fliplr(x(x>th))], [zeros(1,sum(x>th)), normpdf((fliplr(x(x>th))-mu1)/Sigma)], 'k', 'FaceColor', colors{1}, 'LineStyle', 'none', 'FaceAlpha', 0.5);
plot(x, normpdf((x-mu2)/Sigma), '-', 'color', colors{2}, 'lineWidth', 2, 'DisplayName', names{2});
fill([x(x>th) fliplr(x(x>th))], [zeros(1,sum(x>th)), normpdf((fliplr(x(x>th))-mu2)/Sigma)], 'k', 'FaceColor', colors{2}, 'LineStyle', 'none', 'FaceAlpha', 0.5);
ylim(ylim*1.2);
plot([1 1]*th, ylim, 'k--', 'lineWidth', 2);
xlabel('Accumulated response (\eta(t))');
ylabel('Probability');


%%
%%%%%% population mean response %%%%%%
encoder.MeanPopulationActivity(true);
h = findobj(gcf, 'type', 'axes');
delete(h([1:4, 6:9]));
h = h([5 10]);
%%
set(h, 'FontSize', 30, 'XColor', 'k', 'YColor', 'k');
h(1).Title.String = '';
subplot(h(1));
ylabel({'2 cpd', 'Spike rate (s^{-1})'});
set(gca, 'position', [0.07    0.5811    0.2620    0.4100]);
h(1).Legend.FontSize = 20;
h(1).Legend.Position = [0.1911    0.7624    0.1198    0.2102];
h(1).Legend.String = {'Ecc=0, present'  'Ecc=4, present'  'Ecc=8, present'  'Ecc=0, absent'  'Sac. off'  'Mean Sac. On'};
subplot(h(2));
xlabel('Time from saccade off (ms)');
ylabel({'10 cpd', 'Spike rate (^{-1})'});
set(gca, 'position', [0.07    0.125    0.2620    0.4100]);
h = get(gca, 'children');

%% compare model with empirical
root = '../../Data';
load(fullfile(root, 'Simulated Activities', sbj, folder, 'figures - withInternalNoise/thresholding-uni - durOffset=57/PerformanceData.mat'));
load(fullfile(root, 'Data/statsTable.mat'));
colors = {'r', 'g', 'b', [1 0.5 0.5], [0.5 1 0.5], [0.5 0.5 1], [0.5 0 0], [0 0.5 0], [0 0 0.5]};
durs = [0 50 150 500];
SFs = [2 10]; nSFs = size(SFs,2);
Eccs = [0 4 8]; nEccs = size(Eccs,2);
for(iSF = 1:nSFs)
	axes( 'Position', [0.07+0.96/3*2, 0.125+0.92/2*(2-iSF), 0.2620, 0.9/2-0.04], ...
		  'xlim', [0 550], 'xtick', durs(2:end), 'ylim', [0.9 45], 'ytick', [1 10], 'YScale', 'log', 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
    h= [];
	ind = cell(1,nEccs);
	modelSen = zeros(nEccs, 4);
	for(iEcc = 1:nEccs)
		empiricalSen = statsTable.logSens{find(statsTable.spatialFreq == SFs(iSF) & statsTable.eccentricity == Eccs(iEcc), 1, 'last')};
		empiricalSenSD = statsTable.logSensSD{find(statsTable.spatialFreq == SFs(iSF) & statsTable.eccentricity == Eccs(iEcc), 1, 'last')};

        ind{iEcc} = statsTable.Subj{find(statsTable.spatialFreq == SFs(iSF) & statsTable.eccentricity == Eccs(iEcc), 1, 'first')};
		for(iSbj = 1:size(ind{iEcc},1))
			if(strcmpi(sbj, ind{iEcc}.subject{iSbj}))
				matchSbjIdx = iSbj;
				h(2) = plot([50 150 500], [ind{iEcc}.logSens50(iSbj), ind{iEcc}.logSens150(iSbj), ind{iEcc}.logSens500(iSbj)], 'd--', 'color', colors{6+iEcc}, 'MarkerSize', 14, 'LineWidth', 2, 'DisplayName', sprintf('Ecc=%d, human individual', Eccs(iEcc)));
			else
				;% plot([50 150 500], [ind{iEcc}.logSens50(iSbj), ind{iEcc}.logSens150(iSbj), ind{iEcc}.logSens500(iSbj)], 'd--', 'color', colors{3+iEcc}, 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Ecc=%d, Empirical', Eccs(iEcc)));
			end
		end

		h(3) = plot([50 150 500]+20, empiricalSen, 'o-', 'color', colors{3+iEcc}, 'MarkerFaceColor', 'none', 'MarkerSize', 14, 'LineWidth', 2, 'DisplayName', sprintf('Ecc=%d, human population', Eccs(iEcc)));
		plot(reshape([1; 1; NaN]*([50 150 500]+20), 1, []), reshape([empiricalSen-empiricalSenSD; empiricalSen+empiricalSenSD; NaN(1,3)], 1, []), '-', 'color', colors{3+iEcc}, 'LineWidth', 2);

		modelSen(iEcc,:) = 1./squeeze(Thresholds(5,iSF,iEcc,:));
		modelSen(iEcc,:) = modelSen(iEcc,:) / modelSen(iEcc,end) * ind{iEcc}.logSens500(matchSbjIdx);
		% modelSen(iEcc,:) = modelSen(iEcc,:) / (1./squeeze(Thresholds(5,iSF,1,end))) * ind{iEcc}.logSens500(matchSbjIdx);
		h(1) = plot(durs(2:end), modelSen(iEcc, 2:end), 's-', 'color', colors{iEcc}, 'MarkerSize', 16, 'LineWidth', 2, 'DisplayName', sprintf('Ecc=%d, model', Eccs(iEcc)));
	end
	ylabel('Sensitivity');
    set(gca, 'YTickLabel', sprintf('%g\n',(get(gca, 'YTick'))));
    if(iSF == 1)
    	legend(h, 'location', 'best');
    	set(gca, 'xticklabel', []);
    else
    	xlabel('Time from saccade off (ms)');
    end
end

%% model
for(iSF = 1:nSFs)
	axes( 'Position', [0.07+0.96/3*1, 0.125+0.92/2*(2-iSF), 0.2620, 0.9/2-0.04], ...
		  'xlim', [0 550], 'xtick', durs(2:end), 'ylim', [0.9 45], 'ytick', [1 10], 'YScale', 'log', 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
    h = [];
	ind = cell(1,nEccs);
	modelSen = zeros(nEccs, 4);
	for(iEcc = 1:nEccs)
		modelSen(iEcc,:) = 1./squeeze(Thresholds(5,iSF,iEcc,:));
		h(iEcc) = plot(durs(2:end), modelSen(iEcc, 2:end), 's-', 'color', colors{iEcc}, 'MarkerSize', 16, 'LineWidth', 2, 'DisplayName', sprintf('Ecc=%d, model', Eccs(iEcc)));
	end
	ylabel('Sensitivity');
    set(gca, 'YTickLabel', sprintf('%g\n',(get(gca, 'YTick'))));
    if(iSF == 1)
    	legend(h, 'location', 'best');
    	set(gca, 'xticklabel', []);
    else
    	xlabel('Time from saccade off (ms)');
    end
end


%%
%%%%%% Stabilization %%%%%%
driftStab = load(fullfile(root, 'Simulated Activities', sbj, 'Noise & Grating Simulated Separately - Drift Stabilized', 'figures - withInternalNoise/thresholding-uni - durOffset=57/PerformanceData.mat'));
sacStab = load(fullfile(root, 'Simulated Activities', sbj, 'Noise & Grating Simulated Separately - Saccade Stabilized (static+static+drift)', 'figures - withInternalNoise/thresholding-uni - durOffset=57/PerformanceData.mat'));
stabThresholds = {driftStab.Thresholds, sacStab.Thresholds};
stabNames = {'drift stabilized', 'saccade stabilized'};
figure('color', 'w'); pause(0.1); jf = get(handle(gcf),'javaframe'); jf.setMaximized(1); pause(1);
colors = {'r', 'g', 'b', [1 0.5 0.5], [0.5 1 0.5], [0.5 0.5 1], [0.5 0 0], [0 0.5 0], [0 0 0.5]};
durs = [0 50 150 500];
SFs = [2 10]; nSFs = size(SFs,2);
Eccs = [0 4 8]; nEccs = size(Eccs,2);

for(iSF = 1:nSFs)
	for(iStab = 1:2)
		axes( 'Position', [0.27+0.96/2.5*(iStab-1), 0.125+0.92/2*(2-iSF), 0.32, 0.9/2-0.04], ...
			  'xlim', [-25 525], 'xtick', durs,  'YScale', 'log', 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
	    h = [];
		ind = cell(1,nEccs);
		modelSen = zeros(nEccs, 4);
		for(iEcc = 1:nEccs)
			modelSen(iEcc,:) = 1./squeeze(Thresholds(5,iSF,iEcc,:));
			modelSenStab(iEcc,:) = 1./squeeze(stabThresholds{iStab}(5,iSF,iEcc,:));
			h(iEcc) = plot(durs(1:end), modelSen(iEcc, 1:end), 's-', 'color', colors{iEcc}, 'MarkerSize', 16, 'LineWidth', 2, 'DisplayName', sprintf('Ecc=%d, unstabilized', Eccs(iEcc)));
			h(iEcc+nEccs) = plot(durs(1:end), modelSenStab(iEcc, 1:end), 'd--', 'color', colors{iEcc}, 'MarkerSize', 16, 'LineWidth', 2, 'DisplayName', sprintf('Ecc=%d, stabilized', Eccs(iEcc)));
		end
		ylabel('Sensitivity');
	    if(iSF == 1)
	    	if(iStab == 1)
		    	legend(h, 'position', [0.0043 0.6518 0.2047 0.2842]);
		    end
	    	set(gca, 'xticklabel', [], 'ylim', [2 22], 'ytick', [3 10]);
	    else
	    	xlabel('Time from saccade off (ms)');
            set(gca, 'ylim', [0.3 25], 'ytick', [1 10]);
	    end
	    set(gca, 'YTickLabel', sprintf('%g\n',(get(gca, 'YTick'))));
	end
end