%%
load('F:/Post Saccadic Dynamics Modeling/Data/Simulated Activities/SacDB/UG - Noise & Grating Simulated Separately/figures - withInternalNoise/uni-no_bias-fa25_hit75 - durOffset=0/PerformanceData.mat');
iL = 5;
sen2 = shiftdim(Sensitivities(iL,1,:,:), 2)';
sen2SD = shiftdim(SensitivitiesSTD(iL,1,:,:), 2)';
sen10 = shiftdim(Sensitivities(iL,2,:,:), 2)';
sen10SD = shiftdim(SensitivitiesSTD(iL,2,:,:), 2)';

% durs = durs+1; - 7;57;
% durs = [0 durs];
% sen2 = [ones(1,size(sen2,2))*eps; sen2];
% sen2SD = [ones(1,size(sen2SD,2))*eps; sen2SD];
% sen10 = [ones(1,size(sen10,2))*eps; sen10];
% sen10SD = [ones(1,size(sen10SD,2))*eps; sen10SD];

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
ylim([0.5 10]);
YLim1 = ylim;

subplot(2,2,2); hold on;
for(k = [1:8])
	fill([durs, fliplr(durs)], [sen10(:,k) - sen10SD(:,k); flipud(sen10(:,k) + sen10SD(:,k))], 'k', 'LineStyle', 'none', 'FaceColor', colors{k}, 'FaceAlpha', 0.5);
	plot(durs, sen10(:,k), 'color', colors{k}, 'LineWidth', 2, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('10 cpd');
ylim([0.5 10]);
YLim2 = ylim;

subplot(2,2,3); hold on;
for(k = [1:8])
	fill([durs, fliplr(durs)], [sen2(:,k) - sen2SD(:,k); flipud(sen2(:,k) + sen2SD(:,k))] ./ sen2(end,k), 'k', 'LineStyle', 'none', 'FaceColor', colors{k}, 'FaceAlpha', 0.5);
	plot(durs, sen2(:,k) ./ sen2(end,k), 'color', colors{k}, 'LineWidth', 2, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('2 cpd (Normalized)');
ylim([0.1 1.3]);

subplot(2,2,4); hold on; h = [];
for(k = [1:8])
	fill([durs, fliplr(durs)], [sen10(:,k) - sen10SD(:,k); flipud(sen10(:,k) + sen10SD(:,k))] ./ sen10(end,k), 'k', 'LineStyle', 'none', 'FaceColor', colors{k}, 'FaceAlpha', 0.5);
	h(k) = plot(durs, sen10(:,k) ./ sen10(end,k), 'color', colors{k}, 'LineWidth', 2, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('10 cpd (Normalized)');
legend(h, 'location', 'SouthEast');
ylim([0.1 1.3]);

set(findobj(gcf,'type','axes'), 'xlim', [0 550], 'yscale', 'log', 'fontsize', 20, 'linewidth', 2, 'box', 'off');

%% un-normalized results
figure('color', 'w');  pause(0.1); jf = get(handle(gcf),'javaframe'); jf.setMaximized(1); pause(1);
colors = {'r', 'g', 'b', [1 0.5 0.5], [0.5 1 0.5], [0.5 0.5 1], [0.5 0 0], [0 0.5 0], [0 0 0.5]};
tIdx = durs >= 0;
iSF = 1;
h = [];
axes( 'Position', [0.17+0.33*(iSF-1), 0.35, 0.2620, 0.9/2-0.04], ...
	  'xlim', [-25 525], 'xtick', 0:150:500, 'ylim', [0.7 10], 'ytick', [1 2 4 8], 'YScale', 'log', 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
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
	  'xlim', [-25 525], 'xtick', 0:150:500, 'ylim', [0.7 10], 'ytick', [1 2 4 8], 'YScale', 'log', 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
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
	  'xlim', [-25 525], 'xtick', 0:150:500, 'ylim', [0.7 10], 'ytick', [1 2 4 8], 'YScale', 'log', 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
for(k = [1 3 5])
	h(k) = errorbar(durs(tIdx), sen2(tIdx,k), -sen2SD(tIdx,k), sen2SD(tIdx,k), 'color', colors{k}, 'LineWidth', 3, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('2 cpd');
legend(h([1 3 5]), 'location', 'SouthEast');

iSF = 2;
axes( 'Position', [0.17+0.33*(iSF-1), 0.35, 0.2620, 0.9/2-0.04], ...
	  'xlim', [-25 525], 'xtick', 0:150:500, 'ylim', [0.7 10], 'ytick', [1 2 4 8], 'YScale', 'log', 'fontsize', 30, 'LineWidth', 2, 'nextplot', 'add', 'XColor', 'k', 'YColor', 'k' );
for(k = [1 3 5])
	errorbar(durs(tIdx), sen10(tIdx,k), -sen10SD(tIdx,k), sen10SD(tIdx,k), 'color', colors{k}, 'LineWidth', 3, 'DisplayName', sprintf('Ecc = %d', k*2-2));
end
ylabel('Sensitivity');
xlabel('Time from saccade off (ms)');
title('10 cpd');


%% results for ecc=0,4,8 and t=50,150,500, un-normalized, used to compare with empirical results
durs = durs - 7;	% compensate the difference btween online & offline saccade detection 
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
durs = durs + 7;


%% Video showing P On cell activities averaged across trials / example trial
paramFolder = './Parameters/';
HsBe1 = true;
encoder = Encoder(paramFolder, HsBe1);
encoder.LoadExampleCellsActivities(fullfile('../../Data/Simulated Activities', 'SacDB', 'UG - Noise & Grating Simulated Separately', ['All Conditions' '.mat']));
encoder.AddInternalNoise([], fullfile('../../Data/Simulated Activities', 'SacDB', 'UG - Noise & Grating Simulated Separately', ['All Conditions' '.mat']));

%%
trials = encoder.activityParams.trials;
nTrials = size(trials,2);
for(iTrial = [273 158 102 281 156])%randperm(nTrials, 40))% [80 283 194 133 120 126])%)
	figure('NumberTitle', 'off', 'name', sprintf('Demo: P On cell Activity Map - iTrial=%d', iTrial), 'color', 'w'); pause(0.1); jf = get(handle(gcf),'javaframe'); jf.setMaximized(1); pause(1);
	iL = 1;
	ecc = 0; 8;
	iEcc = 1; 5;	% ecc = 8

	cellXRange = [min(encoder.layers(iL).locations(encoder.layers(iL).idxExampleCells{iEcc}, 1)), max(encoder.layers(iL).locations(encoder.layers(iL).idxExampleCells{iEcc}, 1))];
	cellYRange = [min(encoder.layers(iL).locations(encoder.layers(iL).idxExampleCells{iEcc}, 2)), max(encoder.layers(iL).locations(encoder.layers(iL).idxExampleCells{iEcc}, 2))];
	cellXRange = [1.1  -0.1; -0.1 1.1] * cellXRange';
	cellYRange = [1.1  -0.1; -0.1 1.1] * cellYRange';
	
	% iTrial = 47;1;
	egTrial = trials(iTrial);
	tTicks = -150:600;	% aligned to saccade off

	% cell response maps
	r = mean( [encoder.layers(iL).sRFParams(encoder.layers(iL).idxExampleCells{iEcc}).centerRadii] ) / 2.2 / 2;
	frBG = encoder.ExampleCellsActivitiesOnCondition(encoder.layers(iL).name, ecc, 0, 0, 'saccadeOff', tTicks([1 end]), false, true);
	fr2 = encoder.ExampleCellsActivitiesOnCondition(encoder.layers(iL).name, ecc, 2, 0.5, 'saccadeOff', tTicks([1 end]), false, true);
	fr10 = encoder.ExampleCellsActivitiesOnCondition(encoder.layers(iL).name, ecc, 10, 0.5, 'saccadeOff', tTicks([1 end]), false, true);
	fr = {fr2, frBG, fr10};
	% fr = frBG(~isnan(frBG(:))); frMaxBG = std(fr)*5;
	% fr = fr2(~isnan(fr2(:))); frMax2 = std(fr)*5;
	% fr = fr10(~isnan(fr10(:))); frMax10 = std(fr)*5;
	% frMax = [frMax2, frMaxBG, frMax10];
	frMax = [max(max(mean(fr2(:,:,iTrial),3),[],1)), max(max(mean(frBG(:,:,iTrial),3),[],1)), max(max(mean(fr10(:,:,iTrial),3),[],1))];
	frMax(:) = max(frMax);
	names = {'2 cpd', 'Absent', '10 cpd'};
	clear hFR;
	for(k = 1 : 3)
		axes('position', [0.4169+(k-2)*0.2691, 0.5346 0.1661 0.3804], 'nextplot', 'add');
		colormap(gca, 'hot');
		colors = colormap(gca);
		iTick = 1;
		eyeX = 0;
		eyeY = 0;
		for( iCell = length(encoder.layers(iL).idxExampleCells{iEcc}) : -1 : 1 )
			hFR(k,iCell) = fill( r*cosd(0:10:360) + encoder.layers(iL).locations(encoder.layers(iL).idxExampleCells{iEcc}(iCell), 1) + eyeX, r*sind(0:10:360) + encoder.layers(iL).locations(encoder.layers(iL).idxExampleCells{iEcc}(iCell), 2) + eyeY,...
							   'k', 'FaceColor', colors(round(mean(fr{k}(iCell, iTick, iTrial), 3) / frMax(k) * 255) + 1, :), 'FaceAlpha', 1, 'LineStyle', 'none' );
		end
		% xlabel('Horizontal position (\circ)');
		% ylabel('Vertical position (\circ)');
		% title(['Normalized Cell Response | ', names{k}]);
		% if(k == 1)
		% 	ylabel('Normalized Cell Response');
		% end
		axis equal;
		set(gca, 'xlim', cellXRange + eyeX, 'ylim', cellYRange, 'xtick', [], 'ytick', [], 'fontsize', 16, 'LineWidth', 2, 'color', 'k', 'XColor', 'k', 'YColor', 'k');
		colorbar;
	end
	% axes('position', [0 0 1 1], 'visible', 'off');
	% hTxt = text(0.5, 0.95, sprintf('Time from Saccade Offset: %d ms', tTicks(iTick)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 24);

	% cell responses as a function of time
	axes('position', [0.1300 0.3169 0.7750 0.2157], 'nextplot', 'add'); h = [];
	colors2 = {[0.0745    0.6235    1.0000], 'k', [1.0000    0.4118    0.1608]};
	for(k = 1:3)
		m = mean(fr{k}(:,:,iTrial), 1);
		sem = std(fr{k}(:,:,iTrial), [], 1) / sqrt(size(fr{k}, 1));
		fill([tTicks, fliplr(tTicks)], [m-sem, fliplr(m+sem)], 'k', 'FaceColor', colors2{k}, 'LineStyle', 'none', 'FaceAlpha', 0.5 );
		h(k) = plot(tTicks, m, '-', 'color', colors2{k}, 'lineWidth', 2, 'displayname', names{k});
	end
	hTime1 = plot([1 1]*tTicks(iTick), ylim, 'k--', 'LineWidth', 2);
	plot([0 0], ylim, '-', 'color', [0.5 0.5 0.5], 'lineWidth', 2);
	legend(h, 'location', 'northeast');
	ylabel('Firing rate (s^{-1})');
	hTitle = title(sprintf('Eye Trace | t = %d ms', tTicks(iTick)));
	set(gca, 'xlim', tTicks([1 end]) + [-10 10], 'ylim', ylim, 'XTickLabel', [], 'lineWidth', 2, 'fontsize', 16, 'XColor', 'k', 'YColor', 'k');

	% eye trace
	axes('position', [0.1300 0.0734 0.7750 0.2157], 'nextplot', 'add'); h = [];
	h(1) = plot(tTicks, egTrial.x.position(egTrial.saccadeOff + round(tTicks/1000*egTrial.sRate)), 'LineWidth', 2, 'displayname', 'Horizontal');
	h(2) = plot(tTicks, egTrial.y.position(egTrial.saccadeOff + round(tTicks/1000*egTrial.sRate)), 'LineWidth', 2, 'displayname', 'Vertical');
	hTime2 = plot([1 1]*tTicks(iTick), ylim, 'k--', 'LineWidth', 2);
	plot([0 0], ylim, '-', 'color', [0.5 0.5 0.5], 'lineWidth', 2);
	xlabel('Time from saccade off (ms)');
	ylabel('Eye position (arcmin)');
	legend(h, 'location', 'east');
	set(gca, 'xlim', tTicks([1 end]) + [-10 10], 'ylim', ylim, 'LineWidth', 2, 'FontSize', 16, 'XColor', 'k', 'YColor', 'k');

% 	continue;

	% generate movie
	saveFolder = '../../Manuscript/FB Renewal 2021/Videos';
	filename = fullfile(saveFolder, sprintf('Demo_Activity - - P On Cells - Fovea - iTrial=%d', iTrial));
	writerObj = VideoWriter(filename, 'MPEG-4');
	open(writerObj);
	for(iTick = 1 : size(tTicks,2))

		% cell response
		for(k = 1:3)
			for(iCell = 1 : length(encoder.layers(iL).idxExampleCells{iEcc}))
				hFR(k,iCell).FaceColor = colors(min(256, round(mean(fr{k}(iCell, iTick, iTrial), 3) / frMax(k) * 255) + 1), :);
			end
		end
		% hTxt.String = sprintf('Time from Saccade Offset: %d ms', tTicks(iTick));
		hTime1.XData(:) = tTicks(iTick);
		hTime2.XData(:) = tTicks(iTick);
		hTitle.String = sprintf('Eye Trace | t = %d ms', tTicks(iTick));
		drawnow;
		writeVideo(writerObj, getframe(gcf));
		% pause;
	end

	close(writerObj);
end


%% example stimuli
figure('NumberTitle', 'off', 'name', 'Example Stimuli', 'color', 'w'); pause(0.1); jf = get(handle(gcf),'javaframe'); jf.setMaximized(1); pause(1);
[noise, inputX, inputY] = encoder.LoadNoise( fullfile('../../Data', 'noise/noise_1.bin'), 1/60 );
grating2 = encoder.GenerateGrating( 2, 0, 1/60, -1, 45 );	% uniform grating across visual field
grating10 = encoder.GenerateGrating( 10, 0, 1/60, -1, 45 );	% uniform grating across visual field
xIdx = 850:1050;
yIdx = 565:635;
subplot(1,3,1);
imshow((1.2*noise(yIdx,xIdx) + 0.8*grating2(yIdx,xIdx)) * 0.25 + 0.5);
subplot(1,3,2);
imshow(1.2*noise(yIdx,xIdx) * 0.25 + 0.5);
subplot(1,3,3);
imshow((1.2*noise(yIdx,xIdx) + 0.8*grating10(yIdx,xIdx)) * 0.25 + 0.5);



%% video showing dynamics of the sensitivity maps across the visual field out of linear interpolation
for(isNorm2End = [false true])
	for(isEqualAxis = [true false])
		for(isLog = [true false])

			[X, Y] = meshgrid(-14:14, -14:14);
			ecc_ = 0:0.1:20;
			names = {'POn', 'POff', 'MOn', 'MOff'};
			density_ = 0;
			density = 0;
			for(k = 1:4)
				density_ = density_ + WatsonRGCModel.RFSpacingDensityMeridian( ecc_, WatsonRGCModel.enumeratedMeridianNames{1}, names{k} );   % temporal meridian
				density = density + WatsonRGCModel.RFSpacingDensity([X(:), Y(:)], names{k});
			end
			temporalEccDegs = interp1( density_, ecc_, density, 'linear', 'extrap' );
			locIdx = temporalEccDegs <= 14;

			% plot cell locations
			figure('numbertitle', 'off', 'name', sprintf('Temporal Equivalent Locations of Cells'), 'color', 'w');  pause(0.1); jf = get(handle(gcf),'javaframe'); jf.setMaximized(1); pause(1);
			hold on;
			plot(X(locIdx), Y(locIdx), 'bo', 'linewidth', 2, 'displayname', 'original');
			plot(X(locIdx)./(eps+sqrt(X(locIdx).^2+Y(locIdx).^2)).*temporalEccDegs(locIdx), Y(locIdx)./(eps+sqrt(X(locIdx).^2+Y(locIdx).^2)).*temporalEccDegs(locIdx), 'rx', 'linewidth', 2, 'displayname', 'temporal equivalence');
			xlabel('Horizontal position (\circ)');
			xlabel('Vertical position (\circ)');
			legend('location', 'northeastoutside');
			axis equal;
			set(gca, 'fontsize', 20, 'linewidth', 2);

			% generate movie of dynamics of sensitivity map
			filename = fullfile('../../Manuscript/FB Renewal 2021/Videos', sprintf('Dynamics of Sensitivity Map - isNorm2End=%d - isEqualAxis=%d - isLog=%d', isNorm2End, isEqualAxis, isLog));
			writerObj = VideoWriter(filename, 'MPEG-4');
			% writerObj.FrameRate = 15;
			open(writerObj);
			figure('numbertitle', 'off', 'name', sprintf('Movie Showing Dynamics of Sensitivity Maps - isNorm2End=%d - isEqualAxis=%d - isLog=%d', isNorm2End, isEqualAxis, isLog), 'color', 'w');  pause(0.1); jf = get(handle(gcf),'javaframe'); jf.setMaximized(1); pause(1);
			clear h;
			tTicks = 0 : durs(end);
			clear data;
			for(k = 1 : size(sen2,2))
				data{1}(:,k) = interp1(durs, sen2(:,k), tTicks, 'linear');
				data{2}(:,k) = interp1(durs, sen10(:,k), tTicks, 'linear');
			end
			if(isNorm2End)
				data = {data{1} ./ data{1}(end,:), data{2} ./ data{2}(end,:)};
			end
			if(isLog)
				data = {log10(data{1}), log10(data{2})};
			end
			SFs = [2 10];
			for(iTick = 1 : size(tTicks, 2))
				for(iSF = 1 : 2)
					senMap = ones(size(X)) * eps;
					senMap(locIdx) = interp1(unique([conditions.eccentricity]), data{iSF}(iTick,:), temporalEccDegs(locIdx), 'linear');

					subplot(1,2,iSF);
					[~, h(iSF)] = contour( -14:14, -14:14, senMap, 100, 'LineStyle', 'none', 'fill', 'on' );
			        h = handle(h);
					hBar = colorbar;
					hBar.Label.String = 'Sensitivity';
					if(isLog)
						hBar.Label.String = ['Log ' hBar.Label.String];
					end
					if(isNorm2End)
						hBar.Label.String = ['Normalized ' hBar.Label.String];
					end
					hBar.FontSize = 20;
					if(isEqualAxis)
						caxis([min(reshape([data{1}(2:end,:); data{2}(2:end,:)], 1, [])), max([data{1}(:); data{2}(:)])]);
					else
						caxis([min(reshape([data{1}(2:end,:); data{2}(2:end,:)], 1, [])), max(data{iSF}(:))]);
					end
					colormap('hot');
					title(sprintf('SF = %d', SFs(iSF)));
					xlabel('Horizontal position (\circ)');
					ylabel('Vertical position (\circ)');
					axis equal;
					set(gca, 'FontSize', 20, 'lineWidth', 2, 'color', 'k');
				end
				axes('position', [0 0 1 1], 'visible', 'off');
				hTxt = text(0.5, 0.9, sprintf('Post-Saccade Exposure: %d ms', tTicks(iTick)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 24);
			    drawnow;
				writeVideo(writerObj, getframe(gcf));
			end
			close(writerObj);
		end
	end
end



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
