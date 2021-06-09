
%%
% encoder = Encoder([], true);
% encoder = Encoder([], false);   % H_s not 1

%%
sbjs = {'A0HL', 'A092', 'A0NK'};
for(iSbj = 1:nSbj)
	sbj = sbjs{iSbj};
	for(fd = {	...
				'Noise & Grating Simulated Separately', ...
				...%'Noise & Grating Simulated Separately - Hs Not 1', ...
				...%'Noise & Grating Simulated Separately - Drift Stabilized', ...
				% 'Noise & Grating Simulated Separately - Saccade Stabilized (static+static+drift)', ...
				})

		folder = fd{1};

		withInternalNoise = true;
		% durOffset = 0;
		durOffset = 57;
		
		% encoder.SimulateExampleCellsActivities([1 4 7], '../../Data', [], 'none', folder);
		% encoder.SimulateExampleCellsActivities([1 4 7], '../../Data', [], 'drift', folder);
		% encoder.SimulateExampleCellsActivities([1 4 7], '../../Data', [], 'saccade', folder);

		% gather data
		% encoder.GatherSimulatedExampleCellsActivities('../../Data', folder);
	% 	encoder.LoadExampleCellsActivities(fullfile('../../Data/Simulated Activities', folder, [folder '.mat']));

		% activity profiles
		if(withInternalNoise)
			figFolder = fullfile('../../Data/Simulated Activities', sbj, folder, 'figures - withInternalNoise');
		else
			figFolder = fullfile('../../Data/Simulated Activities', sbj, folder, 'figures');
	    end
	    if(~exist(figFolder, 'dir'))
	        mkdir(figFolder);
	    end
		encoder.MeanPopulationActivity(withInternalNoise); drawnow;
		saveas(gcf, fullfile(figFolder, 'Mean+SEM of Population Activity of Example Cells - Contrast=0.5.fig'));
		saveas(gcf, fullfile(figFolder, 'Mean+SEM of Population Activity of Example Cells - Contrast=0.5.png'));
		encoder.AccumulatedPopulationActivity(withInternalNoise); drawnow;
		saveas(gcf, fullfile(figFolder, 'Mean+SEM of Mean Accumulated Population Activity of Example Cells - Contrast=0.5.fig'));
		saveas(gcf, fullfile(figFolder, 'Mean+SEM of Mean Accumulated Population Activity of Example Cells - Contrast=0.5.png'));
		% encoder.DetectionDynamics(durOffset);

		clc;
		encoder.ContrastDetection('thresholding-uni', 0.1, fullfile('../../Data/Simulated Activities', folder), durOffset, withInternalNoise);
		% encoder.ContrastDetection('thresholding-uni-optim_fa', 0.1, fullfile('../../Data/Simulated Activities', folder), durOffset, withInternalNoise);
		% encoder.ContrastDetection('thresholding-sw_sine-tw_uni', 0.1, fullfile('../../Data/Simulated Activities', folder), durOffset, withInternalNoise);
		% encoder.ContrastDetection('thresholding-sw_uni-tw_dprime', 0.1, fullfile('../../Data/Simulated Activities', folder), durOffset, withInternalNoise);
		% encoder.ContrastDetection('thresholding-sw_uni-tw_pval', 0.1, fullfile('../../Data/Simulated Activities', folder), durOffset, withInternalNoise);
	end
end



%% compare model and human
sbj = 'A014';
root = '../../Data';
folder = 'Noise & Grating Simulated Separately';
load(fullfile(root, 'Simulated Activities', sbj, folder, 'figures - withInternalNoise/thresholding-uni - durOffset=57/PerformanceData.mat'));

load(fullfile(root, 'Data/statsTable.mat'));
figure('color', 'w');
colors = {'r', 'g', 'b', [1 0.5 0.5], [0.5 1 0.5], [0.5 0.5 1], [0.5 0 0], [0 0.5 0], [0 0 0.5]};
SFs = [2 10]; nSFs = size(SFs,2);
Eccs = [0 4 8]; nEccs = size(Eccs,2);
for(iSF = 1:nSFs)
	subplot(1,2,iSF); hold on; h = [];
	ind = cell(1,nEccs);
	modelSen = zeros(nEccs, 4);
	for(iEcc = 1:nEccs)
		empiricalSen = statsTable.logSens{find(statsTable.spatialFreq == SFs(iSF) & statsTable.eccentricity == Eccs(iEcc), 1, 'last')};
		empiricalSenSD = statsTable.logSensSD{find(statsTable.spatialFreq == SFs(iSF) & statsTable.eccentricity == Eccs(iEcc), 1, 'last')};

        ind{iEcc} = statsTable.Subj{find(statsTable.spatialFreq == SFs(iSF) & statsTable.eccentricity == Eccs(iEcc), 1, 'first')};
		for(iSbj = 1:size(ind{iEcc},1))
			if(strcmpi(sbj, ind{iEcc}.subject{iSbj}))
				matchSbjIdx = iSbj;
				h(2) = plot([50 150 500], [ind{iEcc}.logSens50(iSbj), ind{iEcc}.logSens150(iSbj), ind{iEcc}.logSens500(iSbj)], 'd-', 'color', colors{3+iEcc}, 'MarkerSize', 14, 'LineWidth', 2, 'DisplayName', sprintf('Ecc=%d, Empirical Individual', Eccs(iEcc)));
			else
				;% plot([50 150 500], [ind{iEcc}.logSens50(iSbj), ind{iEcc}.logSens150(iSbj), ind{iEcc}.logSens500(iSbj)], 'd--', 'color', colors{3+iEcc}, 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', sprintf('Ecc=%d, Empirical', Eccs(iEcc)));
			end
		end

		h(3) = plot([50 150 500]+20, empiricalSen, 'o-', 'color', colors{iEcc}, 'MarkerFaceColor', 'none', 'MarkerSize', 14, 'LineWidth', 2, 'DisplayName', sprintf('Ecc=%d, Empirical Population', Eccs(iEcc)));
		plot(reshape([1; 1; NaN]*([50 150 500]+20), 1, []), reshape([empiricalSen-empiricalSenSD; empiricalSen+empiricalSenSD; NaN(1,3)], 1, []), '-', 'color', colors{iEcc}, 'LineWidth', 2);

		modelSen(iEcc,:) = 1./squeeze(Thresholds(5,iSF,iEcc,:));
		modelSen(iEcc,:) = modelSen(iEcc,:) / modelSen(iEcc,end) * ind{iEcc}.logSens500(matchSbjIdx);
		% modelSen(iEcc,:) = modelSen(iEcc,:) / (1./squeeze(Thresholds(5,iSF,1,end))) * ind{iEcc}.logSens500(matchSbjIdx);
		h(1) = plot(durs(2:end), modelSen(iEcc, 2:end), 'd--', 'color', colors{6+iEcc}, 'MarkerSize', 16, 'LineWidth', 2, 'DisplayName', sprintf('Ecc=%d, Model', Eccs(iEcc)));
	end

	xlabel('Post-saccadic duration (ms)');
	ylabel('Sensitivity');
	title(sprintf('SF = %d', SFs(iSF)));
	set(gca, 'XTick', [0 50 150 500], 'XLim', [0 550], 'YScale', 'log', 'FontSize', 20, 'LineWidth', 2);
    set(gca, 'YTickLabel', sprintf('%g\n',(get(gca, 'YTick'))));
    if(iSF == 1)
    	legend(h, 'location', 'best');
    end
end