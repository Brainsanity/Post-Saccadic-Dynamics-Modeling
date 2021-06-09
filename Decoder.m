classdef Decoder < handle

	properties (SetAccess = public)%private)
		encoder = Encoder();
	end


	methods
		function obj = Decoder(stabilize)
			%% Constructor function
			%   stabilize:			'normal', 'drift' stabilized, or 'saccade' stabilized

			if( ~exist('stabilize', 'var') || isempty(stabilize) )
				stabilize = 'normal';
			end

			paramFolder = './Parameters/';
			HsBe1 = true;
			obj.encoder = Encoder(paramFolder, HsBe1);

			sbj = 'SacDB';
			switch stabilize
				case 'drift'
					folder = 'UG - Noise & Grating Simulated Separately - drift';
				case 'saccade'
					folder = 'UG - Noise & Grating Simulated Separately - saccade';
				otherwise
					folder = 'UG - Noise & Grating Simulated Separately';
			end
			inFolder = 'UG - Noise & Grating Simulated Separately';
			
			obj.encoder.LoadExampleCellsActivities(fullfile('../../Data/Simulated Activities', sbj, folder, ['All Conditions' '.mat']));
			obj.encoder.AddInternalNoise([], fullfile('../../Data/Simulated Activities', sbj, inFolder, ['All Conditions' '.mat']));
		end


		function [tTicks, Thresholds, ThresholdsSTD, Sensitivities, SensitivitiesSTD] = ContrastThreshold(obj, alignEvent, durMax, tStep, dataFolder, durOffset, withInternalNoise)
			
			if(~exist('alignEvent', 'var') || isempty(alignEvent))
				alignEvent = 'saccadeOff';
			end
			if( ~exist('durMax', 'var') || isempty(durMax) )
				durMax = 600;
			end
			if( ~exist('tStep', 'var') || isempty(tStep) )
				tStep = 10;
			end
			if( ~exist('dataFolder', 'var') || isempty(dataFolder) )
				dataFolder = '../../Data/Simulated Activities/SacDB/UG-Noise & Grating Simulated Separately';
			end
			if( ~exist('durOffset', 'var') || isempty(durOffset) )
				durOffset = 0;	%50 + 7;			% resonse delay of 50 ms + online saccade off later by 7 ms
			end
			if( ~exist('withInternalNoise', 'var') )
				withInternalNoise = true;
			end

			classifier = 'uni-no_bias-fa25_hit75';

			if(withInternalNoise)
				dataFolder = fullfile(dataFolder, 'figures - withInternalNoise', sprintf('%s - durOffset=%d', classifier, durOffset));
			else
				dataFolder = fullfile(dataFolder, 'figures', sprintf('%s - durOffset=%d', classifier, durOffset));
			end
			if(~exist(dataFolder, 'dir'))
				mkdir(dataFolder);
			end

			sbj = 'SacDB';
			if(strcmpi(classifier, 'thresholding-uni-optim_fa'))		% optimal false alarm rate from distributions of human experimental conditions
				if(~isfield(obj.encoder.activityParams, 'FA') || ~obj.encoder.activityParams.FA.isKey(sbj) || ~strcmpi(obj.encoder.activityParams.FA(sbj).alignEvent, alignEvent) || obj.encoder.activityParams.FA(sbj).T(1) > durMax(1) || obj.encoder.activityParams.FA(sbj).T(end) < durMax(2))
					obj.encoder.FalseAlarmRate(sbj, [], durOffset, alignEvent, durMax);
				end
			else
				[fpr, nNoPresent, durs, ~] = EmpiricalBox.FalsePositiveRate('A014');
				fpr = nansum(fpr .* nNoPresent, 2) ./ sum(nNoPresent, 2);			% all conditions are mixed within blocks, therefore the non-present condition has no eccentricity
				fprFun = @(dur) max(0, interp1(durs(:,1), fpr, dur, 'linear', 'extrap'));		% get false positive rate according to duration and eccentricity
			end

			nBoots = 4;
			nTrials = size(obj.encoder.activityParams.trials,2);

			conditions = obj.encoder.activityParams.conditions;
			Eccs = unique([conditions.eccentricity]);
			nEccs = size(Eccs,2);

			SFs = unique([conditions.sf]);
			SFs(SFs == 0) = [];
			nSFs = size(SFs,2);

			tTicks = 0 : tStep : durMax;
			durs = tTicks;

			load( fullfile(dataFolder, 'PerformanceData.mat') );
			% Thresholds = zeros(5, nSFs, nEccs, size(tTicks,2));
			% ThresholdsSTD = Thresholds;
			% Sensitivities = Thresholds;
			% SensitivitiesSTD = Thresholds;

			cellNumAmplifier =cat(1, obj.encoder.layers.nAllCells) ./ cat(1,obj.encoder.layers.nExampleCells);		% inverse of proportion of cells used
			cellNumAmplifier(5,:) = mean(cellNumAmplifier,1);

			for(iL = 1 : 5)
				for(iSF = 1 : nSFs)
					for(iEcc = 1 : nEccs)
						fprintf('iL = %d, SF = %d, Ecc = %d...\n', iL, SFs(iSF), Eccs(iEcc));
						if(iL ~= 5)
							frAbsent = obj.encoder.ExampleCellsActivitiesOnCondition(obj.encoder.layers(iL).name, Eccs(iEcc), 0, 0, alignEvent, [0 durMax], false, withInternalNoise);
							lfrPresent = obj.encoder.ExampleCellsActivitiesOnCondition(obj.encoder.layers(iL).name, Eccs(iEcc), SFs(iSF), 1.0, alignEvent, [0 durMax], false, false);
						else
							frAbsent = cat(1, obj.encoder.ExampleCellsActivitiesOnCondition(obj.encoder.layers(1).name, Eccs(iEcc), 0, 0, alignEvent, [0 durMax], false, withInternalNoise),...
											  obj.encoder.ExampleCellsActivitiesOnCondition(obj.encoder.layers(2).name, Eccs(iEcc), 0, 0, alignEvent, [0 durMax], false, withInternalNoise),...
											  obj.encoder.ExampleCellsActivitiesOnCondition(obj.encoder.layers(3).name, Eccs(iEcc), 0, 0, alignEvent, [0 durMax], false, withInternalNoise),...
											  obj.encoder.ExampleCellsActivitiesOnCondition(obj.encoder.layers(4).name, Eccs(iEcc), 0, 0, alignEvent, [0 durMax], false, withInternalNoise));
							lfrPresent = cat(1, obj.encoder.ExampleCellsActivitiesOnCondition(obj.encoder.layers(1).name, Eccs(iEcc), SFs(iSF), 1.0, alignEvent, [0 durMax], false, false),...
												obj.encoder.ExampleCellsActivitiesOnCondition(obj.encoder.layers(2).name, Eccs(iEcc), SFs(iSF), 1.0, alignEvent, [0 durMax], false, false),...
												obj.encoder.ExampleCellsActivitiesOnCondition(obj.encoder.layers(3).name, Eccs(iEcc), SFs(iSF), 1.0, alignEvent, [0 durMax], false, false),...
												obj.encoder.ExampleCellsActivitiesOnCondition(obj.encoder.layers(4).name, Eccs(iEcc), SFs(iSF), 1.0, alignEvent, [0 durMax], false, false));
						end
						for(iTick = 1 : size(tTicks,2))
							if(iTick == 1), tic; end
							fprintf('\tduration = %d...\n', tTicks(iTick));

							Y0_ = frAbsent(:, 1 : tTicks(iTick)+1, :);
							tw = ones(1, size(Y0_,2)) / size(Y0_,2);		% uniform temporal weights
							Y0_ = tw * shiftdim(nanmean(Y0_, 1), 1);

							Y1_ = lfrPresent(:, 1 : tTicks(iTick)+1, :);

							if(iTick == 1)
								c = ones(1, nBoots) * 0.5;	% initial contrast
							else
								c = ones(1, nBoots) * Thresholds(iL,iSF,iEcc,iTick-1);
							end
							for(iBoot = 1 : nBoots)
							% for(k = 1 : 1000)
							% 	parfor(iBoot = (k-1)*nBoots/1000+1 : k*nBoots/1000)
									if(~mod(iBoot-1, round(nBoots/10)))
										fprintf('\t\tiBoot = %d/%d...', iBoot, nBoots);
									end
									Y0 = Y0_(randi(nTrials, 1, nTrials));
									db = prctile(Y0, 75);						% decision boundary at 25% false-alarm rate
									
									cb = [0 Inf];	% initial contrast boundaries
									% c(iBoot) = 0.5;		% initial contrast
									extScale = 2;	% scale for extending over cb(2)
									while(abs(diff(cb)) / c(iBoot) > 0.01)
										Y1 = obj.encoder.AddInternalNoise(c(iBoot)*Y1_(:,:,randi(nTrials, 1, nTrials)));	% add non-linearity and internal noise
										Y1 = tw * shiftdim(nanmean(Y1, 1), 1);
										if(sum(Y1 < db) / size(Y1,2) < 0.25)		% miss rate should be 25%
											% need to decrease the contrast
											cb(2) = c(iBoot);
											c(iBoot) = mean(cb);
										else
											% need to increase the contrast
											cb(1) = c(iBoot);
											if(cb(2) == Inf)
												c(iBoot) = c(iBoot) * extScale;
												extScale = extScale * 2;
											else
												c(iBoot) = mean(cb);
											end
										end
									end

									if(~mod(iBoot-1, round(nBoots/10)))
										fprintf(' | t = %.6f\n', toc); tic;
									end
	                            % end
	                        end
                            
                            c(isoutlier(c)) = [];

							Thresholds(iL,iSF,iEcc,iTick) = mean(c);
							ThresholdsSTD(iL,iSF,iEcc,iTick) = std(c);
							Sensitivities(iL,iSF,iEcc,iTick) = mean(1./c);
							SensitivitiesSTD(iL,iSF,iEcc,iTick) = std(1./c);
						end
					end
				end
			end
				
			save( fullfile(dataFolder, 'PerformanceData.mat'), 'tTicks', 'conditions', 'durs', 'Thresholds', 'ThresholdsSTD', 'Sensitivities', 'SensitivitiesSTD' );

		end
	end
end