classdef EmpiricalBox < handle

	methods (Access = private)
		function obj = EmpiricalBox()
		end
	end

	methods (Static)
		function trials = LoadSingleData(dataFile, isDriftOnly)
			%% load data for one single subject
			if(~exist('isDriftOnly', 'var') || isempty(isDriftOnly))
				isDriftOnly = true;
			end

			% load experiment data
			load( dataFile, 'ppt', 'eminfo', 'counter' );
			trials = [ppt{:}];
			if(isDriftOnly)
				idx = (eminfo.DriftOnly') & [trials.contrast] <= 0.5;% & ~isnan(eminfo.SaccadeStart');
                eminfo.SaccadeStart = eminfo.SaccadeStart(idx);
				trials = trials(idx);
            else
            	idx = [trials.contrast] <= 0.5;% & ~isnan(eminfo.SaccadeStart');
                eminfo.SaccadeStart = eminfo.SaccadeStart(idx);
				trials = trials(idx);
			end

			% set flashOn to the middle of saccade
			for( iTrial = 1 : size(trials,2) )
				trials(iTrial).flashOn = find( trials(iTrial).x.position(eminfo.SaccadeStart(iTrial):end) <= 200, 1, 'first' ) + eminfo.SaccadeStart(iTrial)-1;		% in samples
			end
			
			% redo saccade detection
			trials = SaccadeTool.GetSacs( trials, 'minmsa', 3 );
            flags = false(size(trials));
			for( iTrial = size(trials,2) : -1 : 1 )
                if(isempty(trials(iTrial).saccades.start))
                    flags(iTrial) = true;
                    continue;
                end
				[~, saccadeIdx(iTrial)] = min(abs( trials(iTrial).saccades.start - trials(iTrial).flashOn /trials(iTrial).sRate*1000 ));
				saccadeOn(iTrial) = trials(iTrial).saccades.start(saccadeIdx(iTrial));																		% in samples
				saccadeOff(iTrial) = trials(iTrial).saccades.start(saccadeIdx(iTrial)) + trials(iTrial).saccades.duration(saccadeIdx(iTrial)) - 1;			% in samples
            end
            trials(flags) = [];
            saccadeOn(flags) = [];
            saccadeOff(flags) = [];
			saccadeOn = num2cell(saccadeOn);
			saccadeOff = num2cell(saccadeOff);
			[trials.saccadeOn] = saccadeOn{:};
			[trials.saccadeOff] = saccadeOff{:};
			trials([trials.saccadeOff] - [trials.saccadeOn] > 150./1000*[trials.sRate]) = [];

			% saccade landing time: first sample within 5 arcmin (horizontally) of saccade target
			for( iTrial = 1 : size(trials,2) )
				trials(iTrial).saccadeLand = find( abs( trials(iTrial).x.position( trials(iTrial).saccadeOn : trials(iTrial).saccadeOff ) ) <= 5, 1, 'first' ) + trials(iTrial).saccadeOn-1;		% in samples
				if( isempty(trials(iTrial).saccadeLand) )
					[~, trials(iTrial).saccadeLand] = min( abs( trials(iTrial).x.position( trials(iTrial).saccadeOn : trials(iTrial).saccadeOff ) ) );
					trials(iTrial).saccadeLand = trials(iTrial).saccadeLand + trials(iTrial).saccadeOn - 1;
				end
			end
		end

		function data = LoadMultipleData(dataFolder, subjects)
			%% load data for specified subjects

			if( nargin() < 2 )	% subjects not specified, load all available subjects in dataFolder
				dataFiles = dir( fullfile(dataFolder, '*.mat') );
				subjects = {dataFiles.name};
			else
				subjects = cellfun( @(sbj) {[sbj,'.mat']}, subjects(:)' );
			end

			data = cellfun( @(sbj) struct( 'sbj', sbj(1:end-4), 'trials', EmpiricalBox.LoadSingleData(fullfile(dataFolder, sbj)) ), subjects );
		end

		function [tpr, fpr, nPresent, nNoPresent] = SinglePerformance(dataFolder, sbj)
			if(~exist('dataFolder', 'var') || isempty(dataFolder))
				dataFolder = '../../Data/Data';
			end

			dataFile = fullfile(dataFolder, [sbj, '.mat']);

			durs = [0 100 250 650];
			% SFs = [2 10];
			eccs = [0 4 8 12];

			trials = EmpiricalBox.LoadSingleData(dataFile);
			for( iDur = size(durs,2)-1 : -1 : 1 )
				durations = [trials.stimOff] - [trials.saccOff];
				idxDur = durs(iDur) < durations & durations < durs(iDur+1);

				for( iEcc = size(eccs,2) : -1 : 1 )
					idxEcc = [trials.eccentricity] == eccs(iEcc);

					idxNoPresent = idxDur & idxEcc & ~[trials.present];
					nNoPresent(iDur,iEcc) = sum(idxNoPresent);
					fpr(iDur,iEcc) = sum( idxNoPresent & [trials.responseValue] ) / nNoPresent(iDur,iEcc);

					idxPresent = idxDur & idxEcc & [trials.present];
					nPresent(iDur,iEcc) = sum(idxPresent);
					tpr(iDur,iEcc) = sum( idxPresent & [trials.responseValue] ) / nPresent(iDur,iEcc);
				end
			end
		end

		function [fpr, nNoPresent, durs, eccs] = FalsePositiveRate(sbj)
			dataFolder = 'F:\Post Saccadic Dynamics Modeling\Data\Data';
			dataFiles = dir( fullfile(dataFolder, 'A*.mat') );
			subjects = {dataFiles.name};

			durs = [0 100 250 650];
			% SFs = [2 10];
			eccs = [0 4 8 12];

			if( exist( fullfile(dataFolder, 'FalsePositiveVSDuration.mat') ) )
				load( fullfile(dataFolder, 'FalsePositiveVSDuration.mat') );
			
			else
				for( iSbj = size(subjects,2) : -1 : 1 )
					trials = EmpiricalBox.LoadSingleData( fullfile(dataFolder, subjects{iSbj}), false );
					
					for( iDur = size(durs,2)-1 : -1 : 1 )
						durations = [trials.stimOff] - [trials.saccOff];
						idxDur = durs(iDur) < durations & durations < durs(iDur+1);

						for( iEcc = size(eccs,2) : -1 : 1 )
							idxEcc = [trials.eccentricity] == eccs(iEcc);

							idxNoPresent = idxDur & idxEcc & ~[trials.present];
							nNoPresent(iSbj,iDur,iEcc) = sum(idxNoPresent);
							fpr(iSbj,iDur,iEcc) = sum( idxNoPresent & [trials.responseValue] ) / nNoPresent(iSbj,iDur,iEcc);
						end
					end
				end
				save( fullfile(dataFolder, 'FalsePositiveVSDuration.mat'), 'fpr', 'nNoPresent' );
			end

			for(iSbj = 1 : size(subjects,2))
				if(strcmpi(sbj, subjects{iSbj}(1:end-4)))
					fpr = squeeze(fpr(iSbj,:,:));
					nNoPresent = squeeze(nNoPresent(iSbj,:,:));
				end
			end

			[eccs, durs] = meshgrid(eccs, [50 150 500]);
		end

		function [fpr, nNoPresent] = StimDur2FA(dataFolder)
			%% analyze correlation between false alarm rate and stimulus duration
			if( nargin() < 1 )
				dataFolder = 'F:\Post Saccadic Dynamics Modeling\Data\Data';
			end

			dataFiles = dir( fullfile(dataFolder, 'A*.mat') );
			subjects = {dataFiles.name};

			durs = [0 100 250 650];
			% SFs = [2 10];
			eccs = [0 4 8 12];

			if( exist( fullfile(dataFolder, 'FalsePositiveVSDuration.mat') ) )
				load( fullfile(dataFolder, 'FalsePositiveVSDuration.mat') );
			
			else
				for( iSbj = size(subjects,2) : -1 : 1 )
					trials = EmpiricalBox.LoadSingleData( fullfile(dataFolder, subjects{iSbj}) );
					
					for( iDur = size(durs,2)-1 : -1 : 1 )
						durations = [trials.stimOff] - [trials.saccOff];
						idxDur = durs(iDur) < durations & durations < durs(iDur+1);

						for( iEcc = size(eccs,2) : -1 : 1 )
							idxEcc = [trials.eccentricity] == eccs(iEcc);

							idxNoPresent = idxDur & idxEcc & ~[trials.present];
							nNoPresent(iSbj,iDur,iEcc) = sum(idxNoPresent);
							fpr(iSbj,iDur,iEcc) = sum( idxNoPresent & [trials.responseValue] ) / nNoPresent(iSbj,iDur,iEcc);
						end
					end
				end
			end

			figure('NumberTitle', 'off', 'name', 'False Alarm Rate VS Stimulus Duration | Pooled Across Conditions', 'color', 'w'); hold on;
			colors = { 'r', 'g', 'b', 'm', 'y', 'c', [0.5 0 0], [0 0.5 0], [0 0 0.5], [0.5 0 0.5], [0.5 0.5 0], [0 0.5 0.5], [1 0.5 0.5], [0.5 1 0.5], [0.5 0.5 1], [1 0.5 1], [1 1 0.5], [0.5 1 1] };
			poolFA = nansum( fpr .* nNoPresent, 3 ) ./ nansum( nNoPresent, 3 );
			for( iSbj = 1 : size(poolFA,1) )
				h(iSbj) = plot( [50 150 500], poolFA(iSbj,:), '^', 'color', colors{iSbj}, 'markersize', 8, 'linewidth', 2, 'displayname', subjects{iSbj}(1:end-4) );
				plot( [50 150 500], poolFA(iSbj,:), '--', 'color', colors{iSbj}, 'linewidth', 1 );
			end
			fa = poolFA;
			fa( any(isnan(fa(:,1:2)), 2), : ) = [];
			m = nanmean(fa,1);
			sem = nanstd(fa,1) ./ sqrt( sum(~isnan(fa), 1) );
			h(end+1) = plot( [50 150 500] + 10, m, 'ks-', 'markersize', 12, 'linewidth', 2, 'displayname', 'Average' );
			plot( reshape([50 150  500; 50 150 500; nan nan nan], 1, []) + 5, reshape([m+sem; m-sem; nan nan nan], 1, [] ), 'k-', 'linewidth', 2 );
			legend(h,'location', 'northeast');
			xlabel('Stimulus duration (ms)');
			ylabel('False alarm rate');
			set(gca, 'xtick', [50 150 500], 'xlim', [0 650], 'fontsize', 20, 'linewidth', 2);

			figure('NumberTitle', 'off', 'name', 'False Alarm Rate VS Stimulus Duration', 'color', 'w'); hold on;
			nRows = 1;
			nCols = size(eccs,2);
			h = [];
			for( iEcc = 1 : size(eccs,2) )
				subplot(nRows, nCols, iEcc); hold on;
				for( iSbj = 1 : size(nNoPresent,1) )
					h(iSbj) = plot( [50 150 500], fpr(iSbj,:,iEcc), '^', 'color', colors{iSbj}, 'markersize', 8, 'linewidth', 2, 'displayname', subjects{iSbj}(1:end-4) );
					plot( [50 150 500], fpr(iSbj,:,iEcc), '--', 'color', colors{iSbj}, 'linewidth', 1 );
				end
				fa = fpr(:,:,iEcc);
				fa( any(isnan(fa(:,1:2)), 2), : ) = [];
				m = nanmean(fa,1);
				sem = nanstd(fa,1) ./ sqrt( sum(~isnan(fpr(:,:,iEcc)), 1) );
				h(end+1) = plot( [50 150 500] + 10, m, 'ks-', 'markersize', 12, 'linewidth', 2, 'displayname', 'Average' );
				plot( reshape([50 150  500; 50 150 500; nan nan nan], 1, []) + 5, reshape([m+sem; m-sem; nan nan nan], 1, [] ), 'k-', 'linewidth', 2 );
				if( iEcc == size(eccs,2) )
					legend(h,'location', 'northeast');
				end
				xlabel('Stimulus duration (ms)');
				if( iEcc == 1 )
					ylabel('False alarm rate');
				end
				title(sprintf('Ecc = %d', eccs(iEcc)));
				set(gca, 'xtick', [50 150 500], 'xlim', [0 650], 'fontsize', 20, 'linewidth', 2);
			end
		end
	end



end