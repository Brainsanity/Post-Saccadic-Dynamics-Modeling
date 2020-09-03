classdef Encoder < handle

	properties (SetAccess = private)
		layers;%(4) = struct( 'name', [], 'locations', [], 'sRFRc', [], 'sRFKc', [], 'sRFRs', [], 'sRFKs', [] );
		SpatialModel;
		TemporalModel;
	end


	methods
		
		function obj = Encoder(paramFolder)
			%% Constructor function
			%	paramFolder:		folder containing cell parameters

			if( nargin() < 1 || isempty(paramFolder) )
				paramFolder = './Parameters/';
			end

			obj.TemporalModel = TemporalRF();
			
			obj.layers = struct( 'name', cell(1,4), 'locations', [], 'sRFParams', [], 'tRFParams', [] );
			names = {'POn', 'POff', 'MOn', 'MOff'};
			[obj.layers.name] = names{:};

			obj.SpatialModel = CronerKaplanRGCModel( 'dataSetToFit', 'raw', 'fitIntercept', false, 'randomizeCenterRadii', true, 'randomizeCenterSensitivities', true, 'randomizeSurroundRadii', true, 'randomizeSurroundSensitivities', true );

			for( k =  1 : size(obj.layers,2) )
				fprintf( 'k = %d\ncellType = %s\n\n', k, names{k} );

				% get cell locations
				if( exist( fullfile( paramFolder, ['Mosaic_' names{k} '_Radius20.0deg_maxMovPrctile20.mat'] ), 'file' ) )
					obj.layers(k).locations = load( fullfile( paramFolder, ['Mosaic_' names{k} '_Radius20.0deg_maxMovPrctile20.mat'] ), 'rfPositions' );
					obj.layers(k).locations = obj.layers(k).locations.rfPositions;
				else
					obj.layers(k).locations = WatsonRGCModel.CreateMosaic(20, names{k}, 'saveFolder', paramFolder, 'maxIterations', 3000);
				end

				% get spatial RF parameters
				if( exist( fullfile( paramFolder, ['SpatialRFParams_', names{k}, '.mat'] ), 'file' ) )
					obj.layers(k).sRFParams = load( fullfile( paramFolder, ['SpatialRFParams_', names{k}, '.mat'] ), 'sRFParams' );
					obj.layers(k).sRFParams = obj.layers(k).sRFParams.sRFParams;
				else
					sRFParams = obj.SpatialModel.SynthesizeRFParams( sqrt(sum(obj.layers(k).locations.^2,2))', names{k} );
					save( fullfile( paramFolder, ['SpatialRFParams_', names{k}, '.mat'] ), 'sRFParams' );
					obj.layers(k).sRFParams = sRFParams;
				end

				% get temporal RF parameters
				obj.layers(k).tRFParams = obj.TemporalModel.SynthesizeParams( names{k}, size(obj.layers(k).locations, 1) );
			end
		end


		function [LFR, time, conditions, sFR, sFR_c, sFR_s, trials, trialsIdx, cellIdx, nCells, nCellsUsed] = ExampleCells(obj, dataFolder, alignEvent, contrast, stabilize)
			%% Show example modeled cells at locations (0,0), (0,4), (0,8), (0,12)
			%   dataFolder:			folder containing experiment data and noise background
			%	alignEvent:			'saccadeOn', 'saccadeOff', 'flashOn'
			%   contrast:			contrast of grating
			%   stabilize:			'normal' (no stabilization), 'drift' (only stabilize drift), or 'full' (full stabilization)

			if( nargin() < 2 || isempty(dataFolder) )
				dataFolder = '../../data/';
			end
			if( nargin() < 3 || isempty(alignEvent) )
				alignEvent = 'saccadeOff'; %'saccadeLand';%
			end
			if( nargin() < 5 || isempty(stabilize) )
				stabilize = 'normal';
			end

			trials = obj.LoadExpData(dataFolder);
			
			% get example cell indices
			eccs = [0, 4, 8, 12];
			for( iL = size(obj.layers,2) : -1 : 1 )
				index = obj.layers(iL).locations(:,1) >= 0 & abs( obj.layers(iL).locations(:,2) ) <= 0.5;		% right visual field and vertically within +-0.5 deg
				for( iEcc = size(eccs,2) : -1 : 1 )
					d2 = sum( obj.layers(iL).locations.^2, 2 );
					if( eccs(iEcc) == 0 )
						cellIdx{iL,iEcc} = find( d2 <= trials(1).gratingWidth^2 & index );
						nCells{iL,iEcc} = sum( d2 <= trials(1).gratingWidth^2 );
					else
						w = trials(1).gratingWidth/2;
						cellIdx{iL,iEcc} = find( (eccs(iEcc)-w)^2 <= d2 & d2 <= (eccs(iEcc)+w)^2 & index );
						nCells{iL,iEcc} = sum( (eccs(iEcc)-w)^2 <= d2 & d2 <= (eccs(iEcc)+w)^2 );
					end
					cellIdx{iL,iEcc} = cellIdx{iL,iEcc}( randperm( size(cellIdx{iL,iEcc}, 1) ) );
					nCellsUsed{iL,iEcc} = round( nCells{iL,iEcc} / 100 );
					% [~, cellIdx(iL,iEcc)] = min( sum( (obj.layers(iL).locations - cLocs(iEcc,:)).^2, 2 ) );
				end
			end
			% nCells = 30;
            nTrials = 30;


			% obj.SpatialModel = CronerKaplanRGCModel( 'dataSetToFit', 'raw', 'fitIntercept', false, 'randomizeCenterRadii', true, 'randomizeCenterSensitivities', true, 'randomizeSurroundRadii', true, 'randomizeSurroundSensitivities', true );
			% obj.TemporalModel = TemporalRF();

			%% compute responses
% 			eccs = [0, 4, 8, 12];
% 			SFs = [2,10];
% % 			durs = [50, 150, 500];
%             durs = [500];
% % 			reports = [0,1];
% 			presents = [0,1];	% whether grating displayed or absent
% 			categories(1).name = 'Eccentricity';
% 			categories(1).values = eccs;
% 			categories(2).name = 'Spatial Frequency';
% 			categories(2).values = SFs;
% 			categories(3).name = 'Duration';
% 			categories(3).values = durs;
% 			categories(4).name = 'Present';
% 			categories(4).values = presents;

			conditions = struct( ...
							'eccentricity',	{  0,   4,   8,  12,   0,   4,   8,  12,   0,   4,   8,  12}, ...
							'sf',			{  2,   2,   2,   2,   2,   2,   2,   2,  10,  10,  10,  10}, ...
							'duration',		{500}, ...
							'present',		{  0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1}, ...
							'alignEvent',	{alignEvent} );
			conditions([4:4:end]) = [];
			Eccs = unique([conditions.eccentricity]);

			for( iCond = size(conditions,2) : -1 : 1 )
				iEcc = find( conditions(iCond).eccentricity == Eccs );

				trialsIdx{iCond} = find( ... 
										 ...%[trials.eccentricity] == eccs(1,iEcc) & ...
										 ... %[trials.spatialFreq] == SFs(iSF) & ...
										 abs( [trials.stimOff] - [trials.saccadeOff] - conditions(iCond).duration ) < 50 & ...
										 true );%[trials.present] == presents(iPre) );

				idx = trialsIdx{iCond};    %fprintf('nTrials = %d\n', size(idx,2)); continue;
                idx = idx( 1 : min(nTrials,end) );
				tMax = round(max( [trials(idx).saccadeOff] + round(conditions(iCond).duration/1000*[trials(idx).sRate]) - [trials(idx).(alignEvent)] ));	% in samples
				tMin = round(min( [trials(idx).saccadeOn] - round(0.150*[trials(idx).sRate]) - [trials(idx).(alignEvent)] ));				% in samples
				for( iL = 4 : -1 : 1 )
					sFR{iCond,iL} = zeros( nCellsUsed{iL,iEcc}, tMax-tMin+1, size(idx,2) );		% 1st dim: neurons;		2nd dim: time;		3rd dim: trials;
				end
				sFR_c{iCond,iL} = sFR{iCond,iL};
				sFR_s{iCond,iL} = sFR{iCond,iL};
				LFR{iCond,iL} = sFR{iCond,iL};
				for( k = size(idx,2) : -1 : 1 )
					% tic;
					index = tMin + trials(idx(k)).(alignEvent) : min( size(trials(idx(k)).x.position,2), tMax + trials(idx(k)).(alignEvent) );
					x  = trials(idx(k)).x.position(index) / 60;
					y  = trials(idx(k)).y.position(index) / 60;
                    
					% before saccade
					eyeIdx{1} = 1 : trials(idx(k)).saccadeOn - trials(idx(k)).(alignEvent) - tMin;
					s1 = size(eyeIdx{1},2);

					% during saccade
					eyeIdx{2} = eyeIdx{1}(end) + ( 1 : trials(idx(k)).saccadeOff - trials(idx(k)).saccadeOn + 1 );
					s2 = size(eyeIdx{2},2);

					% after saccade
					eyeIdx{3} = eyeIdx{2}(end)+1 : size(x,2);
					s3 = size(eyeIdx{3},2);

					if( strcmpi( stabilize, 'drift' ) )
						x(eyeIdx{1}) = x(eyeIdx{1}(end));
						y(eyeIdx{1}) = y(eyeIdx{1}(end));
						x(eyeIdx{3}) = x(eyeIdx{3}(1));
						y(eyeIdx{3}) = y(eyeIdx{3}(1));
					elseif( strcmpi( stabilize, 'full' ) )
						x(:) = x(eyeIdx{3}(1));
						y(:) = y(eyeIdx{3}(1));
					end

                    [noise, inputX, inputY] = obj.LoadNoise( fullfile(dataFolder, trials(idx(k)).backgroundImage), trials(idx(k)).pixelAngle/60 );
					grating = obj.GenerateGrating( conditions(iCond).eccentricity, trials(idx(k)).gratingWidth, conditions(iCond).sf, trials(idx(k)).phase, trials(idx(k)).pixelAngle/60 );
					noise = noise * trials(idx(k)).backgroundContrast;
					grating = grating * contrast * conditions(iCond).present;
                    bg = 1;
					stimulus = noise + grating;
                    
					% avContrast = ((s1+s2)*trials(idx(k)).backgroundContrast + s2*contrast*conditions(iCond).present) / (s1+s2);
					avContrast = trials(idx(k)).backgroundContrast;
					% contrastF = [ zeros(1,s1), ones(1,s2) * contrast ] + trials(idx(k)).backgroundContrast;
					% contrastF = obj.ComputeContrasts( noise+grating+bg,  inputX, inputY, obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),1), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),2), [obj.layers(iL).sRFParams(cellIdx{iL,iEcc}(1:nCells)).surroundRadii], x([eyeIdx1,eyeIdx2]), y([eyeIdx1,eyeIdx2]) );
					contrastF = 1;%contrast;

					for( iL = 1:4 )
						tic;
						for( m = 1 : 3 )
							xIdx = min(x(eyeIdx{m})) + min(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),1)) - 4 <= inputX & inputX <= max(x(eyeIdx{m})) + max(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),1)) + 4; %720 : 1460;
		                    yIdx = min(y(eyeIdx{m})) + min(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),2)) - 4 <= inputY & inputY <= max(y(eyeIdx{m})) + max(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),2)) + 4; %300 : 780;
	                    
							[sFR{iCond,iL}(:,eyeIdx{m},k), sFR_c{iCond,iL}(:,eyeIdx{m},k), sFR_s{iCond,iL}(:,eyeIdx{m},k)] = obj.SpatialModel.LinearResponse( stimulus(yIdx,xIdx), inputX(xIdx), inputY(yIdx), x(eyeIdx{m}), y(eyeIdx{m}), obj.layers(iL).sRFParams(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc})), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),1), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),2));
						end
						if( lower(obj.layers(iL).name(1)) == 'm' )
							LFR{iCond,iL}(:,:,k) = obj.TemporalModel.LinearResponse( obj.layers(iL).name, obj.layers(iL).tRFParams(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc})), avContrast, trials(idx(k)).sRate, contrastF .* sFR{iCond,iL}(:,:,k) );
						else
							LFR{iCond,iL}(:,:,k) = obj.TemporalModel.LinearResponse( obj.layers(iL).name, obj.layers(iL).tRFParams(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc})), avContrast, trials(idx(k)).sRate, contrastF .* sFR_c{iCond,iL}(:,:,k), contrastF .* sFR_s{iCond,iL}(:,:,k) );
						end
						fprintf('Ecc=%d, SF=%d, Dur=%d, Present=%d, iL=%d, iTrials=%d/%d, t=%f\n', conditions(iCond).eccentricity, conditions(iCond).sf, conditions(iCond).duration, conditions(iCond).present, iL, k, size(idx,2), toc);
					end
					
				end
				if(isempty(idx))
					time{iCond} = [];
				else
					time{iCond} = (tMin:tMax) / trials(idx(1)).sRate * 1000;
				end
            end
            
            save( fullfile( dataFolder, 'figures', sprintf( 'Example Cells LFR; contrast = %f.mat', contrast ) ), 'LFR', 'time', 'conditions', 'sFR', 'sFR_c', 'sFR_s', 'trials', 'trialsIdx', 'cellIdx' );
		end


		function DisplayExampleCells(obj, FR, time, conditions, alignEvent, trials, trialsIdx, cellIdx, applyRectify)
			if( nargin() < 9 || isempty(applyRectify) )
				applyRectify = true;
			end


			% conditions = struct( ...
			% 				'eccentricity',	{  0,   4,   8,  12,   0,   4,   8,  12,   0,   4,   8,  12}, ...
			% 				'sf',			{  2,   2,   2,   2,   2,   2,   2,   2,  10,  10,  10,  10}, ...
			% 				'duration',		{500}, ...
			% 				'present',		{  0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1}, ...
			% 				'alignEvent',	{'flashOn'} );
			
			% nTrials = size( FR{1}, 1 );
			% nCells = size( FR{1}, 3 );
			nTrials = size( FR{1}, 3 );
			nCells = size( FR{1}, 1 );
			idx = trialsIdx{1}(1:nTrials);

			Eccs = unique([conditions.eccentricity]);

			nRows = size( Eccs, 2 );
			nCols = 1;

			for( iL = 1 : 4 )
				figure( 'NumberTitle', 'off', 'name', ['Example Cells: ', obj.layers(iL).name], 'color', 'w' );
				pause(0.1);
				jf = get(handle(gcf),'javaframe');
				jf.setMaximized(1);
				pause(1);

				colors = {'b', 'g', 'r'};
				lineStyles = {'-', '-.', '--'};
				names = {'Absent', '2 cpd', '10 cpd'};
				for(iCond = 1 : size(conditions,2))
					iEcc = find( conditions(iCond).eccentricity == Eccs );
					iSF = find( conditions(iCond).sf == [2 10] );
					iDur = find( conditions(iCond).duration == [500] );
					iPre = find( conditions(iCond).present == [0 1] );

					subplot( nRows, nCols, (iEcc-1)*nCols + 1 ); hold on;
					k = iSF + iPre - 1;
					% fr = FR{iEcc,iSF,iDur,iPre}(:,:,iL);
					tmpFR = FR{iCond,iL}(:,:,:);
					if(applyRectify)
						tmpFR(tmpFR<0) = 0;
					end

					tMin = max( round( time{iCond}(1)/1000*[trials(idx).sRate] ) + [trials(idx).(conditions(iCond).alignEvent)] - [trials(idx).(alignEvent)] );
					tMax = min( round( time{iCond}(end)/1000*[trials(idx).sRate] ) + [trials(idx).(conditions(iCond).alignEvent)] - [trials(idx).(alignEvent)] );
					fr = zeros( size(tmpFR,1), tMax-tMin+1, size(tmpFR,3) );
					for( iTrial = 1 : nTrials )
						s = find( time{iCond} >= round( ( trials(idx(iTrial)).(alignEvent) - trials(idx(iTrial)).(conditions(iCond).alignEvent) ) / trials((idx(iTrial))).sRate * 1000 ), 1, 'first' );
						fr(:,:,iTrial) = tmpFR(:,s+(tMin:tMax),iTrial);
					end

                    fr = mean( fr, 3 );
					m = mean( fr, 1 );
					sem = std( fr, [], 1 ) / sqrt(size(fr,1));

					t = (tMin:tMax) / trials(idx(1)).sRate * 1000;
					h(k) = plot( t, m, 'color', colors{k}, 'LineStyle', lineStyles{k}, 'lineWidth', 2, 'displayName', names{k} );
					fill( [t, t(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], colors{k}, 'LineStyle', 'none', 'FaceAlpha', 0.2 );
					if( k == 3 )
						plot( [0 0], get(gca, 'ylim'), 'k--', 'lineWidth', 2 );
						h(4) = plot( [1 1] * mean([trials(idx).saccadeOn] - [trials(idx).(alignEvent)]), get(gca, 'ylim'), 'b--', 'lineWidth', 2, 'displayName', 'Mean SacOn' );
						h(5) = plot( [1 1] * mean([trials(idx).saccadeLand] - [trials(idx).(alignEvent)]), get(gca, 'ylim'), 'm--', 'lineWidth', 2, 'displayName', 'Mean SacLand' );
						% h(6) = plot( [1 1] * mean([trials(idx).saccadeOff] - [trials(idx).(alignEvent)]), get(gca, 'ylim'), 'r--', 'lineWidth', 2, 'displayName', 'Mean SacOff' );
					end

					if( iEcc == nRows && k == 3 )
						legend( h, 'location', 'northwest' );
					end
					set( gca, 'xlim', [-150 700], 'fontsize', 18, 'lineWidth', 2 );
					if( k == 3 )
						title( sprintf('Ecc=%d | Dur=%d', conditions(iCond).eccentricity, conditions(iCond).duration) );
						ylabel( sprintf('Firing rate') );
					end
					if( iEcc == nRows && k == 3 )
						xlabel( sprintf('Time aligned to %s (ms)', alignEvent ) );
					end
				end

				drawnow;
			end


			figure( 'NumberTitle', 'off', 'name', 'Eye Movements & Cell Properties', 'color', 'w' );
			pause(0.1);
			jf = get(handle(gcf),'javaframe');
			jf.setMaximized(1);
			pause(1);

			% saccadeOn - flashOn, saccadeOff - flashOn
			subplot(3, 2, 1); hold on;
			centers = -100:5:100;
			data = hist( [trials(idx).saccadeOn] - [trials(idx).flashOn], centers );
			bar( centers, data, 0.9, 'b', 'LineStyle', 'none', 'displayName', 'sacOn - flashOn' );
			data = hist( [trials(idx).saccadeOff] - [trials(idx).flashOn], centers );
			bar( centers, data, 0.9, 'r', 'LineStyle', 'none', 'displayName', 'sacOff - flashOn' );
			data = hist( [trials(idx).saccadeLand] - [trials(idx).flashOn], centers );
			bar( centers, data, 0.9, 'm', 'LineStyle', 'none', 'displayName', 'sacLand - flashOn' );
			set( legend, 'location', 'northwest' );
			xlabel('Time aligned to flashOn');
			ylabel('Frequency');
			set( gca, 'lineWidth', 2, 'fontsize', 16 );

			% eye traces
			subplot(3, 2, 2); hold on;
			for( k = size(idx,2) : -1 : 1 )
				x = trials(idx(k)).x.position( -100+trials(idx(k)).saccadeLand : min(end, 100+trials(idx(k)).saccadeLand) );
				eyeX( k, 1:size(x,2) ) = x;
			end
			m = mean( eyeX, 1 );
			sd = std( eyeX, [], 1 );
			plot( -100:100, m, 'b', 'lineWidth', 2 );
			fill( [-100:100, 100:-1:-100], [m-sd, m(end:-1:1)+sd(end:-1:1)], 'b', 'LineStyle', 'none', 'FaceAlpha', 0.5 );
			plot( [0 0], get(gca,'ylim'), 'k--', 'lineWidth', 2 );
			xlabel('Time aligned to saccadeLand');
			ylabel('Horizontal eye (mean+-std) (arcmin)');
			set( gca, 'lineWidth', 2, 'fontsize', 16 );

			% cell spatial sensitivity function
			colors = {[1 0 0], [0 0 1], 'm', 'c'};
			sf = 0.1:0.1:40;
			for(iEcc = 1:size(Eccs,2))
				subplot(3, size(Eccs,2), size(Eccs,2)+iEcc); hold on; h = [];
				for(iL = 1 : 4)
					nCells = size(FR{iEcc,iL},1);
					idx = cellIdx{iL,iEcc}(1:nCells);
					m = mean( abs( obj.SpatialModel.SpatialSensitivity( obj.layers(iL).sRFParams(idx), sf ) ), 1 );
					sem = std( abs( obj.SpatialModel.SpatialSensitivity( obj.layers(iL).sRFParams(idx), sf ) ), [], 1 ) / sqrt(nCells);
					h(iL) = plot( sf, m, 'color', colors{iL}, 'lineWidth', 2, 'displayName', obj.layers(iL).name );
					fill( [sf, sf(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], 'k', 'FaceColor', colors{iL}, 'LineStyle', 'none', 'FaceAlpha', 0.5 );
				end
				set( gca, 'xscale', 'log', 'yscale', 'log', 'lineWidth', 2, 'fontsize', 16 );
				plot( [2 2], get(gca,'ylim'), 'k--', 'lineWidth', 2 );
				plot( [10 10], get(gca,'ylim'), 'k--', 'lineWidth', 2 );
				if(iEcc==nRows), legend(h); end
				xlabel('Spatial frequency (CPD)');
				if(iEcc==1), ylabel('Sensitivity'); end
				title( sprintf('Ecc = %d', conditions(iEcc).eccentricity) );
			end

			% cell temporal sensitivity function
			tf = 0.1:0.1:60;
			for(iEcc = 1:size(Eccs,2))
				subplot(3, size(Eccs,2), size(Eccs,2)*2+iEcc); hold on; h = [];
				for(iL = 1 : 4)
					nCells = size(FR{iEcc,iL},1);
					idx = cellIdx{iL,iEcc}(1:nCells);
					if( iL <= 2 )
						[C, S] = obj.TemporalModel.TemporalSensitivity( obj.layers(iL).name, obj.layers(iL).tRFParams(idx), tf );
						m = mean(C,1);
						sem = std(C, [], 1) / sqrt(nCells);
						h(end+1) = plot( tf, m, 'color', colors{iL}, 'lineWidth', 2, 'displayName', [obj.layers(iL).name ' Center'] );
						fill( [tf, tf(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], 'k', 'FaceColor', colors{iL}, 'LineStyle', 'none', 'FaceAlpha', 0.5 );

						m = mean(S,1);
						sem = std(S, [], 1) / sqrt(nCells);
						h(end+1) = plot( tf, m, 'color', colors{iL}/2, 'lineWidth', 2, 'displayName', [obj.layers(iL).name ' Surround'] );
						fill( [tf, tf(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], 'k', 'FaceColor', colors{iL}, 'LineStyle', 'none', 'FaceAlpha', 0.5 );
					else
						S = obj.TemporalModel.TemporalSensitivity( obj.layers(iL).name, obj.layers(iL).tRFParams(idx), tf, 1 );
						m = mean(S,1);
						sem = std(S, [], 1) / sqrt(nCells);
						h(end+1) = plot( tf, m, 'color', colors{iL}, 'lineWidth', 2, 'displayName', obj.layers(iL).name );
						fill( [tf, tf(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], 'k', 'FaceColor', colors{iL}, 'LineStyle', 'none', 'FaceAlpha', 0.5 );
					end
				end
				if(iEcc==nRows), legend(h); end
				xlabel('Temporal frequency (Hz)');
				if(iEcc == 1), ylabel('Sensitivity'); end
				set( gca, 'xscale', 'log', 'yscale', 'log', 'lineWidth', 2, 'fontsize', 16 );
			end
		end


		function fr = ExampleCells4TrialPredict(obj, FR, time, conditions, trials, trialsIdx, cellIdx, applyRectify, alignEvent, LBOffset, UBOffset, saveFolder)
			if( nargin() < 9 || isempty(applyRectify) )
				applyRectify = true;
			end
			if( nargin() < 12 )
				saveFolder = [];
			end

			nTrials = size( FR{1}, 3 );
			idx = trialsIdx{1}(1:nTrials);

			Eccs = unique([conditions.eccentricity]);
			nEccs = size(Eccs,2);

			nRows = size( unique([conditions.eccentricity]), 2 );
			nCols = 4;

			fr = cell(size(conditions,2),5);
			for( iL = 1 : 5 )
				if( iL ~= 5 )
					layerName = obj.layers(iL).name;
				else
					layerName = 'Layers Average';
				end
				figure( 'NumberTitle', 'off', 'name', sprintf('Example Cells : Trial-by-Trial Prediction : %s+[%d,%d] : %s', alignEvent, LBOffset, UBOffset, layerName), 'color', 'w' );
				pause(0.1);
				jf = get(handle(gcf),'javaframe');
				jf.setMaximized(1);
				pause(1);

				colors = {'b', 'g', 'r'};
				lineStyles = {'-', '-.', '--'};
				names = {'Absent', '2 cpd', '10 cpd'};
				for(iCond = 1 : size(conditions,2))
					iEcc = find( conditions(iCond).eccentricity == Eccs );
					iSF = find( conditions(iCond).sf == [2 10] );
					iDur = find( conditions(iCond).duration == [500] );
					iPre = find( conditions(iCond).present == [0 1] );
					k = iSF + iPre - 1;

					if( iL ~= 5 )
						tmpFR = FR{iCond,iL}(:,:,:);
						if(applyRectify)
							tmpFR(tmpFR<0) = 0;
                        end
						
                        fr{iCond,iL} = zeros(size(FR{iCond,iL},1), UBOffset - LBOffset + 1, nTrials);
						for( iTrial = 1 : nTrials )
							s = find( time{iCond} >= round( ( trials(idx(iTrial)).(alignEvent) - trials(idx(iTrial)).(conditions(iCond).alignEvent) ) / trials((idx(iTrial))).sRate * 1000 ), 1, 'first' );
							fr{iCond,iL}(:,:,iTrial) = tmpFR(:,s+(LBOffset:UBOffset),iTrial);
						end
					else
						fr{iCond,iL}(:,:,:) = cat( 1, fr{iCond,1}(:,:,:), fr{iCond,2}(:,:,:), fr{iCond,3}(:,:,:), fr{iCond,4}(:,:,:) );
					end

					%% Mean (accumulated) firing rate across cell population for each trial
					subplot( nRows, nCols, (iEcc-1)*nCols + 1 ); hold on;
					h = plot( squeeze( mean( sum(fr{iCond,iL}(:,:,:), 2), 1 ) ), 'lineStyle', lineStyles{k}, 'lineWidth', 2, 'color', colors{k}, 'displayName', names{k} );
					% if( iEcc == 1 )
					% 	hs(k) = h;
					% end
					
					if( k == 3 )
						set( gca, 'fontsize', 18, 'lineWidth', 2 );
						title( sprintf('Ecc=%d', conditions(iCond).eccentricity) );
						if( iEcc == 2 )
							ylabel( sprintf('Population mean FR') );
						end
						% if( iEcc == 1 )
						% 	legend( hs, 'location', 'northwest' );
						% end
						if( iEcc == nRows )
							xlabel('Trial index');
						end
					end

					%% Prabability density of population mean (accumulated) firing rate across trials
					subplot( nRows, nCols, (iEcc-1)*nCols + 2 ); hold on;
					points = squeeze( mean( sum(fr{iCond,iL}(:,:,:), 2), 1 ) )';
					st = max(points) / (round(nTrials/3));
					edges = -st/2+min(points) : st : max(points)+st+st/2;
					cnts = histc( points, edges ); cnts(end) = [];
					cnts = cnts / sum(cnts) / (edges(2)-edges(1));
					h = plot( (edges(1:end-1)+edges(2:end))/2, cnts, 'LineStyle', lineStyles{k}, 'color', colors{k}, 'lineWidth', 2, 'displayName', names{k} );
					if( iEcc == 1 )
						hs(k) = h;
					end

					if( k == 3 )
						set( gca, 'fontsize', 18, 'lineWidth', 2 );
						% title( sprintf('Ecc=%d | %s+[%d,%d]', conditions(iCond).eccentricity, alignEvent, LBOffset, UBOffset) );
						if( iEcc == 2 )
							ylabel( sprintf('Probability density') );
						end
						if( iEcc == 1 )
							set( legend( hs, 'location', 'northeast' ), 'position', [0.4627, 0.8528, 0.0552, 0.0884] );
						end
						if( iEcc == nRows )
							xlabel('Population mean FR');
						end
					end

					%% ROC by thresholding
					if( k==2 || k==3 )
						subplot( nRows, nCols, (iEcc-1)*nCols + 3 ); hold on;
						t = [zeros(1,nTrials), ones(1,nTrials)];
						data = cat( 4, fr{iEcc+([1 k]-1)*nEccs, iL} );
						y = reshape( mean( sum(data,2), 1 ), 1, [] );
						y = (y-min(y)) / (max(y)-min(y));
						[tpr, fpr, th] = roc( t, y );
						plot( fpr, tpr, 'lineWidth', 2, 'color', colors{k} );
						text( 1, (4-k)*0.2-0.1, sprintf( 'AUC: %.4f', AreaUnderROC( unique( [tpr;fpr]', 'rows' ) ) ), 'fontsize', 18, 'color', colors{k}, 'horizontalAlignment', 'right', 'verticalAlignment', 'bottom' );

						if( k == 3 )
							set( gca, 'fontsize', 18, 'lineWidth', 2 );
							if( iEcc == 1 )
								title('Thresholding Classifier');
							end
							if( iEcc == 2 )
								ylabel( sprintf('False alarm rate') );
							end
							if( iEcc == nRows )
								xlabel('Correct detection rate');
							end
						end
					end


					%% ROC with GLM classifier
					if( k == 3 )
						subplot( nRows, nCols, (iEcc-1)*nCols + 4 ); hold on;
						data = cat( 4, fr{iEcc+(0:2)*nEccs, iL} );
						X = reshape( mean( data, 2 ), size(data,1), [] )';
						T = [zeros(1,nTrials), ones(1,nTrials*2)]';
						N = round(nTrials/3*2);
						nBoots = 100;
						for(iBoot = nBoots:-1:1)
							rndIdx = [randperm(nTrials), randperm(nTrials)+nTrials, randperm(nTrials)+nTrials*2];
							trainIdx = [1:N, (1:N)+nTrials, (1:N)+2*nTrials];
							mdl = fitglm( X(rndIdx(trainIdx),:), T(rndIdx(trainIdx)), 'linear', 'distribution', 'binomial' );
							for( m = 1 : 2 )
								testIdx = [N+1:nTrials, (N+1:nTrials)+nTrials*m];
								[tpr, fpr] = roc( T(testIdx)', mdl.predict(X(testIdx,:))' );	% True Positive Rate, False Positive Rate
								AUC(m,iBoot) = AreaUnderROC( unique( [tpr;fpr]', 'rows' ) );
							end
						end

						for( m = 1 : 2 )
							M = mean(AUC(m,:));
							sd = std(AUC(m,:));
							plot( m, M, 's', 'color', colors{m+1}, 'markersize', 8, 'lineWidth', 2 );
							plot( [m m], M+[-1 1]*sd, 'lineWidth', 2, 'color', colors{m+1} );
							text( m, 0.1, sprintf( '%.4f\\pm\n%.4f', M, sd ), 'fontsize', 14, 'color', colors{m+1}, 'horizontalAlignment', 'center', 'verticalAlignment', 'bottom' );
						end

						set( gca, 'xlim', [0 3], 'ylim', [0 1.1], 'xtick', [1 2], 'XTickLabel', {'2 cpd', '10 cpd'}, 'fontsize', 18, 'lineWidth', 2 );
						% if( iEcc == 2 )
						% 	ylabel( sprintf('False alarm rate') );
						% end
						if( iEcc == 1 )
							title('GLM Classifier');
						end
						if( iEcc == nRows )
							xlabel( sprintf('Bootstrap AUC (%d)', nBoots) );
						end
					end
				end
				if( ~isempty(saveFolder) && exist(saveFolder,'dir') )
					saveas( gcf, fullfile( saveFolder, sprintf( 'Trial-by-Trial Prediction, %s+[%d,%d], %s.fig', alignEvent, LBOffset, UBOffset, layerName ) ) );
					saveas( gcf, fullfile( saveFolder, sprintf( 'Trial-by-Trial Prediction, %s+[%d,%d], %s.png', alignEvent, LBOffset, UBOffset, layerName ) ) );
				end
			end

			% [coef, scores, latent] = pca( fr );
			% scores = [xVel;yVel]' * coef;
			% plot( [0 coef(1,1)*sqrt(latent(1))], [0 coef(2,1)*sqrt(latent(1))], '-', 'color', colors{iCond}, 'LineWidth', 2 );
			% plot( [0 coef(1,2)*sqrt(latent(2))], [0 coef(2,2)*sqrt(latent(2))], '--', 'color', colors{iCond}, 'LineWidth', 2 );
		end


		function [tTicks, accAUC, segAUC, frAcc, frSeg] = TrialPredictWithAUC(obj, FR, time, conditions, trials, trialsIdx, cellIdx, applyRectify, alignEvent, LBOffset, UBOffset, tStep, tWin, saveFolder)
			if( nargin() < 9 || isempty(applyRectify) )
				applyRectify = true;
			end
			if( nargin() < 12 || isempty(tStep) )
				tStep = 10;
			end
			if( nargin() < 13 || isempty(tWin) )
				tWin = 50;
			end
			if( nargin() < 14 )
				saveFolder = [];
			end
			tWin = round(tWin/2)*2;

			nTrials = size( FR{1}, 3 );
			trialIdx = trialsIdx{1}(1:nTrials);

			Eccs = unique([conditions.eccentricity]);
			nEccs = size(Eccs,2);

			tTicks = LBOffset : tStep : UBOffset;
			accAUC = nan(5,2,nEccs,size(tTicks,2));		% AUC computed with accumulated firing rate in the segment of [LBOffset, tTicks(k)]; layers X SFs X Eccs X time
			segAUC = accAUC;							% AUC computed with firing rate in the segment [-tWin/2,tWin/2]+tTicks(k)

			frAcc = cell(size(conditions,2),5);
			frSeg = frAcc;


			for(iTick = 1 : size(tTicks,2))
				for( iL = 1 : 5 )
					for(iCond = 1 : size(conditions,2))
						iEcc = find( conditions(iCond).eccentricity == Eccs );
						iSF = find( conditions(iCond).sf == [2 10] );
						iDur = find( conditions(iCond).duration == [500] );
						iPre = find( conditions(iCond).present == [0 1] );
						k = iSF + iPre - 1;

						if( iL ~= 5 )
							tmpFR = FR{iCond,iL}(:,:,:);
							if(applyRectify)
								tmpFR(tmpFR<0) = 0;
	                        end
							
	                        frAcc{iCond,iL} = nan(size(FR{iCond,iL},1), tTicks(iTick) - LBOffset + 1, nTrials);
	                        frSeg{iCond,iL} = nan(size(FR{iCond,iL},1), tWin+1, nTrials);
							for( iTrial = 1 : nTrials )
								s = find( time{iCond} >= round( ( trials(trialIdx(iTrial)).(alignEvent) - trials(trialIdx(iTrial)).(conditions(iCond).alignEvent) ) / trials((trialIdx(iTrial))).sRate * 1000 ), 1, 'first' );
								frAcc{iCond,iL}(:,:,iTrial) = tmpFR(:, s+(LBOffset:tTicks(iTick)), iTrial);

								if( s+tTicks(iTick)-tWin/2 < 1 )
									idx = 1 : s+tTicks(iTick)+tWin/2;
									frSeg{iCond,iL}(:,idx,iTrial) = tmpFR(:,idx,iTrial);
								elseif( s+tTicks(iTick)+tWin/2 > size(tmpFR,2) )
									idx = s+tTicks(iTick)-tWin/2 : size(tmpFR,2);
									frSeg{iCond,iL}(:, idx-idx(end)+tWin+1, iTrial) = tmpFR(:,idx,iTrial);
								else
									frSeg{iCond,iL}(:,:,iTrial) = tmpFR(:, s+tTicks(iTick)+(-tWin/2:tWin/2), iTrial);
								end
							end
						else
							frAcc{iCond,iL} = cat( 1, frAcc{iCond,1}(:,:,:), frAcc{iCond,2}(:,:,:), frAcc{iCond,3}(:,:,:), frAcc{iCond,4}(:,:,:) );
							frSeg{iCond,iL} = cat( 1, frSeg{iCond,1}(:,:,:), frSeg{iCond,2}(:,:,:), frSeg{iCond,3}(:,:,:), frSeg{iCond,4}(:,:,:) );
						end

						if( k > 1 )
							t = [zeros(1,nTrials), ones(1,nTrials)];
							data = cat( 4, frAcc{iEcc+([1 k]-1)*nEccs, iL} );
							y = reshape( mean( sum(data,2), 1 ), 1, [] );
							y = (y-min(y)) / (max(y)-min(y));
							[tpr, fpr, th] = roc( t, y );
							accAUC(iL,k-1,iEcc,iTick) = AreaUnderROC( unique( [tpr;fpr]', 'rows' ) );
							
							data = cat( 4, frSeg{iEcc+([1 k]-1)*nEccs, iL} );
							y = reshape( mean( sum(data,2), 1 ), 1, [] );
							y = (y-min(y)) / (max(y)-min(y));
							[tpr, fpr, th] = roc( t, y );
							segAUC(iL,k-1,iEcc,iTick) = AreaUnderROC( unique( [tpr;fpr]', 'rows' ) );
						end
					end
				end
			end


			%% plot figures
			nRows = size( unique([conditions.eccentricity]), 2 );
			nCols = 4;
			colors = {[0 0 1], [0 1 0], [1 0 0]};
			lineStyles = {'-', '-.', '--'};
			names = {'Absent', '2 cpd', '10 cpd'};
			
			for( iL = 1 : 5 )
				if( iL ~= 5 )
					layerName = obj.layers(iL).name;
				else
					layerName = 'Layers Average';
				end
				figure( 'NumberTitle', 'off', 'name', sprintf('Example Cells : Trial-by-Trial Prediction : %s+[%d,%d] : %s', alignEvent, LBOffset, UBOffset, layerName), 'color', 'w' );
				pause(0.1);
				jf = get(handle(gcf),'javaframe');
				jf.setMaximized(1);
				pause(1);

				for(iCond = 1 : size(conditions,2))
					iEcc = find( conditions(iCond).eccentricity == Eccs );
					iSF = find( conditions(iCond).sf == [2 10] );
					iDur = find( conditions(iCond).duration == [500] );
					iPre = find( conditions(iCond).present == [0 1] );
					k = iSF + iPre - 1;

					%% Prabability density of population mean (accumulated) firing rate across trials
					subplot( nRows, nCols, (iEcc-1)*nCols + 1 ); hold on;
					points = squeeze( mean( sum(frAcc{iCond,iL}(:,:,:), 2), 1 ) )';
					st = max(points) / (round(nTrials/3));
					edges = -st/2+min(points) : st : max(points)+st+st/2;
					cnts = histc( points, edges ); cnts(end) = [];
					cnts = cnts / sum(cnts) / (edges(2)-edges(1));
					h = plot( (edges(1:end-1)+edges(2:end))/2, cnts, 'LineStyle', lineStyles{k}, 'color', colors{k}, 'lineWidth', 2, 'displayName', names{k} );
					if( iEcc == 1 )
						hs(k) = h;
						legend( hs, 'location', 'northeast' );
					end

					if( k == 3 )
						set( gca, 'fontsize', 18, 'lineWidth', 2 );
						title( sprintf('Ecc=%d | %s+[%d,%d]', conditions(iCond).eccentricity, alignEvent, LBOffset, UBOffset) );
						if( iEcc == 2 )
							ylabel( sprintf('Probability density') );
						end
						% if( iEcc == 1 )
						% 	set( legend( hs, 'location', 'northeast' ), 'position', [0.4627, 0.8528, 0.0552, 0.0884] );
						% end
						if( iEcc == nRows )
							xlabel('Population mean accumulated FR');
						end
					end

					%% ROC by thresholding
					if( k==2 || k==3 )
						subplot( nRows, nCols, (iEcc-1)*nCols + 2 ); hold on;
						t = [zeros(1,nTrials), ones(1,nTrials)];
						data = cat( 4, frAcc{iEcc+([1 k]-1)*nEccs, iL} );
						y = reshape( mean( sum(data,2), 1 ), 1, [] );
						y = (y-min(y)) / (max(y)-min(y));
						[tpr, fpr, th] = roc( t, y );
						plot( fpr, tpr, 'lineWidth', 2, 'color', colors{k} );
						% text( 1, (4-k)*0.2-0.1, sprintf( 'AUC: %.4f', AreaUnderROC( unique( [tpr;fpr]', 'rows' ) ) ), 'fontsize', 18, 'color', colors{k}, 'horizontalAlignment', 'right', 'verticalAlignment', 'bottom' );
						text( 1, (4-k)*0.2-0.1, sprintf( 'AUC: %.4f', accAUC(iL,k-1,iEcc,end) ), 'fontsize', 18, 'color', colors{k}, 'horizontalAlignment', 'right', 'verticalAlignment', 'bottom' );

						if( k == 3 )
							set( gca, 'fontsize', 18, 'lineWidth', 2 );
							if( iEcc == 1 )
								title('Thresholding Classifier');
							end
							if( iEcc == 2 )
								ylabel( sprintf('False alarm rate') );
							end
							if( iEcc == nRows )
								xlabel('Correct detection rate');
							end
						end
					end

					%% AUC as a function of time computed with accumulated firing rate in the segment of [LBOffset, tTicks(k)]
					if( k > 1 )
						subplot( 2, 2, 2 ); hold on;
						plot( tTicks, squeeze(accAUC(iL,k-1,iEcc,:)), 'color', colors{k}*(nEccs-iEcc+1)/nEccs, 'lineWidth', 2, 'displayName', sprintf('SF=%d, Ecc=%d', conditions(iCond).sf, conditions(iCond).eccentricity) );
						if( k == 3 )
							set( gca, 'fontsize', 18, 'lineWidth', 2 );
							% set( legend, 'location', 'northwest' );
							title('AUC Evolving over Time');
							ylabel('AUC');
						end
					end

					%% AUC as a function of time computed with firing rate in the segment [-tWin/2,tWin/2]+tTicks(k)
					if( k > 1 )
						subplot( 2, 2, 4 ); hold on;
						plot( tTicks, squeeze(segAUC(iL,k-1,iEcc,:)), 'color', colors{k}*(nEccs-iEcc+1)/nEccs, 'lineWidth', 2, 'displayName', sprintf('SF=%d, Ecc=%d', conditions(iCond).sf, conditions(iCond).eccentricity) );
						if( k == 3 )
							set( gca, 'fontsize', 18, 'lineWidth', 2 );
							set( legend, 'location', 'northwest', 'position', [0.8998, 0.3585, 0.0979, 0.1732] );
							title(sprintf('AUC in Sliding Window of %d ms',tWin));
							xlabel(['Time aligned to ', alignEvent]);
							ylabel('AUC');
						end
					end

					%% ROC with GLM classifier
					% if( k == 3 )
					% 	subplot( nRows, nCols, (iEcc-1)*nCols + 4 ); hold on;
					% 	data = cat( 4, fr{iEcc+(0:2)*nEccs, iL} );
					% 	X = reshape( mean( data, 2 ), size(data,1), [] )';
					% 	T = [zeros(1,nTrials), ones(1,nTrials*2)]';
					% 	N = round(nTrials/3*2);
					% 	nBoots = 100;
					% 	for(iBoot = nBoots:-1:1)
					% 		rndIdx = [randperm(nTrials), randperm(nTrials)+nTrials, randperm(nTrials)+nTrials*2];
					% 		trainIdx = [1:N, (1:N)+nTrials, (1:N)+2*nTrials];
					% 		mdl = fitglm( X(rndIdx(trainIdx),:), T(rndIdx(trainIdx)), 'linear', 'distribution', 'binomial' );
					% 		for( m = 1 : 2 )
					% 			testIdx = [N+1:nTrials, (N+1:nTrials)+nTrials*m];
					% 			[tpr, fpr] = roc( T(testIdx)', mdl.predict(X(testIdx,:))' );	% True Positive Rate, False Positive Rate
					% 			AUC(m,iBoot) = AreaUnderROC( unique( [tpr;fpr]', 'rows' ) );
					% 		end
					% 	end

					% 	for( m = 1 : 2 )
					% 		M = mean(AUC(m,:));
					% 		sd = std(AUC(m,:));
					% 		plot( m, M, 's', 'color', colors{m+1}, 'markersize', 8, 'lineWidth', 2 );
					% 		plot( [m m], M+[-1 1]*sd, 'lineWidth', 2, 'color', colors{m+1} );
					% 		text( m, 0.1, sprintf( '%.4f\\pm\n%.4f', M, sd ), 'fontsize', 14, 'color', colors{m+1}, 'horizontalAlignment', 'center', 'verticalAlignment', 'bottom' );
					% 	end

					% 	set( gca, 'xlim', [0 3], 'ylim', [0 1.1], 'xtick', [1 2], 'XTickLabel', {'2 cpd', '10 cpd'}, 'fontsize', 18, 'lineWidth', 2 );
					% 	% if( iEcc == 2 )
					% 	% 	ylabel( sprintf('False alarm rate') );
					% 	% end
					% 	if( iEcc == 1 )
					% 		title('GLM Classifier');
					% 	end
					% 	if( iEcc == nRows )
					% 		xlabel( sprintf('Bootstrap AUC (%d)', nBoots) );
					% 	end
					% end


				end

				if( ~isempty(saveFolder) && exist(saveFolder,'dir') )
					saveas( gcf, fullfile( saveFolder, sprintf( 'Time Course of Trial-by-Trial Prediction, %s+[%d,%d], %s.fig', alignEvent, LBOffset, UBOffset, layerName ) ) );
					saveas( gcf, fullfile( saveFolder, sprintf( 'Time Course of Trial-by-Trial Prediction, %s+[%d,%d], %s.png', alignEvent, LBOffset, UBOffset, layerName ) ) );
				end
			end

			% [coef, scores, latent] = pca( fr );
			% scores = [xVel;yVel]' * coef;
			% plot( [0 coef(1,1)*sqrt(latent(1))], [0 coef(2,1)*sqrt(latent(1))], '-', 'color', colors{iCond}, 'LineWidth', 2 );
			% plot( [0 coef(1,2)*sqrt(latent(2))], [0 coef(2,2)*sqrt(latent(2))], '--', 'color', colors{iCond}, 'LineWidth', 2 );
		end


		function trials = LoadExpData(obj, dataFolder)

			% load experiment data
			load( fullfile(dataFolder, 'A014.mat'), 'ppt', 'eminfo', 'counter' );
			trials = [ppt{:}];
			eminfo.SaccadeStart = eminfo.SaccadeStart( eminfo.DriftOnly' & [trials.contrast] <= 0.5 );
			trials = trials( eminfo.DriftOnly' & [trials.contrast] <= 0.5 );

			% set flashOn to the middle of saccade
			for( iTrial = 1 : size(trials,2) )
				trials(iTrial).flashOn = find( trials(iTrial).x.position(eminfo.SaccadeStart(iTrial):end) <= 200, 1, 'first' ) + eminfo.SaccadeStart(iTrial)-1;		% in samples
			end
			
			% redo saccade detection
			trials = SaccadeTool.GetSacs( trials, 'minmsa', 3 );
			for( iTrial = size(trials,2) : -1 : 1 )
				[~, saccadeIdx(iTrial)] = min(abs( trials(iTrial).saccades.start - trials(iTrial).flashOn /trials(iTrial).sRate*1000 ));
				saccadeOn(iTrial) = trials(iTrial).saccades.start(saccadeIdx(iTrial));																		% in samples
				saccadeOff(iTrial) = trials(iTrial).saccades.start(saccadeIdx(iTrial)) + trials(iTrial).saccades.duration(saccadeIdx(iTrial)) - 1;			% in samples
			end
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


		function [img, inputX, inputY] = LoadNoise(obj, noiseFile, degPerPix)
			%% load noise background from the file specified by noiseFile
			%   noiseFile:			full file path for the noise
			%	degPerPix:			number of degrees per pixel
			%
			%   img:				image pixel values loaded from the file (range of [-1 1]); 1st dimension for vertical, 2nd for horizontal
			%	inputX:				horizontal coordinates of each pixel in img (deg)
			%	inputY:				vertical coordinates of each pixel in img (deg)

			f = fopen(noiseFile);
			w = 1920; h = 1080;
			img = zeros(w,h);
			img(:) = fread(f, w*h, 'float32');
			img = img(:,end:-1:1)'*2 - 1;
			inputX = ( (1:w) - (1+w)/2 ) * degPerPix;
			inputY = ( (1:h) - (1+h)/2 ) * degPerPix;
            fclose(f);
		end


		function [grating, inputX, inputY] = GenerateGrating(obj, radius, width, sf, phase, degPerPix)
			%% generate a circular grating
			%   radius:				radius of the grating (in the middle)
			%	width:				width of the circle (deg)
			%	sf:					spatial frequency (c/deg)
			%	phase:				phase of the grating (0/1)
			%	degPerPix:			number of degrees per pixel
			%
			%   grating:			pixel values of the grating (range of [-1 1]); 1st dimension for vertical, 2nd for horizontal
			%	inputX:				horizontal coordinates of each pixel in img (deg)
			%	inputY:				vertical coordinates of each pixel in img (deg)

			w = 1920; h = 1080;
			grating = zeros(h,w);
			inputX = ( (1:w) - (1+w)/2 ) * degPerPix;
			inputY = ( (1:h) - (1+h)/2 ) * degPerPix;

			if(radius == 0)
				width = width*2;
			end

			r = sqrt( inputX.^2 + inputY'.^2 );
			idx = abs(r - radius) < width/2;
			grating(idx) = cos( pi * (r(idx)-radius) / width ) .* cos( 2*pi * sf * (r(idx)-radius) + phase*pi );
		end


		function contrasts = ComputeContrasts(obj, img, inputX, inputY, rfX, rfY, rfSurroundRadius, eyeX, eyeY)
			% contrasts:	1st dim for cells, 2nd dim for eye positions

			r = 10 * rfSurroundRadius / sqrt(2);		% 6*std of surround
			nNeurons = length(rfX);
			nEyePos = length(eyeX);
			MAXs = NaN(nNeurons,nEyePos);
			MINs = NaN(nNeurons,nEyePos);
			
			parfor( k = 0 : nNeurons * nEyePos-1 )
		        iEyePos = floor(k/nNeurons) + 1;
		        iNeuron = mod(k,nNeurons) + 1;

		        % square area
		        dX = inputX' - (rfX(iNeuron) + eyeX(iEyePos));
		        dY = inputY - (rfY(iNeuron) + eyeY(iEyePos));
		        idxX = abs(dX) <= r(iNeuron);
		        idxY = abs(dY) <= r(iNeuron);
		        MAX = max( max( img(idxY, idxX) ) );
		        MIN = min( min( img(idxY, idxX) ) );
		        if( ~isempty(MAX) )
		        	MAXs(k+1) = MAX;
		        	MINs(k+1) = MIN;
		        end
		    end

		    contrasts = (MAXs - MINs) ./ (MAXs + MINs);
		    contrasts(isnan(contrasts)) = 0;
		end


		function DisplayMosaics(obj, radius)
			% figure( 'NumberTitle', 'off', 'name', 'Mosaics', 'color', 'k' );
			colors = {'r', 'b', 'm', 'c'};
			for( k = 1 : size(obj.layers,2) )
				figure( 'NumberTitle', 'off', 'name', sprintf( 'Mosaics | %s | [%.1f, %.1f]', obj.layers(k).name, -radius, radius ), 'color', 'k' );
				% subplot(2,2,k);
				plot( obj.layers(k).locations(:,1), obj.layers(k).locations(:,2), '.', 'color', colors{k}, 'markersize', 4 );
				% title( obj.layers(k).name );
				text( 0, radius*1.02, obj.layers(k).name, 'horizontalAlignment', 'center', 'verticalAlignment', 'bottom', 'fontsize', 20, 'color', 'w' );
				axis equal;
				set( gca, 'xlim', [-radius radius], 'ylim', [-radius radius], 'fontsize', 20, 'visible', 'off', 'color', 'k', 'xcolor', 'w', 'ycolor', 'w' );
			end
		end
	end
end