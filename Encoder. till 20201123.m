classdef Encoder < handle

	properties (SetAccess = public)%private)
		layers;%(4) = struct( 'name', [], 'locations', [], 'sRFRc', [], 'sRFKc', [], 'sRFRs', [], 'sRFKs', [] );
		SpatialModel;
		TemporalModel;
	end


	methods
		
		function obj = Encoder(paramFolder, HsBe1)
			%% Constructor function
			%	paramFolder:		folder containing cell parameters

			if( nargin() < 1 || isempty(paramFolder) )
				paramFolder = './Parameters/';
			end
			if( nargin() < 2 || isempty(HsBe1) )
				HsBe1 = false;
			end

			obj.TemporalModel = TemporalRF();
			
			obj.layers = struct( 'name', cell(1,4), 'locations', [], 'sRFParams', [], 'tRFParams', [] );
			names = {'POn', 'POff', 'MOn', 'MOff'};
			[obj.layers.name] = names{:};

			% obj.SpatialModel = CronerKaplanRGCModel( 'dataSetToFit', 'raw', 'fitIntercept', false, 'randomizeCenterRadii', true, 'randomizeCenterSensitivities', true, 'randomizeSurroundRadii', true, 'randomizeSurroundSensitivities', true );
			obj.SpatialModel = CronerKaplanRGCModel( 'dataSetToFit', 'raw', 'fitIntercept', false, 'randomizeCenterRadii', false, 'randomizeCenterSensitivities', false, 'randomizeSurroundRadii', false, 'randomizeSurroundSensitivities', false );

			for( k =  1 : size(obj.layers,2) )
				fprintf( 'k = %d\ncellType = %s\n\n', k, names{k} );

				% get cell locations
				if( exist( fullfile( paramFolder, ['Mosaic_' names{k} '_Radius20.0deg_maxMovPrctile20.mat'] ), 'file' ) )
					mosaic_ = load( fullfile( paramFolder, ['Mosaic_' names{k} '_Radius20.0deg_maxMovPrctile20.mat'] ), 'rfPositions', 'spacing' );
					obj.layers(k).locations = mosaic_.rfPositions;
					obj.layers(k).spacing = mosaic_.spacing;
				else
					obj.layers(k).locations = WatsonRGCModel.CreateMosaic(20, names{k}, 'saveFolder', paramFolder, 'maxIterations', 3000);
					obj.layers(k).spacing = WatsonRGCModel.AveSpacingInMosaic(obj.layers(k).locations);
				end

				% get spatial RF parameters
				if( exist( fullfile( paramFolder, ['SpatialRFParams_', names{k}, '.mat'] ), 'file' ) )
					obj.layers(k).sRFParams = load( fullfile( paramFolder, ['SpatialRFParams_', names{k}, '.mat'] ), 'sRFParams' );
					obj.layers(k).sRFParams = obj.layers(k).sRFParams.sRFParams;
				else
					% sRFParams = obj.SpatialModel.SynthesizeRFParams( sqrt(sum(obj.layers(k).locations.^2,2))', names{k} );
					% compute temporal equivalent eccentricity
					ecc_ = 0:0.001:120;
					[~, spacing_] = WatsonRGCModel.RFSpacingDensityMeridian( ecc_, WatsonRGCModel.enumeratedMeridianNames{1}, [cellType OnOff] );   % temporal meridian
					[~, spacing] = WatsonRGCModel.RFSpacingDensity(obj.layers(k).locations, obj.layers(k).name);
					temporalEccDegs = interp1( spacing_, ecc_, spacing, 'linear', 'extrap' );
					
					sRFParams = obj.SpatialModel.Spacing2RFParams(names{k}, obj.layers(k).spacing', temporalEccDegs);
					obj.layers(k).sRFParams = sRFParams;
					save( fullfile( paramFolder, ['SpatialRFParams_', names{k}, '.mat'] ), 'sRFParams' );
				end

				% get temporal RF parameters
				obj.layers(k).tRFParams = obj.TemporalModel.SynthesizeParams( names{k}, size(obj.layers(k).locations, 1), true, '97', HsBe1 );
			end
		end


		function [LFR, time, conditions, sFR, sFR_c, sFR_s, trials, trialsIdx, cellIdx, nCells, nCellsUsed] = ExampleCells(obj, dataFolder, alignEvent, stabilize, saveFolder)
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
			if( nargin() < 4 || isempty(stabilize) )
				stabilize = 'normal';
			end
			if( nargin() < 5 || isempty(saveFolder) )
				saveFolder = sprintf( 'Example Cells LFR' ); %sprintf( 'Example Cells LFR; contrast = %f', contrast );
			end
			if( saveFolder(end) == '/' || saveFolder(end) == '\' )
				saveFolder(end) = [];
			end
			if( ~exist(fullfile(dataFolder,'figures',saveFolder), 'dir') )
				mkdir( fullfile(dataFolder,'figures',saveFolder) );
			end

			trials = EmpiricalBox.LoadSingleData(fullfile(dataFolder, 'A014.mat'));
			
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
            nTrials = 70;%30;


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
							'sf',			{  0,   0,   0,   0,   2,   2,   2,   2,  10,  10,  10,  10}, ...
							'duration',		{500}, ...
							'present',		{  0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1}, ...
							'contrast', 	{0}, ...
							'alignEvent',	{alignEvent} );
			conditions([4:4:end]) = [];
			
			contrasts = [0.01, 0.03, 0.05, 0.07, 0.09, 0.10, 0.20, 0.30, 0.40, 0.50];
			ecc = [conditions.eccentricity];
			ecc = num2cell([ ecc, repmat( ecc(4:end), 1, size(contrasts,2)-1 ) ]);
			sf = [conditions.sf];
			sf = num2cell([ sf, repmat( sf(4:end), 1, size(contrasts,2)-1 ) ]);
			present = num2cell([ 0 0 0, ones(1, size(sf,2)-3) ]);
			contrasts = num2cell([ 0 0 0, reshape( repmat(contrasts, 6, 1), 1, [] ) ]);

			conditions(size(sf,2)) = conditions(1);
			[conditions.eccentricity] = ecc{:};
			[conditions.sf] = sf{:};
			[conditions.present] = present{:};
			[conditions.contrast] = contrasts{:};
			[conditions.duration] = deal(500);
			[conditions.alignEvent] = deal(alignEvent);

			Eccs = unique([conditions.eccentricity]);

			for( iCond = size(conditions,2) : -1 : 51 )
				iEcc = find( conditions(iCond).eccentricity == Eccs );
				contrast = conditions(iCond).contrast;

				trialsIdx{iCond} = find( ... 
										 ...%[trials.eccentricity] == eccs(1,iEcc) & ...
										 ... %[trials.spatialFreq] == SFs(iSF) & ...
										 abs( [trials.stimOff] - [trials.saccadeOff] - conditions(iCond).duration ) < 50 & ...
										 true );%[trials.present] == presents(iPre) );

				idx = trialsIdx{iCond};    %fprintf('nTrials = %d\n', size(idx,2)); continue;
                idx = idx( (1:nTrials) + 30 );% idx( 1 : min(nTrials,end) );
				tMax = round(max( [trials(idx).saccadeOff] + round(conditions(iCond).duration/1000*[trials(idx).sRate]) - [trials(idx).(alignEvent)] ));	% in samples
				tMin = round(min( [trials(idx).saccadeOn] - round(0.150*[trials(idx).sRate]) - [trials(idx).(alignEvent)] ));				% in samples
				for( iL = 4 : -1 : 1 )
					% sFR{iCond,iL} = zeros( nCellsUsed{iL,iEcc}, tMax-tMin+1, size(idx,2) );		% 1st dim: neurons;		2nd dim: time;		3rd dim: trials;
					sFR{iL} = zeros( nCellsUsed{iL,iEcc}, tMax-tMin+1, size(idx,2) );		% 1st dim: neurons;		2nd dim: time;		3rd dim: trials;
				end
				sFR_c = sFR;	%sFR_c{iCond,iL} = sFR{iCond,iL};
				sFR_s = sFR;	%sFR_s{iCond,iL} = sFR{iCond,iL};
				LFR(iCond,:) = sFR;	%LFR{iCond,iL} = sFR{iCond,iL};
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
					avContrast = trials(idx(k)).backgroundContrast;% + contrast * conditions(iCond).present;
					% contrastF = [ zeros(1,s1), ones(1,s2) * contrast ] + trials(idx(k)).backgroundContrast;
					% contrastF = obj.ComputeContrasts( noise+grating+bg,  inputX, inputY, obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),1), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),2), [obj.layers(iL).sRFParams(cellIdx{iL,iEcc}(1:nCells)).surroundRadii], x([eyeIdx1,eyeIdx2]), y([eyeIdx1,eyeIdx2]) );
					% contrastF = 1;%contrast;

					for( iL = 1:4 )
						tic;
						for( m = 1 : 3 )
							xIdx = min(x(eyeIdx{m})) + min(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),1)) - 4 <= inputX & inputX <= max(x(eyeIdx{m})) + max(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),1)) + 4; %720 : 1460;
		                    yIdx = min(y(eyeIdx{m})) + min(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),2)) - 4 <= inputY & inputY <= max(y(eyeIdx{m})) + max(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),2)) + 4; %300 : 780;
	                    
							% [sFR{iCond,iL}(:,eyeIdx{m},k), sFR_c{iCond,iL}(:,eyeIdx{m},k), sFR_s{iCond,iL}(:,eyeIdx{m},k)] = obj.SpatialModel.LinearResponse( stimulus(yIdx,xIdx), inputX(xIdx), inputY(yIdx), x(eyeIdx{m}), y(eyeIdx{m}), obj.layers(iL).sRFParams(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc})), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),1), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),2));
							[sFR{iL}(:,eyeIdx{m},k), sFR_c{iL}(:,eyeIdx{m},k), sFR_s{iL}(:,eyeIdx{m},k)] = obj.SpatialModel.LinearResponse( stimulus(yIdx,xIdx), inputX(xIdx), inputY(yIdx), x(eyeIdx{m}), y(eyeIdx{m}), obj.layers(iL).sRFParams(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc})), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),1), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),2));
						end
						if( lower(obj.layers(iL).name(1)) == 'm' )
							LFR{iCond,iL}(:,:,k) = obj.TemporalModel.LinearResponse( obj.layers(iL).name, obj.layers(iL).tRFParams(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc})), avContrast, trials(idx(k)).sRate, sFR{iL}(:,:,k) );
						else
							LFR{iCond,iL}(:,:,k) = obj.TemporalModel.LinearResponse( obj.layers(iL).name, obj.layers(iL).tRFParams(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc})), avContrast, trials(idx(k)).sRate, sFR_c{iL}(:,:,k), sFR_s{iL}(:,:,k) );
						end
						fprintf('Contrast = %.2f, Present=%d, SF=%d, Ecc=%d, iTrials=%d/%d, iL=%d, t=%f\n', conditions(iCond).contrast, conditions(iCond).present, conditions(iCond).sf, conditions(iCond).eccentricity, k, size(idx,2), iL, toc);
					end
					
				end
				if(isempty(idx))
					time{iCond} = [];
				else
					time{iCond} = (tMin:tMax) / trials(idx(1)).sRate * 1000;
				end

				if(iCond == 51)% || conditions(iCond-1).contrast ~= contrast)
					save( fullfile( dataFolder, 'figures', saveFolder, sprintf('%s-. 51-end.mat', saveFolder) ), 'LFR', 'time', 'conditions', 'trials', 'trialsIdx', 'cellIdx', 'nCells', 'nCellsUsed' );
					LFR = cellfun( @(x) {[]}, LFR );
				end
            end
            
            % save( fullfile( dataFolder, 'figures', saveFolder, [saveFolder, '.mat'] ), 'LFR', 'time', 'conditions', 'sFR', 'sFR_c', 'sFR_s', 'trials', 'trialsIdx', 'cellIdx', 'nCells', 'nCellsUsed' );
            % save( fullfile( dataFolder, 'figures', saveFolder, [saveFolder, '.mat'] ), 'LFR', 'time', 'conditions', 'trials', 'trialsIdx', 'cellIdx', 'nCells', 'nCellsUsed' );
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
			data = hist( [trials(idx).saccadeOn] - [trials(idx).saccadeOff], centers );
			bar( centers, data, 0.9, 'b', 'LineStyle', 'none', 'displayName', 'sacOn - sacOff' );
			% data = hist( [trials(idx).saccadeOff] - [trials(idx).flashOn], centers );
			% bar( centers, data, 0.9, 'r', 'LineStyle', 'none', 'displayName', 'sacOff - flashOn' );
			data = hist( [trials(idx).saccadeLand] - [trials(idx).saccadeOff], centers );
			bar( centers, data, 0.9, 'm', 'LineStyle', 'none', 'displayName', 'sacLand - sacOff' );
			set( legend, 'location', 'northwest' );
			xlabel('Time aligned to saccadeOff');
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
			colors = {[1 0 0], [0 0 1], [1 0 1], [0 1 1]};
			sf = 0.1:0.1:40;
			for(iL = 1 : 4)
				subplot(3, 4, 4+iL); hold on; h = [];
				for(iEcc = 1:size(Eccs,2))
					nCells = size(FR{iEcc,iL},1);
					idx = cellIdx{iL,iEcc}(1:nCells);
					m = mean( abs( obj.SpatialModel.SpatialSensitivity( obj.layers(iL).sRFParams(idx), sf ) ), 1 );
					sem = std( abs( obj.SpatialModel.SpatialSensitivity( obj.layers(iL).sRFParams(idx), sf ) ), [], 1 ) / sqrt(nCells);
					fill( [sf, sf(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], 'k', 'FaceColor', colors{iL}, 'LineStyle', 'none', 'FaceAlpha', 0.5 );
					h(iEcc) = plot( sf, m, 'color', colors{iL}*(size(Eccs,2)-iEcc+1)/size(Eccs,2), 'lineWidth', 2, 'displayName', sprintf('Ecc = %d', conditions(iEcc).eccentricity) );
				end
				set( gca, 'xscale', 'log', 'yscale', 'log', 'lineWidth', 2, 'fontsize', 16 );
				plot( [2 2], get(gca,'ylim'), 'k--', 'lineWidth', 2 );
				plot( [10 10], get(gca,'ylim'), 'k--', 'lineWidth', 2 );
				set(legend(h), 'location', 'southwest');
				xlabel('Spatial frequency (CPD)');
				if(iL==1), ylabel('Sensitivity'); end
				title(obj.layers(iL).name);
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
								ylabel( sprintf('Correct detection rate') );
							end
							if( iEcc == nRows )
								xlabel('False alarm rate');
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
						nBoots = 3;
						for(iBoot = nBoots:-1:1)
							rndIdx = [randperm(nTrials), randperm(nTrials)+nTrials, randperm(nTrials)+nTrials*2];
							trainIdx = [1:N, (1:N)+nTrials, (1:N)+2*nTrials];
							mdl = fitglm( X(rndIdx(trainIdx),:), T(rndIdx(trainIdx)), 'linear', 'distribution', 'binomial' );
							for( m = 1 : 2 )
								testIdx = [N+1:nTrials, (N+1:nTrials)+nTrials*m];
								[tpr, fpr] = roc( T(rndIdx(testIdx))', mdl.predict(X(rndIdx(testIdx),:))' );	% True Positive Rate, False Positive Rate
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
						if( iEcc == 2 )
							ylabel( sprintf('AUC') );
						end
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


		function [tTicks, accAUCm, accAUCsd, segAUCm, segAUCsd, frAcc, frSeg, Weights] = TrialPredictWithAUC(obj, classifier, FR, time, conditions, trials, trialsIdx, nCells, cellIdx, cell2AnalyzeIdx, applyRectify, alignEvent, LBOffset, UBOffset, tStep, tWin, saveFolder)
			if( nargin() < 10 || isempty(applyRectify) )
				applyRectify = true;
			end
			if( nargin() < 13 || isempty(tStep) )
				tStep = 10;
			end
			if( nargin() < 14 || isempty(tWin) )
				tWin = 50;
			end
			if( nargin() < 15 )
				saveFolder = [];
			end
			tWin = round(tWin/2)*2;

			classifiers = {'threshold-Uni', 'thresholding-TW_Uni', 'thresholding-TW_pVal', 'IdealObserver', 'GLM', 'FLDA', 'FLDA_Time', 'FLDA_Single'};
			% classifier = classifiers{1};

			nTrials = size( FR{1}, 3 );
			trialIdx = trialsIdx{1}(1:nTrials);

			Eccs = unique([conditions.eccentricity]);
			nEccs = size(Eccs,2);

			SFs = unique([conditions.sf]);
			nSFs = size(SFs,2);

			tTicks = LBOffset : tStep : UBOffset;
			accAUC = nan(5,2,nEccs,size(tTicks,2));		% AUC computed with accumulated firing rate in the segment of [LBOffset, tTicks(k)]; layers X SFs X Eccs X time
			segAUC = accAUC;							% AUC computed with firing rate in the segment [-tWin/2,tWin/2]+tTicks(k)

			accAUCm = accAUC;
			segAUCm = accAUC;
			accAUCsd = accAUC;
			segAUCsd = accAUC;

			frAcc = cell(size(conditions,2),5);
			frSeg = frAcc;

			%% previously, data were collected with a vertical range of [-0.5, 0.5], here we select those cells within a sector
			for(iCond = 1 : size(FR,1))
				iEcc = find( conditions(iCond).eccentricity == Eccs );
				for(iL = 1 : size(FR,2))
					FR{iCond,iL} = FR{iCond,iL}( cell2AnalyzeIdx{iL,iEcc}, :, : );
				end
			end

			cellNumAmplifier = cell2mat(nCells) ./ cellfun(@(x) size(x,1), cell2AnalyzeIdx);		% inverse of proportion of cells used
			cellNumAmplifier(5,:) = mean(cellNumAmplifier,1);
			areaAmplifier = [1 1 1];	% it's already accounted for in nCells and nCellsUsed	 % [4.800, 25.039, 50.182];	% only an area of +-0.5 deg vertically in right visual field was sampled for each eccentricity


			for(iTick = 1 : size(tTicks,2))
				for( iL = 1 : 5 )
					for(iCond = 1 : size(conditions,2))
						iEcc = find( conditions(iCond).eccentricity == Eccs );
						iSF = find( conditions(iCond).sf == SFs );
						iPre = find( conditions(iCond).present == [0 1] );
						k = iSF + iPre - 1;

						if( iL ~= 5 )
							tmpFR = FR{iCond,iL}(:,:,:);
							if(applyRectify)
								tmpFR(tmpFR<0) = 0;
	                        end
							
	                        frAcc{iCond,iL} = nan(size(FR{iCond,iL},1), tTicks(iTick) - LBOffset + 1, nTrials, 'single');
	                        frSeg{iCond,iL} = nan(size(FR{iCond,iL},1), tWin+1, nTrials, 'single');
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
							for(iiL = 1 : 5)
                                frAcc{iCond,iiL}(:, any(any(isnan(frAcc{iCond,iiL}), 1), 3), :) = [];
                                frSeg{iCond,iiL}(:, any(any(isnan(frSeg{iCond,iiL}), 1), 3), :) = [];
                            end
						end

						switch(lower(classifier))
							case 'thresholding-uni'
								%% ROC with thresholding, uniform weights across space and time
								if( k > 1 )
									T = [zeros(1,nTrials), ones(1,nTrials)];
									data = cat( 3, frAcc{iEcc+([1 k]-1)*nEccs, iL} );	% catenate in the trials dimension
									w = ones(1,size(data,2));		% uniform temporal weight
									w = w / sum(w);
									y = w * shiftdim(nanmean(data,1),1);
									y(T==0) = (y(T==0) - mean(y(T==0))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==0));	% compensate for the low sampling
									y(T==1) = (y(T==1) - mean(y(T==1))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==1));	% compensate for the low sampling
									% if(tTicks(iTick) <= 50)
										thresholdAcc{iCond,iL} = prctile( y(T==0), 90 + max(0,5*tTicks(iTick)/max(tTicks)) );		% set threshold at the level giving 10%~5% false alarm rate
									% end
									accAUCm(iL,k-1,iEcc,iTick) = sum( (y > thresholdAcc{iCond,iL}) == T & T == 1 ) / sum(T==1);	% true positive rate of prediction
									accAUCsd(iL,k-1,iEcc,iTick) = 0;
									
									data = cat( 3, frSeg{iEcc+([1 k]-1)*nEccs, iL} );	% catenate in the trials dimension
									w = ones(1,size(data,2));		% uniform temporal weight
									w = w / sum(w);
									y = w * shiftdim(nanmean(data,1),1);
									y(T==0) = (y(T==0) - mean(y(T==0))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==0));	% compensate for the low sampling
									y(T==1) = (y(T==1) - mean(y(T==1))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==1));	% compensate for the low sampling
									y = (y-min(y)) / (max(y)-min(y));
									% if(tTicks(iTick) <= 50)
										thresholdSeg{iCond,iL} = prctile( y(T==0), 90 + max(0,5*tTicks(iTick)/max(tTicks)) );		% set threshold at the level giving 10%~5% false alarm rate
									% end
									segAUCm(iL,k-1,iEcc,iTick) = sum( (y > thresholdSeg{iCond,iL}) == T & T == 1 ) / sum(T==1);	% true positive rate of prediction
									segAUCsd(iL,k-1,iEcc,iTick) = 0;
								end

							case 'thresholding-sw_sine-tw_uni'
								%% ROC with thresholding, sine wave spatial weights and uniform temporal weights
								% 	For this approach, we combine ON and OFF cells
								if( k > 1 && any(iL == [2,4,5]) )		% 1-POn, 2-POff, 3-MOn, 4-MOff, 5-All
									if(iL ~= 5)
										dataOn = cat( 3, frAcc{iEcc+([1 k]-1)*nEccs, iL-1} );	% catenate in the trials dimension
										dataOff = cat( 3, frAcc{iEcc+([1 k]-1)*nEccs, iL} );	% catenate in the trials dimension
										rfLocOn = obj.layers(iL-1).locations( cellIdx{iL-1,iEcc}(cell2AnalyzeIdx{iL-1,iEcc}), : );
										rfLocOff = obj.layers(iL).locations( cellIdx{iL,iEcc}(cell2AnalyzeIdx{iL,iEcc}), : );
									else
										dataOn = cat( 1, cat( 3, frAcc{iEcc+([1 k]-1)*nEccs, 1} ), cat( 3, frAcc{iEcc+([1 k]-1)*nEccs, 3} ) );
										dataOff = cat( 1, cat( 3, frAcc{iEcc+([1 k]-1)*nEccs, 2} ), cat( 3, frAcc{iEcc+([1 k]-1)*nEccs, 4} ) );
										rfLocOn = [obj.layers(1).locations( cellIdx{1,iEcc}(cell2AnalyzeIdx{1,iEcc}), : ); obj.layers(3).locations( cellIdx{3,iEcc}(cell2AnalyzeIdx{3,iEcc}), : )];
										rfLocOff = [obj.layers(2).locations( cellIdx{2,iEcc}(cell2AnalyzeIdx{2,iEcc}), : ); obj.layers(4).locations( cellIdx{4,iEcc}(cell2AnalyzeIdx{4,iEcc}), : )];
									end
									wOn = cos( 2*pi*SFs(iSF) * (sqrt(sum(rfLocOn.^2,2))' - Eccs(iEcc)) ) / 2 + 0.5;		% column vector
									wOff = -cos( 2*pi*SFs(iSF) * (sqrt(sum(rfLocOff.^2,2))' - Eccs(iEcc)) ) / 2 + 0.5;
									w = [wOn, wOff];
									w = w / sum(w);
									y = w * squeeze(nanmean([dataOn; dataOff],2));

									T = [zeros(1,nTrials), ones(1,nTrials)];
									y(T==0) = (y(T==0) - mean(y(T==0))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==0));	% compensate for the low sampling
									y(T==1) = (y(T==1) - mean(y(T==1))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==1));	% compensate for the low sampling
									% if(tTicks(iTick) <= 50)
										thresholdAcc{iCond,iL} = prctile( y(T==0), 90 );%+ max(0,5*tTicks(iTick)/max(tTicks)) );		% set threshold at the level giving 10%~5% false alarm rate
									% end
									accAUCm(iL,k-1,iEcc,iTick) = sum( (y > thresholdAcc{iCond,iL}) == T & T == 1 ) / sum(T==1);	% true positive rate of prediction
									accAUCsd(iL,k-1,iEcc,iTick) = 0;
									
									if(iL ~= 5)
										dataOn = cat( 3, frSeg{iEcc+([1 k]-1)*nEccs, iL-1} );	% catenate in the trials dimension
										dataOff = cat( 3, frSeg{iEcc+([1 k]-1)*nEccs, iL} );	% catenate in the trials dimension
										rfLocOn = obj.layers(iL-1).locations( cellIdx{iL-1,iEcc}(cell2AnalyzeIdx{iL-1,iEcc}), : );
										rfLocOff = obj.layers(iL).locations( cellIdx{iL,iEcc}(cell2AnalyzeIdx{iL,iEcc}), : );
									else
										dataOn = cat( 1, cat( 3, frSeg{iEcc+([1 k]-1)*nEccs, 1} ), cat( 3, frSeg{iEcc+([1 k]-1)*nEccs, 3} ) );
										dataOff = cat( 1, cat( 3, frSeg{iEcc+([1 k]-1)*nEccs, 2} ), cat( 3, frSeg{iEcc+([1 k]-1)*nEccs, 4} ) );
										rfLocOn = [obj.layers(1).locations( cellIdx{1,iEcc}(cell2AnalyzeIdx{1,iEcc}), : ); obj.layers(3).locations( cellIdx{3,iEcc}(cell2AnalyzeIdx{3,iEcc}), : )];
										rfLocOff = [obj.layers(2).locations( cellIdx{2,iEcc}(cell2AnalyzeIdx{2,iEcc}), : ); obj.layers(4).locations( cellIdx{4,iEcc}(cell2AnalyzeIdx{4,iEcc}), : )];
									end
									wOn = cos( 2*pi*SFs(iSF) * (sqrt(sum(rfLocOn.^2,2))' - Eccs(iEcc)) ) / 2 + 0.5;		% column vector
									wOff = -cos( 2*pi*SFs(iSF) * (sqrt(sum(rfLocOff.^2,2))' - Eccs(iEcc)) ) / 2 + 0.5;
									w = [wOn, wOff];
									w = w / sum(w);
									y = w * squeeze(nanmean([dataOn; dataOff],2));

									T = [zeros(1,nTrials), ones(1,nTrials)];
									y(T==0) = (y(T==0) - mean(y(T==0))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==0));	% compensate for the low sampling
									y(T==1) = (y(T==1) - mean(y(T==1))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==1));	% compensate for the low sampling
									% if(tTicks(iTick) <= 50)
										thresholdSeg{iCond,iL} = prctile( y(T==0), 90 );%+ max(0,5*tTicks(iTick)/max(tTicks)) );		% set threshold at the level giving 10%~5% false alarm rate
									% end
									segAUCm(iL,k-1,iEcc,iTick) = sum( (y > thresholdSeg{iCond,iL}) == T & T == 1 ) / sum(T==1);	% true positive rate of prediction
									segAUCsd(iL,k-1,iEcc,iTick) = 0;
								end

							case 'thresholding-sw_uni-tw_pval'
								%% ROC with thresholding, uniform spatial weights, p-values as temporal weights
								if( k > 1 )
									T = [zeros(1,nTrials), ones(1,nTrials)];
									data = cat( 3, frAcc{iEcc+([1 k]-1)*nEccs, iL} );	% catenate in the trials dimension
									w = zeros(1,size(data,2));		% temporal weight by seperation of two distributions
									for( iT = 1 : size(data,2) )
										r = squeeze( mean(data(:,iT,:), 1) )';	% average across neurons
										if( any(isnan(r)) )
											continue;
										end
										[~,w(iT)] = ttest2( r(T==0), r(T==1) );
										w(iT) = 1 - w(iT);			% smaller p-value corresponds to greater weight
									end
									w = w / sum(w);
									y = w * shiftdim(nanmean(data,1),1);
									y(T==0) = (y(T==0) - mean(y(T==0))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==0));	% compensate for the low sampling
									y(T==1) = (y(T==1) - mean(y(T==1))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==1));	% compensate for the low sampling
									% if(tTicks(iTick) <= 50)
										thresholdAcc{iCond,iL} = prctile( y(T==0), 90 );%+ max(0,5*tTicks(iTick)/max(tTicks)) );		% set threshold at the level giving 10%~5% false alarm rate
									% end
									accAUCm(iL,k-1,iEcc,iTick) = sum( (y > thresholdAcc{iCond,iL}) == T & T == 1 ) / size(y(T==1),2);	% true positive rate of prediction
									accAUCsd(iL,k-1,iEcc,iTick) = 0;
									
									data = cat( 3, frSeg{iEcc+([1 k]-1)*nEccs, iL} );	% catenate in the trials dimension
									w = zeros(1,size(data,2));		% temporal weight by seperation of two distributions
									for( iT = 1 : size(data,2) )
										r = squeeze( mean(data(:,iT,:), 1) )';	% average across neurons
										if( any(isnan(r)) )
											continue;
										end
										[~,w(iT)] = ttest2( r(T==0), r(T==1) );
										w(iT) = 1 - w(iT);			% smaller p-value corresponds to greater weight
									end
									w = w / sum(w);
									y = w * shiftdim(nanmean(data,1),1);
									y(T==0) = (y(T==0) - mean(y(T==0))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==0));	% compensate for the low sampling
									y(T==1) = (y(T==1) - mean(y(T==1))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==1));	% compensate for the low sampling
									% if(tTicks(iTick) <= 50)
										thresholdSeg{iCond,iL} = prctile( y(T==0), 90 );%+ max(0,5*tTicks(iTick)/max(tTicks)) );		% set threshold at the level giving 10%~5% false alarm rate
									% end
									segAUCm(iL,k-1,iEcc,iTick) = sum( (y > thresholdSeg{iCond,iL}) == T & T == 1 ) / size(y(T==1),2);	% true positive rate of prediction
									segAUCsd(iL,k-1,iEcc,iTick) = 0;
								end

							case 'thresholding-sw_uni-tw_dprime'
								%% ROC with thresholding, uniform spatial weights, d-prime as temporal weights
								if( k > 1 )
									T = [zeros(1,nTrials), ones(1,nTrials)];
									data = cat( 3, frAcc{iEcc+([1 k]-1)*nEccs, iL} );	% catenate in the trials dimension
									w = zeros(1,size(data,2));		% temporal weight by seperation of two distributions
									for( iT = 1 : size(data,2) )
										r = squeeze( mean(data(:,iT,:), 1) )';	% average across neurons
										if( any(isnan(r)) )
											continue;
										end
										w(iT) = max(0, (mean(r(T==1)) - mean(r(T==0))) / sqrt(var(r(T==1) + var(r(T==0)))));
									end
									w = w / sum(w);
									y = w * shiftdim(nanmean(data,1),1);
									y(T==0) = (y(T==0) - mean(y(T==0))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==0));	% compensate for the low sampling
									y(T==1) = (y(T==1) - mean(y(T==1))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==1));	% compensate for the low sampling
									% if(tTicks(iTick) <= 50)
										thresholdAcc{iCond,iL} = prctile( y(T==0), 90 );%+ max(0,5*tTicks(iTick)/max(tTicks)) );		% set threshold at the level giving 10%~5% false alarm rate
									% end
									accAUCm(iL,k-1,iEcc,iTick) = sum( (y > thresholdAcc{iCond,iL}) == T & T == 1 ) / size(y(T==1),2);	% true positive rate of prediction
									accAUCsd(iL,k-1,iEcc,iTick) = 0;
									
									data = cat( 3, frSeg{iEcc+([1 k]-1)*nEccs, iL} );	% catenate in the trials dimension
									w = zeros(1,size(data,2));		% temporal weight by seperation of two distributions
									for( iT = 1 : size(data,2) )
										r = squeeze( mean(data(:,iT,:), 1) )';	% average across neurons
										if( any(isnan(r)) )
											continue;
										end
										[~,w(iT)] = ttest2( r(T==0), r(T==1) );
										w(iT) = 1 - w(iT);			% smaller p-value corresponds to greater weight
									end
									w = w / sum(w);
									y = w * shiftdim(nanmean(data,1),1);
									y(T==0) = (y(T==0) - mean(y(T==0))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==0));	% compensate for the low sampling
									y(T==1) = (y(T==1) - mean(y(T==1))) / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc)) + mean(y(T==1));	% compensate for the low sampling
									% if(tTicks(iTick) <= 50)
										thresholdSeg{iCond,iL} = prctile( y(T==0), 90 );%+ max(0,5*tTicks(iTick)/max(tTicks)) );		% set threshold at the level giving 10%~5% false alarm rate
									% end
									segAUCm(iL,k-1,iEcc,iTick) = sum( (y > thresholdSeg{iCond,iL}) == T & T == 1 ) / size(y(T==1),2);	% true positive rate of prediction
									segAUCsd(iL,k-1,iEcc,iTick) = 0;
								end

							case 'idealobserver'
								%% Bayesian inference as ideal observer
								if( k > 1 )
									T = [zeros(1,nTrials), ones(1,nTrials)];
									data = cat( 3, frAcc{iEcc+([1 k]-1)*nEccs, iL} );	% catenate in the trials dimension
									logP = 0;
									for( iT = 1 : size(data,2) )
										r = squeeze( mean(data(:,iT,:), 1) )';	% average across neurons
										if( any(isnan(r)) )
											continue;
										end
										m0 = mean(r(T==0));
										sd0 = std(r(T==0));
										m1 = mean(r(T==1));
										sd1 = std(r(T==1));
										sd0 = sd0 / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc));	% compensate for the low sampling
										sd1 = sd1 / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc));	% compensate for the low sampling
										logP = logP + log(normpdf(r,m1,sd1)) - log(normpdf(r,m0,sd0));
									end
									accAUCm(iL,k-1,iEcc,iTick) = sum( (logP >= 0) == T  ) / size(T,2);	% correct rate
									accAUCsd(iL,k-1,iEcc,iTick) = 0;
									
									data = cat( 3, frSeg{iEcc+([1 k]-1)*nEccs, iL} );	% catenate in the trials dimension
									logP = 0;
									for( iT = 1 : size(data,2) )
										r = squeeze( mean(data(:,iT,:), 1) )';	% average across neurons
										if( any(isnan(r)) )
											continue;
										end
										m0 = mean(r(T==0));
										sd0 = std(r(T==0));
										m1 = mean(r(T==0));
										sd1 = std(r(T==0));
										sd0 = sd0 / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc));	% compensate for the low sampling
										sd1 = sd1 / sqrt(cellNumAmplifier(iL,iEcc) * areaAmplifier(iEcc));	% compensate for the low sampling
										logP = logP + log(normpdf(r,m1,sd1)) - log(normpdf(r,m0,sd0));
									end
									segAUCm(iL,k-1,iEcc,iTick) = sum( (logP >= 0) == T ) / size(T,2);	% correct rate
									segAUCsd(iL,k-1,iEcc,iTick) = 0;
								end

							case 'glm'
								%% ROC with GLM classifier
								if( k == 3 )
									data = {cat(4, frAcc{iEcc+(0:2)*nEccs, iL}), cat(4, frAcc{iEcc+(0:2)*nEccs, iL})};
									for( ii = 2:-1:1 )
										X = reshape( nanmean(data{ii},2), size(data{ii},1), [] )';
										[coef, scores, latent] = pca(X);
										iLat = 0; s = 0; upperLat = sum(latent) * 0.95;	% 99% of variance
										while( s < upperLat )
											iLat = iLat+1;
											s = s + latent(iLat);
										end
										iLat = 4;
										X = scores(:, 1:iLat);
										T = [zeros(1,nTrials), ones(1,nTrials*2)]';
										N = round(nTrials/3*2);
										nBoots = 3;
										for(iBoot = nBoots:-1:1)
											rndIdx = [randperm(nTrials), randperm(nTrials)+nTrials, randperm(nTrials)+nTrials*2];
											trainIdx = [1:N, (1:N)+nTrials, (1:N)+2*nTrials];
											mdl = fitglm( X(rndIdx(trainIdx),:), T(rndIdx(trainIdx)), 'linear', 'distribution', 'binomial' );
											for( m = 2:-1:1 )
												predicts = mdl.predict(X(rndIdx(trainIdx),:));
												[tpr, fpr] = roc( T(rndIdx(trainIdx))', predicts' );	% True Positive Rate, False Positive Rate
												[~,iThresh] = max(tpr-fpr);
												threshold = sort(unique(predicts), 'descend');
												threshold = threshold(iThresh-1);
												testIdx = [N+1:nTrials, (N+1:nTrials)+nTrials*m];
												AUC(m,iBoot) = sum( (mdl.predict(X(rndIdx(testIdx),:)) > threshold) == T(rndIdx(testIdx)) ) / size(testIdx,2);		% prediction correct rate
												% [tpr, fpr] = roc( T(rndIdx(testIdx))', mdl.predict(X(rndIdx(testIdx),:))' );	% True Positive Rate, False Positive Rate
												% AUC(m,iBoot) = AreaUnderROC( unique( [tpr;fpr]', 'rows' ) );
											end
										end
										auc(ii,1,:) = mean(AUC,2);
										auc(ii,2,:) = std(AUC,[],2);
									end
									accAUCm(iL,:,iEcc,iTick) = auc(1,1,:);
									accAUCsd(iL,:,iEcc,iTick) = auc(1,2,:);
									segAUCm(iL,:,iEcc,iTick) = auc(2,1,:);
									segAUCsd(iL,:,iEcc,iTick) = auc(2,2,:);
								end
								
							case 'flda'
								%% ROC with Fisher's Linear Discriminant Analysis
								% each variable/dimenstion corresponds to a cell's mean firing rate across time, each trial corresponds to one sample
								if( k == 3 )
									data = { frAcc{iEcc,iL}, cat(4, frAcc{iEcc+(1:2)*nEccs,iL}); frSeg{iEcc,iL}, cat(4, frSeg{iEcc+(1:2)*nEccs,iL}) };
									for( ii = 2:-1:1 )
										X1 = reshape( nanmean(data{ii,1},2), size(data{ii,1},1), [] );
										X2 = reshape( nanmean(data{ii,2},2), size(data{ii,2},1), [] );
										[coef, scores, latent] = pca( [X1,X2]' );
										iLat = 0; s = 0; upperLat = sum(latent) * 0.80;	% 99% of variance
										while( s < upperLat )
											iLat = iLat+1;
											s = s + latent(iLat);
										end
										X1 = scores(1:nTrials, 1:iLat-1)';
										X2 = scores(nTrials+1:end, 1:iLat-1)';
										T1 = false(1,nTrials);
										T2 = true(1,nTrials*2);
										N = round(nTrials/3*2);
										trainIdx1 = 1:N;
										trainIdx2 = [1:N, nTrials+(1:N)];
										testIdx1 = N+1:nTrials;
										testIdx2 = (0:1)'*nTrials+(N+1:nTrials);
										nBoots = 3;
										for(iBoot = nBoots:-1:1)
											rndIdx1 = randperm(nTrials);
											rndIdx2 = [randperm(nTrials), randperm(nTrials)+nTrials];
											mu1 = mean( X1(:,rndIdx1(trainIdx1)), 2 );
											mu2 = mean( X2(:,rndIdx2(trainIdx2)), 2 );
											S1 = X1(:,rndIdx1(trainIdx1)) - mu1;	S1 = S1 * S1';
											S2 = X2(:,rndIdx2(trainIdx2)) - mu2;	S2 = S2 * S2';
											W = inv(S1+S2)' * (mu1-mu2);	% Fisher Discriminant vector
											for( m = 2:-1:1 )
												F = W' * [X1(:,rndIdx1(testIdx1)), X2(:,rndIdx2(testIdx2(m,:)))];
												predicts = F < W' * (mu1+mu2)/2;		% false for class 1, true for class 2
												AUC(m,iBoot) = sum( predicts == [T1(rndIdx1(testIdx1)), T2(rndIdx2(testIdx2(m,:)))] ) / size(predicts,2);	% prediction correct rate
												% [tpr, fpr] = roc( [T1(rndIdx1(testIdx1)), T2(rndIdx2(testIdx2(m,:)))], F );	% True Positive Rate, False Positive Rate
												% AUC(m,iBoot) = AreaUnderROC( unique( [tpr;fpr]', 'rows' ) );
											end
										end
										auc(ii,1,:) = mean(AUC,2);
										auc(ii,2,:) = std(AUC,[],2);
									end
									accAUCm(iL,:,iEcc,iTick) = auc(1,1,:);
									accAUCsd(iL,:,iEcc,iTick) = auc(1,2,:);
									segAUCm(iL,:,iEcc,iTick) = auc(2,1,:);
									segAUCsd(iL,:,iEcc,iTick) = auc(2,2,:);
								end

							case 'flda_time'
								%% ROC with Fisher's Linear Discriminant Analysis
								% each variable/dimenstion corresponds to one cell's firing rate at a specific time point, each trial corresponds to one sample
								if( k == 3 )
									data = { frAcc{iEcc,iL}, cat(4, frAcc{iEcc+(1:2)*nEccs,iL}); frSeg{iEcc,iL}, cat(4, frSeg{iEcc+(1:2)*nEccs,iL}) };
									for( ii = 2:-1:1 )
										if( isempty(data{ii,1}) || isempty(data{ii,2}) || ii == 2 )
											auc(ii,1:2,1:2) = 0;
											continue;
										end
										X1 = reshape( data{ii,1}, size(data{ii,1},1), [] );
										X2 = reshape( data{ii,2}, size(data{ii,2},1), [] );
										[coef, scores, latent] = pca( [X1,X2]' );
										iLat = 0; s = 0; upperLat = sum(latent) * 0.80;	% 99% of variance
										while( s < upperLat )
											iLat = iLat+1;
											s = s + latent(iLat);
										end
                                        X1 = reshape( scores(:,1:iLat)', iLat, size(data{ii,1},2), [] );
 										X2 = X1(:,:,nTrials+1:end);
                                        X1 = X1(:,:,1:nTrials);
                                        
										T1 = false(1,nTrials);
										T2 = true(1,nTrials*2);
										N = round(nTrials/3*2);
										trainIdx1 = 1:N;
										trainIdx2 = [1:N, nTrials+(1:N)];
										testIdx1 = N+1:nTrials;
										testIdx2 = (0:1)'*nTrials+(N+1:nTrials);
										nBoots = 3;
										
										if( iTick == size(tTicks,2) && ii == 1 )
											Weights{iEcc,iL} = zeros(size(X1,1),1);
										end
										
										for(iBoot = nBoots:-1:1)
											rndIdx1 = randperm(nTrials);
											rndIdx2 = [randperm(nTrials), randperm(nTrials)+nTrials];
											xx1 = reshape( X1(:,:,rndIdx1(trainIdx1)), size(X1,1), [] );
											xx2 = reshape( X2(:,:,rndIdx2(trainIdx2)), size(X2,1), [] );
											mu1 = mean( xx1, 2 );
											mu2 = mean( xx2, 2 );
											S1 = xx1 - mu1;	S1 = S1 * S1';
											S2 = xx2 - mu2;	S2 = S2 * S2';
											W = inv(S1+S2)' * (mu1-mu2);	% Fisher Discriminant vector
											if( iTick == size(tTicks,2) && ii == 1 )
												Weights{iEcc,iL} = Weights{iEcc,iL} + W;
											end
											for( m = 2:-1:1 )
												F = W' * [reshape( X1(:,:,rndIdx1(testIdx1)), size(X1,1), [] ), reshape( X2(:,:,rndIdx2(testIdx2(m,:))), size(X2,1), [] )];
												predicts = F < W' * (mu1+mu2)/2;		% false for class 1, true for class 2
												predicts = mean( reshape(predicts, size(data{ii,1},2), []), 1 ) > 0.5;
												AUC(m,iBoot) = sum( predicts == [T1(rndIdx1(testIdx1)), T2(rndIdx2(testIdx2(m,:)))] ) / size(predicts,2);	% prediction correct rate
												% [tpr, fpr] = roc( [T1(rndIdx1(testIdx1)), T2(rndIdx2(testIdx2(m,:)))], F );	% True Positive Rate, False Positive Rate
												% AUC(m,iBoot) = AreaUnderROC( unique( [tpr;fpr]', 'rows' ) );
											end
										end
										auc(ii,1,:) = mean(AUC,2);
										auc(ii,2,:) = std(AUC,[],2);

										if( iTick == size(tTicks,2) && ii == 1 )
											Weights{iEcc,iL} = coef(:,1:size(X1,1)) * (Weights{iEcc,iL}/nBoots);
										end
									end
									accAUCm(iL,:,iEcc,iTick) = auc(1,1,:);
									accAUCsd(iL,:,iEcc,iTick) = auc(1,2,:);
									segAUCm(iL,:,iEcc,iTick) = auc(2,1,:);
									segAUCsd(iL,:,iEcc,iTick) = auc(2,2,:);
								end

							case 'flda_single'
								%% PCA + Fisher's Linear Discriminant Analysis in single neuron level; using time points as variables
								sigm = @(x) 1 ./ (1+exp(-x));	% sigmoid function
								if( k == 3 )
									data = { frAcc{iEcc,iL}, cat(4, frAcc{iEcc+(1:2)*nEccs,iL}); frSeg{iEcc,iL}, cat(4, frSeg{iEcc+(1:2)*nEccs,iL}) };
									for( ii = 2:-1:1 )
										N = round(nTrials/3*2);
										T1 = false(1,nTrials);
										T2 = true(1,nTrials*2);
										trainIdx1 = 1:N;
										trainIdx2 = [1:N, nTrials+(1:N)];
										testIdx1 = N+1:nTrials;
										testIdx2 = (0:1)'*nTrials+(N+1:nTrials);

										nBoots = 10;
										for(iBoot = nBoots:-1:1)
											rndIdx1 = randperm(nTrials);
											rndIdx2 = [randperm(nTrials), randperm(nTrials)+nTrials];

											p_present = zeros(2, (nTrials-N)*2, size(data{ii,1},1));
											for(iCell = 1 : size(data{ii,1},1))
												X1 = shiftdim(data{ii,1}(iCell,:,:), 1);
												X2 = shiftdim(data{ii,2}(iCell,:,:), 1);
												[coef, scores, latent] = pca( [X1,X2]' );
												iLat = min(size(scores,2),4);
												X1 = scores(1:nTrials, 1:iLat)';
												X2 = scores(nTrials+1:end, 1:iLat)';

												mu1 = mean( X1(:,rndIdx1(trainIdx1)), 2 );
												mu2 = mean( X2(:,rndIdx2(trainIdx2)), 2 );
												S1 = X1(:,rndIdx1(trainIdx1)) - mu1;	S1 = S1 * S1';
												S2 = X2(:,rndIdx2(trainIdx2)) - mu2;	S2 = S2 * S2';
												W = inv(S1+S2)' * (mu1-mu2);	% Fisher Discriminant vector
												for( m = 2:-1:1 )
													F = W' * [X1(:,rndIdx1(testIdx1)), X2(:,rndIdx2(testIdx2(m,:)))];
													predicts = F < W' * (mu1+mu2)/2;		% false for class 1, true for class 2
													p_present(m,:,iCell) = predicts;
												end
											end
											p_present = mean(p_present,3);
											AUC(:,iBoot) = sum( (p_present > 0.5) == [T1(rndIdx1(testIdx1)), T2(rndIdx2(testIdx2(1,:))); T1(rndIdx1(testIdx1)), T2(rndIdx2(testIdx2(2,:)))], 2 ) / ((nTrials-N)*2);
										end
										auc(ii,1,:) = mean(AUC,2);
										auc(ii,2,:) = std(AUC,[],2);
									end
									accAUCm(iL,:,iEcc,iTick) = auc(1,1,:);
									accAUCsd(iL,:,iEcc,iTick) = auc(1,2,:);
									segAUCm(iL,:,iEcc,iTick) = auc(2,1,:);
									segAUCsd(iL,:,iEcc,iTick) = auc(2,2,:);
								end

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
				hFig = figure( 'NumberTitle', 'off', 'name', sprintf('Example Cells : Trial-by-Trial Prediction : %s+[%d,%d] : %s', alignEvent, LBOffset, UBOffset, layerName), 'color', 'w' );
				pause(0.1);
				jf = get(handle(gcf),'javaframe');
				jf.setMaximized(1);
				pause(1);

				figure(hFig);
				for(iCond = 1 : size(conditions,2))
					iEcc = find( conditions(iCond).eccentricity == Eccs );
					iSF = find( conditions(iCond).sf == [2 10] );
					iPre = find( conditions(iCond).present == [0 1] );
					k = iSF + iPre - 1;

					%% Prabability density of population mean (accumulated) firing rate across trials
					subplot( nRows, nCols, (iEcc-1)*nCols + 1 ); hold on;
					points = squeeze( mean( nansum(frAcc{iCond,iL}(:,:,:), 2), 1 ) )';
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
					drawnow;

					%% ROC by thresholding
					if( k==2 || k==3 )
						subplot( nRows, nCols, (iEcc-1)*nCols + 2 ); hold on;
						t = [zeros(1,nTrials), ones(1,nTrials)];
						data = cat( 4, frAcc{iEcc+([1 k]-1)*nEccs, iL} );
						y = reshape( mean( nansum(data,2), 1 ), 1, [] );
						y = (y-min(y)) / (max(y)-min(y));
						[tpr, fpr, th] = roc( t, y );
						plot( fpr, tpr, 'lineWidth', 2, 'color', colors{k} );
						text( 1, (4-k)*0.2-0.1, sprintf( 'AUC: %.4f', AreaUnderROC( unique( [tpr;fpr]', 'rows' ) ) ), 'fontsize', 18, 'color', colors{k}, 'horizontalAlignment', 'right', 'verticalAlignment', 'bottom' );
						% text( 1, (4-k)*0.2-0.1, sprintf( 'AUC: %.4f', accAUCm(iL,k-1,iEcc,end) ), 'fontsize', 18, 'color', colors{k}, 'horizontalAlignment', 'right', 'verticalAlignment', 'bottom' );

						if( k == 3 )
							set( gca, 'fontsize', 18, 'lineWidth', 2 );
							if( iEcc == 1 )
								title('Thresholding Classifier');
							end
							if( iEcc == 2 )
								ylabel( sprintf('Correct detection rate') );
							end
							if( iEcc == nRows )
								xlabel('False alarm rate');
							end
						end
						drawnow;
					end

					% %% AUC as a function of time computed with accumulated firing rate in the segment of [LBOffset, tTicks(k)]
					% if( k > 1 )
					% 	subplot( 2, 2+2*useGLM, 2 ); hold on;
					% 	plot( tTicks, squeeze(accAUC(iL,k-1,iEcc,:)), 'color', colors{k}*(nEccs-iEcc+1)/nEccs, 'lineWidth', 2, 'displayName', sprintf('SF=%d, Ecc=%d', conditions(iCond).sf, conditions(iCond).eccentricity) );
					% 	if( k == 3 )
					% 		set( gca, 'fontsize', 18, 'lineWidth', 2 );
					% 		% set( legend, 'location', 'northwest' );
					% 		title('AUC Evolving over Time');
					% 		ylabel('AUC');
					% 	end
					% 	drawnow;
					% end

					% %% AUC as a function of time computed with firing rate in the segment [-tWin/2,tWin/2]+tTicks(k)
					% if( k > 1 )
					% 	subplot( 2, 2+2*useGLM, 4 ); hold on;
					% 	plot( tTicks, squeeze(segAUC(iL,k-1,iEcc,:)), 'color', colors{k}*(nEccs-iEcc+1)/nEccs, 'lineWidth', 2, 'displayName', sprintf('SF=%d, Ecc=%d', conditions(iCond).sf, conditions(iCond).eccentricity) );
					% 	if( k == 3 )
					% 		set( gca, 'fontsize', 18, 'lineWidth', 2 );
					% 		set( legend, 'location', 'northwest', 'position', [0.8998, 0.3585, 0.0979, 0.1732] );
					% 		title(sprintf('AUC in Sliding Window of %d ms',tWin));
					% 		xlabel(['Time aligned to ', alignEvent]);
					% 		ylabel('AUC');
					% 	end
					% 	drawnow;
					% end

					%% Classification performance
					if( k == 3 )
						for( m = 1 : 2 )
							subplot( 2, 2, 2 ); hold on;
							M = squeeze(accAUCm(iL,m,iEcc,:))';
							SD = squeeze(accAUCsd(iL,m,iEcc,:))';
							plot( tTicks, M, 'color', colors{m+1}*(nEccs-iEcc+1)/nEccs, 'lineWidth', 2 );
							fill( [tTicks, tTicks(end:-1:1)], [M+SD, fliplr(M)-fliplr(SD)], colors{m}*(nEccs-iEcc+1)/nEccs, 'LineStyle', 'none', 'FaceAlpha', 0.2 );
							set( gca, 'fontsize', 18, 'lineWidth', 2 );
							% title('AUC Evolving over Time');
							% ylabel('AUC');
							title('Prediction Evolving over Time');
							ylabel('Corerct rate');

							subplot( 2, 2, 4 ); hold on;
							M = squeeze(segAUCm(iL,m,iEcc,:))';
							SD = squeeze(segAUCsd(iL,m,iEcc,:))';
							plot( tTicks, M, 'color', colors{m+1}*(nEccs-iEcc+1)/nEccs, 'lineWidth', 2 );
							fill( [tTicks, tTicks(end:-1:1)], [M+SD, fliplr(M)-fliplr(SD)], colors{m}*(nEccs-iEcc+1)/nEccs, 'LineStyle', 'none', 'FaceAlpha', 0.2 );
							set( gca, 'fontsize', 18, 'lineWidth', 2 );
							% title(sprintf('AUC in Sliding Window of %d ms',tWin));
							title(sprintf('Prediction in Sliding Window of %d ms',tWin));
							xlabel(['Time aligned to ', alignEvent]);
							% ylabel('AUC');
							ylabel('Corerct rate');
						end
					end
				end

				if( ~isempty(saveFolder) )
					if( ~exist(saveFolder,'dir') )
						mkdir(saveFolder);
					end
					saveas( gcf, fullfile( saveFolder, sprintf( 'Time Course of Trial-by-Trial Prediction with %s, %s+[%d,%d], %s.fig', classifier, alignEvent, LBOffset, UBOffset, layerName ) ) );
					saveas( gcf, fullfile( saveFolder, sprintf( 'Time Course of Trial-by-Trial Prediction with %s, %s+[%d,%d], %s.png', classifier, alignEvent, LBOffset, UBOffset, layerName ) ) );
					% save(tTicks, accAUC, segAUC, glmAccAUCm, glmAccAUCsd, glmSegAUCm, glmSegAUCsd, frAcc, frSeg);
				end
			end

			% [coef, scores, latent] = pca( fr );
			% scores = [xVel;yVel]' * coef;
			% plot( [0 coef(1,1)*sqrt(latent(1))], [0 coef(2,1)*sqrt(latent(1))], '-', 'color', colors{iCond}, 'LineWidth', 2 );
			% plot( [0 coef(1,2)*sqrt(latent(2))], [0 coef(2,2)*sqrt(latent(2))], '--', 'color', colors{iCond}, 'LineWidth', 2 );
		end


		function [tTicks, accTPRm, accTPRsd, Thresholds, ThresholdsSTD, segTPRm, segTPRsd, frAcc, frSeg] = ContrastDetection(obj, classifier, targetRatio, chance, dataFolder, saveFolder)
			if( nargin() < 2 )
				classifier = 'thresholding-uni';
			end
			if( nargin() < 3 )
				targetRatio = 0.625;
			end
			if( nargin() < 4 )
				change = 0.1;
			end
			if( nargin() < 5 )
				dataFolder = '../../data/figures';
			end
			
			contrasts = [0.01, 0.03, 0.05, 0.07, 0.09, 0.10, 0.20, 0.30];%, 0.40, 0.50];
            layerNames = {'POn', 'POff', 'MOn', 'MOff', 'LayerAverage'};
			if( exist( fullfile(saveFolder, 'PerformanceData.mat') ) )
				load( fullfile(saveFolder, 'PerformanceData.mat') );
			else
				fprintf( 'Loading %s ...\n', fullfile(dataFolder, 'LFR.mat') );
                load( fullfile(dataFolder, 'LFR.mat') );
				% folders = cellfun( @(c) {sprintf('Real ratio cells, avContrast=0.5, contrast=%.2f, same 30 trials; Hs=1', c)}, num2cell(contrasts) );
				% for( iFolder = 1 : size(folders,2) )
				for( iContrast = 1 : size(contrasts,2) )
					% fprintf( 'Processing %s ...\n', fullfile(dataFolder, folders{iFolder}, [folders{iFolder}, '.mat']) );
					fprintf( 'Processing contrast=%.2f ...\n', contrasts(iContrast) );
					% load( fullfile(dataFolder, folders{iFolder}, [folders{iFolder}, '.mat']), 'LFR', 'time', 'conditions', 'trials', 'trialsIdx', "cellIdx" );
					condIdx = [1:3, (4:9)+(iContrast-1)*6];
					tmpConditions = conditions(condIdx);
					[tmpConditions(1:3).sf] = deal(2);
					[tTicks, accTPRm(:,:,:,:,iContrast), accTPRsd(:,:,:,:,iContrast), segTPRm(:,:,:,:,iContrast), segTPRsd(:,:,:,:,iContrast), frAcc(:,:,iContrast), frSeg(:,:,iContrast)] = obj.TrialPredictWithAUC(classifier, LFR(condIdx,:), time(1:9), tmpConditions, trials, trialsIdx, nCells, cellIdx, cell2AnalyzeIdx, true, 'saccadeOff', 0, 500, [], [], fullfile(saveFolder, sprintf('contrast=%.2f',contrasts(iContrast))));
					drawnow;
					pause(3);
					for(iL = 5:-1:1)
% 						saveas( gcf, fullfile( saveFolder, sprintf( 'Time Course of Trial-by-Trial Prediction with %s, %s+[%d,%d], %s.fig', classifier, 'saccadeOff', 0, 500, layerNames{iL} ) ) );
% 						saveas( gcf, fullfile( saveFolder, sprintf( 'Time Course of Trial-by-Trial Prediction with %s, %s+[%d,%d], %s.png', classifier, 'saccadeOff', 0, 500, layerNames{iL} ) ) );
						close(gcf);
					end
					close all;
				end
				if( ~exist( saveFolder ) )
					mkdir(saveFolder);
				end
				save( fullfile(saveFolder, 'PerformanceData.mat'), 'tTicks', 'accTPRm', 'accTPRsd', 'segTPRm', 'segTPRsd', 'frAcc', 'frSeg', 'conditions' );
			end

			Eccs = unique([conditions.eccentricity]);
			nEccs = size(Eccs,2);
			SFs = unique([conditions.sf]);
			SFs(SFs == 0) = [];
			nSFs = size(SFs,2);
			if( ~exist('durs', 'var') )
				durs = [50 150 300 500];
				durs = [1 50 75 150 300 500];
			end
            nDurs = size(durs,2);
			durOffset = 50 + 7;		% resonse delay of 50 ms + online saccade off later by 7 ms 
			layerNames = {'POn', 'POff', 'MOn', 'MOff', 'LayerAverage'};

			%% psychometric curves
			if( ~exist('Thresholds', 'var') )
				Thresholds = zeros(size(layerNames,2),nSFs,nEccs,size(durs,2));
				ThresholdsSTD = Thresholds;
				fitThreshold = true;
			else
				fitThreshold = true;false;
			end
			CHANCE = chance;
			for( iL = 1 : 5 )
				colors = { [0 1 0], [0 0.7 0], [0 0.4 0], [0 0.1 0]; ...	% for 2cpd
						   [1 0 0], [0.7 0 0], [0.4 0 0], [0.1 0 0] };		% for 10cpd
				colors = { [0 1 0], [0 0.775 0], [0 0.55 0], [0 0.325 0], [0 0.1 0]; ...	% for 2cpd
						   [1 0 0], [0.775 0 0], [0.55 0 0], [0.325 0 0], [0 0.1 0] };		% for 10cpd
				figure( 'NumberTitle', 'off', 'name', sprintf('Psychometric Curve - Prediction by AverageThresholding - %s', layerNames{iL}), 'color', 'w' );
				pause(0.1);
				jf = get(handle(gcf),'javaframe');
				jf.setMaximized(1);
				pause(1);
				for(iSF = 1 : nSFs)
					for(iEcc = 1 : nEccs)
						subplot( nSFs, nEccs, (iSF-1)*nEccs + iEcc ); hold on;
						for( iDur = 1 : nDurs )
							chance = CHANCE - max(0, durs(iDur)/500*0.05);
							targetRatio = (1+chance)/2;
							[~,iTick] = min(abs(tTicks - (durs(iDur)+durOffset)));
							h(iDur) = plot( contrasts, squeeze(accTPRm(iL,iSF,iEcc,iTick,:)), 'o', 'color', colors{iSF,1}*(nDurs-iDur+1)/nDurs, 'displayName', sprintf('%d ms', durs(iDur)), 'markersize', 8, 'lineWidth', 3 );
							plot( contrasts, squeeze(accTPRm(iL,iSF,iEcc,iTick,:)), '--', 'color', colors{iSF,1}*(nDurs-iDur+1)/nDurs, 'displayName', sprintf('%d ms', durs(iDur)), 'lineWidth', 1 );

							if(fitThreshold)
								if( ~any(isnan(accTPRm(iL,iSF,iEcc,iTick,:))) )
									nTrials = 100;%size(frSeg{1},3);
									hits = zeros(1, nTrials*size(contrasts,2));
									for( iCont = 1 : size(contrasts,2) )
										hits( (iCont-1)*nTrials + ( 1 : round(nTrials * accTPRm(iL,iSF,iEcc,iTick,iCont)) ) ) = 1;
									end
									[~, ~, nThresh, nPar, g, chisq] = psyfit( reshape(repmat(contrasts,nTrials,1),1,[]), hits, 'Thresh', targetRatio, 'Chance', chance, 'Lapses', 0, 'Log', 'Extra', 'PlotOff', 'Boots', 10, 'disttype', 'normal' );
									par = mean(nPar,2);
									stdPar(2) = std(nPar(2,:));
									stdPar(1) = std(nPar(1,:));
									thresh = mean(nThresh);
									stdThresh = std(nThresh);
									x = linspace( 0, max(contrasts)*1.1, 10000 );
									y = psyfun( x, par(1), par(2), chance, 0, false, true, 'normal' );
									yLow = psyfun( x, par(1) - stdPar(1), par(2) - stdPar(2), 0.1, 0, false, true, 'normal' );
									yUp = psyfun( x, par(1) + stdPar(1), par(2) + stdPar(2), 0.1, 0, false, true, 'normal' );
									plot( x, y, '-', 'LineWidth', 2.5, 'color', colors{iSF,1}*(nDurs-iDur+1)/nDurs, 'DisplayName', sprintf('%d cpd, %d ms', SFs(iSF), durs(iDur)) );
									% set( fill( [x(2:end) x(end:-1:2)], [yLow(2:end) yUp(end:-1:2)], colors{iSF,iDur} ), 'LineStyle', 'none', 'FaceAlpha', 0.5 );
									plot( [1, 1] * thresh, [0, targetRatio], '--', 'color', colors{iSF,1}*(nDurs-iDur+1)/nDurs, 'lineWidth', 2 );
									% set( fill( [-1, 1, 1, -1]*stdThresh + thresh, [0, 0, targetRatio, targetRatio], colors{iSF,iDur} ), 'LineStyle', 'none', 'FaceAlpha', 0.5 );

									Thresholds(iL,iSF,iEcc,iDur) = thresh;
									ThresholdsSTD(iL,iSF,iEcc,iDur) = stdThresh;
								else
									Thresholds(iL,iSF,iEcc,iDur) = nan;
									ThresholdsSTD(iL,iSF,iEcc,iDur) = nan;
								end
							end
						end
						if(iSF == 1 && iEcc == 1)
							legend(h, 'location', 'southeast');
						end
						if(iSF == nSFs)
							xlabel('Contrast');
						end
						if(iEcc == 1)
							ylabel( {sprintf('SF=%d', SFs(iSF)), 'Proportion yes'} );
                        end
                        if(iSF == 1)
                            title(sprintf('Ecc = %d', Eccs(iEcc)));
                        end
						set( gca, 'xlim', [0.01 0.6], 'xscale', 'log', 'ylim', [0 1], 'lineWidth', 2, 'fontsize', 18 );
						drawnow;
					end
				end
				saveas( gcf, fullfile( saveFolder, sprintf('Psychometric Curve - Prediction by AverageThresholding - %s.fig', layerNames{iL}) ) );
				saveas( gcf, fullfile( saveFolder, sprintf('Psychometric Curve - Prediction by AverageThresholding - %s.png', layerNames{iL}) ) );
			end

			if(fitThreshold)
				%save( fullfile(dataFolder, saveFolder, 'PerformanceData.mat'), 'tTicks', 'accTPRm', 'accTPRsd', 'segTPRm', 'segTPRsd', 'frAcc', 'frSeg', 'conditions', 'Thresholds', 'ThresholdsSTD' );
                save( fullfile(saveFolder, 'PerformanceData.mat'), 'tTicks', 'accTPRm', 'accTPRsd', 'segTPRm', 'segTPRsd', 'conditions', 'durs', 'Thresholds', 'ThresholdsSTD' );
			end


			%% contrast threshold as a function of stimulus duration
			figure( 'NumberTitle', 'off', 'name', sprintf('Contrast Threshold - Prediction by AverageThresholding - %s', layerNames{iL}), 'color', 'w' );
			for( iL = 1 : 5 )
				pause(0.1);
				jf = get(handle(gcf),'javaframe');
				jf.setMaximized(1);
				pause(1);
				nRows = nSFs;
				nCols = size(layerNames,2);
				colors = { [0 1 0], [0 0.6 0], [0 0.2 0]; ...	% for 2cpd
						   [1 0 0], [0.6 0 0], [0.2 0 0] };		% for 10cpd
				for(iSF = 1 : nSFs)
					subplot(nRows, nCols, (iSF-1)*nCols + iL); hold on; h = [];
					for(iEcc = 1 : nEccs)
						m = squeeze(Thresholds(iL,iSF,iEcc,:))';
						sd = squeeze(ThresholdsSTD(iL,iSF,iEcc,:))';
						h(iEcc) = plot( durs, m, 's', 'color', colors{iSF,iEcc}, 'markersize', 10, 'lineWidth', 3, 'displayName', sprintf('Ecc=%d', Eccs(iEcc)) );
						plot( durs, m, '--', 'color', colors{iSF,iEcc}, 'lineWidth', 1 );
						plot( reshape([durs; durs; nan(size(durs))], 1, []), reshape([m+sd;m-sd;nan(size(durs))], 1, []), 'color', colors{iSF,iEcc}, 'lineWidth', 2 );
					end
					if(iSF == 1 && iL == 1)
						legend(h, 'location', 'southeast');
					end
					if(iSF == nSFs)
						xlabel('Duration (ms)');
					end
					if(iL == 1)
						ylabel( {sprintf('SF=%d', SFs(iSF)), 'Contrast threshold'} );
                    end
                    if(iSF == 1)
                        title(layerNames{iL});
                    end
					set( gca, 'xlim', [0 600], 'xtick', [1 50 150 500], 'ylim', [0.008 1], 'ytick', [0.01, 0.05, 0.1 0.5], 'xscale', 'log', 'yscale', 'log', 'ydir', 'reverse', 'lineWidth', 2, 'fontsize', 18 );
					if( max(get(gca,'ylim')) > 1 )
						set( gca, 'ylim', [0 1] );
					end
					drawnow;
				end
			end
			saveas( gcf, fullfile( saveFolder, sprintf('Contrast Threshold - Prediction by AverageThresholding - %s.fig', layerNames{iL}) ) );
			saveas( gcf, fullfile( saveFolder, sprintf('Contrast Threshold - Prediction by AverageThresholding - %s.png', layerNames{iL}) ) );
			
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


		%% Video4ExampleCells
		function Video4ExampleCells(arg)
			%%
			Eccs = unique([conditions.eccentricity]);
			nEccs = size(Eccs,2);
			SFs = [0, 2, 10];
			contrasts = [0, 0.5, 0.5];
			colors = { [1 0 0], [0 0 1], [1 0 1], [0 1 1] };

			iTrial = 1;
			fr = reshape( cat( 1, LFR{ [conditions.contrast] == contrasts(1) | [conditions.contrast] == contrasts(2), : } ), 1, [] );
			fr(fr<0) = 0;
			fr(isnan(fr)) = [];
			frMax = std(fr)*5;

			%% get average radii
			for(iL = 4:-1:1)
				for(iEcc = size(Eccs,2) : -1 : 1)
					Radii(iL,iEcc) = mean( [encoder.layers(iL).sRFParams(cellIdx{iL,iEcc}).centerRadii] ) / 2;
				end
			end

			%%
			for( iCont = 1 : size(contrasts,2) )

				filename = sprintf('../../Data/videos/contrast=%.2f, SF=%d', contrasts(iCont), SFs(iCont));
				writerObj = VideoWriter(filename, 'MPEG-4');
				open(writerObj);

				figure( 'NumberTitle', 'off', 'color', 'k', 'name', sprintf('Movie: contrast=%0.2f, SF=%d', contrasts(iCont), SFs(iCont)) );
				pause(0.1);
				jf = get(handle(gcf),'javaframe');
				jf.setMaximized(1);
				pause(0.5);
				for( iEcc = 1 : size(Eccs,2) )
					for( iL = 1 : 4 )
						subplot(nEccs, 4, (iEcc-1)*4 + iL); hold on;
						set( gca, 'fontsize', 20, 'linewidth', 2, 'color', 'k', 'XColor', 'w', 'YColor', 'w' );
						axis equal;
					end
				end

				tStart = -50;
				for( iTick = find(time{1} == tStart) : size(time{1},2) )

					% create graphic objects
					if( iTick == find(time{1} == tStart) )
						for( iEcc = 1 : size(Eccs,2) )
							for( iL = 1 : 4 )
								subplot(nEccs, 4, (iEcc-1)*4 + iL);
								title(sprintf('%s | Ecc=%d | t=%dms', encoder.layers(iL).name, Eccs(iEcc), time{1}(iTick)), 'color', 'w');
								for( iCell = 1 : nCellsUsed{iL,iEcc} )
									r = Radii(iL,iEcc);
									iCond = [conditions.contrast] == contrasts(iCont) & [conditions.eccentricity] == Eccs(iEcc) & [conditions.sf] == SFs(iCont);
									% h{iTick,iL,iEcc}(iCell) = rectangle( 'position', [ encoder.layers(iL).locations(cellIdx{iL,iEcc}(iCell), 1), encoder.layers(iL).locations(cellIdx{iL,iEcc}(iCell), 2), 2*r, 2*r ],...
									% 									 'FaceColor', 'k', 'LineStyle', 'none', 'Curvature', [1 1] );
									h{iL,iEcc}(iCell) = fill( r*cosd(0:60:300)+encoder.layers(iL).locations(cellIdx{iL,iEcc}(iCell), 1), r*sind(0:60:300)+encoder.layers(iL).locations(cellIdx{iL,iEcc}(iCell), 2),...
																	'k', 'FaceColor', 'k', 'FaceAlpha', 0.5, 'LineStyle', 'none' );
								end
								eyeIdx = time{1}(iTick) + trials(trialsIdx{1}(iTrial)).(conditions(iCond).alignEvent) + (-10:0);
								if(iEcc == 1)
									ecc = 0.75;
								else
									ecc = Eccs(iEcc);
								end
								hEye{iL,iEcc} = plot( -trials(trialsIdx{1}(iTrial)).x.position(eyeIdx)/60 + ecc, -trials(trialsIdx{1}(iTrial)).y.position(eyeIdx)/60, 'w-', 'linewidth', 1 );
							end
						end
						drawnow;
					end

					% update graphic objects for each time point
					for( iEcc = 1 : size(Eccs,2) )
			            iCond = [conditions.contrast] == contrasts(iCont) & [conditions.eccentricity] == Eccs(iEcc) & [conditions.sf] == SFs(iCont);
						
						for( iL = 1 : 4 )
							subplot(nEccs, 4, (iEcc-1)*4 + iL);
							title(sprintf('%s | Ecc=%d | t=%dms', encoder.layers(iL).name, Eccs(iEcc), time{1}(iTick)), 'color', 'w');
								
							for( iCell = 1 : nCellsUsed{iL,iEcc} )
								h{iL,iEcc}(iCell).FaceColor = colors{iL} * min( 1, max(0, LFR{iCond,iL}(iCell,iTick,iTrial)) / frMax );
							end
							eyeIdx = time{1}(iTick) + trials(trialsIdx{1}(iTrial)).(conditions(iCond).alignEvent) + (-10:0);
							if(iEcc == 1)
								ecc = 0.75;
							else
								ecc = Eccs(iEcc);
							end
							hEye{iL,iEcc}.XData = -trials(trialsIdx{1}(iTrial)).x.position(eyeIdx)/60 + ecc;
							hEye{iL,iEcc}.YData = -trials(trialsIdx{1}(iTrial)).y.position(eyeIdx)/60;
						end
					end
					drawnow;
					writeVideo(writerObj, getframe(gcf));
					% pause;
				end

				close(writerObj);
			end
		end
		
	end
end