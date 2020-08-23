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


		function [LFR, time, conditions, sFR, sFR_c, sFR_s, trials, trialsIdx, cellIdx] = ExampleCells(obj, dataFolder, alignEvent, contrast)
			%% Show example modeled cells at locations (0,0), (0,4), (0,8), (0,12)
			%   dataFolder:			folder containing experiment data and noise background
			%	alignEvent:			'saccadeOn', 'saccadeOff', 'flashOn'

			if( nargin() < 1 || isempty(dataFolder) )
				dataFolder = '../../data/';
			end
			if( nargin() < 2 || isempty(alignEvent) )
				alignEvent = 'saccadeOff';
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
					else
						w = trials(1).gratingWidth/2;
						cellIdx{iL,iEcc} = find( (eccs(iEcc)-w)^2 <= d2 & d2 <= (eccs(iEcc)+w)^2 & index );
					end
					cellIdx{iL,iEcc} = cellIdx{iL,iEcc}( randperm( size(cellIdx{iL,iEcc}, 1) ) );
					% [~, cellIdx(iL,iEcc)] = min( sum( (obj.layers(iL).locations - cLocs(iEcc,:)).^2, 2 ) );
				end
			end
			nCells = 30;
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
							'sf',			{  0,   0,   0,   0,   2,   2,   2,   2,  10,  10,  10,  10}, ...
							'duration',		{500}, ...
							'present',		{  0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1}, ...
							'alignEvent',	{alignEvent} );

			for( iCond = size(conditions,2) : -1 : 1 )
				iEcc = find( conditions(iCond).eccentricity == [0, 4, 8, 12] );

				trialsIdx{iCond} = find( ... 
										 ...%[trials.eccentricity] == eccs(1,iEcc) & ...
										 ... %[trials.spatialFreq] == SFs(iSF) & ...
										 abs( [trials.stimOff] - [trials.saccadeOff] - conditions(iCond).duration ) < 50 & ...
										 true );%[trials.present] == presents(iPre) );

				idx = trialsIdx{iCond};    %fprintf('nTrials = %d\n', size(idx,2)); continue;
                idx = idx( 1 : min(nTrials,end) );
				tMax = round(max( [trials(idx).saccadeOff] + round(conditions(iCond).duration/1000*[trials(idx).sRate]) - [trials(idx).(alignEvent)] ));	% in samples
				tMin = round(min( [trials(idx).saccadeOn] - round(0.300*[trials(idx).sRate]) - [trials(idx).(alignEvent)] ));				% in samples
				sfr = zeros( nCells, tMax-tMin+1, size(idx,2), 4 );		% 1st dim: neurons;		2nd dim: time;		3rd dim: trials;	4th dim: layers
				sfr_c = sfr;
				sfr_s = sfr;
				tfr = sfr;
				for( k = size(idx,2) : -1 : 1 )
					tic;
					x  = trials(idx(k)).x.position / 60;
					y  = trials(idx(k)).y.position / 60;
                    
					% before grating
					eyeIdx1 = tMin + trials(idx(k)).(alignEvent) : trials(idx(k)).flashOn-1;
					s1 = size(eyeIdx1,2);

					% during grating
					eyeIdx2 = trials(idx(k)).flashOn : min( size(x,2), tMax + trials(idx(k)).(alignEvent) );
					s2 = size(eyeIdx2,2);

                    [noise, inputX, inputY] = obj.LoadNoise( fullfile(dataFolder, trials(idx(k)).backgroundImage), trials(idx(k)).pixelAngle/60 );
					grating = obj.GenerateGrating( conditions(iCond).eccentricity, trials(idx(k)).gratingWidth, conditions(iCond).sf, trials(idx(k)).phase, trials(idx(k)).pixelAngle/60 );
					noise = noise * trials(idx(k)).backgroundContrast;
					grating = grating * contrast * conditions(iCond).present;
                    bg = 1;
                    
					% avContrast = ((s1+s2)*trials(idx(k)).backgroundContrast + s2*contrast*conditions(iCond).present) / (s1+s2);
					avContrast = trials(idx(k)).backgroundContrast;
					% contrastF = [ zeros(1,s1), ones(1,s2) * contrast ] + trials(idx(k)).backgroundContrast;
					% contrastF = obj.ComputeContrasts( noise+grating+bg,  inputX, inputY, obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),1), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),2), [obj.layers(iL).sRFParams(cellIdx{iL,iEcc}(1:nCells)).surroundRadii], x([eyeIdx1,eyeIdx2]), y([eyeIdx1,eyeIdx2]) );
					contrastF = contrast;

					for( iL = 1:4 )
						xIdx = min(x([eyeIdx1,eyeIdx2])) + min(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),1)) - 4 <= inputX & inputX <= max(x([eyeIdx1,eyeIdx2])) + max(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),1)) + 4; %720 : 1460;
	                    yIdx = min(y([eyeIdx1,eyeIdx2])) + min(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),2)) - 4 <= inputY & inputY <= max(y([eyeIdx1,eyeIdx2])) + max(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),2)) + 4; %300 : 780;
	                    
						[sfr(:,s1+(1:s2),k,iL), sfr_c(:,s1+(1:s2),k,iL), sfr_s(:,s1+(1:s2),k,iL)] = obj.SpatialModel.LinearResponse( noise(yIdx,xIdx)+grating(yIdx,xIdx), inputX(xIdx), inputY(yIdx), x(eyeIdx2), y(eyeIdx2), obj.layers(iL).sRFParams(cellIdx{iL,iEcc}(1:nCells)), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),1), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),2));
						[sfr(:,1:s1,k,iL), sfr_c(:,1:s1,k,iL), sfr_s(:,1:s1,k,iL)] = obj.SpatialModel.LinearResponse( noise(yIdx,xIdx), inputX(xIdx), inputY(yIdx), x(eyeIdx1), y(eyeIdx1), obj.layers(iL).sRFParams(cellIdx{iL,iEcc}(1:nCells)), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),1), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),2));

						if( lower(obj.layers(iL).name(1)) == 'm' )
							tfr(:,:,k,iL) = obj.TemporalModel.LinearResponse( obj.layers(iL).name, obj.layers(iL).tRFParams(cellIdx{iL,iEcc}(1:nCells)), avContrast, trials(idx(k)).sRate, contrastF .* sfr(:,:,k,iL) );
						else
							tfr(:,:,k,iL) = obj.TemporalModel.LinearResponse( obj.layers(iL).name, obj.layers(iL).tRFParams(cellIdx{iL,iEcc}(1:nCells)), avContrast, trials(idx(k)).sRate, contrastF .* sfr_c(:,:,k,iL), contrastF .* sfr_s(:,:,k,iL) );
						end
					end
					fprintf('Ecc=%d, SF=%d, Dur=%d, Present=%d, iTrials=%d/%d, t=%f\n', conditions(iCond).eccentricity, conditions(iCond).sf, conditions(iCond).duration, conditions(iCond).present, k, size(idx,2), toc);
				end
				sFR{iCond} = sfr;
				sFR_c{iCond} = sfr_c;
				sFR_s{iCond} = sfr_s;
				LFR{iCond} = tfr;
				if(isempty(idx))
					time{iCond} = [];
				else
					time{iCond} = (tMin:tMax) / trials(idx(1)).sRate * 1000;
				end
			end
		end


		function DisplayExampleCells(obj, FR, time, conditions, alignEvent, trials, trialsIdx, cellIdx)


			conditions = struct( ...
							'eccentricity',	{  0,   4,   8,  12,   0,   4,   8,  12,   0,   4,   8,  12}, ...
							'sf',			{  2,   2,   2,   2,   2,   2,   2,   2,  10,  10,  10,  10}, ...
							'duration',		{500}, ...
							'present',		{  0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1}, ...
							'alignEvent',	{'flashOn'} );
			
			% nTrials = size( FR{1}, 1 );
			% nCells = size( FR{1}, 3 );
			nTrials = size( FR{1}, 3 );
			nCells = size( FR{1}, 1 );
			idx = trialsIdx{1}(1:nTrials);

			nRows = size( unique([conditions.eccentricity]), 2 );
			nCols = 1;

			for( iL = 1 : 4 )
				figure( 'NumberTitle', 'off', 'name', ['Example Cells: ', obj.layers(iL).name], 'color', 'w' );
				pause(0.1);
				jf = get(handle(gcf),'javaframe');
				jf.setMaximized(1);
				pause(1);

				colors = {'b', 'g', 'r'};
				names = {'Absent', '2 cpd', '10 cpd'};
				for(iCond = 1 : size(conditions,2))
					iEcc = find( conditions(iCond).eccentricity == [0 4 8 12] );
					iSF = find( conditions(iCond).sf == [2 10] );
					iDur = find( conditions(iCond).duration == [500] );
					iPre = find( conditions(iCond).present == [0 1] );

					subplot( nRows, nCols, (iEcc-1)*nCols + 1 ); hold on;
					k = iSF + iPre - 1;
					% fr = FR{iEcc,iSF,iDur,iPre}(:,:,iL);
					fr = FR{iCond}(:,:,:,iL);
					fr(fr<0) = 0;
                    fr = mean( fr, 3 );
					m = mean( fr, 1 );
					sem = std( fr, [], 1 ) / sqrt(size(fr,1));

					% t = time{iEcc,iSF,iDur,iPre};
					t = time{iCond};

					h(k) = plot( t, m, 'color', colors{k}, 'lineWidth', 2, 'displayName', names{k} );
					fill( [t, t(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], colors{k}, 'LineStyle', 'none', 'FaceAlpha', 0.2 );
					if( k == 3 )
						plot( [0 0], get(gca, 'ylim'), 'k--', 'lineWidth', 2 );
						plot( [1 1] * mean([trials(idx).saccadeOn] - [trials(idx).flashOn]), get(gca, 'ylim'), 'b--', 'lineWidth', 2 );
						plot( [1 1] * mean([trials(idx).saccadeOff] - [trials(idx).flashOn]), get(gca, 'ylim'), 'r--', 'lineWidth', 2 );
						plot( [1 1] * mean([trials(idx).saccadeLand] - [trials(idx).flashOn]), get(gca, 'ylim'), 'm--', 'lineWidth', 2 );
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
			subplot(3, 4, [1 2]); hold on;
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
			subplot(3, 4, [3 4]); hold on;
			for( k = size(idx,2) : -1 : 1 )
				x = trials(idx(k)).x.position( -100+trials(idx(k)).flashOn : min(end, 150+trials(idx(k)).flashOn) );
				eyeX( k, 1:size(x,2) ) = x;
			end
			m = mean( eyeX, 1 );
			sd = std( eyeX, [], 1 );
			plot( -100:150, m, 'b', 'lineWidth', 2 );
			fill( [-100:150, 150:-1:-100], [m-sd, m(end:-1:1)+sd(end:-1:1)], 'b', 'LineStyle', 'none', 'FaceAlpha', 0.5 );
			plot( [0 0], get(gca,'ylim'), 'k--', 'lineWidth', 2 );
			xlabel('Time aligned to flashOn');
			ylabel('Horizontal eye (mean+-std) (arcmin)');
			set( gca, 'lineWidth', 2, 'fontsize', 16 );

			% cell spatial sensitivity function
			colors = {[1 0 0], [0 0 1], 'm', 'c'};
			sf = 0.1:0.1:20;
			for(iEcc = 1:4)
				subplot(3, 4, 4+iEcc); hold on; h = [];
				for(iL = 1 : 4)
					% idx = cellIdx(iL,iEcc);
					idx = cellIdx{iL,iEcc};
					m = mean( abs( obj.SpatialModel.SpatialSensitivity( obj.layers(iL).sRFParams(idx), sf ) ), 1 );
					sem = std( abs( obj.SpatialModel.SpatialSensitivity( obj.layers(iL).sRFParams(idx), sf ) ), [], 1 ) / sqrt(nCells);
					h(iL) = plot( sf, m, 'color', colors{iL}, 'lineWidth', 2, 'displayName', obj.layers(iL).name );
					fill( [sf, sf(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], 'k', 'FaceColor', colors{iL}, 'LineStyle', 'none', 'FaceAlpha', 0.5 );
				end
				plot( [2 2], get(gca,'ylim'), 'k--', 'lineWidth', 2 );
				plot( [10 10], get(gca,'ylim'), 'k--', 'lineWidth', 2 );
				if(iEcc==4), legend(h); end
				xlabel('Spatial frequency (CPD)');
				if(iEcc==1), ylabel('Sensitivity'); end
				title( sprintf('Ecc = %d', conditions(iEcc).eccentricity) );
				set( gca, 'lineWidth', 2, 'fontsize', 16 );
			end

			% cell temporal sensitivity function
			tf = 0.1:0.1:60;
			for(iEcc = 1:4)
				subplot(3, 4, 8+iEcc); hold on; h = [];
				for(iL = 1 : 4)
					% idx = cellIdx(iL,iEcc);
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
				if(iEcc==4), legend(h); end
				xlabel('Temporal frequency (Hz)');
				if(iEcc == 1), ylabel('Sensitivity'); end
				set( gca, 'lineWidth', 2, 'fontsize', 16 );
			end
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