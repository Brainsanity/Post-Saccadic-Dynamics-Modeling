%% display neurons
% figure( 'NumberTitle', 'off', 'color', 'w');
% SynIdx.POn = [cellIdx{1,1}(cell2AnalyzeIdx{1,1}); cellIdx{1,2}(cell2AnalyzeIdx{1,2}); cellIdx{1,3}(cell2AnalyzeIdx{1,3})];
% SynIdx.POff = [cellIdx{2,1}(cell2AnalyzeIdx{2,1}); cellIdx{2,2}(cell2AnalyzeIdx{2,2}); cellIdx{2,3}(cell2AnalyzeIdx{2,3})];
% SynIdx.MOn = [cellIdx{3,1}(cell2AnalyzeIdx{3,1}); cellIdx{3,2}(cell2AnalyzeIdx{3,2}); cellIdx{3,3}(cell2AnalyzeIdx{3,3})];
% SynIdx.MOff = [cellIdx{4,1}(cell2AnalyzeIdx{4,1}); cellIdx{4,2}(cell2AnalyzeIdx{4,2}); cellIdx{4,3}(cell2AnalyzeIdx{4,3})];
SynIdx.POn = [cellIdx{1,1}; cellIdx{1,2}; cellIdx{1,3}];
SynIdx.POff = [cellIdx{2,1}; cellIdx{2,2}; cellIdx{2,3}];
SynIdx.MOn = [cellIdx{3,1}; cellIdx{3,2}; cellIdx{3,3}];
SynIdx.MOff = [cellIdx{4,1}; cellIdx{4,2}; cellIdx{4,3}];
% encoder.SpatialModel.DisplaySynthesizedRFParams(encoder.layers(1).sRFParams(SynIdx.POn), encoder.layers(3).sRFParams(SynIdx.MOn));
% encoder.SpatialModel.DisplaySynthesizedRFParams(encoder.layers(2).sRFParams(SynIdx.POff), encoder.layers(4).sRFParams(SynIdx.MOff));

%%
figure('color', 'w');
					iL = 1;
					t = encoder.activityParams.timeline;
					for(iCond = 2:-1:1)
						fr{iCond} = encoder.ExampleCellsActivitiesOnCondition(encoder.layers(iL).name, 0, (iCond-1)*2, 0.5, 'saccadeOff', t([1 end]));
						for(iTrial = 1 : 4)
							subplot(2,4,iTrial+(iCond-1)*4); hold on;
							for(iCell = 1:50)
								plot3( t, iCell*ones(size(encoder.activityParams.timeline)),  fr{iCond}(iCell, :, iTrial) );
							end
							set(gca, 'xlim', [-150 550], 'view', [14.7560 19.3337], 'fontsize', 20, 'linewidth', 2);
							if(iTrial==4) ylabel('Neuron'); end
							if(iTrial==1) xlabel('Time (ms)'); end
							zlabel('Firing rate');
						end
					end

					figure('color', 'w'); hold on;
					[data, r] = hist( squeeze(sum(mean(fr{1}(:, 0<=t & t<=500, :), 1), 2)), 50 );
					thresh = prctile( squeeze(sum(mean(fr{1}(:, 0<=t & t<=500, :), 1), 2)), 90 );
					plot(r, data / sum(data) / (r(2)-r(1)), 'b-', 'lineWidth', 2);
					
					[data, r] = hist( squeeze(sum(mean(fr{2}(:, 0<=t & t<=500, :), 1), 2)), 50 );
					plot(r, data / sum(data) / (r(2)-r(1)), 'r-', 'lineWidth', 2);
					plot([1 1]*thresh, get(gca, 'ylim'), 'k-', 'lineWidth', 2);
					xlabel('Accumulate population mean FR (spikes/s)');ylabel('Probability density');
					set(gca, 'fontsize', 20, 'lineWidth', 2);


%%
obj = encoder;
% encoder2 = Encoder([], true);
% encoder2.DisplayExampleCells(LFR, time, conditions, 'saccadeOff', trials, trialsIdx, cellIdx);
iL = 1; plot( [obj.layers(iL).sRFParams.temporalEccDegs], ([obj.layers(iL).sRFParams.surroundRadii] ./ [obj.layers(iL).sRFParams.centerRadii]).^2 .* ([obj.layers(iL).sRFParams.surroundPeakSensitivities] ./ [obj.layers(iL).sRFParams.centerPeakSensitivities]), '.' )
%%
subplot(2,2,1); hold on; cla;
for(iL = [1 3])
	plot( [encoder2.layers(iL).sRFParams.centerRadii].^2 .* [encoder2.layers(iL).sRFParams.centerPeakSensitivities], [encoder2.layers(iL).sRFParams.surroundRadii].^2 .* [encoder2.layers(iL).sRFParams.surroundPeakSensitivities], '.' )
end
plot( [min(get(gca,'xlim'), get(gca,'ylim')), max(get(gca,'xlim'), get(gca,'ylim'))], [min(get(gca,'xlim'), get(gca,'ylim')), max(get(gca,'xlim'), get(gca,'ylim'))], 'k-' );

subplot(2,2,2); hold on; cla;
for(iL = [1 3])
	plot( [encoder2.layers(iL).sRFParams.eccDegs], [encoder2.layers(iL).sRFParams.surroundRadii] ./ [encoder2.layers(iL).sRFParams.centerRadii], '.' );
end
plot( [min(get(gca,'xlim'), get(gca,'ylim')), max(get(gca,'xlim'), get(gca,'ylim'))], [min(get(gca,'xlim'), get(gca,'ylim')), max(get(gca,'xlim'), get(gca,'ylim'))], 'k-' );

subplot(2,2,3); hold on; cla;
for(iL = [1 3])
	plot( [encoder2.layers(iL).sRFParams.centerPeakSensitivities], [encoder2.layers(iL).sRFParams.surroundPeakSensitivities], '.' );
end
plot( [min(get(gca,'xlim'), get(gca,'ylim')), max(get(gca,'xlim'), get(gca,'ylim'))], [min(get(gca,'xlim'), get(gca,'ylim')), max(get(gca,'xlim'), get(gca,'ylim'))], 'k-' );

%%
figure
obj = encoder;
for(iPM = 1 : 2)
	subplot(2,2,iPM); hold on;
	% rCenter = obj.SpatialModel.([PM{iPM} 'CenterRadiusFunction'])( obj.SpatialModel.([PM{iPM} 'CenterRadiusParams']), ecc);
	% rSurround = obj.SpatialModel.([PM{iPM} 'SurroundRadiusFunction'])( obj.SpatialModel.([PM{iPM} 'SurroundRadiusParams']), ecc);
	[~, spacing] = WatsonRGCModel.RFSpacingDensityMeridian( ecc, WatsonRGCModel.enumeratedMeridianNames{iMeridian}, [PM{iPM}, 'On'] );
	% rCenter = fun.(obj.layers(iPM*2-1).name).Center(spacing);
	rCenter = fun.([obj.layers(iPM*2-1).name(1) 'Center'])(spacing);
	if(iPM <= 1)
		% rSurround = fun.(obj.layers(iPM*2-1).name).Surround(spacing);
		rCenter = fun.([obj.layers(iPM*2-1).name(1) 'Surround'])(spacing);
	else
		% rSurround = fun.(obj.layers(iPM*2-3).name).Surround(spacing) ./ fun.(obj.layers(iPM*2-3).name).Center(spacing) .* fun.(obj.layers(iPM*2-1).name).Center(spacing);
		rSurround = fun.([obj.layers(iPM*2-3).name(1), 'Surround'])(spacing) ./ fun.([obj.layers(iPM*2-3).name(1), 'Center'])(spacing) .* fun.([obj.layers(iPM*2-1).name(1), 'Center'])(spacing);
	end
	plot( ecc, rCenter, 'r-' );
	plot( ecc, rSurround, 'b-' );
    plot( ecc, (rSurround./rCenter).^2, 'g.' )
    plot( ecc, obj.SpatialModel.([PM{iPM} 'CenterPeakSensitivityFunction'])( obj.SpatialModel.([PM{iPM} 'CenterPeakSensitivityParams']), rCenter)./ ...
    		   obj.SpatialModel.([PM{iPM} 'SurroundPeakSensitivityFunction'])( obj.SpatialModel.([PM{iPM} 'SurroundPeakSensitivityParams']), rSurround), 'c.' );
	plot( obj.SpatialModel.([PM{iPM} 'CenterData'])('size').eccDegs, obj.SpatialModel.([PM{iPM} 'CenterData'])('size').radiusDegs, 'ro' );
	plot( obj.SpatialModel.([PM{iPM} 'SurroundData'])('size').eccDegs, obj.SpatialModel.([PM{iPM} 'SurroundData'])('size').radiusDegs, 'bo' );

	subplot(2,2,2+iPM); hold on;
	plot( rCenter, obj.SpatialModel.([PM{iPM} 'CenterPeakSensitivityFunction'])( obj.SpatialModel.([PM{iPM} 'CenterPeakSensitivityParams']), rCenter), 'r-' );
	plot( rSurround, obj.SpatialModel.([PM{iPM} 'SurroundPeakSensitivityFunction'])( obj.SpatialModel.([PM{iPM} 'SurroundPeakSensitivityParams']), rSurround), 'b-' );
	plot( obj.SpatialModel.([PM{iPM} 'CenterData'])('sensitivity').radiusDegs, obj.SpatialModel.([PM{iPM} 'CenterData'])('sensitivity').peakSensitivity, 'ro' );
	plot( obj.SpatialModel.([PM{iPM} 'SurroundData'])('sensitivity').radiusDegs, obj.SpatialModel.([PM{iPM} 'SurroundData'])('sensitivity').peakSensitivity, 'bo' );
end

%%
for(iL = 1 : 4)
	[~, spacing] = WatsonRGCModel.RFSpacingDensity(encoder.layers(iL).locations, encoder.layers(iL).name);
	% cR = num2cell(fun.(encoder.layers(iL).name).Center(spacing));
	cR = num2cell(fun.([encoder.layers(iL).name(1) 'Center'])(spacing));
	if(iL <= 2)
		% sR = num2cell(fun.(encoder.layers(iL).name).Surround(spacing));
		sR = num2cell(fun.([encoder.layers(iL).name(1) 'Surround'])(spacing));
	else
		% sR = num2cell(fun.(encoder.layers(iL-2).name).Surround(spacing) ./ fun.(encoder.layers(iL-2).name).Center(spacing) .* fun.(encoder.layers(iL).name).Center(spacing));
		sR = num2cell(fun.([encoder.layers(iL-2).name(1) 'Surround'])(spacing) ./ fun.([encoder.layers(iL-2).name(1) 'Center'])(spacing) .* fun.([encoder.layers(iL).name(1) 'Center'])(spacing));
	end
% 	cR = num2cell(encoder.SpatialModel.([encoder.layers(iL).name(1) 'CenterRadiusFunction'])( encoder.SpatialModel.([encoder.layers(iL).name(1) 'CenterRadiusParams']), [encoder.layers(iL).sRFParams.eccDegs]));
% 	sR = num2cell(encoder.SpatialModel.([encoder.layers(iL).name(1) 'SurroundRadiusFunction'])( encoder.SpatialModel.([encoder.layers(iL).name(1) 'SurroundRadiusParams']), [encoder.layers(iL).sRFParams.eccDegs]));
	
	if(iL <= 4)
		cPS = num2cell(encoder.SpatialModel.([encoder.layers(iL).name(1) 'CenterPeakSensitivityFunction'])(encoder.SpatialModel.([encoder.layers(iL).name(1) 'CenterPeakSensitivityParams']), [cR{:}]));
		sPS = num2cell(encoder.SpatialModel.([encoder.layers(iL).name(1) 'SurroundPeakSensitivityFunction'])(encoder.SpatialModel.([encoder.layers(iL).name(1) 'SurroundPeakSensitivityParams']), [sR{:}]));
	else
		cPS = num2cell(encoder.SpatialModel.([encoder.layers(iL).name(1) 'CenterPeakSensitivityFunction'])([2.4036, -1.8920], [cR{:}]));
		sPS = num2cell(encoder.SpatialModel.([encoder.layers(iL).name(1) 'SurroundPeakSensitivityFunction'])([2.4036, -1.8920], [sR{:}]));
	end
	
	[encoder.layers(iL).sRFParams.centerRadii] = cR{:};
	[encoder.layers(iL).sRFParams.surroundRadii] = sR{:};
	[encoder.layers(iL).sRFParams.centerPeakSensitivities] = cPS{:};
	[encoder.layers(iL).sRFParams.surroundPeakSensitivities] = sPS{:};
end
encoder.DisplayExampleCells(LFR, time, conditions, 'saccadeOff', trials, trialsIdx, cellIdx);

return;



%% run multiple classifiers
folder = '../../Data/figures/Real ratio cells, avContrast=0.5, 10contrasts, same 100 trials; Hs=1';

classifiers = { 'thresholding-uni', ...
				'thresholding-sw_uni-tw_dprime', ...
				'thresholding-sw_sine-tw_uni' };
for(k = 1 : size(classifiers,2)-1)
    fprintf('Processing: %s ...\n', classifiers{k});
    if( exist( fullfile(folder, classifiers{k}, 'PerformanceData.mat') ) )
		delete( fullfile(folder, classifiers{k}, 'PerformanceData.mat') );
	end
	[tTicks, accTPRm, accTPRsd, Thresholds, ThresholdsSTD] = encoder.ContrastDetection(classifiers{k}, 0.625, 0.1, folder, fullfile(folder, classifiers{k}));
end

return;


%%
% Eccs = unique([conditions.eccentricity]);
% nEccs = size(Eccs,2);
% SFs = [0, 2, 10];
% contrasts = [0, 0.5, 0.5];
% colors = { [1 0 0], [0 0 1], [1 0 1], [0 1 1] };

iTrial = 7;

for( iCont = 1 : size(contrasts,2) )
	figure( 'NumberTitle', 'off', 'color', 'w', 'name', sprintf('Spatial Dynamics: Contrast=%0.2f, SF=%d', contrasts(iCont), SFs(iCont)) );
	pause(0.1);
	jf = get(handle(gcf),'javaframe');
	jf.setMaximized(1);
	pause(0.5);

	for( iEcc = 1 : size(Eccs,2) )
		for( iL = 1 : 4 )
			subplot(nEccs, 4, (iEcc-1)*4 + iL); hold on;
			set( gca, 'fontsize', 20, 'linewidth', 2 );

			iCond = find( [conditions.contrast] == contrasts(iCont) & [conditions.eccentricity] == Eccs(iEcc) & [conditions.sf] == SFs(iCont) );

			tStart = 0;

			for( iTick = find(time{1} == tStart) : size(time{1},2)-400 )
				eyeIdx = time{1}(iTick) + trials(trialsIdx{1}(iTrial)).(conditions(iCond).alignEvent);
				plot( encoder.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}), 1) + trials(trialsIdx{1}(iTrial)).x.position(eyeIdx)/60, max(0, LFR{iCond,iL}(:,iTick,iTrial)), '.', 'color', colors{iL} );
            end
            drawnow;
		end
	end

end

return;



%% Video4ExampleCells
% Eccs = unique([conditions.eccentricity]);
% nEccs = size(Eccs,2);
% SFs = [0, 2, 10];
% contrasts = [0, 0.5, 0.5];
% colors = { [1 0 0], [0 0 1], [1 0 1], [0 1 1] };

% iTrial = 1;
% fr = reshape( cat( 1, LFR{ [conditions.contrast] == contrasts(1) | [conditions.contrast] == contrasts(2), : } ), 1, [] );
% fr(fr<0) = 0;
% fr(isnan(fr)) = [];
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

	figure( 'NumberTitle', 'off', 'color', 'k', 'name', sprintf('Movie: Contrast=%0.2f, SF=%d', contrasts(iCont), SFs(iCont)) );
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
						iCond = find( [conditions.contrast] == contrasts(iCont) & [conditions.eccentricity] == Eccs(iEcc) & [conditions.sf] == SFs(iCont) );
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
            iCond = find( [conditions.contrast] == contrasts(iCont) & [conditions.eccentricity] == Eccs(iEcc) & [conditions.sf] == SFs(iCont) );
			
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

return

%%
% eccs = [0, 4, 8, 12];
% for( iL = size(encoder.layers,2) : -1 : 1 )
%     index = encoder.layers(iL).locations(:,1) >= 0 & abs( encoder.layers(iL).locations(:,2) ) > 0.5;		% right visual field and vertically outside +-0.5 deg
%     for( iEcc = size(eccs,2) : -1 : 1 )
%         d2 = sum( encoder.layers(iL).locations.^2, 2 );
%         if( eccs(iEcc) == 0 )
%             idx = find( d2 <= trials(1).gratingWidth^2 & index );
%         else
%             w = trials(1).gratingWidth/2;
%             idx = find( (eccs(iEcc)-w)^2 <= d2 & d2 <= (eccs(iEcc)+w)^2 & index );
%         end
%         [a,r] = cart2pol(encoder.layers(iL).locations(idx,1), encoder.layers(iL).locations(idx,2));
%         aMax = min(a(a>0));
%         aMin = max(a(a<0));
%         [a,r] = cart2pol(encoder.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),1), encoder.layers(iL).locations(cellIdx{iL,iEcc}(1:nCellsUsed{iL,iEcc}),2));
%         cell2AnalyzeIdx{iL,iEcc} = find( aMin < a & a < aMax );
%         nCells2Analyze(iL,iEcc) = size(cell2AnalyzeIdx{iL,iEcc}, 1);
%     end
% end
% save('F:\Post Saccadic Dynamics Modeling\Data\figures\Real ratio cells, avContrast=0.5, 10contrasts, same 100 trials; Hs=1\LFR.mat', '-v7.3', 'cellIdx', 'cell2AnalyzeIdx', 'conditions', 'LFR', 'nCells', 'nCellsUsed', 'nCells2Analyze', 'time', 'trials', 'trialsIdx');
%%
[tTicks, accTPRm, accTPRsd, Thresholds, ThresholdsSTD] = encoder.ContrastDetection('thresholding-uni', 0.625, 0.1, '../../Data/figures/Real ratio cells, avContrast=0.5, 10contrasts, same 100 trials; Hs=1', '../../Data/figures/Real ratio cells, avContrast=0.5, 10contrasts, same 100 trials; Hs=1/thresholding-uni');
[tTicks, accTPRm, accTPRsd, Thresholds, ThresholdsSTD] = encoder.ContrastDetection('thresholding-sw_sine-tw_uni', 0.625, 0.1, '../../Data/figures/Real ratio cells, avContrast=0.5, 10contrasts, same 100 trials; Hs=1', '../../Data/figures/Real ratio cells, avContrast=0.5, 10contrasts, same 100 trials; Hs=1/thresholding-sw_sine-tw_uni');
[tTicks, accTPRm, accTPRsd, Thresholds, ThresholdsSTD] = encoder.ContrastDetection('thresholding-sw_uni-tw_pval', 0.625, 0.1, '../../Data/figures/Real ratio cells, avContrast=0.5, 10contrasts, same 100 trials; Hs=1', '../../Data/figures/Real ratio cells, avContrast=0.5, 10contrasts, same 100 trials; Hs=1/thresholding-sw_uni-tw_pval');
[tTicks, accTPRm, accTPRsd, Thresholds, ThresholdsSTD] = encoder.ContrastDetection('idealobserver', 0.75, 0.5, '../../Data/figures/Real ratio cells, avContrast=0.5, 10contrasts, same 100 trials; Hs=1', '../../Data/figures/Real ratio cells, avContrast=0.5, 10contrasts, same 100 trials; Hs=1/idealobserver');
return;

%%
for k = 1 : 10
    data = load( sprintf( '../../Data/figures/Real ratio cells, avContrast=0.5, contrast=%.2f, same 30 trials; Hs=1/Real ratio cells, avContrast=0.5, contrast=%.2f, same 30 trials; Hs=1.mat', contrasts(k), contrasts(k) ) );
    for m = 1 : size(data.conditions,2)
        if(~data.conditions(m).present), continue; end
        iCond = [conditions.contrast] == contrasts(k) & [conditions.eccentricity] == data.conditions(m).eccentricity & [conditions.sf] == data.conditions(m).sf;
        LFR{iCond,1} = cat(3, [nan(size(data.LFR{m,1},1),20,30), single(data.LFR{m,1})], LFR{iCond,1});
        LFR{iCond,2} = cat(3, [nan(size(data.LFR{m,2},1),20,30), single(data.LFR{m,2})], LFR{iCond,2});
        LFR{iCond,3} = cat(3, [nan(size(data.LFR{m,3},1),20,30), single(data.LFR{m,3})], LFR{iCond,3});
        LFR{iCond,4} = cat(3, [nan(size(data.LFR{m,4},1),20,30), single(data.LFR{m,4})], LFR{iCond,4});
    end
end
% save('F:\Post Saccadic Dynamics Modeling\Data\figures\Real ratio cells, avContrast=0.5, 10contrasts, same 100 trials; Hs=1\LFR.mat', '-v7.3', 'cellIdx', 'cell2AnalyzeIdx', 'conditions', 'LFR', 'nCells', 'nCellsUsed', 'nCells2Analyze', 'time', 'trials', 'trialsIdx' );
return;

%%
absent = max(0,cat(1,LFR{1,:}));
present2 = max(0,cat(1,LFR{4,:}));
present10 = max(0,cat(1,LFR{7,:}));

mAbsent = mean(mean(absent,1),3);
sdAbsent = std(mean(absent,1),0,3);
mPresent2 = mean(mean(present2,1),3);
sdPresent2 = std(mean(present2,1),0,3);
mPresent10 = mean(mean(present10,1),3);
sdPresent10 = std(mean(present10,1),0,3);

%%figure; hold on;
subplot(2,2,1); hold on; cla;
h(1) = plot( time{1}, mAbsent, 'b', 'linewidth', 1, 'displayname', 'Absent' );
h(2) = plot( time{1}, mPresent2, 'g', 'linewidth', 1, 'displayname', '2 cpd' );
h(3) = plot( time{1}, mPresent10, 'r', 'linewidth', 1, 'displayname', '10 cpd' );
m = mAbsent;
sem = sdAbsent / sqrt(30);
fill( [time{1}, time{1}(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], 'b', 'LineStyle', 'none', 'FaceAlpha', 0.4 );
m = mPresent2;
sem = sdPresent2 / sqrt(30);
fill( [time{1}, time{1}(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], 'g', 'LineStyle', 'none', 'FaceAlpha', 0.4 );
m = mPresent10;
sem = sdPresent10 / sqrt(30);
fill( [time{1}, time{1}(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], 'r', 'LineStyle', 'none', 'FaceAlpha', 0.4 );
legend(h, 'location', 'northeast');
set(gca, 'xlim', [-50 500], 'linewidth', 2, 'fontsize', 18)
title('Population Mean');
xlabel('Time from saccade off (ms)');
ylabel('Firing rate');

%%
subplot(2,2,2); hold on; cla;
colors = {'b', 'g', 'r'};
iTick = find(time{1} == 0);
accFR = {absent(:,iTick:end,:), present2(:,iTick:end,:), present10(:,iTick:end,:)};
for(k = 1 : 3)
	for( t = 2 : size(accFR{1},2) )
		accFR{k}(:,t,:) = accFR{k}(:,t-1,:) + accFR{k}(:,t,:);
	end
end
for(k = 1 : 3)
	m = mean(mean(accFR{k},1),3);
	sem = std(mean(accFR{k},1),0,3) / sqrt(size(accFR{k},3));
	plot( time{1}(iTick:end), m, 'color', colors{k}, 'linewidth', 1 );
	fill( [time{1}(iTick:end), time{1}(end:-1:iTick)], [m-sem, m(end:-1:1)+sem(end:-1:1)], colors{k}, 'LineStyle', 'none', 'FaceAlpha', 0.4 );
end
set(gca, 'xlim', [0 500], 'linewidth', 2, 'fontsize', 18)
title('Accumulated Population Mean');
xlabel('Time from saccade off (ms)');
ylabel('Accumulated firing rate');

%%
subplot(2,2,3); hold on; cla;
for(k = 1 : 3)
	meanFR{k} = accFR{k} ./ (1:size(accFR{k},2));
end
for(k = 1 : 3)
	m = mean(mean(meanFR{k},1),3);
	sem = std(mean(meanFR{k},1),0,3) / sqrt(size(meanFR{k},3));
	plot( time{1}(iTick:end), m, 'color', colors{k}, 'linewidth', 1 );
	fill( [time{1}(iTick:end), time{1}(end:-1:iTick)], [m-sem, m(end:-1:1)+sem(end:-1:1)], colors{k}, 'LineStyle', 'none', 'FaceAlpha', 0.4 );
end
set(gca, 'xlim', [0 500], 'linewidth', 2, 'fontsize', 18)
title('Average of Population Mean');
xlabel('Time from saccade off (ms)');
ylabel('Average firing rate');


%%
subplot(2,2,4); hold on; cla;
for(k = 2 : 3)
	mDiffAccFR{k} = mean(mean(accFR{k},1),3) - mean(mean(accFR{1},1),3);
	sdDiffAccFR{k} = sqrt( std(mean(accFR{k},1),0,3).^2 + std(mean(accFR{1},1),0,3).^2 );
end
for(k = 2 : 3)
	m = mDiffAccFR{k};
	sem = sdDiffAccFR{k};
	plot( time{1}(iTick:end), m, 'color', colors{k}, 'linewidth', 1 );
	fill( [time{1}(iTick:end), time{1}(end:-1:iTick)], [m-sem, m(end:-1:1)+sem(end:-1:1)], colors{k}, 'LineStyle', 'none', 'FaceAlpha', 0.4 );
end
set(gca, 'xlim', [0 500], 'linewidth', 2, 'fontsize', 18)
title('Difference from Absent for Accumulated Population Mean');
xlabel('Time from saccade off (ms)');
ylabel('Accumulated firing rate');
return

%%
iCond = 4; 1; 7;
Eccs = [0 4 8];
nTrials = 30;
nNeurons = 3;
obj = encoder;
contrast = 0.5;
iEcc = find( conditions(iCond).eccentricity == Eccs );
idx = trialsIdx{iCond};    %fprintf('nTrials = %d\n', size(idx,2)); continue;
idx = idx( 1 : min(nTrials,end) );
tMax = round(max( [trials(idx).saccadeOff] + round(conditions(iCond).duration/1000*[trials(idx).sRate]) - [trials(idx).(alignEvent)] ));	% in samples
tMin = round(min( [trials(idx).saccadeOn] - round(0.150*[trials(idx).sRate]) - [trials(idx).(alignEvent)] ));				% in samples
k = 1;
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
avContrast = trials(idx(k)).backgroundContrast + contrast * conditions(iCond).present;

iL = 3;
tic;
for( m = 1 : 3 )
	xIdx = min(x(eyeIdx{m})) + min(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nNeurons),1)) - 4 <= inputX & inputX <= max(x(eyeIdx{m})) + max(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nNeurons),1)) + 4; %720 : 1460;
    yIdx = min(y(eyeIdx{m})) + min(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nNeurons),2)) - 4 <= inputY & inputY <= max(y(eyeIdx{m})) + max(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nNeurons),2)) + 4; %300 : 780;

	[sfr{iCond,iL}(:,eyeIdx{m},k), sfr_c{iCond,iL}(:,eyeIdx{m},k), sfr_s{iCond,iL}(:,eyeIdx{m},k)] = obj.SpatialModel.LinearResponse( stimulus(yIdx,xIdx), inputX(xIdx), inputY(yIdx), x(eyeIdx{m}), y(eyeIdx{m}), obj.layers(iL).sRFParams(cellIdx{iL,iEcc}(1:nNeurons)), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nNeurons),1), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nNeurons),2));
end
if( lower(obj.layers(iL).name(1)) == 'm' )
	lfr{iCond,iL}(:,:,k) = obj.TemporalModel.LinearResponse( obj.layers(iL).name, obj.layers(iL).tRFParams(cellIdx{iL,iEcc}(1:nNeurons)), avContrast, trials(idx(k)).sRate, sfr{iCond,iL}(:,:,k) );
else
	lfr{iCond,iL}(:,:,k) = obj.TemporalModel.LinearResponse( obj.layers(iL).name, obj.layers(iL).tRFParams(cellIdx{iL,iEcc}(1:nNeurons)), avContrast, trials(idx(k)).sRate, sfr_c{iCond,iL}(:,:,k), sfr_s{iCond,iL}(:,:,k) );
end
fprintf('Ecc=%d, SF=%d, Dur=%d, Present=%d, iL=%d, iTrials=%d/%d, t=%f\n', conditions(iCond).eccentricity, conditions(iCond).sf, conditions(iCond).duration, conditions(iCond).present, iL, k, size(idx,2), toc);

return;


idx = trialsIdx{1};
tMin = time{1}(1);
tMax = time{1}(end);
contrast = 0.1;
iConds = [1 5 9]+4;
iL = 1;
% iEcc = find( conditions(iCond).eccentricity == [0, 4, 8, 12] );
nCells = 15;
alignEvent = 'flashOn';
dataFolder = '../../Data';
figure('color', 'w');
col = {'b', 'g', 'r'};
for( iEcc = 1 : 4 )
	iConds = [0 4 8] + iEcc;
	subplot(2,2,iEcc); hold on;
	for( iCond = iConds )
		for( k = 1 : 15 )
			x  = trials(idx(k)).x.position / 60;
			y  = trials(idx(k)).y.position / 60;

			% before grating
			eyeIdx1 = tMin + trials(idx(k)).(alignEvent) : trials(idx(k)).flashOn-1;
			s1 = size(eyeIdx1,2);

			% during grating
			eyeIdx2 = trials(idx(k)).flashOn : min( size(x,2), tMax + trials(idx(k)).(alignEvent) );
			s2 = size(eyeIdx2,2);

			[noise, inputX, inputY] = encoder.LoadNoise( fullfile(dataFolder, trials(idx(k)).backgroundImage), trials(idx(k)).pixelAngle/60 );
			grating = encoder.GenerateGrating( conditions(iCond).eccentricity, trials(idx(k)).gratingWidth, conditions(iCond).sf, trials(idx(k)).phase, trials(idx(k)).pixelAngle/60 );
			noise = noise * trials(idx(k)).backgroundContrast;
			grating = grating * contrast * conditions(iCond).present;
			bg = 1;

			avContrast = ((s1+s2)*trials(idx(k)).backgroundContrast + s2*contrast*conditions(iCond).present) / (s1+s2);
			% contrastF = [ zeros(1,s1), ones(1,s2) * contrast ] + trials(idx(k)).backgroundContrast;
			contrastF(:,:,k) = encoder.ComputeContrasts( noise+grating+bg,  inputX, inputY, encoder.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),1), encoder.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),2), [encoder.layers(iL).sRFParams(cellIdx{iL,iEcc}(1:nCells)).surroundRadii], x([eyeIdx1,eyeIdx2]), y([eyeIdx1,eyeIdx2]) );
		end
		plot( tMin:tMax, mean(contrastF-(contrastF(:,1,:)),3)', 'color', col{iCond==iConds} );
	end
	title( sprintf('Ecc = %d', conditions(iCond).eccentricity) );
	if(iEcc >= 3), xlabel('Time aligned to flashOn (ms)'); end
	if(mod(iEcc,2)), ylabel('Contrast within 6*std of cell surround'); end
	set(gca, 'fontsize', 20, 'linewidth', 2);
end
return;


%% prepare data
load(['../../Data/Mosaic_POn_Radius20.0deg_maxMovPrctile20.mat']);
rfPositions = rfPositions( abs(rfPositions(:,1))<=2 & abs(rfPositions(:,2))<=2, : );

ck = CronerKaplanRGCModel( 'dataSetToFit', 'raw', 'fitIntercept', false, 'randomizeCenterRadii', true, 'randomizeCenterSensitivities', true, 'randomizeSurroundRadii', true, 'randomizeSurroundSensitivities', true );
PRFParams = ck.SynthesizeRFParams(sqrt(sum(rfPositions.^2,2))','POn');


img = imread('demo.jpg');
img = single(img(end:-1:1,:,1));
inputX = (-275:275)/60;
inputY = (-310:309)/60;

load('ExampleEyeTrace.mat');
eyeX = eyeX/4;
eyeX = eyeX(1501:2500);
eyeY = eyeY(1501:2500);
eyeX = -0.15;
eyeY = -1.88;


%% spatial RF
fr = zeros( size(rfPositions,1), length(eyeX) );
fr_c = fr;
fr_s = fr;

n = 1000;
T = zeros( 1, ceil(size(rfPositions,1)/n) );
for( k = 0 : ceil(size(rfPositions,1)/n)-1 )
	fprintf('%d/%d\n', k, ceil(size(rfPositions,1)/n)-1 );
	idx = k*n+1 : min(k*n+n, size(rfPositions,1));
    tic;
	[fr(idx,:), fr_c(idx,:), fr_s(idx,:)] = ck.LinearResponse(img, inputX, inputY, eyeX, eyeY, PRFParams(idx), rfPositions(idx,1), rfPositions(idx,2));
    T(k+1) = toc;
end

%%	temporal RF
for( k = size(fr,1) : -1 : 1 )
	tfr_c(k,:) = RGC.TemporalLinearModel( fr_c(k,:), 'P', 'center' );
	tfr_s(k,:) = RGC.TemporalLinearModel( fr_s(k,:), 'P', 'surround' );
end
tfr = tfr_c + tfr_s;

%% Display activity
t = 1;
d = 1/60/2;
[x, y] = meshgrid( -2:d:2, -2:d:2 );
frMap = zeros(size(x));
counts = zeros(size(x));
idx = round( (rfPositions(:,1)' - -2 ) / d ) + 1;
idy = round( (rfPositions(:,2)' - -2 ) / d ) + 1;
for( k = find( idx >= 1 & idx <= size(x,1) & idy >= 1 & idy <= size(y,2) ) )
	frMap(idy(k),idx(k)) = frMap(idy(k),idx(k)) + fr(k,t);%max(0,fr(k,t));
	counts(idy(k),idx(k)) = counts(idy(k),idx(k)) + 1;
end
frMap(:) = max( frMap(:), prctile(frMap(:), 0.5) );
frMap(:) = min( frMap(:), prctile(frMap(:), 99.5) );
frMap = (frMap - min(frMap(:))) ./ (max(frMap(:)) - min(frMap(:)));
imshow( cat(3, frMap, zeros(size(frMap,1),size(frMap,2),2)) );
set( gca, 'ydir', 'normal' );