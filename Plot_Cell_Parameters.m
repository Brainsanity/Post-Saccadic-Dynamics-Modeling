figure('NumberTitle', 'off', 'name', 'Cell Parameters', 'color', 'w'); pause(0.1); jf = get(handle(gcf),'javaframe'); jf.setMaximized(1); pause(1);
colors = {'r', 'b', 'm', 'c'};

% density
subplot(2,2,1); hold on;
ecc_ = 0:0.001:15;
for(k = 1:4)
	[density_, spacing_] = WatsonRGCModel.RFSpacingDensityMeridian( ecc_, WatsonRGCModel.enumeratedMeridianNames{1}, obj.layers(k).name );
	plot(ecc_, density_, '-', 'color', colors{k}, 'linewidth', 2, 'displayname', obj.layers(k).name);
end
legend('location', 'northeast');
% xlabel('Eccentricity (\circ)');
ylabel('Density (deg^{-2})');
title('Density');
set(gca, 'yscale', 'log', 'fontsize', 20, 'linewidth', 2);

% ratio
subplot(2,2,2); hold on;
plot(ecc_, WatsonRGCModel.RFSpacingDensityMeridian( ecc_, WatsonRGCModel.enumeratedMeridianNames{1}, 'M' ) ./ WatsonRGCModel.RFSpacingDensityMeridian( ecc_, WatsonRGCModel.enumeratedMeridianNames{1}, 'P' ), 'k-', 'linewidth', 2, 'displayname', 'M/P ratio');
plot(ecc_, WatsonRGCModel.OnOffRatio(ecc_, 'POn') ./ WatsonRGCModel.OnOffRatio(ecc_, 'POff'), 'r-', 'linewidth', 2, 'displayname', 'P On/Off ratio');
plot(ecc_, WatsonRGCModel.OnOffRatio(ecc_, 'MOn') ./ WatsonRGCModel.OnOffRatio(ecc_, 'MOff'), 'm-', 'linewidth', 2, 'displayname', 'M On/Off ratio');
legend('location', 'southeast');
% xlabel('Eccentricity (\circ)');
ylabel('Ratio');
title('Ratio');
set(gca, 'fontsize', 20, 'linewidth', 2);

% ratius
subplot(2,2,3); hold on;
eccs = unique([obj.activityParams.conditions.eccentricity]);
for(k = 1:4)
	% iid = cat(1, obj.layers(k).idxExampleCells{eccs == 0 | eccs == 4 | eccs == 8});

iid = cat(1, obj.layers(k).idxExampleCells{:});	plot([obj.layers(k).sRFParams(iid).temporalEccDegs], [obj.layers(k).sRFParams(iid).centerRadii], '+', 'color', colors{k}, 'displayname', [obj.layers(k).name ' center']);
	plot([obj.layers(k).sRFParams(iid).temporalEccDegs], [obj.layers(k).sRFParams(iid).surroundRadii], 'o', 'color', colors{k}, 'displayname', [obj.layers(k).name ' surround']);
end
legend('fontsize', 16, 'position', [0.1365 0.3236 0.0979 0.2297]);
xlabel('Eccentricity (\circ)');
ylabel('Radius (\circ)');
title('Radius');
set(gca, 'xlim', [-0.2 15], 'ylim', [-0.05 1.5], 'fontsize', 20, 'linewidth', 2);

% peak sensitivity
subplot(2,2,4); hold on;
eccs = unique([obj.activityParams.conditions.eccentricity]);
for(k = 1:4)
	% iid = cat(1, obj.layers(k).idxExampleCells{eccs == 0 | eccs == 4 | eccs == 8});
	iid = cat(1, obj.layers(k).idxExampleCells{:});
	plot([obj.layers(k).sRFParams(iid).temporalEccDegs], abs([obj.layers(k).sRFParams(iid).centerPeakSensitivities]), '+', 'color', colors{k}, 'displayname', [obj.layers(k).name ' center']);
	plot([obj.layers(k).sRFParams(iid).temporalEccDegs], abs([obj.layers(k).sRFParams(iid).surroundPeakSensitivities]), 'o', 'color', colors{k}, 'displayname', [obj.layers(k).name ' surround']);
end
xlabel('Eccentricity (\circ)');
ylabel('Peak sensitivity');
title('Peak Sensitivity');
set(gca, 'xlim', [-0.2 15], 'yscale', 'log', 'fontsize', 20, 'linewidth', 2);