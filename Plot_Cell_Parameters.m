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
subplot(2,2,3); hold on; h = [];
%eccs = unique([obj.activityParams.conditions.eccentricity]);
for(k = 1:4)
	% iid = cat(1, obj.layers(k).idxExampleCells{eccs == 0 | eccs == 4 | eccs == 8});
	% iid = cat(1, obj.layers(k).idxExampleCells{:});
	iid = cat(1, obj.layers(k).locations(:,1) > 0 & abs(obj.layers(k).locations(:,2)) <= 1);
	x = [obj.layers(k).sRFParams(iid).temporalEccDegs];
	y = [obj.layers(k).sRFParams(iid).centerRadii];
	x2 = x.^2;
	[x, ii] = sort(x);
	y = y(ii);
	x2 = x2(ii);
	N = 51; NL = (N-1) / 2;
	x = [x(NL:-1:1), x, x(end:-1:end-NL+1)];
	y = [y(NL:-1:1), y, y(end:-1:end-NL+1)];
    x2 = [x2(NL:-1:1), x2, x2(end:-1:end-NL+1)];
	x = conv(x, ones(1,N)/N, 'same');
	y = conv(y, ones(1,N)/N, 'same');
	x2 = conv(x2, ones(1,N)/N, 'same');
	x = x(NL+1:end-NL);
	y = y(NL+1:end-NL);
	x2 = x2(NL+1:end-NL);
	sd = x2 - x.^2;
	h(end+1) = plot(x, y, '-', 'color', colors{k}, 'linewidth', 2, 'displayname', [obj.layers(k).name ' center']);
	fill([x x(end:-1:1)], [y-sd y(end:-1:1)+sd(end:-1:1)], 'k', 'LineStyle', 'none', 'FaceColor', colors{k}, 'FaceAlpha', 0.5);

	y = [obj.layers(k).sRFParams(iid).surroundRadii];
	x2 = x.^2;
	[x, ii] = sort(x);
	y = y(ii);
	x2 = x2(ii);
	N = 51; NL = (N-1) / 2;
	x = [x(NL:-1:1), x, x(end:-1:end-NL+1)];
	y = [y(NL:-1:1), y, y(end:-1:end-NL+1)];
    x2 = [x2(NL:-1:1), x2, x2(end:-1:end-NL+1)];
	x = conv(x, ones(1,N)/N, 'same');
	y = conv(y, ones(1,N)/N, 'same');
	x2 = conv(x2, ones(1,N)/N, 'same');
	x = x(NL+1:end-NL);
	y = y(NL+1:end-NL);
	x2 = x2(NL+1:end-NL);
	sd = x2 - x.^2;
	h(end+1) = plot(x, y, '--', 'color', colors{k}, 'linewidth', 2, 'displayname', [obj.layers(k).name ' surround']);
	fill([x x(end:-1:1)], [y-sd y(end:-1:1)+sd(end:-1:1)], 'k', 'LineStyle', 'none', 'FaceColor', colors{k}, 'FaceAlpha', 0.5);
end
legend(h, 'fontsize', 16, 'position', [0.1365 0.3236 0.0979 0.2297]);
xlabel('Eccentricity (\circ)');
ylabel('Radius (\circ)');
title('Radius');
set(gca, 'xlim', [-0.2 15], 'ylim', [-0.05 1.5], 'fontsize', 20, 'linewidth', 2);

% peak sensitivity
subplot(2,2,4); hold on;
%eccs = unique([obj.activityParams.conditions.eccentricity]);
for(k = 1:4)
	% iid = cat(1, obj.layers(k).idxExampleCells{eccs == 0 | eccs == 4 | eccs == 8});
	% iid = cat(1, obj.layers(k).idxExampleCells{:});
	iid = cat(1, obj.layers(k).locations(:,1) > 0 & abs(obj.layers(k).locations(:,2)) <= 1);
	x = [obj.layers(k).sRFParams(iid).temporalEccDegs];
	y = abs([obj.layers(k).sRFParams(iid).centerPeakSensitivities]);
	x2 = x.^2;
	[x, ii] = sort(x);
	y = y(ii);
	x2 = x2(ii);
	N = 51; NL = (N-1) / 2;
	x = [x(NL:-1:1), x, x(end:-1:end-NL+1)];
	y = [y(NL:-1:1), y, y(end:-1:end-NL+1)];
    x2 = [x2(NL:-1:1), x2, x2(end:-1:end-NL+1)];
	x = conv(x, ones(1,N)/N, 'same');
	y = conv(y, ones(1,N)/N, 'same');
	x2 = conv(x2, ones(1,N)/N, 'same');
	x = x(NL+1:end-NL);
	y = y(NL+1:end-NL);
	x2 = x2(NL+1:end-NL);
	sd = x2 - x.^2;
	plot(x, y, '-', 'color', colors{k}, 'linewidth', 2, 'displayname', [obj.layers(k).name ' center']);
	fill([x x(end:-1:1)], [y-sd y(end:-1:1)+sd(end:-1:1)], 'k', 'LineStyle', 'none', 'FaceColor', colors{k}, 'FaceAlpha', 0.5);

	x = [obj.layers(k).sRFParams(iid).temporalEccDegs];
	y = abs(abs([obj.layers(k).sRFParams(iid).surroundPeakSensitivities]));
	x2 = x.^2;
	[x, ii] = sort(x);
	y = y(ii);
	x2 = x2(ii);
	N = 51; NL = (N-1) / 2;
	x = [x(NL:-1:1), x, x(end:-1:end-NL+1)];
	y = [y(NL:-1:1), y, y(end:-1:end-NL+1)];
    x2 = [x2(NL:-1:1), x2, x2(end:-1:end-NL+1)];
	x = conv(x, ones(1,N)/N, 'same');
	y = conv(y, ones(1,N)/N, 'same');
	x2 = conv(x2, ones(1,N)/N, 'same');
	x = x(NL+1:end-NL);
	y = y(NL+1:end-NL);
	x2 = x2(NL+1:end-NL);
	sd = x2 - x.^2;
	plot(x, y, '--', 'color', colors{k}, 'linewidth', 2, 'displayname', [obj.layers(k).name ' surround']);
	fill([x x(end:-1:1)], [y-sd y(end:-1:1)+sd(end:-1:1)], 'k', 'LineStyle', 'none', 'FaceColor', colors{k}, 'FaceAlpha', 0.5);
end
xlabel('Eccentricity (\circ)');
ylabel('Peak sensitivity');
title('Peak Sensitivity');
set(gca, 'xlim', [-0.2 15], 'yscale', 'log', 'fontsize', 20, 'linewidth', 2);