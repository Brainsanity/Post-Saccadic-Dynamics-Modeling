% Compare radius predicted by CronerKaplanRGCModel (fitted with raw data and no intercept) with spacing predicted by WatsonRGCModel

% ck = CronerKaplanRGCModel( 'dataSetToFit', 'raw', 'fitIntercept', false );

figure( 'NumberTitle', 'off', 'name', 'Predicted (Raw noIntercept) Radius VS Spacing Predicted by WatsonModel', 'color', 'w' );
lineWidth = 2;
fontSize = 20;
PM = {'P', 'M'};
ecc = 0:0.01:40;
colors = { 'r', 'm', ...
		  [0 1 0], 			[0 0 1], ...
		  [0.5 1 0.5], 		[0.5 0.5 1], ...
		  [0 0.5 0], 		[0 0 0.5], ...
		  [0.25 0.75 0.25], [0.25 0.25 0.75] };
for( iPlot = 2:-1:1 )
	subplot(1,2,iPlot); hold on;
	h = [];
	h(1) = plot( ecc, ck.([PM{iPlot} 'CenterRadiusFunction'])( ck.([PM{iPlot} 'CenterRadiusParams']), ecc ), 'r--', 'lineWidth', lineWidth, 'displayName', 'CK Center' );
	h(2) = plot( ecc, ck.([PM{iPlot} 'SurroundRadiusFunction'])( ck.([PM{iPlot} 'SurroundRadiusParams']), ecc ), 'm--', 'lineWidth', lineWidth, 'displayName', 'CK Surround' );
	for( k = 1 : 4 )
		[~, spacingOn] = WatsonRGCModel.RFSpacingDensityMeridian( ecc, WatsonRGCModel.enumeratedMeridianNames{k}, [PM{iPlot} 'On'] );
		[~, spacingOff] = WatsonRGCModel.RFSpacingDensityMeridian( ecc, WatsonRGCModel.enumeratedMeridianNames{k}, [PM{iPlot} 'Off'] );
		h(k*2+1) = plot( ecc, spacingOn, '-', 'color', colors{k*2+1}, 'lineWidth', lineWidth, 'displayName', ['On ' WatsonRGCModel.enumeratedMeridianNames{k}(1:end-9)] );
		h(k*2+2) = plot( ecc, spacingOff, '-', 'color', colors{k*2+2}, 'lineWidth', lineWidth, 'displayName', ['Off ' WatsonRGCModel.enumeratedMeridianNames{k}(1:end-9)] );
	end
	k = 1;
	[~, spacingOn] = WatsonRGCModel.RFSpacingDensityMeridian( ecc, WatsonRGCModel.enumeratedMeridianNames{k}, [PM{iPlot} 'On'] );
	[~, spacingOff] = WatsonRGCModel.RFSpacingDensityMeridian( ecc, WatsonRGCModel.enumeratedMeridianNames{k}, [PM{iPlot} 'Off'] );
% 	h(end+1) = plot( ecc, fun.([PM{iPlot} 'On']).Center(spacingOn), 'r-.', 'lineWidth', 2, 'displayName', 'a*Spacing: On Center' );
% 	h(end+1) = plot( ecc, fun.([PM{iPlot} 'Off']).Center(spacingOff), 'm-.', 'lineWidth', 2, 'displayName', 'a*Spacing: Off Center' );
% 	h(end+1) = plot( ecc, fun.([PM{iPlot} 'On']).Surround(spacingOn), 'r:', 'lineWidth', 2, 'displayName', 'a*Spacing: On Surround' );
% 	h(end+1) = plot( ecc, fun.([PM{iPlot} 'Off']).Surround(spacingOff), 'm:', 'lineWidth', 2, 'displayName', 'a*Spacing: Off Surround' );
	h(end+1) = plot( ecc, ck.([PM{iPlot} 'CenterSpacing2Radius']).fun(spacingOn), 'r-.', 'lineWidth', 2, 'displayName', 'Radius ~ Spacing: On Center' );
	h(end+1) = plot( ecc, ck.([PM{iPlot} 'CenterSpacing2Radius']).fun(spacingOff), 'm-.', 'lineWidth', 2, 'displayName', 'Radius ~ Spacing: Off Center' );
	h(end+1) = plot( ecc, ck.([PM{iPlot} 'SurroundSpacing2Radius']).fun(spacingOn), 'r:', 'lineWidth', 2, 'displayName', 'Radius ~ Spacing: On Surround' );
	h(end+1) = plot( ecc, ck.([PM{iPlot} 'SurroundSpacing2Radius']).fun(spacingOff), 'm:', 'lineWidth', 2, 'displayName', 'Radius ~ Spacing: Off Surround' );
	title( [PM{iPlot} ' Cell'] );
	xlabel('Eccentricity (deg)');
	set( gca, 'lineWidth', lineWidth, 'fontSize', fontSize );
end
ylabel('Radius or Spacing (deg)');
legend( h([1:2, 3:2:end, 4:2:end]), 'location', 'northWest' );
