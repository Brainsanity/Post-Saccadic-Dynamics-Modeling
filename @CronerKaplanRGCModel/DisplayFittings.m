function DisplayFittings()

	hFigDigitized = figure( 'NumberTitle', 'off', 'name', 'Digitized Data', 'color', 'w' );
	hFigP = figure( 'NumberTitle', 'off', 'name', 'P Cell Fittings', 'color', 'w' );
	hFigM = figure( 'NumberTitle', 'off', 'name', 'M Cell Fittings', 'color', 'w' );
	markerSize = 8;
	lineWidth = 2;
	fontSize = 20;
	colors = {'r', 'b', 'm', 'c', 'g', 'y'};

	ck = CronerKaplanRGCModel( 'dataSetToFit', 'paperFormulas', 'fitIntercept', false );
	
	%%%% display digitized data points
	figure(hFigDigitized);
	subplot(2,2,1); hold on;	% Figure 4a
	h(1) = plot( ck.PCenterData('size').eccDegs, ck.PCenterData('size').radiusDegs, 'o', 'color', colors{1}, 'markerSize', markerSize, 'lineWidth', lineWidth, 'displayName', 'P' );
	h(2) = plot( ck.MCenterData('size').eccDegs, ck.MCenterData('size').radiusDegs, '^', 'color', colors{2}, 'markerSize', markerSize, 'lineWidth', lineWidth, 'displayName', 'M' );
	legend( h, 'location', 'northWest' );
	xlabel('Temporal equivalent eccentricity (deg)');
	ylabel('Center radius (deg)');
	title('Center (Figure 4a)');
	set( gca, 'xlim', [0 40], 'ylim', [0 0.3], 'fontSize', fontSize, 'lineWidth', lineWidth );

	subplot(2,2,2); hold on;	% Figure 4b
	h(1) = plot( ck.PSurroundData('size').eccDegs, ck.PSurroundData('size').radiusDegs, 'o', 'color', colors{1}, 'markerSize', markerSize, 'lineWidth', lineWidth, 'displayName', 'P' );
	h(2) = plot( ck.MSurroundData('size').eccDegs, ck.MSurroundData('size').radiusDegs, '^', 'color', colors{2}, 'markerSize', markerSize, 'lineWidth', lineWidth, 'displayName', 'M' );
	legend( h, 'location', 'northWest' );
	xlabel('Temporal equivalent eccentricity (deg)');
	ylabel('Surround radius (deg)');
	title('Surround (Figure 4b)');
	set( gca, 'xlim', [0.1 100], 'ylim', [0.01 10], 'xscale', 'log', 'yscale', 'log', 'fontSize', fontSize, 'lineWidth', lineWidth );

	subplot(2,2,3); hold on;	% Figure 5a
	h(1) = plot( ck.PCenterData('sensitivity').radiusDegs, ck.PCenterData('sensitivity').peakSensitivity, 'o', 'color', colors{1}, 'markerSize', markerSize, 'lineWidth', lineWidth, 'displayName', 'P Center' );
	h(2) = plot( ck.MCenterData('sensitivity').radiusDegs, ck.MCenterData('sensitivity').peakSensitivity, '^', 'color', colors{2}, 'markerSize', markerSize, 'lineWidth', lineWidth, 'displayName', 'M Center' );
	h(3) = plot( ck.PSurroundData('sensitivity').radiusDegs, ck.PSurroundData('sensitivity').peakSensitivity, 'o', 'color', colors{3}, 'markerSize', markerSize, 'lineWidth', lineWidth, 'displayName', 'P Surround' );
	h(4) = plot( ck.MSurroundData('sensitivity').radiusDegs, ck.MSurroundData('sensitivity').peakSensitivity, '^', 'color', colors{4}, 'markerSize', markerSize, 'lineWidth', lineWidth, 'displayName', 'M Surround' );
	legend( h, 'location', 'southWest' );
	xlabel('Radius (deg)');
	ylabel('Peak sensitivity (imp/(s %contrast deg^2))');
	title('Figure 5a');
	set( gca, 'xlim', [0.01 10], 'ylim', [0.001 1000], 'xscale', 'log', 'yscale', 'log', 'fontSize', fontSize, 'lineWidth', lineWidth );


	%%%% Fitting with paper formula; not application for M cells
	%% P Cell Fittings
	figure(hFigP);
	subplot(2,2,1); hold on;
	plot( ck.PCenterData('size').eccDegs, ck.PCenterData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	% hLinesP{1} = plot( 0.01:0.01:40, ck.PCenterRadiusFunction( ck.PCenterRadiusParams, 0.01:0.01:40 ), '--', 'color', colors{1}, 'lineWidth', lineWidth, 'displayName', 'PaperFormula' );

	subplot(2,2,2); hold on;
	plot( ck.PSurroundData('size').eccDegs, ck.PSurroundData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesP{2} = plot( 0.1:0.01:60, ck.PSurroundRadiusFunction( ck.PSurroundRadiusParams, 0.1:0.01:60 ), '--', 'color', colors{1}, 'lineWidth', lineWidth, 'displayName', 'PaperFormula' );

	subplot(2,2,3); hold on;
	plot( ck.PCenterData('sensitivity').radiusDegs, ck.PCenterData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesP{3} = plot( 0.01:0.01:10, ck.PCenterPeakSensitivityFunction( ck.PCenterPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{1}, 'lineWidth', lineWidth, 'displayName', 'PaperFormula' );

	subplot(2,2,4); hold on;
	plot( ck.PSurroundData('sensitivity').radiusDegs, ck.PSurroundData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesP{4} = plot( 0.01:0.01:10, ck.PSurroundPeakSensitivityFunction( ck.PSurroundPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{1}, 'lineWidth', lineWidth, 'displayName', 'PaperFormula' );


	%%%% Fitting with median values from the paper with no intercept
	ck = CronerKaplanRGCModel( 'dataSetToFit', 'medians', 'fitIntercept', false );
	%% P Cell Fittins
	figure(hFigP);
	subplot(2,2,1); hold on;
	plot( ck.PCenterData('size').eccDegs, ck.PCenterData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesP{1} = plot( 0.01:0.01:40, ck.PCenterRadiusFunction( ck.PCenterRadiusParams, 0.01:0.01:40 ), '--', 'color', colors{2}, 'lineWidth', lineWidth, 'displayName', 'Medians noIntercept' );

	subplot(2,2,2); hold on;
	plot( ck.PSurroundData('size').eccDegs, ck.PSurroundData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesP{2}(end+1) = plot( 0.1:0.01:60, ck.PSurroundRadiusFunction( ck.PSurroundRadiusParams, 0.1:0.01:60 ), '--', 'color', colors{2}, 'lineWidth', lineWidth, 'displayName', 'Medians noIntercept' );

	subplot(2,2,3); hold on;
	plot( ck.PCenterData('sensitivity').radiusDegs, ck.PCenterData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesP{3}(end+1) = plot( 0.01:0.01:10, ck.PCenterPeakSensitivityFunction( ck.PCenterPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{2}, 'lineWidth', lineWidth, 'displayName', 'Medians noIntercept' );

	subplot(2,2,4); hold on;
	plot( ck.PSurroundData('sensitivity').radiusDegs, ck.PSurroundData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesP{4}(end+1) = plot( 0.01:0.01:10, ck.PSurroundPeakSensitivityFunction( ck.PSurroundPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{2}, 'lineWidth', lineWidth, 'displayName', 'Medians noIntercept' );

	%% M Cell Fittins
	figure(hFigM);
	subplot(2,2,1); hold on;
	plot( ck.MCenterData('size').eccDegs, ck.MCenterData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesM{1} = plot( 0.01:0.01:40, ck.MCenterRadiusFunction( ck.MCenterRadiusParams, 0.01:0.01:40 ), '--', 'color', colors{2}, 'lineWidth', lineWidth, 'displayName', 'Medians noIntercept' );

	subplot(2,2,2); hold on;
	plot( ck.MSurroundData('size').eccDegs, ck.MSurroundData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesM{2} = plot( 0.1:0.01:60, ck.MSurroundRadiusFunction( ck.MSurroundRadiusParams, 0.1:0.01:60 ), '--', 'color', colors{2}, 'lineWidth', lineWidth, 'displayName', 'Medians noIntercept' );

	subplot(2,2,3); hold on;
	plot( ck.MCenterData('sensitivity').radiusDegs, ck.MCenterData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesM{3} = plot( 0.01:0.01:10, ck.MCenterPeakSensitivityFunction( ck.MCenterPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{2}, 'lineWidth', lineWidth, 'displayName', 'Medians noIntercept' );

	subplot(2,2,4); hold on;
	plot( ck.MSurroundData('sensitivity').radiusDegs, ck.MSurroundData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesM{4} = plot( 0.01:0.01:10, ck.MSurroundPeakSensitivityFunction( ck.MSurroundPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{2}, 'lineWidth', lineWidth, 'displayName', 'Medians noIntercept' );


	%%%% Fitting with median values from the paper with intercept
	ck = CronerKaplanRGCModel( 'dataSetToFit', 'medians', 'fitIntercept', true );
	%% P Cell Fittins
	figure(hFigP);
	subplot(2,2,1); hold on;
	plot( ck.PCenterData('size').eccDegs, ck.PCenterData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesP{1}(end+1) = plot( 0.01:0.01:40, ck.PCenterRadiusFunction( ck.PCenterRadiusParams, 0.01:0.01:40 ), '--', 'color', colors{3}, 'lineWidth', lineWidth, 'displayName', 'Medians w/Intercept' );

	subplot(2,2,2); hold on;
	plot( ck.PSurroundData('size').eccDegs, ck.PSurroundData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesP{2}(end+1) = plot( 0.1:0.01:60, ck.PSurroundRadiusFunction( ck.PSurroundRadiusParams, 0.1:0.01:60 ), '--', 'color', colors{3}, 'lineWidth', lineWidth, 'displayName', 'Medians w/Intercept' );

	subplot(2,2,3); hold on;
	plot( ck.PCenterData('sensitivity').radiusDegs, ck.PCenterData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesP{3}(end+1) = plot( 0.01:0.01:10, ck.PCenterPeakSensitivityFunction( ck.PCenterPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{3}, 'lineWidth', lineWidth, 'displayName', 'Medians w/Intercept' );

	subplot(2,2,4); hold on;
	plot( ck.PSurroundData('sensitivity').radiusDegs, ck.PSurroundData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesP{4}(end+1) = plot( 0.01:0.01:10, ck.PSurroundPeakSensitivityFunction( ck.PSurroundPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{3}, 'lineWidth', lineWidth, 'displayName', 'Medians w/Intercept' );

	%% M Cell Fittins
	figure(hFigM);
	subplot(2,2,1); hold on;
	plot( ck.MCenterData('size').eccDegs, ck.MCenterData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesM{1}(end+1) = plot( 0.01:0.01:40, ck.MCenterRadiusFunction( ck.MCenterRadiusParams, 0.01:0.01:40 ), '--', 'color', colors{3}, 'lineWidth', lineWidth, 'displayName', 'Medians w/Intercept' );

	subplot(2,2,2); hold on;
	plot( ck.MSurroundData('size').eccDegs, ck.MSurroundData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesM{2}(end+1) = plot( 0.1:0.01:60, ck.MSurroundRadiusFunction( ck.MSurroundRadiusParams, 0.1:0.01:60 ), '--', 'color', colors{3}, 'lineWidth', lineWidth, 'displayName', 'Medians w/Intercept' );

	subplot(2,2,3); hold on;
	plot( ck.MCenterData('sensitivity').radiusDegs, ck.MCenterData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesM{3}(end+1) = plot( 0.01:0.01:10, ck.MCenterPeakSensitivityFunction( ck.MCenterPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{3}, 'lineWidth', lineWidth, 'displayName', 'Medians w/Intercept' );

	subplot(2,2,4); hold on;
	plot( ck.MSurroundData('sensitivity').radiusDegs, ck.MSurroundData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesM{4}(end+1) = plot( 0.01:0.01:10, ck.MSurroundPeakSensitivityFunction( ck.MSurroundPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{3}, 'lineWidth', lineWidth, 'displayName', 'Medians w/Intercept' );


	%%%% Fitting with raw digitized data with no intercept
	ck = CronerKaplanRGCModel( 'dataSetToFit', 'raw', 'fitIntercept', false );
	%% P Cell Fittins
	figure(hFigP);
	subplot(2,2,1); hold on;
	plot( ck.PCenterData('size').eccDegs, ck.PCenterData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesP{1}(end+1) = plot( 0.01:0.01:40, ck.PCenterRadiusFunction( ck.PCenterRadiusParams, 0.01:0.01:40 ), '--', 'color', colors{4}, 'lineWidth', lineWidth, 'displayName', 'Raw noIntercept' );

	subplot(2,2,2); hold on;
	plot( ck.PSurroundData('size').eccDegs, ck.PSurroundData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesP{2}(end+1) = plot( 0.1:0.01:60, ck.PSurroundRadiusFunction( ck.PSurroundRadiusParams, 0.1:0.01:60 ), '--', 'color', colors{4}, 'lineWidth', lineWidth, 'displayName', 'Raw noIntercept' );

	subplot(2,2,3); hold on;
	plot( ck.PCenterData('sensitivity').radiusDegs, ck.PCenterData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesP{3}(end+1) = plot( 0.01:0.01:10, ck.PCenterPeakSensitivityFunction( ck.PCenterPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{4}, 'lineWidth', lineWidth, 'displayName', 'Raw noIntercept' );

	subplot(2,2,4); hold on;
	plot( ck.PSurroundData('sensitivity').radiusDegs, ck.PSurroundData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesP{4}(end+1) = plot( 0.01:0.01:10, ck.PSurroundPeakSensitivityFunction( ck.PSurroundPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{4}, 'lineWidth', lineWidth, 'displayName', 'Raw noIntercept' );

	%% M Cell Fittins
	figure(hFigM);
	subplot(2,2,1); hold on;
	plot( ck.MCenterData('size').eccDegs, ck.MCenterData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesM{1}(end+1) = plot( 0.01:0.01:40, ck.MCenterRadiusFunction( ck.MCenterRadiusParams, 0.01:0.01:40 ), '--', 'color', colors{4}, 'lineWidth', lineWidth, 'displayName', 'Raw noIntercept' );

	subplot(2,2,2); hold on;
	plot( ck.MSurroundData('size').eccDegs, ck.MSurroundData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesM{2}(end+1) = plot( 0.1:0.01:60, ck.MSurroundRadiusFunction( ck.MSurroundRadiusParams, 0.1:0.01:60 ), '--', 'color', colors{4}, 'lineWidth', lineWidth, 'displayName', 'Raw noIntercept' );

	subplot(2,2,3); hold on;
	plot( ck.MCenterData('sensitivity').radiusDegs, ck.MCenterData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesM{3}(end+1) = plot( 0.01:0.01:10, ck.MCenterPeakSensitivityFunction( ck.MCenterPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{4}, 'lineWidth', lineWidth, 'displayName', 'Raw noIntercept' );

	subplot(2,2,4); hold on;
	plot( ck.MSurroundData('sensitivity').radiusDegs, ck.MSurroundData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesM{4}(end+1) = plot( 0.01:0.01:10, ck.MSurroundPeakSensitivityFunction( ck.MSurroundPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{4}, 'lineWidth', lineWidth, 'displayName', 'Raw noIntercept' );


	%%%% Fitting with raw digitized data with intercept
	ck = CronerKaplanRGCModel( 'dataSetToFit', 'raw', 'fitIntercept', true );
	%% P Cell Fittins
	figure(hFigP);
	subplot(2,2,1); hold on;
	plot( ck.PCenterData('size').eccDegs, ck.PCenterData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesP{1}(end+1) = plot( 0.01:0.01:40, ck.PCenterRadiusFunction( ck.PCenterRadiusParams, 0.01:0.01:40 ), '--', 'color', colors{5}, 'lineWidth', lineWidth, 'displayName', 'Raw w/Intercept' );

	subplot(2,2,2); hold on;
	plot( ck.PSurroundData('size').eccDegs, ck.PSurroundData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesP{2}(end+1) = plot( 0.1:0.01:60, ck.PSurroundRadiusFunction( ck.PSurroundRadiusParams, 0.1:0.01:60 ), '--', 'color', colors{5}, 'lineWidth', lineWidth, 'displayName', 'Raw w/Intercept' );

	subplot(2,2,3); hold on;
	plot( ck.PCenterData('sensitivity').radiusDegs, ck.PCenterData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesP{3}(end+1) = plot( 0.01:0.01:10, ck.PCenterPeakSensitivityFunction( ck.PCenterPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{5}, 'lineWidth', lineWidth, 'displayName', 'Raw w/Intercept' );

	subplot(2,2,4); hold on;
	plot( ck.PSurroundData('sensitivity').radiusDegs, ck.PSurroundData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesP{4}(end+1) = plot( 0.01:0.01:10, ck.PSurroundPeakSensitivityFunction( ck.PSurroundPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{5}, 'lineWidth', lineWidth, 'displayName', 'Raw w/Intercept' );

	%% M Cell Fittins
	figure(hFigM);
	subplot(2,2,1); hold on;
	plot( ck.MCenterData('size').eccDegs, ck.MCenterData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesM{1}(end+1) = plot( 0.01:0.01:40, ck.MCenterRadiusFunction( ck.MCenterRadiusParams, 0.01:0.01:40 ), '--', 'color', colors{5}, 'lineWidth', lineWidth, 'displayName', 'Raw w/Intercept' );

	subplot(2,2,2); hold on;
	plot( ck.MSurroundData('size').eccDegs, ck.MSurroundData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
	hLinesM{2}(end+1) = plot( 0.1:0.01:60, ck.MSurroundRadiusFunction( ck.MSurroundRadiusParams, 0.1:0.01:60 ), '--', 'color', colors{5}, 'lineWidth', lineWidth, 'displayName', 'Raw w/Intercept' );

	subplot(2,2,3); hold on;
	plot( ck.MCenterData('sensitivity').radiusDegs, ck.MCenterData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesM{3}(end+1) = plot( 0.01:0.01:10, ck.MCenterPeakSensitivityFunction( ck.MCenterPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{5}, 'lineWidth', lineWidth, 'displayName', 'Raw w/Intercept' );

	subplot(2,2,4); hold on;
	plot( ck.MSurroundData('sensitivity').radiusDegs, ck.MSurroundData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
	hLinesM{4}(end+1) = plot( 0.01:0.01:10, ck.MSurroundPeakSensitivityFunction( ck.MSurroundPeakSensitivityParams, 0.01:0.01:10 ), '--', 'color', colors{5}, 'lineWidth', lineWidth, 'displayName', 'Raw w/Intercept' );


	%%%% set axis properties
	% P Cell
	figure(hFigP);
	subplot(2,2,1);
	xlabel('Temporal equivalent eccentricity (deg)');
	ylabel('Center radius (deg)');
	title('Center (Figure 4a)');
	% legend( hLinesP{1}, 'location', 'northwest' );
	set( gca, 'xlim', [0 40], 'ylim', [0 0.3], 'lineWidth', lineWidth, 'fontSize', fontSize );

	subplot(2,2,2);
	xlabel('Temporal equivalent eccentricity (deg)');
	ylabel('Surround radius (deg)');
	title('Surround (Figure 4b)');
	% legend( hLinesP{2}, 'location', 'northwest' );
	set( gca, 'xlim', [0.1 100], 'ylim', [0.01 10], 'xscale', 'log', 'yscale', 'log', 'fontSize', fontSize, 'lineWidth', lineWidth );

	subplot(2,2,3);
	xlabel('Radius (deg)');
	ylabel('Peak sensitivity (imp/(s %contrast deg^2))');
	title('Center (Figure 5b)');
	% legend( hLinesP{3}, 'location', 'southwest' );
	set( gca, 'xlim', [0.01 1], 'ylim', [1 1000], 'xscale', 'log', 'yscale', 'log', 'fontSize', fontSize, 'lineWidth', lineWidth );

	subplot(2,2,4);
	xlabel('Radius (deg)');
	ylabel('Peak sensitivity (imp/(s %contrast deg^2))');
	title('Surround (Figure 5c)');
	legend( hLinesP{4}, 'location', 'southwest' );
	set( gca, 'xlim', [0.01 10], 'ylim', [0.001 100], 'xscale', 'log', 'yscale', 'log', 'fontSize', fontSize, 'lineWidth', lineWidth );

	% M Cell
	figure(hFigM);
	subplot(2,2,1);
	xlabel('Temporal equivalent eccentricity (deg)');
	ylabel('Center radius (deg)');
	title('Center (Figure 4a)');
	legend( hLinesM{1}, 'location', 'southeast' );
	set( gca, 'xlim', [0 40], 'ylim', [0 0.3], 'lineWidth', lineWidth, 'fontSize', fontSize );

	subplot(2,2,2);
	xlabel('Temporal equivalent eccentricity (deg)');
	ylabel('Surround radius (deg)');
	title('Surround (Figure 4b)');
	% legend( hLinesM{2}, 'location', 'southwest' );
	set( gca, 'xlim', [0.1 100], 'ylim', [0.01 10], 'xscale', 'log', 'yscale', 'log', 'fontSize', fontSize, 'lineWidth', lineWidth );

	subplot(2,2,3);
	xlabel('Radius (deg)');
	ylabel('Peak sensitivity (imp/(s %contrast deg^2))');
	title('Center (Figure 5b)');
	% legend( hLinesM{3}, 'location', 'southwest' );
	set( gca, 'xlim', [0.01 1], 'ylim', [1 1000], 'xscale', 'log', 'yscale', 'log', 'fontSize', fontSize, 'lineWidth', lineWidth );

	subplot(2,2,4);
	xlabel('Radius (deg)');
	ylabel('Peak sensitivity (imp/(s %contrast deg^2))');
	title('Surround (Figure 5c)');
	% legend( hLinesM{4}, 'location', 'southwest' );
	set( gca, 'xlim', [0.01 10], 'ylim', [0.001 100], 'xscale', 'log', 'yscale', 'log', 'fontSize', fontSize, 'lineWidth', lineWidth );

end