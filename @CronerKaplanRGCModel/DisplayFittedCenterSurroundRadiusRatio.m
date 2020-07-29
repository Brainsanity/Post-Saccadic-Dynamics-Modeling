function DisplayFittedCenterSurroundRadiusRatio(obj)
	% obj = CronerKaplanRGCModel( 'dataSetToFit', 'raw', 'fitIntercept', false );
	ecc = 0:0.0001:0.2;

	figure( 'NumberTitle', 'off', 'name', 'Ratio of Center-Surround Radius Fitted with Raw Data and No InterCept', 'color', 'w' );
	lineWidth = 2;
	fontSize = 20;

	subplot(1,2,1); hold on;
	[ax, h2, h3] = plotyy( ecc, obj.PSurroundRadiusFunction( obj.PSurroundRadiusParams, ecc ), ecc, obj.PSurroundRadiusFunction( obj.PSurroundRadiusParams, ecc ) ./ obj.PCenterRadiusFunction( obj.PCenterRadiusParams, ecc ) );
	set( h2, 'color', 'c', 'lineWidth', lineWidth, 'displayName', 'Surround' );
	set( h3, 'lineStyle', '--', 'color', 'r', 'lineWidth', lineWidth, 'displayName', 'Ratio' );
	h1 = plot( ax(1), ecc, obj.PCenterRadiusFunction( obj.PCenterRadiusParams, ecc ), 'b', 'lineWidth', lineWidth, 'displayName', 'Center' );
	xlabel('Eccentricity (deg)');
	ylabel('Radius (deg)');
	title('P Cell');
	legend( [h1 h2 h3], 'location', 'northwest' );
	set( gca, 'ylim', [0 0.085], 'ytick', 0:0.01:0.085, 'fontSize', fontSize, 'lineWidth', lineWidth );
	set( ax(2), 'YColor', 'r', 'fontSize', fontSize, 'lineWidth', lineWidth );

	subplot(1,2,2); hold on;
	[ax, h2, h3] = plotyy( ecc, obj.MSurroundRadiusFunction( obj.MSurroundRadiusParams, ecc ), ecc, obj.MSurroundRadiusFunction( obj.MSurroundRadiusParams, ecc ) ./ obj.MCenterRadiusFunction( obj.MCenterRadiusParams, ecc ) );
	set( h2, 'color', 'c', 'lineWidth', lineWidth, 'displayName', 'Surround' );
	set( h3, 'lineStyle', '--', 'color', 'r', 'lineWidth', lineWidth, 'displayName', 'Ratio' );
	h1 = plot( ax(1), ecc, obj.MCenterRadiusFunction( obj.MCenterRadiusParams, ecc ), 'b', 'lineWidth', lineWidth, 'displayName', 'Center' );
	xlabel('Eccentricity (deg)');
	ylabel('Radius (deg)');
	title('M Cell');
	% legend( [h1 h2 h3], 'location', 'northwest' );
	set( gca, 'ylim', [0 0.6], 'ytick', 0:0.1:0.6, 'fontSize', fontSize, 'lineWidth', lineWidth );
	set( ax(2), 'YColor', 'r', 'fontSize', fontSize, 'lineWidth', lineWidth );
end