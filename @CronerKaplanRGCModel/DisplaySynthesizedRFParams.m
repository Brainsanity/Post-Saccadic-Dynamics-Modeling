function DisplaySynthesizedRFParams(obj, PCells, MCells)

    %%%% Fitting with raw digitized data with no intercept
    ecc = reshape( repmat(0.01:0.1:40, 1, 100), 1, [] );
    if( nargin() < 2 )
        PCells = obj.SynthesizeRFParams(ecc, 'POn');
    end
    if( nargin() <3 )
        MCells = obj.SynthesizeRFParams(ecc, 'MOn');
    end

    hFigP = figure( 'NumberTitle', 'off', 'name', 'Synthesized P Cells', 'color', 'w' );
    hFigM = figure( 'NumberTitle', 'off', 'name', 'Synthesized M Cells', 'color', 'w' );
    lineWidth = 2;
    markerSize = 8;
    fontSize = 20;

    %% P Cell Fittings
    figure(hFigP);
    subplot(2,2,1); hold on;
    h(2) = plot( [PCells.eccDegs], [PCells.centerRadii], 'r.', 'markerSize', 10, 'displayName', 'Synthesized' );
    h(1) = plot( obj.PCenterData('size').eccDegs, obj.PCenterData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize, 'displayName', 'Digitized' );
    plot( 0.01:0.01:40, obj.PCenterRadiusFunction( obj.PCenterRadiusParams, 0.01:0.01:40 ), 'g--', 'lineWidth', lineWidth );
    legend( h, 'location', 'northwest' );

    subplot(2,2,2); hold on;
    plot( [PCells.eccDegs], [PCells.surroundRadii], 'r.', 'markerSize', 10 );
    plot( obj.PSurroundData('size').eccDegs, obj.PSurroundData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
    plot( 0.1:0.01:60, obj.PSurroundRadiusFunction( obj.PSurroundRadiusParams, 0.1:0.01:60 ), 'g--', 'lineWidth', lineWidth );

    subplot(2,2,3); hold on;
    plot( [PCells.centerRadii], [PCells.centerPeakSensitivities], 'r.', 'markerSize', 10 );
    plot( obj.PCenterData('sensitivity').radiusDegs, obj.PCenterData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
    plot( 0.01:0.01:10, obj.PCenterPeakSensitivityFunction( obj.PCenterPeakSensitivityParams, 0.01:0.01:10 ), 'g--', 'lineWidth', lineWidth );

    subplot(2,2,4); hold on;
    plot( [PCells.surroundRadii], [PCells.surroundPeakSensitivities], 'r.', 'markerSize', 10 );
    plot( obj.PSurroundData('sensitivity').radiusDegs, obj.PSurroundData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
    plot( 0.01:0.01:10, obj.PSurroundPeakSensitivityFunction( obj.PSurroundPeakSensitivityParams, 0.01:0.01:10 ), 'g--', 'lineWidth', lineWidth );

    %% M Cell Fittings
    figure(hFigM);
    subplot(2,2,1); hold on;
    h(2) = plot( [MCells.eccDegs], [MCells.centerRadii], 'r.', 'markerSize', 10, 'displayName', 'Fitted' );
    h(1) = plot( obj.MCenterData('size').eccDegs, obj.MCenterData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize, 'displayName', 'Digitized' );
    plot( 0.01:0.01:40, obj.MCenterRadiusFunction( obj.MCenterRadiusParams, 0.01:0.01:40 ), 'g--', 'lineWidth', lineWidth );
    legend( h, 'location', 'northwest' );

    subplot(2,2,2); hold on;
    plot( [MCells.eccDegs], [MCells.surroundRadii], 'r.', 'markerSize', 10 );
    plot( obj.MSurroundData('size').eccDegs, obj.MSurroundData('size').radiusDegs, 'ko', 'lineWidth', lineWidth, 'markerSize', markerSize );
    plot( 0.1:0.01:60, obj.MSurroundRadiusFunction( obj.MSurroundRadiusParams, 0.1:0.01:60 ), 'g--', 'lineWidth', lineWidth );

    subplot(2,2,3); hold on;
    plot( [MCells.centerRadii], [MCells.centerPeakSensitivities], 'r.', 'markerSize', 10 );
    plot( obj.MCenterData('sensitivity').radiusDegs, obj.MCenterData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
    plot( 0.01:0.01:10, obj.MCenterPeakSensitivityFunction( obj.MCenterPeakSensitivityParams, 0.01:0.01:10 ), 'g--', 'lineWidth', lineWidth );

    subplot(2,2,4); hold on;
    plot( [MCells.surroundRadii], [MCells.surroundPeakSensitivities], 'r.', 'markerSize', 10 );
    plot( obj.MSurroundData('sensitivity').radiusDegs, obj.MSurroundData('sensitivity').peakSensitivity, 'ko', 'markerSize', markerSize, 'lineWidth', lineWidth );
    plot( 0.01:0.01:10, obj.MSurroundPeakSensitivityFunction( obj.MSurroundPeakSensitivityParams, 0.01:0.01:10 ), 'g--', 'lineWidth', lineWidth );


    %% set axis properties
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
    set( gca, 'xlim', [0 40], 'ylim', [0 3], 'xscale', 'linear', 'yscale', 'linear', 'fontSize', fontSize, 'lineWidth', lineWidth );

    subplot(2,2,3);
    xlabel('Radius (deg)');
    ylabel('Peak sensitivity (imp/(s %contrast deg^2))');
    title('Center (Figure 5b)');
    % legend( hLinesP{3}, 'location', 'southwest' );
    set( gca, 'xlim', [2e-4 3e-1], 'ylim', [1e0 1.5e6], 'xscale', 'log', 'yscale', 'log', 'fontSize', fontSize, 'lineWidth', lineWidth );

    subplot(2,2,4);
    xlabel('Radius (deg)');
    ylabel('Peak sensitivity (imp/(s %contrast deg^2))');
    title('Surround (Figure 5c)');
%     legend( hLinesP{4}, 'location', 'southwest' );
    set( gca, 'xlim', [6e-5 3e0], 'ylim', [1e-3 2e8], 'xscale', 'log', 'yscale', 'log', 'fontSize', fontSize, 'lineWidth', lineWidth );

    % M Cell
    figure(hFigM);
    subplot(2,2,1);
    xlabel('Temporal equivalent eccentricity (deg)');
    ylabel('Center radius (deg)');
    title('Center (Figure 4a)');
%     legend( hLinesM{1}, 'location', 'southeast' );
    set( gca, 'xlim', [0 40], 'ylim', [0 0.4], 'lineWidth', lineWidth, 'fontSize', fontSize );

    subplot(2,2,2);
    xlabel('Temporal equivalent eccentricity (deg)');
    ylabel('Surround radius (deg)');
    title('Surround (Figure 4b)');
    % legend( hLinesM{2}, 'location', 'southwest' );
    set( gca, 'xlim', [0 40], 'ylim', [0 5], 'xscale', 'linear', 'yscale', 'linear', 'fontSize', fontSize, 'lineWidth', lineWidth );

    subplot(2,2,3);
    xlabel('Radius (deg)');
    ylabel('Peak sensitivity (imp/(s %contrast deg^2))');
    title('Center (Figure 5b)');
    % legend( hLinesM{3}, 'location', 'southwest' );
    set( gca, 'xlim', [2e-3 4e-1], 'ylim', [4e0 2e4], 'xscale', 'log', 'yscale', 'log', 'fontSize', fontSize, 'lineWidth', lineWidth );

    subplot(2,2,4);
    xlabel('Radius (deg)');
    ylabel('Peak sensitivity (imp/(s %contrast deg^2))');
    title('Surround (Figure 5c)');
    % legend( hLinesM{4}, 'location', 'southwest' );
    set( gca, 'xlim', [1e-2 7e0], 'ylim', [9e-2 3e2], 'xscale', 'log', 'yscale', 'log', 'fontSize', fontSize, 'lineWidth', lineWidth );
end
