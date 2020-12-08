function DisplaySpacing2RadiusFitting(obj)

    iMeridian = 1; 	  % temporal meridian
    CS = {'Center', 'Surround'};
    PM = {'P', 'M'};
    lineWidth = 2;
    fontSize = 20;
    colors = {'r', 'b'};
    figure( 'NumberTitle', 'off', 'name', 'Radius VS Spacing Predicted by WatsonModel', 'color', 'w' );
    pause(0.1);
    jf = get(handle(gcf),'javaframe');
    jf.setMaximized(1);
    pause(0.5);
    for( iCS = 1:2 )
    	for( iPM = 1:2 )
    		subplot(2,2,iPM+(iCS-1)*2); hold on;
            h = [];
    		% iL = (iPM-1)*2+iOnOff;
    		
    		[~, spacingOn] = WatsonRGCModel.RFSpacingDensityMeridian( obj.([PM{iPM} CS{iCS} 'Data'])('size').eccDegs, WatsonRGCModel.enumeratedMeridianNames{iMeridian}, [PM{iPM} 'On'] );
    		[~, spacingOff] = WatsonRGCModel.RFSpacingDensityMeridian( obj.([PM{iPM} CS{iCS} 'Data'])('size').eccDegs, WatsonRGCModel.enumeratedMeridianNames{iMeridian}, [PM{iPM} 'Off'] );
    		spacing = mean([spacingOn, spacingOff], 2);			% the CronerKaplanRGCModel does not differentiate On/Off cells, so we use the average
    		
    		h(end+1) = plot( spacing, obj.([PM{iPM} CS{iCS} 'Data'])('size').radiusDegs, 'o', 'color', colors{iCS}, 'lineWidth', lineWidth, 'displayName', 'CK Paper & WatsonModel' );
    		h(end+1) = plot( spacing, obj.([PM{iPM} CS{iCS} 'Spacing2Radius']).fun(spacing), '-', 'color', 'k', 'lineWidth', 2, 'displayName', 'Fitting Radius ~ Spacing' );
            
    		title( [PM{iPM} CS{iCS} ' @' WatsonRGCModel.enumeratedMeridianNames{iMeridian}] );
    		if(iCS == 2) xlabel('Average spacing between On/Off (deg)'); end
    		if(iPM == 1) ylabel('Radius (deg)'); end
    		if(iCS == 1 && iPM == 1) legend( h, 'location', 'northWest' ); end
    		set( gca, 'XScale', 'log', 'YScale', 'log', 'lineWidth', lineWidth, 'fontSize', fontSize );
            if(iPM == 1)
                set(gca, 'XLim', [0.01 0.21], 'XTick', [0.01 0.1], 'XTickLabel', {'0.01' '0.1'});
                if(iCS == 1)
                    set(gca, 'YLim', [0.015 0.3], 'YTick', [0.02 0.1 0.3], 'YTickLabel', {'0.02', '0.1', '0.3'});
                else
                    set(gca, 'YLim', [0.07 2], 'YTick', [0.1 1], 'YTickLabel', {'0.1', '1'});
                end
            else
                set(gca, 'XLim', [0.07 0.21], 'XTick', [0.1 0.2], 'XTickLabel', {'0.1', '0.2'});
                if(iCS == 1)
                    set(gca, 'YLim', [0.06 0.3], 'YTick', [0.1 0.3], 'YTickLabel', {'0.1', '0.3'});
                else
                    set(gca, 'YLim', [0.25 6], 'YTick', [1 5], 'YTickLabel', {'1', '5'});
                end
            end
            if(iPM == 1 || iCS == 1)
                text( max(get(gca,'XLim')), get(gca,'YLim')*[0.99; 0.01], ...
                	 { sprintf('$r=%.3f$', obj.([PM{iPM} CS{iCS} 'Spacing2Radius']).r), ...
                	   sprintf('$p=%.3f$', obj.([PM{iPM} CS{iCS} 'Spacing2Radius']).pVal)}, ...
                	 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', fontSize, 'interpreter', 'latex' );
            end
    	end
    end
end