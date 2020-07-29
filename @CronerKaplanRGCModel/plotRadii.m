function plotRadii(theAxes, d, model, pointSize, color, displayYLabel, theLabel)

    d = d('size');
    
    % Face color
    c = color + 0.5*[1 1 1];
    c(c>1)=1;
    
    hold(theAxes, 'on');
    set(theAxes, 'XScale', 'log', 'YScale', 'log');
    
    % Plot the digitized data
    if (strcmp(theLabel, 'center')) || (strcmp(theLabel, 'surround'))
        scatter(theAxes, d.eccDegs, d.radiusDegs, pointSize, 'MarkerFaceColor', c, 'MarkerEdgeColor', color);
    elseif (strcmp(theLabel, 'retinal center'))
        scatter(theAxes, d.retinalEccDegs, d.retinalRadiusDegs, pointSize, 'MarkerFaceColor', c, 'MarkerEdgeColor', color);
    else
        error('Label must b either ''center'', ''surround'', or ''retinal center''.');
    end
    
    if (isfield(d, 'eccDegsTable')) &&  (~(strcmp(theLabel, 'retinal center')))
        scatter(theAxes, d.eccDegsTable, d.radiusDegsMedianTable, 's', 'MarkerFaceColor', [0.2 0.2 0.2], 'MarkerEdgeColor', [0.2 0.2 0.2]);
    end
    
    % Plot the fitted model
    plot(theAxes, model.eccDegs, model.function(model.params, model.eccDegs), 'k-');
    
    % Add text with model equation
    t=['$\displaystyle R = ', sprintf('%2.4f',model.params(1)),' \times {E}^{', sprintf('%2.4f',model.params(2)), '}$'];
    x = get( theAxes, 'XLim' );
    y = get( theAxes, 'YLim' );
    theTextHandle = text(theAxes, 10^( log10(x) * [0.975;0.025] ), 10^( log10(y) * [0.025;0.975] ), t, 'Interpreter', 'latex');
    set(theTextHandle,'FontSize', 20, 'Color', [0.3 0.3 0.3], 'BackgroundColor', 'none');

    % Finish plot
    % set(theAxes, 'XLim', [0.1 100], 'YLim', [0.01 10]);
    % set(theAxes, 'XTick', [0.01 0.03 0.1 0.3 1 3 10 30 100], 'YTick', [0.01 0.03 0.1 0.3 1 3 10]);
    % set(theAxes, 'XTickLabel', {'0.01' '0.03' '0.1' '0.3' '1' '3' '10' '30' '100'}, ...
    %               'YTickLabel', {'0.01' '0.03' '0.1' '0.3' '1' '3' '10'});
    xlabel(theAxes, 'eccentricity (degs)');
    if (displayYLabel)
        if (strcmp(theLabel, 'retinal center'))
            ylabel(theAxes, 'retinal radius (degs)');
        else
            ylabel(theAxes, 'radius (degs)');
        end
    end
    if (~isempty(theLabel))
        if (strcmp(theLabel, 'retinal center'))
            title(theAxes, 'center', 'Color', color, 'FontSize', 30);
        else
            title(theAxes, theLabel, 'Color', color, 'FontSize', 30);
        end
    end
    grid(theAxes, 'on');
    box(theAxes, 'off');

end
