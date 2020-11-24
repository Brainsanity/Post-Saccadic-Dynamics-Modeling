function rfParams = Spacing2RFParams(obj, cellType, spacing, temporalEccDegs)
    % cellType:     POn/POff (midget) or MOn/MOff (parasol) cell

    OnOff = cellType(2:end);
    cellType = upper(cellType(1));

    if(exist('temporalEccDegs', 'var'))
        rfParams(numel(temporalEccDegs)).temporalEccDegs = [];
        temporalEccDegs = num2cell(temporalEccDegs);
        [rfParams.temporalEccDegs] = temporalEccDegs{:};
    end

    centerRadii = obj.([cellType 'CenterSpacing2Radius']).fun(spacing);
    surroundRadii = obj.([cellType 'SurroundSpacing2Radius']).fun(spacing);
    centerPeakSensitivities = obj.([cellType 'CenterPeakSensitivityFunction'])(obj.([cellType 'CenterPeakSensitivityParams']), centerRadii);
    surroundPeakSensitivities = obj.([cellType 'SurroundPeakSensitivityFunction'])(obj.([cellType 'SurroundPeakSensitivityParams']), surroundRadii);

    assert(all(centerPeakSensitivities>0), 'Found center peak sensitivities < 0');
    assert(all(surroundPeakSensitivities>0), 'Found surround peak sensitivities < 0');
    assert(all(surroundRadii>0), 'Found surround radii < 0');
    assert(all(centerRadii>0), 'Found center radii < 0');

    if( strcmpi( OnOff, 'off' ) )
        centerPeakSensitivities = -centerPeakSensitivities;
        surroundPeakSensitivities = -surroundPeakSensitivities;
    end
    
    centerRadii = num2cell(centerRadii);
    surroundRadii = num2cell(surroundRadii);
    centerPeakSensitivities = num2cell(centerPeakSensitivities);
    surroundPeakSensitivities = num2cell(surroundPeakSensitivities);

    [rfParams.centerRadii] = centerRadii{:};
    [rfParams.surroundRadii] = surroundRadii{:};
    [rfParams.centerPeakSensitivities] = centerPeakSensitivities{:};
    [rfParams.surroundPeakSensitivities] = surroundPeakSensitivities{:};
end