function synthesizedData = synthesizeData(obj, eccDegs, cellType)
    % cellType:     P (midget) or M (parasol) cell

    cellType = upper(cellType);

    rng('shuffle');
    
    % Noise-free parameters
    centerRadii = obj.([cellType 'CenterRadiusFunction'])(obj.([cellType 'CenterRadiusParams']), eccDegs);
    surroundRadii = obj.([cellType 'SurroundRadiusFunction'])(obj.([cellType 'SurroundRadiusParams']), eccDegs);
    centerSurroundRadiusRatios = centerRadii ./ surroundRadii;
    
    centerPeakSensitivities = obj.([cellType 'CenterPeakSensitivityFunction'])(obj.([cellType 'CenterPeakSensitivityParams']), centerRadii);
    surroundPeakSensitivities = obj.([cellType 'SurroundPeakSensitivityFunction'])(obj.([cellType 'SurroundPeakSensitivityParams']), surroundRadii);
    surroundCenterPeakSensitivityRatios = surroundPeakSensitivities./centerPeakSensitivities;
    
    maxAttemptsNo = 1000;
    % pW = 0.15;                      %%%%%% YB. Why 0.15? %%%%%%
    pW = 0.164;                   %%%%%% YB. 3*std of regression error in Figure 11 %%%%%%
    
    % Synthesize params for each rfUnit
    parfor rfUnit = 1:numel(eccDegs)
        % From Cronner and Kaplan '94 Figure 11
        meanIntegratedSurroundToCenterSensitivityRatio = 0.466 + eccDegs(rfUnit)*0.007;%*0.5;         %%%%%% YB. Correced according to the paper %%%%%%
        integratedSurroundToCenterSensitivityRatio = Inf;
        attemptsNo = 0;
        p = abs(randn*pW);
        centerRadii(rfUnit) = 0;
        centerPeakSensitivities(rfUnit) = 0;
        surroundRadii(rfUnit) = 0;
        surroundPeakSensitivities(rfUnit) = 0;
        
        % Adjust surround radii/sensitivity so as to bring the integrated
        % surround/center sensitivity within the distribution of Figure 11
        while (attemptsNo < maxAttemptsNo) && (...
              ((abs(integratedSurroundToCenterSensitivityRatio-meanIntegratedSurroundToCenterSensitivityRatio) > p) || ...
              (integratedSurroundToCenterSensitivityRatio<0.1) || ...
              (integratedSurroundToCenterSensitivityRatio>0.9) || ...
              (centerRadii(rfUnit) < 0) || (centerPeakSensitivities(rfUnit) < 0) || ...
              (surroundRadii(rfUnit)< 0) || (surroundPeakSensitivities(rfUnit) < 0)))
              
            
            [centerRadii(rfUnit), centerPeakSensitivities(rfUnit), ...
                surroundRadii(rfUnit), surroundPeakSensitivities(rfUnit)] = drawSurroundAgain(obj, ...
                        eccDegs(rfUnit), centerSurroundRadiusRatios(rfUnit), surroundCenterPeakSensitivityRatios(rfUnit), cellType );
            
            integratedSurroundToCenterSensitivityRatio = ...
                (surroundRadii(rfUnit)/centerRadii(rfUnit))^2 * (surroundPeakSensitivities(rfUnit)/centerPeakSensitivities(rfUnit));
            
            p = abs(randn*pW);
            attemptsNo = attemptsNo + 1;
            
        end
        
        if (attemptsNo == maxAttemptsNo)
            fprintf('Failed to meet integrated sensitivity ratio after %d attempts for eccentricity %f\n', maxAttemptsNo, eccDegs(rfUnit));
        end
    end % parfor
    
    assert(all(centerPeakSensitivities>=0), 'Found center peak sensitivities < 0');
    assert(all(surroundPeakSensitivities>=0), 'Found surround peak sensitivities < 0');
    assert(all(surroundRadii>=0), 'Found surround radii < 0');
    assert(all(centerRadii>=0), 'Found center radii < 0');
    
    synthesizedData = struct(...
        'eccDegs', eccDegs, ...
        'centerRadii', centerRadii, ...
        'surroundRadii', surroundRadii, ...
        'centerPeakSensitivities', centerPeakSensitivities, ...
        'surroundPeakSensitivities', surroundPeakSensitivities ...
        );
end

     
function [centerRadius,  centerPeakSensitivity, surroundRadius, surroundPeakSensitivity] = ...
        drawSurroundAgain(obj, eccDegs, centerSurroundRadiusRatio, surroundCenterPeakSensitivityRatio, cellType)
    
     if (obj.synthesisOptions.randomizeSurroundRadii)
        % Derive from model of surround size
        surroundRadiusParamsNoisy = normrnd(obj.([cellType 'SurroundRadiusParams']), obj.([cellType 'SurroundRadiusParamsSE']));
        surroundRadius = obj.([cellType 'SurroundRadiusFunction'])(surroundRadiusParamsNoisy, eccDegs);
        % Derive from center radius and stats of surround/center radius ratio
        %stochasticRatio = normrnd(obj.([cellType 'CenterData'])('size').radiusRatioToSurroundStats(1), obj.([cellType 'CenterData'])('size').radiusRatioToSurroundStats(2));
        %surroundRadius = centerRadius / stochasticRatio;
     else
        surroundRadius = obj.([cellType 'SurroundRadiusFunction'])(obj.([cellType 'SurroundRadiusParams']), eccDegs);         %%%%%% YB. This is overwritten by line 92 %%%%%%
     end
    
     if (obj.synthesisOptions.randomizeCenterRadii)
        centerRadiusParamsNoisy = normrnd(obj.([cellType 'CenterRadiusParams']), obj.([cellType 'CenterRadiusParamsSE']));
        centerRadius = obj.([cellType 'CenterRadiusFunction'])(centerRadiusParamsNoisy, eccDegs);
     else
        %centerRadius = surroundRadius * centerSurroundRadiusRatio;
        centerRadius = obj.([cellType 'CenterRadiusFunction'])(obj.([cellType 'CenterRadiusParams']), eccDegs);
     end

     if (~obj.synthesisOptions.randomizeSurroundRadii) && (obj.synthesisOptions.randomizeCenterRadii)
          surroundRadius = centerRadius / centerSurroundRadiusRatio;
     end
    
    if (obj.synthesisOptions.randomizeCenterSensitivities)
        % Noisy center peak sensitivities
        centerPeakSensitivityParamsNoisy = normrnd(obj.([cellType 'CenterPeakSensitivityParams']), obj.([cellType 'CenterPeakSensitivityParamsSE']));
        centerPeakSensitivity = obj.([cellType 'CenterPeakSensitivityFunction'])(centerPeakSensitivityParamsNoisy, centerRadius);
    else
        centerPeakSensitivity = obj.([cellType 'CenterPeakSensitivityFunction'])(obj.([cellType 'CenterPeakSensitivityParams']), centerRadius);
    end
            
    

    if (obj.synthesisOptions.randomizeSurroundSensitivities)
        surroundPeakSensitivityParamsNoisy = normrnd(obj.([cellType 'SurroundPeakSensitivityParams']), obj.([cellType 'SurroundPeakSensitivityParamsSE']));
        surroundPeakSensitivity = obj.([cellType 'CenterPeakSensitivityFunction'])(surroundPeakSensitivityParamsNoisy, surroundRadius);
    else
        % surround sensitivities from noisy center sensitivities based on noise-free ratios
        surroundPeakSensitivity = centerPeakSensitivity * surroundCenterPeakSensitivityRatio;
    end
    
end
