function synthesizedParams = SynthesizeRFParams(obj, eccDegs, cellType)
    % cellType:     P (midget) or M (parasol) cell

    cellType = upper(cellType);

    rng('shuffle');
    
    % Noise-free parameters
    centerRadii = obj.([cellType 'CenterRadiusFunction'])(obj.([cellType 'CenterRadiusParams']), eccDegs);
    surroundRadii = obj.([cellType 'SurroundRadiusFunction'])(obj.([cellType 'SurroundRadiusParams']), eccDegs);
    centerRadiusFitStd = std( obj.([cellType 'CenterRadiusFunction'])(obj.([cellType 'CenterRadiusParams']), obj.([cellType 'CenterData'])('size').eccDegs) - obj.([cellType 'CenterData'])('size').radiusDegs );
    surroundRadiusFitStd = std( obj.([cellType 'SurroundRadiusFunction'])(obj.([cellType 'SurroundRadiusParams']), obj.([cellType 'SurroundData'])('size').eccDegs) - obj.([cellType 'SurroundData'])('size').radiusDegs );
    centerSurroundRadiusRatios = centerRadii ./ surroundRadii;
    
    centerPeakSensitivities = obj.([cellType 'CenterPeakSensitivityFunction'])(obj.([cellType 'CenterPeakSensitivityParams']), centerRadii);
    surroundPeakSensitivities = obj.([cellType 'SurroundPeakSensitivityFunction'])(obj.([cellType 'SurroundPeakSensitivityParams']), surroundRadii);
    centerPeakSensitivityFitStd = std( log10( obj.([cellType 'CenterPeakSensitivityFunction'])(obj.([cellType 'CenterPeakSensitivityParams']), obj.([cellType 'CenterData'])('sensitivity').radiusDegs) ) - log10( obj.([cellType 'CenterData'])('sensitivity').peakSensitivity ) );
    surroundPeakSensitivityFitStd = std( log10( obj.([cellType 'SurroundPeakSensitivityFunction'])(obj.([cellType 'SurroundPeakSensitivityParams']), obj.([cellType 'SurroundData'])('sensitivity').radiusDegs) ) - log10( obj.([cellType 'SurroundData'])('sensitivity').peakSensitivity ) );
    surroundCenterPeakSensitivityRatios = surroundPeakSensitivities./centerPeakSensitivities;

    % center/surround radius ratio according to Fig. 4C;
    % surround/center peak sensitivity ratio according to Fig. 6
    if( cellType == 'P' )
        maxCenterSurroundRadiusRatio = 0.38;
        maxSurroundCenterPeakSensitivityRatio = 0.12;
    else
        maxCenterSurroundRadiusRatio = 0.75;
        maxSurroundCenterPeakSensitivityRatio = 0.44;
    end
    
    maxAttemptsNo = 1000;
    % pW = 0.15;                      %%%%%% YB. Why 0.15? %%%%%%
    pW = 0.164;                   %%%%%% YB. 3*std of regression error in Figure 11 %%%%%%
    
    % Synthesize params for each rfUnit
    parfor rfUnit = 1:numel(eccDegs)
        % From Cronner and Kaplan '94 Figure 11
        meanIntegratedSurroundToCenterSensitivityRatio = 0.466 + eccDegs(rfUnit)*0.007;%*0.5;         %%%%%% YB. Correced according to the paper %%%%%%
        integratedSurroundToCenterSensitivityRatio = Inf;
        attemptsNo = 0;
        p = pW;%abs(randn*pW);
        centerRadii(rfUnit) = 0;
        centerPeakSensitivities(rfUnit) = 0;
        surroundRadii(rfUnit) = 0;
        surroundPeakSensitivities(rfUnit) = 0;
        
        % Adjust surround radii/sensitivity so as to bring the integrated
        % surround/center sensitivity within the distribution of Figure 11
        while (attemptsNo < maxAttemptsNo) && (...
              ((abs(integratedSurroundToCenterSensitivityRatio-meanIntegratedSurroundToCenterSensitivityRatio) > p) || ...
              (integratedSurroundToCenterSensitivityRatio<0.1) || ...
              (integratedSurroundToCenterSensitivityRatio>0.9)) || ...
              (centerRadii(rfUnit) <= 0) || (centerPeakSensitivities(rfUnit) <= 0) || ...
              (surroundRadii(rfUnit) <= 0) || (surroundPeakSensitivities(rfUnit) <= 0) || ...
              centerRadii(rfUnit)/surroundRadii(rfUnit) > maxCenterSurroundRadiusRatio || ...
              surroundPeakSensitivities(rfUnit)/centerPeakSensitivities(rfUnit) > maxSurroundCenterPeakSensitivityRatio )
              
            
            [centerRadii(rfUnit), centerPeakSensitivities(rfUnit), surroundRadii(rfUnit), surroundPeakSensitivities(rfUnit)] = drawSurroundAgain(obj, cellType, centerRadiusFitStd, surroundRadiusFitStd, centerPeakSensitivityFitStd, surroundPeakSensitivityFitStd, eccDegs(rfUnit), centerSurroundRadiusRatios(rfUnit), surroundCenterPeakSensitivityRatios(rfUnit));
            
            integratedSurroundToCenterSensitivityRatio = (surroundRadii(rfUnit)/centerRadii(rfUnit))^2 * (surroundPeakSensitivities(rfUnit)/centerPeakSensitivities(rfUnit));
            
            p = abs(randn*pW);
            attemptsNo = attemptsNo + 1;
            
        end
        
        if (attemptsNo == maxAttemptsNo)
            fprintf('Failed to meet integrated sensitivity ratio after %d attempts for eccentricity %f\n', maxAttemptsNo, eccDegs(rfUnit));
        end
    end % parfor
    
    assert(all(centerPeakSensitivities>0), 'Found center peak sensitivities < 0');
    assert(all(surroundPeakSensitivities>0), 'Found surround peak sensitivities < 0');
    assert(all(surroundRadii>0), 'Found surround radii < 0');
    assert(all(centerRadii>0), 'Found center radii < 0');
    
    synthesizedParams = struct(...
        'eccDegs', eccDegs, ...
        'centerRadii', centerRadii, ...
        'surroundRadii', surroundRadii, ...
        'centerPeakSensitivities', centerPeakSensitivities, ...
        'surroundPeakSensitivities', surroundPeakSensitivities ...
        );

end


function [centerRadius,  centerPeakSensitivity, surroundRadius, surroundPeakSensitivity] = drawSurroundAgain(obj, cellType, centerRadiusFitStd, surroundRadiusFitStd, centerPeakSensitivityFitStd, surroundPeakSensitivityFitStd, ecc, centerSurroundRadiusRatio, surroundCenterPeakSensitivityRatio)
    
     if (obj.synthesisOptions.randomizeSurroundRadii)
        % Derive from model of surround size
        surroundRadius = normrnd( obj.([cellType 'SurroundRadiusFunction'])(obj.([cellType 'SurroundRadiusParams']), ecc), surroundRadiusFitStd );
        % surroundRadiusParamsNoisy = normrnd(obj.([cellType 'SurroundRadiusParams']), obj.([cellType 'SurroundRadiusParamsSE']));
        % surroundRadius = obj.([cellType 'SurroundRadiusFunction'])(surroundRadiusParamsNoisy, ecc);

        % Derive from center radius and stats of surround/center radius ratio
        %stochasticRatio = normrnd(obj.([cellType 'CenterData'])('size').radiusRatioToSurroundStats(1), obj.([cellType 'CenterData'])('size').radiusRatioToSurroundStats(2));
        %surroundRadius = centerRadius / stochasticRatio;
     else
        surroundRadius = obj.([cellType 'SurroundRadiusFunction'])(obj.([cellType 'SurroundRadiusParams']), ecc);         %%%%%% YB. This is overwritten by line 92 %%%%%%
     end
    
     if (obj.synthesisOptions.randomizeCenterRadii)
        centerRadius = normrnd( obj.([cellType 'CenterRadiusFunction'])(obj.([cellType 'CenterRadiusParams']), ecc), centerRadiusFitStd );
        % centerRadiusParamsNoisy = normrnd(obj.([cellType 'CenterRadiusParams']), obj.([cellType 'CenterRadiusParamsSE']));
        % centerRadius = obj.([cellType 'CenterRadiusFunction'])(centerRadiusParamsNoisy, ecc);
     else
        %centerRadius = surroundRadius * centerSurroundRadiusRatio;
        centerRadius = obj.([cellType 'CenterRadiusFunction'])(obj.([cellType 'CenterRadiusParams']), ecc);
     end

     if (~obj.synthesisOptions.randomizeSurroundRadii) && (obj.synthesisOptions.randomizeCenterRadii)
          surroundRadius = centerRadius / centerSurroundRadiusRatio;
     end
    
    if (obj.synthesisOptions.randomizeCenterSensitivities)
        % Noisy center peak sensitivities
        % centerPeakSensitivity = normrnd( obj.([cellType 'CenterPeakSensitivityFunction'])(obj.([cellType 'CenterPeakSensitivityParams']), centerRadius), centerPeakSensitivityFitStd );
        centerPeakSensitivity = obj.([cellType 'CenterPeakSensitivityFunction'])(obj.([cellType 'CenterPeakSensitivityParams']), centerRadius) * 10^normrnd( 0, centerPeakSensitivityFitStd );
    else
        centerPeakSensitivity = obj.([cellType 'CenterPeakSensitivityFunction'])(obj.([cellType 'CenterPeakSensitivityParams']), centerRadius);
    end
            
    

    if (obj.synthesisOptions.randomizeSurroundSensitivities)
        % surroundPeakSensitivity = normrnd( obj.([cellType 'CenterPeakSensitivityFunction'])(obj.([cellType 'SurroundPeakSensitivityParams']), surroundRadius), surroundPeakSensitivityFitStd );
        surroundPeakSensitivity = obj.([cellType 'CenterPeakSensitivityFunction'])(obj.([cellType 'SurroundPeakSensitivityParams']), surroundRadius) * 10^normrnd( 0, surroundPeakSensitivityFitStd );
    else
        % surround sensitivities from noisy center sensitivities based on noise-free ratios
        surroundPeakSensitivity = centerPeakSensitivity * surroundCenterPeakSensitivityRatio;
    end
    
end