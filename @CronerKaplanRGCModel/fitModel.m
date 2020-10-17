function fitModel(obj)
    
    rng(1);
    % rng('shuffle');
    
    switch (lower(obj.dataSetToFit))
        case 'raw'
            fitRawData(obj);
        case 'medians'
            fitMedianData(obj);
        case 'paperformulas'
            fitRawData(obj);    % for M cells
            usePaperFits(obj);
        otherwise
            error('Unknown dataSet: ''%s''.', dataSet)
    end

end


function fitRawData(obj)

    Cell = {'P', 'M'};
    minY = [0.0044, 11.6141*0.0044, eps; 0.0112, 38.8495*0.0112, eps];     % minimal radius as half the minimal cell spacing

    for( k = 1 : 2 )
        fprintf('%s Cell - Fitting Raw Data: ecc => center radius...\n', Cell{k});
        x = obj.([Cell{k} 'CenterData'])('size').eccDegs;
        y = obj.([Cell{k} 'CenterData'])('size').radiusDegs;
        [obj.([Cell{k} 'CenterRadiusFunction']), obj.([Cell{k} 'CenterRadiusParams']), obj.([Cell{k} 'CenterRadiusParamsSE']), isConverged] = nonLinearFitData(x, y, [], [], obj.([Cell{k} 'CenterData'])('size').initialParams, obj.fitIntercept, minY(k,1));

        fprintf('%s Cell - Fitting Raw Data: ecc => surround radius...\n', Cell{k});
        x = obj.([Cell{k} 'SurroundData'])('size').eccDegs;
        y = obj.([Cell{k} 'SurroundData'])('size').radiusDegs;
        [obj.([Cell{k} 'SurroundRadiusFunction']), obj.([Cell{k} 'SurroundRadiusParams']), obj.([Cell{k} 'SurroundRadiusParamsSE']), isConverged] = nonLinearFitData(x, y, [], [], obj.([Cell{k} 'SurroundData'])('size').initialParams, obj.fitIntercept, minY(k,2));
        
        fprintf('%s Cell - Fitting Raw Data: center radius => center peak sensitivity...\n', Cell{k});
        x = obj.([Cell{k} 'CenterData'])('sensitivity').radiusDegs;
        y = obj.([Cell{k} 'CenterData'])('sensitivity').peakSensitivity;
        [obj.([Cell{k} 'CenterPeakSensitivityFunction']), obj.([Cell{k} 'CenterPeakSensitivityParams']), obj.([Cell{k} 'CenterPeakSensitivityParamsSE']), isConverged] = nonLinearFitData(x, y, [], [], obj.([Cell{k} 'CenterData'])('sensitivity').initialParams, obj.fitIntercept, minY(k,3));

        fprintf('%s Cell - Fitting Raw Data: surround radius => surround peak sensitivity...\n', Cell{k});
        x = obj.([Cell{k} 'SurroundData'])('sensitivity').radiusDegs;
        y = obj.([Cell{k} 'SurroundData'])('sensitivity').peakSensitivity;
        [obj.([Cell{k} 'SurroundPeakSensitivityFunction']), obj.([Cell{k} 'SurroundPeakSensitivityParams']), obj.([Cell{k} 'SurroundPeakSensitivityParamsSE']), isConverged] = nonLinearFitData(x, y, [], [], obj.([Cell{k} 'SurroundData'])('sensitivity').initialParams, obj.fitIntercept, minY(k,3));
    end
end

function usePaperFits(obj)
    %% For P cell only

    minY = [0.0044, 11.6141*0.0044, eps];     % minimal radius as half the minimal cell spacing

    % Fit the ecc - center radius data
    fprintf('P Cell - Fitting Median Data: ecc => center radius...\n');
    x = obj.PCenterData('size').eccDegsTable;
    yMedian = obj.PCenterData('size').radiusDegsMedianTable;
    yIQR = obj.PCenterData('size').radiusDegsIQRTable;
    ySamplesNum = obj.PCenterData('size').samplesTable;
    [obj.PCenterRadiusFunction, obj.PCenterRadiusParams, obj.PCenterRadiusParamsSE, isConverged] = nonLinearFitData(x, yMedian, yIQR, ySamplesNum, obj.PCenterData('size').initialParams, obj.fitIntercept, minY(1));
    
    
    % The ecc - surround radius equation from the paper (Figure 4 caption)
    obj.PSurroundRadiusFunction = @(p,x)(p(1)*x.^p(2));
    obj.PSurroundRadiusParams = [0.203 0.472];
    obj.PSurroundRadiusParamsSE = [0 0];
    
    
    % The radius - center sensitivity equation from the paper (Figure 5 caption)
    obj.PCenterPeakSensitivityFunction = @(p,x)(p(1)*x.^p(2));
    obj.PCenterPeakSensitivityParams = [0.391 -1.850];
    obj.PCenterPeakSensitivityParamsSE = [0 0];
    
    % The radius - surround sensitivity equation from the paper (Figure 5 caption)
    obj.PSurroundPeakSensitivityFunction = @(p,x)(p(1)*x.^p(2));
    obj.PSurroundPeakSensitivityParams = [0.128 -2.147];
    obj.PSurroundPeakSensitivityParamsSE = [0 0];
end

function fitMedianData(obj)

    Cell = {'P', 'M'};
    minY = [0.0044, 11.6141*0.0044, eps; 0.0112, 38.8495*0.0112, eps];     % minimal radius as half the minimal cell spacing

    for( k = 1 : 2 )
        % Fit the ecc - center radius data
        fprintf('%s Cell - Fitting Median Data: ecc => center radius...\n', Cell{k});
        x = obj.([Cell{k} 'CenterData'])('size').eccDegsTable;
        yMedian = obj.([Cell{k} 'CenterData'])('size').radiusDegsMedianTable;
        yIQR = obj.([Cell{k} 'CenterData'])('size').radiusDegsIQRTable;
        ySamplesNum = obj.([Cell{k} 'CenterData'])('size').samplesTable;
        [obj.([Cell{k} 'CenterRadiusFunction']), obj.([Cell{k} 'CenterRadiusParams']), obj.([Cell{k} 'CenterRadiusParamsSE']), isConverged] = nonLinearFitData(x, yMedian, yIQR, ySamplesNum, obj.([Cell{k} 'CenterData'])('size').initialParams, obj.fitIntercept, minY(k,1));
        
        % Fit the ecc - surround radius data
        fprintf('%s Cell - Fitting Median Data: ecc => surround radius...\n', Cell{k});
        x = obj.([Cell{k} 'SurroundData'])('size').eccDegsTable;
        yMedian = obj.([Cell{k} 'SurroundData'])('size').radiusDegsMedianTable;
        yIQR = obj.([Cell{k} 'SurroundData'])('size').radiusDegsIQRTable;
        ySamplesNum = obj.([Cell{k} 'SurroundData'])('size').samplesTable;
        [obj.([Cell{k} 'SurroundRadiusFunction']), obj.([Cell{k} 'SurroundRadiusParams']), obj.([Cell{k} 'SurroundRadiusParamsSE']), isConverged] = nonLinearFitData(x, yMedian, yIQR, ySamplesNum, obj.([Cell{k} 'SurroundData'])('size').initialParams, obj.fitIntercept, minY(k,2));

        
        % Fit the radius - center sensitivity data
        fprintf('%s Cell - Fitting Median Data: center radius => center peak sensitivity...\n', Cell{k});
        x = obj.([Cell{k} 'CenterData'])('size').radiusDegsMedianTable;
        yMedian = obj.([Cell{k} 'CenterData'])('sensitivity').peakSensitivityMedianTable;
        yIQR = obj.([Cell{k} 'CenterData'])('sensitivity').peakSensitivityIQRTable;
        ySamplesNum = obj.([Cell{k} 'CenterData'])('sensitivity').samplesTable;
        [obj.([Cell{k} 'CenterPeakSensitivityFunction']), obj.([Cell{k} 'CenterPeakSensitivityParams']), obj.([Cell{k} 'CenterPeakSensitivityParamsSE']), isConverged] = nonLinearFitData(x, yMedian, yIQR, ySamplesNum, obj.([Cell{k} 'CenterData'])('sensitivity').initialParams, obj.fitIntercept,minY(k,3));
        
        
        % Fit the radius - surround sensitivity data
        fprintf('%s Cell - Fitting Median Data: surround radius => surround peak sensitivity...\n', Cell{k});
        x = obj.([Cell{k} 'SurroundData'])('size').radiusDegsMedianTable;
        yMedian = obj.([Cell{k} 'SurroundData'])('sensitivity').peakSensitivityMedianTable;
        yIQR = obj.([Cell{k} 'SurroundData'])('sensitivity').peakSensitivityIQRTable;
        ySamplesNum = obj.([Cell{k} 'SurroundData'])('sensitivity').samplesTable;
        [obj.([Cell{k} 'SurroundPeakSensitivityFunction']), obj.([Cell{k} 'SurroundPeakSensitivityParams']), obj.([Cell{k} 'SurroundPeakSensitivityParamsSE']), isConverged] = nonLinearFitData(x, yMedian, yIQR, ySamplesNum, obj.([Cell{k} 'SurroundData'])('sensitivity').initialParams, obj.fitIntercept, minY(k,3));
    end
end


function [powerFunction, fittedParams, fittedParamsSE, isConverged] = nonLinearFitData(x, y, yIQR, ySamplesNum, initialParams, fitIntercept, minY)
    
    % Objective Function
    if(fitIntercept)
        powerFunction = @(p,x) max(minY, (p(1)*x.^p(2)+p(3)));
    else
        powerFunction = @(p,x) max(minY, (p(1)*x.^p(2)));
    end

    opts.RobustWgtFun = 'talwar';
    opts.MaxIter = 3000;
    if (~isempty(yIQR))
        xx = [];
        yy = [];
        for k = 1:numel(y)
            ySamples = normrnd(y(k), yIQR(k)/1.35, [1 ySamplesNum(k)]);
            yy = cat(2, yy, ySamples);
            xx = cat(2, xx, repmat(x(k), [1 ySamplesNum(k)]));
        end
        y = yy;
        x = xx;
    end

    isConverged = false;
    for(itr = 1 : 100)
        lastwarn('');
        [fittedParams,~,~,varCovarianceMatrix,~] = nlinfit(x, y, powerFunction, initialParams(1:2+fitIntercept), opts);
        msg = lastwarn;
        if( any( abs( fittedParams - initialParams(1:2+fitIntercept) ) > 1e-6 | abs( fittedParams - initialParams(1:2+fitIntercept) ) ./ fittedParams > 1e-6 ) || ~isempty(msg) )
            initialParams(1:2+fitIntercept) = fittedParams;
        else
            isConverged = true;
            fprintf('\tFitting converged after %d iterations of nlinfit()\n\n', itr);
            break;
        end
    end
    if(~isConverged)
        fprintf('\tWarning: Fitting not converged after %d iterations of nlinfit()\n\n', 100);
    end

    % standard error of the mean
    fittedParamsSE = sqrt(diag(varCovarianceMatrix));
    fittedParamsSE = fittedParamsSE';
   
    % make it standard deviation
    %fittedParamsSE = fittedParamsSE * sqrt(mean(ySamplesNum));
end
