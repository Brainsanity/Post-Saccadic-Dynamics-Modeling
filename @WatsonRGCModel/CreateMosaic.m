function [rfPositions, maxMovements, reTriangulationIterations, terminationReason, histogramDataHistory, histogramWidths] = CreateMosaic(radius, cellType, varargin)
    % Create a mosaic of the specified RGC cell type, and return the visual locations of all cells
    %   radius:         radius of the mosaic in degrees
    %   cellType:       P, POn, POff, M, MOn, MOff
    %
    %   rfPositions:    visual locations of all cells in the created mosaic
    %                   1st row for horizontal, 2nd row for vertical
    %                   negative for nasal & inferior, and positive for temporal & superior

    args = inputParser;       
    args.addParameter('saveFolder', './', @ischar);
    args.addParameter('saveName', '', @ischar);
    args.addParameter('ellipseAxes', [1 1.2247], @isnumeric);
    args.addParameter('maxIterations', 30, @isnumeric);
    args.addParameter('saveHistory', false, @islogical);    % whether save history of rfPositions
    args.parse(varargin{:});
    
    % Visualize mosaic and progress
    visualizeProgress = false;

    % Termination conditions
    % 1. Stop if cells move less than this positional tolerance (x gridParams.lambdaMin) in microns
    dTolerance = 1.0e-4;
    
    % 2. Stop if we exceed this many iterations
    maxIterations = args.Results.maxIterations;
    
    % 3. Trigger Delayun triangularization if rfmovement exceeds this number
    maxMovementPercentile = 20;    
    
    % 4. Do not trigger Delayun triangularization if less than minIterationsBeforeRetriangulation have passed since last one
    minIterationsBeforeRetriangulation = 5;
    
    % 5. Trigger Delayun triangularization if more than maxIterationsBeforeRetriangulation have passed since last one
    maxIterationsBeforeRetriangulation = 30;
    
    % 6. Interval to query user whether he/she wants to terminate
    queryUserIntervalMinutes = 60*12;
    
    % Save filename
    saveFileName = fullfile(args.Results.saveFolder, sprintf('Mosaic_%s_Radius%2.1fdeg_maxMovPrctile%d.mat', cellType, radius, maxMovementPercentile));
    
    gridParams.ellipseAxes = args.Results.ellipseAxes;
    [~,gridParams.lambdaMin] = WatsonRGCModel.RFSpacingDensityMeridian(0, WatsonRGCModel.enumeratedMeridianNames{1}, cellType);%2;     %%%%%% YB. 2020.07.10. To match maximal density of P cell. Not Perfect. %%%%%%
    gridParams.borderTolerance = 0.001 * gridParams.lambdaMin;
    gridParams.dTolerance = gridParams.lambdaMin * dTolerance;
    gridParams.maxMovementPercentile = maxMovementPercentile;

    gridParams.saveHistory = args.Results.saveHistory;
    
    tStart = tic;

    %% Generate initial RF positions and downsample according to the density
    % rfPositions = generateInitialRFpositions(radius*1.07, gridParams.lambdaMin, gridParams.micronsPerDegree);
    %%%%%% YB. 2020.07.06. %%%%%%
    scaleF = sqrt(3) / 2;
    margin = 0.1;
    eccX = 0 : scaleF * gridParams.lambdaMin : radius*1.07*1.2+margin;
    eccY = 0 : gridParams.lambdaMin : radius*1.07*1.2+margin;
    [X, Y] = meshgrid( [-eccX(end:-1:2), eccX], [-eccY(end:-1:2), eccY] );
    Y(:,2:2:end) = Y(:,2:2:end) - 0.5*gridParams.lambdaMin;
    rfPositions = ([X(:), Y(:)]);
    clear eccX eccY X Y;


    %% Down sample rfPositions
    rng('shuffle');    
    rfsNum = size(rfPositions,1);
    if (rfsNum > 1000*1000)
        fprintf('Started with %2.1f million RFs, time lapsed: %f minutes\n', rfsNum/1000000, toc(tStart)/60);
    else
        fprintf('Started with %2.1f thousand RFs, time lapsed: %f minutes\n', rfsNum/1000, toc(tStart)/60);
    end

    % Remove cells outside the desired elliptical region
    fprintf('Removing cells outside the ellipse ...');
    gridParams.radius = max(abs(rfPositions(:)));
    d = ellipticalDomainFunction(rfPositions, gridParams.radius, gridParams.ellipseAxes);
    rfPositions = rfPositions(d < gridParams.borderTolerance, :);
    fprintf('... time lapsed: %f minutes.\n', toc(tStart)/60);

    % sample probabilistically according to cellspacingFunction
    rfsNum = size(rfPositions,1);
    if (rfsNum > 1000*1000)
        fprintf('Computing separations for %2.1f million nodes ...', rfsNum/1000000);
    else
        fprintf('Computing separations for %2.1f thousand nodes ...', rfsNum/1000);
    end
    [~, spacing] = WatsonRGCModel.RFSpacingDensity(rfPositions, cellType);
    fprintf('... time lapsed: %f minutes.',  toc(tStart)/60);

    fprintf('\nProbabilistic sampling ...');
    normalizedSpacing = spacing / gridParams.lambdaMin;
    densityP = WatsonRGCModel.densityFromSpacing(normalizedSpacing);

    % Remove cells accordingly
    fixedRFPositionsRadius = 1;     % central region not to remove
    radii = sqrt(sum(rfPositions.^2,2));
    keptRFIndices = (rand(size(rfPositions, 1), 1) < densityP) | ((radii < fixedRFPositionsRadius*gridParams.lambdaMin));
    rfPositions = rfPositions(keptRFIndices, :);
    fprintf(' ... done ! After %f minutes.\n', toc(tStart)/60);


    
    %% Adjust cell positions
    rfsNum = size(rfPositions,1);
    if (rfsNum > 1000*1000)
        fprintf('Iteration: 0, Adusting %2.1f million nodes, time lapsed: %f minutes\n', size(rfPositions,1)/1000000, toc(tStart)/60);
    else
        fprintf('Iteration: 0, Adusting %2.1f thousand nodes, time lapsed: %f minutes\n', size(rfPositions,1)/1000, toc(tStart)/60);
    end
    
    % Do it
    [rfPositions, rfPositionsHistory, iterationsHistory, rfPositionsHistory2, iterationsHistory2, maxMovements, reTriangulationIterations, terminationReason, histogramDataHistory, histogramWidths] = ...
        smoothGrid(gridParams, rfPositions,  minIterationsBeforeRetriangulation, maxIterationsBeforeRetriangulation, maxIterations, queryUserIntervalMinutes, visualizeProgress, radius, cellType, tStart, saveFileName);

    % normalize according to model predicted number of cells within a central region
    r = radius*0.1;
    d = 0.001;
    [x,y] = meshgrid(-r:d:r, -r:d:r);
    density = WatsonRGCModel.RFSpacingDensity([x(:),y(:)], cellType);
    nGT = sum( density * d^2 );
    n = sum( abs(rfPositions(:,1)) <= r+d/2 & abs(rfPositions(:,2)) <= r+d/2 );
    rfPositions = rfPositions * sqrt(n/nGT);

    % Save results
    save(saveFileName, 'rfPositions', 'rfPositionsHistory', 'iterationsHistory', 'rfPositionsHistory2', 'iterationsHistory2', 'maxMovements', 'reTriangulationIterations', 'terminationReason', 'histogramDataHistory', 'histogramWidths', '-v7.3');
    fprintf('History saved in %s\n', saveFileName);

end

    
function [rfPositions, rfPositionsHistory, iterationsHistory, rfPositionsHistory2, iterationsHistory2, maxMovements, reTriangulationIterations, terminationReason, histogramDataHistory, histogramWidths] = ...
    smoothGrid(gridParams, rfPositions,  minIterationsBeforeRetriangulation, maxIterationsBeforeRetriangulation, maxIterations, queryUserIntervalMinutes, visualizeProgress, radius, cellType, tStart, saveFileName)

    gridParams.maxIterations = maxIterations;
    deps = gridParams.lambdaMin * sqrt(eps);
    deltaT = 0.2;

    % Initialize convergence
    forceMagnitudes = [];

    % Turn off Delaunay triangularization warning
    warning('off', 'MATLAB:qhullmx:InternalWarning');
    
    % Number of cells
    rfsNum = size(rfPositions, 1);
    
    % Iteratively adjust the cell positions until the forces between nodes
    % (rfPositions) reach equilibrium.
    notConverged = true;
    terminateNowDueToReductionInLatticeQuality = false;
    oldRFPositions = inf;
    
    iteration = 0;
    maxMovements = [];
    rfPositionsHistory = [];
    rfPositionsHistory2 = [];           %%%%%%  YB. Add. 2020.07.02  %%%%%%
    iterationsHistory = [];
    iterationsHistory2 = [];
    
    lastTriangularizationAtIteration = iteration;
    minimalIterationsPerformedAfterLastTriangularization = 0;
    histogramWidths = [];
    histogramDataHistory = [];
    reTriangulationIterations = [];
    timePrevious = clock;
    userRequestTerminationAtIteration = [];
    terminateNow = false;
    
    
    while all(~isnan(rfPositions(:))) && (~terminateNow) && (~terminateNowDueToReductionInLatticeQuality) && (notConverged) && (iteration <= gridParams.maxIterations) || ...
            ((lastTriangularizationAtIteration > iteration-minimalIterationsPerformedAfterLastTriangularization)&&(iteration > gridParams.maxIterations))
        
        if ((lastTriangularizationAtIteration > iteration-minimalIterationsPerformedAfterLastTriangularization)&&(iteration > gridParams.maxIterations))
            fprintf('Exceed max iterations (%d), but last triangularization was less than %d iterations before so we will do one more iteration\n', gridParams.maxIterations,minimalIterationsPerformedAfterLastTriangularization);
        end
        
        iteration = iteration + 1;

        % compute cell positional diffs
        positionalDiffs = sqrt(sum((rfPositions-oldRFPositions).^ 2,2)); 
        
        % Check if we need to re-triangulate
        %positionalDiffsMetric = max(positionalDiffs);
        %positionalDiffsMetric = median(positionalDiffs);
       % positionalDiffsMetric = prctile(positionalDiffs, 99);
        
        % We need to triangulate again if the positionalDiff is above the set tolerance
%         if ((positionalDiffsMetric > gridParams.positionalDiffToleranceForTriangularization))
%             reTriangulationIsNeeded = true;
%             triangularizationTriggerEvent = 'movement > tolerance';
%         end
        
        % We need to triangulate again if the movement in the current iteration was > the average movement in the last 2 iterations 
        if (numel(maxMovements)>3) && (maxMovements(iteration-1) > 1.02*0.5*(maxMovements(iteration-2)+maxMovements(iteration-3)))
            reTriangulationIsNeeded = true;
             triangularizationTriggerEvent = 'movement stopped decreasing';
        end
        
        % We need to triangulate again if we went for maxIterationsToRetriangulate + some more since last triangularization
        if ((abs(lastTriangularizationAtIteration-iteration-1)) > maxIterationsBeforeRetriangulation+min([10 round(iteration/20)]))
            reTriangulationIsNeeded = true;
            triangularizationTriggerEvent = 'maxIterations passed';
        end
        
        % Do not triangulare if we did one less than minIterationsBeforeRetriangulation before
        if ((abs(lastTriangularizationAtIteration-iteration)) < minIterationsBeforeRetriangulation+min([5 round(iteration/50)]))
            reTriangulationIsNeeded = false;
        end
        
        %
        if (iteration==1)
            reTriangulationIsNeeded = true;
            triangularizationTriggerEvent = '1st iteration';
        end
        
        if (reTriangulationIsNeeded)
            lastTriangularizationAtIteration = iteration;
            % save old cell positions
            oldRFPositions = rfPositions;
            
            % Perform new Delaunay triangulation to determine the updated
            % topology of the truss.
            triangleIndices = delaunayn(rfPositions);
            % Compute the centroids of all triangles
            centroidPositions = 1.0/3.0 * (...
                    rfPositions(triangleIndices(:, 1), :) + ...
                    rfPositions(triangleIndices(:, 2), :) + ...
                    rfPositions(triangleIndices(:, 3), :));
            
            % Remove centroids outside the desired region by applying the
            % signed distance function
            d = ellipticalDomainFunction(centroidPositions, gridParams.radius, gridParams.ellipseAxes);
            triangleIndices = triangleIndices(d < gridParams.borderTolerance, :);
            
           % Create a list of the unique springs (each spring connecting 2 cells)
           springs = [...
                    triangleIndices(:, [1, 2]); ...
                    triangleIndices(:, [1, 3]); ...
                    triangleIndices(:, [2, 3]) ...
           ];
           springs = unique(sort(springs, 2), 'rows');
            
           % find all springs connected to this cell
           springIndices = cell(1,rfsNum);
           for rfIndex = 1:rfsNum
               springIndices{rfIndex} = find((springs(:, 1) == rfIndex) | (springs(:, 2) == rfIndex));
           end
        end % reTriangulationIsNeeded
        
        % Compute spring vectors
        springVectors =  rfPositions(springs(:, 1), :) - rfPositions(springs(:, 2), :);
        % their centers
        springCenters = (rfPositions(springs(:, 1), :) + rfPositions(springs(:, 2), :)) / 2.0;
        % and their lengths
        springLengths = sqrt(sum(springVectors.^2, 2));

        if (reTriangulationIsNeeded)
            % Compute desired spring lengths
            [~, desiredSpringLengths] = WatsonRGCModel.RFSpacingDensity(springCenters, cellType);
        end
        
        % Normalize spring lengths
        normalizingFactor = sqrt(sum(springLengths .^ 2) / ...
            sum(desiredSpringLengths .^ 2));
        % normalizingFactor = gridParams.lambdaMin / min(desiredSpringLengths);   %%%%%% YB. 2020.07.10. %%%%%%
        desiredSpringLengths = desiredSpringLengths * normalizingFactor;
        
        gain = 1;1.1;
        springForces = max(gain * desiredSpringLengths - springLengths, 0);

        % compute x, y-components of forces on each of the springs
        springForceXYcomponents = abs(springForces ./ springLengths * [1, 1] .* springVectors);

        % Compute net forces on each cell
        netForceVectors = zeros(rfsNum, 2);
        
        parfor rfIndex = 1:rfsNum
           % compute net force from all connected springs
           deltaPos = -bsxfun(@minus, springCenters(springIndices{rfIndex}, :), rfPositions(rfIndex, :));
           netForceVectors(rfIndex, :) = sum(sign(deltaPos) .* springForceXYcomponents(springIndices{rfIndex}, :), 1);
        end
            
        % update cell positions according to netForceVectors
        oldRFPositions = rfPositions;
        rfPositions = rfPositions + deltaT * netForceVectors;

        % And project them back to the domain
        d = ellipticalDomainFunction(rfPositions, gridParams.radius, gridParams.ellipseAxes);
        outIdx = d > 0;
        if (sum(outIdx) > 1)
                % Compute numerical gradient along x-positions
                dXgradient = (ellipticalDomainFunction( [rfPositions(outIdx, 1) + deps, rfPositions(outIdx, 2)], gridParams.radius, gridParams.ellipseAxes) - d(outIdx)) / deps;
                dYgradient = (ellipticalDomainFunction( [rfPositions(outIdx, 1), rfPositions(outIdx, 2) + deps], gridParams.radius, gridParams.ellipseAxes) - d(outIdx)) / deps;

                % Project these points back to boundary
                delta = [d(outIdx) .* dXgradient, d(outIdx) .* dYgradient];
                delta( abs(delta) == Inf ) = 0;
                rfPositions(outIdx, :) = rfPositions(outIdx, :) - delta;       %%%%%% move back by d*gradient ******
        end
            
        % Check if all interior nodes move less than dTolerance
        movementAmplitudes = sqrt(sum(deltaT * netForceVectors(d < -gridParams.borderTolerance, :) .^2 , 2));
        maxMovement = prctile(movementAmplitudes, gridParams.maxMovementPercentile);
        maxMovements(iteration) = maxMovement;
        
        if maxMovement < gridParams.dTolerance
            notConverged = false; 
        end
          
        % Check for early termination due to decrease in hex lattice quality
        if (true||reTriangulationIsNeeded)
            reTriangulationIterations = cat(2,reTriangulationIterations, iteration);
            [terminateNowDueToReductionInLatticeQuality, histogramData, histogramWidths, histogramDiffWidths, checkedBins] = ...
                checkForEarlyTerminationDueToHexLatticeQualityDecrease(rfPositions, triangleIndices, histogramWidths);
            histogramDataHistory = [histogramDataHistory, histogramData];
        end
        
        if  ( reTriangulationIsNeeded || terminateNowDueToReductionInLatticeQuality)  

            % See if we need to query the user about terminating
            timeLapsedMinutes = etime(clock, timePrevious)/60;
            
            if (timeLapsedMinutes > queryUserIntervalMinutes)
                queryUserWhetherToTerminateSoon = true;
            else
                queryUserWhetherToTerminateSoon = false;
            end
            
            fprintf('\t>Triangularization at iteration: %d/%d (%s) - medianMov: %2.6f, tolerance: %2.3f, time lapsed: %2.1f minutes\n', ...
                iteration, gridParams.maxIterations, triangularizationTriggerEvent, maxMovement, gridParams.dTolerance, toc(tStart)/60);
            
            if(gridParams.saveHistory)
                if (isempty(rfPositionsHistory))
                    rfPositionsHistory(1,:,:) = single(oldRFPositions);
                    iterationsHistory = iteration;
                else
                    rfPositionsHistory = cat(1, rfPositionsHistory, reshape(single(oldRFPositions), [1 size(oldRFPositions,1) size(oldRFPositions,2)]));
                    iterationsHistory = cat(2, iterationsHistory, iteration);
                end
                
                if (visualizeProgress)
                    plotMosaic([], rfPositions, triangleIndices, maxMovements, reTriangulationIterations, histogramDiffWidths, histogramData, checkedBins, gridParams.dTolerance, radius*2, gridParams.micronsPerDegree);
                else
                    plotMovementSequence([],maxMovements, gridParams.dTolerance)
                    plotMeshQuality([],histogramData, checkedBins, iterationsHistory);
                end
            end
        else
            % Store positions at intermediate iterations (every 5 iterations)
            if(gridParams.saveHistory)
                if (mod(iteration,3)==0)
                    if (isempty(rfPositionsHistory2))       %%%%%%  YB. Change from rfPositionsHistory to rfPositionHistory2. 2020.07.02 %%%%%%
                        rfPositionsHistory2(1,:,:) = single(rfPositions);
                        iterationsHistory2 = iteration;
                    else
                        rfPositionsHistory2 = cat(1, rfPositionsHistory2, reshape(single(rfPositions), [1 size(rfPositions,1) size(rfPositions,2)]));
                        iterationsHistory2 = cat(2, iterationsHistory2, iteration);
                    end
                end
            end
            
        end
        
        
        if (queryUserWhetherToTerminateSoon)
            fprintf('Another %d minute period has passed. Terminate soon?', queryUserIntervalMinutes);
            userTermination = GetWithDefault(' If so enter # of iteration to terminate on. Otherwise hit enter to continue', 'continue');
            if (~strcmp(userTermination, 'continue'))
                userRequestTerminationAtIteration = str2double(userTermination);
                if (isnan(userRequestTerminationAtIteration))
                    userRequestTerminationAtIteration = [];
                end
            else
                fprintf('OK, will ask again in %d minutes.', queryUserIntervalMinutes);
            end
            timePrevious = clock;
        end
        queryUserWhetherToTerminateSoon = false;
        
        
        if (terminateNowDueToReductionInLatticeQuality)
            % Return the last cell positions
            rfPositions = oldRFPositions;
        end
        
        if (~isempty(userRequestTerminationAtIteration)) && (iteration >= userRequestTerminationAtIteration)
            rfPositionsHistory = cat(1, rfPositionsHistory, reshape(single(oldRFPositions), [1 size(oldRFPositions,1) size(oldRFPositions,2)]));
            iterationsHistory = cat(2, iterationsHistory, iteration);
            reTriangulationIterations = cat(2,reTriangulationIterations, iteration);
            fprintf('Current iteration: %d, user request stop iteration: %d\n', iteration,userRequestTerminationAtIteration)
            terminateNow = true;
        end

        if (any(isnan(rfPositions(:))) || any(abs(rfPositions(:)) == Inf))
            rfPositions = oldRFPositions;
            terminationReason = 'Get NaN values for rfPositions';
        end
        
        save(sprintf('%s_Iteration%d.mat', saveFileName(1:end-4), iteration), 'rfPositions', 'iterationsHistory', 'iterationsHistory2', 'maxMovements', 'reTriangulationIterations', 'histogramDataHistory', 'histogramWidths', '-v7.3');
        
    end
    
    if (terminateNow)
            terminationReason = sprintf('User requested termination at iteration %d', userRequestTerminationAtIteration);
    else
        if (notConverged)
            if (terminateNowDueToReductionInLatticeQuality)
                terminationReason = 'Decrease in hex lattice quality.';
            else
                terminationReason = 'Exceeded max number of iterations.';
            end
        else
            terminationReason = 'Converged.';
        end
    end

    fprintf('Hex lattice adjustment ended. Reason: %s\n', terminationReason);
end

function distances = ellipticalDomainFunction(rfPositions, radius, ellipseAxes)
    xx = rfPositions(:, 1);
    yy = rfPositions(:, 2);
    radii = sqrt((xx / ellipseAxes(1)) .^ 2 + (yy / ellipseAxes(2)) .^ 2);
    distances = radii - radius;
end


function [terminateNow, histogramData, widths, diffWidths, bin1Percent] = checkForEarlyTerminationDueToHexLatticeQualityDecrease(currentRFPositions, triangleIndices, widths)
    
    qDist = computeQuality(currentRFPositions, triangleIndices);
    qBins = [0.5:0.01:1.0];
    [counts,centers] = hist(qDist, qBins);
    bin1Percent = prctile(qDist,[0.8 3 7 15 99.8]);
    [~, idx1] = min(abs(centers-bin1Percent(2)));
    [~, idx2] = min(abs(centers-bin1Percent(3)));
    [~, idx3] = min(abs(centers-bin1Percent(4)));
    [~, idxEnd] = min(abs(centers-bin1Percent(end)));
    if (isempty(widths))
        k = 1;
    else
        k = size(widths,1)+1;
    end
    widths(k,:) = centers(idxEnd)-[centers(idx1) centers(idx2) centers(idx3)];
    if (k == 1)
        diffWidths = nan;
    else
        diffWidths = diff(widths,1)./(widths(end,:));
    end

    histogramData.x = centers;
    histogramData.y = counts;

    % Termination condition
    cond1 = bin1Percent(1) > 0.85;
    cond2 = (any(diffWidths(:) > 0.1)) && (~any((isnan(diffWidths(:)))));
    if (cond1 && cond2)
        fprintf(2,'Should terminate here\n');
        terminateNow = true;
    else
        terminateNow = false;
    end
        
end

function plotMeshQuality(figNo,histogramData, bin1Percent, iterationsHistory)
    if (isempty(figNo))
        hFig = figure(10); 
        subplotIndex = mod(numel(iterationsHistory)-1,12)+1;
        if (subplotIndex == 1)
            clf;
            set(hFig, 'Position', [10 10 820 930]);
        end
        
        rows = 5; cols = 3;
        posVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', rows, ...
           'colsNum', cols, ...
           'heightMargin',  0.07, ...
           'widthMargin',    0.03, ...
           'leftMargin',     0.03, ...
           'rightMargin',    0.01, ...
           'bottomMargin',   0.03, ...
           'topMargin',      0.02);
        row = 1+floor((subplotIndex-1)/cols);
        col = 1+mod((subplotIndex-1),cols);
        subplot('Position', posVectors(row,col).v);
    end
 
    qLims = [0.6 1.005]; 
    bar(histogramData.x,histogramData.y,1); hold on;
    plot(bin1Percent(1)*[1 1], [0 max(histogramData.y)], 'r-', 'LineWidth', 1.5);
    plot(bin1Percent(end)*[1 1], [0 max(histogramData.y)], 'c-', 'LineWidth', 1.5);
    plot(bin1Percent(2)*[1 1], [0 max(histogramData.y)], 'k-',  'LineWidth', 1.5);
    plot(bin1Percent(3)*[1 1], [0 max(histogramData.y)], 'k-', 'LineWidth', 1.5);
    plot(bin1Percent(4)*[1 1], [0 max(histogramData.y)], 'k-', 'LineWidth', 1.5);
    set(gca, 'XLim', qLims, 'YLim', [0 max(histogramData.y)], ...
        'XTick', [0.6:0.05:1.0],  'XTickLabel', {'.6', '', '.7', '', '.8', '', '.9', '', '1.'}, ...
        'YTickLabel', {}, 'FontSize', 12);
    grid on
    xlabel('hex-index $\left(\displaystyle 2 r_{ins} / r_{cir} \right)$', 'Interpreter', 'latex', 'FontSize', 12);
    if (isempty(figNo))
        title(sprintf('iteration:%d', iterationsHistory(end)))
        drawnow;
    end
    
    if (isempty(figNo))
        figure(11); hold on;
    end
    
end

function plotMovementSequence(figNo, maxMovements, dTolerance)
    if (isempty(figNo))
        figure(11); clf;
    end
    
    if (numel(maxMovements) < 10) 
        markerSize = 12;
    elseif (numel(maxMovements) < 50)
        markerSize = 10;
    elseif (numel(maxMovements) < 100)
        markerSize = 8;
    elseif (numel(maxMovements) < 500)
        markerSize = 6;
    else
        markerSize = 4;
    end
    
    plot(1:numel(maxMovements), maxMovements, 'ko-', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', markerSize);
    hold on;
    plot([1 numel(maxMovements)], dTolerance*[1 1], 'r-', 'LineWidth', 1.5);
    set(gca, 'YLim', [dTolerance*0.5 max(maxMovements)], 'YScale', 'log', 'FontSize', 16);
    xlabel('iteration');
    ylabel('median movement', 'FontSize', 16)
end


function plotMosaic(hFig, rfPositions, triangleIndices, maxMovements,  reTriangulationIterations, widths, histogramData, bin1Percent,  dTolerance, fovDegs, micronsPerDeg)

    eccDegs = (sqrt(sum(rfPositions.^2, 2)))/micronsPerDeg;
    idx = find(eccDegs <= min([1 fovDegs])/2);
    %idx = 1:size(rfPositions,1);
    
    if (isempty(hFig))
        hFig = figure(1);
        set(hFig, 'Position', [10 10 1596 1076]);
    end
    
    clf;
    subplot(2,3,[1 2 4 5]);
    
    if (1==1)
        plotTriangularizationGrid = true;
        if (plotTriangularizationGrid)
            visualizeLatticeState(rfPositions, triangleIndices);
        end

        plot(rfPositions(idx,1), rfPositions(idx,2), 'r.');
        maxPos = max(max(abs(rfPositions(idx,:))));
        set(gca, 'XLim', maxPos*[-1 1], 'YLim', maxPos*[-1 1], 'FontSize', 16);
        axis 'square'
    end
      
    subplot(2,3,3);
    yyaxis left
    plotMovementSequence(hFig, maxMovements, dTolerance);
    title(sprintf('Iteration: %d', reTriangulationIterations(end)));
    yyaxis right
    if (~isnan(widths))
        plot(reTriangulationIterations(2:end), widths(:,1), 'rs-', 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.0, 'MarkerSize', 10); hold on
        plot(reTriangulationIterations(2:end), widths(:,2), 'rs-', 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.0, 'MarkerSize', 10);
        plot(reTriangulationIterations(2:end), widths(:,3), 'rs-', 'MarkerFaceColor', [1 0.5 0.5], 'LineWidth', 1.0, 'MarkerSize', 10);
    end
    set(gca, 'YLim', [-1.5 0.1]);
     
    subplot(2,3,6);
    plotMeshQuality(hFig,histogramData, bin1Percent, []);
    drawnow
end

function q = computeQuality(rfLocs, triangles)
    %% q:   (a+b-c)(a+c-b)(b+c-a) / (abc), which is maximal to 1 when a=b=c

    trianglesNum = size(triangles,1);
    X = rfLocs(:,1);
    Y = rfLocs(:,2);
    
    q = zeros(1,trianglesNum);
    for triangleIndex = 1:trianglesNum
        for node = 1:3
            x(node) = X(triangles(triangleIndex,node));
            y(node) = Y(triangles(triangleIndex,node));
        end 
        aLength = sqrt((x(1)-x(2))^2 + (y(1)-y(2))^2);
        bLength = sqrt((x(1)-x(3))^2 + (y(1)-y(3))^2);
        cLength = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
        q(triangleIndex) = (bLength+cLength-aLength)*(cLength+aLength-bLength)*(aLength+bLength-cLength)/(aLength*bLength*cLength);
    end
end


function visualizeLatticeState(rfPositions, triangleIndices)
    x = rfPositions(:,1);
    y = rfPositions(:,2);
    
    xx = []; yy = [];
    for triangleIndex = 1:size(triangleIndices, 1)
        coneIndices = triangleIndices(triangleIndex, :);
        xCoords = x(coneIndices);
        yCoords = y(coneIndices);
        for k = 1:numel(coneIndices)
            xx = cat(2, xx, xCoords);
            yy = cat(2, yy, yCoords);
        end
    end
    
    patch(xx, yy, [0 0 1], 'EdgeColor', [1 0 0], ...
        'EdgeAlpha', 0.8, 'FaceAlpha', 0.0, ...
        'FaceColor', [0.99 0.99 0.99], 'LineWidth', 1.0, ...
        'LineStyle', '-', 'Parent', gca); 
    hold on;
end


function inputVal = GetWithDefault(prompt,defaultVal)
    if (ischar(defaultVal))
        inputVal = input(sprintf([prompt ' [%s]: '],defaultVal),'s');
    else
        inputVal = input(sprintf([prompt ' [%g]: '],defaultVal));
    end
    if (isempty(inputVal))
        inputVal = defaultVal;
    end

end
