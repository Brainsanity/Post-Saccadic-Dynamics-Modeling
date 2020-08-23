idx = trialsIdx{1};
tMin = time{1}(1);
tMax = time{1}(end);
contrast = 0.1;
iConds = [1 5 9]+4;
iL = 1;
% iEcc = find( conditions(iCond).eccentricity == [0, 4, 8, 12] );
nCells = 15;
alignEvent = 'flashOn';
dataFolder = '../../Data';
figure('color', 'w');
col = {'b', 'g', 'r'};
for( iEcc = 1 : 4 )
	iConds = [0 4 8] + iEcc;
	subplot(2,2,iEcc); hold on;
	for( iCond = iConds )
		for( k = 1 : 15 )
			x  = trials(idx(k)).x.position / 60;
			y  = trials(idx(k)).y.position / 60;

			% before grating
			eyeIdx1 = tMin + trials(idx(k)).(alignEvent) : trials(idx(k)).flashOn-1;
			s1 = size(eyeIdx1,2);

			% during grating
			eyeIdx2 = trials(idx(k)).flashOn : min( size(x,2), tMax + trials(idx(k)).(alignEvent) );
			s2 = size(eyeIdx2,2);

			[noise, inputX, inputY] = encoder.LoadNoise( fullfile(dataFolder, trials(idx(k)).backgroundImage), trials(idx(k)).pixelAngle/60 );
			grating = encoder.GenerateGrating( conditions(iCond).eccentricity, trials(idx(k)).gratingWidth, conditions(iCond).sf, trials(idx(k)).phase, trials(idx(k)).pixelAngle/60 );
			noise = noise * trials(idx(k)).backgroundContrast;
			grating = grating * contrast * conditions(iCond).present;
			bg = 1;

			avContrast = ((s1+s2)*trials(idx(k)).backgroundContrast + s2*contrast*conditions(iCond).present) / (s1+s2);
			% contrastF = [ zeros(1,s1), ones(1,s2) * contrast ] + trials(idx(k)).backgroundContrast;
			contrastF(:,:,k) = encoder.ComputeContrasts( noise+grating+bg,  inputX, inputY, encoder.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),1), encoder.layers(iL).locations(cellIdx{iL,iEcc}(1:nCells),2), [encoder.layers(iL).sRFParams(cellIdx{iL,iEcc}(1:nCells)).surroundRadii], x([eyeIdx1,eyeIdx2]), y([eyeIdx1,eyeIdx2]) );
		end
		plot( tMin:tMax, mean(contrastF-(contrastF(:,1,:)),3)', 'color', col{iCond==iConds} );
	end
	title( sprintf('Ecc = %d', conditions(iCond).eccentricity) );
	if(iEcc >= 3), xlabel('Time aligned to flashOn (ms)'); end
	if(mod(iEcc,2)), ylabel('Contrast within 6*std of cell surround'); end
	set(gca, 'fontsize', 20, 'linewidth', 2);
end
return;


%% prepare data
load(['../../Data/Mosaic_POn_Radius20.0deg_maxMovPrctile20.mat']);
rfPositions = rfPositions( abs(rfPositions(:,1))<=2 & abs(rfPositions(:,2))<=2, : );

ck = CronerKaplanRGCModel( 'dataSetToFit', 'raw', 'fitIntercept', false, 'randomizeCenterRadii', true, 'randomizeCenterSensitivities', true, 'randomizeSurroundRadii', true, 'randomizeSurroundSensitivities', true );
PRFParams = ck.SynthesizeRFParams(sqrt(sum(rfPositions.^2,2))','POn');


img = imread('demo.jpg');
img = single(img(end:-1:1,:,1));
inputX = (-275:275)/60;
inputY = (-310:309)/60;

load('ExampleEyeTrace.mat');
eyeX = eyeX/4;
eyeX = eyeX(1501:2500);
eyeY = eyeY(1501:2500);
eyeX = -0.15;
eyeY = -1.88;


%% spatial RF
fr = zeros( size(rfPositions,1), length(eyeX) );
fr_c = fr;
fr_s = fr;

n = 1000;
T = zeros( 1, ceil(size(rfPositions,1)/n) );
for( k = 0 : ceil(size(rfPositions,1)/n)-1 )
	fprintf('%d/%d\n', k, ceil(size(rfPositions,1)/n)-1 );
	idx = k*n+1 : min(k*n+n, size(rfPositions,1));
    tic;
	[fr(idx,:), fr_c(idx,:), fr_s(idx,:)] = ck.LinearResponse(img, inputX, inputY, eyeX, eyeY, PRFParams(idx), rfPositions(idx,1), rfPositions(idx,2));
    T(k+1) = toc;
end

%%	temporal RF
for( k = size(fr,1) : -1 : 1 )
	tfr_c(k,:) = RGC.TemporalLinearModel( fr_c(k,:), 'P', 'center' );
	tfr_s(k,:) = RGC.TemporalLinearModel( fr_s(k,:), 'P', 'surround' );
end
tfr = tfr_c + tfr_s;

%% Display activity
t = 1;
d = 1/60/2;
[x, y] = meshgrid( -2:d:2, -2:d:2 );
frMap = zeros(size(x));
counts = zeros(size(x));
idx = round( (rfPositions(:,1)' - -2 ) / d ) + 1;
idy = round( (rfPositions(:,2)' - -2 ) / d ) + 1;
for( k = find( idx >= 1 & idx <= size(x,1) & idy >= 1 & idy <= size(y,2) ) )
	frMap(idy(k),idx(k)) = frMap(idy(k),idx(k)) + fr(k,t);%max(0,fr(k,t));
	counts(idy(k),idx(k)) = counts(idy(k),idx(k)) + 1;
end
frMap(:) = max( frMap(:), prctile(frMap(:), 0.5) );
frMap(:) = min( frMap(:), prctile(frMap(:), 99.5) );
frMap = (frMap - min(frMap(:))) ./ (max(frMap(:)) - min(frMap(:)));
imshow( cat(3, frMap, zeros(size(frMap,1),size(frMap,2),2)) );
set( gca, 'ydir', 'normal' );