
%% prepare data
load(['../../Data/Mosaic_POn_Radius20.0deg_maxMovPrctile20.mat']);
rfPositions = rfPositions( abs(rfPositions(:,1))<=2 & abs(rfPositions(:,2))<=2, : );

ck = CronerKaplanRGCModel( 'dataSetToFit', 'raw', 'fitIntercept', false, 'randomizeCenterRadii', true, 'randomizeCenterSensitivities', true, 'randomizeSurroundRadii', true, 'randomizeSurroundSensitivities', true );
PRFParams = ck.SynthesizeRFParams(sqrt(sum(rfPositions.^2,2))','P');


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
	PRFParams2 = PRFParams;
	PRFParams2.centerRadii = PRFParams2.centerRadii(idx);
	PRFParams2.centerPeakSensitivities = PRFParams2.centerPeakSensitivities(idx);
	PRFParams2.surroundRadii = PRFParams2.surroundRadii(idx);
	PRFParams2.surroundPeakSensitivities = PRFParams2.surroundPeakSensitivities(idx);
    tic;
	[fr(idx,:), fr_c(idx,:), fr_s(idx,:)] = ck.LinearResponse(img, inputX, inputY, eyeX, eyeY, PRFParams2, rfPositions(idx,1), rfPositions(idx,2));
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
frMap(:) = max( frMap(:), prctile(frMap(:), 0.1) );
frMap(:) = min( frMap(:), prctile(frMap(:), 99.9) );
frMap = (frMap - min(frMap(:))) ./ (max(frMap(:)) - min(frMap(:)));
imshow( cat(3, frMap, zeros(size(frMap,1),size(frMap,2),2)) );
set( gca, 'ydir', 'normal' );