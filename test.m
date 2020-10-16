%%
absent = max(0,cat(1,LFR{1,:}));
present2 = max(0,cat(1,LFR{4,:}));
present10 = max(0,cat(1,LFR{7,:}));

mAbsent = mean(mean(absent,1),3);
sdAbsent = std(mean(absent,1),0,3);
mPresent2 = mean(mean(present2,1),3);
sdPresent2 = std(mean(present2,1),0,3);
mPresent10 = mean(mean(present10,1),3);
sdPresent10 = std(mean(present10,1),0,3);

%%figure; hold on;
subplot(2,2,1); hold on; cla;
h(1) = plot( time{1}, mAbsent, 'b', 'linewidth', 1, 'displayname', 'Absent' );
h(2) = plot( time{1}, mPresent2, 'g', 'linewidth', 1, 'displayname', '2 cpd' );
h(3) = plot( time{1}, mPresent10, 'r', 'linewidth', 1, 'displayname', '10 cpd' );
m = mAbsent;
sem = sdAbsent / sqrt(30);
fill( [time{1}, time{1}(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], 'b', 'LineStyle', 'none', 'FaceAlpha', 0.4 );
m = mPresent2;
sem = sdPresent2 / sqrt(30);
fill( [time{1}, time{1}(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], 'g', 'LineStyle', 'none', 'FaceAlpha', 0.4 );
m = mPresent10;
sem = sdPresent10 / sqrt(30);
fill( [time{1}, time{1}(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], 'r', 'LineStyle', 'none', 'FaceAlpha', 0.4 );
legend(h, 'location', 'northeast');
set(gca, 'xlim', [-50 500], 'linewidth', 2, 'fontsize', 18)
title('Population Mean');
xlabel('Time from saccade off (ms)');
ylabel('Firing rate');

%%
subplot(2,2,2); hold on; cla;
colors = {'b', 'g', 'r'};
iTick = find(time{1} == 0);
accFR = {absent(:,iTick:end,:), present2(:,iTick:end,:), present10(:,iTick:end,:)};
for(k = 1 : 3)
	for( t = 2 : size(accFR{1},2) )
		accFR{k}(:,t,:) = accFR{k}(:,t-1,:) + accFR{k}(:,t,:);
	end
end
for(k = 1 : 3)
	m = mean(mean(accFR{k},1),3);
	sem = std(mean(accFR{k},1),0,3) / sqrt(size(accFR{k},3));
	plot( time{1}(iTick:end), m, 'color', colors{k}, 'linewidth', 1 );
	fill( [time{1}(iTick:end), time{1}(end:-1:iTick)], [m-sem, m(end:-1:1)+sem(end:-1:1)], colors{k}, 'LineStyle', 'none', 'FaceAlpha', 0.4 );
end
set(gca, 'xlim', [0 500], 'linewidth', 2, 'fontsize', 18)
title('Accumulated Population Mean');
xlabel('Time from saccade off (ms)');
ylabel('Accumulated firing rate');

%%
subplot(2,2,3); hold on; cla;
for(k = 1 : 3)
	meanFR{k} = accFR{k} ./ (1:size(accFR{k},2));
end
for(k = 1 : 3)
	m = mean(mean(meanFR{k},1),3);
	sem = std(mean(meanFR{k},1),0,3) / sqrt(size(meanFR{k},3));
	plot( time{1}(iTick:end), m, 'color', colors{k}, 'linewidth', 1 );
	fill( [time{1}(iTick:end), time{1}(end:-1:iTick)], [m-sem, m(end:-1:1)+sem(end:-1:1)], colors{k}, 'LineStyle', 'none', 'FaceAlpha', 0.4 );
end
set(gca, 'xlim', [0 500], 'linewidth', 2, 'fontsize', 18)
title('Average of Population Mean');
xlabel('Time from saccade off (ms)');
ylabel('Average firing rate');


%%
subplot(2,2,4); hold on; cla;
for(k = 2 : 3)
	mDiffAccFR{k} = mean(mean(accFR{k},1),3) - mean(mean(accFR{1},1),3);
	sdDiffAccFR{k} = sqrt( std(mean(accFR{k},1),0,3).^2 + std(mean(accFR{1},1),0,3).^2 );
end
for(k = 2 : 3)
	m = mDiffAccFR{k};
	sem = sdDiffAccFR{k};
	plot( time{1}(iTick:end), m, 'color', colors{k}, 'linewidth', 1 );
	fill( [time{1}(iTick:end), time{1}(end:-1:iTick)], [m-sem, m(end:-1:1)+sem(end:-1:1)], colors{k}, 'LineStyle', 'none', 'FaceAlpha', 0.4 );
end
set(gca, 'xlim', [0 500], 'linewidth', 2, 'fontsize', 18)
title('Difference from Absent for Accumulated Population Mean');
xlabel('Time from saccade off (ms)');
ylabel('Accumulated firing rate');
return

%%
iCond = 4; 1; 7;
Eccs = [0 4 8];
nTrials = 30;
nNeurons = 3;
obj = encoder;
contrast = 0.5;
iEcc = find( conditions(iCond).eccentricity == Eccs );
idx = trialsIdx{iCond};    %fprintf('nTrials = %d\n', size(idx,2)); continue;
idx = idx( 1 : min(nTrials,end) );
tMax = round(max( [trials(idx).saccadeOff] + round(conditions(iCond).duration/1000*[trials(idx).sRate]) - [trials(idx).(alignEvent)] ));	% in samples
tMin = round(min( [trials(idx).saccadeOn] - round(0.150*[trials(idx).sRate]) - [trials(idx).(alignEvent)] ));				% in samples
k = 1;
% tic;
index = tMin + trials(idx(k)).(alignEvent) : min( size(trials(idx(k)).x.position,2), tMax + trials(idx(k)).(alignEvent) );
x  = trials(idx(k)).x.position(index) / 60;
y  = trials(idx(k)).y.position(index) / 60;

% before saccade
eyeIdx{1} = 1 : trials(idx(k)).saccadeOn - trials(idx(k)).(alignEvent) - tMin;
s1 = size(eyeIdx{1},2);

% during saccade
eyeIdx{2} = eyeIdx{1}(end) + ( 1 : trials(idx(k)).saccadeOff - trials(idx(k)).saccadeOn + 1 );
s2 = size(eyeIdx{2},2);

% after saccade
eyeIdx{3} = eyeIdx{2}(end)+1 : size(x,2);
s3 = size(eyeIdx{3},2);

if( strcmpi( stabilize, 'drift' ) )
	x(eyeIdx{1}) = x(eyeIdx{1}(end));
	y(eyeIdx{1}) = y(eyeIdx{1}(end));
	x(eyeIdx{3}) = x(eyeIdx{3}(1));
	y(eyeIdx{3}) = y(eyeIdx{3}(1));
elseif( strcmpi( stabilize, 'full' ) )
	x(:) = x(eyeIdx{3}(1));
	y(:) = y(eyeIdx{3}(1));
end

[noise, inputX, inputY] = obj.LoadNoise( fullfile(dataFolder, trials(idx(k)).backgroundImage), trials(idx(k)).pixelAngle/60 );
grating = obj.GenerateGrating( conditions(iCond).eccentricity, trials(idx(k)).gratingWidth, conditions(iCond).sf, trials(idx(k)).phase, trials(idx(k)).pixelAngle/60 );
noise = noise * trials(idx(k)).backgroundContrast;
grating = grating * contrast * conditions(iCond).present;
bg = 1;
stimulus = noise + grating;
avContrast = trials(idx(k)).backgroundContrast + contrast * conditions(iCond).present;

iL = 3;
tic;
for( m = 1 : 3 )
	xIdx = min(x(eyeIdx{m})) + min(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nNeurons),1)) - 4 <= inputX & inputX <= max(x(eyeIdx{m})) + max(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nNeurons),1)) + 4; %720 : 1460;
    yIdx = min(y(eyeIdx{m})) + min(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nNeurons),2)) - 4 <= inputY & inputY <= max(y(eyeIdx{m})) + max(obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nNeurons),2)) + 4; %300 : 780;

	[sfr{iCond,iL}(:,eyeIdx{m},k), sfr_c{iCond,iL}(:,eyeIdx{m},k), sfr_s{iCond,iL}(:,eyeIdx{m},k)] = obj.SpatialModel.LinearResponse( stimulus(yIdx,xIdx), inputX(xIdx), inputY(yIdx), x(eyeIdx{m}), y(eyeIdx{m}), obj.layers(iL).sRFParams(cellIdx{iL,iEcc}(1:nNeurons)), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nNeurons),1), obj.layers(iL).locations(cellIdx{iL,iEcc}(1:nNeurons),2));
end
if( lower(obj.layers(iL).name(1)) == 'm' )
	lfr{iCond,iL}(:,:,k) = obj.TemporalModel.LinearResponse( obj.layers(iL).name, obj.layers(iL).tRFParams(cellIdx{iL,iEcc}(1:nNeurons)), avContrast, trials(idx(k)).sRate, sfr{iCond,iL}(:,:,k) );
else
	lfr{iCond,iL}(:,:,k) = obj.TemporalModel.LinearResponse( obj.layers(iL).name, obj.layers(iL).tRFParams(cellIdx{iL,iEcc}(1:nNeurons)), avContrast, trials(idx(k)).sRate, sfr_c{iCond,iL}(:,:,k), sfr_s{iCond,iL}(:,:,k) );
end
fprintf('Ecc=%d, SF=%d, Dur=%d, Present=%d, iL=%d, iTrials=%d/%d, t=%f\n', conditions(iCond).eccentricity, conditions(iCond).sf, conditions(iCond).duration, conditions(iCond).present, iL, k, size(idx,2), toc);

return;


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