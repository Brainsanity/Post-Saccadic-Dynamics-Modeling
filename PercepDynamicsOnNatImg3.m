%%
load('../../Data/Simulated Activities/SacDB/UG - Noise & Grating Simulated Separately/figures - withInternalNoise/uni-no_bias-fa25_hit75 - durOffset=0/PerformanceData.mat');
iL = 5;
sen2 = shiftdim(Sensitivities(iL,1,:,:), 2)';
sen2SD = shiftdim(SensitivitiesSTD(iL,1,:,:), 2)';
sen10 = shiftdim(Sensitivities(iL,2,:,:), 2)';
sen10SD = shiftdim(SensitivitiesSTD(iL,2,:,:), 2)';


%%%% Video showing perceptual dynamics over natural scenes
%%
% imgName = 'DinoShower'; dx = 3.242; dy = -1.458; blurSD = 20; min2 = 0.2;	min10 = 0;
imgName = 'Crab'; dx = 0; dy = 2.4; blurSD = 40; min2 = 0.4; min10 = 0;
imgName = 'Elephants';  dx = 1.2; dy = 0.4; blurSD = 20; min2 = 0.4; min10 = 0;
imgFile = sprintf('../../manuscript/FB Renewal 2021/Perceptual Dynamics over Natural Scenes/%s.jpg', imgName);
img = imread(imgFile);
img = double(img) / 255;
img = img(:,:,1) * 4/13 + img(:,:,2) * 8/13 + img(:,:,3) * 1/13;		% change to gray scale
pixPerArcmin = 1;
figure;
imshow(img);

%% assemble an eye trace
imgX = (-size(img,2)/2 + 0.5 : size(img,2)/2 - 0.5) / pixPerArcmin  / 60;
imgY = (-size(img,1)/2 + 0.5 : size(img,1)/2 - 0.5) / pixPerArcmin  / 60;
cutSF = 3;		% cutoff spatial frequency (cpd)
figure;
imshow(img, 'XData', imgX, 'YData', imgY); hold on;
set(gca, 'visible', 'on', 'nextplot', 'add');

trials = EmpiricalBox.LoadSacDB();
iTrials = 1:3;
eyeX = trials(iTrials(1)).x.position;
eyeY = trials(iTrials(1)).y.position;
sacsOn = trials(iTrials(1)).saccadeOn;% - 50;
sacsOff = trials(iTrials(1)).saccadeOff;
for(iTrial = iTrials(2:end))
	% sacsOn = [sacsOn, trials(iTrial).saccadeOn - 50 + size(eyeX,2)];
	sacsOn = [sacsOn, trials(iTrial).saccadeOn + size(eyeX,2)];
	sacsOff = [sacsOff, trials(iTrial).saccadeOff + size(eyeX,2)];
	eyeX = [eyeX, trials(iTrial).x.position - trials(iTrial).x.position(1) + eyeX(end)];
	eyeY = [eyeY, trials(iTrial).y.position - trials(iTrial).y.position(1) + eyeY(end)];
end
eyeX = eyeX / 60;
eyeY = eyeY / 60;

eyeX = eyeX + dx;
eyeY = eyeY + dy;

plot(eyeX, eyeY, 'r.');

figure; hold on;
plot([eyeX; eyeY]');
plot([sacsOn; sacsOn], ylim' * ones(size(sacsOn)), 'g--');
plot([sacsOff; sacsOff], ylim' * ones(size(sacsOff)), 'r--');


%% get dynamics of sensitivity
Eccs = unique([conditions.eccentricity, 0.5]);
nEccs = size(Eccs,2);
sen2 = [sen2(:,1), sen2(:,1), sen2(:,2:end)];	% add ecc = 0.5 as equal to ecc =0
sen10 = [sen10(:,1), sen10(:,1), sen10(:,2:end)];	% add ecc = 0.5 as equal to ecc =0
sen2(1,:) = sen2(end,:) * min2;
sen10(1,:) = sen10(end,:) * min10;
sen2 = sen2 / mean(sen2(end-4:end, 1));			% normalize to the late period at the fovea
sen10 = sen10 / mean(sen10(end-4:end, 1));		% normalize to the late period at the fovea
% for(iEcc = 1:nEccs)
% 	sen10(:,iEcc) = linspace(sen10(1,iEcc), sen10(end,iEcc), size(sen10,1));		% replace the modeled dynamics with linear dynamics
% end
figure;
subplot(1,2,1);
plot(durs, sen2);
xlabel('Time (ms)');
ylabel('Sensitivity');
title('2 CPD');
subplot(1,2,2);
plot(durs, sen10);
title('10 CPD');
xlabel('Time (ms)');
ylabel('Sensitivity');

%%
tTicksL = 0 : durs(end);
nTicksL = size(tTicksL,2);
sL = interp2(Eccs, durs, sen2, Eccs', tTicksL)';

slowRate = 1;
tTicksH = 0 : durs(end) * slowRate;
nTicksH = size(tTicksH,2);
sH = interp2(Eccs, durs * slowRate, sen10, Eccs', tTicksH)';		% slowdown the dynamics

senLow = ones(nEccs, size(eyeX,2));
senHigh = senLow;
for(iEcc = 1 : nEccs)
	% fixation before 1st saccade
	senLow(iEcc, 1 : sacsOn(1)-1) = sL(iEcc,end);
	senHigh(iEcc, 1 : sacsOn(1)-1) = sH(iEcc,end);

	% saccade then post-saccade
	for(iSac = 1 : size(sacsOn,2))
		% during saccade, roughly following Fig. 1 of Idrees et al., Nature Communications, 2020, but with different minima
		downT = 30;	% 70;
		delay = 10;	% 20;
		senLow(iEcc, sacsOn(iSac) + (0:downT-1)) = senLow(iEcc, sacsOn(iSac)-1) * ((downT-1:-1:0)/(downT-1) * (1-min2) + min2);
		senHigh(iEcc, sacsOn(iSac) + (delay:downT-1)) = senHigh(iEcc, sacsOn(iSac)-1) * ((downT-1-delay:-1:0)/(downT-1-delay) * (1-min10) + min10);
		senLow(iEcc, sacsOn(iSac) + downT : sacsOff(iSac)) = min2 * senLow(iEcc, sacsOn(iSac)-1);
		senHigh(iEcc, sacsOn(iSac) + downT : sacsOff(iSac)) = min10 * senHigh(iEcc, sacsOn(iSac)-1);
		senHigh(iEcc, sacsOn(iSac) + (0:delay-1)) = senHigh(iEcc, sacsOn(iSac)-1);

		senLow(iEcc, sacsOff(iSac) + (1:nTicksL)) = sL(iEcc,:);
		senHigh(iEcc, sacsOff(iSac) + (1:nTicksH)) = sH(iEcc,:);
		if(iSac < size(sacsOn,2))
			senLow(iEcc, sacsOff(iSac)+nTicksL+1 : sacsOn(iSac+1)-1) = senLow(iEcc, sacsOff(iSac)+nTicksL);
			senHigh(iEcc, sacsOff(iSac)+nTicksH+1 : sacsOn(iSac+1)-1) = senHigh(iEcc, sacsOff(iSac)+nTicksH);
		else
			senLow(iEcc, sacsOff(iSac)+nTicksL+1 : end) = senLow(iEcc, sacsOff(iSac)+nTicksL);
			senHigh(iEcc, sacsOff(iSac)+nTicksH+1 : end) = senHigh(iEcc, sacsOff(iSac)+nTicksH);
		end

	end
end

figure; hold on;
colors = {'r', 'g', 'b'};
for(iEcc = 1:nEccs)
	plot(squeeze(senLow(iEcc, :)), '-', 'color', colors{mod(iEcc-1, size(colors,2))+1}, 'displayname', sprintf('Low | Ecc = %d', Eccs(iEcc)));
end
for(iEcc = 1:nEccs)
	plot(squeeze(senHigh(iEcc, :)), '--', 'color', colors{mod(iEcc-1, size(colors,2))+1}, 'displayname', sprintf('Low | Ecc = %d', Eccs(iEcc)));
end
legend;


%% filtering of the whole image
m = mean(mean(img,1),2);
imLow = imgaussfilt(img, blurSD) - m;
imHigh = img - imLow - m;
figure; imshow(imLow + m); title('Low SF');
figure; imshow(imHigh + m); title('High SF');

%% generate movie showing the perceptual dynamics over a natural scene
figure('NumberTitle', 'off', 'name', sprintf('Perceptual Dynamics over a Natural Image - SD = %d px', blurSD));
imshow(img, 'XData', imgX, 'YData', imgY); hold on;
set(gca, 'position', [0 0 1 1]);
hImg = get(gca, 'children');
hEye = plot(eyeX(1), eyeY(1), 'ro', 'linewidth', 1, 'markersize', 20);
hText = text(0, min(imgY), sprintf('t = %d ms', 1), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 20, 'color', 'r');
tFrames = 1 : 4 : size(eyeX,2);
nFrames = size(tFrames,2);

saveFolder = '../../Manuscript/FB Renewal 2021/Videos';
filename = fullfile(saveFolder, sprintf('Gray. Perceptual Dynamics over %s. SD = %d. nFrames = %d.px', imgName, blurSD, nFrames));
writerObj = VideoWriter(filename, 'MPEG-4');
open(writerObj);

for(iTick = 1:nFrames)
	t = tFrames(iTick);
	hEye.XData = eyeX(t);
	hEye.YData = eyeY(t);
	hImg.CData(:) = 0.5;
	hText.String = sprintf('t = %d ms', t);
	eccs = sqrt((imgX - eyeX(t)).^2 + (imgY' - eyeY(t)).^2);
	wLow = reshape(interp1(Eccs, senLow(:,t), eccs(:)', 'linear', NaN), size(eccs,1), []);
	wHigh = reshape(interp1(Eccs, senHigh(:,t), eccs(:)', 'linear', NaN), size(eccs,1), []);
	wLow(isnan(wLow)) = senLow(end,t);
	wHigh(isnan(wHigh)) = senHigh(end,t);

	hImg.CData = imLow .* wLow + imHigh .* wHigh + m;

	drawnow;
	writeVideo(writerObj, getframe(gcf));
end
close(writerObj);