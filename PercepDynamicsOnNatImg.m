%%%% Video showing perceptual dynamics over natural scenes
%%
imgFile = '../../2021.06. FB Grants Renewal/Perceptual Dynamics over Natural Scenes/DinoShower.jpg';
img = imread(imgFile);
img = double(img) / 255;
pixPerArcmin = 1;
figure;
imshow(img);

%% filtering of the whole image
cutSF = 3;		% cutoff spatial frequency (cpd)
[nY, nX, ~] = size(img);
sfX = (-nX/2 : nX/2-1) / nX * (60/pixPerArcmin);
sfY = (-nY/2 : nY/2-1) / nY * (60/pixPerArcmin);
sf2 = sfX.^2 + sfY'.^2;
imHigh = zeros(size(img));
imLow = zeros(size(img));
for(k = 1:3)
	f = fftshift(fft2(img(:,:,k)));
	imHigh(:,:,k) = real(ifft2(ifftshift(f .* (sf2 > cutSF^2))));
	imLow(:,:,k) = real(ifft2(ifftshift(f .* (sf2 < cutSF^2))));
end
figure; imshow(imLow); title(sprintf('[0 %d]', cutSF));
figure; imshow(imHigh); title(sprintf('[%d Inf)', cutSF));

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
sacsOn = trials(iTrials(1)).saccadeOn - 50;
sacsOff = trials(iTrials(1)).saccadeOff;
for(iTrial = iTrials(2:end))
	sacsOn = [sacsOn, trials(iTrial).saccadeOn - 50 + size(eyeX,2)];
	sacsOff = [sacsOff, trials(iTrial).saccadeOff + size(eyeX,2)];
	eyeX = [eyeX, trials(iTrial).x.position - trials(iTrial).x.position(1) + eyeX(end)];
	eyeY = [eyeY, trials(iTrial).y.position - trials(iTrial).y.position(1) + eyeY(end)];
end
eyeX = eyeX / 60;
eyeY = eyeY / 60;
plot(eyeX, eyeY, 'r.');

figure; hold on;
plot([eyeX; eyeY]');
plot([sacsOn; sacsOn], ylim' * ones(size(sacsOn)), 'g--');
plot([sacsOff; sacsOff], ylim' * ones(size(sacsOff)), 'r--');


%% get dynamics of sensitivity
Eccs = [0 2 4]; nEccs = size(Eccs,2);
senLow = ones(size(Eccs,2), size(eyeX,2));
senHigh = senLow;

% post-saccadic dynamics
tTicks = 0 : durs(end);
nTicks = size(tTicks,2);
min2 = 0;	% 0.2
min10 = 0;	% 0.25
sen2(1,:) = sen2(end,:) * min2;
sen10(1,:) = sen10(end,:) * min10;
for(iEcc = nEccs:-1:1)
	sL(iEcc,:) = interp1(durs, sen2(:, Eccs(iEcc) == unique([conditions.eccentricity])), tTicks, 'linear');
	sH(iEcc,:) = interp1(durs, sen10(:, Eccs(iEcc) == unique([conditions.eccentricity])), tTicks, 'linear');
end
m = max([sL(:); sH(:)]);
sL = sL / m;
sH = sH / m;

for(iEcc = 1:nEccs)
	% fixation before 1st saccade
	senLow(iEcc, 1 : sacsOn(1)-1) = sL(iEcc,end);
	senHigh(iEcc, 1 : sacsOn(1)-1) = sH(iEcc,end);

	% saccade then post-saccade
	for(iSac = 1 : size(sacsOn,2))
		% during saccade, roughly following Fig. 1 of Idrees et al., Nature Communications, 2020
		senLow(iEcc, sacsOn(iSac) + (0:69)) = senLow(iEcc, sacsOn(iSac)-1) * ((69:-1:0)/69 * (1-min2) + min2);
		senHigh(iEcc, sacsOn(iSac) + (20:69)) = senHigh(iEcc, sacsOn(iSac)-1) * ((49:-1:0)/49 * (1-min10) + min10);
		senLow(iEcc, sacsOn(iSac) + 70 : sacsOff(iSac)) = min2 * senLow(iEcc, sacsOn(iSac)-1);
		senHigh(iEcc, sacsOn(iSac) + 70 : sacsOff(iSac)) = min10 * senHigh(iEcc, sacsOn(iSac)-1);
		senHigh(iEcc, sacsOn(iSac) + (0:19)) = senHigh(iEcc, sacsOn(iSac)-1);

		senLow(iEcc, sacsOff(iSac) + (1:nTicks)) = sL(iEcc,:);
		senHigh(iEcc, sacsOff(iSac) + (1:nTicks)) = sH(iEcc,:);
		if(iSac < size(sacsOn,2))
			senLow(iEcc, sacsOff(iSac)+nTicks+1 : sacsOn(iSac+1)-1) = senLow(iEcc, sacsOff(iSac)+nTicks);
			senHigh(iEcc, sacsOff(iSac)+nTicks+1 : sacsOn(iSac+1)-1) = senHigh(iEcc, sacsOff(iSac)+nTicks);
		else
			senLow(iEcc, sacsOff(iSac)+nTicks+1 : end) = senLow(iEcc, sacsOff(iSac)+nTicks);
			senHigh(iEcc, sacsOff(iSac)+nTicks+1 : end) = senHigh(iEcc, sacsOff(iSac)+nTicks);
		end

	end
end

figure; hold on;
colors = {'r', 'g', 'b'};
for(iEcc = 1:nEccs)
	plot(senLow(iEcc, :), '-', 'color', colors{iEcc}, 'displayname', sprintf('Low | Ecc = %d', Eccs(iEcc)));
end
for(iEcc = 1:nEccs)
	plot(senHigh(iEcc, :), '--', 'color', colors{iEcc}, 'displayname', sprintf('High | Ecc = %d', Eccs(iEcc)));
end
legend;


%% generate movie showing the perceptual dynamics over a natural scene
rApt = 1;	% radius of aperture (deg)

saveFolder = '../../Manuscript/FB Renewal 2021/Videos';
filename = fullfile(saveFolder, 'Perceptual Dynamics over Natural Image');
writerObj = VideoWriter(filename, 'MPEG-4');
open(writerObj);

figure('NumberTitle', 'off', 'name', 'Perceptual Dynamics over a Natural Image');
imshow(img, 'XData', imgX, 'YData', imgY); hold on;
set(gca, 'position', [0 0 1 1]);
hImg = get(gca, 'children');
hText = text(0, min(imgY), sprintf('t = %d ms', t), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 20, 'color', 'r');
for(t = 1:size(eyeX,2))
	hImg.CData(:) = 0.5;
	hText.String = sprintf('t = %d ms', t);
	for(iEcc = 1:nEccs)
		idx = abs(eyeX(t) + Eccs(iEcc) - imgX) <= rApt;
		idy = abs(eyeY(t) - imgY) <= rApt;
		imLocal = img(idy,idx,:);
		[nY, nX, ~] = size(imLocal);
		sfX = (-nX/2 : nX/2-1) / nX * (60/pixPerArcmin);
		sfY = (-nY/2 : nY/2-1) / nY * (60/pixPerArcmin);
		sf2 = sfX.^2 + sfY'.^2;
		for(k = 1:3)
			f = fftshift(fft2(imLocal(:,:,k)));
			f(sf2 > cutSF^2) = f(sf2 > cutSF^2) * senHigh(iEcc,t);
			f(sf2 <= cutSF^2 & sf2 > 0) = f(sf2 <= cutSF^2 & sf2 > 0) * senLow(iEcc,t);
			hImg.CData(idy,idx,k) = real(ifft2(ifftshift(f)));
		end
	end
	drawnow;
	writeVideo(writerObj, getframe(gcf));
end
close(writerObj);