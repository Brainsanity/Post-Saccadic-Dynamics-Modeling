%%%% Video showing perceptual dynamics over natural scenes
%%
imgFile = '../../manuscript/FB Renewal 2021/Perceptual Dynamics over Natural Scenes/DinoShower.jpg';
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
	imLow(:,:,k) = real(ifft2(ifftshift(f .* (sf2 < cutSF^2 & sf2 > 0))));
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
patchW = 1;	% width of each patch (deg)
patchLocX = -4 : patchW : 4;
patchLocY = (-4 : patchW : 4)';

% post-saccadic dynamics
tTicks = 0 : durs(end);
nTicks = size(tTicks,2);
min2 = 0;	% 0.2
min10 = 0;	% 0.25
sen2(1,:) = sen2(end,:) * min2;
sen10(1,:) = sen10(end,:) * min10;
% for(iEcc = size(sen2,2) : -1 : 1)
% 	sen2_(iEcc,:) = interp1(durs, sen2(:,iEcc), tTicks, 'linear');
% 	sen10_(iEcc,:) = interp1(durs, sen10(:,iEcc), tTicks, 'linear');
% end
% for(iTick = nTicks:-1:1)
% 	sL()
% end
sL = interp2(unique([conditions.eccentricity]), durs, sen2, reshape(sqrt(patchLocY.^2 + patchLocX.^2), [], 1), tTicks);
sH = interp2(unique([conditions.eccentricity]), durs, sen10, reshape(sqrt(patchLocY.^2 + patchLocX.^2), [], 1), tTicks);
sL = reshape(shiftdim(sL,1), size(patchLocY,1), size(patchLocX, 2), []);
sH = reshape(shiftdim(sH,1), size(patchLocY,1), size(patchLocX, 2), []);

m = max([sL(:); sH(:)]);
sL = sL / m;
sH = sH / m;

senLow = ones(size(patchLocY,1), size(patchLocX,2), size(eyeX,2));
senHigh = senLow;
for(iX = 1 : size(patchLocX,2))
	for(iY = 1 : size(patchLocY,1))
		% fixation before 1st saccade
		senLow(iY, iX, 1 : sacsOn(1)-1) = sL(iY, iX,end);
		senHigh(iY, iX, 1 : sacsOn(1)-1) = sH(iY, iX,end);

		% saccade then post-saccade
		for(iSac = 1 : size(sacsOn,2))
			% during saccade, roughly following Fig. 1 of Idrees et al., Nature Communications, 2020
			senLow(iY, iX, sacsOn(iSac) + (0:69)) = senLow(iY, iX, sacsOn(iSac)-1) * ((69:-1:0)/69 * (1-min2) + min2);
			senHigh(iY, iX, sacsOn(iSac) + (20:69)) = senHigh(iY, iX, sacsOn(iSac)-1) * ((49:-1:0)/49 * (1-min10) + min10);
			senLow(iY, iX, sacsOn(iSac) + 70 : sacsOff(iSac)) = min2 * senLow(iY, iX, sacsOn(iSac)-1);
			senHigh(iY, iX, sacsOn(iSac) + 70 : sacsOff(iSac)) = min10 * senHigh(iY, iX, sacsOn(iSac)-1);
			senHigh(iY, iX, sacsOn(iSac) + (0:19)) = senHigh(iY, iX, sacsOn(iSac)-1);

			senLow(iY, iX, sacsOff(iSac) + (1:nTicks)) = sL(iY, iX,:);
			senHigh(iY, iX, sacsOff(iSac) + (1:nTicks)) = sH(iY, iX,:);
			if(iSac < size(sacsOn,2))
				senLow(iY, iX, sacsOff(iSac)+nTicks+1 : sacsOn(iSac+1)-1) = senLow(iY, iX, sacsOff(iSac)+nTicks);
				senHigh(iY, iX, sacsOff(iSac)+nTicks+1 : sacsOn(iSac+1)-1) = senHigh(iY, iX, sacsOff(iSac)+nTicks);
			else
				senLow(iY, iX, sacsOff(iSac)+nTicks+1 : end) = senLow(iY, iX, sacsOff(iSac)+nTicks);
				senHigh(iY, iX, sacsOff(iSac)+nTicks+1 : end) = senHigh(iY, iX, sacsOff(iSac)+nTicks);
			end

		end
	end
end

figure; hold on;
colors = {'r', 'g', 'b'};
for(iX = 1 : size(patchLocX,2))
	for(iY = 1 : size(patchLocY,1))
		plot(squeeze(senLow(iY, iX, :)), '-', 'color', colors{mod((iX-1)*size(patchLocY,1)+iY, nEccs)+1}, 'displayname', sprintf('Low | X = %d, Y = %d', patchLocX(iX), patchLocY(iY)));
	end
end
for(iX = 1 : size(patchLocX,2))
	for(iY = 1 : size(patchLocY,1))
		plot(squeeze(senHigh(iY, iX, :)), '--', 'color', colors{mod((iX-1)*size(patchLocY,1)+iY, nEccs)+1}, 'displayname', sprintf('Low | X = %d, Y = %d', patchLocX(iX), patchLocY(iY)));
	end
end
legend;


%% generate movie showing the perceptual dynamics over a natural scene
% rApt = 1;	% radius of aperture (deg)


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
	for(iX = 1 : size(patchLocX,2))
		for(iY = 1 : size(patchLocY,1))
			idx = abs(eyeX(t) + patchLocX(iX) - imgX) <= patchW/2;
			idy = abs(eyeY(t) + patchLocY(iY) - imgY) <= patchW/2;
			imLocal = img(idy,idx,:);
			[nY, nX, ~] = size(imLocal);
			sfX = (-nX/2 : nX/2-1) / nX * (60/pixPerArcmin);
			sfY = (-nY/2 : nY/2-1) / nY * (60/pixPerArcmin);
			sf2 = sfX.^2 + sfY'.^2;
			for(k = 1:3)
				% f = fftshift(fft2(imLocal(:,:,k)));
				% f(sf2 > cutSF^2) = f(sf2 > cutSF^2) * senHigh(iY,iX,t);
				% f(sf2 <= cutSF^2 & sf2 > 0) = f(sf2 <= cutSF^2 & sf2 > 0) * senLow(iY,iX,t);
				% hImg.CData(idy,idx,k) = real(ifft2(ifftshift(f)));
				hImg.CData(idy,idx,k) = imLow(idy,idx,k) * senLow(iY,iX,t) + imHigh(idy,idx,k) * senHigh(iY,iX,t);
			end
		end
	end
% 	mask = 
	drawnow;
	writeVideo(writerObj, getframe(gcf));
end
close(writerObj);