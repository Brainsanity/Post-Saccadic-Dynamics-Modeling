load(['../../Data/Mosaic_POn_Radius20.0deg_maxMovPrctile20.mat']);
rfPositions = rfPositions( abs(rfPositions(:,1))<=2 & abs(rfPositions(:,2))<=2, : );

ck = CronerKaplanRGCModel( 'dataSetToFit', 'raw', 'fitIntercept', false, 'randomizeCenterRadii', true, 'randomizeCenterSensitivities', true, 'randomizeSurroundRadii', true, 'randomizeSurroundSensitivities', true );
PRFParams = ck.SynthesizeRFParams(sqrt(sum(rfPositions.^2,2))','POn');
POFffRFParams = ck.SynthesizeRFParams(sqrt(sum(rfPositions.^2,2))','POff');

%%
iRF = 100;

f = 0.001:0.001:20;
fr_c = PRFParams(iRF).centerPeakSensitivities * PRFParams(iRF).centerRadii^2 * exp(-(pi*PRFParams(iRF).centerRadii*f).^2);
fr_s = - PRFParams(iRF).surroundPeakSensitivities * PRFParams(iRF).surroundRadii^2 * exp(-(pi*PRFParams(iRF).surroundRadii*f).^2);
[~, k] = max(fr_c + fr_s);

img = ToolKit.Gabor( 1/f(k)*60,0, 0, 601, 601 );

%%
sFR = []; sFR_c = []; sFR_s = [];
for( k = 100:-1:1 )
	if(~mod(k-1,10))
		fprintf('%d/100\n',k);
	end
	% [sFR(k,:), sFR_c(k,:), sFR_s(k,:)] = ck.LinearResponse( img*(k/100)/2+0.5, (-300:300)/60, (-300:300)/60, (1:300)/1000*2, zeros(1,300), PRFParams(iRF), 0, 0);
	% [sFR(k,:), sFR_c(k,:), sFR_s(k,:)] = ck.LinearResponse( img*(k/100), (-300:300)/60, (-300:300)/60, (1:300)/1000*2, zeros(1,300), PRFParams(iRF), 0, 0);
	tModulate = [-ones(1,300), ones(1, round(1000/1.06/2)), -ones(1, round(1000/1.06/2))] * k/100;
	for( t = 1 : 1200 )
		[sFR(k,t), sFR_c(k,t), sFR_s(k,t)] = ck.LinearResponse( img*tModulate(t)/2, (-300:300)/60, (-300:300)/60, 0, 0, PRFParams(iRF), 0, 0);
	end
end

%%
tFR = []; tRF = [];
for( k = 100:-1:1 )
	if(~mod(k-1,10))
		fprintf('%d/100\n',k);
	end
	[tFR(k,:), tRF(k,:)] = RGC.TemporalLinearModel( sFR(k,:), 'm', 'on', k/100 );
    tFR(k,:) = tFR(k,:)*k/100;
%     tRF(k,:) = tRF(k,:)*k/100;
end

%%
subplot(2,2,1); hold on;
for( k = 1:100 )
	if(~mod(k-1,10))
		fprintf('%d/100\n',k);
	end
	plot( sFR(k,:) );
	pause
end

%%
subplot(2,2,2); hold on;
for( k = 1:100 )
	if(~mod(k-1,10))
		fprintf('%d/100\n',k);
	end
	plot( tFR(k,:) );
	pause
end

%%
subplot(2,2,3); hold on;
for( k = 1:100 )
	if(~mod(k-1,10))
		fprintf('%d/100\n',k);
	end
	plot( tRF(k,:) );
	pause
end

%% temporal frequency sensitivity profile
subplot(2,2,4); hold on;
tFreqs = 0.1 : 0.1 : 32;
for(c = 0.01:0.01:1)
	plot( log(tFreqs)/log(2), log( c*RGC.TemporalFreqGainProfile( tFreqs, 'm', 'on', c ) ) / log(2), 'displayname', num2str(c) );
    pause;
end
legend