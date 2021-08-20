if true

	if false
		encoder = Encoder([], true);
		sbj = 'SacDB';
		inFolder = 'UG - Noise & Grating Simulated Separately';
		folder = 'UG - Noise & Grating Simulated Separately';
		encoder.LoadExampleCellsActivities(fullfile('../../Data/Simulated Activities', sbj, folder, ['All Conditions' '.mat']));
	end

	%% Plot the profile of mean activity of a population of example cells
	withInternalNoise = false;
	alignEvent = 'saccadeOff';

	conditions = encoder.activityParams.conditions;

	SFs = unique([conditions.sf]);
	SFs(SFs == 0) = [];
	Eccs = unique([conditions.eccentricity]);

	figure( 'NumberTitle', 'off', 'name', 'Mean FR of Example Cells', 'color', 'w' );
	pause(0.1);
	jf = get(handle(gcf),'javaframe');
	jf.setMaximized(1);
	pause(1);
	colors = {'r', 'g', 'b'};
	for(k = 10:-1:1)
	    hAxes(k) = axes( 'Position', [0.06+0.96/5*(mod(k-1,5)), 0.08+0.92/2*(k<6), 0.96/5-0.03, 0.9/2-0.04], ...
	                     'xlim', [-150 550], 'xtick', -100:200:500, 'fontsize', 18, 'LineWidth', 2, 'nextplot', 'add' );
	    if(k > 5)
	        xlabel('Time (ms)');
	    else
	        set(gca, 'XTickLabel', []);
	        if(k < 5)
	        	title(encoder.layers(k).name);
	        else
	        	title('Layers Average');
	        end
	    end
	    if(k == 1 || k == 6)
	        ylabel({sprintf('%d cpd', SFs((k>1)+1)), 'Firing rate (spikes/s)'});
	    end

	    iSF = (k>5) + 1;
        iL = mod(k-1, 5) + 1;
	    h = [];
	    for(iEcc = size(Eccs,2) : -1 : 1)
	    	iCond = [conditions.eccentricity] == Eccs(iEcc) & [conditions.sf] == SFs(iSF);
	    	tMin = max( round( encoder.activityParams.timeline(1)/1000*[encoder.activityParams.trials.sRate] ) + [encoder.activityParams.trials.(conditions(iCond).alignEvent)] - [encoder.activityParams.trials.(alignEvent)] );
			tMax = min( round( encoder.activityParams.timeline(end)/1000*[encoder.activityParams.trials.sRate] ) + [encoder.activityParams.trials.(conditions(iCond).alignEvent)] - [encoder.activityParams.trials.(alignEvent)] );
			t = (tMin:tMax) / encoder.activityParams.trials(1).sRate * 1000;
			
			if(iL ~= 5)
				frBG = encoder.ExampleCellsActivitiesOnCondition(encoder.layers(iL).name, Eccs(iEcc), 0, 0, alignEvent, t([1 end]), false, withInternalNoise);
				fr = max(0, encoder.ExampleCellsActivitiesOnCondition(encoder.layers(iL).name, Eccs(iEcc), SFs(iSF), 0.5, alignEvent, t([1 end]), true, withInternalNoise));
			else
				frBG = cat(1, encoder.ExampleCellsActivitiesOnCondition(encoder.layers(1).name, Eccs(iEcc), 0, 0, alignEvent, t([1 end]), false, withInternalNoise), ...
							  encoder.ExampleCellsActivitiesOnCondition(encoder.layers(2).name, Eccs(iEcc), 0, 0, alignEvent, t([1 end]), false, withInternalNoise), ...
							  encoder.ExampleCellsActivitiesOnCondition(encoder.layers(3).name, Eccs(iEcc), 0, 0, alignEvent, t([1 end]), false, withInternalNoise), ...
							  encoder.ExampleCellsActivitiesOnCondition(encoder.layers(4).name, Eccs(iEcc), 0, 0, alignEvent, t([1 end]), false, withInternalNoise));
				fr = cat(1, max(0, encoder.ExampleCellsActivitiesOnCondition(encoder.layers(1).name, Eccs(iEcc), SFs(iSF), 0.5, alignEvent, t([1 end]), true, withInternalNoise)), ...
						    max(0, encoder.ExampleCellsActivitiesOnCondition(encoder.layers(2).name, Eccs(iEcc), SFs(iSF), 0.5, alignEvent, t([1 end]), true, withInternalNoise)), ...
						    max(0, encoder.ExampleCellsActivitiesOnCondition(encoder.layers(3).name, Eccs(iEcc), SFs(iSF), 0.5, alignEvent, t([1 end]), true, withInternalNoise)), ...
						    max(0, encoder.ExampleCellsActivitiesOnCondition(encoder.layers(4).name, Eccs(iEcc), SFs(iSF), 0.5, alignEvent, t([1 end]), true, withInternalNoise)));
			end

			frBG = mean(frBG, 1);	% mean population activity
            mBG = mean(frBG, 3);
            semBG = std(frBG, [], 3) / sqrt(size(frBG,3));
            qua25BG = prctile(frBG, 25, 3);
            qua75BG = prctile(frBG, 75, 3);

            fr = mean(fr, 1);		% mean population activity
			m = mean(fr, 3);
			sem = std(fr, [], 3) / sqrt(size(fr,3));
            qua25 = prctile(fr, 25, 3);
            qua75 = prctile(fr, 75, 3);

            if(any(iEcc == 1:3))
	            fill( [t, t(end:-1:1)], [m-sem, m(end:-1:1)+sem(end:-1:1)], colors{iEcc}, 'LineStyle', 'none', 'FaceAlpha', 0.2 );
				h(iEcc) = plot( t, m, '-', 'color', colors{iEcc}, 'lineWidth', 2, 'displayName', sprintf('Grating ecc=%d', Eccs(iEcc)) );

	            % fill( [t, t(end:-1:1)], [mBG-semBG, mBG(end:-1:1)+semBG(end:-1:1)], colors{iEcc}, 'LineStyle', 'none', 'FaceAlpha', 0.2 );
				% h(size(Eccs,2)+1) = plot( t, mBG, '-.', 'color', colors{iEcc}, 'lineWidth', 2, 'displayName', 'Natural noise' );
			
            	drawnow;
            end
        end
        h(3+1) = plot( [0 0], [0, max(get(gca, 'ylim'))], 'k--', 'lineWidth', 2, 'displayName', 'Saccade off' );
        h(3+2) = plot( [1 1] * mean(([encoder.activityParams.trials.saccadeOn] - [encoder.activityParams.trials.(alignEvent)]) ./ [encoder.activityParams.trials.sRate]*1000), get(gca, 'ylim'), '--', 'color', [0.5 0.5 0.5], 'lineWidth', 2, 'displayName', 'Mean SacOn' );
        drawnow;
        
        if(k == 5)
            legend(h, 'location', 'northeast');
        end

        set(gca, 'ylim', [0 max(get(gca, 'ylim'))]);
	end

	% save mean population firing rates
	tSacOn = mean(([encoder.activityParams.trials.saccadeOn] - [encoder.activityParams.trials.(alignEvent)]) ./ [encoder.activityParams.trials.sRate]*1000);
	t(t < tSacOn) = [];
	for(k = 2:-1:1)		% P cells, M cells
		for(iSF = 2:-1:1)	% 2 cpd, 10 cpd
			for(iEcc = size(Eccs,2) : -1 : 1)
				fr = cat(1, max(0, encoder.ExampleCellsActivitiesOnCondition(encoder.layers((k-1)*2+1).name, Eccs(iEcc), SFs(iSF), 0.5, alignEvent, t([1 end]), true, withInternalNoise)), ...
						    max(0, encoder.ExampleCellsActivitiesOnCondition(encoder.layers((k-1)*2+2).name, Eccs(iEcc), SFs(iSF), 0.5, alignEvent, t([1 end]), true, withInternalNoise)));
				popFR{k}{iSF} = mean(mean(fr, 1), 3);
			end
		end
	end
	tPopFR = t;
	save('../../Data/Simulated Activities/SacDB/UG - Noise & Grating Simulated Separately/figures - withInternalNoise/uni-no_bias-fa25_hit75 - durOffset=0/PopFR.mat', 'popFR', 'tPopFR');
end

load('../../Data/Simulated Activities/SacDB/UG - Noise & Grating Simulated Separately/figures - withInternalNoise/uni-no_bias-fa25_hit75 - durOffset=0/PopFR.mat', 'popFR', 'tPopFR');