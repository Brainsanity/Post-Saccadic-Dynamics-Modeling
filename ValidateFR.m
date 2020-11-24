%%
isLoad = true; false;
if(isLoad)
	contrasts = [0    0.0100    0.0300    0.0500    0.0700    0.0900    0.1000    0.2000    0.3000    0.4000    0.5000];
	for(k = size(contrasts,2) : -1 : 2)
		load(sprintf('F:/Post Saccadic Dynamics Modeling/Data/figures/Real ratio cells, avContrast=0.5, contrast=%.2f, same 30 trials; Hs=1/Real ratio cells, avContrast=0.5, contrast=%.2f, same 30 trials; Hs=1.mat', contrasts(k), contrasts(k)));
		idx = find([conditions.eccentricity] == 0 & [conditions.sf] == 2);
		lfr(:,:,k) = [mean(LFR{idx(1),1}(1,:,:),3); mean(LFR{idx(2),1}(1,:,:),3)];		% 1st POn cell & 1st trial; 	1st row: absent;	2nd row: ecc=0 & sf=2
        NeuronIdx{k} = cellIdx{1,1};
	end
end

%%
figure; hold on;
for(k = 2 : size(contrasts,2))
	plot( time{1}, lfr(1,:,k)', 'linewidth', 2, 'displayname', sprintf('%.2f', contrasts(k)) );
end
legend;