encoder = Encoder([], true);
% encoder = Encoder([], false);   % H_s not 1

sbj = 'SacDB'; 'A014'; %'A0NK';

inFolder = 'UG - Noise & Grating Simulated Separately';
folder = 'UG - Noise & Grating Simulated Separately';
% folder = 'UG - Noise & Grating Simulated Separately - Drift Stabilized';
% folder = 'UG - Noise & Grating Simulated Separately - Saccade Stabilized';
withInternalNoise = true;
% durOffset = 0;
durOffset = 57;

% gather data
encoder.LoadExampleCellsActivities(fullfile('../../Data/Simulated Activities', sbj, folder, ['All Conditions' '.mat']));
encoder.AddInternalNoise([], fullfile('../../Data/Simulated Activities', sbj, inFolder, ['All Conditions' '.mat']));