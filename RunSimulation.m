function RunSimulation(iConds, nHours)
	if(~exist('nHours', 'var') || isempty(nHours))
		nHours = Inf;
	end
	addpath('../../EMAT-2.1.0.0');

	encoder = Encoder([], true);
	sbj = 'SacDB';
	folder = 'UG - Noise & Grating Simulated Separately';
	encoder.SimulateExampleCellsActivities_SacDB(iConds, '../../Data', [], 'none', folder, nHours);
end