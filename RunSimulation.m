function RunSimulation(iConds, nHours, stabilize)
	if(~exist('nHours', 'var') || isempty(nHours))
		nHours = Inf;
	end
	if( ~exist('stabilize', 'var') || isempty(stabilize) )
		stabilize = 'normal';
	end

	addpath('../../EMAT-2.1.0.0');

	encoder = Encoder([], true);
	sbj = 'SacDB';
	if(strcmpi(stabilize, 'drift'))
		folder = 'UG - Noise & Grating Simulated Separately - Drift Stabilized';
	elseif(strcmpi(stabilize, 'saccade'))
		folder = 'UG - Noise & Grating Simulated Separately - Saccade Stabilized';
	else
		folder = 'UG - Noise & Grating Simulated Separately';
	end
	encoder.SimulateExampleCellsActivities_SacDB(iConds, '../../Data', [], stabilize, folder, nHours);
end