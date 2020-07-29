classdef TemporalRF < handle

	properties (Static)
		paperParams;
	end

	methods
		
		function obj = TemporalRF()
			paperParams = containers.Map();
			paperParams('POnCenter') = struct( 	'A',		[13.98], ...	% spikes / (s * unit contrast); mean +- std
									            'D',		4.0/1000, ...	% s, delay of transmission from the optic chiasm to the LGN
									            'N_L',		27, ...
									            'Tau_L',	2.05/1000, ...  % s
									            'H_S',		0.68, ...
									            'N_H',		1, ...
									            'Tau_H',	21.25/1000 );	% s

			paperParams('POffCenter') = struct(	'A',		10.42, ...
									        	'D',		4.0/1000, ...
									        	'N_L',		37, ...
									        	'Tau_L',	1.54/1000, ...
									        	'H_S',		0.61, ...
									        	'N_H',		1, ...
									        	'Tau_H',	24.41/1000 );

			paperParams('MOn') = struct(	'A', 		566.92, ...
								            'D', 		2.2/1000, ...   % s
								            'N_L', 		30.30, ...
								            'Tau_L', 	1.41/1000, ...  % s
								            'H_S', 		0.98, ...        % 
								            'Tau_0', 	54.6/1000, ...  % s
								            'C_half', 	0.056 );

			paperParams('MOn') = struct(	'A',		550.14, ...
								            'D',		2.31/1000, ...   % s
								            'N_L',		22.60, ...
								            'Tau_L',	1.98/1000, ...  % s
								            'H_S',		0.93, ...        % 
								            'Tau_0',	153.34/1000, ...  % s
								            'C_half',	0.051 );

			obj.paperParams = paperParams;
		end

	end

end