classdef TemporalRF < handle

	properties (SetAccess = private)
		paperParams;
		cutOff = 2;			% interpolate frequencies below this, so that sensitivity is 0 when frequency is 0 Hz
	end

	methods
		
		function obj = TemporalRF()
			%% Constructor function: loading parameters from papers

			paperParams = containers.Map();
			

			% Benardete & Kaplan, Visual Neuroscience, 1997, I
			%	min, max, median, mean, s.d., coefficient of variation (std divided by the mean)
			paperParams('POnCenter97')	= struct(	'A',		[19.41,		128.10,		67.59,		63.68,		33.06,		0.52], ...	% spikes / (s * unit contrast); mean +- std
										            'D',		[2.5,		4,			3.5,		3,47,		0.51,		0.15] / 1000, ...	% s, delay of transmission from the optic chiasm to the LGN
										            'N_L',		[20,		46,			38,			36.38,		6.63,		0.18], ...
										            'Tau_L',	[39.64/20,	55.06/46,	48.15/38,	48.05/36.38,4.181/6.63,	0.09/0.18] / 1000, ...  % s
										            'H_S',		[0.57		0.86,		0.69,		0.69,		0.078,		0.11], ...
										            'Tau_H',	[11.9,		52.66,		29.36,		31.24,		13.93,		0.45] / 1000 );	% s

			paperParams('POnSurround97')	= struct( 	'A',		[29.67,		519.63,		49.98,		126.46,		193.23,		1.53], ...	% spikes / (s * unit contrast); mean +- std
										            'D',		[2.5,		4,			3.5,		3,47,		0.51,		0.15] / 1000, ...	% s, delay of transmission from the optic chiasm to the LGN
										            'N_L',		[98,		133,		111,		111.5,		13.19,		0.12], ...
										            'Tau_L',	[47.79/98,	63.99/133,	55.02/111,	55.84/111.5,5.37/13.19,	0.10/0.12] / 1000, ...  % s
										            'H_S',		[0,			0.92,		0.48,		0.43,		0.32,		0.74], ...
										            'Tau_H',	[0,			38.71,		18.62,		18.93,		20.39,		1.08] / 1000 );	% s

			paperParams('POffCenter97') 	= struct(	'A',		[29.58,		114.05,		52.95,		68.96,		25.31,		0.37], ...
										        	'D',		[2.5,		4,			3.5,		3.5,		0.61,		0.17] / 1000, ...
										        	'N_L',		[21,		41,			27,			27.64,		5.48,		0.20], ...
										        	'Tau_L',	[42.39/21,	55.40/41,	50.38/27,	50.24/27.64,4.32/5.48,	0.09/0.20] / 1000, ...
										        	'H_S',		[0.50,		0.82,		0.75,		0.72,		0.10,		0.14], ...
										        	'Tau_H',	[20.85,		101.77,		35.68,		45.35,		28.29,		0.62] / 1000 );

			paperParams('POffSurround97')	= struct( 	'A',		[12.54,		420.20,		23.91,		76.07,		139.97,		1.84], ...	% spikes / (s * unit contrast); mean +- std
										            'D',		[2.5,		4,			3.5,		3.5,		0.61,		0.17] / 1000, ...	% s, delay of transmission from the optic chiasm to the LGN
										            'N_L',		[45,		168,		67,			78.88,		41.46,		0.53], ...
										            'Tau_L',	[53.00/45,	64.20/168,	60.33/67,	59.89/78.88,3.61/41.46,	0.06/0.53] / 1000, ...  % s
										            'H_S',		[0,			0.65,		0.42,		0.38,		0.22,		0.59], ...
										            'Tau_H',	[0,			14.05,		52.00,		60.15,		42.79,		0.71] / 1000 );	% s

			paperParams('nPOnCenter97') = 13;
			paperParams('nPOffCenter97') = 11;
			paperParams('nPOnSurround97') = 6;
			paperParams('nPOffSurround97') = 8;



			% Benardete & Kaplan, JOP, 1999		(values for D taken from Benardete & Kaplan, Visual Neuroscience, 1997, I)
			%	min, max, median, mean, s.d.
			paperParams('POnCenter99')	= struct(	'A',		[7.48,		111.42,		13.7,		22.87,		27.11], ...	% spikes / (s * unit contrast); mean +- std
										            'D',		[2.5,		4,			3.5,		3,47,		0.51] / 1000, ...	% s, delay of transmission from the optic chiasm to the LGN
										            'N_L',		[16,		36,			26,			26.92,		4.65], ...
										            'Tau_L',	[44.28/16,	64.53/36,	55.33/26,	53.80/26.92,5.25/4.65] / 1000, ...  % s
										            'H_S',		[0.41,		0.94,		0.64,		0.67,		0.15], ...
										            'Tau_H',	[1.78,		51.85,		16.31,		21.10,		14.91] / 1000 );	% s

			paperParams('POnSurround99')	= struct( 	'A',		[2.41,		43.12,		10.44,		12.63,		9.69], ...	% spikes / (s * unit contrast); mean +- std
										            'D',		[2.5,		4,			3.5,		3,47,		0.51] / 1000, ...	% s, delay of transmission from the optic chiasm to the LGN
										            'N_L',		[27,		101,		39,			46.15,		20.64], ...
										            'Tau_L',	[46.82/27,	65.31/101,	55.88/39,	56.35/46.15,5.38/20.64] / 1000, ...  % s
										            'H_S',		[0.48,		0.87,		0.57,		0.62,		0.11], ...
										            'Tau_H',	[2.04,		40.92,		29.59,		25.88,		10.78] / 1000 );	% s

			paperParams('POffCenter99') 	= struct(	'A',		[10.15,		20.29,		14.20,		14.27,		3.34], ...
										        	'D',		[2.5,		4,			3.5,		3.5,		0.61] / 1000, ...
										        	'N_L',		[19,		36,			27,			26.71,		5.21], ...
										        	'Tau_L',	[45.40/19,	60.76/36,	55.28/27,	52.76/26.71,6.19/5.21] / 1000, ...
										        	'H_S',		[0.56,		0.79,		0.74,		0.71,		0.09], ...
										        	'Tau_H',	[14.60,		30.86,		26.65,		24.22,		5.63] / 1000 );

			paperParams('POffSurround99')	= struct( 	'A',		[8.27,		65.13,		10.11,		18.04,		20.85], ...	% spikes / (s * unit contrast); mean +- std
										            'D',		[2.5,		4,			3.5,		3.5,		0.61] / 1000, ...	% s, delay of transmission from the optic chiasm to the LGN
										            'N_L',		[20,		45,			32,			32,			8.52], ...
										            'Tau_L',	[48.66/20,	64.71/45,	58.70/32,	57.78/32,	6.65/8.52] / 1000, ...  % s
										            'H_S',		[0.51,		0.90,		0.72,		0.72,		0.15], ...
										            'Tau_H',	[2.94,		38.21,		16.33,		20.16,		14.07] / 1000 );	% s

			paperParams('nPOnCenter99') = 13;
			paperParams('nPOffCenter99') = 7;
			paperParams('nPOnSurround99') = 13;
			paperParams('nPOffSurround99') = 7;




			% Benardete & Kaplan, Visual Neuroscience, 1999
			%	min, max, median, mean, SEM
			paperParams('MOn') = struct(	'A', 		[249.03,	1137.04,	499.77,		566.92, 	89.25], ...
								            'D', 		[2,			2.5,		2,			2.2,		0.095] / 1000, ...   % s
								            'N_L', 		[21,		41,			30,			30.30,		1.86], ...
								            'Tau_L', 	[0.91,		2.03,		1.44,		1.41,		0.11] / 1000, ...  % s
								            'H_S', 		[0.91,		1.0,		1.0,		0.98, 		0.01], ...        % 
								            'Tau_0', 	[27.89,		213.36,		37.34,		54.60,		17.76] / 1000, ...  % s
								            'C_half', 	[0.015,		0.087,		0.051,		0.056,		0.026] );

			paperParams('MOff') = struct(	'A',		[149.38,	856.43,		517.40,		550.14, 	70.92], ...
								            'D',		[2,			2.5,		2.4,		2.31,		0.086] / 1000, ...   % s
								            'N_L',		[18,		31,			22,			22.60,		1.19], ...
								            'Tau_L',	[1.39,		2.47,		2.05,		1.98,		0.09] / 1000, ...  % s
								            'H_S',		[0.80,		1.00,		0.97,		0.93,		0.03], ...        % 
								            'Tau_0',	[12.04,		611.42,		57.62,		153.34,		60.41] / 1000, ...  % s
								            'C_half',	[0.018,		0.116,		0.039,		0.051,		0.035] );

			paperParams('nMOnCell') = 10;
			paperParams('nMOffCell') = 10;


			obj.paperParams = paperParams;
		end


		function params = SynthesizeParams( obj, cellType, nCells, isMedian, whichPaper, HsBe1 )
			%% Generate paparameters according to paper reports
			%   isMedian:		true for median values; false for mean values
			%   whichPaper:		used for P cells: '97' refers to Benardete & Kaplan, Visual Neuroscience, 1997, I; '99' refers to Benardete & Kaplan, JOP, 1999

			if( nargin() < 4 || isempty(isMedian) )
				isMedian = true;
			end
			if( nargin() < 5 || isempty(whichPaper) )
				whichPaper = '97';
			end
			if( nargin() < 6 || isempty(HsBe1) )
				HsBe1 = false;
			end

			cellType = [upper(cellType(1:2)), lower(cellType(3:end))];

			params = [];

			switch cellType
				case {'POn', 'POff'}

					params(nCells).centerA = [];
					[params.centerA] = deal(obj.paperParams([cellType, 'Center', whichPaper ]).A(3+(~isMedian)));
					[params.centerD] = deal(obj.paperParams([cellType, 'Center', whichPaper ]).D(3+(~isMedian)));
					[params.centerN_L] = deal(obj.paperParams([cellType, 'Center', whichPaper ]).N_L(3+(~isMedian)));
					[params.centerTau_L] = deal(obj.paperParams([cellType, 'Center', whichPaper ]).Tau_L(3+(~isMedian)));
					if(HsBe1)
						[params.centerH_S] = deal(1);
					else
						[params.centerH_S] = deal(obj.paperParams([cellType, 'Center', whichPaper ]).H_S(3+(~isMedian)));
					end
					[params.centerTau_H] = deal(obj.paperParams([cellType, 'Center', whichPaper ]).Tau_H(3+(~isMedian)));

					[params.surroundA] = deal(obj.paperParams([cellType, 'Surround', whichPaper ]).A(3+(~isMedian)));
					[params.surroundD] = deal(obj.paperParams([cellType, 'Surround', whichPaper ]).D(3+(~isMedian)));
					[params.surroundN_L] = deal(obj.paperParams([cellType, 'Surround', whichPaper ]).N_L(3+(~isMedian)));
					[params.surroundTau_L] = deal(obj.paperParams([cellType, 'Surround', whichPaper ]).Tau_L(3+(~isMedian)));
					if(HsBe1)
						[params.surroundH_S] = deal(1);
					else
						[params.surroundH_S] = deal(obj.paperParams([cellType, 'Surround', whichPaper ]).H_S(3+(~isMedian)));
					end
					[params.surroundTau_H] = deal(obj.paperParams([cellType, 'Surround', whichPaper ]).Tau_H(3+(~isMedian)));
				
				case {'MOn', 'MOff'}
					params(nCells).A = [];
					[params.A] = deal(obj.paperParams(cellType).A(3+(~isMedian)));
					[params.D] = deal(obj.paperParams(cellType).D(3+(~isMedian)));
					[params.N_L] = deal(obj.paperParams(cellType).N_L(3+(~isMedian)));
					[params.Tau_L] = deal(obj.paperParams(cellType).Tau_L(3+(~isMedian)));
					if(HsBe1)
						[params.H_S] = deal(1);
					else
						[params.H_S] = deal(obj.paperParams(cellType).H_S(3+(~isMedian)));
					end
					[params.Tau_0] = deal(obj.paperParams(cellType).Tau_0(3+(~isMedian)));
					[params.C_half] = deal(obj.paperParams(cellType).C_half(3+(~isMedian)));

				otherwise
					fprintf('Cell type %s not recognized!!!\n', cellType);

			end
			
		end


		function [fr, fr_c, fr_s] = LinearResponse( obj, cellType, rfParams, contrasts, samRate, input1, input2 )
			%% Compute linear response with specified temporal parameters
			%	rfParams:		Structure array of temporal RF parameters generated with obj.SynthesizeParams(cellType, nCells)
			%	contrasts:		Contrast for each neuron
			%	samRate:		Sampling rate of input signals
			%	input1:			Input signal, which is linear spatial responses obtained with CronerKaplanRGCModel.LinearResponse(obj, stimulus, inputX, inputY, eyeX, eyeY, rfParams, rfX, rfY);
			%					each row for one neuron;
			%					center responses if P cell
			%	input2:			Input signal for the surround if P cell
			%	
			%	fr: 	        Linear response; 1st dimension for neurons, 2nd dimension for time
		    %	fr_c:        	Center linear response if P Cell
    		%	fr_s:        	Surround linear response if P Cell

    		cellType = [upper(cellType(1:2)), lower(cellType(3:end))];

    		contrasts = contrasts(:)';

    		switch cellType
				case {'POn', 'POff'}
					fr_c = linearFR( input1, [rfParams.centerA],	[rfParams.centerD],		[rfParams.centerN_L],	[rfParams.centerTau_L],		[rfParams.centerH_S], 	[rfParams.centerTau_H] );
					fr_s = linearFR( input2, [rfParams.surroundA],	[rfParams.surroundD],	[rfParams.surroundN_L],	[rfParams.surroundTau_L], 	[rfParams.surroundH_S],	[rfParams.surroundTau_H] );
					fr = fr_c + fr_s;

				case {'MOn', 'MOff'}
					Tau_H = [rfParams.Tau_0] ./ ( 1 + (contrasts./[rfParams.C_half]).^2 );	% s
					fr = linearFR( input1, [rfParams.A], [rfParams.D], [rfParams.N_L], [rfParams.Tau_L], [rfParams.H_S], Tau_H );
					fr_c = [];
					fr_s = [];

				otherwise
					fprintf('Cell type %s not recognized!!!\n', cellType);

			end

			function FR = linearFR( input, A, D, N_L, Tau_L, H_S, Tau_H )
				% cutOff = 0.5;	% interpolate frequencies below this, so that sensitivity is 0 when frequency is 0 Hz

				if(true)% size(input,2) < samRate )	% apply convolution in time domain if too few samples, apply multiplication in frequency domain otherwise
		            nSamples = samRate;	% number of samples for 1 s
		            isOdd = 0;
		        else
		        	nSamples = floor( size(input,2) / 2 ) * 2;
		        	isOdd = size(input,2) - nSamples;
		        end
	            t = (1:nSamples) / samRate;
	            f = (0:nSamples/2) / (nSamples/samRate);
	            w = 2*pi * f( 1 : nSamples/2+1 );
	            K = A' .* exp( -i*D'*w ) .* ( 1 - H_S' ./ (1 + i*Tau_H'*w) ) .* ( 1 ./ (1 + i*Tau_L'*w) ).^(N_L');
	            % K(:,1) = 0;
	            % idx = f==0|f>=obj.cutOff;
	            % K = interp1( f(idx), K(:,idx), f, 'pchip' );
	            K(:, end+1:end*2-2) = conj( K(:, end-1:-1:2) );

	            if(false)% nSamples == size(input,2) - isOdd )
	            	s = fft(input(:,1:nSamples), [], 2);		% along 2nd dimension
	            	FR = real( ifft(s.*K, [], 2) );				%%%%%% this method seems to produce error at the initial period
	            	if(isOdd)
	            		FR = [FR, FR(:,end)];
	            	end
	            else
		        	tRF = real( ifft(K, [], 2) );		% along 2nd dimension
		        	tRF = tRF(:, t < 0.5);
                    FR = zeros( size(tRF,1), size(input,2) );
                    parfor( iCell = 1 : size(tRF,1) )
                        FR(iCell,:) = filter( tRF(iCell,:), 1, input(iCell,:), [], 2 );		% along 2nd dimension
                    end
						% FR = FR(:,1:size(input,2));
		        end

		        FR = FR .* contrasts';
			end

		end

		function [sensitivity, sensitivity_s] = TemporalSensitivity(obj, cellType, rfParams, tf, contrasts)
			% cutOff = 0.5;	% interpolate frequencies below this, so that sensitivity is 0 when frequency is 0 Hz
			
			tf = tf(:);
			w = 2*pi * tf(:)';
			
			switch cellType
				case {'POn', 'POff'}
					sensitivity = abs( [rfParams.centerA]' .* exp( -i * [rfParams.centerD]' * w ) .* ( 1 - [rfParams.centerH_S]' ./ (1 + i * [rfParams.centerTau_H]' * w) ) .* ( 1 ./ (1 + i * [rfParams.centerTau_L]' * w) ).^([rfParams.centerN_L]') );
					sensitivity_s = abs( [rfParams.surroundA]' .* exp( -i * [rfParams.surroundD]' * w ) .* ( 1 - [rfParams.surroundH_S]' ./ (1 + i * [rfParams.surroundTau_H]' * w) ) .* ( 1 ./ (1 + i * [rfParams.surroundTau_L]' * w) ).^([rfParams.surroundN_L]') );
					
				case {'MOn', 'MOff'}
					Tau_H = [rfParams.Tau_0] ./ ( 1 + (contrasts./[rfParams.C_half]).^2 );	% s
					sensitivity = abs( [rfParams.A]' .* exp( -i * [rfParams.D]' * w ) .* ( 1 - [rfParams.H_S]' ./ (1 + i * Tau_H' * w) ) .* ( 1 ./ (1 + i * [rfParams.Tau_L]' * w) ).^([rfParams.N_L]') );
					
				otherwise
					fprintf('Cell type %s not recognized!!!\n', cellType);

			end
		end

	end

end