classdef RGC < handle
    methods (Static)
    	function [FRs, tFRs, movs] = main( version, modulator )
    		set( figure, 'NumberTitle', 'off', 'name', ['Firing Rate | ' modulator] ); hold on;

    		swPix = 1366;
    		swMm = 600;
    		shPix = 768;
    		shMm = 335;
    		sDist = 1620;

    		% grating simulus
    		contrast = 0.03;
    		sf = 3;
    		orientation = 0;
    		phase = 0;
    		wlPix = swPix ./ ( atand(swMm/2/sDist) * 2 * sf );
    		img = ToolKit.Gabor( wlPix, orientation, phase, swPix, shPix, 'grating' )' * contrast;
            % img = ToolKit.Gabor( wlPix, orientation, phase, swPix, shPix, 'grating' )' * 1 + 128;   % in intensity
    		
    		% stimulus gain as function of time: ramp + plateau
            tRamp = 1024 + 512; % ms
            tPlateau = 2048;    % ms
    		T = tRamp + tPlateau;	% ms
    		Gt = ones(1,T);
    		Gt(1:tRamp) = (0:tRamp-1)/(tRamp-1);

    		load('Trials.mat');
    		Trials(51) = [];
            Trials = Trials(1:5);
    		nValidTrials = size(Trials,2);
    		iTrial = 1;

    		% spatial receptive field
    		rfWDeg = 5;%10;
    		rfHDeg = 5;%10;
    		rfWPix = ceil( tand(rfWDeg/2) * 2 * sDist / swMm * swPix );
    		rfHPix = ceil( tand(rfHDeg/2) * 2 * sDist / shMm * shPix );
    		precision = rfWDeg / rfWPix;
    		rfLoc = round( ( [size(img)] - [rfWPix, rfHPix] ) / 2 ) - 1;	% bottom left point [x,y]
    		rfImg = zeros(rfHPix, rfWPix, T);

    		FRs = ones(size(Trials,2),T) * NaN;
            tFRs = ones(size(Trials,2),T) * NaN;

    		for( iTrial = 1 : size(Trials,2) )
	    		%% simulate blink as covering the whole visual field briefly (150 ms) 
	    		% get eye drift
				x = zeros( 1, sum([Trials(iTrial).drifts.duration]) );
				y = zeros(size(x));
				index = 0;
				Trials(iTrial).blinks.start = [Trials(iTrial).blinks.start, Trials(iTrial).notracks.start];
				Trials(iTrial).blinks.duration = [Trials(iTrial).blinks.duration, Trials(iTrial).notracks.duration];
				for( i = 1 : size( Trials(iTrial).drifts.start, 2 ) )
					idx = find( Trials(iTrial).drifts.start(i) == Trials(iTrial).blinks.start + Trials(iTrial).blinks.duration );
					if( ~isempty(idx) )
						st = Trials(iTrial).drifts.start(i) + 250;
						dur = Trials(iTrial).drifts.duration(i) - 250;
						if( dur < 1 )
							continue;
						else
							Trials(iTrial).drifts.start(i) = st;
							Trials(iTrial).drifts.duration(i) = dur;
						end
					end
					idx = find( Trials(iTrial).drifts.start(i) + Trials(iTrial).drifts.duration(i) == Trials(iTrial).blinks.start );
					if( ~isempty(idx) )
						dur = Trials(iTrial).drifts.duration(i) - 50;
						if( dur < 1 )
							continue;
						else
							Trials(iTrial).drifts.duration(i) = dur;
						end
					end
					x( (1 : Trials(iTrial).drifts.duration(i)) + index ) = Trials(iTrial).x.position( (0 : Trials(iTrial).drifts.duration(i)-1) + Trials(iTrial).drifts.start(i) );
					y( (1 : Trials(iTrial).drifts.duration(i)) + index ) = Trials(iTrial).y.position( (0 : Trials(iTrial).drifts.duration(i)-1) + Trials(iTrial).drifts.start(i) );
					if( index > 0 )
						x( (1 : Trials(iTrial).drifts.duration(i)) + index ) = x( (1 : Trials(iTrial).drifts.duration(i)) + index ) - x(index+1) + x(index);
						y( (1 : Trials(iTrial).drifts.duration(i)) + index ) = y( (1 : Trials(iTrial).drifts.duration(i)) + index ) - y(index+1) + y(index);
					end
					index = index + Trials(iTrial).drifts.duration(i);
				end
				x( index+1 : end ) = [];
				y( index+1 : end ) = [];
				if( size(x,2) < T )
					nValidTrials = nValidTrials - 1;
					continue;
				end
				x = x(1:T);
				y = y(1:T);

				% rotate by -orientation
				xx = x * cosd(orientation) + y * sind(orientation);
				yy = y * cosd(orientation) - x * sind(orientation);

				% center eye trace
				x = x - mean(x);
				y = y - mean(y);

				% convert to pixels
				x = round( sDist * tand(x/60) / swMm * swPix );
				y = round( sDist * tand(y/60) / shMm * shPix );

                % x(:) = 0;
                % y(:) = 0;

				if( strcmpi( version, 'covering' ) )

					% blink trace simulated as covering input image from top to bottom and then back to top
					[bx, by] = meshgrid( 1:swPix, 1:shPix );
					bx = bx';
					by = by';
					a = sind(-orientation);
					b = cosd(-orientation);
					dc = max( a * bx(:) - b * by(:) ) - min( a * bx(:) - b * by(:) ) + 1;
					bTrace = ones(1,T) * min( a * bx(:) - b * by(:) ) - 1;
					
					if( ~strcmp( modulator, 'drift' ) )
						tFullCover = 200;			% 150 ms of complete covering
						d = - round(tFullCover/2);% + randi(50) - 26;		% blink central time relative to stimulus center
						bTrace( (1:tFullCover) + (end-tFullCover)/2 + d ) = max( a * bx(:) - b * by(:) );
						dur = round(dc/swPix*20);	% move by PSX in 20 ms
						% bTrace( (1-dur:0) + (end-tFullCover)/2 + d ) = min( a * bx(:) - b * by(:) ) - 1 + round( (1:dur)/dur * dc );
						% bTrace( (dur:-1:1) + tFullCover + (end-tFullCover)/2 + d ) = min( a * bx(:) - b * by(:) ) - 1 + round( (1:dur)/dur * dc );
						
						% freeze eye trace during blink
						iStart = find( bTrace > min( a * bx(:) - b * by(:) ) - 1, 1, 'first' );	% start of blink
						iEnd = find( bTrace > min( a * bx(:) - b * by(:) ) - 1, 1, 'last' );		% end of blink
						if( ~isempty(iStart) && ~isempty(iEnd) )
							x(iStart:iEnd) = x(iStart-1);
							x(iEnd+1:end) = x(iEnd+1:end) - ( x(iEnd+1) - x(iStart-1) );
							y(iStart:iEnd) = y(iStart-1);
							y(iEnd+1:end) = y(iEnd+1:end) - ( y(iEnd+1) - y(iStart-1) );
						end
					end

					% create 3D input
					for( t = 1 : T )
						tmpMov = zeros(swPix,shPix);
						tmpMov( max(1,1-x(t)) : min(swPix,swPix-x(t)), max(1,1-y(t)) : min(shPix,shPix-y(t)) ) = img( max(1,1+x(t)) : min(swPix,swPix+x(t)), max(1,1+y(t)) : min(shPix,shPix+y(t)) ) * Gt(t) + 1;
						tmpMov( a * bx - b * by <= bTrace(t) ) = 0;
						rfImg(:,:,t) = tmpMov( (1:rfWPix)+rfLoc(1), (1:rfHPix)+rfLoc(2) )';
						% imshow( tmpMov', [-1 1]*contrast ); pause(0.1);
						% imshow( rfImg, [-1 1]*contrast ); pause(0.1);
					end
				end

                movs{iTrial} = rfImg;
                % continue;
	    		
	    		% Linearity
	    		% [fRate, fRate_C, fRate_S] = RGC.SpatioLinearModel( reshape( rfImg(:) * Gt, rfHPix, rfWPix, [] ), precision, round(rfWPix/2), round(rfHPix/2), 'p', 40 );
	    		[fRate, fRate_C, fRate_S] = RGC.SpatioLinearModel( rfImg, precision, round(rfWPix/2), round(rfHPix/2), 'm', 25 );
	    		% fRate = RGC.TemporalLinearModel( fRate_C, 'p', 'center' ) + RGC.TemporalLinearModel( fRate_S, 'p', 'surround' );
	    		fRate = RGC.TemporalLinearModel( fRate, 'm', 'on' );
	    		plot( 1:T, fRate, 'r' );

                tFRs(iTrial,:) = fRate;

	    		% Non-linearity: sigmoid
	    		K = 30;	% peak firing rate
	    		g = 1;	% gain, which modulates the slope of the S-curve
	    		thresh = 0;%60;%0.2;
	    		fRate = K ./ ( 1 + exp( -g * (fRate - thresh) ) );
	    		plot( 1:T, fRate, 'g' );

	    		FRs(iTrial,:) = fRate;

	    		% Noise

	    		
	    		% plot( 1:T, fRate, 'b' );

	    		% pause(0.5);
	    	end

            plot( 1:T, nanmean(tFRs,1), 'color', [0.5 0 0], 'LineWidth', 2 );
	    	plot( 1:T, nanmean(FRs,1), 'color', [0 0.5 0], 'LineWidth', 2 );


    	end

        function [FR, FR_C, FR_S] = SpatioLinearModel( input, precision, xCenter, yCenter, cellType, eccentricity, isInterp )
        	%% From Croner & Kaplan, Vision Research, 1995
            %% input:       2nd dimension as horizontal (x), 1st dimension as vertical (y); 3rd dimension as time
            %  precision:   degrees per pixel in input 
            %  xCenter:     horizontal location of the receptive field center on the input image, which is along the 2nd dimension of the input matrix
            %  yCenter:     vertical location of the receptive field cente on the input image, which is along the 1st dimension of the input matrx

            if( nargin() < 6 || eccentricity < 0 || eccentricity > 40 )
                index = 1;
                isInterp = false;   % interpolation for eccentricity infeasible in this case
            else
                index = 0;
            end
            if( nargin() < 7 )
                isInterp = false;   % no interpolation for eccentricity by default
            end

            if( lower(cellType) == 'p' )    % P cell
                ectRange = [ 0, 40; 0, 5; 5, 10; 10, 20; 20, 30; 30, 40 ];  % eccentricity ranges; degrees
                r_c = [ 0.05;   0.03;   0.05;   0.07;   0.09;   0.15 ]; % center radius, which gives 1/e of the peak sensitivity;   degrees
                K_c = [ 106.3;  325.2;  114.7;  77.8;   57.2;   18.6 ]; % center peak sensitivity;     spikes / ( s * %contrast * degree^2 )
                r_s = [ 0.42;   0.18;   0.43;   0.54;   0.73;   0.65 ]; % surround radius, which gives 1/e of the peak sensitivity;   degrees
                K_s = [ 1.1;    4.4;    0.7;    0.6;    0.8;    1.1 ];  % surround peak sensitivity;     spikes / ( s * %contrast * degree^2 )
            elseif( lower(cellType) == 'm' )    % M cell
                ectRange = [ 0, 40; 0, 10; 10, 20; 20, 30 ];  % eccentricity ranges; degrees
                r_c = [ 0.17;   0.10;   0.18;   0.23 ]; % center radius, which gives 1/e of the peak sensitivity;   degrees
                K_c = [ 84.7;   148.0;  115.0;  63.8 ]; % center peak sensitivity;     spikes / ( s * %contrast * degree^2 )
                r_s = [ 0.80;   0.72;   1.19;   0.58 ]; % surround radius, which gives 1/e of the peak sensitivity;   degrees
                K_s = [ 1.2;    1.1;    2.0;    1.6 ];  % surround peak sensitivity;     spikes / ( s * %contrast * degree^2 )
            else
                FR = [];
                FR_C = [];
                FR_S = [];
                return;
            end

            if( ~isInterp )
                if( index ~= 1 )
                    index = find( ectRange(1:end,1) <= eccentricity & eccentricity <= ectRange(1:end,2), 1, 'first' );
                    if( isempty(index) )
                        return;
                    else
                        index = index(1);
                    end
                end
                r_c = r_c(index);
                K_c = K_c(index);
                r_s = r_s(index);
                K_s = K_s(index);
            end

            x = ( (1 : size(input,2)) - xCenter ) * precision;
            y = ( (1 : size(input,1)) - yCenter ) * precision;
            [x,y] = meshgrid(x,y);
            x = x(:)';
            y = y(:)';
            input = reshape( input, [], size(input,3) );
            FR_C = K_c * exp( -(x.^2+y.^2) / r_c^2 ) * input * precision^2;
            FR_S = - K_s * exp( -(x.^2+y.^2) / r_s^2 ) * input * precision^2;
            FR = FR_C + FR_S;
        end


        function [FR, tRF ] = TemporalLinearModel( input, cellType, specifier, c )
        	%% From Benardete & Kaplan, JPhy, 1999 (P cells); and Benardete & Kaplan, VisNeuro, 1999 (M cells)
        	%  input:			firing rate resulted from spatial receptive field
        	%  cellType:		P cell or M cell
        	%  specifier:		ON cell or OFF cell when cellType is M; center or surround when cellType is P

           	FR = [];
           	tRF = [];
        	
        	if( lower(cellType) == 'p' )    % P cell
        		if( strcmpi( specifier, 'center' ) )
	                A = 13.98;     % spikes / (s * unit contrast)
		            D = 4.0/1000;	% s, delay of transmission from the optic chiasm to the LGN
		            N_L = 27;
		            Tau_L = 2.05/1000;  % s
		            H_S = 0.68;
		            N_H = 1;
		            Tau_H = 21.25/1000;	% s
		        elseif( strcmpi( specifier, 'surround' ) )
		        	A = 10.42;
		        	D = 4.0/1000;
		        	N_L = 37;
		        	Tau_L = 1.54/1000;
		        	H_S = 0.61;
		        	N_H = 1;
		        	Tau_H = 24.41/1000;
		        else
		        	return;
		        end

            elseif( lower(cellType) == 'm' )    % M cell
            	if( strcmpi( specifier, 'on' ) )
            		A = 566.92;
		            D = 2.2/1000;   % s
		            N_L = 30.30;
		            Tau_L = 1.41/1000;  % s
		            H_S = 0.98;        % H_S = 1 makes the summation of the kernel be zero
		            Tau_0 = 54.6/1000;  % s
		            C_half = 0.056;
		        elseif( strcmpi( specifier, 'off' ) )
            		A = 550.14;
		            D = 2.31/1000;   % s
		            N_L = 22.60;
		            Tau_L = 1.98/1000;  % s
		            H_S = 1;%0.93;        % H_S = 1 makes the summation of the kernel be zero
		            Tau_0 = 153.34/1000;  % s
		            C_half = 0.051;
            	else
            		return;
            	end
            	N_H = 1;
            	% c = 0.04;	% contrast
            	Tau_H = Tau_0 / ( 1 + (c/C_half)^2 );	% s
            
            else
            	return;
            end

        	sRate = 1000;   % sampling rate of 1000 Hz
            nSamples = 1000;%round( .350 * sRate );   % number of samples for 350 ms; even number
            t = (1:nSamples) / sRate;    % 350 ms
            f = (0:nSamples/2) / (nSamples/sRate);
            w = 2*pi * f( 1 : nSamples/2+1 );
            K = A * exp( -i*w*D ) .* ( 1 - H_S ./ (1 + i*w*Tau_H) ).^N_H .* ( 1 ./ (1 + i*w*Tau_L) ).^N_L;
            K(end+1:end*2-2) = conj(K(end-1:-1:2));
            tRF = real(ifft(K));

            if(isempty(input)) return; end

            FR = conv( input, tRF ) * (1/sRate);
            FR = FR(1:size(input,2));

        end


        function senseProfile = TemporalFreqGainProfile( tFreqs, cellType, specifier, c )
            %% From Benardete & Kaplan, JPhy, 1999 (P cells); and Benardete & Kaplan, VisNeuro, 1999 (M cells)
            %  tFreqs:          Temporal frequencies
            %  cellType:        P cell or M cell
            %  specifier:       ON cell or OFF cell when cellType is M; center or surround when cellType is P

            FR = [];
            tRF = [];
            
            if( lower(cellType) == 'p' )    % P cell
                if( strcmpi( specifier, 'center' ) )
                    A = 13.98;     % spikes / (s * unit contrast)
                    D = 4.0/1000;   % s, delay of transmission from the optic chiasm to the LGN
                    N_L = 27;
                    Tau_L = 2.05/1000;  % s
                    H_S = 0.68;
                    N_H = 1;
                    Tau_H = 21.25/1000; % s
                elseif( strcmpi( specifier, 'surround' ) )
                    A = 10.42;
                    D = 4.0/1000;
                    N_L = 37;
                    Tau_L = 1.54/1000;
                    H_S = 0.61;
                    N_H = 1;
                    Tau_H = 24.41/1000;
                else
                    return;
                end

            elseif( lower(cellType) == 'm' )    % M cell
                if( strcmpi( specifier, 'on' ) )
                    A = 566.92;
                    D = 2.2/1000;   % s
                    N_L = 30.30;
                    Tau_L = 1.41/1000;  % s
                    H_S = 1;%0.98;
                    Tau_0 = 54.6/1000;  % s
                    C_half = 0.056;
                elseif( strcmpi( specifier, 'off' ) )
                    A = 550.14;
                    D = 2.31/1000;   % s
                    N_L = 22.60;
                    Tau_L = 1.98/1000;  % s
                    H_S = 0.93;
                    Tau_0 = 153.34/1000;  % s
                    C_half = 0.051;
                else
                    return;
                end
                N_H = 1;
                % c = 0.4;    % contrast
                Tau_H = Tau_0 / ( 1 + (c/C_half)^2 );   % s
            
            else
                return;
            end

            w = 2*pi * tFreqs;
            senseProfile = abs( A * exp( -i*w*D ) .* ( 1 - H_S ./ (1 + i*w*Tau_H) ).^N_H .* ( 1 ./ (1 + i*w*Tau_L) ).^N_L );
        end


        function Test()
            t = 0:0.001:0.999;
            x = sind( 360*10*t );
            x(1:300) = 0;
            freqs = (0:500)/0.99;
            figure; hold on;
            h(1) = plot( t, x, '--', 'color', [0.5 0.5 0.5], 'DisplayName', 'x' );
            h(2) = plot( t, RGC.TemporalLinearModel(x,'m','on'), 'r', 'DisplayName', 'FR-t_domain' );

            
            K = RGC.TemporalFreqGainProfile(freqs,'m','on');
            K(end+1:end*2-2) = conj(K(end-1:-1:2));
            fProfile = (fft(x));
            % fProfile(2:500) = fProfile(2:500)*2;
            % fProfile = fProfile(1:501);
            FR = real(ifft(K.*fProfile));
            h(3) = plot( t, FR, 'g', 'DisplayName', 'FR-f_profile' );

            [~,tRF] = RGC.TemporalLinearModel([],'m','on');
            K = fft(tRF);
            FR = (ifft(K.*fProfile));
            sum(FR == real(FR))
            h(4) = plot( t, real(FR), 'b', 'DisplayName', 'FR-f_reverse' );

            legend(h);

            figure;
            x = [zeros(1,500),ones(1,500)];
            y = abs(fft(x));
            y(2:500) = y(2:500) * 2;
            y = y(1:501);
            plot(freqs,y);



        end

        
        function FR = Demo_CSAntagnism()
            [img, MAP ] = imread('Demo_CSAntagnism.png');
            OriImg = reshape( MAP(img+1,1), size(img) );
            
            figure;
            subplot(1,2,1);
            imshow(OriImg);
            
            OriImg = ( OriImg - mean(OriImg(:)) ) / mean(OriImg(:)) * 100;      % convert to contrast
            precision = 0.12;   % degrees/pixel
            distance = 0.12;    % distance between two adjacent RGCs (degrees)
            r = 1.8;              % radius of visual patch presented in RGC receptive field (degrees)
            FR = zeros( floor( size(img,1)/(distance/precision) ), floor( size(img,2)/(distance/precision) ) );
            
            for( i = 1 : size(FR,1) )
                for( j = 1 : size(FR,2) )
                    iCenter = (i-1) * distance / precision + 1;
                    jCenter = (j-1) * distance / precision + 1;
                    iLow = max( [1, round(iCenter - r/precision)] );
                    iHigh = min( [size(img,1), round(iCenter + r/precision)] );
                    jLow = max( [1, round(jCenter - r/precision)] );
                    jHigh = min( [size(img,2), round(jCenter + r/precision)] );
                    patch = OriImg( iLow : iHigh, jLow : jHigh );
                    FR(i,j) = SenseProfile( patch, precision, iCenter - iLow + 1, jCenter - jLow + 1 );
                end
            end
            
            subplot(1,2,2);
            imshow( FR + 1 - max(FR(:)) );
            
            
            function FR = SenseProfile( input, precision, iCenter, jCenter )
                K_c = 148;
                r_c = 0.1;
                K_s = 1.1*3;%1.1*1.7;
                r_s = 0.72;
                x = ( (1 : size(input,2)) - jCenter ) * precision;
                y = ( (1 : size(input,1)) - iCenter ) * precision;
                [x,y] = meshgrid(x,y);
                FR = sum( sum( ( K_c * exp( -(x.^2+y.^2) / r_c^2 ) - K_s * exp( -(x.^2+y.^2) / r_s^2 ) ) .* input * precision^2 ) );
            end
        end

        function Croner_Kaplan_1995()
            
            %% P cells
            ectRange = [ 0, 40; 0, 5; 5, 10; 10, 20; 20, 30; 30, 40 ];  % eccentricity ranges; degrees
            r_c = [ 0.05;   0.03;   0.05;   0.07;   0.09;   0.15 ]; % center radius, which gives 1/e of the peak sensitivity;   degrees
            K_c = [ 106.3;  325.2;  114.7;  77.8;   57.2;   18.6 ]; % center peak sensitivity;     spikes / ( s * %contrast * degree^2 )
            r_s = [ 0.42;   0.18;   0.43;   0.54;   0.73;   0.65 ]; % surround radius, which gives 1/e of the peak sensitivity;   degrees
            K_s = [ 1.1;    4.4;    0.7;    0.6;    0.8;    1.1 ];  % surround peak sensitivity;     spikes / ( s * %contrast * degree^2 )
            
            C = 50;  % contrast, in percent
            sf = 0.1:0.1:100;
            FontSize  = 20;

            % direct computation using Equation 1 in the paper
            rCenter = C * repmat(K_c,1,size(sf,2)) * pi .* repmat(r_c,1,size(sf,2)).^2 .* exp( -(pi * r_c * sf).^2 );    % center
            rSurround = C * repmat(K_s,1,size(sf,2)) * pi .* repmat(r_s,1,size(sf,2)).^2 .* exp( -(pi * r_s * sf).^2 );    % surround
            set( figure, 'NumberTitle', 'off', 'name', 'Croner&Kaplan.1995. Fig.1(a). Equation 1. P Cells' );
            for( i = 1 : size(r_c,1) )
                subplot(2,3,i); hold on;
                h = [];
                h(end+1) = plot( sf, rCenter(i,:), '--', 'LineWidth', 2, 'DisplayName', 'Center' );    % center
                h(end+1) = plot( sf, rSurround(i,:), '--', 'LineWidth', 2, 'DisplayName', 'Surround');    % surround
                h(end+1) = plot( sf, rCenter(i,:) - rSurround(i,:), 'LineWidth', 2, 'DisplayName', 'C - S' );    % center - surround
                title( sprintf( 'Eccentricity: [%d %d]', ectRange(i,:) ) );
                xlabel( 'Spatial frequency (cpd)' );
                ylabel( 'Firing rate (spikes/s)' );
                set( gca, 'XLim', [0.1 100], 'XTick', [0.1 1 10 100], 'XTickLabel', [0.1 1.0 10.0 100.0], 'YLim', [0.1 100], 'YTick', [0.1 1 10 100], 'YTickLabel', [0.1 1.0 10.0 100.0], 'LineWidth', 2, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log' );
            end
            legend(h);

            % convolution with gratings
            precision = 1/120/2;  % half half arc min
            w = round(8/precision);  % 4 degrees
            for( iSF = 1:size(sf,2) )
                img = C * repmat( cos( 2 * pi * sf(iSF) * ((1:w)-round(w/2)) * precision ), w, 1 );
                for( iEccentricity = 1:6 )
                    eccentricity = ectRange(iEccentricity,1);
                    if( iEccentricity == 1 )
                        eccentricity = -1;
                    end
                    [~, rCenter(iEccentricity,iSF), rSurround(iEccentricity,iSF)] = RGC.SpatioLinearModel( img, precision, round(w/2), round(w/2), 'P', eccentricity );
                end
            end
            set( figure, 'NumberTitle', 'off', 'name', 'Croner&Kaplan.1995. Fig.1(a). Convolution. P Cells' );
            for( i = 1 : size(r_c,1) )
                subplot(2,3,i); hold on;
                h = [];
                h(end+1) = plot( sf, rCenter(i,:), '--', 'LineWidth', 2, 'DisplayName', 'Center' );    % center
                h(end+1) = plot( sf, -rSurround(i,:), '--', 'LineWidth', 2, 'DisplayName', 'Surround');    % surround
                h(end+1) = plot( sf, rCenter(i,:) + rSurround(i,:), 'LineWidth', 2, 'DisplayName', 'C - S' );    % center - surround
                title( sprintf( 'Eccentricity: [%d %d]', ectRange(i,:) ) );
                xlabel( 'Spatial frequency (cpd)' );
                ylabel( 'Firing rate (spikes/s)' );
                set( gca, 'XLim', [0.1 100], 'XTick', [0.1 1 10 100], 'XTickLabel', [0.1 1.0 10.0 100.0], 'YLim', [0.1 100], 'YTick', [0.1 1 10 100], 'YTickLabel', [0.1 1.0 10.0 100.0], 'LineWidth', 2, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log' );
            end
            legend(h);


            %% M cells
            ectRange = [ 0, 40; 0, 10; 10, 20; 20, 30 ];  % eccentricity ranges; degrees
            r_c = [ 0.17;   0.10;   0.18;   0.23 ]; % center radius, which gives 1/e of the peak sensitivity;   degrees
            K_c = [ 84.7;   148.0;  115.0;  63.8 ]; % center peak sensitivity;     spikes / ( s * %contrast * degree^2 )
            r_s = [ 0.80;   0.72;   1.19;   0.58 ]; % surround radius, which gives 1/e of the peak sensitivity;   degrees
            K_s = [ 1.2;    1.1;    2.0;    1.6 ];  % surround peak sensitivity;     spikes / ( s * %contrast * degree^2 )
            
            C = 25;  % contrast, in percent
            sf = 0.1:0.1:100;
            FontSize  = 20;

            % direct computation using Equation 1 in the paper
            rCenter = C * repmat(K_c,1,size(sf,2)) * pi .* repmat(r_c,1,size(sf,2)).^2 .* exp( -(pi * r_c * sf).^2 );    % center
            rSurround = C * repmat(K_s,1,size(sf,2)) * pi .* repmat(r_s,1,size(sf,2)).^2 .* exp( -(pi * r_s * sf).^2 );    % surround
            set( figure, 'NumberTitle', 'off', 'name', 'Croner&Kaplan.1995. Fig.1(a). Equation 1. M Cells' );
            for( i = 1 : size(r_c,1) )
                subplot(2,2,i); hold on;
                h = [];
                h(end+1) = plot( sf, rCenter(i,:), '--', 'LineWidth', 2, 'DisplayName', 'Center' );    % center
                h(end+1) = plot( sf, rSurround(i,:), '--', 'LineWidth', 2, 'DisplayName', 'Surround');    % surround
                h(end+1) = plot( sf, rCenter(i,:) - rSurround(i,:), 'LineWidth', 2, 'DisplayName', 'C - S' );    % center - surround
                title( sprintf( 'Eccentricity: [%d %d]', ectRange(i,:) ) );
                xlabel( 'Spatial frequency (cpd)' );
                ylabel( 'Firing rate (spikes/s)' );
                set( gca, 'XLim', [0.1 100], 'XTick', [0.1 1 10 100], 'XTickLabel', [0.1 1.0 10.0 100.0], 'YLim', [0.1 300], 'YTick', [0.1 1 10 100], 'YTickLabel', [0.1 1.0 10.0 100.0], 'LineWidth', 2, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log' );
            end
            legend(h);

            % convolution with gratings
            precision = 1/120/2;  % half half arc min
            w = round(8/precision);  % 8 degrees
            for( iSF = 1:size(sf,2) )
                img = C * repmat( cos( 2 * pi * sf(iSF) * ((1:w)-round(w/2)) * precision ), w, 1 );
                for( iEccentricity = 1 : size(r_c,1) )
                    eccentricity = ectRange(iEccentricity,1);
                    if( iEccentricity == 1 )
                        eccentricity = -1;
                    end
                    [~, rCenter(iEccentricity,iSF), rSurround(iEccentricity,iSF)] = RGC.SpatioLinearModel( img, precision, round(w/2), round(w/2), 'M', eccentricity );
                end
            end
            set( figure, 'NumberTitle', 'off', 'name', 'Croner&Kaplan.1995. Fig.1(a). Convolution. M Cells' );
            for( i = 1 : size(r_c,1) )
                subplot(2,2,i); hold on;
                h = [];
                h(end+1) = plot( sf, rCenter(i,:), '--', 'LineWidth', 2, 'DisplayName', 'Center' );    % center
                h(end+1) = plot( sf, -rSurround(i,:), '--', 'LineWidth', 2, 'DisplayName', 'Surround');    % surround
                h(end+1) = plot( sf, rCenter(i,:) + rSurround(i,:), 'LineWidth', 2, 'DisplayName', 'C - S' );    % center - surround
                title( sprintf( 'Eccentricity: [%d %d]', ectRange(i,:) ) );
                xlabel( 'Spatial frequency (cpd)' );
                ylabel( 'Firing rate (spikes/s)' );
                set( gca, 'XLim', [0.1 100], 'XTick', [0.1 1 100], 'XTickLabel', [0.1 1.0 100.0], 'YLim', [0.1 300], 'YTick', [0.1 1 10 100], 'YTickLabel', [0.1 1.0 10.0 100.0], 'LineWidth', 2, 'FontSize', FontSize, 'XScale', 'log', 'YScale', 'log' );
            end
            legend(h);
        end


        function [K R] = Benardete_Kaplan_JPhy_1999(arg)
            A = 199.92;     % spikes / (s * unit contrast)
            N_L = 66;
            Tau_L = 0.74/1000;  % s
            H = 0.60;
            N_H = 1;
            Tau_H = 32.25/1000;
            D = 4.0/1000;

            % A = 345.24;
            % N_L = 29;
            % Tau_L = 1.65/1000;
            % H = 0.68;
            % N_H = 1;
            % Tau_H = 10.68/1000;
            % D = 4.0/1000;

            % A = 601.48;
            % N_L = 51;
            % Tau_L = 0.87/1000;
            % H = 0.77;
            % N_H = 1;
            % Tau_H = 31.73/1000;
            % D = 4.0/1000;

            % A = 102.28;
            % N_L = 121;
            % Tau_L = 0.33/1000;
            % H = 0.99;
            % N_H = 1;
            % Tau_H = 19.12/1000;
            % D = 4.0/1000;

            sRate = 1000;   % sampling rate of 1000 Hz
            nSamples = round( 1.350 * sRate );   % number of samples for 1350 ms; even number
            t = (1:nSamples) / sRate;    % 350 ms
            f = (0:nSamples/2) / (nSamples/sRate);
            w = 2*pi * f( 1 : nSamples/2+1 );
            K = A * exp( -i*w*D ) .* ( 1 - H ./ (1 + i*w*Tau_H) ).^N_H .* ( 1 ./ (1 + i*w*Tau_L) ).^N_L;
            % K = exp( real(K) ) + i * exp( real(-i*K) );
            figure;
            FontSize = 20;
            LineWidth = 2;
            subplot(2,2,1); hold on;
            plot( log2(f), log2(abs(K)), 'LineWidth', LineWidth );
            set( gca, 'XLim', [0 log2(64)], 'XTick', 0 : nextpow2(f(end)), 'XTickLabel', 2.^(0 : nextpow2(f(end))), 'YLim', [0 8], 'YTick', 0:2:8, 'YTickLabel', 2.^(0:2:8), 'LineWidth', 2, 'FontSize', FontSize );
            subplot(2,2,3); hold on;
            plot( log2(f), log2( angle(K) ), 'LineWidth', LineWidth );
            set( gca, 'XLim', [0 log2(64)], 'XTick', 0 : nextpow2(f(end)), 'XTickLabel', 2.^(0 : nextpow2(f(end))), 'LineWidth', 2, 'FontSize', FontSize );


            % K(2:end-1) = K(2:end-1) / 2;
            K(end+1:end*2-2) = conj(K(end-1:-1:2));
            R = ifft(K);

            subplot(2,2,[2 4]); hold on;
            plot( t([1 end]), [0 0], ':', 'LineWidth', LineWidth );
            plot( t, real(R), 'LineWidth', LineWidth );
            set( gca, 'XLim', [0 0.2], 'YLim', [-20 20], 'LineWidth', 2, 'FontSize', FontSize );

        end
        
        
        function Benardete_Kaplan_1999()
            
            tStep = 0.001;  % s
            f = 2 : 0.01 : 1/tStep;
            w = 2*pi*f;
            figure;
            FontSize = 20;
            LineWidth = 2;
            

            %% On cells
            A = 566.92;
            D = 2.2/1000;   % s
            N_L = 30.30;
            Tau_L = 1.41/1000;  % s
            H_S = 0.98;
            Tau_0 = 54.6/1000;  % s
            C_half = 0.056;
            
            count = 5;
            for( c = [0.01 0.02 0.04 0.08 0.12] )
                K = A * exp( -1i*w*D ) .* ( 1 - H_S ./ ( 1 + 1i*w * Tau_0 / ( 1 + (c/C_half)^2 ) ) ) ./ ( 1 + 1i*w*Tau_L ).^N_L;
                % contrast gain
                subplot(2,3,1); hold on;
                h(count) = plot( log2(f), log2(abs(K)), 'DisplayName', sprintf( '%.2f', c ), 'LineWidth', LineWidth );
                
                % phase
                subplot(2,3,2); hold on;
                angs = angle(K);
                index = find( angs(2:end) - angs(1:end-1) > 0 ) + 1;
                for( ii = index )
                    angs(ii:end) = angs(ii:end) - 2*pi;     % realign by -2*pi, so that the phase goes monotonically
                end
                plot( log2(f), angs/pi, 'color', get(h(count),'color'), 'LineWidth', LineWidth );

                count = count - 1;
            end
            subplot(2,3,1);
            set( gca, 'XLim', [0 6], 'XTick', 0:6, 'XTickLabel', 2.^(0:6), 'YLim', [3 10], 'YTick', 3:10, 'YTickLabel', 2.^(3:10), 'FontSize', FontSize, 'LineWidth', LineWidth );
            ylabel('Gain (spikes/second-u.c.)');
            subplot(2,3,2);
            set( gca, 'XLim', [0 6], 'XTick', 0:6, 'XTickLabel', 2.^(0:6), 'FontSize', FontSize, 'LineWidth', LineWidth );
            ylabel('Phase (\pi)');
            legend(h);

            subplot(2,3,3);
            c = 0:0.01:0.14;
            plot( c, Tau_0 ./ ( 1 + (c/C_half).^2 ), 'LineWidth', LineWidth );
            set( gca, 'LineWidth', LineWidth, 'FontSize', FontSize );
            ylabel('$\tau_s$ (ms)', 'interpreter', 'latex');

            
            %% Off cells
            A = 550.14;
            D = 2.31/1000;   % s
            N_L = 22.60;
            Tau_L = 1.98/1000;  % s
            H_S = 0.93;
            Tau_0 = 153.34/1000;  % s
            C_half = 0.051;

            count = 3;
            h = [];
            for( c = [0.03125 0.0625 0.125] )
                K = A * exp( -1i*w*D ) .* ( 1 - H_S ./ ( 1 + 1i*w * Tau_0 / ( 1 + (c/C_half)^2 ) ) ) ./ ( 1 + 1i*w*Tau_L ).^N_L;
                % contrast gain
                subplot(2,3,4); hold on;
                h(count) = plot( log2(f), log2(abs(K)), 'DisplayName', sprintf( '%.4f', c ), 'LineWidth', LineWidth );
                
                % phase
                subplot(2,3,5); hold on;
                angs = angle(K);
                index = find( angs(2:end) - angs(1:end-1) > 0 ) + 1;
                for( ii = index )
                    angs(ii:end) = angs(ii:end) - 2*pi;     % realign by -2*pi, so that the phase goes monotonically
                end
                plot( log2(f), angs/pi, 'color', get(h(count),'color'), 'LineWidth', LineWidth );

                count = count - 1;
            end
            subplot(2,3,4);
            set( gca, 'XLim', [0 6], 'XTick', 0:6, 'XTickLabel', 2.^(0:6), 'YLim', [3 10], 'YTick', 3:10, 'YTickLabel', 2.^(3:10), 'FontSize', FontSize, 'LineWidth', LineWidth );
            xlabel('Frequency (Hz)');
            ylabel('Gain (spikes/second-u.c.)');
            subplot(2,3,5);
            set( gca, 'XLim', [0 6], 'XTick', 0:6, 'XTickLabel', 2.^(0:6), 'FontSize', FontSize, 'LineWidth', LineWidth );
            xlabel('Frequency (Hz)');
            ylabel('Phase (\pi)');
            legend(h);
            
            subplot(2,3,6);
            c = 0:0.01:0.14;
            plot( c, Tau_0 ./ ( 1 + (c/C_half).^2 ) * 1000, 'LineWidth', LineWidth );
            set( gca, 'LineWidth', LineWidth, 'FontSize', FontSize );
            xlabel('Contrast');
            ylabel('$\tau_s$ (ms)', 'interpreter', 'latex');


            axes( 'position', [0 0 1 1], 'visible', 'off' );
            text( 0.04, 0.75, 'On Cells', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'rotation', 90, 'FontSize', FontSize+5 );
            text( 0.04, 0.25, 'Off Cells', 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'rotation', 90, 'FontSize', FontSize+5 );

        end
    end
end