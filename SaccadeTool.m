classdef SaccadeTool < handle
	
	methods ( Access = private )
		function obj = SaccadeTool()
		end
	end

	methods ( Static )
		function sac = Saccade()
			sac.start		= [];
			sac.duration	= [];
			sac.angle		= [];	% in degrees
			sac.amplitude	= [];
			sac.termiPoints	= single([]);	% 2-by-2: 1st column for start point; 2nd for end point
			sac.velocity	= single([]);
			sac.speed		= single([]);
			sac.peakVel	= single([]);
		end

		function [outData, debugSacs] = GetSacs( inData, varargin )
			%% outData = GetSacs( inData, varargin )
			%  analysis given eye traces and return all micro/saccades detected.
			%	inData: 		either structure array of trials, or
			%					a 2-row matrix in which the first row referring to horizontal eye position and the second row designating vertical eye position (in arcmin).
			%	sRate:  		sampling rate of eye traces; by default, 1000 Hz, or use 'sRate' in trial structure
			%	minMSA:			minimum microsaccade amplitude (arcmin); 3 by default
			%	maxMSA:			maximum microsaccade amplitude (arcmin); 30 by default
			%	minInterval:	minimum interval between micro/saccades (ms); 15 by default
			%	minDur:			minimum duration of micro/saccades (ms); 15 by default
			%   maxDur:			maximum duration of micro/saccades (ms); 350 by default
			%   minVel:			minimum peak velocity (velocity threshold) (arcmin/s); 180 by default
			%   minAcl:			minimum acceleration (acceleration threshold) (deg/s^2); 80 by default
			%   PLOT:			whether plot results; false by default
			%	
			%
			%	outData: 	structure array of trials with microsaccades and saccades, or
			%				a structure array containing microsaccades and saccades detected.
			%   debugSacs:	structure array containing all microsaccades and saccades (not classified), before minMSA, maxMSA, minInterval & minDur applied


			% Initialize all variables
			sRate = 1000;	% Hz
			minMSacAmplitude = 3;	% arcmin
			maxMSacAmplitude = 30;	% arcmin
			minInterval = 15;	% ms
			minDur = 15;		% ms
			maxDur = 350;		% ms
			minVelocity = 3*60;	% arcmin/s
			minAcl = 80;		% deg/s^2
			PLOT_FIGURE = false;

			if isfield( inData, 'sRate' )
				sRate = [];			% will be passed trial by trial later
			end

			% Interpret the user parameters
			k = 1;
			while k <= length(varargin) && ischar(varargin{k})
				switch (lower(varargin{k}))
					case 'srate'
						sRate = varargin{k+1};
					case 'minmsa'
						minMSacAmplitude = varargin{k+1};
					case 'maxmsa'
						maxMSacAmplitude = varargin{k+1};
					case 'mininterval'
						minInterval = varargin{k+1};
					case 'mindur'
						minDur = varargin{k+1};
					case 'maxdur'
						maxDur = varargin{k+1};
					case 'minvel'
						minVelocity = varargin{k+1};
					case 'minacl'
						minAcl = varargin{k+1};
					case 'plot'
						PLOT_FIGURE = varargin{k+1};
					otherwise
						error( 'outData = GetSacs( inData, varargin )' );
				end
				k = k + 2;
			end

			if isstruct(inData)		% array of trials
				for iTrial = length(inData) : -1 : 1
					if isempty(sRate)	% pass sRate trial by trial
						[sac, dbSacs] = SaccadeTool.GetSacs( [ inData(iTrial).x.position; inData(iTrial).y.position ], 'sRate', inData(iTrial).sRate, 'minMSA', minMSacAmplitude, 'maxMSA', maxMSacAmplitude, 'minInterval', minInterval, 'minDur', minDur, 'maxDur', maxDur, 'minVel', minVelocity, 'minAcl', minAcl );
					else
						[sac, dbSacs] = SaccadeTool.GetSacs( [ inData(iTrial).x.position; inData(iTrial).y.position ], 'sRate', sRate, 'minMSA', minMSacAmplitude, 'maxMSA', maxMSacAmplitude, 'minInterval', minInterval, 'minDur', minDur, 'maxDur', maxDur, 'minVel', minVelocity, 'minAcl', minAcl );
					end
					outData(iTrial) = inData(iTrial);
					fields = { 'microsaccades', 'saccades', 'start', 'duration', 'amplitude', 'angle', 'peakVel' };
					for i = 1 : 2
						for k = 3 : size(fields,2)
							outData(iTrial).(fields{i}).(fields{k}) = [ sac.(fields{i}).(fields{k}) ];
						end

						%% filter m/sacs with analysis, blinks, notracks and invalid
						index = false( size( outData(iTrial).(fields{i}).start ) );
						for k = 1 : size(index,2)
							index(k) = ~isIncludedIn( outData(iTrial).(fields{i}).start(k), outData(iTrial).(fields{i}).duration(k), outData(iTrial).analysis );
							for field = { 'blinks', 'notracks', 'invalid' }
								if ~index(k)
									index(k) = isIntersectedIn( outData(iTrial).(fields{i}).start(k), outData(iTrial).(fields{i}).duration(k), outData(iTrial).(field{1}) );
								end
							end
						end
						for k = 3 : size(fields,2)
							outData(iTrial).(fields{i}).(fields{k})(index) = [];
						end
					end

					debugSacs{iTrial} = dbSacs;
				end
				return;
			end

			minInterval = minInterval / 1000;	% convert from ms to s
			minDur = minDur / 1000;			% convert from ms to s
			maxDur = maxDur / 1000;		% convert from ms to s
			minVelocity = minVelocity / 60;		% convert from arcmin to degree

			saccades = struct( 'start', [], 'duration', [], 'angle', [], 'amplitude', [], 'peakVel', [], 'velocity', [], 'points', [] );
			saccades(1) = [];
			debugSacs = saccades;
			outData.saccades = saccades;
			outData.microsaccades = saccades;

		
			eyeTrace = inData / 60;	% input is eye trace; convert from arcmin to degrees
			if isempty(eyeTrace) || size(eyeTrace,1) ~= 2
				saccades = [];
				disp('EyeTrace must be a 2-row matrix!');
				return;
			end

			nDots = size(eyeTrace,2);				% number of sample dots
			if nDots / sRate < minDur				% too short sample duration
				return;
			end

			smoothEyeTrace = SaccadeTool.SmoothTrace( eyeTrace, round( 0.009 * sRate ) );	% sd of 0.009 s

			if PLOT_FIGURE
				FontSize = 16;
				LineWidth = 1.5;
				figure;
				subplot(2,2,1); hold on;
				plot( (0:nDots-1)/sRate*1000, eyeTrace(1,:), 'b-', 'DisplayName', 'X', 'LineWidth', 1 );
				plot( (0:nDots-1)/sRate*1000, eyeTrace(2,:), 'r-', 'DisplayName', 'Y', 'LineWidth', 1 );
				plot( (0:nDots-1)/sRate*1000, smoothEyeTrace(1,:), 'c-', 'DisplayName', 'Smoothed X', 'LineWidth', 1 );
				plot( (0:nDots-1)/sRate*1000, smoothEyeTrace(2,:), 'm-', 'DisplayName', 'Smoothed Y', 'LineWidth', 1 );
				xlabel('Time (ms)');
				ylabel('Eye position (deg)');
				set( gca, 'FontSize', FontSize, 'LineWidth', LineWidth );
				set( legend(), 'position', [0.0033 0.6038 0.0964 0.2919] );
			end

			%% calculate velocity
			velocity = SaccadeTool.SmoothTrace( gradient( smoothEyeTrace, 1/sRate ), round( 0.009 * sRate ) );	% vecctor
			speed = sqrt( velocity(1,:).^2 + velocity(2,:).^2 );	% quantity

			%% calculate acceleration
			acceleration = SaccadeTool.SmoothTrace( gradient( velocity, 1/sRate ), round( 0.009 * sRate ) );	% vector
			aclQuantity = sqrt( acceleration(1,:).^2 + acceleration(2,:).^2 );	% quantity
			aclOfSpeed = SaccadeTool.SmoothTrace( gradient( speed, 1/sRate ), round( 0.009 * sRate ) );	% acceleration of speed
			
			if PLOT_FIGURE
				subplot(2,2,2); hold on;
				ax = plotyy( (0:nDots-1)/sRate*1000, velocity(1,:), (0:nDots-1)/sRate*1000, acceleration(1,:) );
				xlabel('Time (ms)');
				ax(1).YLabel.String = 'X velocity (deg/s)';
				ax(2).YLabel.String = 'X acceleration (deg/s^2)';
				set( ax, 'FontSize', FontSize, 'LineWidth', LineWidth );
				
				subplot(2,2,3);
				ax = plotyy( (0:nDots-1)/sRate*1000, velocity(2,:), (0:nDots-1)/sRate*1000, acceleration(2,:) );
				xlabel('Time (ms)');
				ax(1).YLabel.String = 'Y velocity (deg/s)';
				ax(2).YLabel.String = 'Y acceleration (deg/s^2)';
				set( ax, 'FontSize', FontSize, 'LineWidth', LineWidth );

				subplot(2,2,4);
				[ax, h1, h2] = plotyy(  (0:nDots-1)/sRate*1000, speed, (0:nDots-1)/sRate*1000, aclOfSpeed );
				set( h1, 'DisplayName', 'Speed (deg/s)' );
				set( h2, 'DisplayName', 'Speed gradient (deg/s^2)' );

				set( ax, 'NextPlot', 'add' );
				plot( ax(2), (0:nDots-1)/sRate*1000, aclQuantity, 'c', 'DisplayName', 'Acceleration quantity (deg/s^2' );
				xlabel('Time (ms)');
				set( ax, 'FontSize', FontSize, 'LineWidth', LineWidth );
			end


			%% thresholding: look for saccades candidates
			mark = aclQuantity >= minAcl;	% pick out saccade candidates using an acceleration threshold
			mark(1) = 0;
			mark(end) = 0;

			if PLOT_FIGURE
				subplot(2,2,1);
				plot( (0:nDots-1)/sRate*1000, mark, '-', 'color', [0.5 0.5 0.5], 'DisplayName', 'Thresholding' );
			end

			mark = mark(2:end) - mark(1:end-1);
			iBounds = [ find( mark == 1 ); find( mark == -1 ) + 1 ];	% first row for start positions and second for end ones	
			if isempty(iBounds)
				saccades = [];
				return;
			end


			%% adjust candidate bounds according to peak acceleration, so that adjacent candidates from the same saccade have a higher chance to be merged together, while adjacent saccades have a lower chance to be merged together
			for k = 1 : size(iBounds,2)
				threshold = max( minAcl, max( aclQuantity( iBounds(1,k) : iBounds(2,k) ) ) * 0.05 );
				if aclQuantity(iBounds(1,k)) < threshold
					while aclQuantity(iBounds(1,k)) < threshold && iBounds(1,k) < nDots && aclQuantity(iBounds(1,k)+1) > aclQuantity(iBounds(1,k))
						iBounds(1,k) = iBounds(1,k) + 1;
					end
				else
					while aclQuantity(iBounds(1,k)) > threshold && iBounds(1,k) > 1 && aclQuantity(iBounds(1,k)-1) < aclQuantity(iBounds(1,k))
						iBounds(1,k) = iBounds(1,k) - 1;
					end
				end

				if aclQuantity(iBounds(2,k)) < threshold
					while aclQuantity(iBounds(2,k)) < threshold && iBounds(2,k) > 1 && aclQuantity(iBounds(2,k)-1) > aclQuantity(iBounds(2,k))
						iBounds(2,k) = iBounds(2,k) - 1;
					end
				else
					while aclQuantity(iBounds(2,k)) > threshold && iBounds(2,k) < nDots && aclQuantity(iBounds(2,k)+1) < aclQuantity(iBounds(2,k))
						iBounds(2,k) = iBounds(2,k) + 1;
					end
				end
			end

			if PLOT_FIGURE
				plot( ( [iBounds(1,:), nan, iBounds(2,:)] - 1 ) / sRate * 1000, [ones(1,size(iBounds,2))*2, nan, ones(1,size(iBounds,2))*2+0.1], 'r.', 'MarkerSize', 8, 'DisplayName', 'Adaptive thresh' );
			end

			
			%% split or merge candidates according to pre-defined profiles of aclQuantity & aclOfSpeed
			mark = ( gradient( aclQuantity, 1/sRate ) >= 0 );
			iTroughs_Q = find( mark(2:end) - mark(1:end-1) > 0 );	% troughs index for aclQuantity
			iPeaks_Q = find( mark(2:end) - mark(1:end-1) < 0 );		% peaks index for aclQuantity
			
			mark = ( gradient( aclOfSpeed, 1/sRate ) >= 0 );
			iTroughs_S = find( mark(2:end) - mark(1:end-1) > 0 );	% troughs index for aclOfSpeed
			iPeaks_S = find( mark(2:end) - mark(1:end-1) < 0 );		% peaks index for aclOfSpeed

			m = 1; k = 2;
			while k <= size(iBounds,2)

				%% split
				iTroughs = iTroughs_Q( iBounds(1,m) < iTroughs_Q & iTroughs_Q < iBounds(2,m) );		% troughs index for aclQuantity in candiate m
				if size(iTroughs,2) > 1		% presumably there should be only 1 trough in aclQuantity
					iCut = iTroughs(2);		% cut at the 2nd trough
				else		% check is there are actually extra peaks and troughs in aclOfSpeed; ideally there should be 1 troughs + 1 peak in aclOfSpeed (special case might exist)
					iPeaks = iPeaks_S( iBounds(1,m) < iPeaks_S & iPeaks_S < iBounds(2,m) );
					iPeaks = iPeaks( aclOfSpeed(iPeaks) > minAcl );
					iTroughs = iTroughs_S( iBounds(1,m) < iTroughs_S & iTroughs_S < iBounds(2,m) );
					iTroughs = iTroughs( -aclOfSpeed(iTroughs) > minAcl );
					if size(iPeaks,2) > 0 && size(iTroughs,2) > 0 && size(iPeaks,2) + size(iTroughs,2) > 2	% more than 1 trough + 1 peak in aclOfSpeed
						iCut = max( [iPeaks(1), iTroughs(1)] );
						mark = ( aclOfSpeed( iCut : iBounds(2,m) ) >= 0 );
						ex = find( mark(2:end) - mark(1:end-1), 1, 'first' ) + iCut;	% 1st zero crossing of aclOfSpeed
						if isempty(ex)
							ex = iBounds(2,m);
						end
						ex1 = iPeaks_S( find( iCut < iPeaks_S & iPeaks_S <= iBounds(2,m), 1, 'first' ) );	% 1st peak after iCut
						if ~isempty(ex1) && ex1 < ex
							ex = ex1;
						end
						ex1 = iTroughs_S( find( iCut < iTroughs_S & iTroughs_S <= iBounds(2,m), 1, 'first' ) );		% 1st trough after iCut
						if ~isempty(ex1) && ex1 < ex
							ex = ex1;
						end

						if ex < iBounds(2,m)
							iCut = ex;
						else
							iCut = nan;
						end

					else
						iCut = nan;
					end
				end
				
				if ~isnan(iCut)
					% left part
					iPeak_L = iPeaks_S( iBounds(1,m) <= iPeaks_S & iPeaks_S <= iCut );	% might be multiple ones
					iPeak_L = iPeak_L( aclOfSpeed(iPeak_L) > minAcl );
					iTrough_L = iTroughs_S( iBounds(1,m) <= iTroughs_S & iTroughs_S <= iCut );	% might be multiple ones
					iTrough_L = iTrough_L( -aclOfSpeed(iTrough_L) > minAcl );
					
					% right part
					iPeak_R = iPeaks_S( iCut <= iPeaks_S & iPeaks_S <= iBounds(2,m) );
					iPeak_R = iPeak_R( aclOfSpeed(iPeak_R) > minAcl );
					iTrough_R = iTroughs_S( iCut <= iTroughs_S & iTroughs_S <= iBounds(2,m) );	% might be multiple ones
					iTrough_R = iTrough_R( -aclOfSpeed(iTrough_R) > minAcl );

					if ~isempty(iPeak_L) && ~isempty(iTrough_L) &&...
					   iPeak_L(1) < iTrough_L(end) &&...											% left: peak at left side of trough, for the left part
					   max(abs(aclOfSpeed(iPeak_L))) / max(abs(aclOfSpeed(iTrough_L))) < 2.75 &&...	% left: peak/trough ratio small enough
					   max(abs(aclOfSpeed(iTrough_L))) / max(abs(aclOfSpeed(iPeak_L))) < 2.75 &&... % left: trough/peak ratio small enough
					   ~isempty(iPeak_R) &&...																				% right: contains at least 1 peak
					   ( isempty(iTrough_R) || max(abs(aclOfSpeed(iTrough_R))) / max(abs(aclOfSpeed(iPeak_R))) < 2.75 )		% right: contains no extreme low trough
					   																											% (if it does contain, then the right part is not likely to be a valid candidate, therefore not cut out)
						
						% cut out the left part
						iBounds = [ iBounds(:,1:k-1), [iCut; iBounds(2,m)], iBounds(:,k:end) ];	
						iBounds(2,m) = iCut;
						continue;
					end
				end

				%% merge if too close, and with some additional criteria (current candidate has no low trough or next candidate has no high peak)
				% m: current candidate
				iPeak_m = iPeaks_S( iBounds(1,m) <= iPeaks_S & iPeaks_S <= iBounds(2,m) );	% might be multiple ones
				iPeak_m = iPeak_m( aclOfSpeed(iPeak_m) > minAcl*0.4 );		% lower the threshold to 40%
				iTrough_m = iTroughs_S( iBounds(1,m) <= iTroughs_S & iTroughs_S <= iBounds(2,m) );	% might be multiple ones
				iTrough_m = iTrough_m( -aclOfSpeed(iTrough_m) > minAcl*0.4 );	% lower the threshold to 40%
				
				% k: next candidate
				iPeak_k = iPeaks_S( iBounds(1,k) <= iPeaks_S & iPeaks_S <= iBounds(2,k) );
				iPeak_k = iPeak_k( aclOfSpeed(iPeak_k) > minAcl*0.4 );	% lower the threshold to 40%
				iTrough_k = iTroughs_S( iBounds(1,k) <= iTroughs_S & iTroughs_S <= iBounds(2,k) );	% might be multiple ones
				iTrough_k = iTrough_k( -aclOfSpeed(iTrough_k) > minAcl*0.4 );	% lower the threshold to 40%

				if iBounds(1,k) - iBounds(2,m) <= minInterval*sRate &&...													% closer than minInterval
				   ( isempty(iPeak_m) || isempty(iTrough_m) ||...																% m: no peak or trough (even with a lower threshold)
				     iPeak_m(1) > iTrough_m(end) ||...																			% m: all peaks are at right side of all troughs
				     max(abs(aclOfSpeed(iPeak_m))) / max(abs(aclOfSpeed(iTrough_m))) > 3 ||...									% m: peak too high relative to trough
				     max(abs(aclOfSpeed(iTrough_m))) / max(abs(aclOfSpeed(iPeak_m))) > 3 ||...									% m: trough too low relative to peak
				     isempty(iPeak_k) ||...																						% k: no peak
				     ( ~isempty(iTrough_k) && max(abs(aclOfSpeed(iTrough_k))) / max(abs(aclOfSpeed(iPeak_k))) > 3 ) ) ||...		% k: extreme low trough exists
				   ...
				   iBounds(1,k) - iBounds(2,m) <= minInterval*2*sRate &&...													% closer than minInterval*2
				   ( size(iPeak_m,2) == 1 &&...																					% m: only 1 peak
				     isempty(iTrough_m) &&...																					% m: no trough
				     aclOfSpeed(iPeak_m) > minAcl*0.6 &&...																		% m: more strict threshold
				     isempty(iPeak_k) &&...																						% k: no peak
				     size(iTrough_k,2) == 1 &&...																				% k: only 1 trough
				     -aclOfSpeed(iTrough_k) > minAcl*0.6 )																		% k: more strict threshold

					% merge k into m
					iBounds(2,m) = iBounds(2,k);
				else
					% copy k into m+1
					m = m + 1;
					iBounds(:,m) = iBounds(:,k);
				end

				k = k + 1;
			end
			iBounds(:,m+1:end) = [];

			if PLOT_FIGURE
				plot( ( [iBounds(1,:), nan, iBounds(2,:)] - 1 ) / sRate * 1000, [ones(1,size(iBounds,2))*2.5, nan, ones(1,size(iBounds,2))*2.6], 'cs', 'DisplayName', 'Split & Merge' );
			end


			%% exclude invalid candidates
			for k = 1 : size(iBounds,2)
				iPeaks_q = iPeaks_Q( iBounds(1,k) <= iPeaks_Q & iPeaks_Q <= iBounds(2,k) );		% peaks in aclQuantity
				peaks_q = sort(aclQuantity(iPeaks_q));		% sort peaks
				if size(peaks_q,2) < 2 || peaks_q(end) / peaks_q(end-1) > 2.75		% presumably, there should be two roughly equal peaks in aclQuantity
					iBounds(:,k) = nan;
					continue;
				end

				iPeaks_s = iPeaks_S( iBounds(1,k) <= iPeaks_S & iPeaks_S <= iBounds(2,k) );		% peaks in aclOfSpeed
				iPeaks_s = iPeaks_s( aclOfSpeed(iPeaks_s) >= minAcl*0.4 );						% lower threshold to 40%
				iTroughs_s = iTroughs_S( iBounds(1,k) <= iTroughs_S & iTroughs_S <= iBounds(2,k) );	% troughs in aclOfSpeed
				iTroughs_s = iTroughs_s( -aclOfSpeed(iTroughs_s) >= minAcl*0.4 );					% lower threshold to 40%
				flag = true;
				%% check cases not to exclude
				if ~isempty(iPeaks_s) && ~isempty(iTroughs_s)		% no peak, no trough
					for iPeak = iPeaks_s
						if any( abs(aclOfSpeed(iPeak)) ./ abs(aclOfSpeed(iTroughs_s)) <= 2.75 &...	% peak not too high relative to trough
								abs(aclOfSpeed(iTroughs_s)) / abs(aclOfSpeed(iPeak)) <= 2.75 &...	% trough not too low relative to peak
								( iPeak < iTroughs_s &...																						% peak at left side of trough
								  max( speed( iBounds(1,k) : iBounds(2,k) ) ) - min( speed( iBounds(1,k) : iBounds(2,k) ) ) > minVelocity |...	% maximum speed exceeds minVelocity relative to minimum speed
								  aclOfSpeed(iPeak) >= minAcl & -aclOfSpeed(iTroughs_s) >= minAcl ) )											% peak at right side of trough, but peak & trough exceed normal threshold
							
							% should not be excluded
							flag = false;
							break;
						end
					end
				end

				if flag
					%% use maximum & minimum values for peak & trough, and check again
					[peak, iPeak] = max( aclOfSpeed( iBounds(1,k) : iBounds(2,k) ) );		% peak
					[trough, iTrough] = min( aclOfSpeed( iBounds(1,k) : iBounds(2,k) ) );	% trough
					if -trough < minAcl*0.4 ||...				% lower threshold to 40%
						peak < minAcl*0.4 ||...					% lower threshold to 40%
						abs(trough) / abs(peak) > 2.75 ||...	% trough too low relative to peak
						abs(peak) / abs(trough) > 2.75 ||...	% peak too high relative to trough
						( iPeak > iTrough ||...																									% peak at right side of trough
						  max( speed( iBounds(1,k) : iBounds(2,k) ) ) - min( speed( iBounds(1,k) : iBounds(2,k) ) ) < minVelocity ) &&...		% maximum speed is smaller than minVelocity relative to minimum speed
						( -trough < minAcl || peak < minAcl  ||...																				% trough or peak not exceed normal threshold
						  abs(trough) / abs(peak) > 2.75 || abs(peak) / abs(trough) > 2.75 )													% trough & peak too unequal
						
						% exclude candidate
						iBounds(:,k) = nan;
					end
				end
			end
			iBounds( :, isnan(iBounds(1,:)) ) = [];

			if PLOT_FIGURE
				plot( ( [iBounds(1,:), nan, iBounds(2,:)] - 1 ) / sRate * 1000, [ones(1,size(iBounds,2))*3, nan, ones(1,size(iBounds,2))*3.1], 'mo', 'DisplayName', 'Exlcuding' );
			end


			%% refine onset and offset time points
			velocity = SaccadeTool.SmoothTrace( gradient( eyeTrace, 1/sRate ), round( 0.004 * sRate ) );	% smaller SD for smoothing
			speed = sqrt( velocity(1,:).^2 + velocity(2,:).^2 );
			
			n1 = ceil( 0.01 * sRate );		% number of samples of 10ms
			n2 = ceil( 0.02 * sRate );		% number of samples of 20ms
			valueOn = abs( conv( [speed(1)*ones(1,n1), speed, speed(end)*ones(1,n1)], [zeros(1,n1), 1, -ones(1,n1)/n1], 'same' ) );		% onset: current speed - mean of previous 10ms
			valueOn = valueOn(n1+1:end-n1);
			valueOff = abs( conv( [speed(1)*ones(1,n2), speed, speed(end)*ones(1,n2)], [ones(1,n2)/n2, 0, -ones(1,n2)/n2], 'same' ) );	% offset: mean of successive 20ms - mean of previous 20ms
			valueOff = valueOff(n2+1:end-n2);
			if PLOT_FIGURE
				plot( ax(2), (0:nDots-1)/sRate*1000, valueOn, 'g', 'DisplayName', 'Onset criteria' );
				plot( ax(2), (0:nDots-1)/sRate*1000, valueOff, 'm', 'DisplayName', 'Offset criteria' );
				legend( [ get( ax(1), 'children' )', get( ax(2), 'children' )' ] );
			end
			
			for i = 1 : size(iBounds,2)
				bound = iBounds(:,i);

				%% for saccade onset
				x = ( gradient( valueOn( iBounds(1,i) : iBounds(2,i) ), 1/sRate) >= 0 );
				iPeaks = find( x(2:end) - x(1:end-1) < 0 );		% indices of peaks
				if ~isempty(iPeaks)
					iPeak = iPeaks( find( valueOn( iBounds(1,i) + iPeaks ) > max( valueOn( iBounds(1,i) + iPeaks ) ) * 0.25, 1, 'first' ) );		% leftmost of high peaks
				else
					iPeak = 1;
				end
				index = find( valueOn( iBounds(1,i) : iBounds(1,i)+iPeak-1 ) <= max( 2, 0.05 * max( valueOn( iBounds(1,i) : iBounds(1,i)+iPeak-1 ) ) ), 1, 'last' );		% last point below the threshold before iPeak (greater of 2 and 0.05*maxima)
				if ~isempty(index)
					bound(1) = iBounds(1,i) + index-1;
				end
				
				%% for saccade offset
				x = ( gradient( valueOff( iBounds(1,i) : iBounds(2,i) ), 1/sRate) >= 0 ) * 2 - 1;
				iPeaks = find( x(2:end) - x(1:end-1) < 0 );		% indices of peaks
				if ~isempty(iPeaks)
					iPeak = iPeaks( find( valueOff( iBounds(1,i) + iPeaks ) > max( valueOff( iBounds(1,i) + iPeaks ) ) * 0.25, 1, 'last' ) );		% rightmost of high peaks
				else
					iPeak = 1;
				end
				index = find( valueOff( iBounds(1,i)+iPeak-1 : iBounds(2,i) ) >= max( 2, 0.05 * max( valueOff( iBounds(1,i)+iPeak-1 : iBounds(2,i) ) ) ), 1, 'last' );	% last point above the threshold after iPeak (greater of 2 and 0.05*maxima)
				if ~isempty(index)
					bound(2) = iBounds(1,i) + iPeak-1 + index-1;
				end

				iBounds(:,i) = bound;

				% get sac properties
				saccades(i).start = iBounds(1,i);
				saccades(i).duration = iBounds(2,i) - iBounds(1,i);
				[ saccades(i).angle, saccades(i).amplitude ] = cart2pol( eyeTrace(1,iBounds(2,i)) - eyeTrace(1,iBounds(1,i)), eyeTrace(2,iBounds(2,i)) - eyeTrace(2,iBounds(1,i)) );
				saccades(i).angle = saccades(i).angle * 180 / pi;
				saccades(i).amplitude = saccades(i).amplitude * 60;		% convert degree to arcmin
				saccades(i).peakVel = max( speed( iBounds(1,i) : iBounds(2,i) ) );	% deg/s
				saccades(i).velocity = velocity( : , iBounds(1,i) : iBounds(2,i) );		% deg/s
				saccades(i).points = [ eyeTrace(1,iBounds(1,i)), eyeTrace(1,iBounds(2,i)); eyeTrace(2,iBounds(1,i)), eyeTrace(2,iBounds(2,i)) ];	% 1st row: start point; 2nd row: end point (degrees)
			end

			if( size(iBounds,2) ) < 1
				saccades = struct( 'start', [], 'duration', [], 'angle', [], 'amplitude', [], 'peakVel', [], 'velocity', [], 'points', [] );
				saccades(1) = [];
			end
			
			debugSacs = saccades;

			%% merge saccades with too short intervals
			m = 1;
			for k = 2 : size(saccades,2)
				if  saccades(k).start - ( saccades(m).start + saccades(m).duration ) >= minInterval * sRate
					m = m + 1;
					saccades(m) = saccades(k);
				else
					saccades(m).duration = saccades(k).duration + saccades(k).start - saccades(m).start;
					bounds = [ saccades(m).start, saccades(m).start + saccades(m).duration ];
					saccades(m).velocity = velocity( :, bounds(1) : bounds(2) );
					saccades(m).speed = speed( bounds(1) : bounds(2) );
					saccades(m).peakVel = max( saccades(m).speed );
					[ saccades(m).angle, saccades(m).amplitude ] = cart2pol( eyeTrace(1,bounds(2)) - eyeTrace(1,bounds(1)), eyeTrace(2,bounds(2)) - eyeTrace(2,bounds(1)) );
					saccades(m).angle = saccades(m).angle * 180 / pi;
					saccades(m).amplitude = saccades(m).amplitude * 60;
					saccades(m).points = [ eyeTrace(1,bounds(1)), eyeTrace(1,bounds(2)); eyeTrace(2,bounds(1)), eyeTrace(2,bounds(2)) ];
				end
			end
			saccades(m+1:end) = [];

			%% delete saccades with improper durations
			saccades( [saccades.duration] / sRate >= maxDur | [saccades.duration] / sRate <= minDur ) = [];

			if PLOT_FIGURE
				plot( ( reshape( [1,1,nan]' * iBounds(1,:), 1, [] ) - 1 ) / sRate * 1000, reshape( [get(gca,'YLim'), nan]' * ones(1,size(iBounds,2)), 1, [] ), 'g', 'DisplayName', 'Sac onset' );
				plot( ( reshape( [1,1,nan]' * iBounds(2,:), 1, [] ) - 1 ) / sRate * 1000, reshape( [get(gca,'YLim'), nan]' * ones(1,size(iBounds,2)), 1, [] ), '-', 'color', [0 0.7 0], 'DisplayName', 'Sac offset' );
				plot( (0:nDots-1)/sRate*1000, speed./50, 'DisplayName', 'Speed/50 (deg/s)' );
			end

			outData.saccades = saccades( [saccades.amplitude] > maxMSacAmplitude );
			outData.microsaccades = saccades( minMSacAmplitude <= [saccades.amplitude] & [saccades.amplitude] <= maxMSacAmplitude );
		end


		function trials = FilterNoise( trials, samples, sd_thresh, isPlot )
			if nargin() < 2 || isempty(samples)
				samples = trials;
			end
			if nargin() < 3 || isempty(sd_thresh)
				sd_thresh = 3;
			end
			if nargin() < 4 || isempty(isPlot)
				isPlot = false;
			end

			sacs = [trials.microsaccades, trials.saccades];
			data = [ [sacs.amplitude]'/60, [sacs.peakVel]', [sacs.duration]'/trials(1).sRate ];
			logData = log10(data);

			index = true( size(data,1), 1 );
			idx = [ 1, 2; 1, 3; 2, 3 ];
			for i = 3 : -1 : 1
				x = logData( :, idx(i,1) );
				y = logData( :, idx(i,2) );
				b(i,1:2) = polyfit( x, y, 1 );
				d = b(i,1)*x + b(i,2) - y;
				b(i,3) = sd_thresh * std(d);

				if i == 1
					index( d <= -b(i,3) | d >= b(i,3) ) = false;
				else
					index( d <= -b(i,3)/2 | d >= b(i,3) ) = false;
				end
			end

			if isPlot
				figure( 'NumberTitle', 'off', 'name', 'Main Sequence' );
				plot3( data(index,1)', data(index,2)', data(index,3)', 'k.', 'MarkerSize', 10 ); hold on;
				plot3( data(~index,1)', data(~index,2)', data(~index,3)', 'r.', 'MarkerSize', 10 );
				x = reshape( repmat( [min(logData(:,1))/2, max(logData(:,1))*2, nan]', 1, 3 ), 1, [] );
				plot3( 10.^x, 10.^( b(1,1)*x + b(1,2) + reshape( [1; 1; nan] * [-b(1,3), 0, b(1,3)], 1, []) ), repmat( max(data(:,3))*2, 1, 9 ) );
				plot3( 10.^x, repmat( max(data(:,2))*2, 1, 9 ), 10.^( b(2,1)*x + b(2,2) + reshape( [1; 1; nan] * [-b(2,3), 0, b(2,3)/2], 1, []) ) );
				x = reshape( repmat( [min(logData(:,2))/2, max(logData(:,2))*2, nan]', 3, 1 ), 1, [] );
				plot3( repmat( max(data(:,1))*2, 1, 9 ), 10.^x, 10.^( b(3,1)*x + b(3,2) + reshape( [1; 1; nan] * [-b(2,3), 0, b(2,3)/2], 1, []) ) );
				set( gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'log' );
				xlabel('Amplitude (deg)');
				ylabel('Peak velocity (deg/s)');
				zlabel('Duration (s)');
			end

			for iTrial = 1 : max(size(trials))
				fields1 = {'microsaccades', 'saccades'};
				for iField = 1 : 2
					sacs = trials(iTrial).(fields1{iField});
					if isempty(sacs.amplitude)
						continue;
					end
					data = log10( [ sacs.amplitude'/60, sacs.peakVel', sacs.duration'/trials(iTrial).sRate ] );
					index = true( size(sacs.amplitude) );
					for i = 1 : 3
						x = data( :, idx(i,1) );
						y = data( :, idx(i,2) );
						d = b(i,1)*x + b(i,2) - y;
						if i == 1
							index( d <= -b(i,3) | d >= b(i,3) ) = false;
						else
							index( d <= -b(i,3)/2 | d >= b(i,3) ) = false;
						end
					end
					fields2 = fields(sacs);
					for iF = 1 : size(fields2,1)
						trials(iTrial).(fields1{iField}).(fields2{iF})(~index) = [];
					end
				end
			end
		end


		function trace = SmoothTrace( trace, sd )
			%% smooth the trace with a guassian window with a standard deviation of sd
			convRadius = sd * 6;
			convFunctor = normpdf( -convRadius:convRadius, 0, sd );
			for i = 1 : size(trace,1)
				x = conv( [ ones(1,convRadius) * trace(i,1), trace(i,:), ones(1,convRadius) * trace(i,end) ], convFunctor, 'same' );
				trace(i,:) = x( convRadius+1 : end-convRadius );
			end
		end


		function Debug()
			global trials;
			global badTrials;
			global goodTrials;
			global badIndex;
			if(isempty(trials))
				folders = dir('F:\SmoothPursuit\DPI\RAW\A025');
				trials = [];
				for( i = 3 : size(folders,1) )
					load( fullfile( folders(i).folder, folders(i).name, 'Trials.mat' ) )
					trials = [trials, Trials];
				end
				% trials = SaccadeTool.GetSacs( trials, 'minmsa', 3 );
				% for iTrial = 1 : size(trials,2)
				% 	trial = trials(iTrial);
				% 	trial.x = trial.xRelative;
				% 	trial.y = trial.yRelative;
				% 	trial = SaccadeTool.GetSacs( trial, 'minmsa', 3 );
				% 	trial = findDrifts(trial);
				% 	trials(iTrial).msacsRelative = trial.microsaccades;
				% 	trials(iTrial).sacsRelative = trial.saccades;
				% 	trials(iTrial).driftsRelative = trial.drifts;
				% end

				badIndex = [251, 252, 253, 254, 256, 257, 260, 262, 269, 286, 290, 291, 321, 341, 346, 360, 381, 398, 401, 484, 527, 528, 536, 541, 543, 580, 612, 614, 619, 628, 630, 631, 633, 634, 635, 638, 640, 641, 645, 646, 648, 649, 650, 656, 657, 658, 664, 669, 670, 671, 681, 682, 684, 685, 686, 689, 696, 698, 699, 700, 702, 703, 705, 706, 708, 710, 711, 713, 716, 719, 750, 756, 759, 763, 765, 770, 776, 787, 800, 806, 808, 818, 820, 822, 823, 824, 825, 827, 828, 847, 849, 850, 851, 861, 864, 866, 868, 869, 873, 877, 878, 879, 880, 882, 883, 884, 886, 887, 889, 890, 891, 893, 894, 895, 900, 907, 908, 911, 914, 915, 916, 921, 922, 923, 924, 931, 932, 933, 934, 937, 939, 943, 945, 946, 947, 948, 950, 951, 952, 953, 954, 956, 957, 958, 961, 963, 964, 965, 966, 967, 968, 969, 970, 972, 973, 974, 979, 983, 985, 989, 992, 995, 996, 1004, 1005, 1007, 1008, 1009, 1010, 1015, 1016, 1017, 1019, 1020, 1021, 1022, 1023, 1026, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1038, 1039, 1040, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1059, 1060, 1063, 1064, 1065, 1066, 1067, 1068, 1070, 1072, 1073, 1076, 1078, 1079, 1081, 1086, 1087, 1088, 1089, 1090, 1094, 1096, 1097, 1103, 1107, 1110, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1120, 1157, 1158, 1159, 1164, 1167, 1172, 1185, 1194, 1198, 1199, 1200, 1201, 1202, 1204, 1217, 1218, 1219, 1220, 1221, 1222, 1226, 1230, 1231, 1234, 1235, 1246, 1251, 1253, 1254, 1256, 1257, 1258, 1259, 1260, 1262, 1263, 1264, 1265, 1270, 1271, 1272, 1274, 1275, 1277, 1278, 1289, 1293, 1295, 1296, 1298, 1299, 1302, 1303, 1304, 1305, 1306, 1309, 1310, 1311, 1315, 1316, 1317, 1319, 1320, 1326, 1327, 1338, 1347, 1357, 1363, 1373, 1374, 1376, 1377, 1381, 1389, 1392, 1398, 1399, 1400, 1402, 1405, 1411, 1412, 1413, 1417, 1418, 1420, 1422, 1424, 1425, 1426, 1427, 1445, 1460, 1472, 1474, 1480, 1481, 1482, 1487, 1488, 1490, 1501, 1504, 1512, 1518, 1519, 1520, 1523, 1526, 1530, 1531, 1532, 1538, 1541, 1542, 1544, 1545, 1546, 1547, 1551, 1552, 1555, 1557, 1559, 1561, 1563];
				badTrials = trials(badIndex);
				goodTrials = trials;
				goodTrials(badIndex) = [];
			end
			return;

			SmoothPursuit.CheckEyeTrace(badTrials);
			SmoothPursuit.CheckEyeTrace(goodTrials);

			set( figure, 'NumberTitle', 'off', 'name', 'Main Sequence | Population' ); hold on;
			sacs = [goodTrials.microsaccades, goodTrials.saccades];
			plot3( [sacs.amplitude]/60, [sacs.peakVel], [sacs.duration]/1000, '.', 'markersize', 5, 'markeredgecolor', 'b' );
			sacs = [badTrials.microsaccades, badTrials.saccades];
			plot3( [sacs.amplitude]/60, [sacs.peakVel], [sacs.duration]/1000, '.', 'markersize', 5, 'markeredgecolor', 'r' );
			set( gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'log' );
			xlabel('Amplitude (deg)');
			ylabel('Peak velocity (deg/s)');
			zlabel('Duration (s)');



			sacs = [goodTrials.microsaccades, goodTrials.saccades];
			data = [[sacs.amplitude]/60; sacs.peakVel; [sacs.duration]/1000];
			data = data - mean(data,2);
			data * data';

return;
			set( figure, 'NumberTitle', 'off', 'name', 'Main Sequence | Individual dots' ); hold on;
			for( iTrial = 1 : size(goodTrials,2) )
				for( iSac = 1 : size(goodTrials(iTrial).microsaccades.start,2) )
					plot3( [goodTrials(iTrial).microsaccades.amplitude(iSac)]/60, [goodTrials(iTrial).microsaccades.peakVel(iSac)], [goodTrials(iTrial).microsaccades.duration(iSac)]/1000, '.', 'markersize', 5, 'markeredgecolor', 'b', 'tag', sprintf('iT: %d; iM: %d', iTrial, iSac) );
				end
				for( iSac = 1 : size(goodTrials(iTrial).saccades.start,2) )
					plot3( [goodTrials(iTrial).saccades.amplitude(iSac)]/60, [goodTrials(iTrial).saccades.peakVel(iSac)], [goodTrials(iTrial).saccades.duration(iSac)]/1000, '.', 'markersize', 5, 'markeredgecolor', 'b', 'tag', sprintf('iT: %d; iS: %d', iTrial, iSac) );
				end
			end
			for( iTrial = 1 : size(badTrials,2) )
				for( iSac = 1 : size(badTrials(iTrial).microsaccades.start,2) )
					plot3( [badTrials(iTrial).microsaccades.amplitude(iSac)]/60, [badTrials(iTrial).microsaccades.peakVel(iSac)], [badTrials(iTrial).microsaccades.duration(iSac)]/1000, '.', 'markersize', 5, 'markeredgecolor', 'r', 'tag', sprintf('iT: %d; iM: %d', iTrial, iSac) );
				end
				for( iSac = 1 : size(badTrials(iTrial).saccades.start,2) )
					plot3( [badTrials(iTrial).saccades.amplitude(iSac)]/60, [badTrials(iTrial).saccades.peakVel(iSac)], [badTrials(iTrial).saccades.duration(iSac)]/1000, '.', 'markersize', 5, 'markeredgecolor', 'r', 'tag', sprintf('iT: %d; iS: %d', iTrial, iSac) );
				end
			end
			set( gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'log' );
			xlabel('Amplitude (\deg)');
			ylabel('Peak velocity (\deg/s)');
			zlabel('Duration (s)');
		end

		function DBPlotBT( trials, iTrial, M_L, iSac, marker, width, size, color )
			
			if( strcmpi( M_L, 'm' ) )
				% plot3( [trials(iTrial).microsaccades.amplitude(iSac)]/60, [trials(iTrial).microsaccades.peakVel(iSac)], [trials(iTrial).microsaccades.duration(iSac)]/1000, 'LineWidth', width, 'Marker', marker, 'markersize', size, 'markeredgecolor', color );
				plot( [trials(iTrial).microsaccades.amplitude(iSac)]/60, [trials(iTrial).microsaccades.peakVel(iSac)], 'LineWidth', width, 'Marker', marker, 'markersize', size, 'markeredgecolor', color );
			else
				% plot3( [trials(iTrial).saccades.amplitude(iSac)]/60, [trials(iTrial).saccades.peakVel(iSac)], [trials(iTrial).saccades.duration(iSac)]/1000, 'LineWidth', width, 'Marker', marker, 'markersize', size, 'markeredgecolor', color );
				plot( [trials(iTrial).saccades.amplitude(iSac)]/60, [trials(iTrial).saccades.peakVel(iSac)], 'LineWidth', width, 'Marker', marker, 'markersize', size, 'markeredgecolor', color );
			end

		end


	end

end