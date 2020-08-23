function [fr, fr_c, fr_s] = LinearResponse(obj, stimulus, inputX, inputY, eyeX, eyeY, rfParams, rfX, rfY)
	% Method to calculate response to given stimulus of spatial RFs with specified parameters and RF locations
	%  stimulus:	2nd dimension as horizontal (x), 1st dimension as vertical (y)
	%  inputX:   	horizontal coordinates of each pixel in the stimulus (degree)
	%  inputX:   	vertical coordinates of each pixel in the stimulus (degree)
    %  eyeX:        horizontal eye position (degree)
    %  eyeY:        vertical eye position (degree)
	%  rfParams:	structure array containing RF parameters
    %  rfX:			horizontal locations of RF centers (degree)
    %  rfY:			vertical locations of RF centers (degree)
    %
    %  fr:          net linear response; 1st dimension for neurons, 2nd dimension for eye positions
    %  fr_c:        center linear response; 1st dimension for neurons, 2nd dimension for eye positions
    %  fr_s:        surround linear response; 1st dimension for neurons, 2nd dimension for eye positions
    
    inputX = inputX(:);     % column vector
    inputY = inputY(:);

    r_c_sq = cat(1, rfParams.centerRadii).^2;
    K_c = cat(1, rfParams.centerPeakSensitivities);
    r_s_sq = cat(1, rfParams.surroundRadii).^2;
    K_s = cat(1, rfParams.surroundPeakSensitivities);

    %% method 1
    % [x, y] = meshgrid(inputX, inputY);
    % x = reshape(x,1,1,[]) - ( rfX(:) + eyeX(:)' );		% 1st dimension for RFs, 2nd for eyeX, 3rd for horizontal coordinates in stimulus
    % y = reshape(y,1,1,[]) - ( rfY(:) + eyeY(:)' );		% 1st dimension for RFs, 2nd for eyeY, 3rd for vertical coordinates in stimulus
    
    % [ax, ay] = meshgrid(gradient(inputX,1), gradient(inputY,1));
    % areas = ax(:) .* ay(:);

    % fr_c = reshape( reshape(K_c .* exp( -(x.^2+y.^2) ./ r_c_sq ), [], numel(stimulus)) * (stimulus(:) .* areas), length(rfX), [] );
    % fr_s = - reshape( reshape(K_s .* exp( -(x.^2+y.^2) ./ r_s_sq ), [], numel(stimulus)) * (stimulus(:) .* areas), length(rfX), [] );
    %% end method 1



    %% par-for-loop
    nNeurons = length(rfX);
    nEyePos = length(eyeX);
    fr_c = zeros(nNeurons,nEyePos);
    fr_s = fr_c;
    areas = gradient(inputX',1) .* gradient(inputY,1);

    %% method 2
    % rfX = repmat( rfX(:), nEyePos, 1 );     % make sliced variables for the par-for-loop
    % rfY = repmat( rfY(:), nEyePos, 1 );
    % K_c = repmat( K_c, nEyePos, 1 );
    % K_s = repmat( K_s, nEyePos, 1 );
    % r_c_sq = repmat( r_c_sq, nEyePos, 1 );
    % r_s_sq = repmat( r_s_sq, nEyePos, 1 );
    % eyeX = repmat( eyeX(:)', nNeurons, 1 );
    % eyeY = repmat( eyeY(:)', nNeurons, 1 );

    % n = 10;
    % r_s_sq_n = 10^2 * r_s_sq/2;      % +-n std
    % parfor( k = 0 : nNeurons * nEyePos-1 )
    %     distSq = (inputX' - (rfX(k) + eyeX(k))).^2 + (inputY - (rfY(k) + eyeY(k))).^2;
    %     idx = distSq(:) <= r_s_sq_n(k);     % within +-5 std; column vector
    %     fr_c(k) = K_c(k) * exp( -distSq(idx') / r_c_sq(k) ) * (stimulus(idx) .* areas(idx));
    %     fr_s(k) = - K_s(k) * exp( -distSq(idx') / r_s_sq(k) ) * (stimulus(idx) .* areas(idx));
    % end
    %% end method 2

    
    %% method 3
    n = 6;
    r_s_n = n * cat(1, rfParams.surroundRadii) / sqrt(2);
    r_s_sq_n = n^2 * r_s_sq/2;      % +-n std
    parfor( k = 0 : nNeurons * nEyePos-1 )
        iEyePos = floor(k/nNeurons) + 1;
        iNeuron = mod(k,nNeurons) + 1;

        % % circular area
        % distSq = (inputX' - (rfX(iNeuron) + eyeX(iEyePos))).^2 + (inputY - (rfY(iNeuron) + eyeY(iEyePos))).^2;
        % idx = distSq(:) <= r_s_sq_n(iNeuron);     % within +-n std; column vector
        % fr_c(k+1) = K_c(iNeuron) * exp( -distSq(idx') / r_c_sq(iNeuron) ) * (stimulus(idx) .* areas(idx));
        % fr_s(k+1) = - K_s(iNeuron) * exp( -distSq(idx') / r_s_sq(iNeuron) ) * (stimulus(idx) .* areas(idx));

        % square area
        dX = inputX' - (rfX(iNeuron) + eyeX(iEyePos));
        dY = inputY - (rfY(iNeuron) + eyeY(iEyePos));
        idxX = abs(dX) <= r_s_n(iNeuron);
        idxY = abs(dY) <= r_s_n(iNeuron);
        idx = idxX & idxY;
        distSq = dX(idxX).^2 + dY(idxY).^2;
        fr_c(k+1) = K_c(iNeuron) * exp( -distSq(:)' / r_c_sq(iNeuron) ) * (stimulus(idx) .* areas(idx));
        fr_s(k+1) = - K_s(iNeuron) * exp( -distSq(:)' / r_s_sq(iNeuron) ) * (stimulus(idx) .* areas(idx));
    end
    %% end method 3

    fr = fr_c + fr_s;
end