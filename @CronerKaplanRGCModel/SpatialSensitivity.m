function [sensitivity, C, S] = SpatialSensitivity(obj, rfParams, SFs)
	% sensitivity:		1st dim for rfParams, 2nd dim for SFs

	C = [rfParams.centerPeakSensitivities] .* pi .* [rfParams.centerRadii].^2 .* exp( -(pi * SFs(:) * [rfParams.centerRadii]).^2 );
	S = -[rfParams.surroundPeakSensitivities] .* pi .* [rfParams.surroundRadii].^2 .* exp( -(pi * SFs(:) * [rfParams.surroundRadii]).^2 );
	C = C';
	S = S';
	sensitivity = C + S;
end