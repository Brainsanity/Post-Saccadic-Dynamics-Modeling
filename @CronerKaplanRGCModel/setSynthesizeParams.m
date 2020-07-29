function setSynthesizeParams(obj, varargin)
	p = inputParser;
    p.addParameter('randomizeCenterRadii', true, @islogical);
    p.addParameter('randomizeCenterSensitivities', true, @islogical);
    p.addParameter('randomizeSurroundRadii', true, @islogical);
    p.addParameter('randomizeSurroundSensitivities', true, @islogical);
    p.parse(varargin{:});
    
    obj.synthesisOptions = struct( ...
        'randomizeCenterRadii', p.Results.randomizeCenterRadii, ...
        'randomizeCenterSensitivities', p.Results.randomizeCenterSensitivities, ...
        'randomizeSurroundRadii', p.Results.randomizeSurroundRadii, ...
        'randomizeSurroundSensitivities', p.Results.randomizeSurroundSensitivities);
end