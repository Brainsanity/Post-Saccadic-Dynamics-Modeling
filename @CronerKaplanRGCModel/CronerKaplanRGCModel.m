classdef CronerKaplanRGCModel < handle
    % Create a CronerKaplan RGC Model
    
    % References:
    %    Croner&Kaplan (1994). 'Receptive fields of P and M Ganglion cells 
    %    across the primate retina.',Vis. Res, (35)(1), pp.7-24
    % History:
    %    11/8/19  NPC, ISETBIO Team     Wrote it.
    
    properties (SetAccess = private)
        % Digitized data from Figure 4 & 5
        PCenterData;
        PSurroundData;
        MCenterData;
        MSurroundData;
        
        % Synthesized data
        PSynthesizedData;
        MSynthesizedData;
        
        % Model of center radius with eccentricity
        PCenterRadiusFunction;
        PCenterRadiusParams;
        PCenterRadiusParamsSE;
        PCenterSpacing2Radius;
        MCenterRadiusFunction;
        MCenterRadiusParams;
        MCenterRadiusParamsSE;
        MCenterSpacing2Radius;
        
        % Model of surround radius with eccentricity
        PSurroundRadiusFunction;
        PSurroundRadiusParams;
        PSurroundRadiusParamsSE;
        PSurroundSpacing2Radius;
        MSurroundRadiusFunction;
        MSurroundRadiusParams;
        MSurroundRadiusParamsSE;
        MSurroundSpacing2Radius;
        
        % Model of center sensitivity with center radius
        PCenterPeakSensitivityFunction;
        PCenterPeakSensitivityParams;
        PCenterPeakSensitivityParamsSE;
        MCenterPeakSensitivityFunction;
        MCenterPeakSensitivityParams;
        MCenterPeakSensitivityParamsSE;
        
        % Model of surround sensitivity with surround radius
        PSurroundPeakSensitivityFunction;
        PSurroundPeakSensitivityParams;
        PSurroundPeakSensitivityParamsSE;
        MSurroundPeakSensitivityFunction;
        MSurroundPeakSensitivityParams;
        MSurroundPeakSensitivityParamsSE;
        
        % Directory with psf deconvolution results
        psfDeconvolutionDir;
        
        synthesisOptions;
        
        plotlabOBJ;

        % fiiting parameters
        fitIntercept;
        dataSetToFit;
    end
    
    methods
        % Constructor
        function obj = CronerKaplanRGCModel(varargin) 
            % Parse input
            p = inputParser;
            p.addParameter('generateAllFigures', false, @islogical);
            p.addParameter('instantiatePlotLab', false, @islogical);
            p.addParameter('dataSetToFit', 'medians', @(x)(ismember(x, {'medians', 'raw', 'paperFormulas'})));
            p.addParameter('fitIntercept', false, @islogical);
            p.addParameter('randomizeCenterRadii', true, @islogical);
            p.addParameter('randomizeCenterSensitivities', true, @islogical);
            p.addParameter('randomizeSurroundRadii', true, @islogical);
            p.addParameter('randomizeSurroundSensitivities', true, @islogical);
            p.parse(varargin{:});
            
            obj.psfDeconvolutionDir = strrep(fileparts(which(mfilename())), ...
                '@CronerKaplanRGCModel', 'VisualToRetinalCorrectionData');

            obj.dataSetToFit = p.Results.dataSetToFit;
            obj.fitIntercept = p.Results.fitIntercept;
            
            obj.loadRawData();
            obj.fitModel();
            
            obj.synthesisOptions = struct( ...
                'randomizeCenterRadii', p.Results.randomizeCenterRadii, ...
                'randomizeCenterSensitivities', p.Results.randomizeCenterSensitivities, ...
                'randomizeSurroundRadii', p.Results.randomizeSurroundRadii, ...
                'randomizeSurroundSensitivities', p.Results.randomizeSurroundSensitivities);
            
            if (p.Results.instantiatePlotLab)
                obj.setupPlotLab();
            end
            
            if (p.Results.generateAllFigures)
                obj.plotDigitizedData();
            end
        end
        
        % Fit the model to a data set, either 'medians', or 'raw'
        fitModel(obj);
        
        % Method to synthesize data for a sample of eccentricity values
        synthesizedParams = SynthesizeRFParams(obj, cellType, eccDegs);

        % Method to generate RF params for given spacing
        rfParams = Spacing2RFParams(obj, cellType, spacing, temporalEccDegs);

        % Method to calculate response to given input of spatial RFs with specified parameters and RF locations
        [fr, fr_c, fr_s] = LinearResponse(obj, stimulus, inputX, inputY, eyeX, eyeY, rfParams, rfX, rfY);
        [sensitivity, centerSensitivity, surroundSensitivity] = SpatialSensitivity(obj, rfParams, SFs);

        % Method to display synthesized parameters
        DisplaySynthesizedRFParams(obj, PCells, MCells);

        DisplayFittedCenterSurroundRadiusRatio(obj);

        DisplaySpacing2RadiusFitting(obj);

        % Method to set synthesize parameters
        setSynthesizeParams(obj, varargin);
        
        % Method to simulate the Croner&Kaplan results
        simulateCronerKaplanResults(obj, varargin);
        
        % Compute the Gaussian-PSF convolution data
        performGaussianConvolutionWithPolansPSFanalysis(obj, varargin);
        
        % Method to generate retinal RF params given the retinal center radius
        % and eccentricity as inputs. This uses (via computeDeconvolutionModel()),
        % the Gaussian-PSF convolution data generated by performGaussianConvolutionWithPolansPSFanalysis().
        % It is to be used with mRGC mosaics whose centers are determined by connectivity to an underlying
        % cone mosaic.
        synthesizedRFParams = synthesizeRetinalRFparamsConsistentWithVisualRFparams(obj, retinalCenterRadii, retinalCenterMicrons, rgcIndices);
    end
    
    methods (Static)
        % plotSensitivities(theAxes, d, model, pointSize, color,displayYLabel, theLabel);
        % plotRadii(theAxes, d, model, pointSize, color, displayYLabel, theLabel);
        % [hEcc, vEcc, thePSFs, thePSFsupportDegs] = psfAtEccentricity(goodSubjects, imposedRefractionErrorDiopters, wavefrontSpatialSamples, eccXrange, eccYrange, deltaEcc);

        DisplayFittings();
        % DisplayFittedCenterSurroundRadiusRatio();
    end
    
    methods (Access=private)
        setupPlotLab(obj);
        
        % Method to compute the visual->retinal deconvolution model
        % based on the convolution results with the Polans data
        deconvolutionModel = computeDeconvolutionModel(obj);
    end
    
end

