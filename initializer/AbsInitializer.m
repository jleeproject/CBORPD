classdef AbsInitializer < matlab.mixin.Copyable
    properties (Abstract)
        type
    end
    
    properties
        simulator
        xDomain

        nInitSample
        samplesizeForVariance
        nConstraints
    
    end
    
    methods (Abstract)
        [yy_arr, xx_arr, samplesize_arr] = getSamples(this);
    end
    
    methods
        function this = AbsInitializer(simulator, xDomain, nInitSample, nConstraints, samplesizeForVariance)
            this.simulator = simulator;
            this.xDomain = xDomain;
            
            this.nInitSample = nInitSample;
            this.nConstraints = nConstraints;
            if(nargin>4)
                this.samplesizeForVariance = samplesizeForVariance;
            else
                this.samplesizeForVariance = [];
            end
            
        end
        
    end
end

