classdef AbsEstimatorGp < matlab.mixin.Copyable

    properties
        xDomain        
    end
    
    properties (Abstract)
        param_BO
        type
    end
    
    methods (Abstract)
        estimateGp(this, x_train, y_train, cumul_samplesizes)
    end
end

