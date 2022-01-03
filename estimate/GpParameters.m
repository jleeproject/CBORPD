classdef GpParameters
    properties
        sigmaGp
        hyperparameter
        % Used For individual GP Package
        parameterObject
    end
    methods
        function this = GpParameters(parameterObject, sigmaGp, hyperparameter)
            if(nargin>1)
                this.sigmaGp = sigmaGp;
                this.hyperparameter = hyperparameter;
            end
            this.parameterObject = parameterObject;
        end
    end
end

