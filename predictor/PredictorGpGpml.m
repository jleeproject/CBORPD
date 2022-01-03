classdef PredictorGpGpml < AbsPredictorGpWrapper
    
    properties
    end
    
    methods
       
        function [mu,sig,covMat] = predict(this, xx)
            [mu, var] = this.fnPredict(xx);
            sig = sqrt(var);
            covMat = [];
            
        end
    end
end

