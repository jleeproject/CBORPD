classdef PredictorGpDirect < AbsPredictorGpWrapper
    
    properties
    end
    
    methods
       
        function [mu,sig,covMat] = predict(this, xx)
            [mu, covMat, sig] = this.fnPredict(xx);
        end
    end
end

