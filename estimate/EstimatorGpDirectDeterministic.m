classdef EstimatorGpDirectDeterministic < AbsEstimatorGp
    %ESTIMATORGP_GPML Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        param_BO
        gpParam
        type = TypeEstimatorGp.DirectDeterministic;
        
    end
    
    methods
        function this = EstimatorGpDirectDeterministic(param_BO, gpSetting, xDomain)
            this.param_BO = param_BO;
            this.gpParam = gpSetting.gpParam;
%             this.xDomain = xDomain;
        end
        
        function [predictor, gpParam, truePredictor] = estimateGp(this, x_train, y_train, cumul_samplesizes, nouse)
            [funcHs, learnedHPs, funcCovModifiaiblePost] = gpDirectDeterministic(x_train, y_train, this.gpParam, this.param_BO);
            predictor = PredictorGpDirect(funcHs, funcCovModifiaiblePost, this.param_BO);
            gpParam = GpParameters(learnedHPs, sqrt(learnedHPs.scale), learnedHPs.bw);
            truePredictor = predictor;
        end
    end
end

