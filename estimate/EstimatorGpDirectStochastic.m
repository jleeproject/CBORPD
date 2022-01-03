classdef EstimatorGpDirectStochastic < AbsEstimatorGp
    % GP with noise variance estimation
    
    properties
        param_BO
        gpParam
        type = TypeEstimatorGp.DirectStochastic;
    end
    
    methods
        function this = EstimatorGpDirectStochastic(param_BO, gpSetting, xDomain)
            this.param_BO = param_BO;
            this.gpParam = gpSetting.gpParam;
%             this.xDomain = xDomain;
        end
        
        function [predictor, gpParam, truePredictor] = estimateGp(this, x_train, y_train, cumul_samplesizes, nouse)
            [funcHs, learnedHPs, funcCovModifiaiblePost] = gpDirectStochastic(x_train, y_train, this.gpParam, this.param_BO);
            predictor = PredictorGpDirect(funcHs, funcCovModifiaiblePost, this.param_BO);
            gpParam = GpParameters(learnedHPs, sqrt(learnedHPs.scale), learnedHPs.bw);
            truePredictor = predictor;
        end
    end
end

