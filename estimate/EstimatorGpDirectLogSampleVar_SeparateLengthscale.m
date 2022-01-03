classdef EstimatorGpDirectLogSampleVar_SeparateLengthscale < AbsEstimatorGp
    %ESTIMATORGP_GPML Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        param_BO
        gpParam
        type = TypeEstimatorGp.DirectLogSampleVar_SeparateLengthscale;
    end
    
    methods
        function this = EstimatorGpDirectLogSampleVar_SeparateLengthscale(param_BO, gpSetting)
            this.param_BO = param_BO;
            this.gpParam = gpSetting.gpParam;
        end
        
        function [predictor, gpParam, truePredictor] = estimateGp(this, x_train, y_train, cumul_samplesizes)
            if(sum(isnan(y_train))>0)
                idx = find(~isnan(y_train));
                x_train = x_train(idx,:);
                y_train = y_train(idx,:);
                fprintf('[WARNING] [EstimatorGpDirectLogSampleVar] y_train includes nan. forced to omit nan value.');
            end
            [funcHs, learnedHPs, funcCovModifiaiblePost] = gpDirectLogSampleVar_SeparateLengthscale(x_train, y_train, cumul_samplesizes, this.gpParam, this.param_BO);
            predictor = PredictorGpDirect(funcHs, funcCovModifiaiblePost, this.param_BO);
            gpParam = GpParameters(learnedHPs, sqrt(learnedHPs.scale), learnedHPs.bw);
            truePredictor = predictor;
        end
    end
end

