classdef EstimatorGpDirectStochasticNoiseVar_SeparateLengthscale < AbsEstimatorGp
    % GP using given noise variance
    
    properties
        param_BO
        gpParam
        type = TypeEstimatorGp.DirectLogSampleVar_SeparateLengthscale;
    end
    
    methods
        function this = EstimatorGpDirectStochasticNoiseVar_SeparateLengthscale(param_BO, gpSetting, xDomain)
            this.param_BO = param_BO;
            this.gpParam = gpSetting.gpParam;
            this.xDomain = xDomain;
        end
        
        function [predictor, gpParam, truePredictor] = estimateGp(this, x_train, y_train, cumul_samplesizes, fn_mu_sigma_predict)
            if(nargin<5)
                throwError('Noise Variance is not specified.');
            end
            noiseVars = exp(2.*fn_mu_sigma_predict(x_train))./cumul_samplesizes;
%             noiseVars = exp(2.*fn_mu_sigma_predict(x_train));
            
            
            [funcHs, learnedHPs, funcCovModifiaiblePost] = gpDirectGivenNoiseVariance_SeparateLengthscale(x_train, y_train, cumul_samplesizes, this.gpParam, this.param_BO, noiseVars);
            predictor = PredictorGpDirect(funcHs, funcCovModifiaiblePost, this.param_BO);
            gpParam = GpParameters(learnedHPs, sqrt(learnedHPs.scale), learnedHPs.bw);
            truePredictor = predictor;
        end
    end
end

