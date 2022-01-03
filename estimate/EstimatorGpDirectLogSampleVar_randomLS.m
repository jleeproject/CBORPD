% classdef EstimatorGpDirectLogSampleVar_randomLS < AbsEstimatorGp
%     %ESTIMATORGP_GPML Summary of this class goes here
%     %   Detailed explanation goes here
%     
%     properties
%         param_BO
%         gpParam
%         type = TypeEstimatorGp.DirectLogSampleVar;
%         sigLS = 0.05;
%     end
%     
%     methods
%         function this = EstimatorGpDirectLogSampleVar_randomLS(param_BO, gpSetting, xDomain, sigLS)
%             this.param_BO = param_BO;
%             this.gpParam = gpSetting.gpParam;
%             this.xDomain = xDomain;
%             if(nargin>3)
%                 this.sigLS = sigLS;
%             end
%         end
%         
%         function [infillPredictor, gpParam, truePostPredictor] = estimateGp(this, x_train, y_train, cumul_samplesizes)
%             [funcHs, learnedHPs, truePostGp] = gpDirectLogSampleVar_randomLS(x_train, y_train, cumul_samplesizes, this.gpParam, this.param_BO);
%             infillPredictor = PredictorGpDirect(funcHs);
%             truePostPredictor = PredictorGpDirect(truePostGp);
%             gpParam = GpParameters(learnedHPs, sqrt(learnedHPs.scale), learnedHPs.bw);
%         end
%     end
% end
% 
