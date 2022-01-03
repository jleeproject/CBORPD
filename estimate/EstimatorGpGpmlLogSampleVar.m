classdef EstimatorGpGpmlLogSampleVar < AbsEstimatorGp
    %ESTIMATORGP_GPML Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        param_BO
        type = TypeEstimatorGp.GpmlLogSampleVar;
        %
        max_iter_num
        inf
        meanfunc
        covfunc
        likfunc
        hyp_input
        gpParam
    end
    
    methods
%         function this = EstimatorGpGpmlLogSampleVar(param_BO, meanfunc, covfunc, likfunc, inf, hyp_input, max_iter_num)
        function this = EstimatorGpGpmlLogSampleVar(param_BO, gpSetting, xDomain)
            this.gpParam = gpSetting.gpParam;
            gpParam = gpSetting.gpParam;
            this.param_BO =param_BO;
            this.max_iter_num = gpParam.max_iter_num;
            this.inf = gpParam.inf;
            this.meanfunc = gpParam.meanfunc;
            this.covfunc = gpParam.covfunc;
            this.likfunc = gpParam.likfunc;
            this.hyp_input = gpParam.hyp_input;
%             this.xDomain = xDomain;
        end
        
        function [predictor, gpParam, truePredictor] = estimateGp(this, x_train, y_train, cumul_samplesizes)
        %%%%% --------------------- ESTIMATE SURROGATE MODEL ---------------------------------------
            hyp = setfield(this.hyp_input, 'samplesize',cumul_samplesizes);
        %     hyp = struct('mean', [], 'cov', [0 0], 'samplesize',samplesize_arr);
        %                     hyp2 = minimize_jaesung(hyp, @gp, -max_iter_num, inf, meanfunc, covfunc, likfunc, x, y);
        %                     hyp2 = minimize_jaesung_UCB(hyp, @gp, -max_iter_num, inf, meanfunc, covfunc, likfunc, x, y);
            if(this.param_BO.isConservative())
                hyp2 = minimize_jaesung_UCB(hyp, @gp, -this.max_iter_num, this.inf, this.meanfunc, this.covfunc, this.likfunc, x_train, y_train);
            else
                hyp2 = minimize_jaesung(hyp, @gp, -this.max_iter_num, this.inf, this.meanfunc, this.covfunc, this.likfunc, x_train, y_train);
            end
            % add samplesize
            hyp = setfield(hyp2, 'samplesize',cumul_samplesizes);
        %                     hyp.mean = meanfuncV;
        %%%%% --------------------- EVALUATION POINTS ---------------------------------------
        %     xs;

        %%%%% --------------------- PREDICT BY SURROGATE MODEL ---------------------------------------
        funcGpmlPredict = @(xs) gp(hyp, this.inf, this.meanfunc, this.covfunc, this.likfunc, x_train, y_train, xs);
%             if(this.param_BO.isSearchDirect)
% %                     funcHs = @(xs) gp(hyp, inf, meanfunc, covfunc, likfunc, x, y, xs);
%                 funcGpmlPredict = @(xs) gp(hyp, this.inf, this.meanfunc, this.covfunc, this.likfunc, x, y, xs);
%                 funcHs = @(xs) transFunctionReturnVar2Std(funcGpmlPredict, xs);
%             elseif(this.param_BO.isSearchGridGpml)
%                 [mu sigma2] = gp(hyp, this.inf, this.meanfunc, this.covfunc, this.likfunc, x, y, xs);
%                 sigma = sqrt(sigma2);
%             end
            gpParam = GpParameters(hyp, exp(hyp.cov(2)), exp(hyp.cov(1)));
            predictor = PredictorGpGpml(funcGpmlPredict);
            truePredictor = predictor;

        end
    end
end

