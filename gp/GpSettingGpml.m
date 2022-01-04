classdef GpSettingGpml
    properties
        gpParam
    end
    methods

        function this = GpSettingGpml(type_gp_fit, conservative, type_gp_kernel, yy_arr, meanFuncConservative, xDomain, max_opt_iter)
            % -------------------------------------------------------------------------
            if(conservative)
                meanfunc = @meanConst;    
%                 switch type_acquisition
%                     case TypeAcquisition.EI
%                          meanfuncV = [mean(yy_arr)];                    % empty: don't use a mean function
%                     case TypeAcquisition.LCB
%                         minV = min(yy_arr);
%                         rangeV = range(yy_arr);
%                         maxV = max(yy_arr);
%                         meanV = mean(yy_arr);
%                         meanfuncV = meanFuncConservatives.LCB(meanV, minV, maxV, rangeV);
%                     otherwise
%                         throwError()
%                 end
                meanfuncV = meanFuncConservative(meanV, minV, maxV, rangeV);
            end
            switch type_gp_kernel 
                case TypeGpKernel.Matern52;
                    covfunc = @(arg) covMaterniso(5,arg);
                case TypeGpKernel.Matern;
                    covfunc = @covMaterniso;
                case TypeGpKernel.Se;
                    covfunc = @covSEiso;              % Squared Exponental covariance function
                otherwise
                    throwError(sprintf('Undefined covfunc:%s',(covfunc)))
            end

            switch type_gp_fit 
                case TypeEstimatorGp.GpmlLogSampleVar
                    likfunc = @likGauss_taylored_log_sample_std;              % Gaussian likelihood
                    % Inference Likelihood Function
                    if(conservative)
                        inf =  @infGaussLik_taylored_log_sample_std_noCalcMean;
                    else
                        meanfunc = [];  
                        meanfuncV = [];
                        inf =  @infGaussLik_taylored_log_sample_std;
                    end
                case TypeEstimatorGp.GpmlStochastic
                    likfunc = @likGauss;
                    inf =  @infGaussLik;
                    if(~param_BO.isConservative)
                        meanfunc = [];  
                        meanfuncV = [];
                    end
                otherwise
                    throwError(sprintf('Undefined type_gp_fit:%s',(type_gp_fit)))
            end
            hyp = struct('mean', meanfuncV, 'cov', [0 0]);
            hyp_input = hyp; % hyp_input: without samplesize

            gpParam.meanfunc = meanfunc;
            gpParam.covfunc = covfunc;
            gpParam.likfunc = likfunc;
            gpParam.inf = inf;
            gpParam.hyp_input = hyp_input;
            gpParam.max_iter_num = max_opt_iter;
            
            this.gpParam = gpParam;
        end
        
        function out = getGpParam(this)
            out = this.setting.gpParam;
        end

    end
    
end

