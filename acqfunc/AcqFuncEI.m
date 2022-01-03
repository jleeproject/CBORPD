classdef AcqFuncEI < AbsAcquisitionFunc & handle
    %_CLASS_TEMPLATE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Gaussian_CDF = @(x) normcdf(x, 0,1);
        Gaussian_PDF = @(x) normpdf(x, 0,1);
        name = 'EI';
        f_opt;
        type = TypeAcquisition.EI;
        isMaximizeAcq = true;
        isFirst = true;
    end
    
    methods
        function validate(this)
            if(this.isMaximizeObj)
                disp('EI:maximization, need to double check')
            end
        end
        
        function init(this)
%             if(this.isMaximizeObj)
%                 this.f_opt = - Inf;
%             else
%                 this.f_opt = Inf;
%             end
        end
        
        function [acq, mu, sigma] = acquire(this, predictor, x_pred, iter, f_opt)
            [mu, sigma] = predictor.predict(x_pred);
            if(this.isFirst)
                this.f_opt = mu(1);
                this.isFirst = false;
            end
            if(nargin<5)
                f_opt = this.f_opt; 
            end
%             disp(this.f_opt);
            if(this.isMaximizeObj)
                f_max = f_opt;
%                 disp('EI:maximization, need to double check')
                delta = (mu-f_max);
                stdDelta = delta./sigma;
                acq=(delta).*this.Gaussian_CDF(stdDelta)+sigma.*this.Gaussian_PDF(stdDelta);
            else
                f_min = f_opt;
                delta = (f_min-mu);
                stdDelta = delta./sigma;
                acq=(delta).*this.Gaussian_CDF(stdDelta)+sigma.*this.Gaussian_PDF(stdDelta);
            end
%             [optAcqV, idxOptAcq] = max(acq);
%             optAcqX = x_pred(idxOptAcq);
%            [optAcqV, optAcqX, idxOptAcq]
        end
        
        function updateFOpt(this, f_opt)
            this.f_opt = f_opt;
            this.isFirst =false;
%             if(this.param_BO.debug)
%                 fprintf('[AcqFuncEI] f_opt = %.4g\n',f_opt);
%             end
        end
    end
end

