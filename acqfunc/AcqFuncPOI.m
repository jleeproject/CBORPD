classdef AcqFuncPOI < AbsAcquisitionFunc & handle
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
        
        function [acq, mu, sigma] = acquire(this, predictor, x_pred, iter)
            [mu, sigma] = predictor.predict(x_pred);
            if(this.isFirst)
                this.f_opt = mu(1);
                this.isFirst = false;
            end
            if(this.isMaximizeObj)
                f_max = this.f_opt;
                delta = -(mu-f_max);
                stdDelta = delta./sigma;
                acq=this.Gaussian_CDF(stdDelta);
            else
                f_min = this.f_opt;
                delta = (f_min-mu);
                stdDelta = delta./sigma;
                acq=this.Gaussian_CDF(stdDelta);
            end
%             [optAcqV, idxOptAcq] = max(acq);
%             optAcqX = x_pred(idxOptAcq);
%            [optAcqV, optAcqX, idxOptAcq]
        end
        
        function updateFOpt(this, f_opt)
            this.f_opt = f_opt;
            this.isFirst = false;
            if(this.param_BO.debug)
                fprintf('[AcqFuncEI] f_opt = %.4g\n',f_opt);
            end
        end
    end
end

