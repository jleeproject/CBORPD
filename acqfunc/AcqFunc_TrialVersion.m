classdef AcqFuncTesting < AbsAcquisitionFunc & handle
    %_CLASS_TEMPLATE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Gaussian_CDF = @(x) normcdf(x, 0,1);
        Gaussian_PDF = @(x) normpdf(x, 0,1);
        name = 'Testing';
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
        
        function [acq, mu, sigma] = acquire(this, objPredictor, x_pred, data)
            [mu, sigma] = objPredictor.predict(x_pred);
            if(numel(data)==0)
                data = mu-1.5.*sigma;
                delta = (data);
                stdDelta = delta./sigma;
%                 acq=(delta).*this.Gaussian_CDF(stdDelta)+sigma.*this.Gaussian_PDF(stdDelta);
                acq=this.Gaussian_CDF(stdDelta)+this.Gaussian_PDF(stdDelta);
                fprintf('method 1 : coded\n');
            else
                stdDelta = data./sigma;
                acq=this.Gaussian_CDF(stdDelta)+this.Gaussian_PDF(stdDelta);
                fprintf('method 1 : argument\n');
            end
            if(this.isFirst)
                this.f_opt = mu(1);
                this.isFirst = false;
            end
        end
        
        function updateFOpt(this, f_opt)
            this.f_opt = f_opt;
            this.isFirst =false;
            if(this.param_BO.debug)
                fprintf('[AcqFunc:Testing] f_opt = %.4g\n',f_opt);
            end
        end
    end
end

