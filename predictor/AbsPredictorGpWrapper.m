classdef AbsPredictorGpWrapper < handle
    
    properties
        fnPredict
        fnSigModifiablePredict = [];
        param_BO
    end
    
    methods (Abstract)
        predict(this)
        %[mu,sig,covMat] = predict(this, xx)
    end
    methods
        function this = AbsPredictorGpWrapper(fnPredict, fnSigModifiablePredict, param_BO)
            this.fnPredict = fnPredict;
            this.param_BO  = param_BO;
            %
            if(nargin>1)
                this.fnSigModifiablePredict = fnSigModifiablePredict; % fnSigModifiablePredict(x,cov)
            end
        end
        
        function this = addNugget(this, nugget_proportion)
            if( numel(this.fnSigModifiablePredict)>0)
                this.fnPredict = this.fnSigModifiablePredict(nugget_proportion);
%             else
            end
        end
        
        function [A, mu] = getCholDecSig(this, xx)
            [mu,~,covMat] = predict(this, xx);
            try
                A  = chol(covMat,'lower');
            catch
                succeed = false;
%                 prevPredic = this.fnPredict;
                for i = 1:15
                    try
                        if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget 1e%d.*mean.\n', (-15+i));end;
%                         this.fnPredict = prevPredic;
                        nugget_proportion = 10^(-15+i);
                    
                        addNugget(this, nugget_proportion);
                        [mu,~,covMat] = predict(this, xx);
                        A  = chol(covMat,'lower');
%                         cov_new2 = covMat + eye(size(covMat,1)).*nugget;
%                         A  = chol(cov_new2 ,'lower');
                        succeed = true;
%                         addNugget(this, cov_new2);
                        break;
                    catch err
                        showErrors(err)
                    end
                end
                if(succeed)
                    if(this.param_BO.debug);fprintf('Succeed.\n');end
                else
                    fprintf('Failed in Cholesky Decomposition.\n');
                    showErrors(err);
                    error('cov is not positive definite (InfillOptimizerDirectConstNEI).');
                end
            end
            
            
        end
        
%         function [mu,sig,covMat] = method1(obj,inputArg)
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

