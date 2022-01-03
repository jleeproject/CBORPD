classdef AcqFuncLCB < AbsAcquisitionFunc 
    %_CLASS_TEMPLATE Summary of this class goes here
    %   Detailed explanation goes here
    
%     properties
% %         beta_t = @(numDims, t) 0.2 * numDims * log(2*numDims*t);
%     end
    
    properties
        name = 'LCB';
        beta_t;
        type = TypeAcquisition.LCB;
        isMaximizeAcq = false;
        f_opt
    end
    
    methods
        function validate(this)
            if(this.isMaximizeObj)
                throwError('[ERROR] LCB:maximization, LCB is to minimize');
                return;
            end
        end
        
        function init(this)
            if( isa(this.objFunc,'AbsFunction'))
%                 this.xDomain = optProb.objFunc.getXDomain();
                numDims = this.objFunc.getDim();
            elseif( isa(this.objFunc,'GrapheneModelSolver'))
                numDims = size(this.objFunc.decisionVariables.domains.all,1);
%                 this.xDomain = optProb.objFunc.decisionVariables.domains.all;
            else
                error('[AbsInfillOptimizer] Undefined Type.');
            end
            this.beta_t = @(iter) 0.2 * numDims * log(2*numDims*iter);
        end
        
        function [acq, mu, sigma] = acquire(this, predictor, x_pred, iter)
%             [mu, sigma] = predictor(x_pred);
%             numDims = size(x_pred, 2);

            [mu, sigma] = predictor.predict(x_pred);
            uncerts = sqrt(this.beta_t(iter)) * sigma;
%             acq = mu - uncerts - zetas(1);
            acq = mu - uncerts;
%             [optAcqV, idxOptAcq] = max(acq);
%             optAcqX = x_pred(idxOptAcq);
%            [optAcqV, optAcqX, idxOptAcq]
        end
    end
end
