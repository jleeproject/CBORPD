classdef AbsInfillOptimizer < handle & matlab.mixin.Copyable
    
    properties (Abstract)
        name
        type
    end
    
    properties
        acqusitionFunc
        param_BO
        objFunc
        
        acqOptMethod % Direct / Grid
        
        acqFuncName
        isMaximizeAcq % If it needs to find the maximum of acquisition function.
        isMaximizeObj
        % For EI, UCB, LCB, this is the same as those
        
        isSetIfMaximize = false;
        init_y = [];
        meanFuncConservative
        isReady = false;
        constraints = [];
        
        strStep = '';

        xDomain
        xDim  % number of parameters
    end
    
    methods (Abstract)
        init(this);
        find(this, predictor, iter)
    end
    
    methods
        function this = AbsInfillOptimizer(param_BO, acqusitionFunc, optProb)
            this.objFunc = optProb.objFunc;
            if(optProb.nCon>0)
                this.constraints = optProb.con;
            else
            end
            % Construct an instance of this class
            this.acqusitionFunc = acqusitionFunc;
            this.param_BO = param_BO;
            
            this.acqFuncName = acqusitionFunc.name;
            this.isMaximizeAcq = acqusitionFunc.isMaximizeAcq; 
            this.isMaximizeObj = acqusitionFunc.isMaximizeObj; 
            
            if( isa(optProb.objFunc,'AbsFunction'))
                this.xDomain = optProb.objFunc.getXDomain();
            elseif( isa(optProb.objFunc,'GrapheneModelSolver'))
                this.xDomain = optProb.objFunc.decisionVariables.domains.all;
            else
                error('[AbsInfillOptimizer] Undefined Type.');
            end
            this.xDim = size(this.xDomain,1);
%             this.init_y = initial_ys;
%             this.meanFuncConservative = meanFuncConservative;
            init(this)
        end
        
%         function prepare(this, initial_ys, meanFuncConservative)
%             this.init_y = initial_ys;
%             this.meanFuncConservative = meanFuncConservative;
%             this.isReady = true;
%         end
    end
end

% -----------------------------------------------------------
%         infillSampler.sample(predictor, f_min, iter)
% 
%         if(param_BO.isSearchGridGpml())
%         %% Acquisition Function
%             if(param_BO.isTypeAcquisitionEI())
%                 % calcuate the EI value
%                 acq=(f_min-mu).*Gaussian_CDF((f_min-mu)./sigma)+sigma.*Gaussian_PDF((f_min-mu)./sigma);
%                 [optAcqV, idxOptAcq] = max(acq);
%             elseif(param_BO.isTypeAcquisitionUCB());
%                 beta_t = 0.2 * fn_dim_eval * log(2*fn_dim_eval*iter);
%                 uncerts = sqrt(beta_t) * sigma;
%                 acq =  mu - uncerts;
% %                     acq = -acq;
%                 [optAcqV, idxOptAcq] = min(acq);
% 
%                 hist_beta_t(nEvaluations+1,1) = beta_t;
% 
%             else
%                 fprintf('[ERROR] Wrong Acquisition input\n');
%             end
%             xnew = xs(idxOptAcq,:);
%             sigma_xnew = sigma(idxOptAcq,1);
%             % ### --- Grid Search Result ---
%             % xmin_history = argmin E[f|y]
%             % fmin_history = min E[f|y])
%             [f_min, idx_min_mu] = min(mu);
%             fmin_history(nEvaluations+1,1) = f_min;
%             xmin_history(nEvaluations+1,:) = xs(idx_min_mu,:);
% %         elseif(strcmp('direct',gp_package))
%         elseif(param_BO.isSearchDirect())
%     %         [funcHs, learnedHPs] = gp_direct(x, y, cumul_samplesizes, gpHPs);
%             % First maximise the MF-GP-UCB acquisition function.
%             acquisition = @(arg) acqGP_UCB_min(arg, funcHs, iter, gpHPs.zetas);
%             % acquisition returns: [acq, uncerts, mu, sigma]
%             [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq, queryVals_uncert, queryVals_mu, queryVals_sigma] = ...
%             diRectWrap_min_func_nargout_4( acquisition, fn_x_domain, gpHPs.diRectParams);
%             [nextPtAcq, uncerts] = acquisition(nextPt);
%             xnew = nextPt;
% %                 final_est_idx = xnew;
% 
%             [mu, ~, sigma_xnew] = funcHs(xnew);
%             % [TODO] May need visualization over Acquisition Function Search
%             [f_min, idx_min_mu] = min(queryVals_mu);
%             fmin_history(nEvaluations+1,1) = f_min;
%             xmin_history(nEvaluations+1,:) = queries_all(idx_min_mu,:);
%         else
%             fprintf('[ERROR] wrong package name\n');
%             return;
%         end