classdef InfillOptimizerDirectConstWEI < AbsInfillOptimizer
    
    
    properties
        name = 'Weighted EI';
        optParam
        type = TypeInfillOptimizer.DirectConstWEI;
        pofCalculator
        isFirstRun = true;
        feasibilitySamplingConstraint
    end
    
%     methods (Abstract)
%         method2(this)
%     end
    methods
        function init(this)
%             this.gpHPs.diRectParams = 0
            this.optParam = gpDirectInit(this);
            
            if(~isa(this.acqusitionFunc,'AcqFuncEI'))
                fprintf('[WARNING] EI should be used for the acquisition function. (InfillOptimizerDirectConstWEI)\n');
                fprintf('[WARNING] Replace to EI\n');
                prevAcqFunc = this.acqusitionFunc;
                this.acqusitionFunc = AcqFuncEI(prevAcqFunc.param_BO, prevAcqFunc.objFunc, prevAcqFunc.isMaximizeObj);
                this.isMaximizeAcq = true;
            end
            this.pofCalculator = PoFCalculator(this.param_BO);
            this.feasibilitySamplingConstraint = FeasibilitySamplingConstraint(this.param_BO);
            this.acqFuncName = 'Weighted EI';
        end
        
%         function setFeasibilitySamplingConstraint(this, feasibilitySamplingConstraint)
%             this.feasibilitySamplingConstraint = feasibilitySamplingConstraint;
%         end
        
        function [wei, mu, sig] = weiWrapper(this, predictor, iter, multiConstrPredictors, arg)
            [ei, mu, sig] =   this.acqusitionFunc.acquire(predictor, arg, iter);
            wei = ei.*this.pofCalculator.getPof( multiConstrPredictors, arg) ;
            
        end
        
        
        function [optAcqX, optAcqV, idxOptAcq, histories, optMuX, optMuV, acquisition, funcFeas, problem, urnd ] = find(this, predictor, iter, multiConstrPredictors)
            if(this.isFirstRun)
                acquisition =  @(arg) this.pofCalculator.getPof( multiConstrPredictors, arg);
                [nextPtAcq, nextPt, history_imprv, queries_all] = ...
                diRectWrapMax_module( acquisition, this.objFunc.getXDomain(), this.optParam);
            
                optMuV = this.acqusitionFunc.acquire(predictor, nextPt, iter);
                this.acqusitionFunc.updateFOpt(optMuV);
                this.isFirstRun = false;
            end
            
%             acquisition =  @(arg) (this.acqusitionFunc.acquire(predictor, arg, iter).* this.pofCalculator.getPof( multiConstrPredictors, arg) );
            acquisition =  @(arg) weiWrapper(this, predictor, iter, multiConstrPredictors, arg);
    %         [funcHs, learnedHPs] = gp_direct(x, y, cumul_samplesizes, gpHPs);
            % First maximise the MF-GP-UCB acquisition function.
%             acquisition = @(arg) acqGP_UCB_min(arg, funcHs, iter, gpHPs.zetas);
            % acquisition returns: [acq, uncerts, mu, sigma]
%             if(this.isMaximizeAcq)
                [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq, queryVals_mu, queryVals_sigma] = ...
                DiRectWrap_max_module_nargout_3( acquisition, this.objFunc.getXDomain(), this.optParam);
%                 [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq] = ...
%                 diRectWrapMax_module( acquisition, this.objFunc.getXDomain(), this.optParam);
%             else
%                 [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq, queryVals_mu, queryVals_sigma] = ...
%                 DiRectWrap_min_module_nargout_3( acquisition, this.objFunc.getXDomain(), this.optParam);
% %                 [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq] = ...
% %                 diRectWrapMin_module( acquisition, this.objFunc.getXDomain(), this.optParam);
%             end
%             [optAcqV] = acquisition(nextPt);
            optAcqV = nextPtAcq;
            optAcqX = nextPt;
            idxOptAcq = [];
            
            histories.improv = history_imprv;
            histories.all  = queries_all;
            histories.acq  = queryVals_acq;
            histories.mu   = queryVals_mu;
            histories.sigma = queryVals_sigma;
            
            funcFeas = @(arg, urnd) this.feasibilitySamplingConstraint.isInFeasible(multiConstrPredictors, arg, urnd);
            problem = struct();
            urnd = [];

%             xxs = queries_all(history_imprv(:,2),:);
%             if(this.isMaximizeObj)
%                 [optMuV,idxMax] = max(queryVals_mu(history_imprv(:,2)));
%             else
%                 [optMuV,idxMax] = min(queryVals_mu(history_imprv(:,2)));
%             end
%             optMuX = xxs(idxMax,:);
            optMuX = optAcqX;
            optMuV = optAcqV;

            %% Not implement maximizing mu. 
            % For this, the acqusition function needs to be defined into three outputs.
            % For that, pof calculator needs to be rows of [prob, 1 1].
            % This will make useless overhead of calculations
            % If this is implemented, diRectWrap needs to be modified into 3 outs.
%             if(this.isMaximizeObj)
%                 [ optMuV ,idxOpt] = max(histories.mu);
%             else
%                 [ optMuV ,idxOpt] = min(histories.mu);
%             end
%             optMuX = queries_all(idxOpt,:);
%             optMuX = optAcqX;
%             optMuV = nextPtAcq;
%             [mu_xnew, ~] = predictor.predict(nextPt);
%             f_opt = mu_xnew; %% SEE IF THIS CAN GET into infillSampler
            this.acqusitionFunc.updateFOpt(optMuV);
%             [mu, ~, sigma_xnew] = funcHs(xnew);
%             % [TODO] May need visualization over Acquisition Function Search
%             [f_min, idx_min_mu] = min(queryVals_mu);
%             fmin_history(nEvaluations+1,1) = f_min;
%             xmin_history(nEvaluations+1,:) = queries_all(idx_min_mu,:);
        end
%         function ready = isReady(this,inputArg)
%             if(contains(this.acqFuncName,'EI','IgnoreCase',true))
%                 ready = this.Property1 + inputArg;
%             elseif(contains(this.acqFuncName,'UCB','IgnoreCase',true))
%             end
%         
        function optParam = gpDirectInit(this)
        %     diRectParams.maxevals = ceil(7 * min(5,numDims)^2 * sqrt(min(iter, 1000)));
        %     diRectParams.maxevals = 200;
            optParam.maxevals = this.param_BO.maxOptEval;
            optParam.maxits = this.param_BO.maxOptIter;
        %     fprintf('t = %d, diREctEvals: %d\n', boIter, diRectParams.maxevals);
        end

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