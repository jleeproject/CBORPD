classdef InfillOptimizerDirectConstGPMean < AbsInfillOptimizer
    
    
    properties
        name = 'GP Mean';
        optParam
        type = TypeInfillOptimizer.DirectConstSobol;
        pofCalculator
        isFirstRun = true;
        feasibilitySamplingConstraint
        history_x
        nEvals = 0;
    end
    
%     methods (Abstract)
%         method2(this)
%     end
    methods
        function init(this)
%             this.gpHPs.diRectParams = 0
            this.optParam = gpDirectInit(this);
            
%             if(~isa(this.acqusitionFunc,'AcqFuncEI'))
%                 fprintf('[WARNING] EI should be used for the acquisition function. (InfillOptimizerDirectConstWEI)\n');
%                 fprintf('[WARNING] Replace to EI\n');
%                 prevAcqFunc = this.acqusitionFunc;
%                 this.acqusitionFunc = AcqFuncEI(prevAcqFunc.param_BO, prevAcqFunc.objFunc, prevAcqFunc.isMaximizeObj);
%                 this.isMaximizeAcq = true;
%             end
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
            urnd = 0
            idxOptAcq = [];
            %% NEXT POINT TO EVAL
            nextPt = net(sobolset(size(this.xDomain,1),'skip',iter),1);
            nextPtVal = predictor.predict(nextPt);
            

            optAcqV = nextPtVal;
            optAcqX = nextPt;
                       
            %% Select Current Optimal Value  ------------------------------------------------------------

            % if feasible then 0, infeasible then 1
            funcFeas = @(arg) this.feasibilitySamplingConstraint.feasGap(multiConstrPredictors.predict(arg));
            problem.constraint.func = funcFeas;
            problem.numconstraints = 1;
            problem.constraint.penalty = 10000;

            acquisition =  @(arg) predictor.predict(arg);

%             if(this.isMaximizeAcq)
%                 [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq, queryVals_mu, queryVals_sigma] = ...
%                 DiRectWrap_max_constrained_nargout_3( problem, acquisition, this.xDomain, this.optParam);
%             else
                [minVal, minPt, history_imprv, queries_all, queryVals_acq] = ...
                diRectWrapMin_module( problem, acquisition, this.xDomain, this.optParam);
%                 diRectWrapMin_given_constraint( problem, acquisition, this.xDomain, this.optParam);
%                 DiRectWrap_min_constrained_nargout_3( problem, acquisition, this.xDomain, this.optParam);
%             end
            [optAcqV] = acquisition(nextPt);
            optAcqX = nextPt;
            
            histories.improv = history_imprv;
            histories.all  = queries_all;
            histories.acq  = queryVals_acq;
            histories.mu   = minVal;
            histories.sigma = 0;
            
            
%             if sum(feas_x_test>0)>0
%             if(funcFeas(minPt))
%                 idx_feas = (feas_x_test>0);
%                 xx_sel_train = xx_train(idx_feas,:);
%                 v_train = acquisition(xx_sel_train);
%                 [v_min, idx_max ] = max(v_train);
            optMuX = minPt;
            optMuV = predictor.predict(minPt);
%             else
%                 feasgap_test = this.feasibilitySamplingConstraint.feasGap(mu_g);
%                 [~, idx_min ] = min(feasgap_test);
%                 optMuX = xx_train(idx_min,:);
%             end

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