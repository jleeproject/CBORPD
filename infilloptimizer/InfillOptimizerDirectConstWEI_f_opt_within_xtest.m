classdef InfillOptimizerDirectConstWEI_f_opt_within_xtest < AbsInfillOptimizer
    
    
    properties
        name = 'Weighted EI (f_opt w/in xtest)';
        optParam
        type = TypeInfillOptimizer.DirectConstWEI;
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
            
            if(~isa(this.acqusitionFunc,'AcqFuncEI'))
                fprintf('[WARNING] EI should be used for the acquisition function. (InfillOptimizerDirectConstWEI)\n');
                fprintf('[WARNING] Replace to EI\n');
                prevAcqFunc = this.acqusitionFunc;
                this.acqusitionFunc = AcqFuncEI(prevAcqFunc.param_BO, prevAcqFunc.objFunc, prevAcqFunc.isMaximizeObj);
                this.isMaximizeAcq = true;
            end
            this.pofCalculator = PoFCalculator(this.param_BO);
            this.feasibilitySamplingConstraint = FeasibilitySamplingConstraint(this.param_BO);
            this.acqFuncName = 'Weighted EI (f_opt w/in xtest)';
            
            this.history_x = zeros(100,this.xDim);
            
            gpSettingMaker = this.param_BO.gpSettingMaker;
            addInitValues(this, gpSettingMaker.xx_arr)
        end
        
        function addInitValues(this, xx)
            nEvals = size(xx,1);
            this.history_x(1:nEvals,:) = xx;
            this.nEvals = nEvals;
        end
        
        function checkHistoryX(this)
            if(this.nEvals >= size(this.history_x,1))
                newMat = zeros(floor(this.nEvals*1.5),this.xDim);
                newMat(1:size(this.history_x,1),:) = this.history_x;
                this.history_x = newMat;
            end
        end
        
        function saveX(this,x)
            if(size(x,1)>1);error('[InfillOptimizerDirectConstNEI] only one point can be saved.');end;
            if(size(x,2)~=this.xDim);error('[InfillOptimizerDirectConstNEI] Dimension mismatch.');end;
            
            checkHistoryX(this)
            this.history_x(this.nEvals+1,:) = x;
            this.nEvals = this.nEvals + 1;
        end
        
        function [xx, nEvals ]= getHistoryX(this)
            nEvals = this.nEvals;
            xx = this.history_x(1:nEvals,:);
        end
        
%         function setFeasibilitySamplingConstraint(this, feasibilitySamplingConstraint)
%             this.feasibilitySamplingConstraint = feasibilitySamplingConstraint;
%         end
        
        function [wei, mu, sig] = weiWrapper(this, predictor, iter, multiConstrPredictors, arg)
            [ei, mu, sig] =   this.acqusitionFunc.acquire(predictor, arg, iter);
            wei = ei.*this.pofCalculator.getPof( multiConstrPredictors, arg) ;
            
        end
        
        
        function [optAcqX, optAcqV, idxOptAcq, histories, optMuX, optMuV, acquisition, funcFeas, problem, urnd ] = find(this, predictor, iter, multiConstrPredictors)
%             if(this.isFirstRun)
%                 acquisition =  @(arg) this.pofCalculator.getPof( multiConstrPredictors, arg);
%                 [nextPtAcq, nextPt, history_imprv, queries_all] = ...
%                 diRectWrapMax_module( acquisition, this.xDomain, this.optParam);
%             
%                 optMuV = this.acqusitionFunc.acquire(predictor, nextPt, iter);
%                 this.acqusitionFunc.updateFOpt(optMuV);
%                 this.isFirstRun = false;
%             end
            
            acquisition =  @(arg) weiWrapper(this, predictor, iter, multiConstrPredictors, arg);
            
            [ xx_train, ~] = getHistoryX(this);
            v_train = acquisition(xx_train);
            [~, idx_max ] = max(v_train);
            v_min = predictor.predict(xx_train(idx_max,:));
            this.acqusitionFunc.updateFOpt(v_min);
            
            acquisition =  @(arg) weiWrapper(this, predictor, iter, multiConstrPredictors, arg);

            [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq, queryVals_mu, queryVals_sigma] = ...
                DiRectWrap_max_module_nargout_3( acquisition, this.xDomain, this.optParam);

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

            saveX(this, nextPt);
            
            [ xx_train, ~] = getHistoryX(this);
            v_train = acquisition(xx_train);
%             [v_min, idx_min ] = max(v_train);
            [~, idx_max ] = max(v_train);
            optMuX = xx_train(idx_max,:);
            optMuV = predictor.predict(optMuX);
            this.acqusitionFunc.updateFOpt(optMuV);
          

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