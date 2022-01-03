classdef InfillOptimizerDirectConst_Testing3 < AbsInfillOptimizer
    %% Testing class. Will be transfer to a new class
    
    properties
        name = 'Feas:POF';
        optParam
        type = TypeInfillOptimizer.DirectConstTesting3;
        pofCalculator
        isFirstRun = true;
        feasibilitySamplingConstraint
%         optStep = true;
        optStep = false;
%         strStep = 'opt';
        
        acqEI
        acqLCB
        acqPOI
        cnt = 1;
        reqCnt = 1000;
    end
    
%     methods (Abstract)
%         method2(this)
%     end
    methods
        function init(this)
            this.acqFuncName = 'Testing';
%             this.gpHPs.diRectParams = 0
            this.optParam = gpDirectInit(this);
            
%             if(~isa(this.acqEI,'AcqFuncEI'))
%                 fprintf('[WARNING] EI should be used for the acquisition function. (InfillOptimizerDirectConstWEI)\n');
%                 fprintf('[WARNING] Replace to EI\n');
%                 prevAcqFunc = this.acqEI;
                this.acqEI = AcqFuncEI(this.param_BO, this.objFunc, false);
                this.acqLCB = AcqFuncLCB(this.param_BO, this.objFunc, false);
                this.acqPOI = AcqFuncPOI(this.param_BO, this.objFunc, false);
                this.isMaximizeAcq = true;
%             end
            this.pofCalculator = PoFCalculator(this.param_BO);
            this.feasibilitySamplingConstraint = FeasibilitySamplingConstraint(this.param_BO);
        end
        
%         function setFeasibilitySamplingConstraint(this, feasibilitySamplingConstraint)
%             this.feasibilitySamplingConstraint = feasibilitySamplingConstraint;
%         end
        
        function [wei, mu, sig] = objoptWrapper(this, objPredictor, iter, multiConstrPredictors, arg)
            [ei, mu, sig] =   this.acqEI.acquire(objPredictor, arg, iter);
%             [ei, mu, sig] =  this.acqLCB.acquire(objPredictor, arg, iter);
            wei = ei.*this.pofCalculator.getPof( multiConstrPredictors, arg) ;
            
        end
        
        function [acq, mu, sig] = feasWrapper(this, objPredictor, iter, multiConstrPredictors, xx)
                %% Candidate 1.
% %             [ei, mu, sig] =   this.acqEI.acquire(objPredictor, arg, iter);
% %             wei = ei.*this.pofCalculator.getPof( multiConstrPredictors, arg) ;
%             [poi_obj, mu, sig] =   this.acqPOI.acquire(objPredictor, arg, iter);
%             
%             pof_con = this.pofCalculator.getPof( multiConstrPredictors, arg) ;
%             [lcb_con, mu, sig] =   this.acqLCB.acquire(multiConstrPredictors, arg, iter);
% 
%             acq = poi_obj .* pof_con.*lcb_con;
                %% Candidate 2.
                pofs = this.feasibilitySamplingConstraint.pof( multiConstrPredictors, xx);
                
                
%                 [lcb_obj, mu, sig] =   this.acqLCB.acquire(objPredictor, arg, iter);
                [ei_obj, mu, sig] =   this.acqEI.acquire(objPredictor, xx, iter);
%                 [mu_con, sigma_con] = multiConstrPredictors.predict(xx);
% %                 data = mu_con;
%                 data = mu_con - sigma_con;
% 
%                 stdDelta = -data./sigma_con;
%                 acq=normcdf(stdDelta)+normpdf(stdDelta);
%                 mu = [];
%                 sig = [];
                acq = pofs;
%             
        end
        
        
        function [optAcqX, optAcqV, idxOptAcq, histories, optMuX, optMuV, acquisition, funcFeas, problem, urnd ] = find(this, objPredictor, iter, multiConstrPredictors)
            if(this.isFirstRun)
                acquisition =  @(arg) this.pofCalculator.getPof( multiConstrPredictors, arg);
                [nextPtAcq, nextPt, history_imprv, queries_all] = ...
                diRectWrapMax_module( acquisition, this.objFunc.getXDomain(), this.optParam);
            
                optMuV = this.acqEI.acquire(objPredictor, nextPt, iter);
%                 fprintf('optMuV : %s\n',optMuV);
                this.acqEI.updateFOpt(optMuV);
                this.acqusitionFunc.updateFOpt(optMuV);
                this.acqEI.updateFOpt(optMuV);
                this.isFirstRun =false;
                this.strStep = 'opt';
            end
            
            if(this.optStep)
                %% optimize step
                acquisition =  @(arg) objoptWrapper(this, objPredictor, iter, multiConstrPredictors, arg);
%                 acq = EI(obj).*PoF(con);
                [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq, queryVals_mu, queryVals_sigma] = ...
                DiRectWrap_max_module_nargout_3( acquisition, this.objFunc.getXDomain(), this.optParam);
%                 DiRectWrap_min_module_nargout_3( acquisition, this.objFunc.getXDomain(), this.optParam);
%                 
                this.cnt = this.cnt + 1;
                if(this.cnt >= this.reqCnt)
                    this.optStep = false;
                    this.cnt = 1;
                    this.strStep = 'feas';
                end
            else
                %% feasibility step
                %% Candidate 1.
                acquisition =  @(arg) feasWrapper(this, objPredictor, iter, multiConstrPredictors, arg);
                [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq, queryVals_mu, queryVals_sigma] = ...
                DiRectWrap_max_module_nargout_3( acquisition, this.objFunc.getXDomain(), this.optParam);
%                 DiRectWrap_min_module_nargout_3( acquisition, this.objFunc.getXDomain(), this.optParam);
                %% Candidate 2.
                
                
                this.cnt = this.cnt + 1;
                if(this.cnt >= this.reqCnt)
                    this.optStep = true;
                    this.cnt = 1;
                    this.strStep = 'opt';
                end
            end
            
            
            optAcqV = objPredictor.predict(nextPt);
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

            optMuX = optAcqX;
            optMuV = optAcqV;

            %% Not implement maximizing mu. 
            this.acqusitionFunc.updateFOpt(optMuV);
            this.acqPOI.updateFOpt(optMuV);
            this.acqEI.updateFOpt(optMuV);
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