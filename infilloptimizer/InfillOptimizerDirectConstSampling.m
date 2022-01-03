classdef InfillOptimizerDirectConstSampling < AbsInfillOptimizer
    
    
    properties
        name = 'ConstSampling';
        optParam
        type = TypeInfillOptimizer.DirectConstSampling;
        feasibilitySamplingConstraint
    end
    
%     methods (Abstract)
%         method2(this)
%     end
    methods
        function init(this)
%             this.gpHPs.diRectParams = 0
            this.optParam = gpDirectInit(this);
            this.feasibilitySamplingConstraint = FeasibilitySamplingConstraint(this.param_BO);
        end
        
%         function setFeasibilitySamplingConstraint(this, feasibilitySamplingConstraint)
%             this.feasibilitySamplingConstraint = feasibilitySamplingConstraint;
%         end
        
        
        
        function [optAcqX, optAcqV, idxOptAcq, histories, optMuX, optMuV, acquisition, funcFeas, problem, urnd ] = find(this, predictor, iter, multiConstrPredictors)
%             if(~this.isReady)
%                 throwError('InfillSamplerGrid not ready. Need prepare() in advance.');
%             end


            % at xx
            % Return 1 if xx is feasible.
            % feasible
%             predict = @(xx) multiConstrPredictors.predict(xx)
%             urnd = unifrnd();
            
            % g2<UB
            % g1>LB
            % OUTPUT : [row vector of UB-g2 or g1-LB]
            % INPUT: Each Constraint's : function, lb, ub (in param_BO), predictor, 
%             feasibilitySamplingConstraint
%             feasibilitySamplingConstraint = FeasibilitySamplingConstraint
            urnd = rand();
            problem.constraint.func = @(arg) this.feasibilitySamplingConstraint.isInFeasible(multiConstrPredictors, arg, urnd);
            problem.numconstraints = 1;
            problem.constraint.penalty = 10;
            
            funcFeas = @(arg, urnd) this.feasibilitySamplingConstraint.isInFeasible(multiConstrPredictors, arg, urnd);
%                    Problem.numconstraints = number of constraints
%                    Problem.constraint(i).func    = i-th constraint handle
%                    Problem.constraint(i).penalty = penalty parameter for
%                                                    i-th constraint

            acquisition =  @(arg) this.acqusitionFunc.acquire(predictor, arg, iter);
    %         [funcHs, learnedHPs] = gp_direct(x, y, cumul_samplesizes, gpHPs);
            % First maximise the MF-GP-UCB acquisition function.
%             acquisition = @(arg) acqGP_UCB_min(arg, funcHs, iter, gpHPs.zetas);
            % acquisition returns: [acq, uncerts, mu, sigma]
            if(this.isMaximizeAcq)
                [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq, queryVals_mu, queryVals_sigma] = ...
                DiRectWrap_max_constrained_nargout_3( problem, acquisition, this.xDomain, this.optParam);
            else
                [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq, queryVals_mu, queryVals_sigma] = ...
                DiRectWrap_min_constrained_nargout_3( problem, acquisition, this.xDomain, this.optParam);
            end
%             [nextPtAcq, optAcqV] = acquisition(nextPt);
            [optAcqV] = acquisition(nextPt);
            optAcqX = nextPt;
            idxOptAcq = [];

            if(funcFeas(nextPt,urnd)==1)
%                 disp('infeasible');
                prevPt = nextPt;
                probPred =  @(arg) this.feasibilitySamplingConstraint.pof(multiConstrPredictors, arg, urnd);
                [nextPtAcq, nextPt, history_imprv, queries_all] = ...
                diRectWrapMax_module( probPred, this.xDomain, this.optParam);
                if(this.param_BO.debug)
                    fprintf('Infeasible Problem with given %.4g, maximize PoF.\n', urnd);
                    fprintf('PoF : %.4g -> %.4g' ,probPred(prevPt), probPred(nextPt) );
                    fprintf(' (x: ', prevPt, nextPt );
                    fprintf('%.4g ', prevPt );
                    fprintf(' -> ');
                    fprintf('%.4g ', nextPt );
                    fprintf(')\n');
                end
            end
            
            idx = history_imprv(:,2)<this.optParam.maxevals;
            history_imprv = history_imprv(idx,:);
            
            histories.improv = history_imprv;
            histories.all  = queries_all;
            histories.acq  = queryVals_acq;
            histories.mu   = queryVals_mu;
            histories.sigma = queryVals_sigma;
%             histories.uncert = queryVals_uncert;

            
%             if(this.isMaximizeObj)
%                 [ optMuV ,idxOpt] = max(histories.mu);
%             else
%                 [ optMuV ,idxOpt] = min(histories.mu);
%             end
%             optMuX = queries_all(idxOpt,:);
            xxs = queries_all(history_imprv(:,2),:);
            if(this.isMaximizeObj)
                [optMuV,idxMax] = max(queryVals_mu(history_imprv(:,2)));
            else
                [optMuV,idxMax] = min(queryVals_mu(history_imprv(:,2)));
            end
            optMuX = xxs(idxMax,:);
            
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