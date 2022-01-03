classdef InfillOptimizerGridConstSampling < AbsInfillOptimizer
    
    properties
        name = 'Grid';
        type = TypeInfillOptimizer.GridConstSampling;
    end
    
    properties        
        gridCoord
    end
    
%     methods (Abstract)
%         method2(this)
%     end
    methods
%         function this = InfillSamplerGrid(param_BO, acqusitionFunc)
%             % Construct an instance of this class
%             this.Property1 = inputArg1 + inputArg2;
%         end
        function init(this)
%             gridCoord = param_BO.
%             if(nargin==1)
%                 counts = this.param_BO.getGridBin();
%             end
%             domain = this.param_BO.getGridDomains(counts);
%             % COMBINATIONS
%             xs = domain(1,:);
%             for i = 2:this.param_BO.fn_dim_eval
%                 xs = combvec(xs, domain(i,:));
%             end
%             this.gridCoord = xs';
        end
        
        function [optAcqX, optAcqV, idxOptAcq, histories, optMuX, optMuV, acquisition] = find(this, predictor, iter)
%             if(~this.isReady)
%                 throwError('InfillSamplerGrid not ready. Need prepare() in advance.');
%             end
            acquisition = @(x) this.acqusitionFunc.acquire(predictor, x, iter);
            [acqs, mu, sig] = this.acqusitionFunc.acquire(predictor, this.gridCoord, iter);
            if(this.isMaximizeAcq)
                [valOpt, idxOpt] = max(acqs);
            else
                [valOpt, idxOpt] = min(acqs);
            end
            optAcqX = this.gridCoord(idxOpt);
            optAcqV = valOpt;
            idxOptAcq = idxOpt;

            histories.acq  = acqs;
            histories.mu   = mu;
            histories.sigma = sig;
            histories.uncert = uncert;
            
            if(this.isMaximizeObj)
                [ optMuV ,idxOptMu] = max(histories.mu);
            else
                [ optMuV ,idxOptMu] = min(histories.mu);
            end
            optMuX = queries_all(idxOptMu,:);
            
%             [mu_xnew, ~] = predictor.predict(nextPt);
%             f_opt = mu_xnew; %% SEE IF THIS CAN GET into infillSampler
            this.acqusitionFunc.updateFOpt(optMuV)
%             if(param_BO.isSearchGridGpml())
%             %% Acquisition Function
%                 if(param_BO.isTypeAcquisitionEI())
%                     % calcuate the EI value
%                     acq=(f_min-mu).*Gaussian_CDF((f_min-mu)./sigma)+sigma.*Gaussian_PDF((f_min-mu)./sigma);
%                     [optAcqV, idxOptAcq] = max(acq);
%                 elseif(param_BO.isTypeAcquisitionUCB());
%                     beta_t = 0.2 * fn_dim_eval * log(2*fn_dim_eval*iter);
%                     uncerts = sqrt(beta_t) * sigma;
%                     acq =  mu - uncerts;
%     %                     acq = -acq;
%                     [optAcqV, idxOptAcq] = min(acq);
%     
%                     hist_beta_t(nEvaluations+1,1) = beta_t;
%     
%                 else
%                     fprintf('[ERROR] Wrong Acquisition input\n');
%                 end
%                 xnew = xs(idxOptAcq,:);
%                 sigma_xnew = sigma(idxOptAcq,1);
%                 % ### --- Grid Search Result ---
%                 % xmin_history = argmin E[f|y]
%                 % fmin_history = min E[f|y])
%                 [f_min, idx_min_mu] = min(mu);
%                 fmin_history(nEvaluations+1,1) = f_min;
%                 xmin_history(nEvaluations+1,:) = xs(idx_min_mu,:);
%     %         elseif(strcmp('direct',gp_package))
%             elseif(param_BO.isSearchDirect())
%         %         [funcHs, learnedHPs] = gp_direct(x, y, cumul_samplesizes, gpHPs);
%                 % First maximise the MF-GP-UCB acquisition function.
%                 acquisition = @(arg) acqGP_UCB_min(arg, funcHs, iter, gpHPs.zetas);
%                 % acquisition returns: [acq, uncerts, mu, sigma]
%                 [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq, queryVals_uncert, queryVals_mu, queryVals_sigma] = ...
%                 diRectWrap_min_func_nargout_4( acquisition, fn_x_domain, gpHPs.diRectParams);
%                 [nextPtAcq, uncerts] = acquisition(nextPt);
%                 xnew = nextPt;
%     %                 final_est_idx = xnew;
%     
%                 [mu, ~, sigma_xnew] = funcHs(xnew);
%                 % [TODO] May need visualization over Acquisition Function Search
%                 [f_min, idx_min_mu] = min(queryVals_mu);
%                 fmin_history(nEvaluations+1,1) = f_min;
%                 xmin_history(nEvaluations+1,:) = queries_all(idx_min_mu,:);
%             else
%             end
        end
%         function ready = isReady(this,inputArg)
%             if(contains(this.acqFuncName,'EI','IgnoreCase',true))
%                 ready = this.Property1 + inputArg;
%             elseif(contains(this.acqFuncName,'UCB','IgnoreCase',true))
%             end
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