classdef InfillOptimizerDirectConstNEI < AbsInfillOptimizer
    
    
    properties
        name = 'Noisy EI (Constrained)';
        optParam
        type = TypeInfillOptimizer.DirectConstNEI;
        pofCalculator
        isFirstRun = true;
        feasibilitySamplingConstraint
        quasi_sample_size = 5;
        penalty = 1000;
        
        history_x
        nEvals = 0;
        gpEstimatorObj 
        gpEstimatorCon
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
            this.acqFuncName = this.name;
            this.history_x = zeros(100,this.xDim);
            

            gpSettingMaker = this.param_BO.gpSettingMaker;
            gpSettingMaker.typeGpFit = TypeEstimatorGp.DirectDeterministic;
            
            gpSettingObj = gpSettingMaker.getGpSetting4ObjFunc;
            this.gpEstimatorObj = EstimatorGpDirectDeterministic(this.param_BO, gpSettingObj, this.xDomain) ;
            
            gpSettingsCon = gpSettingMaker.getGpSetting4Constraints;
            this.gpEstimatorCon = MultipleConstraintsEstimatorGP(this.param_BO, gpSettingsCon);
            this.gpEstimatorCon.addAllConstraints(this.param_BO.setting.optProb.con, TypeEstimatorGp.DirectDeterministic);
            
            addInitValues(this, gpSettingMaker.xx_arr);
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

        
        
        
        function [objPredictors_sampled, conPredictors_sampled, f_opts, hasFeasibleValue] = sampleModels(this, predictor, multiConstrPredictors)
            covMats = cell(multiConstrPredictors.nPredictor+1,1);
            mus = cell(multiConstrPredictors.nPredictor+1,1);
            
            [ xxs, n_xxs] = getHistoryX(this);
             
            [mu_obj,~,covMats_obj] = predictor.predict(xxs);
            for i=1:multiConstrPredictors.nPredictor
                [mus{i}, ~, covMats{i}] = multiConstrPredictors.cell_predictors{i}.predict(xxs);
            end
            mu_all = vertcat(mus{:});
            cov_all = blkdiag(covMats{:});
            [A ~] = chol(cov_all);
            
            objPredictors_sampled = cell(this.quasi_sample_size, 1);
            conPredictors_sampled = cell(this.quasi_sample_size, 1);
            f_opts = zeros(this.quasi_sample_size,1);
            hasFeasibleValue = zeros(this.quasi_sample_size,1);
            
%             tt = net(sobolset(1, 'Skip',1), n_xxs * (multiConstrPredictors.nPredictor + 1));
            cnt = 1;
            
            cs = zeros(n_xxs, multiConstrPredictors.nPredictor);
            for i = 1: this.quasi_sample_size
%                 tt = unifrnd(n_xxs * (multiConstrPredictors.nPredictor + 1) ,1)

%                 out = getNormalSample(this, tt, mu, cov);

                tt = net(sobolset(1, 'Skip',cnt), n_xxs); cnt = cnt + n_xxs;
                f = getNormalSample(this, tt, mu_obj, covMats_obj);

                for j=1:multiConstrPredictors.nPredictor
                    tt = net(sobolset(1, 'Skip',cnt), n_xxs); cnt = cnt + n_xxs;
                    cs(:,j) = getNormalSample(this, tt, mus{j}, covMats{j});
                end

%                 cs_vec = samples(n_xxs+1:end);
%                 cs = reshape(cs_vec, n_xxs, multiConstrPredictors.nPredictor);
                
                [~, ~, objPredictors_sampled{i}] = this.gpEstimatorObj.estimateGp( xxs, f, []);
                [~, ~, conPredictors_sampled{i}] = this.gpEstimatorCon.estimateGp( xxs, cs, [], []);
                
                feas = this.feasibilitySamplingConstraint.isThisValueFeasible(cs);
                if(sum(feas)>0)
                    idx_feas = find(feas);
% try
                    f_feas = f(idx_feas);
% catch err;
%     showErrors(err);
% end
%                     x_feas = xxs(idx_feas);
                    f_opts(i) = min(f_feas);
%                     x_min = 
                    hasFeasibleValue(i) = true;
                else
                    hasFeasibleValue(i) = false;
                end
            end
            
%             if(this.nEvals>10)
%                 disp();
%                 figure(11);clf;
%                 
% %                 xxxd1 = make2dRangesForNumericalStudy(0,6, 10, 0,6, 10);
% %                 xxxd2 = xxxd1';
%                 nBins = 50
%                 range = make1dRangesForNumericalStudy(0,6, nBins);
%                 
%                 xxx = combvec(range,range)';
%                 
%                 subplot(3,6,fnGetIdxOfSubplotWithRowCol([3, 6], 1,1));cla;
%                 imagesc(range, range, reshape(predictor.predict(xxx), nBins, nBins));
% 
%                 for i = 1:this.quasi_sample_size
%                 subplot(3,6,fnGetIdxOfSubplotWithRowCol([3, 6], 1,1+i));cla;
%                 imagesc(range, range, reshape(objPredictors_sampled{i}.predict(xxx), nBins, nBins));
%                 end
%                                 
%                 
% 
%                 for j=1: multiConstrPredictors.nPredictor
%                     subplot(3,6,fnGetIdxOfSubplotWithRowCol([3, 6], 1+j,1));cla;
%                     val = multiConstrPredictors.cell_predictors{j}.predict(xxx);
%                     imagesc(range, range, reshape(val, nBins, nBins));
%                     for i = 1:this.quasi_sample_size
%                     subplot(3,6,fnGetIdxOfSubplotWithRowCol([3, 6], 1+j,1+i));cla;
%                     val = conPredictors_sampled{i}.cell_predictors{j}.predict(xxx);
%                     imagesc(range, range, reshape(val, nBins, nBins));
%                     end
%                     
%                 end
%                 
%                 
%                 subplot(3,6,fnGetIdxOfSubplotWithRowCol([3, 6], 3,1));cla;
%                 ei = this.acqusitionFunc.acquire(predictor, xxx);
%                 pof = this.pofCalculator.getPof(multiConstrPredictors.cell_predictors{j}, xxx);
%                 val = ei.*pof;
%                 imagesc(range, range, reshape(val, nBins, nBins));
%                 for i = 1:this.quasi_sample_size
%                     subplot(3,6,fnGetIdxOfSubplotWithRowCol([3, 6], 3,1+i));cla;
%                     ei = this.acqusitionFunc.acquire(objPredictors_sampled{i}, xxx);
%                     pof = this.pofCalculator.getPof(conPredictors_sampled{i}.cell_predictors{j}, xxx);
%                     val = ei.*pof;
%                     imagesc(range, range, reshape(val, nBins, nBins));
%                 end
%                 
%             end
            
        end

        
        function [eix] = eiWrapper(this, arg, objPredictor_sampled, conPredictor_sampled, f_opt, hasFeasibleValue)

            pofs = this.pofCalculator.getPof( conPredictor_sampled, arg) ;
            if(hasFeasibleValue)
                this.acqusitionFunc.updateFOpt(f_opt)
                [ei, ~, ~] =   this.acqusitionFunc.acquire(objPredictor_sampled, arg, []);
                eix = ei .* pofs;
            else
                [mu,~,~] = objPredictor_sampled.predict(arg);
                eix = (this.penalty - mu) .* pofs;
            end
%             [ei, mu, sig] =   this.acqusitionFunc.acquire(predictor, arg, iter);
%             wei = ei.*this.pofCalculator.getPof( multiConstrPredictors, arg) ;
            
        end
        
        function [nei, mu, sig] = neiWrapper(this, arg, objPredictors_sampled, conPredictors_sampled, f_opts, hasFeasibleValue)

            nei = zeros(size(arg,1),1);
            for i = 1:this.quasi_sample_size
                eix = eiWrapper(this, arg, objPredictors_sampled{i}, conPredictors_sampled{i}, f_opts(i), hasFeasibleValue(i));
                nei = nei + eix./this.quasi_sample_size;
            end
            mu = 0;
            sig = 0;
            
%             [ei, mu, sig] =   this.acqusitionFunc.acquire(predictor, arg, iter);
%             wei = ei.*this.pofCalculator.getPof( multiConstrPredictors, arg) ;
            
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
            [objPredictors_sampled, conPredictors_sampled, f_opts, hasFeasibleValue] = sampleModels(this, predictor, multiConstrPredictors);
            
%             acquisition =  @(arg) (this.acqusitionFunc.acquire(predictor, arg, iter).* this.pofCalculator.getPof( multiConstrPredictors, arg) );
            acquisition =  @(arg) neiWrapper(this, arg, objPredictors_sampled, conPredictors_sampled, f_opts, hasFeasibleValue);
    %         [funcHs, learnedHPs] = gp_direct(x, y, cumul_samplesizes, gpHPs);
            % First maximise the MF-GP-UCB acquisition function.
%             acquisition = @(arg) acqGP_UCB_min(arg, funcHs, iter, gpHPs.zetas);
            % acquisition returns: [acq, uncerts, mu, sigma]
%             if(this.isMaximizeAcq)
                [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq, queryVals_mu, queryVals_sigma] = ...
                DiRectWrap_max_module_nargout_3( acquisition, this.xDomain, this.optParam);
            
            saveX(this, nextPt);
%                 [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq] = ...
%                 diRectWrapMax_module( acquisition, this.xDomain, this.optParam);
%             else
%                 [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq, queryVals_mu, queryVals_sigma] = ...
%                 DiRectWrap_min_module_nargout_3( acquisition, this.xDomain, this.optParam);
% %                 [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq] = ...
% %                 diRectWrapMin_module( acquisition, this.xDomain, this.optParam);
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
            optMuV = predictor.predict(optMuX);

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
        
        function out = getNormalSample(this, tt, mu, cov)
            try
                A  = chol(cov,'lower');
            catch 
                try
                    if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget 1e-15.*mean.\n');end;
                    A  = chol(cov + eye(size(cov,1)).*1e-15.*mean(diag(cov),'all'),'lower');
                catch 
                    try
                        if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget 1e-10.*mean.\n');end;
                        A  = chol(cov + eye(size(cov,1)).*1e-10.*mean(diag(cov),'all'),'lower');
                    catch 
                        try
                            if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget 1e-5.*mean.\n');end;
                            A  = chol(cov + eye(size(cov,1)).*1e-5.*mean(diag(cov),'all'),'lower');
                        catch 
                            try
                                if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget 1e-3.*mean.\n');end;
                                A = chol(cov + eye(size(cov,1)).*1e-3.*mean(diag(cov),'all'),'lower');
                            catch 
                                try
                                    if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget 1e-2.*mean.\n');end;
                                    A = chol(cov + eye(size(cov,1)).*1e-2.*mean(diag(cov),'all'),'lower');
                                catch 
                                    try
                                        if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget 1e-1.*mean.\n');end;
                                        A = chol(cov + eye(size(cov,1)).*1e-1.*mean(diag(cov),'all'),'lower');
                                    catch err
                                        try
                                            if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget mean.\n');end;
                                            A = chol(cov + eye(size(cov,1)).*mean(diag(cov),'all'),'lower');
                                        catch err
                                            fprintf('Failed in Cholesky Decomposition.\n');
                                            showErrors(err);
                                            error('cov is not positive definite (InfillOptimizerDirectConstNEI).');
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                if(this.param_BO.debug);fprintf('Succeed.\n');end
            end
            z = norminv(tt);
            try
                [out] = A*z + mu;            
            catch err
                showErrors(err);
            end
            
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