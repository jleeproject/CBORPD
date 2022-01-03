classdef InfillOptimizerDirectConstRPDWEI < AbsInfillOptimizer
    
    
    properties
        name = 'Quadratic(WEI)';
        optParam
        type = TypeInfillOptimizer.DirectConstNEI;
        pofCalculator
        isFirstRun = true;
        feasibilitySamplingConstraint
%         quasi_sample_size = 5;
        quasi_sample_size = 10e3;
%         quasi_sample_size = 1000;
        penalty = 1000;
        
        history_x
        nEvals = 0;
        gpEstimatorObj 
        gpEstimatorCon
        
        acqEI 
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
            
            
            this.acqEI = AcqFuncEI(this.param_BO, this.objFunc, false);

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
        
        
        
        function [f_at_trains, c_at_trains, f_opts, hasFeasibleValue] = sampleModels(this, predictor, multiConstrPredictors)
            
            [ xx_train, n_xxs] = getHistoryX(this);
             
            c_at_trains = zeros(n_xxs, this.quasi_sample_size, multiConstrPredictors.nPredictor);
            c_vec_at_trains = zeros(n_xxs*this.quasi_sample_size, multiConstrPredictors.nPredictor);
                
%             for i=1:n_xxs
            f_at_trains = getNormalSample(this, predictor, xx_train, n_xxs); % dim [n_xxs, this.quasi_sample_size]
            f_vec_at_trains = reshape(f_at_trains,[],1);

            for j=1:multiConstrPredictors.nPredictor
                c_at_trains(:,:,j) = getNormalSample(this, multiConstrPredictors.cell_predictors{j} , xx_train, n_xxs);
                c_vec_at_trains(:,j) = reshape(c_at_trains(:,:,j),[],1);
            end
%             end

            feas_vec = this.feasibilitySamplingConstraint.isThisValueFeasible(c_vec_at_trains);
            if(sum(feas_vec)>0)
                % f_opts
                f_feas_vec = f_vec_at_trains.*(feas_vec)+ max(f_vec_at_trains).*(1-feas_vec); % feasible point has original values & infeasible Maximum value; this will be minimized.
                f_feas_to_be_minimized = reshape(f_feas_vec, n_xxs, this.quasi_sample_size)'; % dim [this.quasi_sample_size, n_xxs]
                f_opts = min(f_feas_to_be_minimized, [], 2); % minimum among training
                
                % hasFeasibleValue
                feas = reshape(feas_vec, n_xxs, this.quasi_sample_size)'; % dim [this.quasi_sample_size, n_xxs]
                hasFeasibleValue = sum(feas,2)>1; % if at least one of the training points is feasible at each sobol sequence.
            else
                hasFeasibleValue = zeros(this.quasi_sample_size,1);
                f_opts = max(f_vec_at_trains).* ones(this.quasi_sample_size, 1);
            end
            
%             if(size(f_opts,1)~= this.quasi_sample_size)
%                 fprintf(2,'  \n');
% %                 fprintf(2,'  \n');
%             end

            
        end

        
        function [nei, mus, sigs] = neiWrapper(this, xx, objPredictor, multiConstrPredictors, f_at_trains, c_at_trains, f_opts, hasFeasibleValue)

            % matrix init            
            n_xxs = size(xx,1);
            nei = zeros(n_xxs , 1);
%             % jslee
%             c_at_test = zeros(n_xxs, this.quasi_sample_size, multiConstrPredictors.nPredictor);
%             c_vec_at_test = zeros(n_xxs*this.quasi_sample_size, multiConstrPredictors.nPredictor);
            
            % sampling at x_test
%             f_at_test = getNormalSample(this, objPredictor, xx, n_xxs); % dim: [n_xxs, this.quasi_sample_size]
%             mu = mean(f_at_test,2);
%             sig = std(f_at_test,[],2);
            

%             %jslee
%             for j=1:multiConstrPredictors.nPredictor
%                 c_at_test(:,:,j) = getNormalSample(this, multiConstrPredictors.cell_predictors{j} , xx, n_xxs); % dim: [n_xxs, this.quasi_sample_size, multiConstrPredictors.nPredictor]
%             end
            
%             % Version 1: EI & FOP are Independent. Each row is calculated.
%             for j=1:n_xxs
% %                 f_at_test; %% dim [ n_xxs , quasi_sample_size ]
%                 f_obj = f_at_test(j,:)'; %% f_obj: dim [ quasi_sample_size , 1 ]
%                 feas = this.feasibilitySamplingConstraint.isThisValueFeasible(  reshape(c_at_test(j,:,:), this.quasi_sample_size, multiConstrPredictors.nPredictor )  );
%                 % f_obj: dim [ quasi_sample_size , 1 ]
%                 
%                 try
%                 ei = hasFeasibleValue.* (f_opts - f_obj) + (1-hasFeasibleValue).*(this.penalty - mu);
%                 catch err
%                     showErrors(err)
%                 end
%                 nei(j) = mean(ei).* mean(feas); % NEI = EI(x|f)p(f) *POF
% %                 nei(j) = (sum(ei)./this.quasi_sample_size) .* (sum(feas)./this.quasi_sample_size); % NEI = EI(x|f)p(f) *POF
% 
%                 % Version 2: EI & FOP are Dependent.
% %                 if(sum(feas)>0)
% %                     idx_feas = find(feas);
% %                     f_feas = f_obj(idx_feas);
% %                     f_opts_feas = f_opts(idx_feas); % For each value in the sobol sequence, optimal point is calculated.
% %                     if(hasFeasibleValue )
% %                         sum_ei_feas = sum(f_opts_feas - f_feas ); % EIx = E( f_opt - f(x) | f(x) ~GP )PoF(x) = (f_opt - f)/n (#feas(x)/n)
% %                         % EI = E( f_opt - f(x) | f(x) ~GP ) = 1/n (f_opt - f) where f is samples in n-rows column vector.
% %                     else
% %                         [mu, ~, ~] = objPredictor.predict(xx(j,:));
% %                         sum_ei_feas = (this.penalty - mu).*sum(feas); % EI 
% %                     end
% %                 else
% %                     sum_ei_feas = 0;
% %                 end
% %                 nei(j) = sum_ei_feas./this.quasi_sample_size;
%                 % NEI = EI(x|f)P(c is feasible)
%             end
            
            % Final Version: Matrix Operation
%             %jslee
%             feas = zeros(n_xxs,this.quasi_sample_size);
%             for j=1:n_xxs
%                 feas(j,:) = this.feasibilitySamplingConstraint.isThisValueFeasible(  reshape(c_at_test(j,:,:), this.quasi_sample_size, multiConstrPredictors.nPredictor )  )';
% 
%             end

            mus = zeros(n_xxs, 1);
            sigs = zeros(n_xxs, 1);
            for j =1:n_xxs
%                 mat_hasFeasibleValue = (ones(n_xxs,1).*hasFeasibleValue');
                ei = zeros(this.quasi_sample_size, 1);
                try
    %                 impv = max( mat_hasFeasibleValue.*( ones(n_xxs,1).*f_opts' - f_at_test ) + (1-mat_hasFeasibleValue).*(this.penalty - mu), 0);
                    idx_feas = find(hasFeasibleValue);
                    idx_infeas = find(1-hasFeasibleValue);
                    ei(idx_feas) = this.acqEI.acquire( objPredictor, xx(j,:), 1, f_opts(idx_feas,:));
                    
                    [mus(j,1), sigs(j,1)] = objPredictor.predict(xx(j,:));
                    ei(idx_infeas) = this.penalty - mus(j,1);

                catch err
                    fprintf(2, 'Size(mat_hasFeasibleValue) = [%d, %d]\n', size(mat_hasFeasibleValue,1), size(mat_hasFeasibleValue,2))
                    fprintf(2, 'Size(f_opts) = [%d, %d]\n', size(f_opts,1), size(f_opts,2))
                    fprintf(2, 'Size(f_at_test) = [%d, %d]\n', size(f_at_test,1), size(f_at_test,2))
                    fprintf(2, 'Size(this.penalty) = [%d, %d]\n', size(this.penalty,1), size(this.penalty,2))
                    fprintf(2, 'Size(mu) = [%d, %d]\n', size(mu,1), size(mu,2))
                    showErrors(err);
                end

                pof = this.pofCalculator.getPof(multiConstrPredictors, xx(j,:));
    %             %jslee
    %             nei = mean(impv,2) .*  mean(feas,2)  ;
    %             nei = mean(impv,2).*pof;
                nei(j) = mean(ei).*pof;
            end
            
%             %% DEBUG TESTING
%             % COMPARE:
%             disp(this.acqEI.acquire( objPredictor, xx, 1, 0));
%             f_opts_test = 0; current_ei = mean(max( f_opts_test' - f_at_test,0));disp(current_ei);
%             acqEI = AcqFuncEI(this.param_BO, this.objFunc, false);acqEI.updateFOpt(0);test_ei=acqEI.acquire(objPredictor, xx, 1);disp(test_ei);
        end
        
        
        function [optAcqX, optAcqV, idxOptAcq, histories, optMuX, optMuV, acquisition, funcFeas, problem, urnd ] = find(this, predictor, iter, multiConstrPredictors)
            
            [f_at_trains, c_at_trains, f_opts, hasFeasibleValue] = sampleModels(this, predictor, multiConstrPredictors);
            acquisition = @(xx) neiWrapper(this, xx, predictor, multiConstrPredictors, f_at_trains, c_at_trains, f_opts, hasFeasibleValue);
            
            % old
%             [objPredictors_sampled, conPredictors_sampled, f_opts, hasFeasibleValue] = sampleModels(this, predictor, multiConstrPredictors);
%             acquisition =  @(arg) neiWrapper(this, arg, objPredictors_sampled, conPredictors_sampled, f_opts, hasFeasibleValue);
            [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq, queryVals_mu, queryVals_sigma] = ...
                DiRectWrap_max_module_nargout_3( acquisition, this.xDomain, this.optParam);
            
            saveX(this, nextPt);

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

            optMuX = optAcqX;
            optMuV = predictor.predict(optMuX);


            this.acqusitionFunc.updateFOpt(optMuV);
        end
        
        function optParam = gpDirectInit(this)
            optParam.maxevals = this.param_BO.maxOptEval;
            optParam.maxits = this.param_BO.maxOptIter;
        end
        
        
        
        
        function out = getNormalSample(this, predictor, xx,  n_xxs)
             
            tt = net(scramble(sobolset(n_xxs),'MatousekAffineOwen'),this.quasi_sample_size)';
            out = zeros(n_xxs,this.quasi_sample_size);
            % matrix: dim [n_xxs, quasi_size] 
            for i=1:n_xxs
                [A, mu] = predictor.getCholDecSig(xx(i,:));
                z = norminv(tt(i,:));
                try
                    out(i,:) = A*z + mu;           % dim [1 , this.quasi_sample_size]  
                catch err
                    showErrors(err);
                end
            end
            
        end

    end
end        
%         function out = getNormalSample(this, tt, mu, cov)
%             try
%                 A  = chol(cov,'lower');
%             catch 
%                 try
%                     if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget 1e-15.*mean.\n');end;
%                     A  = chol(cov + eye(size(cov,1)).*1e-15.*mean(diag(cov),'all'),'lower');
%                 catch 
%                     try
%                         if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget 1e-10.*mean.\n');end;
%                         A  = chol(cov + eye(size(cov,1)).*1e-10.*mean(diag(cov),'all'),'lower');
%                     catch 
%                         try
%                             if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget 1e-5.*mean.\n');end;
%                             A  = chol(cov + eye(size(cov,1)).*1e-5.*mean(diag(cov),'all'),'lower');
%                         catch 
%                             try
%                                 if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget 1e-3.*mean.\n');end;
%                                 A = chol(cov + eye(size(cov,1)).*1e-3.*mean(diag(cov),'all'),'lower');
%                             catch 
%                                 try
%                                     if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget 1e-2.*mean.\n');end;
%                                     A = chol(cov + eye(size(cov,1)).*1e-2.*mean(diag(cov),'all'),'lower');
%                                 catch 
%                                     try
%                                         if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget 1e-1.*mean.\n');end;
%                                         A = chol(cov + eye(size(cov,1)).*1e-1.*mean(diag(cov),'all'),'lower');
%                                     catch err
%                                         try
%                                             if(this.param_BO.debug);fprintf('[Warning] cov is not positive definite. Trying adding nugget mean.\n');end;
%                                             A = chol(cov + eye(size(cov,1)).*mean(diag(cov),'all'),'lower');
%                                         catch err
%                                             fprintf('Failed in Cholesky Decomposition.\n');
%                                             showErrors(err);
%                                             error('cov is not positive definite (InfillOptimizerDirectConstNEI).');
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%                 if(this.param_BO.debug);fprintf('Succeed.\n');end
%             end
%             z = norminv(tt);
%             try
%                 [out] = A*z + mu;            
%             catch err
%                 showErrors(err);
%             end
%             
%         end
% 
%     end
% end

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