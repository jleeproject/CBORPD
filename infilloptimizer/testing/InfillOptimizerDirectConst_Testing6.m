classdef InfillOptimizerDirectConst_Testing6 < AbsInfillOptimizer
    %% Testing class. Will be transfer to a new class
    
    properties
        name = 'Exploit(mu)_f_opt_before';
        optParam
        type = TypeInfillOptimizer.DirectConstTesting9;
        pofCalculator
        isFirstRun = true;
        feasibilitySamplingConstraint
%         optStep = true;
        optStep = false;
%         strStep = 'feas';
        
        acqEI
        acqUCB
        acqLCB
        acqPOI
        cnt = 1;
        reqCnt = 1000;
        penalty = 1000;
        
        
        history_x
        nEvals = 0;
    end
    
%     methods (Abstract)
%         method2(this)
%     end
    methods
        function init(this)
            this.acqFuncName = 'Proposed:Exploit(mu)_f_opt_before';
%             this.gpHPs.diRectParams = 0
            this.optParam = gpDirectInit(this);
            
%             if(~isa(this.acqEI,'AcqFuncEI'))
%                 fprintf('[WARNING] EI should be used for the acquisition function. (InfillOptimizerDirectConstWEI)\n');
%                 fprintf('[WARNING] Replace to EI\n');
%                 prevAcqFunc = this.acqEI;
                this.acqEI = AcqFuncEI(this.param_BO, this.objFunc, false);
                this.acqLCB = AcqFuncLCB(this.param_BO, this.objFunc, false);
                this.acqUCB = AcqFuncUCB(this.param_BO, this.objFunc, true);
                this.acqPOI = AcqFuncPOI(this.param_BO, this.objFunc, false);
                this.isMaximizeAcq = true;
%             end
            this.pofCalculator = PoFCalculator(this.param_BO);
            this.feasibilitySamplingConstraint = FeasibilitySamplingConstraint(this.param_BO);
            %--------------------------
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
        
%         function [acq, mu, sig, poi_obj, pof_con] = objoptWrapper(this, objPredictor, iter, multiConstrPredictors, arg)
        function [acq, mu, sig, empty, empty_] = objoptWrapper(this, objPredictor, iter, multiConstrPredictors, arg)
%             [poi_obj, mu, sig] =   this.acqPOI.acquire(objPredictor, arg, iter);
            empty = 0;
            empty_ = 0;
            
            [mu, sig] = objPredictor.predict(arg);
%             pof_con = this.pofCalculator.getPof( multiConstrPredictors, arg) ;
            [mu_con, sigma_con] = multiConstrPredictors.predict(arg);
                
%             % Expected POI
%             acq = poi_obj .* pof_con;
            % New acquisition for optimal
            feasgap_test = this.feasibilitySamplingConstraint.feasGap(mu_con);
            acq = mu + this.penalty .* feasgap_test;
%             pof_con = this.pofCalculator.getPof( multiConstrPredictors, arg) ;
%                 
%             [mu, sig] = objPredictor.predict(arg);
%             f_opt = this.acqPOI.f_opt;
%             poi_obj = normcdf( (f_opt - mu +(1-pof_con).*1*abs(f_opt) )./sig  );
%             % Expected POI
%             acq = poi_obj .* pof_con;
        end
        
        function [acq, mu, sig, poi_obj, pof_con, sigma_con] = feasWrapper(this, objPredictor, iter, multiConstrPredictors, arg)
                %% Candidate 1.
% %             [ei, mu, sig] =   this.acqEI.acquire(objPredictor, arg, iter);
% %             wei = ei.*this.pofCalculator.getPof( multiConstrPredictors, arg) ;
            [poi_obj, mu, sig] =   this.acqPOI.acquire(objPredictor, arg, iter);
%             
%             pof_con = this.pofCalculator.getPof( multiConstrPredictors, arg, -0.01) ;
            pof_con = this.pofCalculator.getPof( multiConstrPredictors, arg) ;
            
            [mu_con, sigma_con] = multiConstrPredictors.predict(arg);
%             pdf_con = normpdf(-mu_con./sigma_con);
            %% Original way
%             acq = poi_obj .* pof_con.* sigma_con;
%             acq = poi_obj .* pof_con.* prod(sigma_con,2).^(1./size(sigma_con,2))  - this.penalty .*(1-pof_con) ;
            acq = poi_obj .* pof_con.* prod(sigma_con,2).^(1./size(sigma_con,2))  ;

%             %% NEW WAY
%             [mu, sig] = objPredictor.predict(arg);
%             f_opt = this.acqPOI.f_opt;
%             poi_obj = normcdf( (f_opt - mu + (1-pof_con).*1*abs(f_opt) )./sig  );
%             acq = poi_obj .* pof_con.* sigma_con;
            %             
        end
        
        
        function [optAcqX, optAcqV, idxOptAcq, histories, optMuX, optMuV, acquisition, funcFeas, problem, urnd ] = find(this, objPredictor, iter, multiConstrPredictors)

            acquisition_optsol =  @(arg) objoptWrapper(this, objPredictor, iter, multiConstrPredictors, arg);
            
            [ xx_train, ~] = getHistoryX(this);
            v_train = acquisition_optsol(xx_train);
            [~, idx_max ] = max(v_train);
            v_min = objPredictor.predict(xx_train(idx_max,:));
            this.acqusitionFunc.updateFOpt(v_min);
            this.acqPOI.updateFOpt(v_min);
            this.acqEI.updateFOpt(v_min);
            
            
            %% feasibility step
                %% Candidate 1.
                acquisition =  @(arg) feasWrapper(this, objPredictor, iter, multiConstrPredictors, arg);
                [nextPtAcq, nextPt, history_imprv, queries_all, queryVals_acq, queryVals_mu, queryVals_sigma] = ...
                DiRectWrap_max_module_nargout_3( acquisition, this.xDomain, this.optParam);
%                 DiRectWrap_min_module_nargout_3( acquisition, this.objFunc.getXDomain(), this.optParam);
                %% Candidate 2.
                
                
                this.cnt = this.cnt + 1;
                if(this.cnt >= this.reqCnt)
                    this.optStep = true;
                    this.cnt = 1;
                    this.strStep = 'opt';
                end
%             end
            

%                 acq = EI(obj).*PoF(con);
                [optMuV, optMuX, history_imprv, queries_all, queryVals_acq, queryVals_mu, queryVals_sigma] = ...
                DiRectWrap_min_module_nargout_3( acquisition_optsol, this.xDomain, this.optParam);
%                 DiRectWrap_max_module_nargout_3( acquisition_optsol, this.xDomain, this.optParam);
%                 DiRectWrap_min_module_nargout_3( acquisition, this.objFunc.getXDomain(), this.optParam);

% 
%                 pof_fn =  @(arg) this.pofCalculator.getPof( multiConstrPredictors, arg);
%                 [optpofV, optpofX] =  diRectWrapMax_module( pof_fn, this.xDomain, this.optParam);

%                 
            if(this.param_BO.debug)
                [acq, mu, sig, poi_obj, pof_con, sigma_con] = acquisition(nextPt);
                fprintf('[infill feas] poi = %g \t pof = %g \t sig = %.2g \t (acq = %.2g)\n', poi_obj, pof_con, sigma_con, acq);
                
                [acq, mu, sig, poi_obj, pof_con] =  acquisition_optsol(optMuX);
                fprintf('[infill opt] poi = %g \t pof = %g (acq = %.2g)\n', poi_obj, pof_con, acq);
                
                
                acquisition_pof =  @(arg) this.pofCalculator.getPof( multiConstrPredictors, arg);
                [~, mostFeasPoint, ~, ~] = ...
                diRectWrapMax_module( acquisition_pof, this.xDomain, this.optParam);
                val = acquisition_pof(mostFeasPoint);
                fprintf('[most feasible P(C(x))=%g]\n',val);
            

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

%             optMuX = optAcqX;
            optMuV = objPredictor.predict(optMuX);

%             %% Not implement maximizing mu. 
%             this.acqusitionFunc.updateFOpt(optMuV);
%             this.acqPOI.updateFOpt(optMuV);
%             this.acqEI.updateFOpt(optMuV);
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