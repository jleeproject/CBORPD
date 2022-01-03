nEvaluations = nEvaluations+1;
if(param_BO.isDebug())
    if(mod(iter,10)==1)
        fprintf('iter: %d/%d\n', iter, niter);
        fprintf('# Samples: %d/%d\n', sum(storage.getSamplesizeHistoryTill(nEvaluations)), budget );
    end
end
%% ------------------------------------------------------------------------
%% ----- Obj Func : GP Fit Estimation ----------------------------------------------------
    x_train = storage.getXHistoryTill(nEvaluations-1);
    y_train = storage.getYHistoryTill(nEvaluations-1);
    cumul_samplesizes = storage.getSamplesizeHistoryTill(nEvaluations-1);
    [infillPredictor, param, exactPredictor] = estimator.estimateGp( x_train, y_train, cumul_samplesizes);
    if(param_BO.hasConstraints)
        y_con_train = storage.getYConsHistoryTill(nEvaluations-1);
        [multiPredictors, gpParam, multiTruePredictors] = estimatorCon.estimateGp(x_train, y_con_train, cumul_samplesizes, @exactPredictor.predict);
    end
%% ----- Constraints : GP Fit Estimation ----------------------------------------    
    % [TODO] Gp Fit. Ideally GP with different GP settings
    % May combine with infill predictor

%% ----- SELECT INFILL SAMPLE ----------------------------------------------------
    % [TODO]: infilloptimizer will deal with Constrainted BO.
if(param_BO.hasConstraints)
    [xnew, optAcqV, idxOptAcq, queries_hist, optMuX, optMuV, acquisition, funcFeas, problem, urnd_feasiblity] = infillOptimizer.find(infillPredictor, iter, multiPredictors);
else
    [xnew, optAcqV, idxOptAcq, queries_hist, optMuX, optMuV, acquisition] = infillOptimizer.find(infillPredictor, iter);
    problem = struct();
end

if(setting.typeProblem == TypeProblem.RobustDesignOptimization)
%% ----- Determine Sample Size ----------------------------------------------------
    samplesize = samplesizeSelector.selectSamplesize(xnew, infillPredictor);
    if(sum(cumul_samplesizes) + samplesize>=budget)
        samplesize = budget - sum(cumul_samplesizes);
        stop = true;
    end
else
    samplesize = 1;
%     if(n_init_sample + nEvaluations>=budget)
    if(nEvaluations>=budget)
%         samplesize = sum(cumul_samplesizes) - budget;
        stop = true;
    end
end

if(~stop)
    %% ----- EVALUATE NEW POINT ----------------------------------------------------
    % ynew = simulator.evaluate(xnew, samplesize); % Only sample if samplesize >1
    sttEval = tic();
    [ynew, y_cons, y_orig] = MultipleSamplerObjCon.evaluate(simulator, xnew, samplesize);
    elapsedEval = toc(sttEval);
    
    storage.addElapsedEval(elapsedEval);
    elapsedIter = toc(sttIter);
    storage.addElapsedExp(elapsedIter);
    storage.addElapsedBoWoEval(elapsedIter - elapsedEval);

    if(param_BO.typeProblem == TypeProblem.ConventionalBoDeterministic)
        infillOptimizer.acqusitionFunc.updateFOpt( ynew );
    end

    if(numel(ynew)>0)
        % --- Save Intermediate Results ----
        %
        storage.saveXHistory(xnew, nEvaluations);
        storage.saveYHistory(ynew, nEvaluations);
        storage.saveObsHistory(y_orig, nEvaluations);
        if(param_BO.hasConstraints)
            storage.saveYConsHistory(y_cons, nEvaluations);
        end
        storage.saveSamplesizeHistory(samplesize, nEvaluations);
        % True value at xnew

    %     if(param_BO.getTypeGpFit == TypeEstimatorGp.DirectLogSampleVarRandomLS)
            % New definition of optimal value.
    %     if(wantDetailedOptimalValues && ~type_infill_opt.strcmp('DirectConstWEI'))
%         if(wantDetailedOptimalValues && ~type_infill_opt.strcmp('DirectConstWEI') && ~type_infill_opt.strcmp('DirectConstNEI') && ~contains((type_infill_opt.char),'Testing'))
% 
%             if(param_BO.hasConstraints)
% 
%                 problem.constraint.func = @(arg) funcFeas(arg, .95);
%                 if(setting.optProb.objFunc.isMaximize)
%                     [f_opt, x_optMu, history_imprv, queries_all, queryVals_mu] = ...
%                     diRectWrapMax_given_constraint( @exactPredictor.predict, setting.optProb.objFunc.getXDomain() , gpSettingsObj.gpParam.diRectParams, problem);
%                 else
%                     [f_opt, x_optMu, history_imprv, queries_all, queryVals_mu] = ...
%                     diRectWrapMin_given_constraint( @exactPredictor.predict, setting.optProb.objFunc.getXDomain() , gpSettingsObj.gpParam.diRectParams, problem);
%                 end
% 
%                 storage.saveXMuMinHistory(x_optMu,nEvaluations);
%                 storage.saveYMuMinHistory(f_opt,nEvaluations);
%                 storage.saveEvalFAtXnewHistory( (setting.optProb.objFunc.evalWithVecX(x_optMu)) , nEvaluations)
%                 problem.constraint.func = @(arg) funcFeas(arg, .05);
%                 if(setting.optProb.objFunc.isMaximize)
%                     [f_opt, x_optMu, history_imprv, queries_all, queryVals_mu] = ...
%                     diRectWrapMax_given_constraint( @exactPredictor.predict, setting.optProb.objFunc.getXDomain() , gpSettingsObj.gpParam.diRectParams, problem);
%                 else
%                     [f_opt, x_optMu, history_imprv, queries_all, queryVals_mu] = ...
%                     diRectWrapMin_given_constraint( @exactPredictor.predict, setting.optProb.objFunc.getXDomain() , gpSettingsObj.gpParam.diRectParams, problem);
%                 end
%                 if(~storage.hasLiberal); storage.hasLiberal = true;end;
%                 storage.saveXMuMinLiberalHistory(x_optMu,nEvaluations);
%                 storage.saveYMuMinLiberalHistory(f_opt,nEvaluations);
%                 storage.saveEvalFAtXnewLiberalHistory( (setting.optProb.objFunc.evalWithVecX(x_optMu)) , nEvaluations)
% 
%                 feas = 1;
%                 for i = 1:simul.nCon
%                     gval = simul.con{i}.meanFunc.evalWithVecX(x_optMu);
%                     if(optProb.con.constraints{i}.hasLb && optProb.con.constraints{i}.hasUb )
%                         lb = optProb.con.constraints{i}.lb;
%                         ub = optProb.con.constraints{i}.ub;
%                         feas = feas.*( gval<=ub && lb<=gval );
%                     elseif(optProb.con.constraints{i}.hasLb && ~optProb.con.constraints{i}.hasUb )
%                         lb = optProb.con.constraints{i}.lb;
%                         feas = feas.*( lb<=gval ) ;
%                     elseif(~optProb.con.constraints{i}.hasLb && optProb.con.constraints{i}.hasUb )
%                         ub = optProb.con.constraints{i}.ub;
%                         feas = feas.*( gval<=ub ) ;
%                     elseif(~optProb.con.constraints{i}.hasLb && ~optProb.con.constraints{i}.hasUb )
%                         error('Both ub and lb are not specified.');
%                     else
%                     end
%                 end
%                 storage.saveFeasAtXnewHistory( feas , nEvaluations)
%             else
%                 if(setting.optProb.objFunc.isMaximize)
%                     [f_opt, x_optMu, history_imprv, queries_all, queryVals_mu] = ...
%                     diRectWrapMax_module( @exactPredictor.predict, setting.optProb.objFunc.getXDomain() , gpSettingsObj.gpParam.diRectParams);
%                 else
%                     [f_opt, x_optMu, history_imprv, queries_all, queryVals_mu] = ...
%                     diRectWrapMin_module( @exactPredictor.predict, setting.optProb.objFunc.getXDomain() , gpSettingsObj.gpParam.diRectParams);
%                 end
%                 storage.saveXMuMinHistory(x_optMu,nEvaluations);
%                 storage.saveYMuMinHistory(f_opt,nEvaluations);
%                 storage.saveEvalFAtXnewHistory( (setting.optProb.objFunc.evalWithVecX(x_optMu)) , nEvaluations)
%             end
%         else
            storage.saveXMuMinHistory(optMuX,nEvaluations);
            storage.saveYMuMinHistory(optMuV,nEvaluations);
                if( isa(optProb.objFunc,'AbsFunction'))
                    storage.saveEvalFAtXnewHistory( (setting.optProb.objFunc.evalWithVecX(optMuX)) , nEvaluations)
                    feas = 1;
                    feasgap = 0;
                    for i = 1:optProb.nCon
                        gval = optProb.conFunc{i}.evalWithVecX(optMuX);
                        if(optProb.con.constraints{i}.hasLb && optProb.con.constraints{i}.hasUb )
                            lb = optProb.con.constraints{i}.lb;
                            ub = optProb.con.constraints{i}.ub;
                            feas = feas.*( gval<=ub && lb<=gval );
                            feasgap = max( max( [gval - ub , lb - gval, 0]), feasgap);
                        elseif(optProb.con.constraints{i}.hasLb && ~optProb.con.constraints{i}.hasUb )
                            lb = optProb.con.constraints{i}.lb;
                            feas = feas.*( lb<=gval ) ;
                            feasgap = max( max( [lb - gval , 0]), feasgap) ;
                        elseif(~optProb.con.constraints{i}.hasLb && optProb.con.constraints{i}.hasUb )
                            ub = optProb.con.constraints{i}.ub;
                            feas = feas.*( gval<=ub ) ;
                            feasgap = max( max( [gval - ub, 0]), feasgap) ;
                        elseif(~optProb.con.constraints{i}.hasLb && ~optProb.con.constraints{i}.hasUb )
                            error('Both ub and lb are not specified.');
                        else
                        end
                        storage.saveFeasAtXnewHistory( feas , nEvaluations)
                        storage.saveEvalGXnewHistory( gval , nEvaluations, i )
                        storage.saveFeasgapAtXnewHistory( feasgap , nEvaluations, i )
                    end
                    storage.saveFeasAtXnewHistory( feas , nEvaluations)
                    
                    
                elseif( isa(optProb.objFunc,'GrapheneModelSolver'))
                    data_samples_for_validation = simulator.evaluate(optMuX, setting.casestudy.validation.samplesize);
                    samples_for_validation = data_samples_for_validation.samples;
                    storage.saveEvalFAtXnewHistory( std(samples_for_validation,0,'all') , nEvaluations)
                    
                    gval = mean(samples_for_validation,'all');
                    if(optProb.con.constraints{1}.hasLb && optProb.con.constraints{1}.hasUb )
                        lb = optProb.con.constraints{1}.lb;
                        ub = optProb.con.constraints{1}.ub;
                        feas = gval<=ub && lb<=gval ;
                        feasgap = max( [gval - ub , lb - gval, 0]);
                    elseif(optProb.con.constraints{1}.hasLb && ~optProb.con.constraints{1}.hasUb )
                        lb = optProb.con.constraints{1}.lb;
                        feas = lb<=gval ;
                        feasgap = max( [lb - gval , 0]) ;
                    elseif(~optProb.con.constraints{1}.hasLb && optProb.con.constraints{1}.hasUb )
                        ub = optProb.con.constraints{1}.ub;
                        feas = gval<=ub ;
                        feasgap = max( [gval - ub, 0]) ;
                    elseif(~optProb.con.constraints{1}.hasLb && ~optProb.con.constraints{1}.hasUb )
                        error('Both ub and lb are not specified.');
                    else
                    end
                    storage.saveFeasAtXnewHistory( feas , nEvaluations)
                    storage.saveEvalGXnewHistory( gval , nEvaluations, i )
                    storage.saveFeasgapAtXnewHistory( feasgap , nEvaluations, i )
                else
                    error('[AbsInfillOptimizer] Undefined Type.');
                end
%         end
    end
else
%     if(stop)
        nEvaluations = nEvaluations - 1; % The last ynew failed in observation
%     end
end
        
        
%% ------------------------------------------------------------------------
%% ----- VISUALIZATION ----------------------------------------------------
%% ------------------------------------------------------------------------
% if(param_BO.isPlotOnline())
%     scr_visualize_BO_online_multid_v4
% end
% if(param_BO.isDebug())
%     finIter = toc(sttIter);
%     fprintf('[Iter After Evaluation] :%s\n', showPrettyElapsedTime(finIter));
% end

%                 hyp_input = hyp2;
%% ------------------------------------------------------------------------

%% Record time
        storage.saveCumElapsedTime(toc(stt_experiment), nEvaluations);
        storage.saveIterElapsedTime(toc(sttIter), nEvaluations);
%% Termination Condition               
% if(sum(samplesize_history(1:nEvaluations,1))>budget)
if(stop)
%     if(wantDetailedOptimalValues && ~type_infill_opt.strcmp('DirectConstWEI')&& ~type_infill_opt.strcmp('DirectConstNEI') && ~contains((type_infill_opt.char),'Testing'))
%         final_opt_x = x_optMu;
%         final_opt_f_mu = f_opt;
% 
%         final_opt_f_true = (setting.optProb.objFunc.evalWithVecX(x_optMu) );
%     else
        final_opt_x = optMuX;
        final_opt_f_mu = optMuV;

        final_opt_f_true = (setting.optProb.objFunc.evalWithVecX(optMuX) );
%     end
end

% -------------------------------------------------------------------------

