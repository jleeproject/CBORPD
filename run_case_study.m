resume = false;
% resume = true;
if(resume); if (~askYesNo('Are you sure you want RESUME?'));return;end;end
clc;fprintf('BO starts at %s\n', showPrettyDateTime(now()));setVersion; fprintf('%s\n',versionBO);if(~resume);clear all; resume=false;end

scr_run_on_server
% scr_run_on_my_pc

scr_Casestudy;

scr_common_v2;

stt_all = tic();cnt_exp = 1;; time_exp_start = now;
elapsedTime_settings = zeros(nTypeSettings,1);
tot_exp = numel(cell_obj_func)* nRepeats* nTypeSettings* numel(type_cell_samplesize);
stt_periodic = tic();
for repeat = 1:nRepeats;
    stt_repeat = tic();
    for idx_function = 1:numel(cell_obj_func)
        scr_setting_1_in_functions_sameObj;

        stt_function = tic();
        warning('off');
        if(~resume);scr_setting_2_in_repeats_same_initPts;end;
        for idx_setting = 1:nTypeSettings
            stt_setting = tic;
            scr_setting_3_in_settings; 
    
            for idx_samplesize = 1:numel(type_cell_samplesize)
                fprintf('[Exp] func:%d/%d, rep:%d/%d, set:%d/%d, samp:%d/%d (%d/%d)\n', ...
                    idx_function, numel(cell_obj_func), repeat, nRepeats, idx_setting, nTypeSettings,...
                    idx_samplesize, numel(type_cell_samplesize),...
                    cnt_exp, tot_exp);
                
                scr_setting_4_in_samples;
                
                scr_set_settings_pre_run;
                
                stt_experiment = tic();
                
                if(param_BO.isInfo()); show_BO_parameters_module; end
                
                if(setting.type_infill_opt == TypeInfillOptimizer.DirectConstAdmm)
                    is_original_admm = true;
                    scr_external_admmbo;
                    scr_set_admmbo_results
                elseif(setting.type_infill_opt == TypeInfillOptimizer.DirectConstAdmmNew)
                    is_original_admm = false;
                    scr_external_admmbo;
                    scr_set_admmbo_results
                else
                    niter = param_BO.getEstUpperLimitNumIter(typeSamplesize, givenSamplesize);

%                     if(param_BO.isInfo()); show_BO_parameters_module; end

                    try
                        if(setting.type_infill_opt == TypeInfillOptimizer.DirectConstPatternSearch)
                            scr_pattern_search;
                        else
                            for iter = 1:niter;
                                sttIter = tic();
                                scr_bayesian_optimization

                                if(param_BO.debug)
                                    finIter = toc(sttIter);
    %                                 fprintf('[BO: Eval = %d/%d (Iter = %d/%d)] :%s\n', nEvaluations, budget, iter, niter, showPrettyElapsedTime(finIter));
                                    fprintf('[%s: Eval = %d/%d (Iter = %d/%d)], iter time: %s\n', showPrettyDateTime(now), nEvaluations, budget, iter, niter, showPrettyElapsedTime(finIter));
                                elseif(param_BO.info)
                                    if(toc(stt_periodic)>print_interval) %% Every 60 second;
                                        estRemainedT = getEstimatedRemainedTime(toc(stt_experiment), iter, niter);
                                        fprintf('[%s: Eval = %d/%d (Iter = %d/%d)], Elapsed Time (This Exp.): %s, Est. Remained T: %s\n', showPrettyDateTime(now), nEvaluations, budget, iter, niter, showPrettyElapsedTime(toc(stt_experiment)), showPrettyElapsedTime(estRemainedT) );
                                        stt_periodic = tic();
                                    end
                                end
    %                             pause();
                                if(stop)
                                    break;
                                else
                                    if(param_BO.isPlotOnline())
                                        sttVisual = tic ;
                                        if( setting.type_infill_opt ~= TypeInfillOptimizer.DirectConstSobol)
                                            scr_visualize_BO_online_multid_v4
                                            if( ~isa(optProb.objFunc,'GrapheneModelSolver') && ~setting.show_only_performance )
                                                scr_test_subplot;
                                            end
                                            if(param_BO.info)
                                                finVinsual = toc(sttVisual);
                                                fprintf('[Visual] :%s\n', showPrettyElapsedTime(finVinsual));
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    catch err
                        nEvaluations = nEvaluations-1;
                        showErrors(err)
%                         disp(err)
                        error(err.message);
                    end
                    scr_set_results;
                end
                fprintf('> Single Experiment Time: %s\n', showPrettyElapsedTime(toc(stt_experiment)));
%                 fprintf('Expected End Time:%s\n', showPrettyDateTime(now + seconds(toc(stt_all)/cnt_exp* (tot_exp-cnt_exp))))
                showEstimatedEndRemainedTime(toc(stt_all), cnt_exp, tot_exp);
                saveresults_intermediate;
                fprintf('==========================================================\n');
                cnt_exp = cnt_exp+1;
            end
            elapsedTime_settings(idx_setting)=elapsedTime_settings(idx_setting)+ toc(stt_setting);
        end
        fin_func = toc(stt_function);
        fprintf('Total time for a function ''%s'' : %s\n',setting.optProb.objFunc.name, showPrettyElapsedTime(fin_func));
    end
    fin_rep = toc(stt_repeat);
    fprintf('Total time for a repeat ''%s'' : %s\n',setting.optProb.objFunc.name, showPrettyElapsedTime(fin_rep));
    saveresults_intermediate;
end
fin_all = toc(stt_all);
fprintf('Total time: %s\n', showPrettyElapsedTime(fin_all));
fprintf('Current Time: %s\n',showPrettyDateTime(now()));

fprintf('Elapsed time for each setting: \n');
fprintf(' > %s\n', showPrettyElapsedTime( elapsedTime_settings./nRepeats./numel(cell_obj_func)./numel(type_cell_samplesize) ));
save('result_itermediate');
fprintf('%s\n',A_RUN_PURPOSE)
Visualize_results_figures
