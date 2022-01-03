isFirstRun = true;
% isFirstRun = false;
% summarize = false;
% showIndividual=true;
showIndividual=false;
% showDetail = true;
showDetail = false;
showOnlySummary = true;
show_fmin_diff =false;
clc;
% clearvars -except results_all;
    close all;
    
    
% mysubplot = @(x1,x2,x3) subtightplot(x1,x2,x3, [.06, 0.02]);
mysubplot = @(x1,x2,x3) subplot(x1,x2,x3);
marks = {'r-','b-','c-','g-','k-','m-','y-', 'r--','b--','c--','g--','k--','m--','y--', 'r:','b:','c:','g:','k:','m:','y:'};
marks_liberal = {'r:','b:','c:','g:','k:','m:','y:', 'r-.','b-.','c-.','g-.','k-.','m-.','y-.'};
%     figure(99)
if(isFirstRun)
%     warning('off','all')
%     load save/result_210309
%     load result_210310_power_0_3.mat
%     load save/results_210312_rep100_power1_budget500.mat
%     load save/result_210315_rep_100_power_0_3_budget_1000
%     load save/intermediate_result_210315_rep50_100_budget1000_power1.mat
%     load results_210320_random_lengthscale.mat
%     load results_210321_rep_100_random_lengthscale.mat
%     load save/results_210320_power_0_1_v0_1_0.mat
%     load save/results_2010322_constraints_only_first_exp_meaningful;
%     load results_210323_constraints_trial;
%     load save/result_210317_random_ls;
    [dim_cases_samplesize, dim_repeats, dim_func, num_settings] = size(results_all);
%     dim_func = 1;
%     warning('on','all')
end
% return;
% fig = figure(80);clf;
% fig.Position = [500 0 1500 1000]; % For 3x3
% fig.Position = [100 0 2000 1000]; % For 3 x 4
% fig.Position = [100 100 1500 400]; % For 1 x 3

if(showDetail)
    mysubplot_frame = [dim_cases_samplesize,4];
else
    mysubplot_frame = [dim_cases_samplesize,num_settings*dim_func];
end
idx_subfigure_row = 1;


%%
%%
    idxSumCon1 = 21;
    idxSumCon2 = 22;
    idxSumCon3 = 23;
    idxSumCon4 = 24;
    idxSumUncon1 = 3;
    idxSummary = 4;
%    
    idxIndCon1 = 11;
    idxIndCon2 = 12;
    idxIndCon3 = 13;
    idxIndCon4 = 14;
    idxIndUncon1 = 1;
%     idxIndUncon2 = 2;
% if(~showDetail)
%     fig1 = figure(1);clf;
%     figUncon1 = figure(2);clf;
%     fig3 = figure(3);clf;
%     figCon1 = figure(11);clf;
%     figCon2 = figure(12);clf;
% end
    % visualize_online_BO_2d_acq_search
    
predicted_cell = cell(dim_cases_samplesize, dim_func, num_settings);

dim_all  = zeros(dim_cases_samplesize, dim_func, num_settings);
yys_cell = cell(dim_cases_samplesize, dim_func, num_settings);
xxs_cell = cell(dim_cases_samplesize, dim_func, num_settings);
ypred_cell = cell(dim_cases_samplesize, dim_func, num_settings);

if(param_BO.hasConstraints)
    % xxs_cell{i,k,l}
    feasible_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
    xxsFeasible_cell  = cell(dim_cases_samplesize, dim_func, num_settings);
    optFeasible_cell = cell(dim_cases_samplesize, dim_func, num_settings);
    feasgap_cell     = cell(dim_cases_samplesize, dim_func, num_settings);
    probCons_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
    

    % pred
    feasible_pred_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
    xxsFeasible_pred_cell = cell(dim_cases_samplesize, dim_func, num_settings);
    optFeasible_pred_cell = cell(dim_cases_samplesize, dim_func, num_settings);
    feasgap_pred_cell     = cell(dim_cases_samplesize, dim_func, num_settings);
    probCons_pred_cell    = cell(dim_cases_samplesize, dim_func, num_settings);

%     if(hasLiberal )
            % 2) liberal
        feasible_liberal_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
        xxsFeasible_liberal_cell  = cell(dim_cases_samplesize, dim_func, num_settings);
        optFeasible_liberal_cell = cell(dim_cases_samplesize, dim_func, num_settings);
        feasgap_liberal_cell     = cell(dim_cases_samplesize, dim_func, num_settings);
        probCons_liberal_cell    = cell(dim_cases_samplesize, dim_func, num_settings);

        feasible_liberal_pred_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
        xxsFeasible_liberal_pred_cell = cell(dim_cases_samplesize, dim_func, num_settings);
        optFeasible_liberal_pred_cell = cell(dim_cases_samplesize, dim_func, num_settings);
        feasgap_liberal_pred_cell     = cell(dim_cases_samplesize, dim_func, num_settings);
        probCons_liberal_pred_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
%     end
    
end

minObj_over_funcs = Inf * ones(dim_func,1);
maxObj_over_funcs = -Inf * ones(dim_func,1);
minMeanObj_over_funcs = Inf * ones(dim_func,1);
maxMeanObj_over_funcs = -Inf * ones(dim_func,1);
                

% 
%% Loop: Settings
for l=1:num_settings
    setting = type_cell_settings{l};
    if(iscell(setting))
        str_setting = sprintf('%s,',setting{:} );
    else
        str_setting = sprintf('%s,',setting );
    end
%% Loop: Functions 
    for k=1:dim_func;
%     for k=2:dim_func;
% %% ----------------------------------------------------------------------
        if(showDetail)
            fig1 = figure(50*(l-1)+k);clf;
        end
%% Loop: Sample sizes
        for i=1:dim_cases_samplesize
            idx_pred_stt = 1;            idx_pred_end = 0;
            idx_feas_stt = 1;            idx_feas_end = 0;
            idx_feas_liberal_stt = 1;            idx_feas_liberal_end = 0;
            %% Loop: Repeats
            hasValue = false;
            for j=1:dim_repeats
                result = results_all{i,j,k,l};
%                 result.setting.optProb.name
%                 continue;
%             results_all{idx_samplesize,repeat, idx_function, idx_setting} = struct_result;
                if ([0,0]==size(result))
%                     return;
                    continue;
                end
                hasValue = true;
%     %             mysubplot(dim_cases_samplesize,1,i);
%                 x_history = result.x_history;
%                 y_history = result.y_history;
%                 eval_f_at_xnew_history = result.eval_f_at_xnew_history;
%                 samplesize_history = result.samplesize_history;
%                 yMuMinHistory = result.yMuMinHistory;
%                 xMuMinHistory = result.xMuMinHistory;
%                 hist_beta_t = result.hist_beta_t;
%                 nEvaluations = result.nEvaluations;
%                 n_init_sample = result.n_init_sample;
%                 param_BO = result.param_BO;
%                 visualizer = result.visualizer;
%                 budget = result.budget;
%                 givenSamplesize = result.givenSamplesize;
%                 typeSamplesize = result.typeSamplesize;
%                 elapsedTime = result.elapsedTime;
%                 target_function = result.objFunc;
%                 iterElapsedTime = result.iterElapsedTime;
% %                 scr_true_functions;
%                 objFunc = result.objFunc;
%                 opt_sol_fn = result.opt_sol_fn;
%                 opt_val_fn = result.opt_val_fn;
%                 power = result.power;
                
                scr_setResultValuesToVariables
                
                cum_samplesize =  cumsum(samplesize_history(n_init_sample+1:nEvaluations,1));

                sttX = sum(samplesize_history(1:n_init_sample));
%                 domainX = [sttX:budget]';

                %% Prediction
                idx_pred_end = idx_pred_end + nEvaluations - n_init_sample;
                % ----- Assign values --------------------------------------
                xxs_cell{i,k,l}(idx_pred_stt:idx_pred_end)= ...
                cumsum(samplesize_history(n_init_sample+1:nEvaluations,1))  + sum(samplesize_history(1:n_init_sample));
            
                yys_cell{i,k,l}(idx_pred_stt:idx_pred_end)= ...
                    eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1);

                show_exact_fmin = ~show_fmin_diff;
                if(param_BO.hasConstraints)
                    constraints = param_BO.constraints;
                    % ----- Assign Calculated values --------------------------------------
                    cum_samplesize =  cumsum(samplesize_history(n_init_sample+1:nEvaluations,1));
                    %% 1. feasible
                    [feasible, strCon] =  isFeasible(constraints.constraints{1}, meanFunc.fnEval , xMuMinHistory(n_init_sample+1:nEvaluations,:));
                    feasible_cell{i,k,l}(idx_pred_stt:idx_pred_end)  = feasible;
                    
                    %% 2. optimal among feasible (val(idxFeas))
                    idxFeas = find(feasible);
                    idx_feas_end = idx_feas_end + numel(idxFeas);
                    if(show_exact_fmin)
                        eval_f =   eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1) ;
                        val = eval_f;
                    else
                        diff_f_min = (  eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1)  - log(objFunc.getOptVal() ));
                        val = diff_f_min;
                    end
                    optfeas = val(idxFeas);
                    optFeasible_cell{i,k,l}(idx_feas_stt:idx_feas_end)= optfeas;
                    xxoptfeas = cum_samplesize(idxFeas);
                    xxsFeasible_cell{i,k,l}(idx_feas_stt:idx_feas_end) = xxoptfeas;

                    %% 3. feasiblity gap
                    feasgap = feasibilityGap(constraints.constraints{1}, meanFunc.fnEval, xMuMinHistory(n_init_sample+1:nEvaluations,:));
                    feasgap_cell{i,k,l}(idx_pred_stt:idx_pred_end)= feasgap;

                    if(isfield(result,'funcFeas') && exist('funcFeas','var'))
                        %% 4. probability of the constraints being satisfied
                        [~, pCons] = funcFeas(xMuMinHistory(n_init_sample+1:nEvaluations,:), 0);
                        probCons_cell{i,k,l}(idx_pred_stt:idx_pred_end)= pCons;
                    end
                    idx_feas_stt = idx_feas_end + 1;

          %% ----- LIBERALS --------------------------------------
                    if(hasLiberal )
                        %% 1. feasible
                        [feasible, strCon] =  isFeasible(constraints.constraints{1}, meanFunc.fnEval , xMuMinLiberalHistory(n_init_sample+1:nEvaluations,:));
                        feasible_liberal_cell{i,k,l}(idx_pred_stt:idx_pred_end)  = feasible;

                        %% 2. optimal among feasible (val(idxFeas))
                        idxFeas = find(feasible);
                        idx_feas_liberal_end = idx_feas_liberal_end + numel(idxFeas);
                        if(show_exact_fmin)
                            eval_f =   eval_f_at_xnew_liberal_history(n_init_sample+1:nEvaluations,1) ;
                            val = eval_f;
                        else
                            diff_f_min = (  eval_f_at_xnew_liberal_history(n_init_sample+1:nEvaluations,1)  - log(objFunc.getOptVal() ));
                            val = diff_f_min;
                        end
                        optfeas = val(idxFeas);
                        optFeasible_liberal_cell{i,k,l}(idx_feas_liberal_stt:idx_feas_liberal_end)= optfeas;
                        xxoptfeas = cum_samplesize(idxFeas);
                        xxsFeasible_liberal_cell{i,k,l}(idx_feas_liberal_stt:idx_feas_liberal_end) = xxoptfeas;

                        %% 3. feasiblity gap
                        feasgap = feasibilityGap(constraints.constraints{1}, meanFunc.fnEval, xMuMinLiberalHistory(n_init_sample+1:nEvaluations,:));
                        feasgap_liberal_cell{i,k,l}(idx_pred_stt:idx_pred_end)= feasgap;

                        if(isfield(result,'funcFeas') && exist('funcFeas','var'))
                            %% 4. probability of the constraints being satisfied
                            [~, pCons] = funcFeas(xMuMinLiberalHistory(n_init_sample+1:nEvaluations,:), 0);
                            probCons_liberal_cell{i,k,l}(idx_pred_stt:idx_pred_end)= pCons;
                        end
                        idx_feas_liberal_stt = idx_feas_liberal_end + 1;
                    end

                end
                % ----------------------------------------------------
                idx_pred_stt = idx_pred_end + 1;

                xxs_cell{i,k,l} = xxs_cell{i,k,l}(1:idx_pred_stt-1);
                yys_cell{i,k,l} = yys_cell{i,k,l}(1:idx_pred_stt-1);
                
                if(param_BO.hasConstraints)
                    feasible_cell{i,k,l} = feasible_cell{i,k,l}(1:idx_pred_stt-1);

                    xxsFeasible_cell{i,k,l} = xxsFeasible_cell{i,k,l}(1:idx_feas_stt-1);
                    optFeasible_cell{i,k,l} = optFeasible_cell{i,k,l}(1:idx_feas_stt-1);
                    feasgap_cell{i,k,l} = feasgap_cell{i,k,l}(1:idx_pred_stt-1);
                    if(isfield(result,'funcFeas') && exist('funcFeas','var'))
                        probCons_cell{i,k,l} = probCons_cell{i,k,l}(1:idx_pred_stt-1);
                    end
                    
                    if(hasLiberal )
                        xxsFeasible_liberal_cell{i,k,l} = xxsFeasible_liberal_cell{i,k,l}(1:idx_feas_liberal_stt-1);
                        optFeasible_liberal_cell{i,k,l} = optFeasible_liberal_cell{i,k,l}(1:idx_feas_liberal_stt-1);
                        feasgap_liberal_cell{i,k,l} = feasgap_liberal_cell{i,k,l}(1:idx_pred_stt-1);
                        if(isfield(result,'funcFeas') && exist('funcFeas','var'))
                            probCons_liberal_cell{i,k,l} = probCons_liberal_cell{i,k,l}(1:idx_pred_stt-1);
                        end
                    end
                end

                ypred = smooth(xxs_cell{i,k,l},yys_cell{i,k,l},.1,'loess');
                ypred_cell{i,k,l} = ypred;
               

                if(param_BO.hasConstraints)
                    feasible_pred_cell{i,k,l} = smooth(xxs_cell{i,k,l},feasible_cell{i,k,l},.1,'loess');
                    optFeasible_pred_cell{i,k,l} = smooth(xxsFeasible_cell{i,k,l},optFeasible_cell{i,k,l},.1,'loess');
                    feasgap_pred_cell{i,k,l} = smooth(xxs_cell{i,k,l},feasgap_cell{i,k,l},.1,'loess');
                    if(isfield(result,'funcFeas') && exist('funcFeas','var')) 
                        probCons_pred_cell{i,k,l} = smooth(xxs_cell{i,k,l},probCons_cell{i,k,l},.1,'loess');
                    end
                    
                    if(hasLiberal )
                        feasible_liberal_pred_cell{i,k,l} = smooth(xxs_cell{i,k,l},feasible_liberal_cell{i,k,l},.1,'loess');
                        optFeasible_liberal_pred_cell{i,k,l} = smooth(xxsFeasible_liberal_cell{i,k,l},optFeasible_liberal_cell{i,k,l},.1,'loess');
                        feasgap_liberal_pred_cell{i,k,l} = smooth(xxs_cell{i,k,l},feasgap_liberal_cell{i,k,l},.1,'loess');
                        if(isfield(result,'funcFeas') && exist('funcFeas','var')) 
                            probCons_liberal_pred_cell{i,k,l} = smooth(xxs_cell{i,k,l},probCons_liberal_cell{i,k,l},.1,'loess');
                        end
                    end
                end
                
                
%                 showIndividual
                if(showIndividual)
                    if(param_BO.hasConstraints)

                        figIndCon1 = figure(idxIndCon1);
                        sp = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));

                        if(show_exact_fmin)
                            if(exist('objFunc','var')&& numel(objFunc)>0)
                                hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
                            elseif(exist('opt_val_fn','var')&& numel(opt_val_fn)>0)
                                hideLegend(plot([ sttX,budget]', log([opt_val_fn, opt_val_fn]'),'k--' ));
                            end
                        else
                            hideLegend(plot([ sttX,budget]', [0, 0]','k--' ));
                        end
                        hold on;
                        pp = plot(xxoptfeas, optfeas , marks{i});
                        
                        pp.Color(4)=.2;
                        hideLegend( pp ,arr_samplesize{i});
                        mm = plot(xxoptfeas(end), optfeas(end),'bo');
                        mm.Color(4)=.2;
                        hideLegend( mm ,arr_samplesize{i});
                        showTitles(objFunc, str_setting, k, l);

                        if(minObj_over_funcs(k) > min(optfeas))
                            minObj_over_funcs(k) = min(optfeas);
                        end
                        if(maxObj_over_funcs(k) < max(optfeas))
                            maxObj_over_funcs(k) = max(optfeas);
                        end
                        sp.YLim = [minObj_over_funcs(k), maxObj_over_funcs(k)];
                    else
                        figure(idxIndUncon1);
                        sp = mysubplot(mysubplot_frame(1),mysubplot_frame(2),fnGetIdxOfSubplotWithRowCol(mysubplot_frame, i, fnGetIdxOfSubplotWithRowCol([dim_func,num_settings],k,l)));
                        hold on;
                        cum_samplesize =  cumsum(samplesize_history(n_init_sample+1:nEvaluations,1));

                        obj_f_min = (  eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1)  );
                        showOneLegend(    plot(  cum_samplesize, obj_f_min , marks{i}, [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), diff_f_min(end), 'ob') ,'Diff(y)');
                        showTitles(objFunc, str_setting, k, l);
                        
                        if(minObj_over_funcs(k) > min(obj_f_min))
                            minObj_over_funcs(k) = min(obj_f_min);
                        end
                        if(maxObj_over_funcs(k) < max(obj_f_min))
                            maxObj_over_funcs(k) = max(obj_f_min);
                        end
                        
                        sp.YLim = [minObj_over_funcs(k), maxObj_over_funcs(k)];

                    end

                    
                end

            end


            arr_samplesize = {'adj','3','10','50'};
            % names = {'adj',}
%% ----------------------------------------------------------------------
            
            
            
            if(~hasValue);continue;end;
            upperBound = cum_samplesize;
            [domainX, idxDomain] = sort(xxs_cell{i,k,l});
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            if(param_BO.hasConstraints)
                    % 2. optimal among feasible (val(idxFeas))
                if(~showOnlySummary)
                    figCon1 = figure(idxSumCon1);
                    sp = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));

                    if(show_exact_fmin)
                        if(exist('objFunc','var')&& numel(objFunc)>0)
                            hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
                        elseif(exist('opt_val_fn','var')&& numel(opt_val_fn)>0)
                            hideLegend(plot([ sttX,budget]', log([opt_val_fn, opt_val_fn]'),'k--' ));
                        end
                    else
                        hideLegend(plot([ sttX,budget]', [0, 0]','k--' ));
                    end
                    hold on;

                    [domainFeasX, idxFeasDomain] = sort(xxsFeasible_cell{i,k,l});
                    yPred = optFeasible_pred_cell{i,k,l};
                    showOneLegend( plot(domainFeasX, yPred(idxFeasDomain),marks{i}),arr_samplesize{i} )

                    if(hasLiberal )
                        [domainFeasX, idxFeasDomain] = sort(xxsFeasible_liberal_cell{i,k,l});
                        yPred2 = optFeasible_liberal_pred_cell{i,k,l};
                        showOneLegend( plot(domainFeasX, yPred2(idxFeasDomain),marks_liberal{i}),arr_samplesize{i} )
                    end

                    showTitles(objFunc, str_setting, k, l);

                    if(minMeanObj_over_funcs(k) > min(yPred))
                        minMeanObj_over_funcs(k) = min(yPred);
                    end
                    if(maxMeanObj_over_funcs(k) < max(yPred))
                        maxMeanObj_over_funcs(k) = max(yPred);
                    end
                    if(minMeanObj_over_funcs(k) > min(yPred2))
                        minMeanObj_over_funcs(k) = min(yPred2);
                    end
                    if(maxMeanObj_over_funcs(k) < max(yPred2))
                        maxMeanObj_over_funcs(k) = max(yPred2);
                    end
                    if( minMeanObj_over_funcs(k) == maxMeanObj_over_funcs(k))
                        minMeanObj_over_funcs(k) = minMeanObj_over_funcs(k) *.95;
                        maxMeanObj_over_funcs(k) = maxMeanObj_over_funcs(k) .*1.05;
                    end
                    sp.YLim = [minMeanObj_over_funcs(k), maxMeanObj_over_funcs(k)];



                    legend();
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        % 1. feasible
                    figCon2 = figure(idxSumCon2);
                    sp12 = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));

                        hideLegend(plot([ sttX,budget]', log([0, 0]'),'k--' ));
                        hideLegend(plot([ sttX,budget]', log([1, 1]'),'k--' ));
                    hold on;
                    yPred = feasible_pred_cell{i,k,l};
                    showOneLegend( plot(domainX, yPred(idxDomain),marks{i}),arr_samplesize{i} )
                    if(hasLiberal )
                        yPred = feasible_liberal_pred_cell{i,k,l};
                        showOneLegend( plot(domainX, yPred(idxDomain),marks_liberal{i}),arr_samplesize{i} )
                        showTitles(objFunc, str_setting, k, l);
                    end

                    legend();
                    sp12.YLim = [-0.05 1.05];
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        % 3. feasiblity gap
                    figCon3 = figure(idxSumCon3);
                    sp13 = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));

                    hold on;
                    yPred = feasgap_pred_cell{i,k,l};
                    showOneLegend( plot(domainX, yPred(idxDomain),marks{i}),arr_samplesize{i} )

                    if(hasLiberal )
                        yPred = feasgap_liberal_pred_cell{i,k,l};
                        showOneLegend( plot(domainX, yPred(idxDomain),marks_liberal{i}),arr_samplesize{i} )
                        showTitles(objFunc, str_setting, k, l);
                    end

                    legend();
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                    if(isfield(result,'funcFeas') && exist('funcFeas','var')) 

                            % 4. probability of the constraints being satisfied
                        figCon4 = figure(idxSumCon4);
                        sp14 = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));

                            hideLegend(plot([ sttX,budget]', log([0, 0]'),'k--' ));
                            hideLegend(plot([ sttX,budget]', log([1, 1]'),'k--' ));
                        hold on;
                        yPred = probCons_pred_cell{i,k,l};
                        showOneLegend( plot(domainX, yPred(idxDomain),marks{i}),arr_samplesize{i} )

                        if(hasLiberal )
                            yPred = probCons_liberal_pred_cell{i,k,l};
                            showOneLegend( plot(domainX, yPred(idxDomain),marks_liberal{i}),arr_samplesize{i} )
                            showTitles(objFunc, str_setting, k, l);
                        end

                        legend();
                        sp14.YLim = [-0.05, 1.05];
                    end
                else
                    figUncon1 =figure(idxSumUncon1);
                    sp2 = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));

                    if(show_exact_fmin)
                        if(exist('objFunc','var')&& numel(objFunc)>0)
                            hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
                        elseif(exist('opt_val_fn','var')&& numel(opt_val_fn)>0)
                            hideLegend(plot([ sttX,budget]', log([opt_val_fn, opt_val_fn]'),'k--' ));
                        end
                    else
                        hideLegend(plot([ sttX,budget]', [0, 0]','k--' ));
                        sp2.YLim(1) = -.05;
                    end
                    hold on;
                    [domainX, idxDomain] = sort(xxs_cell{i,k,l});
                    yPred = ypred_cell{i,k,l};
                    showOneLegend( plot(domainX, yPred(idxDomain),marks{i}),arr_samplesize{i} )
                    showTitles(objFunc, str_setting, k, l);

                    if(minMeanObj_over_funcs(k) > min(yPred))
                        minMeanObj_over_funcs(k) = min(yPred);
                    end
                    if(maxMeanObj_over_funcs(k) < max(yPred))
                        maxMeanObj_over_funcs(k) = max(yPred);
                    end
                    sp2.YLim = [minMeanObj_over_funcs(k), maxMeanObj_over_funcs(k)];



                    legend();
                    end
                end
        %% ----------------------------------------------------------------------
            
            
%             annotation(fig2,'textbox',...
%                     [0.0125 0.934714286920571 0.808928568714431 0.0642857130794299],...
%                     'LineStyle','none',...
%                     'FontWeight','bold',...
%                     'String',{sprintf('%s (%s) power = %.2g',objFunc.name, str_setting , power )});

            ;
            % names = {'adj',}
    %% ----------------------------------------------------------------------
    %         fig3 =figure(50*(l-1)+20+k);clf;
%             fig3 =figure(3);
%             mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));
% %             hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
%             if(exist('objFunc','var')&& numel(objFunc)>0)
%                 hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
%             elseif(exist('opt_val_fn','var')&& numel(opt_val_fn)>0)
%                 hideLegend(plot([ sttX,budget]', log([opt_val_fn, opt_val_fn]'),'k--' ));
%             end
% %                 for i = 2:dim_cases_samplesize
%             if(i>1)
%                 hold on;
%                 res = mean(eval_f_at_xnew_history,2);
%                 yy = res(n_init_sample+1:nEvaluations);
%                 xx = cumsum(samplesize_history(n_init_sample+1:nEvaluations,1)) +  sum(samplesize_history(1:n_init_sample));
%                 showOneLegend( plot( xx,yy,marks{i}), arr_samplesize{i});
%             end
% %                 end
%             showTitles(objFunc, str_setting, k, l);
% %             if(k==1&&l==1)
% %                 title(sprintf('%s(%s)',objFunc.name, str_setting ))
% %             elseif(l==1)
% %                 title(sprintf('%s',objFunc.name ))
% %             elseif(k==1)
% %                 title(sprintf('%s', str_setting ))
% %             else
% %                 title(sprintf('%s(%s)',objFunc.name, str_setting ))
% %             end
%             legend();
%             
            figSummary =figure(idxSummary);
            mysubplot(dim_func,1,k);
            if(show_exact_fmin)
                if(exist('objFunc','var')&& numel(objFunc)>0)
                    hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
                elseif(exist('opt_val_fn','var')&& numel(opt_val_fn)>0)
                    hideLegend(plot([ sttX,budget]', log([opt_val_fn, opt_val_fn]'),'k--' ));
                end
            else
                hideLegend(plot([ sttX,budget]', [0, 0]','k--' ));
            end
            if(param_BO.hasConstraints)
                [domainFeasX, idxFeasDomain] = sort(xxsFeasible_cell{i,k,l});
                yPred = optFeasible_pred_cell{i,k,l};
                showOneLegend( plot(domainFeasX, yPred(idxFeasDomain),marks{dim_cases_samplesize*(l-1)+i}),sprintf('%s:%s', arr_samplesize{i}, str_setting) )
                
                if(hasLiberal)
                    [domainFeasX, idxFeasDomain] = sort(xxsFeasible_liberal_cell{i,k,l});
                    yPred = optFeasible_liberal_pred_cell{i,k,l};
                    showOneLegend( plot(domainFeasX, yPred(idxFeasDomain),marks_liberal{dim_cases_samplesize*(l-1)+i}),sprintf('%s:%s (5%%)', arr_samplesize{i}, str_setting) )
                end
            else
                [domainX, idxDomain] = sort(xxs_cell{i,k,l});
                yPred = ypred_cell{i,k,l};
                showOneLegend( plot(domainX, yPred(idxDomain),marks{dim_cases_samplesize*(l-1)+i}),sprintf('%s:%s', arr_samplesize{i}, str_setting) )
            end
%             hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
%                 for i = 1:dim_cases_samplesize
            hold on;
            showTitles(objFunc, str_setting, k, l);

            legend()




% 
%             if(~showDetail)
% 
%             else
% 
%                 idx_subfigure_row = i;
%                 mysubplot(mysubplot_frame(1),mysubplot_frame(2),fnGetIdxOfSubplotWithRowCol(mysubplot_frame, idx_subfigure_row, 1));
%                 hold on;
%                 cum_samplesize =  cumsum(samplesize_history(n_init_sample+1:nEvaluations,1));
% 
%                 diff_f_min = (  eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1)  - log(opt_val_fn ));
% %                     idxNInf = find(eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1) == -inf);
% %                     minDiff = min(diff_f_min(find(eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1) > -inf)));
%                     showOneLegend(    plot(  cum_samplesize, diff_f_min , 'r-', [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), diff_f_min(end), 'ob') ,'Diff(y)');
%                     title(sprintf('|f_{min} - min(mu)|: %.3g',diff_f_min(end)))
% 
%                     mysubplot(mysubplot_frame(1),mysubplot_frame(2),fnGetIdxOfSubplotWithRowCol(mysubplot_frame, idx_subfigure_row, 2));
%                 hold on;
%                 %     diff_x_f_min = vecnorm(xMuMinHistory(n_init_sample+1:nEvaluations,:) - opt_sol_fn, 2, 2);
%                 %     diff_x_f_min = min(vecnorm(xMuMinHistory(n_init_sample+1:nEvaluations,:) - opt_sol_fn, 2, 2));
%                     diff_x_f_min = getMinimumDistance(xMuMinHistory(n_init_sample+1:nEvaluations,:), opt_sol_fn);
%                     plot( cum_samplesize, diff_x_f_min, 'r-' , [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), diff_x_f_min(end), 'ob');
%                     title(sprintf('|x_{min} - argmin(mu)|: %.3g',diff_x_f_min(end)))
%                 % mysubplot(mysubplot_frame(1),mysubplot_frame(2),9);
%                 % idxFigStatus = 3;
%                 mysubplot(mysubplot_frame(1),mysubplot_frame(2),fnGetIdxOfSubplotWithRowCol(mysubplot_frame, idx_subfigure_row, 3));
%                 hold on;
%                 plot( cum_samplesize, samplesize_history(n_init_sample+1:nEvaluations,1), 'r-' , [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), samplesize_history(nEvaluations,1), 'ob');
%                 title(sprintf('Sample Size: %d',samplesize_history(nEvaluations,1) ))
% 
%                 mysubplot(mysubplot_frame(1),mysubplot_frame(2),fnGetIdxOfSubplotWithRowCol(mysubplot_frame, idx_subfigure_row, 4));
% 
%                 hold on;
%                 plot( cum_samplesize, iterElapsedTime(n_init_sample+1:nEvaluations,1), 'r-' , [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), iterElapsedTime(nEvaluations,1), 'ob');
%                 title(sprintf('Elapsed Time for iteration: %.1f sec',iterElapsedTime(nEvaluations,1)))
%             end


        end
        if(showDetail)
            annotation(fig1,'textbox',...
                    [0.0125 0.934714286920571 0.808928568714431 0.0642857130794299],...
                    'LineStyle','none',...
                    'FontWeight','bold',...
                    'String',{sprintf('%s (%s)',objFunc.name, str_setting )});
        end
    % end
    % 


        if(showDetail)
            figUncon1 =figure(50*(l-1)+10+k);clf;
            mysubplot(1,2,1);
            hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
            for i = 1:dim_cases_samplesize
                hold on;
                showOneLegend( plot(domainX, ypred_cell{i},marks{i}),arr_samplesize{i} )
            end
    %         title(sprintf('(%s,%s)',objFunc.name, str_setting ))
            legend()
            annotation(figUncon1,'textbox',...
                    [0.0125 0.934714286920571 0.808928568714431 0.0642857130794299],...
                    'LineStyle','none',...
                    'FontWeight','bold',...
                    'String',{sprintf('%s (%s) power = %.2g',objFunc.name, str_setting , power )});

            ;
            % names = {'adj',}
    %% ----------------------------------------------------------------------
    %         fig3 =figure(50*(l-1)+20+k);clf;
            mysubplot(1,2,2);
            hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
            for i = 2:dim_cases_samplesize
                hold on;
                res = mean(eval_f_at_xnew_history,2);
                yy = res(n_init_sample+1:nEvaluations);
                xx = cumsum(samplesize_history(n_init_sample+1:nEvaluations,1)) +  sum(samplesize_history(1:n_init_sample));
                showOneLegend( plot( xx,yy,marks{i}), arr_samplesize{i});
            end
            legend();
    %         annotation(fig3,'textbox',...
    %                 [0.0125 0.934714286920571 0.808928568714431 0.0642857130794299],...
    %                 'LineStyle','none',...
    %                 'FontWeight','bold',...
    %                 'String',{sprintf('%s (%s)',objFunc.name, str_setting )});
        end
        
        
        

    end
end


for l=1:num_settings
	for k=1:dim_func
		if(param_BO.hasConstraints)
            if(showIndividual)
                figCon1 = figure(idxIndCon1);
                sp = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));
                sp.YLim = [minObj_over_funcs(k), maxObj_over_funcs(k)];
            end
			
			figCon1 = figure(idxSumCon1);
			sp = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));
			sp.YLim = [minMeanObj_over_funcs(k), maxMeanObj_over_funcs(k)];

		else
            if(showIndividual)
                figure(idxIndUncon1);
                sp = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));
                sp.YLim = [minObj_over_funcs(k), maxObj_over_funcs(k)];
            end
			
			figure(idxSumUncon1);
			sp = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));
			sp.YLim = [minMeanObj_over_funcs(k), maxMeanObj_over_funcs(k)];
		end

	end
end
% if(param_BO.hasConstraints)

% - - - ANNOTATION - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    nRow = dim_func;
    nCol = num_settings;
    strRows = cell(nRow,1);
    strCols = cell(nCol,1);
    for idx_row=1:nRow;
        strRows{idx_row} = cell_obj_func{idx_row}.name;
    end
%     for idx_col=1:nCol;
%         setting = type_cell_settings{l};
%         if(iscell(setting))
%             str_setting = sprintf('%s,',setting{:} );
%         else
%             str_setting = sprintf('%s,',setting );
%         end
%         strCols{idx_col} =  sprintf('%s\n',str_setting );
%     end
    if(exist('figCon1','var'));
        if(exist('show_exact_fmin','var'))
            if(show_exact_fmin)
                strTitleObjFunc = sprintf('f_{min}' );
            else
                strTitleObjFunc = sprintf('|f_{min} - min(mu)|' );
            end
            strTitle = sprintf('Objective Optimal Value (%s) (Only Feasible)', strTitleObjFunc);
        else
            strTitle = sprintf('Objective Optimal Value (Only Feasible)');
        end
        annotateTextOnAxes(figCon1, nRow, nCol, strRows, strCols, strTitle)        
        figCon1.Name = 'Constrainted Optimization: Objective Value(Feasible)';
    end;
    if(exist('figCon2','var'));
        strTitle = sprintf('Fraction of Feasible evaluations');
        annotateTextOnAxes(figCon2, nRow, nCol, strRows, strCols, strTitle)        
        figCon2.Name = 'Constrainted Optimization: Fraction of Feasible evaluations';
    end;
    if(exist('figCon3','var'));
        strTitle = sprintf('Feasibility Gap');
        annotateTextOnAxes(figCon3, nRow, nCol, strRows, strCols, strTitle)        
        figCon3.Name = 'Constrainted Optimization: Feasibility Gap';
    end;
    if(exist('figCon4','var'));
        strTitle = sprintf('P(C(x))' );
        annotateTextOnAxes(figCon4, nRow, nCol, strRows, strCols, strTitle)        
        figCon4.Name = 'Constrainted Optimization: P(C(x))';
    end;
    if(exist('figSummary','var'));
        figSummary.Name = 'Summary: Optimization: Objective Value';
    end;

    if(exist('figUncon1','var'));
        strTitle = sprintf('GP: |f_{min} - min(mu)| (power = %.2g)', power );
        annotateTextOnAxes(figUncon1, nRow, nCol, strRows, strCols, strTitle)        
        figUncon1.Name = 'Unconstrainted Optimization: Objective Value';
    end;
% end

if(~showDetail)
%     annotation(fig1,'textbox',...
%             [0.0125 0.934714286920571 0.808928568714431 0.0642857130794299],...
%             'LineStyle','none',...
%             'FontWeight','bold',...
%             'String',{sprintf('%s','|f_{min} - min(mu)| (power = %.2g)', power )});
    nRow = dim_cases_samplesize;
    nCol = dim_func;
    strRows = arr_samplesize;
    strCols = cell(nCol,1);
    for i=1:nCol;
        strCols{i} = cell_obj_func{i}.name;
    end
%     for i=1:dim_cases_samplesize;
%         strCols(i) = cell_obj_func{i}.name;
%     end
    
%     strRows = {'r1','r2','r3','r4'};
%     strCols = {'c1','c2','c3'};
    if(exist('fig1','var'));
        strTitle = sprintf('%s (power = %.2g)','|f_{min} - min(mu)|', power );
        annotateTextOnAxes(fig1, nRow, nCol, strRows, strCols, strTitle)        
    end
    if(exist('fig3','var'));
        strTitle = sprintf('|f_{min} - min(mu)| (power = %.2g)', power );
        annotateTextOnAxes(fig3, nRow, nCol, strRows, strCols, strTitle)        
    end
end
        
%     annotation(fig2,'textbox',...
%         [0.0125 0.934714286920571 0.808928568714431 0.0642857130794299],...
%         'LineStyle','none',...
%         'FontWeight','bold',...
%         'String',{sprintf('GP: |f_{min} - min(mu)| (power = %.2g)', power )});
%     nRow = dim_func;
%     nCol = num_settings;
%     strRows = cell(nRow,1);
%     strCols = cell(nCol,1);
%     for i=1:nRow;
%         strRows{i} = cell_obj_func{i}.name;
%     end
%     for i=1:nCol;
%         setting = type_cell_settings{l};
%         if(iscell(setting))
%             str_setting = sprintf('%s,',setting{:} );
%         else
%             str_setting = sprintf('%s,',setting );
%         end
%         strCols{i} =  sprintf('%s\n',str_setting );
%     end
%     strRows = {'r1','r2','r3','r4'};
%     strCols = {'c1','c2','c3'};

%     annotation(fig3,'textbox',...
%         [0.0125 0.934714286920571 0.808928568714431 0.0642857130794299],...
%         'LineStyle','none',...
%         'FontWeight','bold',...
%         'String',{sprintf('|f_{min} - min(mu)| (power = %.2g)', power )});
%     nRow = 4;
%     nCol = 3;
%     strRows = {'r1','r2','r3','r4'};
%     strCols = {'c1','c2','c3'};
% end
%     mean(merged_results{1}.final_opt_f_true)
%     mean(merged_results{2}.final_opt_f_true)
%     mean(merged_results{3}.final_opt_f_true)
%     mean(merged_results{4}.final_opt_f_true)
% 
%     median(merged_results{1}.final_opt_f_true)
%     median(merged_results{2}.final_opt_f_true)
%     median(merged_results{3}.final_opt_f_true)
%     median(merged_results{4}.final_opt_f_true)
% 
%     var(merged_results{1}.final_opt_f_true)
%     var(merged_results{2}.final_opt_f_true)
%     var(merged_results{3}.final_opt_f_true)
%     var(merged_results{4}.final_opt_f_true)


function showTitles(objFunc, str_setting, k, l)
    if(exist('objFunc','var')&& numel(objFunc)>0)
        if(k==1&&l==1)
            title(sprintf('%s(%s)',objFunc.name, str_setting ))
        elseif(l==1)
            title(sprintf('%s',objFunc.name ))
        elseif(k==1)
            title(sprintf('%s', str_setting ))
        else
            title(sprintf('%s(%s)',objFunc.name, str_setting ))
        end
    end
end