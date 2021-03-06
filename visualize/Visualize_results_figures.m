MEDIAN = 2;

pred_method = MEDIAN;
ignoreLastEval = true;
isFirstRun = true;
% isFirstRun = false;
% summarize = false;
% showIndividual=true;
showIndividual=false;
% showDetail = true;
showDetail = false;

show_only_ours = true;

show_pcon = false;

samplesize_eval = 2;

% showOnlyFeasible = false;
showOnlyFeasible = true;

show_fmin_diff =false;
clc;
% clearvars -except results_all;
    close all;

%     show_median = true;
    show_median = false;
show_errorbar = false;
% show_errorbar = true;
    show_mean = true;
    
mysubplot = @(x1,x2,x3) subplot(x1,x2,x3);
marks = {'r-','g--','m-.','k:','b:x','ro','go', 'r--','b--','c--','g--','k--','m--','y--', 'r:','b:','c:','g:','k:','m:','y:'};
% marks = {'c-','g--','k--','b:','r:','m-','y-', 'r--','b--','c--','g--','k--','m--','y--', 'r:','b:','c:','g:','k:','m:','y:'};
marks_liberal = {'r:','b:','c:','g:','k:','m:','y:', 'r-.','b-.','c-.','g-.','k-.','m-.','y-.'};
if(isFirstRun)
    [dim_cases_samplesize, dim_repeats, dim_func, num_settings] = size(results_all);
end

if(showDetail)
    mysubplot_frame = [dim_cases_samplesize,4];
else
    mysubplot_frame = [dim_cases_samplesize,num_settings*dim_func];
end
idx_subfigure_row = 1;


%%
%%
    idxSumCon3 = 23;
%     idxSumCon4 = 24;
%     idxSumUncon1 = 3;
    idxSummary = 4;
    
predicted_cell = cell(dim_cases_samplesize, dim_func, num_settings);
max_rep_done = zeros(dim_cases_samplesize, dim_func, num_settings);

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
    optAll_cell = cell(dim_cases_samplesize, dim_func, num_settings);

    % pred
    feasible_pred_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
    xxsFeasible_pred_cell = cell(dim_cases_samplesize, dim_func, num_settings);
    optFeasible_pred_cell = cell(dim_cases_samplesize, dim_func, num_settings);
    feasgap_pred_cell     = cell(dim_cases_samplesize, dim_func, num_settings);
    probCons_pred_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
    optAll_pred_cell = cell(dim_cases_samplesize, dim_func, num_settings);

%     if(hasLiberal )
            % 2) liberal
        feasible_liberal_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
        xxsFeasible_liberal_cell  = cell(dim_cases_samplesize, dim_func, num_settings);
        optFeasible_liberal_cell = cell(dim_cases_samplesize, dim_func, num_settings);
        feasgap_liberal_cell     = cell(dim_cases_samplesize, dim_func, num_settings);
        probCons_liberal_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
        optAll_liberal_cell = cell(dim_cases_samplesize, dim_func, num_settings);

        feasible_liberal_pred_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
        xxsFeasible_liberal_pred_cell = cell(dim_cases_samplesize, dim_func, num_settings);
        optFeasible_liberal_pred_cell = cell(dim_cases_samplesize, dim_func, num_settings);
        feasgap_liberal_pred_cell     = cell(dim_cases_samplesize, dim_func, num_settings);
        probCons_liberal_pred_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
        optAll_liberal_pred_cell = cell(dim_cases_samplesize, dim_func, num_settings);
%     end
    
end

minObj_over_funcs = Inf * ones(dim_func,1);
maxObj_over_funcs = -Inf * ones(dim_func,1);
minMeanObj_over_funcs = Inf * ones(dim_func,1);
maxMeanObj_over_funcs = -Inf * ones(dim_func,1);
                
tot_exp = dim_func* dim_repeats* num_settings* dim_cases_samplesize;
cnt_exp = 1;
sttVis = tic();
sttLastIter = tic();
% 
%% Loop: Settings
%% Loop: Functions 
for k=1:dim_func;
    clear yall
    for l=1:num_settings
        setting = type_cell_settings{l};
        if(iscell(setting))
            name = 0;
            for i=1:numel(setting);
                if(isa(setting{i},'TypeInfillOptimizer'))
                    name = setting{i}.char;
                    break;
                end
            end
            if(name==0)
                name = 'Unknown';
            else
                if(numel(name)>12)
                    name = name(12:end);
                end
            end
            str_setting = sprintf('%s',name);
        else
            str_setting = sprintf('%s,',setting );
        end
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
            isfirst  = true;
            for j=1:dim_repeats
                if(toc(sttLastIter)>60*1e3) % show every 60 sec.
                    fprintf('[Visualization: Results] func:%d/%d, rep:%d/%d, set:%d/%d, samp:%d/%d (%d/%d)\n', ...
                        k, dim_func, j, dim_repeats, l, num_settings,...
                        i, dim_cases_samplesize,...
                        cnt_exp, tot_exp);
                    showEstimatedEndRemainedTime(toc(sttVis), cnt_exp-1, tot_exp);
                    sttLastIter = tic();
                end
                clearvars -except i j k l idx_pred_stt idx_pred_end idx_feas_stt idx_feas_end idx_feas_liberal_stt idx_feas_liberal_end hasValue fig1 str_setting setting minObj_over_funcs  maxObj_over_funcs  minMeanObj_over_funcs  maxMeanObj_over_funcs  predicted_cell  dim_all   yys_cell  xxs_cell  ypred_cell  feasible_cell     xxsFeasible_cell   optFeasible_cell  feasgap_cell      probCons_cell     optAll_cell  feasible_pred_cell     xxsFeasible_pred_cell  optFeasible_pred_cell  feasgap_pred_cell      probCons_pred_cell     optAll_pred_cell  feasible_liberal_cell     xxsFeasible_liberal_cell   optFeasible_liberal_cell  feasgap_liberal_cell      probCons_liberal_cell     optAll_liberal_cell  feasible_liberal_pred_cell     xxsFeasible_liberal_pred_cell  optFeasible_liberal_pred_cell  feasgap_liberal_pred_cell      probCons_liberal_pred_cell     optAll_liberal_pred_cell  dim_cases_samplesize dim_repeats dim_func num_settings showDetail mysubplot_frame  idx_subfigure_row  idxSumCon1  idxSumCon2  idxSumCon3  idxSumCon4  idxSumUncon1  idxSummary  idxIndCon1  idxIndCon2  idxIndCon3  idxIndCon4  idxIndUncon1  isFirstRun  showIndividual showDetail  showOnlyFeasible  show_fmin_diff  mysubplot  marks  marks_liberal  results_all type_cell_settings optProb simul param_BO cell_obj_func A_RUN_PURPOSE show_median show_exact_fmin show_fmin_diff hasLiberal sttX budget ignoreLastEval samplesize_eval sttVis cnt_exp tot_exp sttLastIter cell_feas_mat cell_obj_mat cell_feasgap_mat isfirst n_init_sample show_errorbar max_rep_done show_only_ours prev_ylim pred_method SMOOTH SMOOTH_MEDIAN MEDIAN typeSamplesize givenSamplesize yall show_pcon;  
                result = results_all{i,j,k,l};

                if ([0,0]==size(result))
                    continue;
                end
                max_rep_done(i,k,l) = max_rep_done(i,k,l) + 1;

                hasValue = true;
                scr_setResultValuesToVariables
                
                budget = param_BO.budget;
                
                if(~show_pcon);clear funcFeas;end;
                if(isfirst)
                    cell_feas_mat{i,k,l} = zeros(budget ,dim_repeats);
                    cell_obj_mat{i,k,l} = zeros(budget ,dim_repeats);
                    cell_feasgap_mat{i,k,l} = zeros(budget ,dim_repeats);
                    isfirst = false;
                end
                
                if(ignoreLastEval); nEvaluations = nEvaluations - 1;end;
                [yout, extended] = extendifshort(eval_f_at_xnew_history, budget);
                if(extended); eval_f_at_xnew_history = yout;end;
                [xout, extended] = extendifshort(xMuMinHistory, budget);
                if(extended); xMuMinHistory = xout;end;
                [xlout, extended] = extendifshort(xMuMinLiberalHistory, budget);
                if(extended); xMuMinLiberalHistory = xlout;end;
                [ssize, extended] = extendifshort(samplesize_history, budget);
                if(extended); samplesize_history = ssize;end;
                [ylout, extended] = extendifshort(eval_f_at_xnew_liberal_history, budget);
                if(extended); eval_f_at_xnew_liberal_history = ylout;end;
                
                if(setting.typeProblem ~= TypeProblem.RobustDesignOptimization)
                    for idx_value = 2:budget;
                        if(yout(idx_value)==0)
                            yout(idx_value) = yout(idx_value-1);
                        end
                        if(sum(xout(idx_value,:),2)==0 && nEvaluations<=idx_value)
                            xout(idx_value,:) = xout(idx_value-1,:);
                        end
%                         if(samplesize_history(idx_value)==0)
%                             (idx_value) = 1;
%                         end
                    end
                    samplesize_history = ones(budget,1);
                    xMuMinHistory = xout ;
                    eval_f_at_xnew_history = yout;
                    results_all{i,j,k,l}.eval_f_at_xnew_history = yout;
                    nEvaluations = budget;
                else
                    if(nEvaluations<=budget)
                        samplesize_history(nEvaluations+1) = budget - sum(samplesize_history);
                        xMuMinHistory(nEvaluations+1) = xMuMinHistory(nEvaluations);
                        eval_f_at_xnew_history(nEvaluations+1) =   eval_f_at_xnew_history(nEvaluations);
                        results_all{i,j,k,l}.eval_f_at_xnew_history = yout;eval_f_at_xnew_history;
%                         nEvaluations = budget;
                        nEvaluations = nEvaluations +1;
                    end
                end



                sttX = sum(samplesize_history(1:n_init_sample));

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
                    cum_samplesize =  cumsum(samplesize_history(n_init_sample+1:nEvaluations,1)) + sum(samplesize_history(1:n_init_sample,1));
                    %% 1. feasible
                    
                    
                    % This may be deleted for the future version since it would be saved automatically.
                    eval_g_at_xnew_history  = storage.evalGXnewHistory;
                    feas_history  = storage.feasAtXnewHistory;
                    feasgap_history  = storage.feasgapAtXnewHistory;

                    
%                     [feasible, strCon] =  isFeasible(constraints.constraints{1}, conFunc, xMuMinHistory(n_init_sample+1:nEvaluations,:));
                    feasible = feas_history(n_init_sample+1:nEvaluations);
                    feasible_cell{i,k,l}(idx_pred_stt:idx_pred_end)  = feasible;
                    cell_feas_mat{i,k,l}(n_init_sample+1:nEvaluations,j) = feasible;
                    
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
                    cell_obj_mat{i,k,l}(n_init_sample+1:nEvaluations,j) = val;
                    
                    optAll_cell{i,k,l}(idx_pred_stt:idx_pred_end) = val;
                    optfeas = val(idxFeas);
                    optFeasible_cell{i,k,l}(idx_feas_stt:idx_feas_end)= optfeas;
                    xxoptfeas = cum_samplesize(idxFeas);
                    xxsFeasible_cell{i,k,l}(idx_feas_stt:idx_feas_end) = xxoptfeas;

                    %% 3. feasiblity gap
%                     feasgap = feasibilityGap(constraints.constraints{1}, conFunc, xMuMinHistory(n_init_sample+1:nEvaluations,:));
                    feasgap = max(feasgap_history(n_init_sample +1:nEvaluations,:),[],2);

                    feasgap_cell{i,k,l}(idx_pred_stt:idx_pred_end)= feasgap;
                    cell_feasgap_mat{i,k,l}(n_init_sample+1:nEvaluations,j) = feasgap;

                    idx_feas_stt = idx_feas_end + 1;


                end
                % ----------------------------------------------------
                idx_pred_stt = idx_pred_end + 1;

                xxs_cell{i,k,l} = xxs_cell{i,k,l}(1:idx_pred_stt-1);
                yys_cell{i,k,l} = yys_cell{i,k,l}(1:idx_pred_stt-1);
                
                if(param_BO.hasConstraints)
                    feasible_cell{i,k,l} = feasible_cell{i,k,l}(1:idx_pred_stt-1);

                    xxsFeasible_cell{i,k,l} = xxsFeasible_cell{i,k,l}(1:idx_feas_stt-1);
                    optFeasible_cell{i,k,l} = optFeasible_cell{i,k,l}(1:idx_feas_stt-1);
                    optAll_cell{i,k,l} = optAll_cell{i,k,l}(1:idx_pred_stt-1);
                    feasgap_cell{i,k,l} = feasgap_cell{i,k,l}(1:idx_pred_stt-1);
                    if(isfield(result,'funcFeas') && exist('funcFeas','var'))
                        probCons_cell{i,k,l} = probCons_cell{i,k,l}(1:idx_pred_stt-1);
                    end
                    
                    if(hasLiberal )
                        xxsFeasible_liberal_cell{i,k,l} = xxsFeasible_liberal_cell{i,k,l}(1:idx_feas_liberal_stt-1);
                        optFeasible_liberal_cell{i,k,l} = optFeasible_liberal_cell{i,k,l}(1:idx_feas_liberal_stt-1);
                        feasgap_liberal_cell{i,k,l} = feasgap_liberal_cell{i,k,l}(1:idx_pred_stt-1);
                        optAll_liberal_cell{i,k,l} = optAll_liberal_cell{i,k,l}(1:idx_pred_stt-1);
                        if(isfield(result,'funcFeas') && exist('funcFeas','var'))
                            probCons_liberal_cell{i,k,l} = probCons_liberal_cell{i,k,l}(1:idx_pred_stt-1);
                        end
                    end
                end
                cnt_exp = cnt_exp +1;

            end

                ypred = smooth(xxs_cell{i,k,l},yys_cell{i,k,l},.1,'loess');
                ypred_cell{i,k,l} = ypred;
               


            arr_samplesize = {'adj','3','10','50'};
            % names = {'adj',}
%% ----------------------------------------------------------------------
            
            
            
            if(~hasValue);continue;end;
            [domainX, idxDomain] = sort(xxs_cell{i,k,l});
%             if(param_BO.hasConstraints)
%                     % 3. feasiblity gap
%                 
%                 if(isfield(result,'funcFeas') && exist('funcFeas','var')) 
% 
%                         % 4. probability of the constraints being satisfied
%                     figCon4 = figure(idxSumCon4);
%                     sp14 = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));
% 
%                         hideLegend(plot([ sttX,budget]', log([0, 0]'),'k--' ));
%                         hideLegend(plot([ sttX,budget]', log([1, 1]'),'k--' ));
%                     hold on;
%                     yPred = probCons_pred_cell{i,k,l};
%                     if(numel(yPred)>0)
%                         showOneLegend( plot(domainX, yPred(idxDomain),marks{i}),arr_samplesize{i} )
%                     end
% 
%                     if(hasLiberal )
%                         yPred = probCons_liberal_pred_cell{i,k,l};
%                         showOneLegend( plot(domainX, yPred(idxDomain),marks_liberal{i}),arr_samplesize{i} )
%                         showTitles(objFunc, str_setting, k, l);
%                     end
% 
%                     legend();
%                     sp14.YLim = [-0.05, 1.05];
%                 end
%             else
%                 figUncon1 =figure(idxSumUncon1);
%                 sp2 = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));
% 
%                 if(show_exact_fmin)
%                     if(exist('objFunc','var')&& numel(objFunc)>0)
%                         hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
%                     elseif(exist('opt_val_fn','var')&& numel(opt_val_fn)>0)
%                         hideLegend(plot([ sttX,budget]', log([opt_val_fn, opt_val_fn]'),'k--' ));
%                     end
%                 else
%                     hideLegend(plot([ sttX,budget]', [0, 0]','k--' ));
%                     sp2.YLim(1) = -.05;
%                 end
%                 hold on;
%                 [domainX, idxDomain] = sort(xxs_cell{i,k,l});
%                 yPred = ypred_cell{i,k,l};
%                 showOneLegend( plot(domainX, yPred(idxDomain),marks{i}),arr_samplesize{i} )
%                 showTitles(objFunc, str_setting, k, l);
% 
%                 if(minMeanObj_over_funcs(k) > min(yPred))
%                     minMeanObj_over_funcs(k) = min(yPred);
%                 end
%                 if(maxMeanObj_over_funcs(k) < max(yPred))
%                     maxMeanObj_over_funcs(k) = max(yPred);
%                 end
%                 sp2.YLim = [minMeanObj_over_funcs(k), maxMeanObj_over_funcs(k)];
%                 
% 
%                 
%                 legend();
%             end

            
            figSummary =figure(idxSummary);
            sp = mysubplot(dim_func,1,k);
            if(show_exact_fmin)
                if(exist('objFunc','var')&& numel(objFunc)>0)
                    if isfield(objFunc,'optVal')
                        hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
                    end
                elseif(exist('opt_val_fn','var')&& numel(opt_val_fn)>0)
                    hideLegend(plot([ sttX,budget]', log([opt_val_fn, opt_val_fn]'),'k--' ));
                end
            else
                hideLegend(plot([ sttX,budget]', [0, 0]','k--' ));
            end
            
            if(typeSamplesize == TypeSampleSize.Fixed)
                samplesize = givenSamplesize;
            else
                samplesize = 1;
            end
            % For each evaluation, I need to check if 
            n_actually_repeated = 0;
            for j=1:dim_repeats
                if( numel(results_all{i,j,k,l})>0)
                    n_actually_repeated = n_actually_repeated+1;
                end
            end

            yyy = zeros(floor(budget/samplesize) - n_init_sample, 1);
            nfeas = zeros(floor(budget/samplesize) - n_init_sample, 1);
            fgap = zeros(floor(budget/samplesize) - n_init_sample, 1);
            for idx_value = n_init_sample+1:floor(budget/samplesize);
                summed = 0;
                num_feas = 0;
                list = [];
                fgaplist = [];
                for j=1:dim_repeats
                    if( numel(results_all{i,j,k,l})>0)
                        if( idx_value <= size(results_all{i,j,k,l}.feas_history,1) )
                            feasgap = max(results_all{i,j,k,l}.feasgap_history(idx_value,:));
                            isfeasible = results_all{i,j,k,l}.feas_history(idx_value);
                            val = results_all{i,j,k,l}.eval_f_at_xnew_history(idx_value);

                            if (isfeasible)
                                summed = summed + val;
                                list = [list;val];
                                num_feas = num_feas + isfeasible;
                            end
                            fgaplist = [fgaplist;feasgap];
                        end

                    end
                end
                    yyy(idx_value-n_init_sample,1) = median(list);
                    nfeas(idx_value-n_init_sample,1) = numel(list);
                    fgap(idx_value-n_init_sample,1) = mean(fgaplist);
            end

            idx_feas = find(nfeas>0);
            xxx_orig = n_init_sample+1:floor(budget/samplesize);
            xxx = xxx_orig(idx_feas);
            yyy = yyy(idx_feas);

            if ~exist('yall','var');
                yall = yyy;
            else
                yall = [yall;yyy];
            end
            hold on;
            
            if contains(class(param_BO.objFunc),'Graphene')
                zzz = movmedian(yyy,50);                        
            else
                zzz=yyy; 
            end
            eb = plot([xxx]'*samplesize, zzz, marks{dim_cases_samplesize*(l-1)+i}, 'LineWidth',2, 'MarkerSize',3);
%             eb = plot([xxx]'*samplesize, yyy, marks{dim_cases_samplesize*(l-1)+i}, 'LineWidth',2, 'MarkerSize',3);
            showOneLegend(eb, sprintf('%s',  str_setting));
            if(numel(yall)>1 && min(yall)< quantile(yall,.90))
                ylim([min(yall), quantile(yall,.90)] )
            end
            
            if(l==num_settings);
                legend('Location','bestoutside');
            end

            showTitles(optProb);


            %% FIGURE 3

            figCon3 = figure(idxSumCon3);
            sp13 = mysubplot(dim_func,1,k);

            hold on;
            if contains(class(param_BO.objFunc),'Graphene')
                zzz = movmedian(fgap,50);                        
%                             zzz=fgap; 
            else
                zzz = fgap;
            end

            eb = plot(xxx_orig'*samplesize, zzz , marks{dim_cases_samplesize*(l-1)+i}, 'LineWidth',2, 'MarkerSize',3);
%             eb = plot(xxx_orig'*samplesize, fgap , marks{dim_cases_samplesize*(l-1)+i}, 'LineWidth',2, 'MarkerSize',3);
            showOneLegend( eb , sprintf('%s', str_setting));


            showTitles(optProb);

            if(l==num_settings);
                legend('Location','bestoutside');
            end





        end

    end
end


% fig = figure(idxSummary);
% annotateTextTitle(fig, strrep(A_RUN_PURPOSE,'_',' '));

for l=1:num_settings
	for k=1:dim_func
		if(optProb.nCon>0)
            if(showIndividual)
                if(sum(isinf(minObj_over_funcs(k))) + sum(isinf(maxObj_over_funcs(k)))==0)
                    figCon1 = figure(idxIndCon1);
                    sp = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));
                    sp.YLim = [minObj_over_funcs(k), maxObj_over_funcs(k)];
                end
            end
			
            if(sum(isinf(minObj_over_funcs(k))) + sum(isinf(maxObj_over_funcs(k)))==0)
                figCon1 = figure(idxSumCon1);
                sp = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));
                sp.YLim = [minMeanObj_over_funcs(k), maxMeanObj_over_funcs(k)];
            end

		else
            if(showIndividual)
                if(sum(isinf(minObj_over_funcs(k))) + sum(isinf(maxObj_over_funcs(k)))==0)
                    figure(idxIndUncon1);
                    sp = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));
                    sp.YLim = [minObj_over_funcs(k), maxObj_over_funcs(k)];
                end
            end
			
            if(sum(isinf(minObj_over_funcs(k))) + sum(isinf(maxObj_over_funcs(k)))==0)
                figure(idxSumUncon1);
                sp = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));
                sp.YLim = [minMeanObj_over_funcs(k), maxMeanObj_over_funcs(k)];
            end
		end

	end
end

% - - - ANNOTATION - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    nRow = dim_func;
    nCol = num_settings;
    strRows = cell(nRow,1);
    strCols = cell(nCol,1);
    for idx_row=1:nRow;
        strRows{idx_row} = cell_obj_func{idx_row}.name;
    end

    if(exist('figCon3','var'));
%         strTitle = sprintf('Feasibility Gap');
%         annotateTextOnAxes(figCon3, nRow, nCol, strRows, strCols, strTitle)        
        figCon3.Name = 'Constrainted Optimization: Feasibility Gap';
    end;
    if(exist('figSummary','var'));
        figSummary.Name = 'Optimization: Best Feasible Objective Value';
    end;


        

function showTitles(optProb)
%     if(exist('optProb','var')&& numel(optProb)>0)
        title(sprintf('%s',optProb.name ))
%     end
end
function [out, extended] = extendifshort(array, budget)
    if( size(array,1)<budget)
        out = zeros(budget,size(array,2));
        out(1:size(array,1),:) = array(1:end,:);
        extended = 1;
    else
        out = array;
        extended = 0;
    end
end

function y_cons = conFuncWrapper(fn, simulator, xnew, samplesize)
    out = zeros(size(xnew,1),1);
    for i=1:size(xnew,1)
        [~, y_cons, ~] = MultipleSamplerObjCon.evaluate(simulator, xnew(i,:), samplesize);
        out(i) = y_cons;
    end
end