isFirstRun = true;
% isFirstRun = false;
% summarize = false;
% showIndividual=true;
showIndividual=false;
% showDetail = true;
showDetail = false;
clc;
% clearvars -except results_all;
    close all;
% mysubplot = @(x1,x2,x3) subtightplot(x1,x2,x3, [.06, 0.02]);
mysubplot = @(x1,x2,x3) subplot(x1,x2,x3);
    
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
    load results_210323_constraints_trial;
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
if(~showDetail)
    fig1 = figure(1);clf;
    fig2 = figure(2);clf;
    fig3 = figure(3);clf;
    fig11 = figure(11);clf;
    fig12 = figure(12);clf;
end
    % visualize_online_BO_2d_acq_search
    
predicted_cell = cell(dim_cases_samplesize, dim_func, num_settings);

dim_all  = zeros(dim_cases_samplesize, dim_func, num_settings);
yys_cell = cell(dim_cases_samplesize, dim_func, num_settings);
xxs_cell = cell(dim_cases_samplesize, dim_func, num_settings);
ypred_cell = cell(dim_cases_samplesize, dim_func, num_settings);

if(param_BO.hasConstraints)
    % xxs_cell{i,k,l}
    feasible_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
    optFeasible_cell = cell(dim_cases_samplesize, dim_func, num_settings);
    feasgap_cell     = cell(dim_cases_samplesize, dim_func, num_settings);
    probCons_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
    xxsFeasible_cell  = cell(dim_cases_samplesize, dim_func, num_settings);

    feasible_pred_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
    xxsFeasible_pred_cell = cell(dim_cases_samplesize, dim_func, num_settings);
    optFeasible_pred_cell = cell(dim_cases_samplesize, dim_func, num_settings);
    feasgap_pred_cell     = cell(dim_cases_samplesize, dim_func, num_settings);
    probCons_pred_cell    = cell(dim_cases_samplesize, dim_func, num_settings);
end
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
% %% ----------------------------------------------------------------------
        if(showDetail)
            fig1 = figure(50*(l-1)+k);clf;
        end
%% Loop: Sample sizes
        for i=1:dim_cases_samplesize
            idx_pred_stt = 1;            idx_pred_end = 0;
            idx_feas_stt = 1;            idx_feas_end = 0;
            %% Loop: Repeats
            for j=1:dim_repeats
                result = results_all{i,j,k,l};
                if ([0,0]==size(result))
                    return;
                end
    %             mysubplot(dim_cases_samplesize,1,i);
                x_history = result.x_history;
                y_history = result.y_history;
                eval_f_at_xnew_history = result.eval_f_at_xnew_history;
                samplesize_history = result.samplesize_history;
                y_mu_history = result.yMuMinHistory;
                x_mu_history = result.xMuMinHistory;
                hist_beta_t = result.hist_beta_t;
                nEvaluations = result.nEvaluations;
                n_init_sample = result.n_init_sample;
                param_BO = result.param_BO;
                visualizer = result.visualizer;
                budget = result.budget;
                givenSamplesize = result.givenSamplesize;
                typeSamplesize = result.typeSamplesize;
                elapsedTime = result.elapsedTime;
                target_function = result.objFunc;
                iterElapsedTime = result.iterElapsedTime;
%                 scr_true_functions;
                objFunc = result.objFunc;
                opt_sol_fn = result.opt_sol_fn;
                opt_val_fn = result.opt_val_fn;
                power = result.power;
                
%                 scr_setResultValuesToVariables
                
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

                if(param_BO.hasConstraints)

                    % ----- Assign Calculated values --------------------------------------
                    cum_samplesize =  cumsum(samplesize_history(n_init_sample+1:nEvaluations,1));
                    % 1. feasible
                    [feasible, strCon] =  isFeasible(constraints.constraints{1}, meanFunc.fnEval , x_mu_history(n_init_sample+1:nEvaluations,:));
    %                 cum_samplesize, feasible
                    feasible_cell{i,k,l}(idx_pred_stt:idx_pred_end)  = feasible;
                    % 2. optimal among feasible (val(idxFeas))
                    idxFeas = find(feasible);
    %                 cum_samplesize(idxFeas), val(idxFeas)
                    idx_feas_end = idx_feas_end + numel(idxFeas);
                    xxsFeasible_cell{i,k,l}(idx_feas_stt:idx_feas_end) = cum_samplesize(idxFeas);
                        if(exist('objFunc','var')&& numel(objFunc)>0)
                            opt_val_fn= objFunc.optVal;
                            show_exact_fmin = false;
                        elseif(exist('opt_val_fn','var')&& numel(opt_val_fn)>0)
            %                 opt_val_fn = opt_val_fn;
                            show_exact_fmin = false;
                        else
                            show_exact_fmin = true;
                        end
                        if(show_exact_fmin)
                            eval_f =   eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1) ;
    %                         showOneLegend(    plot(  cum_samplesize, eval_f , 'r-',  [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  cum_samplesize(end), eval_f(end), 'ob') ,'Diff(y)');
%                             title(sprintf('f_{min} : %.3g',eval_f(end)))
                            val = eval_f;
    %                         subp15.YLim(2) = max(eval_f)+1;
                        else
                            diff_f_min = (  eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1)  - log(objFunc.getOptVal() ));
    %                         showOneLegend(    plot(  cum_samplesize, diff_f_min , 'r-',  [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  cum_samplesize(end), diff_f_min(end), 'ob') ,'Diff(y)');
%                             title(sprintf('|f_{min} - min(mu)|: %.3g',diff_f_min(end)))
                            val = diff_f_min;
    %                         subp15.YLim(1) = 0;
                        end
                    optFeasible_cell{i,k,l}(idx_feas_stt:idx_feas_end)= val(idxFeas);

                    % 3. feasiblity gap
                    feasgap = feasibilityGap(constraints.constraints{1}, meanFunc.fnEval, x_mu_history(n_init_sample+1:nEvaluations,:));
    %                 cum_samplesize, feasgap
                    feasgap_cell{i,k,l}(idx_pred_stt:idx_pred_end)= feasgap;


                    % 4. probability of the constraints being satisfied
                    [~, pCons] = funcFeas(x_mu_history(n_init_sample+1:nEvaluations,:), 0);
    %                 cum_samplesize, pCons
                    probCons_cell{i,k,l}(idx_pred_stt:idx_pred_end)= pCons;
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

                    feasgap_cell{i,k,l} = feasgap_cell{i,k,l}(1:idx_pred_stt-1);
                    probCons_cell{i,k,l} = probCons_cell{i,k,l}(1:idx_pred_stt-1);
                end

                ypred = smooth(xxs_cell{i,k,l},yys_cell{i,k,l},.1,'loess');
                ypred_cell{i,k,l} = ypred;
               

                if(param_BO.hasConstraints)
                    feasible_pred_cell{i,k,l} = smooth(xxs_cell{i,k,l},feasible_cell{i,k,l},.1,'loess');
    %                 xxsFeasible_pred_cell  {i,k,l} = smooth(xxs_cell{i,k,l},yys_cell{i,k,l},.1,'loess');
                    optFeasible_pred_cell{i,k,l} = smooth(xxsFeasible_cell{i,k,l},optFeasible_cell{i,k,l},.1,'loess');
                    feasgap_pred_cell{i,k,l} = smooth(xxs_cell{i,k,l},feasgap_cell{i,k,l},.1,'loess');
                    probCons_pred_cell{i,k,l} = smooth(xxs_cell{i,k,l},probCons_cell{i,k,l},.1,'loess');
                end
                
                
%                 showIndividual
                if(showIndividual)
%                     figure(1);
%                     mysubplot(mysubplot_frame(1),mysubplot_frame(2),fnGetIdxOfSubplotWithRowCol(mysubplot_frame, i, fnGetIdxOfSubplotWithRowCol([dim_func,num_settings],k,l)));
%                     hold on;
%                     cum_samplesize =  cumsum(samplesize_history(n_init_sample+1:nEvaluations,1));
% 
%                     diff_f_min = (  eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1)  - log(opt_val_fn ));
%                     showOneLegend(    plot(  cum_samplesize, diff_f_min , 'r-', [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), diff_f_min(end), 'ob') ,'Diff(y)');
%                     showTitles(objFunc, str_setting, k, l);
                    
                    
                    if(param_BO.hasConstraints)
                        figure(11);
                        mysubplot(mysubplot_frame(1),mysubplot_frame(2),fnGetIdxOfSubplotWithRowCol(mysubplot_frame, i, fnGetIdxOfSubplotWithRowCol([dim_func,num_settings],k,l)));
                        hold on;
                        % cum_samplesize =  cumsum(samplesize_history(n_init_sample+1:nEvaluations,1));
                        cum_samplesize =  cumsum(samplesize_history(n_init_sample+1:nEvaluations,1)) + sum(samplesize_history(1:n_init_sample,1));

                        xLimLeft = sum(samplesize_history(1:n_init_sample+1,1));
                        xLimRight = sum(samplesize_history)+1;



                        if(show_exact_fmin)
                            eval_f =   eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1) ;
                            showOneLegend(    plot(  cum_samplesize, eval_f , 'r-',  [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  cum_samplesize(end), eval_f(end), 'ob') ,'Diff(y)');
                            title(sprintf('f_{min} : %.3g',eval_f(end)))
                            val = eval_f;
                            subp15.YLim(2) = max(eval_f)+1;
                            strTitleObjFunc = sprintf('f_{min}' );
                        else
                            diff_f_min = (  eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1)  - log(objFunc.getOptVal() ));
                            showOneLegend(    plot(  cum_samplesize, diff_f_min , 'r-',  [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  cum_samplesize(end), diff_f_min(end), 'ob') ,'Diff(y)');
                            title(sprintf('|f_{min} - min(mu)|: %.3g',diff_f_min(end)))
                            val = diff_f_min;
                            subp15.YLim(1) = 0;
                            strTitleObjFunc = sprintf('|f_{min} - min(mu)|' );
                        end
                        % show feasibility
                        hold on;

                        % 1. feasible
                        [feasible, strCon] =  isFeasible(constraints.constraints{1}, meanFunc.fnEval , x_mu_history(n_init_sample+1:nEvaluations,:));
                        idxFeas = find(feasible);
                        idxInfeas = find(1-feasible);
                        % hideLegend( plot(  cum_samplesize, val, '-r'));
                        hideLegend( plot(  cum_samplesize(idxInfeas), val(idxInfeas), 'xr', 'MarkerSize',12, 'LineWidth',3 ));
                        % 2. optimal among feasible (val(idxFeas))
                        hideLegend( plot(  cum_samplesize(idxFeas), val(idxFeas), 'og', 'MarkerSize',12, 'LineWidth',3) );

                        subp15.XLim(1) = xLimLeft;
                        subp15.XLim(2) = xLimRight;



                        % 3. feasiblity gap
                        feasgap = feasibilityGap(constraints.constraints{1}, meanFunc.fnEval, x_mu_history(n_init_sample+1:nEvaluations,:));
                        showOneLegend( plot(  cum_samplesize, feasgap, '-g') , 'Feas Gap' );    
                        [feasible, strCon] =  isFeasible(constraints.constraints{1}, meanFunc.fnEval , x_mu_history(n_init_sample+1:nEvaluations,:));
                        idxInfeas = find(1-feasible);
                        hideLegend( plot(  cum_samplesize(idxInfeas), feasgap(idxInfeas), 'xr', 'MarkerSize',10, 'LineWidth',3) );    
                        hideLegend( plot(  cum_samplesize(idxFeas), feasgap(idxFeas), 'og', 'MarkerSize',12, 'LineWidth',3) );

                        % 4. probability of the constraints being satisfied
                        [~, pCons] = funcFeas(x_mu_history(n_init_sample+1:nEvaluations,:), 0);
                    %     feasgap = feasibilityGap(cons{1}, meanFunc.fnEval, x_mu_history(n_init_sample+1:nEvaluations,:));
                        showOneLegend( plot(  cum_samplesize, pCons, '-m'), 'P(C(x))' );    
                        hideLegend(plot([cum_samplesize(1),cum_samplesize(end)], [1,1],'k--'));
                        hideLegend(plot([cum_samplesize(1),cum_samplesize(end)], [0,0],'k--'));
                        xLimLeft = sum(samplesize_history(1:n_init_sample+1,1));
                        xLimRight = sum(samplesize_history)+1;
                        subp25.XLim = [xLimLeft, xLimRight];
                        subp25.YLim = [-.05, 1.05];

                        xx = cum_samplesize(idxInfeas);
                    %     tt = text( xx , 0.1.*ones(size(xx)), num2str( floor(100*feasgap(find(xx)))/100));
                        tt = text( xx , 0.1.*ones(size(xx)), num2str( floor(100*feasgap(idxInfeas))/100), 'FontSize', 12, 'Color','r');

                        title('Feasibility');
    %                     legend();
                    end

                    
                end

            end


            arr_samplesize = {'adj','3','10','50'};
            % names = {'adj',}
            marks = {'r-','b-','c-','g-','k-','m-','y-', 'r--','b--','c--','g--','k--','m--','y--', 'r:','b:','c:','g:','k:','m:','y:'};
            upperBound = cum_samplesize;
%% ----------------------------------------------------------------------
            
            
            
            
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            if(param_BO.hasConstraints)
                    % 2. optimal among feasible (val(idxFeas))
                fig11 = figure(12);
                mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));

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
                showTitles(objFunc, str_setting, k, l);

                legend();
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                    % 1. feasible
                fig12 = figure(13);
                sp12 = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));

                    hideLegend(plot([ sttX,budget]', log([0, 0]'),'k--' ));
                    hideLegend(plot([ sttX,budget]', log([1, 1]'),'k--' ));
                hold on;
                yPred = feasible_pred_cell{i,k,l};
                showOneLegend( plot(domainX, yPred(idxDomain),marks{i}),arr_samplesize{i} )
                showTitles(objFunc, str_setting, k, l);

                legend();
                sp12.YLim = [-0.05 1.05];
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                    % 3. feasiblity gap
                fig13 = figure(14);
                sp13 = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));

                hold on;
                yPred = feasgap_pred_cell{i,k,l};
                showOneLegend( plot(domainX, yPred(idxDomain),marks{i}),arr_samplesize{i} )
                showTitles(objFunc, str_setting, k, l);

                legend();
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                    % 4. probability of the constraints being satisfied
                fig14 = figure(15);
                sp14 = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));

                    hideLegend(plot([ sttX,budget]', log([0, 0]'),'k--' ));
                    hideLegend(plot([ sttX,budget]', log([1, 1]'),'k--' ));
                hold on;
                yPred = probCons_pred_cell{i,k,l};
                showOneLegend( plot(domainX, yPred(idxDomain),marks{i}),arr_samplesize{i} )
                showTitles(objFunc, str_setting, k, l);

                legend();
                sp14.YLim = [-0.05, 1.05];
            else
                fig2 =figure(2);
                sp2 = mysubplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));
                
                if(exist('objFunc','var')&& numel(objFunc)>0)
                    opt_val_fn= objFunc.optVal;
                    show_exact_fmin = false;
                elseif(exist('opt_val_fn','var')&& numel(opt_val_fn)>0)
    %                 opt_val_fn = opt_val_fn;
                    show_exact_fmin = false;
                else
                    show_exact_fmin = true;
                end

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
%                 if(exist('objFunc','var')&& numel(objFunc)>0)
%                     hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
%                 elseif(exist('opt_val_fn','var')&& numel(opt_val_fn)>0)
%                     hideLegend(plot([ sttX,budget]', log([opt_val_fn, opt_val_fn]'),'k--' ));
%                 end
    %                 for i = 1:dim_cases_samplesize
                hold on;
                [domainX, idxDomain] = sort(xxs_cell{i,k,l});
                yPred = ypred_cell{i,k,l};
                showOneLegend( plot(domainX, yPred(idxDomain),marks{i}),arr_samplesize{i} )
                showTitles(objFunc, str_setting, k, l);

                legend();
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
            fig4 =figure(4);
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
%                 %     diff_x_f_min = vecnorm(x_mu_history(n_init_sample+1:nEvaluations,:) - opt_sol_fn, 2, 2);
%                 %     diff_x_f_min = min(vecnorm(x_mu_history(n_init_sample+1:nEvaluations,:) - opt_sol_fn, 2, 2));
%                     diff_x_f_min = getMinimumDistance(x_mu_history(n_init_sample+1:nEvaluations,:), opt_sol_fn);
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
            fig2 =figure(50*(l-1)+10+k);clf;
            mysubplot(1,2,1);
            hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
            for i = 1:dim_cases_samplesize
                hold on;
                showOneLegend( plot(domainX, ypred_cell{i},marks{i}),arr_samplesize{i} )
            end
    %         title(sprintf('(%s,%s)',objFunc.name, str_setting ))
            legend()
            annotation(fig2,'textbox',...
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
if(param_BO.hasConstraints)

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


    strTitle = sprintf('%s (Only Feasible)', strTitleObjFunc);
    annotateTextOnAxes(fig11, nRow, nCol, strRows, strCols, strTitle)        
    strTitle = sprintf('Fraction of Feasible evaluations');
    annotateTextOnAxes(fig12, nRow, nCol, strRows, strCols, strTitle)        
    strTitle = sprintf('Feasibility Gap');
    annotateTextOnAxes(fig13, nRow, nCol, strRows, strCols, strTitle)        
    strTitle = sprintf('P(C(x))' );
    annotateTextOnAxes(fig14, nRow, nCol, strRows, strCols, strTitle)        
end

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
    strTitle = sprintf('%s (power = %.2g)','|f_{min} - min(mu)|', power );
annotateTextOnAxes(fig1, nRow, nCol, strRows, strCols, strTitle)        
        
%     annotation(fig2,'textbox',...
%         [0.0125 0.934714286920571 0.808928568714431 0.0642857130794299],...
%         'LineStyle','none',...
%         'FontWeight','bold',...
%         'String',{sprintf('GP: |f_{min} - min(mu)| (power = %.2g)', power )});
    nRow = dim_func;
    nCol = num_settings;
    strRows = cell(nRow,1);
    strCols = cell(nCol,1);
    for i=1:nRow;
        strRows{i} = cell_obj_func{i}.name;
    end
    for i=1:nCol;
        setting = type_cell_settings{l};
        if(iscell(setting))
            str_setting = sprintf('%s,',setting{:} );
        else
            str_setting = sprintf('%s,',setting );
        end
        strCols{i} =  sprintf('%s\n',str_setting );
    end
%     strRows = {'r1','r2','r3','r4'};
%     strCols = {'c1','c2','c3'};
    strTitle = sprintf('GP: |f_{min} - min(mu)| (power = %.2g)', power );
annotateTextOnAxes(fig2, nRow, nCol, strRows, strCols, strTitle)        

%     annotation(fig3,'textbox',...
%         [0.0125 0.934714286920571 0.808928568714431 0.0642857130794299],...
%         'LineStyle','none',...
%         'FontWeight','bold',...
%         'String',{sprintf('|f_{min} - min(mu)| (power = %.2g)', power )});
%     nRow = 4;
%     nCol = 3;
%     strRows = {'r1','r2','r3','r4'};
%     strCols = {'c1','c2','c3'};
    strTitle = sprintf('|f_{min} - min(mu)| (power = %.2g)', power );
annotateTextOnAxes(fig3, nRow, nCol, strRows, strCols, strTitle)        

end
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