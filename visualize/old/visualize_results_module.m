isFirstRun = true;
% isFirstRun = false;
% summarize = false;
showIndividual=true;
% showIndividual=false;
% showDetail = true;
showDetail = false;
clc;
% clearvars -except results_all;
    close all;
    
if(isFirstRun)
%     load save/result_210309
%     load result_210310_power_0_3.mat
%     load save/results_210312_rep100_power1_budget500.mat
%     load save/result_210315_rep_100_power_0_3_budget_1000
%     load save/intermediate_result_210315_rep50_100_budget1000_power1.mat
%     load save/result_210317_random_ls.mat
    

    [dim_cases_samplesize, dim_repeats, dim_func, num_settings] = size(results_all);
    
    results_summary = struct();
    % dim_func=1;
    % dim_func=1;
    for l = 1:num_settings
        for k=1:dim_func;
            for j=1:dim_repeats
                for i=1:dim_cases_samplesize
                    result = results_all{i,j,k,l};
                    if ([0,0]==size(result))
                        return;
                    end
        %             results_summary = putin_array_3args(result, results_summary, i,j,k, dim_func, dim2, dim3);
                    results_summary = putin_array_4dim(result, results_summary, i,j,k,l, dim_cases_samplesize, dim_repeats, dim_func, num_settings);


        %             result
                end
            end
        end
    end
    
    arr_samplesize = cell(numel(type_cell_samplesize),1);
    for i=1:numel(type_cell_samplesize)
        sel_type_samplesize = type_cell_samplesize{i};
        switch sel_type_samplesize{1}
            case TypeSampleSize.Adjust
                arr_samplesize{i} = 'adj';
            case TypeSampleSize.Fixed
                arr_samplesize{i} = num2str(sel_type_samplesize{2});
            otherwise
        end
    end
    
end
% return;
% fig = figure(80);clf;
% fig.Position = [500 0 1500 1000]; % For 3x3
% fig.Position = [100 0 2000 1000]; % For 3 x 4
% fig.Position = [100 100 1500 400]; % For 1 x 3

if(showDetail)
    subplot_frame = [dim_cases_samplesize,4];
else
    subplot_frame = [dim_cases_samplesize,num_settings*dim_func];
end
idx_subfigure_row = 1;


%%
%%
if(~showDetail)
    fig1 = figure(1);clf;
    fig2 = figure(2);clf;
    fig3 = figure(3);clf;
end
    % visualize_online_BO_2d_acq_search
for l=1:num_settings
    setting = type_cell_settings{l};
    if(iscell(setting))
        str_setting = sprintf('%s,',setting{:} );
    else
        str_setting = sprintf('%s,',setting );
    end
    for k=1:dim_func;

        merged_results = cell(dim_cases_samplesize,1);
        fnames = fieldnames(results_summary);
        dim_fields = max(size(fnames));
        for i=1:dim_cases_samplesize
            for h= 1:dim_fields
                fname = fnames{h};
                if(strcmp(fname,'nEvaluations'))
                    disp('');
                end
                sel_field = getfield(results_summary,fname);
                if(~isa(sel_field,'cell'))
                    continue;
                end
                if(~isnumeric(sel_field{1}))
                    continue;
                end
                if(size(sel_field{1,1},2)>1)
                    continue;
                end

%                 maxdim = max(cell2mat(results_summary.nEvaluations(i,:,k,l)));
                zz = getfield(results_summary,fname);
                ll = zz(i,:,k,l);
                nRows = numel(ll);


                isSame= true;
                maxdim = 1;
                for m=1:numel(ll)
                    if(size(ll{1},1) ~= size(ll{m},1))
                        isSame = false;
                    end
                    maxdim = max(maxdim, size(ll{m},1));
                end

                if(isSame)
                    str_cmd = sprintf('merged_results{i}.%s = cell2mat(results_summary.%s(i,:,k,l));',fname, fname);
                else
                    tempRes = zeros(maxdim, nRows );
                    for m=1:numel(ll)
                        tempRes(1:size(ll{m},1),m) = ll{m};
                    end
                    str_cmd = sprintf('merged_results{i}.%s = tempRes;',fname);
                end
        %         merged_results.y_history = cell2mat(results_summary.y_history(i,:,k))
                eval(str_cmd);
            end
    %         merged_results{i}.mean_eval_f_at_xnew_history = mean(merged_results{i}.eval_f_at_xnew_history,2);
    %         subplot(dim_cases_samplesize,1,1);
    %         
        end
%% ----------------------------------------------------------------------
        if(showDetail)
            fig1 = figure(50*(l-1)+k);clf;
        end
        for j=1:dim_repeats
            for i=1:dim_cases_samplesize
                result = results_all{i,j,k,l};
                if ([0,0]==size(result))
                    return;
                end
    %             subplot(dim_cases_samplesize,1,i);
                x_history = result.x_history;
                y_history = result.y_history;
                eval_f_at_xnew_history = result.eval_f_at_xnew_history;
                samplesize_history = result.samplesize_history;
                fmin_history = result.yMuMinHistory;
                xmin_history = result.xMuMinHistory;
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


                if(~showDetail)
                    if(showIndividual)
                        figure(1);
                        subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, i, fnGetIdxOfSubplotWithRowCol([dim_func,num_settings],k,l)));
                        hold on;
                        cum_samplesize =  cumsum(samplesize_history(n_init_sample+1:nEvaluations,1));

                        diff_f_min = (  eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1)  - log(opt_val_fn ));
                        showOneLegend(    plot(  cum_samplesize, diff_f_min , 'r-', [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), diff_f_min(end), 'ob') ,'Diff(y)');
    %                     title(sprintf('%s(%s)', objFunc.name, str_setting  ))


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


                else

                    idx_subfigure_row = i;
                    subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idx_subfigure_row, 1));
                    hold on;
                    cum_samplesize =  cumsum(samplesize_history(n_init_sample+1:nEvaluations,1));

                    diff_f_min = (  eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1)  - log(opt_val_fn ));
    %                     idxNInf = find(eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1) == -inf);
    %                     minDiff = min(diff_f_min(find(eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1) > -inf)));
                        showOneLegend(    plot(  cum_samplesize, diff_f_min , 'r-', [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), diff_f_min(end), 'ob') ,'Diff(y)');
                        title(sprintf('|f_{min} - min(mu)|: %.3g',diff_f_min(end)))

                        subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idx_subfigure_row, 2));
                    hold on;
                    %     diff_x_f_min = vecnorm(xmin_history(n_init_sample+1:nEvaluations,:) - opt_sol_fn, 2, 2);
                    %     diff_x_f_min = min(vecnorm(xmin_history(n_init_sample+1:nEvaluations,:) - opt_sol_fn, 2, 2));
                        diff_x_f_min = getMinimumDistance(xmin_history(n_init_sample+1:nEvaluations,:), opt_sol_fn);
                        plot( cum_samplesize, diff_x_f_min, 'r-' , [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), diff_x_f_min(end), 'ob');
                        title(sprintf('|x_{min} - argmin(mu)|: %.3g',diff_x_f_min(end)))
                    % subplot(subplot_frame(1),subplot_frame(2),9);
                    % idxFigStatus = 3;
                    subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idx_subfigure_row, 3));
                    hold on;
                    plot( cum_samplesize, samplesize_history(n_init_sample+1:nEvaluations,1), 'r-' , [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), samplesize_history(nEvaluations,1), 'ob');
                    title(sprintf('Sample Size: %d',samplesize_history(nEvaluations,1) ))

                    subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idx_subfigure_row, 4));

                    hold on;
                    plot( cum_samplesize, iterElapsedTime(n_init_sample+1:nEvaluations,1), 'r-' , [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), iterElapsedTime(nEvaluations,1), 'ob');
                    title(sprintf('Elapsed Time for iteration: %.1f sec',iterElapsedTime(nEvaluations,1)))
                end


            end
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
    % for k = 1:dim_func
        sttX = sum(merged_results{i}.samplesize_history(1:merged_results{i}.n_init_sample(1)));
        domainX = [sttX:budget]';
        dim_all  = zeros(dim_cases_samplesize,1);
        yys_cell = cell(dim_cases_samplesize,1);
        xxs_cell = cell(dim_cases_samplesize,1);
        for i = 1:dim_cases_samplesize
            dim_all(i) = sum(merged_results{i}.nEvaluations);
            yys_cell{i} = zeros(dim_all(i),1);
            xxs_cell{i} = zeros(dim_all(i),1);
        end

        ypred_cell = cell(dim_cases_samplesize,1);
        for i = 1:dim_cases_samplesize
            cnt = 1;
            cnt_new = 0;
            for j = 1:10
                cnt_new = cnt_new + merged_results{i}.nEvaluations(j) - merged_results{i}.n_init_sample(j);
        %         cnt_new = merged_results{i}.n_init_sample(j);
                yys_cell{i}(cnt:cnt_new)= ...
                    merged_results{i}.eval_f_at_xnew_history(merged_results{i}.n_init_sample(j)+1:merged_results{i}.nEvaluations(j),j);
        %         xx=cumsum(merged_results{1}.samplesize_history(merged_results{1}.n_init_sample:merged_results{1}.nEvaluations+1,j))
                xxs_cell{i}(cnt:cnt_new)= ...
                    cumsum(merged_results{i}.samplesize_history(merged_results{i}.n_init_sample(j)+1:merged_results{i}.nEvaluations(j),j))  + sum(merged_results{i}.samplesize_history(1:merged_results{i}.n_init_sample(j)));
                cnt = cnt_new + 1;
            end
            xxs_cell{i} = xxs_cell{i}(1:cnt-1);
            yys_cell{i} = yys_cell{i}(1:cnt-1);
%             gprMd = fitrgp(xxs_cell{i},yys_cell{i},'PredictMethod', 'exact');
%             [ypred,ysd,yint] = gprMd.predict(domainX);
            ypred = smooth(xxs_cell{i},yys_cell{i},.1,'loess');
            ypred_cell{i} = ypred;
        end


        arr_samplesize = {'adj','3','10','50'};
        % names = {'adj',}
%         marks = {'r-','b-','c-','g-'};
        marks = {'r-','b-','c-','g-','k-','m-','y-', 'r--','b--','c--','g--','k--','m--','y--', 'r:','b:','c:','g:','k:','m:','y:'};
        upperBound = cum_samplesize;
%% ----------------------------------------------------------------------
        if(~showDetail)
            fig2 =figure(2);
            subplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));
            hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
            for i = 1:dim_cases_samplesize
                hold on;
                [domainX, idxDomain] = sort(xxs_cell{i});
                yPred = ypred_cell{i};
                showOneLegend( plot(domainX, yPred(idxDomain),marks{i}),arr_samplesize{i} )
            end
            if(k==1&&l==1)
                title(sprintf('%s(%s)',objFunc.name, str_setting ))
            elseif(l==1)
                title(sprintf('%s',objFunc.name ))
            elseif(k==1)
                title(sprintf('%s', str_setting ))
            else
                title(sprintf('%s(%s)',objFunc.name, str_setting ))
            end
            
            legend()
%             annotation(fig2,'textbox',...
%                     [0.0125 0.934714286920571 0.808928568714431 0.0642857130794299],...
%                     'LineStyle','none',...
%                     'FontWeight','bold',...
%                     'String',{sprintf('%s (%s) power = %.2g',objFunc.name, str_setting , power )});

            ;
            % names = {'adj',}
    %% ----------------------------------------------------------------------
    %         fig3 =figure(50*(l-1)+20+k);clf;
            fig3 =figure(3);
            subplot(dim_func,num_settings,fnGetIdxOfSubplotWithRowCol([dim_func,num_settings], k, l));
            hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
            for i = 2:dim_cases_samplesize
                hold on;
                res = mean(merged_results{i}.eval_f_at_xnew_history,2);
                yy = res(merged_results{i}.n_init_sample(1)+1:merged_results{i}.nEvaluations(1));
                xx = cumsum(merged_results{i}.samplesize_history(merged_results{i}.n_init_sample(1)+1:merged_results{i}.nEvaluations(1),1)) +  sum(merged_results{i}.samplesize_history(1:merged_results{i}.n_init_sample(1)));
                showOneLegend( plot( xx,yy,marks{i}), arr_samplesize{i});
            end
            if(k==1&&l==1)
                title(sprintf('%s(%s)',objFunc.name, str_setting ))
            elseif(l==1)
                title(sprintf('%s',objFunc.name ))
            elseif(k==1)
                title(sprintf('%s', str_setting ))
            else
                title(sprintf('%s(%s)',objFunc.name, str_setting ))
            end
            legend();
            
            
    %% ----------------------------------------------------------------------
            fig4 =figure(4);
            subplot(dim_func,1,k);
            hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
            for i = 1:dim_cases_samplesize
                hold on;
                [domainX, idxDomain] = sort(xxs_cell{i});
                yPred = ypred_cell{i};
                showOneLegend( plot(domainX, yPred(idxDomain),marks{dim_cases_samplesize*(l-1)+i}),sprintf('%s:%s', arr_samplesize{i}, str_setting) )
            end
            if(k==1&&l==1)
                title(sprintf('%s(%s)',objFunc.name, str_setting ))
            elseif(l==1)
                title(sprintf('%s',objFunc.name ))
            elseif(k==1)
                title(sprintf('%s', str_setting ))
            else
                title(sprintf('%s(%s)',objFunc.name, str_setting ))
            end
            
            legend()
            
        else
            fig2 =figure(50*(l-1)+10+k);clf;
            subplot(1,2,1);
            hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
            for i = 1:dim_cases_samplesize
                hold on;
                showOneLegend( plot(xxs_cell{i}, ypred_cell{i},marks{i}),arr_samplesize{i} )
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
            subplot(1,2,2);
            hideLegend(plot([ sttX,budget]', log([objFunc.optVal, objFunc.optVal]'),'k--' ));
            for i = 2:dim_cases_samplesize
                hold on;
                res = mean(merged_results{i}.eval_f_at_xnew_history,2);
                yy = res(merged_results{i}.n_init_sample(1)+1:merged_results{i}.nEvaluations(1));
                xx = cumsum(merged_results{i}.samplesize_history(merged_results{i}.n_init_sample(1)+1:merged_results{i}.nEvaluations(1),1)) +  sum(merged_results{i}.samplesize_history(1:merged_results{i}.n_init_sample(1)));
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
