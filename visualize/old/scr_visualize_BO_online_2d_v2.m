show_xdiff = false;
sttVisual = tic ;
gridCoord = visualizer.getGridCoord();
gridSquareD1=visualizer.getGridSquareD1();
gridSquareD2=visualizer.getGridSquareD2();
gridDomainD1=visualizer.getGridDomainD1();
gridDomainD2= visualizer.getGridDomainD2();
% true_minX = visualizer.getTrueGlobalMinX();
true_fmin = visualizer.getTrueGlobalMinF();
true_minX = objFunc.getOptSol();
% TO VISUALIZE
% mysubplot = @(x1,x2,x3) mysubplot(x1,x2,x3);

mysubplot = @(x1,x2,x3) subtightplot(x1,x2,x3, [.06, 0.02]);
str_optimal= 'Optimal';
% str_min_mu = 'opt(mu)';
% str_max_acq = 'opt(acq)';
str_last = 'last';
subplot_frame = [3,5];
    if(infillOptimizer.isMaximizeAcq)
        str_max_acq = 'max(acq)';
    else
        str_max_acq = 'min(acq)';
    end
    if(infillOptimizer.isMaximizeObj)
        str_min_mu = 'max(mu)';
    else
        str_min_mu = 'min(mu)';
    end

%% VISUALIZE
fig = figure(80);clf;
% fig.Position = [500 0 1500 1000]; % For 3x3
fig.Position = [50 50 1850 950]; % For 3 x 4

color_lightgrey = [.7, .7, .7];
color_darkgrey = [.3, .3, .3];

if(~param_BO.isInfillSearchGrid)
    
    queryVals_mu = queries_hist.mu;
    queryVals_acq = queries_hist.acq;
    queries_all = queries_hist.all;
    history_imprv = queries_hist.improv;
    queryVals_sigma = queries_hist.sigma;

    if(infillOptimizer.isMaximizeAcq)
        [~, idx_max_acq] = max(queryVals_acq);
    else
        [~, idx_max_acq] = min(queryVals_acq);
    end
    if(infillOptimizer.isMaximizeObj)
        [~ , idx_min_mu] = max(queryVals_mu);
    else
        [~ , idx_min_mu] = min(queryVals_mu);
    end

    x_history = storage.xHistory;
    y_history = storage.yHistory;
    y_cons_history = storage.yConsHistory;
    samplesize_history = storage.samplesizeHistory;
    x_mu_history = storage.xMuMinHistory;
    eval_f_at_xnew_history = storage.evalFAtXnewHistory;
    xmin_history = storage.xMuMinHistory;

    % First Column : GRID
    if(visualizer.isShowAllGrid())
    %% ======================================================================================
    % ----- Evaluated points ---------------------------------------------------------------------------------
        idxFig = 1;
        subp13 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 1, 3));
        fn_true_mean=objFunc.getFnEval();
        flippedImagesc(gridDomainD1, gridDomainD2, (reshape(log(fn_true_mean(gridCoord(:,1), gridCoord(:,2))),visualizer.getDim(),visualizer.getDim()))');
        hold on;
        hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
        showOneLegend( scatter(true_minX(:,1), true_minX(:,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
%         showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
%         showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
%         showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);
            showOneLegend( scatter(x_mu_history(nEvaluations,1), x_mu_history(nEvaluations,2),150, 'x', 'LineWidth', 2, 'MarkerEdgeColor', color_darkgrey), 'Est. Opt' ); % Estimated Optimal
            showOneLegend( scatter(x_history(nEvaluations,1), x_history(nEvaluations,2),30, 'LineWidth', 2, 'MarkerEdgeColor', 'r','MarkerEdgeAlpha',.7), 'Eval X' ); % SELECTED POINT
        colorbar();
        title(sprintf('fn:%s,n=%d', objFunc.getName(), sum(samplesize_history(n_init_sample+1:nEvaluations,1))));
        xlabel('x');ylabel('y');
        % --------------------------------------------------------------------------------------

    % ======================================================================================
        % --------------------------------------------------------------------------------------
        idxFigGrid = 1;
        subp11 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 1, 1));

        [mu_est, sig_est] = exactPredictor.predict(gridCoord);
        flippedImagesc(gridDomainD1, gridDomainD2, (reshape(mu_est, visualizer.getDim(), visualizer.getDim()))');

        hold on;
        hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
        showOneLegend( scatter(true_minX(:,1), true_minX(:,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
%         showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
%         showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
%         showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);
            showOneLegend( scatter(x_mu_history(nEvaluations,1), x_mu_history(nEvaluations,2),150, 'x', 'LineWidth', 2, 'MarkerEdgeColor', color_darkgrey), 'Est. Opt' ); % Estimated Optimal
            showOneLegend( scatter(x_history(nEvaluations,1), x_history(nEvaluations,2),30, 'LineWidth', 2, 'MarkerEdgeColor', 'r','MarkerEdgeAlpha',.7), 'Eval X' ); % SELECTED POINT
        colorbar();
        title('Mean(log(Sigma)):\mu_{log\sigma}');
        xlabel('x');ylabel('y');
        % --------------------------------------------------------------------------------------
%         idx_col = 3;
        subp31 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 3, 1));
%         acq_est_grid = acqFunc.acquire(exactPredictor, gridCoord, iter);
        acq_est_grid = acquisition(gridCoord);

        flippedImagesc(gridDomainD1, gridDomainD2, (reshape(acq_est_grid, visualizer.getDim(), visualizer.getDim()))');

        hold on;
        hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
        showOneLegend( scatter(true_minX(:,1), true_minX(:,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
%         showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
%         showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
            showOneLegend( scatter(x_mu_history(nEvaluations,1), x_mu_history(nEvaluations,2),150, 'x', 'LineWidth', 2, 'MarkerEdgeColor', color_darkgrey), 'Est. Opt' ); % Estimated Optimal
            showOneLegend( scatter(x_history(nEvaluations,1), x_history(nEvaluations,2),30, 'LineWidth', 2, 'MarkerEdgeColor', 'r','MarkerEdgeAlpha',.7), 'Eval X' ); % SELECTED POINT
%         showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);
        colorbar();
        title(infillOptimizer.acqFuncName);
        xlabel('x');ylabel('y');
        % --------------------------------------------------------------------------------------
        idx_col = 2;
        subp21 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 2, 1));
        flippedImagesc(gridDomainD1, gridDomainD2, (reshape(sig_est, visualizer.getDim(), visualizer.getDim()))');

        hold on;
        hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
        showOneLegend( scatter(true_minX(:,1), true_minX(:,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
%         showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
%         showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
            showOneLegend( scatter(x_mu_history(nEvaluations,1), x_mu_history(nEvaluations,2),150, 'x', 'LineWidth', 2, 'MarkerEdgeColor', color_darkgrey), 'Est. Opt' ); % Estimated Optimal
            showOneLegend( scatter(x_history(nEvaluations,1), x_history(nEvaluations,2),30, 'LineWidth', 2, 'MarkerEdgeColor', 'r','MarkerEdgeAlpha',.7), 'Eval X' ); % SELECTED POINT
%         showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);
        colorbar();
        title('Sd(Sigma):s(\sigma)');
        xlabel('x');ylabel('y');
    end
    % ======================================================================================
        
        idxFig = 2;
        subp33 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 3, 3));
%         fn_eval_true=objFunc.getFnEval();
        hold on;
        hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),30, y_history(1:nEvaluations,:),'filled') );
        showOneLegend( scatter(true_minX(:,1), true_minX(:,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
%         showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
%         showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
            showOneLegend( scatter(x_mu_history(nEvaluations,1), x_mu_history(nEvaluations,2),150, 'x', 'LineWidth', 2, 'MarkerEdgeColor', color_darkgrey), 'Est. Opt' ); % Estimated Optimal
            showOneLegend( scatter(x_history(nEvaluations,1), x_history(nEvaluations,2),30, 'LineWidth', 2, 'MarkerEdgeColor', 'r','MarkerEdgeAlpha',.7), 'Eval X' ); % SELECTED POINT
        legend('Location','best' )
        colorbar();
        title(sprintf('True fn:%s,n=%d', objFunc.getName(), sum(samplesize_history(n_init_sample+1:nEvaluations,1))));
        xlabel('x');ylabel('y');
        % --------------------------------------------------------------------------------------
    % ======================================================================================
    if(visualizer.isShowSearchTrails())
    idxFigSearch = 2;
    % --------------------------------------------------------------------------------------
    % FIGURE 2-1
    idx_col = 1;
    subp12 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 1, 2));
    hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
    hold on;
    scatter(queries_all(:,1), queries_all(:,2),10, queryVals_mu,'filled');
    if(visualizer.isShowAllGrid())
        showOneLegend( scatter(true_minX(:,1), true_minX(:,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
    end
%     showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
%     showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
            showOneLegend( scatter(x_mu_history(nEvaluations,1), x_mu_history(nEvaluations,2),150, 'x', 'LineWidth', 2, 'MarkerEdgeColor', color_darkgrey), 'Est. Opt' ); % Estimated Optimal
            showOneLegend( scatter(x_history(nEvaluations,1), x_history(nEvaluations,2),30, 'LineWidth', 2, 'MarkerEdgeColor', 'r','MarkerEdgeAlpha',.7), 'Eval X' ); % SELECTED POINT
    colorbar();
    title('Mean');
    %                     contour(zz1, zz2, reshape(mu, dim, dim));
    % --------------------------------------------------------------------------------------
    % FIGURE 2-2
    idx_col = 3;
    subp32 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 3, 2));

    hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
    hold on;
    scatter(queries_all(:,1), queries_all(:,2),10, queryVals_acq,'filled');
    if(visualizer.isShowAllGrid())
        showOneLegend( scatter(true_minX(:,1), true_minX(:,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
    end
%     showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
%     showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);

    plot(queries_all(history_imprv(:,2),1), queries_all(history_imprv(:,2),2),'r-');
            showOneLegend( scatter(x_mu_history(nEvaluations,1), x_mu_history(nEvaluations,2),150, 'x', 'LineWidth', 2, 'MarkerEdgeColor', color_darkgrey), 'Est. Opt' ); % Estimated Optimal
            showOneLegend( scatter(x_history(nEvaluations,1), x_history(nEvaluations,2),30, 'LineWidth', 2, 'MarkerEdgeColor', 'r','MarkerEdgeAlpha',.7), 'Eval X' ); % SELECTED POINT

    colorbar();
    title(sprintf('%s',infillOptimizer.acqFuncName));
    % --------------------------------------------------------------------------------------
    % FIGURE 2-3
    idx_col = 2;
    % mysubplot(subplot_frame(1),subplot_frame(2),6);
    % idxFigSearch = 2;
    subp22 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 2, 2));

    hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
    hold on;
    scatter(queries_all(:,1), queries_all(:,2),10, queryVals_sigma,'filled');
    if(visualizer.isShowAllGrid())
        showOneLegend( scatter(true_minX(:,1), true_minX(:,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
    end
%     showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
%     showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
            showOneLegend( scatter(x_mu_history(nEvaluations,1), x_mu_history(nEvaluations,2),150, 'x', 'LineWidth', 2, 'MarkerEdgeColor', color_darkgrey), 'Est. Opt' ); % Estimated Optimal
            showOneLegend( scatter(x_history(nEvaluations,1), x_history(nEvaluations,2),30, 'LineWidth', 2, 'MarkerEdgeColor', 'r','MarkerEdgeAlpha',.7), 'Eval X' ); % SELECTED POINT
    colorbar();
    title('Sigma');
    % --------------------------------------------------------------------------------------


    end
% acquisition
end
% ======================================================================================
% mysubplot(subplot_frame(1),subplot_frame(2),7);
idxFigStatus = 3;
subp15 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 1, 5));
cla;
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
else
    diff_f_min = (  eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1)  - log(objFunc.getOptVal() ));
    showOneLegend(    plot(  cum_samplesize, diff_f_min , 'r-',  [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  cum_samplesize(end), diff_f_min(end), 'ob') ,'Diff(y)');
    title(sprintf('|f_{min} - min(mu)|: %.3g',diff_f_min(end)))
    val = diff_f_min;
    subp15.YLim(1) = 0;
end
% show feasibility
hold on;
[feasible, strCon] =  isFeasible(constraints.constraints{1}, meanFunc.fnEval , x_mu_history(n_init_sample+1:nEvaluations,:));
idxFeas = find(feasible);
idxInfeas = find(1-feasible);
% hideLegend( plot(  cum_samplesize, val, '-r'));
hideLegend( plot(  cum_samplesize(idxInfeas), val(idxInfeas), 'xr', 'MarkerSize',12, 'LineWidth',3 ));
hideLegend( plot(  cum_samplesize(idxFeas), val(idxFeas), 'og', 'MarkerSize',12, 'LineWidth',3) );

legend('Location','best');;


subp15.XLim(1) = xLimLeft;
subp15.XLim(2) = xLimRight;

% if(param_BO.hasConstraints)
%     feasgap = feasibilityGap(constraints.constraints{1}, meanFunc.fnEval, x_mu_history(n_init_sample+1:nEvaluations,:));
%     showOneLegend( plot(  cum_samplesize, feasgap, '-g') , 'Feas Gap' );    
%     hideLegend( plot(  cum_samplesize(idxInfeas), feasgap(idxInfeas), 'xr', 'MarkerSize',10, 'LineWidth',3) );    
% %     if(sum(feasgap)>10);
% %         disp(sum(feasgap));
% %     end
% %     feasgapXnew = feasibilityGap(constraints.constraints{1}, meanFunc.fnEval, x_history(nEvaluations,:));
% %     if(sum(feasgapXnew)>10);
% %         disp(sum(feasgapXnew));
% %     end
% end

% if(show_xdiff)
% subp25 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 1, 5));
%     diff_x_f_min = getMinimumDistance(x_mu_history(n_init_sample+1:nEvaluations,:), objFunc.getOptSol() );
%     plot( cum_samplesize, diff_x_f_min, 'r-' , [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--', cum_samplesize(end), diff_x_f_min(end), 'ob');
%     title(sprintf('|x_{min} - argmin(mu)|: %.3g',diff_x_f_min(end)))
%     subp25.XLim(1) = xLimLeft;
%     subp25.XLim(2) = xLimRight;
%     subp25.YLim(1) = 0;
% else
% end


%% -----------------------------------------------------
% cla;
if(param_BO.hasConstraints)
    subp25 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 2, 5));cla;
    hold on;
    feasgap = feasibilityGap(constraints.constraints{1}, meanFunc.fnEval, x_mu_history(n_init_sample+1:nEvaluations,:));
    showOneLegend( plot(  cum_samplesize, feasgap, '-g') , 'Feas Gap' );    
    [feasible, strCon] =  isFeasible(constraints.constraints{1}, meanFunc.fnEval , x_mu_history(n_init_sample+1:nEvaluations,:));
    idxInfeas = find(1-feasible);
    hideLegend( plot(  cum_samplesize(idxInfeas), feasgap(idxInfeas), 'xr', 'MarkerSize',10, 'LineWidth',3) );    
    hideLegend( plot(  cum_samplesize(idxFeas), feasgap(idxFeas), 'og', 'MarkerSize',12, 'LineWidth',3) );
% end
% if(param_BO.hasConstraints)
%     subplot(subp25);
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
    legend();
end

%% -----------------------------------------------------
    idxFigStatus = 3;
subp35 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 3, 5));
if(typeSamplesize == TypeSampleSize.Fixed)
    plot( cum_samplesize, samplesize_history(n_init_sample+1:nEvaluations,1), 'r-' , 'LineWidth',2);
%     text( 0,givenSamplesize, sprintf('%d',givenSamplesize));
    annotation(fig,'textbox',[0.5845 0.142 0.071 0.119],...
        'String',{sprintf('%d',givenSamplesize)},...
        'FontSize',60,...
        'LineStyle','none',...
        'FitBoxToText','off');
else
    plot( cum_samplesize, samplesize_history(n_init_sample+1:nEvaluations,1), 'r-' , [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  cum_samplesize(end), samplesize_history(nEvaluations,1), 'ob');
end
title(sprintf('Sample Size: %d (iter=%d)',samplesize_history(nEvaluations,1), iter))
subp35.XLim(1) = xLimLeft;
subp35.XLim(2) = xLimRight;
subp35.YLim(1) = 0;

% --------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------

if(param_BO.hasConstraints)
    cons = constraints.constraints;
    if(cons{1}.typeConstraint  == TypeConstraint.UnknownMeanConstraint)
%         if(cons{1}.hasUb && cons{1}.hasLb)
%             ub = cons{1}.ub;
%             lb = cons{1}.lb;
%         elseif(cons{1}.hasUb)
%             ub = cons{1}.ub;
%         elseif(cons{1}.hasLb)
%             lb = cons{1}.lb;
%         else
%             throwError('scr_visualize_BO_online... constraints bounds are not specified.')
%         end
%         meanFunc

%         if(param_BO.isInfillSearchGrid)
%             subp15 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 1, 5));
% %             flippedImagesc(gridDomainD1, gridDomainD2, (reshape(acq, dim, dim))');
% %             hold on;
% %             scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled');
% %             scatter(gridCoord(idx_min_mu,1), gridCoord(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[1 1 1]);
% %             scatter(gridCoord(maxIdx,1), gridCoord(maxIdx, 2), 20, maxV, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]);
% %             colorbar();
% %             title(sprintf('%s',param_BO.getStrTypeAcquisition()));
%         end
        if(visualizer.isShowAllGrid())
        % ----- Evaluated points ---------------------------------------------------------------------------------
            idxRow = 1;
            subp14 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 1, 4));
            cla;
            fn_true_mean=meanFunc.getFnEval();
            flippedImagesc(gridDomainD1, gridDomainD2, reshape(fn_true_mean(gridCoord(:,1), gridCoord(:,2)),visualizer.getDim,visualizer.getDim)');
            hold on;

%             if(cons{1}.hasUb && cons{1}.hasLb)
%                 strCon = sprintf('ub=%.2g, lb=%.2g', ub, lb);
%                 feas = ...
%                     (  reshape(fn_true_mean(gridCoord(:,1), gridCoord(:,2)),visualizer.getDim,visualizer.getDim)' >lb   ) .*...
%                     (  reshape(fn_true_mean(gridCoord(:,1), gridCoord(:,2)),visualizer.getDim,visualizer.getDim)' <ub  )  ;
%             elseif(cons{1}.hasUb)
%                 strCon = sprintf('ub=%.2g', ub);
%                 feas = ...
%                     (  reshape(fn_true_mean(gridCoord(:,1), gridCoord(:,2)),visualizer.getDim,visualizer.getDim)' <ub  )  ;
%             elseif(cons{1}.hasLb)
%                 strCon = sprintf('lb=%.2g', lb);
%                 feas = ...
%                     (  reshape(fn_true_mean(gridCoord(:,1), gridCoord(:,2)),visualizer.getDim,visualizer.getDim)' >lb   ) ;
%             else
%                 throwError('scr_visualize_BO_online... constraints bounds are not specified.')
%             end
            [feas_rows, strCon] =  isFeasible(cons{1}, meanFunc.fnEval, gridCoord(:,1:2));
            feas_eval = reshape(feas_rows,visualizer.getDim,visualizer.getDim)';
            boundaries_true = bwboundaries(feas_eval);
%             numberOfBoundaries = size(boundaries, 1);
%             for k = 1 : numberOfBoundaries
%                 thisBoundary = boundaries{k};
% %               plot(gridDomainD1(thisBoundary(:,2)), gridDomainD2(thisBoundary(:,1)), 'r', 'LineWidth', 2);
%                 h = fill(gridDomainD1(thisBoundary(:,2)), gridDomainD2(thisBoundary(:,1)), 'w', 'LineWidth', 2, 'EdgeColor', 'r');
%                 set(h,'facealpha',.2)
% %               scatter(thisBoundary(:,2), thisBoundary(:,1), 10, 'r', 'LineWidth', 2);
%             end
%             mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 2, 5));
%             boundaries = bwboundaries(feas);
%             numberOfBoundaries = size(boundaries, 1);
%             for k = 1 : numberOfBoundaries
%               thisBoundary = boundaries{k};
%               plot(thisBoundary(:,2), thisBoundary(:,1), 'r', 'LineWidth', 2);
%             end

%             reshape(fn_true_mean(gridCoord(:,1), gridCoord(:,2)),visualizer.getDim,visualizer.getDim)'
            hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
            showOneLegend( scatter(true_minX(:,1), true_minX(:,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
            showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
            showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
    %         showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);
            hideLegend( scatter(x_history(nEvaluations,1), x_history(nEvaluations,2),30, 'LineWidth', 2, 'MarkerEdgeColor', 'r') );
            colorbar();
            title(sprintf('True Cons: fn:%s,%s', meanFunc.getName(), strCon  ));
            xlabel('x');ylabel('y');
            % --------------------------------------------------------------------------------------
            
            subplot(subp14); hold on;
            drawBoundaries(boundaries_true, gridDomainD1, gridDomainD2, 'w', 'w', 0.5, 0.2)
            subplot(subp13); hold on;
            drawBoundaries(boundaries_true, gridDomainD1, gridDomainD2, 'w', 'w', 0.5, 0.2)


        end
        if(visualizer.isShowSearchTrails())
%------------------------------------------------------------------------------------------------
            subp24 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 2, 4)); cla;
            hold on;
%             hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),30, y_history(1:nEvaluations,:),'filled') );
            fn_true_mean=meanFunc.getFnEval();
            pred = multiPredictors.cell_predictors{1}.predict([gridCoord(:,1) gridCoord(:,2)]);
            % COMMENTED: Estimated GP shape
            flippedImagesc(gridDomainD1, gridDomainD2, reshape(pred,visualizer.getDim,visualizer.getDim)');
            
%             feas_eval = funcFeas([gridCoord(:,1) gridCoord(:,2)], urnd_feasiblity);
%             img=flippedImagesc(gridDomainD1, gridDomainD2, reshape(feas_eval,visualizer.getDim,visualizer.getDim)');
%             img.AlphaData = .2;
            
            hold on;
            
            feas_eval = funcFeas([gridCoord(:,1) gridCoord(:,2)], .05);
            feasible = reshape(feas_eval,visualizer.getDim,visualizer.getDim)' ==-1;
            boundaries_ub = bwboundaries(feasible);
            
            feas_eval = funcFeas([gridCoord(:,1) gridCoord(:,2)], .95);
            feasible = reshape(feas_eval,visualizer.getDim,visualizer.getDim)' ==-1;
            boundaries_lb = bwboundaries(feasible);
%------------------------------------------------------------------------------------------------
            showOneLegend( scatter(true_minX(:,1), true_minX(:,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
            hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
            showOneLegend( scatter(x_mu_history(nEvaluations,1), x_mu_history(nEvaluations,2),150, 'x', 'LineWidth', 2, 'MarkerEdgeColor', color_darkgrey), 'Est. Opt' ); % Estimated Optimal
            showOneLegend( scatter(x_history(nEvaluations,1), x_history(nEvaluations,2),30, 'LineWidth', 2, 'MarkerEdgeColor', 'r'), 'Eval X' ); % SELECTED POINT
            dom = objFunc.getXDomain;
            axis([dom(1,:), dom(2,:)]);
            colorbar();
            title(sprintf('Est Cons: fn:%s,%s (u=%.3g)', meanFunc.getName(), strCon, urnd_feasiblity  ));
            xlabel('x');ylabel('y');

            subplot(subp24); hold on;
            drawBoundaries(boundaries_true, gridDomainD1, gridDomainD2, 'w', 'w', 0.3, 0)
%             drawBoundaries(boundaries_ub, gridDomainD1, gridDomainD2, 'm', 'w', 0.3, 0)
%             drawBoundaries(boundaries_lb, gridDomainD1, gridDomainD2, 'm', 'w', 0.3, 0)
                drawBoundaries(boundaries_ub, gridDomainD1, gridDomainD2, 'm', 'w', 0.3, 0)
                drawBoundaries(boundaries_lb, gridDomainD1, gridDomainD2, 'g', 'w', 0.7, 0)

            if( type_infill_opt.strcmpi('DirectConstSampling'))
                [feas_eval, pCons] = funcFeas([gridCoord(:,1) gridCoord(:,2)], urnd_feasiblity);
                feasible = reshape(feas_eval,visualizer.getDim,visualizer.getDim)' ==-1;
                boundaries_est = bwboundaries(feasible);
                subplot(subp24); hold on;
                drawBoundaries(boundaries_est, gridDomainD1, gridDomainD2);

                subplot(subp33);hold on;
                drawBoundaries(boundaries_est, gridDomainD1, gridDomainD2);
                subplot(subp11);hold on;
                drawBoundaries(boundaries_est, gridDomainD1, gridDomainD2);
                subplot(subp21);hold on;
                drawBoundaries(boundaries_est, gridDomainD1, gridDomainD2);
                subplot(subp22);hold on;
                drawBoundaries(boundaries_est, gridDomainD1, gridDomainD2);
                subplot(subp31);hold on;
                drawBoundaries(boundaries_est, gridDomainD1, gridDomainD2);
                subplot(subp32);hold on;
                drawBoundaries(boundaries_est, gridDomainD1, gridDomainD2);
                subplot(subp12);hold on;
                drawBoundaries(boundaries_est, gridDomainD1, gridDomainD2);
            else
                    subplot(subp33);hold on;
                drawBoundaries(boundaries_ub, gridDomainD1, gridDomainD2, 'r', 'w', 0.5, 0)
                drawBoundaries(boundaries_lb, gridDomainD1, gridDomainD2, 'g', 'w', 0.7, 0)
                    subplot(subp11);hold on;
                drawBoundaries(boundaries_ub, gridDomainD1, gridDomainD2, 'r', 'w', 0.5, 0)
                drawBoundaries(boundaries_lb, gridDomainD1, gridDomainD2, 'g', 'w', 0.7, 0)
                    subplot(subp21);hold on;
                drawBoundaries(boundaries_ub, gridDomainD1, gridDomainD2, 'r', 'w', 0.5, 0)
                drawBoundaries(boundaries_lb, gridDomainD1, gridDomainD2, 'g', 'w', 0.7, 0)
                    subplot(subp22);hold on;
                drawBoundaries(boundaries_ub, gridDomainD1, gridDomainD2, 'r', 'w', 0.5, 0)
                drawBoundaries(boundaries_lb, gridDomainD1, gridDomainD2, 'g', 'w', 0.7, 0)
                    subplot(subp31);hold on;
                drawBoundaries(boundaries_ub, gridDomainD1, gridDomainD2, 'r', 'w', 0.5, 0)
                drawBoundaries(boundaries_lb, gridDomainD1, gridDomainD2, 'g', 'w', 0.7, 0)
                    subplot(subp32);hold on;
                drawBoundaries(boundaries_ub, gridDomainD1, gridDomainD2, 'r', 'w', 0.5, 0)
                drawBoundaries(boundaries_lb, gridDomainD1, gridDomainD2, 'g', 'w', 0.7, 0)
                    subplot(subp12);hold on;
                drawBoundaries(boundaries_ub, gridDomainD1, gridDomainD2, 'r', 'w', 0.5, 0)
                drawBoundaries(boundaries_lb, gridDomainD1, gridDomainD2, 'g', 'w', 0.7, 0)
            end
%             for k = 1 : numberOfBoundaries
%                 thisBoundary = boundaries{k};
%                 h = fill(gridDomainD1(thisBoundary(:,2)), gridDomainD2(thisBoundary(:,1)), 'w', 'LineWidth', 2, 'EdgeColor', 'r');
%                 hideLegend(h);
%                 set(h,'facealpha',.2)
%             end

        end
        
        
%------------------------------------------------------------------------------------------------
            subp34 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 3, 4)); cla;
%             fn_true_mean=meanFunc.getFnEval();
            hold on;
            hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),30, y_cons_history(1:nEvaluations,:),'filled') );
            hideLegend( scatter(x_history(nEvaluations,1), x_history(nEvaluations,2),30, 'LineWidth', 2, 'MarkerEdgeColor', 'r') );

            dom = objFunc.getXDomain;
            axis([dom(1,:), dom(2,:)]);
            colorbar();
            title(sprintf('Cons: fn:%s', meanFunc.getName()  ));
            xlabel('x');ylabel('y');


    end
end
% ======================================================================================
if(param_BO.isInfillSearchGrid)
    mysubplot(subplot_frame(1),subplot_frame(2),1);
    dim = sqrt(size(gridCoord,1));
    flippedImagesc(gridDomainD1, gridDomainD2, (log(objFunc.evalWithVecX(gridSquareD1, gridSquareD2)))');
    hold on;
    scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]);
    scatter(gridCoord(idx_min_mu,1), gridCoord(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[1 1 1]);
    scatter(gridCoord(maxIdx,1), gridCoord(maxIdx, 2), 20, 'k', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]);
    scatter(gridCoord(true_min_idx,1), gridCoord(true_min_idx,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[0 0 1]);
    colorbar();
    title(sprintf('fn:%s,n=%d', objFunc.getName(), sum(samplesize_history(n_init_sample+1:nEvaluations,1))));
    xlabel('x');ylabel('y');
    
    % ---------------------------------------------------------
    mysubplot(subplot_frame(1),subplot_frame(2),2);
    flippedImagesc(gridDomainD1, gridDomainD2, (reshape(mu, dim, dim))');
    hold on;
    scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled');
    scatter(gridCoord(idx_min_mu,1), gridCoord(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[1 1 1]);
    scatter(gridCoord(maxIdx,1), gridCoord(maxIdx, 2), 20, 'k', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]);
    colorbar();
    title('Mean');
    % ---------------------------------------------------------
    mysubplot(subplot_frame(1),subplot_frame(2),3);
    flippedImagesc(gridDomainD1, gridDomainD2, (reshape(acq, dim, dim))');
    hold on;
    scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled');
    scatter(gridCoord(idx_min_mu,1), gridCoord(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[1 1 1]);
    scatter(gridCoord(maxIdx,1), gridCoord(maxIdx, 2), 20, maxV, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]);
    colorbar();
    title(sprintf('%s',param_BO.getStrTypeAcquisition()));
end

% --------------------------------------------------------------------------------------

if(param_BO.info)
    finVinsual = toc(sttVisual);
    fprintf('[Visual] :%s\n', showPrettyElapsedTime(finVinsual));
end
return
