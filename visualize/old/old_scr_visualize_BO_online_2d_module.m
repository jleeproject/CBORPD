sttVisual = tic ;
% if(dim-floor(dim)>0)
%     fprintf('ERROR; or Numerical ERROR?\n');
% end
gridCoord = visualizer.getGridCoord();
gridSquareD1=visualizer.getGridSquareD1();
gridSquareD2=visualizer.getGridSquareD2();
gridDomainD1=visualizer.getGridDomainD1();
gridDomainD2= visualizer.getGridDomainD2();
true_minX = visualizer.getTrueGlobalMinX();
true_fmin = visualizer.getTrueGlobalMinF();
% TO VISUALIZE
% queries_all : queryVals_acq, queryVals_uncert, queryVals_mu, queryVals_sigma

str_optimal= 'Optimal';
str_min_mu = 'min(mu)';
str_max_acq = 'min(acq)';
str_last = 'last';
subplot_frame = [3,4];
%% VISUALIZE
fig = figure(80);clf;
% fig.Position = [500 0 1500 1000]; % For 3x3
fig.Position = [100 0 2000 1000]; % For 3 x 4


if(param_BO.isInfillSearchGrid)
    subplot(subplot_frame(1),subplot_frame(2),1);
%     subplot(subplot_frame(1),subplot_frame(2),2);
%     subplot(subplot_frame(1),subplot_frame(2),3);
%     subplot(subplot_frame(1),subplot_frame(2),4);
%     subplot(subplot_frame(1),subplot_frame(2),5);
%     subplot(subplot_frame(1),subplot_frame(2),6);
%     subplot(subplot_frame(1),subplot_frame(2),7);
%     subplot(subplot_frame(1),subplot_frame(2),8);
%     subplot(subplot_frame(1),subplot_frame(2),9);
    dim = sqrt(size(gridCoord,1));
%     zz1sq=reshape(gridCoord(:,1),dim,dim);
%     zz2sq=reshape(gridCoord(:,2),dim,dim);
    flippedImagesc(gridDomainD1, gridDomainD2, (log(fn_eval_true(gridSquareD1, gridSquareD2)))');
    hold on;
    scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]);
    scatter(gridCoord(idx_min_mu,1), gridCoord(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[1 1 1]);
    scatter(gridCoord(maxIdx,1), gridCoord(maxIdx, 2), 20, 'k', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]);
    scatter(gridCoord(true_min_idx,1), gridCoord(true_min_idx,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[0 0 1]);
    colorbar();
    title(sprintf('fn:%s,n=%d', objFunc.getName(), sum(samplesize_history(n_init_sample+1:nEvaluations,1))));
    xlabel('x');ylabel('y');
    
    % ---------------------------------------------------------
    subplot(subplot_frame(1),subplot_frame(2),2);
    flippedImagesc(gridDomainD1, gridDomainD2, (reshape(mu, dim, dim))');
    hold on;
    scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled');
    scatter(gridCoord(idx_min_mu,1), gridCoord(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[1 1 1]);
    scatter(gridCoord(maxIdx,1), gridCoord(maxIdx, 2), 20, 'k', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]);
    colorbar();
    title('Mean');
    % ---------------------------------------------------------
    subplot(subplot_frame(1),subplot_frame(2),3);
    flippedImagesc(gridDomainD1, gridDomainD2, (reshape(acq, dim, dim))');
    hold on;
    scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled');
    scatter(gridCoord(idx_min_mu,1), gridCoord(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[1 1 1]);
    scatter(gridCoord(maxIdx,1), gridCoord(maxIdx, 2), 20, maxV, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]);
    colorbar();
    title(sprintf('%s',param_BO.getStrTypeAcquisition()));
    % ---------------------------------------------------------
%     subplot(subplot_frame(1),subplot_frame(2),4);
%     imagesc(gridDomainD1, gridDomainD2, (reshape(sigma, dim, dim))');
%     hold on;
%     scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled');
%     scatter(gridCoord(idx_min_mu,1), gridCoord(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[1 1 1]);
%     scatter(gridCoord(maxIdx,1), gridCoord(maxIdx, 2), 20, 'k', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]);
%     colorbar();
%     title('Sigma');
    % ---------------------------------------------------------


else
queryVals_mu = queries_hist.mu;
queryVals_acq = queries_hist.acq;
queries_all = queries_hist.all;
history_imprv = queries_hist.improv;
queryVals_sigma = queries_hist.sigma;


[~ , idx_min_mu] = min(queryVals_mu);
[~, idx_max_acq] = min(queryVals_acq);
    
x_history = storage.xHistory;
samplesize_history = storage.samplesizeHistory;
eval_f_at_xnew_history = storage.evalFAtXnewHistory;
xmin_history = storage.xminHistory;

% First Column : GRID
if(visualizer.isShowAllGrid())
    %                     subplot(subplot_frame(1),subplot_frame(2),[1 2]);
    idxFig = 1;
    subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idxFig, 4));
    % subplot(subplot_frame(1),subplot_frame(2),3);
    fn_eval_true=objFunc.getFnEval();
    flippedImagesc(gridDomainD1, gridDomainD2, (reshape(log(fn_eval_true(gridCoord(:,1), gridCoord(:,2))),visualizer.getDim(),visualizer.getDim()))');
%     hideLegend(    scatter(gridCoord(:,1), gridCoord(:,2),10, (log(fn_eval_true(gridCoord(:,1), gridCoord(:,2)))),'filled') ) ;
    hold on;
    hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
    showOneLegend( scatter(true_minX(1), true_minX(2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
    showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
    showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
    showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);
    % legend('Location','best' )
    colorbar();
    title(sprintf('fn:%s,n=%d', objFunc.getName(), sum(samplesize_history(n_init_sample+1:nEvaluations,1))));
    xlabel('x');ylabel('y');
    % --------------------------------------------------------------------------------------
    idxFig = 2;
    subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idxFig, 4));
    % subplot(subplot_frame(1),subplot_frame(2),3);
    fn_eval_true=objFunc.getFnEval();
    % imagesc(gridDomainD1, gridDomainD2, (log(fn_eval_true(gridSquareD1, gridSquareD2)))');
%     hideLegend(    scatter(gridCoord(:,1), gridCoord(:,2),10, (log(fn_eval_true(gridCoord(:,1), gridCoord(:,2)))),'filled') ) ;
    hold on;
    hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
    showOneLegend( scatter(true_minX(1), true_minX(2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
    showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
    showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
    showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);
    legend('Location','best' )
    colorbar();
    title(sprintf('True fn:%s,n=%d', objFunc.getName(), sum(samplesize_history(n_init_sample+1:nEvaluations,1))));
    xlabel('x');ylabel('y');
    % --------------------------------------------------------------------------------------

% ======================================================================================
% ======================================================================================
    % --------------------------------------------------------------------------------------
    idxFigGrid = 1;
    subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idxFigGrid, 1));
    
    [mu_est, sig_est] = predictor.predict(gridCoord);
%     sig_est = sqrt(var_est);
    flippedImagesc(gridDomainD1, gridDomainD2, (reshape(mu_est, visualizer.getDim(), visualizer.getDim()))');
    
%     hideLegend(    scatter(gridCoord(:,1), gridCoord(:,2),10, res,'filled') ) ;
    hold on;
    hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
    showOneLegend( scatter(true_minX(1), true_minX(2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
    showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
    showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
    showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);
    colorbar();
%     title(sprintf('fn:%s,n=%d', objFunc.getName(), sum(samplesize_history(n_init_sample+1:nEvaluations,1))));
    title('Mean');
    xlabel('x');ylabel('y');
    % --------------------------------------------------------------------------------------
%     idxFigGrid = 1;
    idx_col = 3;
    subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idxFigGrid, idx_col));
    acq_est_grid = acqFunc.acquire(predictor, gridCoord, iter);
    
    flippedImagesc(gridDomainD1, gridDomainD2, (reshape(acq_est_grid, visualizer.getDim(), visualizer.getDim()))');
    
%     hideLegend(    scatter(gridCoord(:,1), gridCoord(:,2),10, res,'filled') ) ;
    hold on;
    hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
    showOneLegend( scatter(true_minX(1), true_minX(2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
    showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
    showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
    showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);
    colorbar();
%     title(sprintf('fn:%s,n=%d', objFunc.getName(), sum(samplesize_history(n_init_sample+1:nEvaluations,1))));
    title('UCB');
    xlabel('x');ylabel('y');
    % --------------------------------------------------------------------------------------
    idx_col = 2;
    subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idxFigGrid, idx_col));
%     acq_est = acquisition(gridCoord);
    flippedImagesc(gridDomainD1, gridDomainD2, (reshape(sig_est, visualizer.getDim(), visualizer.getDim()))');
    
%     hideLegend(    scatter(gridCoord(:,1), gridCoord(:,2),10, res,'filled') ) ;
    hold on;
    hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
    showOneLegend( scatter(true_minX(1), true_minX(2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
    showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
    showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
    showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);
    colorbar();
%     title(sprintf('fn:%s,n=%d', objFunc.getName(), sum(samplesize_history(n_init_sample+1:nEvaluations,1))));
    title('Sigma');
    xlabel('x');ylabel('y');
% else
%     idxFig = 1;
%     subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idxFig, 4));
%     % subplot(subplot_frame(1),subplot_frame(2),3);
%     fn_eval_true=param_BO.getFnEval();
%     % imagesc(gridDomainD1, gridDomainD2, (log(fn_eval_true(gridSquareD1, gridSquareD2)))');
%     hideLegend(    scatter(gridCoord(:,1), gridCoord(:,2),10, (log(fn_eval_true(gridCoord(:,1), gridCoord(:,2)))),'filled') ) ;
%     hold on;
%     hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
%     showOneLegend( scatter(true_minX(1), true_minX(2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
%     showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
%     showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
%     showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);
%     % legend('Location','best' )
%     colorbar();
%     title(sprintf('fn:%s,n=%d', objFunc.getName(), sum(samplesize_history(n_init_sample+1:nEvaluations,1))));
%     xlabel('x');ylabel('y');
end
% ======================================================================================
% ======================================================================================
if(visualizer.isShowSearchTrails())
idxFigSearch = 2;
% --------------------------------------------------------------------------------------
% FIGURE 2-1
idx_col = 1;
subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idxFigSearch, idx_col));
% subplot(subplot_frame(1),subplot_frame(2),4);
% imagesc(gridDomainD1, gridDomainD2, (reshape(mu, dim, dim))');
% [X,Y] = meshgrid(queries_all(:,1), queries_all(:,2));
% contour(X,Y, queryVals_mu);

% contour(queries_all(:,1), queries_all(:,2), queryVals_mu);
hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
hold on;
scatter(queries_all(:,1), queries_all(:,2),10, queryVals_mu,'filled');
if(visualizer.isShowAllGrid())
    showOneLegend( scatter(true_minX(1), true_minX(2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
end
showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);
colorbar();
title('Mean');
%                     contour(zz1, zz2, reshape(mu, dim, dim));
% --------------------------------------------------------------------------------------
% FIGURE 2-2
idx_col = 3;
subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idxFigSearch, idx_col));

hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
hold on;
scatter(queries_all(:,1), queries_all(:,2),10, queryVals_acq,'filled');
if(visualizer.isShowAllGrid())
    showOneLegend( scatter(true_minX(1), true_minX(2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
end
showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);

plot(queries_all(history_imprv(:,2),1), queries_all(history_imprv(:,2),2),'r-');

colorbar();
title(sprintf('%s',param_BO.getStrTypeAcquisition()));
% --------------------------------------------------------------------------------------
% FIGURE 2-3
idx_col = 2;
% subplot(subplot_frame(1),subplot_frame(2),6);
% idxFigSearch = 2;
subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idxFigSearch, idx_col));

hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
hold on;
scatter(queries_all(:,1), queries_all(:,2),10, queryVals_sigma,'filled');
if(visualizer.isShowAllGrid())
    showOneLegend( scatter(true_minX(1), true_minX(2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
end
showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);
colorbar();
title('Sigma');
% --------------------------------------------------------------------------------------
% subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idxFigSearch, 4));
% scatter(queries_all(:,1), queries_all(:,2),10, acquisition(queries_all),'filled');
% hold on;
% hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
% showOneLegend( scatter(true_minX(1), true_minX(2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
% showOneLegend( scatter(queries_all(idx_min_mu,1), queries_all(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[0 1 0]) , str_min_mu);
% showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
% showOneLegend( scatter(queries_all(end,1), queries_all(end,2),20,'>', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5]) , str_last);
% colorbar();
% title('TRUE UCB');


end
% acquisition
end
% ======================================================================================
% ======================================================================================
% subplot(subplot_frame(1),subplot_frame(2),7);
idxFigStatus = 3;
subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idxFigStatus, 1));
cum_samplesize =  cumsum(samplesize_history(n_init_sample+1:nEvaluations,1));
% if(visualizer.isShowAllGrid())
%     diff_f_min = (  eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1)  - true_fmin );
%     idxNInf = find(eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1) == -inf);
%     minDiff = min(diff_f_min(find(eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1) > -inf)));
%     showOneLegend(    plot(  cum_samplesize, diff_f_min , 'r-', cum_samplesize(idxNInf),minDiff*ones(size(idxNInf)),'ro', [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), abs(  eval_f_at_xnew_history(nEvaluations,1)  - true_fmin ), 'ob') ,'Diff(y)');
%     title(sprintf('|f_{min} - min(mu)|: %.3g',abs(  eval_f_at_xnew_history(nEvaluations,1)  - true_fmin )))
% else
%     f_x = eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1);
%     showOneLegend(    plot(  cum_samplesize, f_x , 'r-', [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)),  eval_f_at_xnew_history(nEvaluations,1), 'ob') ,'Diff(y)');
%     title(sprintf('f_(x): %.3g',(  eval_f_at_xnew_history(nEvaluations,1)  )))
% end
    diff_f_min = (  eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1)  - log(objFunc.getOptVal() ));
%     idxNInf = find(eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1) == -inf);
%     minDiff = min(diff_f_min(find(eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1) > -inf)));
%     showOneLegend(    plot(  cum_samplesize, diff_f_min , 'r-', cum_samplesize(idxNInf),minDiff*ones(size(idxNInf)),'ro', [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), diff_f_min(end), 'ob') ,'Diff(y)');
    showOneLegend(    plot(  cum_samplesize, diff_f_min , 'r-',  [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), diff_f_min(end), 'ob') ,'Diff(y)');
    title(sprintf('|f_{min} - min(mu)|: %.3g',diff_f_min(end)))
% opt_sol_fn
% opt_val_fn
% idxFigStatus = 3;
% if(visualizer.isShowAllGrid())
%     subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idxFigStatus, 2));
%     diff_x_f_min = vecnorm(xmin_history(n_init_sample+1:nEvaluations,:) - true_minX, 2, 2);
%     plot( cum_samplesize, diff_x_f_min, 'r-' , [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), norm(xmin_history(nEvaluations,:) - true_minX, 2), 'ob');
%     title(sprintf('|x_{min} - argmin(mu)|: %.3g',sum(abs(xmin_history(nEvaluations,:) - true_minX))))
% end
    subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idxFigStatus, 2));
%     diff_x_f_min = vecnorm(xmin_history(n_init_sample+1:nEvaluations,:) - opt_sol_fn, 2, 2);
%     diff_x_f_min = min(vecnorm(xmin_history(n_init_sample+1:nEvaluations,:) - opt_sol_fn, 2, 2));
    diff_x_f_min = fn_get_minimum_distance(xmin_history(n_init_sample+1:nEvaluations,:), objFunc.getOptSol() );
    plot( cum_samplesize, diff_x_f_min, 'r-' , [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), diff_x_f_min(end), 'ob');
    title(sprintf('|x_{min} - argmin(mu)|: %.3g',diff_x_f_min(end)))
% subplot(subplot_frame(1),subplot_frame(2),9);
idxFigStatus = 3;
subplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, idxFigStatus, 3));
plot( cum_samplesize, samplesize_history(n_init_sample+1:nEvaluations,1), 'r-' , [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  sum(samplesize_history(n_init_sample+1:nEvaluations,1)), samplesize_history(nEvaluations,1), 'ob');
title(sprintf('Sample Size: %d (iter=%d)',samplesize_history(nEvaluations,1), iter))

% --------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------
% --------------------------------------------------------------------------------------

% figure(89);clf;
% 
% subplot(2,2,1);
% imagesc(gridDomainD1, gridDomainD2, (log(fn_eval_true(gridSquareD1, gridSquareD2)))');
% hold on;
% scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),20, y_history(1:nEvaluations,1),'filled','MarkerEdgeColor',[0 0 0]);
% % idx_min_mu : index of mean function (mu) at the last iteration (Grid)
% % maxIdx : index of maxima of acquisition function (Grid)
% % scatter(gridCoord(maxIdx,1), gridCoord(maxIdx, 2), 20, 'k', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]);
% % scatter(gridCoord(true_min_idx,1), gridCoord(true_min_idx,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[0 0 1]);
% colorbar();
% 
% 
% if(param_BO.isTypeAcquisitionUCB)
% 
%     subplot(2,2,2);
%     imagesc(gridDomainD1,gridDomainD2, reshape( (f_min-mu).*Gaussian_CDF((f_min-mu)./sigma)+sigma.*Gaussian_PDF((f_min-mu)./sigma),dim,dim)')
%     title('EI');
%     colorbar();
% 
%     subplot(2,2,4);
%     zz = norminv(.99,0,1);
%     [v_ucb,idx_ucb]=min( mu-zz*sigma );
% 
%     idgridCoord = (mu(idx_ucb) + zz*sigma(idx_ucb) < mu);
%     imagesc(gridDomainD1, gridDomainD2, reshape(idgridCoord,dim,dim)');
%     colorbar();
% else
%     subplot(2,2,2);
%     imagesc(gridDomainD1,gridDomainD2, reshape( -mu + sqrt(0.2 * fn_dim_eval * log(2*dim_eval*iter)) * sigma ,dim,dim)')
%     title('UCB');
%     colorbar();
% end
% subplot(2,2,3);
% plot(cumsum(samplesize_history(n_init_sample+1:nEvaluations,1)), samplesize_history(n_init_sample+1:nEvaluations,1), '-' , [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--');
% title(sprintf('sample size=%d',samplesize_history(nEvaluations,1)));
% colorbar();

finVinsual = toc(sttVisual);
fprintf('[Visual] :%s\n', showPrettyElapsedTime(finVinsual));

return
%                     %% EXAMPLE
%                     gridCoord2 = [gridCoord' fliplr(gridCoord')]';
%                     ys2 = [(mu+2*sigma)' fliplr((mu-2*sigma)')];
%                     patch=fill(gridCoord2,ys2,options.plot.color);
%                     set(patch, 'FaceAlpha', options.plot.alpha);
%                     hold on;
%                     plot(gridCoord, log(fn_eval_true(gridCoord)), 'k--', gridCoord, mu, 'r-', x, y, 'ob', xnew, ynew, 'or', 'LineWidth', 2);
%                     dim = sqrt(size(gridCoord,1));
%                     gridSquareD1=reshape(gridCoord(:,1),dim,dim);
%                     gridSquareD2=reshape(gridCoord(:,2),dim,dim);
%                     contour(zz1, zz2, fn_eval_true(zz1, zz2));
%                     History
%                     scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),20, y_history(1:nEvaluations,1));
%                     scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled');
% %                     Est MINIMUM POINT
% %                     scatter(gridCoord(idx_min_mu,1), gridCoord(idx_min_mu,2),20, fmin_history(nEvaluations+1,1), 'filled','d');
%                     scatter(gridCoord(idx_min_mu,1), gridCoord(idx_min_mu,2),20, 'LineWidth',1.5,'MarkerEdgeColor',[0 0 0]);
%                     Next Point Evaluated 
%                     scatter(gridCoord(maxIdx,1), gridCoord(maxIdx, 2), 20, maxV,'r', 'filled', 'LineWidth',1.5,'MarkerEdgeColor',[0 .5 .5],...
%                         'MarkerFaceColor',[1 0 0]);
%                     TRUE MINIMUM POINT
%                     scatter(gridCoord(idx_min_mu,1), gridCoord(idx_min_mu,2),20, fmin_history(nEvaluations+1,1), 'filled','d');
% true_fmin,true_min_idx


    hideLegend(    scatter(gridCoord(:,1), gridCoord(:,2), 7, (log(fn_eval_true(gridCoord(:,1), gridCoord(:,2)))), '+') ) ;

%     scatter(gridCoord(:,1), gridCoord(:,2), 7, (log(fn_eval_true(gridCoord(:,1), gridCoord(:,2)))), '+')
fig.Position = [100 0 700 700]; % For 3 x 4

clf;
[ ~, idx_min] = min(acquisition(gridCoord));
scatter(gridCoord(:,1), gridCoord(:,2), 7, ((acquisition(gridCoord))),'filled');
hold on;
plot(gridCoord(idx_min,1), gridCoord(idx_min,2),'ro');
colorbar()
% scatter(queries_all(:,1), queries_all(:,2),10, queryVals_acq,'filled');

idx_show = 20;
idxs = [1:idx_show];
clf;
scatter(queries_all(idxs,1), queries_all(idxs,2),10, queryVals_acq(idxs),'filled');
hold on;
showOneLegend( scatter(queries_all(idx_max_acq,1), queries_all(idx_max_acq, 2), 20, 'v', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]) , str_max_acq);
plot(queries_all(history_imprv(:,2),1), queries_all(history_imprv(:,2),2),'r-');
colorbar()
