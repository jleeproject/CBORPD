% optProb.objFunc = optProb.objFunc;
% plotIdx = [ 3 4];
show_xdiff = false;
sttVisual = tic ;
% gridCoord = visualizer.getGridCoord();
x_history = storage.xHistory;
current_x = x_history(nEvaluations,:);
if(numel(current_x)>=4);plotIdx = [ 3 4]; else; plotIdx = [ 1 2];  end;
% gridCoord = visualizer.getGridCoordAtX(current_x);
gridCoord = visualizer.getGridCoord(current_x, plotIdx);

gridSquareD1=visualizer.getGridSquareD1();
gridSquareD2=visualizer.getGridSquareD2();
gridDomainD1=visualizer.getGridDomainD1();
gridDomainD2= visualizer.getGridDomainD2();
% true_minX = visualizer.getTrueGlobalMinX();
true_fmin = visualizer.getTrueGlobalMinF();
true_minX = optProb.objFunc.getOptSol();
% TO VISUALIZE
% mysubplot = @(x1,x2,x3) mysubplot(x1,x2,x3);
xargs = mat2arg(gridCoord);

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

    if(visualizer.isShowAllGrid())
    % [First Column] : Objective Function (Grid)
    %% ======================================================================================
        % ----- True Objective Function (Grid) ---------------------------------------------------------------------------------
        subp13 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 1, 3));
        fn_true_mean=optProb.objFunc.getFnEval();
%         xargs=mat2arg(gridCoord);
        
        flippedImagesc(gridDomainD1, gridDomainD2, (reshape(log(fn_true_mean(xargs{:})),visualizer.getDim(),visualizer.getDim()))');
        hold on; 
        showAllOtherPoints(x_history, nEvaluations, x_mu_history, true_minX, str_optimal, color_darkgrey);
    
        title(sprintf('fn:%s,n=%d', optProb.objFunc.getName(), sum(samplesize_history(n_init_sample+1:nEvaluations,1))));
        % --------------------------------------------------------------------------------------

        subp23 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 2, 3));
        [~, pof] = funcFeas(gridCoord, 1);

%         fn_true_mean=optProb.objFunc.getFnEval();
%         xargs=mat2arg(gridCoord);
        pof_in_shape = (reshape((pof),visualizer.getDim(),visualizer.getDim()))';
        flippedImagesc(gridDomainD1, gridDomainD2, pof_in_shape);
        hold on; 
        boundaries_pof =  bwboundaries(pof_in_shape<0.5);
            drawBoundaries(boundaries_pof,gridDomainD1, gridDomainD2);
        showAllOtherPoints(x_history, nEvaluations, x_mu_history, true_minX, str_optimal, color_darkgrey);
    
%         title(sprintf('exp(f):%s,n=%d', optProb.objFunc.getName(), sum(samplesize_history(n_init_sample+1:nEvaluations,1))));
        title('PoF');
        % --------------------------------------------------------------------------------------

    % ======================================================================================
        % ----- Estimated Objective Function (Grid) ---------------------------------------------------------------------------------
        % --------------------------------------------------------------------------------------
        subp11 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 1, 1));

        [mu_est, sig_est] = exactPredictor.predict(gridCoord);
        flippedImagesc(gridDomainD1, gridDomainD2, (reshape(mu_est, visualizer.getDim(), visualizer.getDim()))');
        
        showAllOtherPoints(x_history, nEvaluations, x_mu_history, true_minX, str_optimal, color_darkgrey);
        title('Obj:Mean(log(Sigma)):\mu_{log\sigma}');
        % --------------------------------------------------------------------------------------
        subp21 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 2, 1));
        flippedImagesc(gridDomainD1, gridDomainD2, (reshape(sig_est, visualizer.getDim(), visualizer.getDim()))');

        showAllOtherPoints(x_history, nEvaluations, x_mu_history, true_minX, str_optimal, color_darkgrey);
        title('Obj:Sd(Sigma):s(\sigma)');
        % --------------------------------------------------------------------------------------
        subp31 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 3, 1));
        acq_est_grid = acquisition(gridCoord);

        flippedImagesc(gridDomainD1, gridDomainD2, (reshape(acq_est_grid, visualizer.getDim(), visualizer.getDim()))');

        showAllOtherPoints(x_history, nEvaluations, x_mu_history, true_minX, str_optimal, color_darkgrey);
        title(sprintf('Acq Func:%s',infillOptimizer.acqFuncName));
    end
    % ======================================================================================
        % ----- Evaluated Points from Objective Function ---------------------------------------------------------------------------------
        subp33 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 3, 3));cla
%         hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),30, y_history(1:nEvaluations,:),'filled') );
    
            [pred_con_mu, sigma_con_mu] = multiPredictors.cell_predictors{1}.predict(gridCoord);
            ucb_g = pred_con_mu + 1.5.* sigma_con_mu;
            ucb_g_inshape = reshape( ucb_g ,visualizer.getDim,visualizer.getDim)';
            boundaries_ucb_g_0 = bwboundaries(ucb_g_inshape<0);
            
            flippedImagesc(gridDomainD1, gridDomainD2, ucb_g_inshape);
            hold on; 
            drawBoundaries(boundaries_ucb_g_0,gridDomainD1, gridDomainD2);
            dom = optProb.objFunc.getXDomain;
            axis([dom(1,:), dom(2,:)]);
            colorbar();
%             title(sprintf('Cons: fn:%s', optProb.conFunc{1}.getName()  ));
            xlabel('x');ylabel('y');
        
        showEstimatedOptimals(x_history, nEvaluations, x_mu_history, true_minX, str_optimal, color_darkgrey);
        legend('Location','best' )
        title(sprintf('Objective Function: Evaluated Points' ));
        % --------------------------------------------------------------------------------------
    % ======================================================================================
    % [Second Column] : Objective Function (Search Trails)
    if(visualizer.isShowSearchTrails())
        subp12 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 1, 2));
        scatter(queries_all(:,1), queries_all(:,2),10, queryVals_mu,'filled');

        showAllOtherPoints(x_history, nEvaluations, x_mu_history, true_minX, str_optimal, color_darkgrey);
        title('Obj:Mean(log(Sigma)):\mu_{log\sigma}');
        % --------------------------------------------------------------------------------------
        subp22 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 2, 2));
        scatter(queries_all(:,1), queries_all(:,2),10, queryVals_sigma,'filled');

        showAllOtherPoints(x_history, nEvaluations, x_mu_history, true_minX, str_optimal, color_darkgrey);
        title('Obj:Sd(Sigma):s(\sigma)');
        % --------------------------------------------------------------------------------------
        subp32 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 3, 2));
        scatter(queries_all(:,1), queries_all(:,2),10, queryVals_acq,'filled');

        showAllOtherPoints(x_history, nEvaluations, x_mu_history, true_minX, str_optimal, color_darkgrey);
        title(sprintf('%s',infillOptimizer.acqFuncName));
        % --------------------------------------------------------------------------------------

    end
end
% ======================================================================================

%% Performance Summary
subp15 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 1, 5));
cla;
hold on;
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
    diff_f_min = (  eval_f_at_xnew_history(n_init_sample+1:nEvaluations,1)  - log(optProb.objFunc.getOptVal() ));
    showOneLegend(    plot(  cum_samplesize, diff_f_min , 'r-',  [0 sum(samplesize_history(n_init_sample+1:nEvaluations,1))], [0 0],'k--',  cum_samplesize(end), diff_f_min(end), 'ob') ,'Diff(y)');
    title(sprintf('|f_{min} - min(mu)|: %.3g',diff_f_min(end)))
    val = diff_f_min;
    subp15.YLim(1) = 0;
end
% show feasibility
hold on;



%% -----------------------------------------------------
if(param_BO.hasConstraints)
    [feasible, strCon] =  isFeasible(optProb.con.constraints{1}, optProb.conFunc{1}.fnEval , x_mu_history(n_init_sample+1:nEvaluations,:));
    idxFeas = find(feasible);
    idxInfeas = find(1-feasible);
    hideLegend( plot(  cum_samplesize(idxInfeas), val(idxInfeas), 'xr', 'MarkerSize',12, 'LineWidth',3 ));
    hideLegend( plot(  cum_samplesize(idxFeas), val(idxFeas), 'og', 'MarkerSize',12, 'LineWidth',3) );

    legend('Location','best');;

    subp15.XLim(1) = xLimLeft;
    subp15.XLim(2) = xLimRight;
    
    

    subp25 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 2, 5));cla;
    hold on;
    feasgap = feasibilityGap(optProb.con.constraints{1}, optProb.conFunc{1}.fnEval, x_mu_history(n_init_sample+1:nEvaluations,:));
    showOneLegend( plot(  cum_samplesize, feasgap, '-g') , 'Feas Gap' );    
    [feasible, strCon] =  isFeasible(optProb.con.constraints{1}, optProb.conFunc{1}.fnEval , x_mu_history(n_init_sample+1:nEvaluations,:));
    idxInfeas = find(1-feasible);
    hideLegend( plot(  cum_samplesize(idxInfeas), feasgap(idxInfeas), 'xr', 'MarkerSize',10, 'LineWidth',3) );    
    hideLegend( plot(  cum_samplesize(idxFeas), feasgap(idxFeas), 'og', 'MarkerSize',12, 'LineWidth',3) );


    [~, pCons] = funcFeas(x_mu_history(n_init_sample+1:nEvaluations,:), 0);
%     feasgap = feasibilityGap(cons{1}, optProb.conFunc{1}.fnEval, x_mu_history(n_init_sample+1:nEvaluations,:));
    showOneLegend( plot(  cum_samplesize, pCons, '-m'), 'P(C(x))' );    
    hideLegend(plot([cum_samplesize(1),cum_samplesize(end)], [1,1],'k--'));
    hideLegend(plot([cum_samplesize(1),cum_samplesize(end)], [0,0],'k--'));
    xLimLeft = sum(samplesize_history(1:n_init_sample+1,1));
    xLimRight = sum(samplesize_history)+1;
    subp25.XLim = [xLimLeft, xLimRight];
    subp25.YLim = [-.05, 1.05];
    
    xx = cum_samplesize(idxInfeas);
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
%% Regarding Constraints
if(param_BO.hasConstraints)
    cons = optProb.con.constraints;
    % --------------------------------------------------------------------------------------
    [feas_rows, strCon] =  isFeasible(cons{1}, optProb.conFunc{1}.fnEval, gridCoord);
    feas_eval = reshape(feas_rows,visualizer.getDim,visualizer.getDim)';
    boundaries_true = bwboundaries(feas_eval);
    if(cons{1}.typeConstraint  == TypeConstraint.UnknownMeanConstraint)

        if(visualizer.isShowAllGrid())
        % ----- Evaluated points ---------------------------------------------------------------------------------
            subp14 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 1, 4)); cla;
            fn_true_mean=optProb.conFunc{1}.getFnEval();
            
            vals = fn_true_mean(xargs{:});
            flippedImagesc(gridDomainD1, gridDomainD2, reshape( vals,visualizer.getDim,visualizer.getDim)');
            hold on;

            showAllOtherPoints(x_history, nEvaluations, x_mu_history, true_minX, str_optimal, color_darkgrey);
            title(sprintf('True Cons: fn:%s,%s', optProb.conFunc{1}.getName(), strCon  ));


           
            subplot(subp14); hold on;
            drawBoundaries(boundaries_true, gridDomainD1, gridDomainD2, 'w', 'w', 0.5, 0.2)
            subplot(subp13); hold on;
            drawBoundaries(boundaries_true, gridDomainD1, gridDomainD2, 'w', 'w', 0.5, 0.2)


        end
        if(visualizer.isShowSearchTrails())
%------------------------------------------------------------------------------------------------
            subp24 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 2, 4)); cla;
            hold on;
            fn_true_mean=optProb.conFunc{1}.getFnEval();
            
%             [pred_con_mu, sigma_con_mu] = multiPredictors.cell_predictors{1}.predict(gridCoord);
            % COMMENTED: Estimated GP shape
            flippedImagesc(gridDomainD1, gridDomainD2, reshape(pred_con_mu,visualizer.getDim,visualizer.getDim)');
            
            hold on;
            
            showAllOtherPoints(x_history, nEvaluations, x_mu_history, true_minX, str_optimal, color_darkgrey);
            dom = optProb.objFunc.getXDomain; axis([dom(1,:), dom(2,:)]);
            title(sprintf('Est Cons: fn:%s,%s (u=%.3g)', optProb.conFunc{1}.getName(), strCon, urnd_feasiblity  ));

%------------------------------------------------------------------------------------------------

            feas_eval = funcFeas(gridCoord, .05);
            feasible = reshape(feas_eval,visualizer.getDim,visualizer.getDim)' ==-1;
            boundaries_ub = bwboundaries(feasible);
            
            feas_eval = funcFeas(gridCoord, .95);
            feasible = reshape(feas_eval,visualizer.getDim,visualizer.getDim)' ==-1;
            boundaries_lb = bwboundaries(feasible);
%------------------------------------------------------------------------------------------------
%             xlabel('x');ylabel('y');

            subplot(subp24); hold on;
            drawBoundaries(boundaries_true, gridDomainD1, gridDomainD2, 'w', 'w', 0.3, 0)
%             drawBoundaries(boundaries_ub, gridDomainD1, gridDomainD2, 'm', 'w', 0.3, 0)
%             drawBoundaries(boundaries_lb, gridDomainD1, gridDomainD2, 'm', 'w', 0.3, 0)
                drawBoundaries(boundaries_ub, gridDomainD1, gridDomainD2, 'm', 'w', 0.3, 0)
                drawBoundaries(boundaries_lb, gridDomainD1, gridDomainD2, 'g', 'w', 0.7, 0)

            if( type_infill_opt.strcmpi('DirectConstSampling'))
                [feas_eval, pCons] = funcFeas(gridCoord, urnd_feasiblity);
                feasible = reshape(feas_eval,visualizer.getDim,visualizer.getDim)' ==-1;
                boundaries_est = bwboundaries(feasible);
                
                drawEstFeasibleBoundaries(subp24, boundaries_est, gridDomainD1, gridDomainD2);
%                 drawEstFeasibleBoundaries(subp33, boundaries_est, gridDomainD1, gridDomainD2);
                drawEstFeasibleBoundaries(subp22, boundaries_est, gridDomainD1, gridDomainD2);
                drawEstFeasibleBoundaries(subp32, boundaries_est, gridDomainD1, gridDomainD2);
                drawEstFeasibleBoundaries(subp12, boundaries_est, gridDomainD1, gridDomainD2);
                
                if(visualizer.isShowAllGrid())
                    drawEstFeasibleBoundaries(subp11, boundaries_est, gridDomainD1, gridDomainD2);
                    drawEstFeasibleBoundaries(subp21, boundaries_est, gridDomainD1, gridDomainD2);
                    drawEstFeasibleBoundaries(subp31, boundaries_est, gridDomainD1, gridDomainD2);
                end
            else
%                 drawLUFeasibleBoundaries(subp33, boundaries_ub, boundaries_lb, gridDomainD1, gridDomainD2);
                drawLUFeasibleBoundaries(subp22, boundaries_ub, boundaries_lb, gridDomainD1, gridDomainD2);
                drawLUFeasibleBoundaries(subp32, boundaries_ub, boundaries_lb, gridDomainD1, gridDomainD2);
                drawLUFeasibleBoundaries(subp12, boundaries_ub, boundaries_lb, gridDomainD1, gridDomainD2);
                if(visualizer.isShowAllGrid())
                    drawLUFeasibleBoundaries(subp11, boundaries_ub, boundaries_lb, gridDomainD1, gridDomainD2);
                    drawLUFeasibleBoundaries(subp21, boundaries_ub, boundaries_lb, gridDomainD1, gridDomainD2);
                    drawLUFeasibleBoundaries(subp31, boundaries_ub, boundaries_lb, gridDomainD1, gridDomainD2);
                end
            end

        end
        
        
%------------------------------------------------------------------------------------------------
            subp34 = mysubplot(subplot_frame(1),subplot_frame(2),fnGetIdxOfSubplotWithRowCol(subplot_frame, 3, 4)); cla;
            hold on;
%             hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),30, y_cons_history(1:nEvaluations,:),'filled') );
%             hideLegend( scatter(x_history(nEvaluations,1), x_history(nEvaluations,2),30, 'LineWidth', 2, 'MarkerEdgeColor', 'r') );
            
            flippedImagesc(gridDomainD1, gridDomainD2, reshape(sigma_con_mu,visualizer.getDim,visualizer.getDim)');
            dom = optProb.objFunc.getXDomain;
            axis([dom(1,:), dom(2,:)]);
            colorbar();
            title(sprintf('Cons: fn:%s', optProb.conFunc{1}.getName()  ));
            xlabel('x');ylabel('y');


    end
end
% ======================================================================================
if(param_BO.isInfillSearchGrid)
    mysubplot(subplot_frame(1),subplot_frame(2),1);
    dim = sqrt(size(gridCoord,1));
    flippedImagesc(gridDomainD1, gridDomainD2, (log(optProb.objFunc.evalWithVecX(gridSquareD1, gridSquareD2)))');
    hold on;
    scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]);
    scatter(gridCoord(idx_min_mu,1), gridCoord(idx_min_mu,2),20,'d', 'LineWidth',1.5,'MarkerEdgeColor',[1 1 1]);
    scatter(gridCoord(maxIdx,1), gridCoord(maxIdx, 2), 20, 'k', 'LineWidth',1.5,'MarkerEdgeColor',[1 0 0]);
    scatter(gridCoord(true_min_idx,1), gridCoord(true_min_idx,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[0 0 1]);
    colorbar();
    title(sprintf('fn:%s,n=%d', optProb.objFunc.getName(), sum(samplesize_history(n_init_sample+1:nEvaluations,1))));
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
annotateTextTitle(fig, infillOptimizer.strStep);

if(param_BO.info)
    finVinsual = toc(sttVisual);
    fprintf('[Visual] :%s\n', showPrettyElapsedTime(finVinsual));
end

function showAllOtherPoints(x_history, nEvaluations, x_mu_history, true_minX, str_optimal, color_darkgrey);
    hold on; 
    showEvaluatedPoints(x_history, nEvaluations);
    showEstimatedOptimals(x_history, nEvaluations, x_mu_history, true_minX, str_optimal, color_darkgrey);
    colorbar();
    xlabel('x');ylabel('y');
end


function showEvaluatedPoints(x_history, nEvaluations);
    hold on;
    hideLegend( scatter(x_history(1:nEvaluations,1), x_history(1:nEvaluations,2),10, 'k','filled','MarkerEdgeColor',[0 0 0]) );
end
function showEstimatedOptimals(x_history, nEvaluations, x_mu_history, true_minX, str_optimal, color_darkgrey)
    hold on;
    showOneLegend( scatter(true_minX(:,1), true_minX(:,2),50, 'LineWidth',1.5,'MarkerEdgeColor',[1 0 1]), str_optimal );
    showOneLegend( scatter(x_mu_history(nEvaluations,1), x_mu_history(nEvaluations,2),150, 'x', 'LineWidth', 2, 'MarkerEdgeColor', color_darkgrey), 'Est. Opt' ); % Estimated Optimal
    showOneLegend( scatter(x_history(nEvaluations,1), x_history(nEvaluations,2),30, 'LineWidth', 2, 'MarkerEdgeColor', 'r','MarkerEdgeAlpha',.7), 'Eval X' ); % SELECTED POINT
end

function drawEstFeasibleBoundaries(subp, boundaries_est, gridDomainD1, gridDomainD2)
    subplot(subp); hold on;
    drawBoundaries(boundaries_est, gridDomainD1, gridDomainD2);
end


function drawLUFeasibleBoundaries(subp, boundaries_ub, boundaries_lb, gridDomainD1, gridDomainD2)
    subplot(subp); hold on;
    drawBoundaries(boundaries_ub, gridDomainD1, gridDomainD2, 'r', 'w', 0.5, 0)
    drawBoundaries(boundaries_lb, gridDomainD1, gridDomainD2, 'g', 'w', 0.7, 0)
end
