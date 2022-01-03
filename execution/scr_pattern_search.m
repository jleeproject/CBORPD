%% Pattern Search algorithm

%% SAMPLE SIZE
if(setting.typeProblem == TypeProblem.RobustDesignOptimization)
%% ----- Determine Sample Size ----------------------------------------------------

    switch samplesizeSelector.type
        case TypeSampleSize.Fixed
            samplesize = samplesizeSelector.fixedSamplesize;
        case TypeSampleSize.Adjust
            throwError('[SampleSizeSelector] Not defined type');
        otherwise
            throwError('[SampleSizeSelector] Not defined type');
    end

%     samplesize = samplesizeSelector.selectSamplesize(xnew, infillPredictor);
%     if(sum(cumul_samplesizes) + samplesize>=budget)
%         samplesize = budget - sum(cumul_samplesizes);
%         stop = true;
%     end
else
    samplesize = 1;
    setting.samplesize_for_variance_init = 1;
%     if(n_init_sample + nEvaluations>=budget)
    if(nEvaluations>=budget)
%         samplesize = sum(cumul_samplesizes) - budget;
        stop = true;
    end
end

maxEval = floor( (budget - n_init_sample*setting.samplesize_for_variance_init)/samplesize );
% maxEval = 10000;

    x_train = storage.getXHistoryTill(nEvaluations);
    y_train = storage.getYHistoryTill(nEvaluations);
    y_con_train = storage.getYConsHistoryTill(nEvaluations);
    
    feasibility = FeasibilitySamplingConstraint(param_BO);
    feasgap = @feasibility.feasGap;
    isfeas = feasgap(y_con_train)==0;
    if sum(isfeas)>0
        yy = y_train(isfeas);
        [minval,minidx] = min(yy);
        X0 = x_train(minidx, :);
    else
        [minval,minidx] = min(y_train);
        X0 = x_train(minidx, :);
    end
    
clear pscls;
pscls = PatternSearchData(simulator, samplesize, feasgap);
% objfcn = @(x) sum(x);
% confcn = @(x) 0;
LB = optProb.xDomain(:,1)';
UB = optProb.xDomain(:,2)';
% PSoptions = optimoptions(@patternsearch,'Display','iter', 'OutputFcn', @pscls.myoutput, 'MaxFunctionEvaluations', maxEval);
PSoptions = optimoptions(@patternsearch,'Display','off', 'OutputFcn', @pscls.myoutput, 'MaxFunctionEvaluations', maxEval);
[Xps,Fps] = patternsearch(@pscls.objfcn,X0,[],[],[],[],LB,UB,@pscls.confcn, PSoptions);
% [Xps,Fps] = patternsearch(objfcn,X0,[],[],[],[],LB,UB,confcn, PSoptions);


% disp(pscls.x_history)
% % disp(pscls.y_con_history)
% disp(pscls.feasgap_history)
% % disp(pscls.y_obj_history)
% % disp(pscls.x_iter_history)
% % disp(pscls.y_iter_history)
% disp(pscls.cnt_eval)
% disp(pscls.iter_history)
%     options = optimset();
%     [x fval] = fminsearch(@objfun, x0,options);
    

if maxEval<pscls.cnt_eval
    n_save = maxEval + 1;
else
    n_save = pscls.cnt_eval;
end

%% ------------------------------------------------------
for i=1:n_save
    nEvaluations = nEvaluations + 1;
    xnew = pscls.x_history(i,:);
    ynew = pscls.y_obj_history(i,:);
    y_orig = pscls.y_obs_history(i,:);
    y_cons = pscls.y_con_history(i,:);
    optMuX = xnew;
    optMuV = ynew;
    
    storage.saveXHistory(xnew, nEvaluations);
    storage.saveYHistory(ynew, nEvaluations);
    storage.saveObsHistory(y_orig, nEvaluations);
    if(param_BO.hasConstraints)
        storage.saveYConsHistory(y_cons, nEvaluations);
    end
    storage.saveSamplesizeHistory(samplesize, nEvaluations);
    % True value at xnew

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

%% ------------------------------------------------------
% storage.xHistory = pscls.x_history(n_init_sample+1:nEvaluations,:);
% storage.yHistory = pscls.y_obj_history(1:nEvaluations,:);
% storage.obsHistory = pscls.y_obs_history(1:nEvaluations,:);
% %         storage.saveXHistory(xnew, nEvaluations);
% %         storage.saveYHistory(ynew, nEvaluations);
% %         storage.saveObsHistory(y_orig, nEvaluations);
%         if(param_BO.hasConstraints)
% %             storage.saveYConsHistory(y_cons, nEvaluations);
%             storage.yConsHistory = pscls.y_con_history(1:nEvaluations,:);
%         end
% storage.samplesizeHistory = samplesize*ones(nEvaluations,1);
% %             storage.saveSamplesizeHistory(samplesize, nEvaluations);
% storage.xMuMinHistory = storage.xHistory;
% storage.yMuMinHistory = storage.yHistory;
% %             storage.saveXMuMinHistory(optMuX,nEvaluations);
% %             storage.saveYMuMinHistory(optMuV,nEvaluations);
%                 if( isa(optProb.objFunc,'AbsFunction'))
% storage.evalFAtXnewHistory = setting.optProb.objFunc.evalWithVecX(storage.xMuMinHistory);
% 
% %                     storage.saveEvalFAtXnewHistory( (setting.optProb.objFunc.evalWithVecX(optMuX)) , nEvaluations)
%                     feas = 1;
%                     feasgap = 0;
%                     for i = 1:optProb.nCon
%                         gval = optProb.conFunc{i}.evalWithVecX(optMuX);
%                         if(optProb.con.constraints{i}.hasLb && optProb.con.constraints{i}.hasUb )
%                             lb = optProb.con.constraints{i}.lb;
%                             ub = optProb.con.constraints{i}.ub;
%                             feas = feas.*( gval<=ub && lb<=gval );
%                             feasgap = max( max( [gval - ub , lb - gval, 0]), feasgap);
%                         elseif(optProb.con.constraints{i}.hasLb && ~optProb.con.constraints{i}.hasUb )
%                             lb = optProb.con.constraints{i}.lb;
%                             feas = feas.*( lb<=gval ) ;
%                             feasgap = max( max( [lb - gval , 0]), feasgap) ;
%                         elseif(~optProb.con.constraints{i}.hasLb && optProb.con.constraints{i}.hasUb )
%                             ub = optProb.con.constraints{i}.ub;
%                             feas = feas.*( gval<=ub ) ;
%                             feasgap = max( max( [gval - ub, 0]), feasgap) ;
%                         elseif(~optProb.con.constraints{i}.hasLb && ~optProb.con.constraints{i}.hasUb )
%                             error('Both ub and lb are not specified.');
%                         else
%                         end
%                         storage.saveFeasAtXnewHistory( feas , nEvaluations)
%                         storage.saveEvalGXnewHistory( gval , nEvaluations, i )
%                         storage.saveFeasgapAtXnewHistory( feasgap , nEvaluations, i )
%                     end
%                     storage.saveFeasAtXnewHistory( feas , nEvaluations)
%                     
%                     
%                 elseif( isa(optProb.objFunc,'GrapheneModelSolver'))
%                     data_samples_for_validation = simulator.evaluate(optMuX, setting.casestudy.validation.samplesize);
%                     samples_for_validation = data_samples_for_validation.samples;
%                     storage.saveEvalFAtXnewHistory( std(samples_for_validation,0,'all') , nEvaluations)
%                     
%                     gval = mean(samples_for_validation,'all');
%                     if(optProb.con.constraints{1}.hasLb && optProb.con.constraints{1}.hasUb )
%                         lb = optProb.con.constraints{1}.lb;
%                         ub = optProb.con.constraints{1}.ub;
%                         feas = gval<=ub && lb<=gval ;
%                         feasgap = max( [gval - ub , lb - gval, 0]);
%                     elseif(optProb.con.constraints{1}.hasLb && ~optProb.con.constraints{1}.hasUb )
%                         lb = optProb.con.constraints{1}.lb;
%                         feas = lb<=gval ;
%                         feasgap = max( [lb - gval , 0]) ;
%                     elseif(~optProb.con.constraints{1}.hasLb && optProb.con.constraints{1}.hasUb )
%                         ub = optProb.con.constraints{1}.ub;
%                         feas = gval<=ub ;
%                         feasgap = max( [gval - ub, 0]) ;
%                     elseif(~optProb.con.constraints{1}.hasLb && ~optProb.con.constraints{1}.hasUb )
%                         error('Both ub and lb are not specified.');
%                     else
%                     end
%                     storage.saveFeasAtXnewHistory( feas , nEvaluations)
%                     storage.saveEvalGXnewHistory( gval , nEvaluations, i )
%                     storage.saveFeasgapAtXnewHistory( feasgap , nEvaluations, i )
%                 else
%                     error('[AbsInfillOptimizer] Undefined Type.');
%                 end
% %         end
%     end
