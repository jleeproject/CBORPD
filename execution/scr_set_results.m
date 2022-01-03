    %         if( param_BO.isInitSmallSamplesize())
    %             fprintf('iter: %d/%d\n', iter, floor((budget - n_init_sample*3)/samplesize)+1);
    %         else
    %             fprintf('iter: %d/%d\n', iter, floor((budget - n_init_sample*samplesize)/samplesize)+1);
    %         end
            fprintf('[[ FINISHED ]] # Samples: %d/%d\n', sum(storage.samplesizeHistory), budget );
%             fprintf('[[ FINISHED ]]\n');
            
%             struct_result.x_history  = x_history ;
%             struct_result.y_history  = y_history ;
%             struct_result.eval_f_at_xnew_history  = eval_f_at_xnew_history ;
%             struct_result.samplesize_history  = samplesize_history ;
%             struct_result.fmin_history  = fmin_history ;
%             struct_result.xmin_history  = xmin_history ;
%             struct_result.hist_beta_t  = hist_beta_t ;
            struct_result.storage = storage;
            struct_result.setting = setting;
            struct_result.x_history  = storage.getXHistoryTill(nEvaluations) ;
            struct_result.y_history  = storage.getYHistoryTill(nEvaluations);
            struct_result.eval_f_at_xnew_history  = storage.evalFAtXnewHistory;

            struct_result.eval_g_at_xnew_history  = storage.evalGXnewHistory;

            struct_result.feas_history  = storage.feasAtXnewHistory;
            struct_result.feasgap_history  = storage.feasgapAtXnewHistory;
            
            
            struct_result.samplesize_history  = storage.samplesizeHistory;
%             struct_result.fmin_history  = storage.
            struct_result.xMuMinHistory  = storage.xMuMinHistory;
            struct_result.yMuMinHistory  = storage.yMuMinHistory;
            struct_result.hist_beta_t  = storage.histBetaT;
% -------------------------------------------------------
            struct_result.typeSamplesize  = storage.typeSamplesize;
            struct_result.givenSamplesize = storage.givenSamplesize;
% -------------------------------------------------------

            struct_result.optProb = optProb;
            struct_result.simul = simul;
        
            struct_result.nEvaluations  = nEvaluations ;        
            struct_result.n_init_sample = n_init_sample;
            struct_result.param_BO = param_BO;
            struct_result.objFunc = optProb.objFunc;
            
            struct_result.type_cell_samplesize = type_cell_samplesize;
            struct_result.cell_obj_func = cell_obj_func;
            struct_result.type_cell_settings = type_cell_settings;
            struct_result.nRepeats = nRepeats;
% -------------------------------------------------------
            if(show_plot_online)
                struct_result.visualizer = visualizer;
            else
                struct_result.visualizer = [];
            end
            struct_result.budget = budget;
%             struct_result.samplesize = samplesize;
            struct_result.elapsedTime = toc(stt_experiment);
            
            struct_result.elapsedTime2 = storage.elapsedTotalEval; % For validation
            struct_result.elapsedTimeExp = storage.elapsedTotalExp;
            struct_result.elapsedTimeBoWoEval = storage.elapsedTotalBoWoEval;
%             struct_result.final_opt_x = final_opt_x;
%             struct_result.final_opt_f_mu = final_opt_f_mu;
%             struct_result.final_opt_f_true = final_opt_f_true;
            struct_result.cumElapsedTime = storage.cumElapsedTime;
            struct_result.iterElapsedTime = storage.iterElapsedTime;
%             struct_result.adj_samplesize = adj_samplesize;


            if( isa(optProb.objFunc,'AbsFunction'))
                struct_result.opt_val_fn = optProb.objFunc.getOptVal();
                struct_result.opt_sol_fn = optProb.objFunc.getOptSol();
                struct_result.obj = simul.obj;
                if(isfield(simul,'con'))
                    struct_result.con = simul.con;
                end
%                 this.xDomain = optProb.objFunc.getXDomain();
            elseif( isa(optProb.objFunc,'GrapheneModelSolver'))
                struct_result.opt_val_fn = [];
                struct_result.opt_sol_fn = [];
%                 this.xDomain = optProb.objFunc.decisionVariables.domains.all;
                struct_result.obj = optProb.objFunc;
%                 if(isfield(optProb,'con'))
%                     struct_result.con = optProb.con;
%                 end
            else
                error('[AbsInfillOptimizer] Undefined Type.');
            end


            if(exist('power','var'))
                struct_result.power = power;
            end
            
            if(exist('funcFeas','var'))
                struct_result.funcFeas = funcFeas;
            end

            struct_result.hasLiberal = storage.hasLiberal;
            struct_result.eval_f_at_xnew_liberal_history = storage.evalFAtXnewLiberalHistory;
            struct_result.xMuMinLiberalHistory = storage.xMuMinLiberalHistory;
            struct_result.yMuMinLiberalHistory = storage.yMuMinLiberalHistory;

            struct_result.cell_mean_func = cell_mean_func;
            struct_result.nTypeSettings = nTypeSettings;
            struct_result.xDomain = xDomain;
            struct_result.simulator = simulator;
            struct_result.initializer = initializer;
            struct_result.acqFunc = acqFunc;
            struct_result.infillOptimizer = infillOptimizer;
%             struct_result.obj.sigmaFunc = simul.obj.sigmaFunc;
%             struct_result.obj.meanFunc = simul.obj.meanFunc;
            struct_result.gpSettingMaker = gpSettingMaker;
            struct_result.gpSettingsObj =gpSettingsObj ;
            struct_result.gpSettingsCon =gpSettingsCon;
            struct_result.estimator =estimator;
            if(param_BO.hasConstraints)
                struct_result.estimatorCon = estimatorCon;
            end
            struct_result.samplesize_cell = samplesize_cell;
            struct_result.samplesizeSelector = samplesizeSelector;
%             struct_result.
            
%             fprintf('> Optimal f Gap : %.4f\n', final_opt_f_true - log(optProb.objFunc.getOptVal() ));
%             fprintf('> Optimal x Gap : %.4f\n', getMinimumDistance(final_opt_x, optProb.objFunc.getOptSol()));
            results_all{idx_samplesize,repeat, idx_function, idx_setting} = struct_result;
    %         results_all
