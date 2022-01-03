%% ----- Initial Samples ----------------------------------------------------
samplesize_cell = type_cell_samplesize{idx_samplesize};
typeSamplesize = samplesize_cell{1};
givenSamplesize = samplesize_cell{2};
if(setting.typeProblem == TypeProblem.RobustDesignOptimization)
    if(typeSamplesize == TypeSampleSize.Fixed && givenSamplesize==1)
        error('[ERROR : scr_setting_4_in_samples] RDO cannot have fixed sample size 1.');
    end
end

% --- Save Intermediate Results ----
storage = ResultStorageBO(param_BO, setting.optProb.objFunc, typeSamplesize, givenSamplesize);
storage.saveXHistory(xx_arr, 1:n_init_sample);
storage.saveYHistory(yy_arr, 1:n_init_sample);
if(param_BO.hasConstraints)
    storage.saveYConsHistory(yy_cons_arr, 1:n_init_sample);
end
storage.saveObsHistory(y_orig_cell, 1:n_init_sample);

storage.saveSamplesizeHistory(samplesize_arr, 1:n_init_sample);
nEvaluations = n_init_sample;



samplesizeSelector = SampleSizeSelector(typeSamplesize, givenSamplesize);

% --- For algorithm. Not user specific ----
stop = false;





param_BO = ObjectCopier.copy(param_BO);

%% SAVE Class Objects
simulator = ObjectCopier.copy(simulator);
initializer = ObjectCopier.copy(initializer);
acqFunc = ObjectCopier.copy(acqFunc);
infillOptimizer = ObjectCopier.copy(infillOptimizer);
estimator = ObjectCopier.copy(estimator);
gpSettingMaker = ObjectCopier.copy(gpSettingMaker);
% if(setting.show_plot_online)
%     visualizer = ObjectCopier.copy(visualizer);
% end
% objFunc = ObjectCopier.copy(objFunc);
% sigmaFunc = ObjectCopier.copy(sigmaFunc);
% meanFunc = ObjectCopier.copy(meanFunc);

if(param_BO.hasConstraints);
	gpSettingsObj = ObjectCopier.copy(gpSettingsObj);
	gpSettingsCon = ObjectCopier.copy(gpSettingsCon);
else
    clear gpSettingsObj gpSettingsCon;
end





%% Visualizer.
if(setting.show_plot_online )
    visualizer = Visualizer_BO();
    if(setting.show_default)
        switch type_infill_opt
            case TypeInfillOptimizer.DirectConstSampling
                visualizer.setShowAllGrid(false);
                visualizer.setShowSearchTrails(true);
            case TypeInfillOptimizer.DirectUnconst
                visualizer.setShowAllGrid(false);
                visualizer.setShowSearchTrails(true);
            case TypeInfillOptimizer.GridConstSampling
                visualizer.setShowAllGrid(true);
                visualizer.setShowSearchTrails(false);
            case TypeInfillOptimizer.GridUnconst
                visualizer.setShowAllGrid(true);
                visualizer.setShowSearchTrails(false);
            otherwise
                disp('[WARNING] Unknown Settings: Show no grid, no trail.');
                visualizer.setShowAllGrid(false);
                visualizer.setShowSearchTrails(false);
%                 show_default = false;
        end
    end
    if(~setting.show_default )
        visualizer.setShowSearchTrails(setting.show_all_search_trails);
        visualizer.setShowAllGrid(setting.show_all_grid_visualization);
    end
end


if(setting.show_plot_online)
    visualizer.setFnEval(setting.optProb.objFunc)
    if(setting.show_all_search_trails)
        visualizer.setGrid(param_BO.getGridBin(), setting.optProb.objFunc)
    end
    if(~visualizer.isReady)
        throwError('Visualizer NOT READY.\n');
        return;
    end
end
%%