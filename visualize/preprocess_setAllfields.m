
param_BO = Param_BO_module(type_gp_fit, type_infill_opt, budget, nRepeats);
param_BO.setConservative(is_conservative);

param_BO.setInitSmallSamplesize();
param_BO.setPlotOnline(show_plot_online);
param_BO.setMaxEvalDiRect(max_opt_iter);
param_BO.setDebug(isDebug);
param_BO.setInfo(isInfo);

param_BO.setGridBin(grid_bin);
param_BO.setTypeAcquisition(type_acquisition)
if(show_plot_online)
    visualizer = Visualizer_BO(param_BO.options_plot);
    if(show_default)
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
                show_default = false;
        end
    end
    if(~show_default)
        if(show_all_search_trails)
            visualizer.setShowSearchTrails(true);
        else
            visualizer.setShowSearchTrails(false);
        end
        if(show_all_grid_visualization)
            visualizer.setShowAllGrid(true);
        else
            visualizer.setShowAllGrid(false);
        end
    end
end


% - For algorithm. Not user pecified------------------------------------------------------------
if(~exist('type_cell_settings','var') || numel(type_cell_settings) ==0)
    nTypeSettings = 1;
    type_cell_settings = [];
else
    nTypeSettings = numel(type_cell_settings);
end

% results_all = cell(numel(type_cell_samplesize), nRepeats, numel(cell_obj_func), nTypeSettings);
% fprintf('Result Dimension:');
% fprintf('%d ',size(results_all));
% fprintf('\n');


% max_iter_num = 10000;
param_BO_proto = param_BO;

cell_mean_func = cell(size(cell_type_mean_func));
for i = 1:size(cell_type_mean_func,1)
    for j = 1:size(cell_type_mean_func,2)
        setting_mean_func = cell_type_mean_func{i,j};
        if(numel(setting_mean_func)>1)
            type_mean_func = setting_mean_func{1};
            mod_mean_func = setting_mean_func{2};
            cell_mean_func{i,j} = FunctionFactory.getFunction(type_mean_func,[],[], mod_mean_func);
        else
            type_mean_func = [];
            mod_mean_func = [];
            cell_mean_func{i,j} = [];
        end
    end
end

if(exist('type_cell_const_func','var') && numel(cell_const_func)>0)
    if(numel(cell_const_func) ~= numel(cell_obj_func))
        throwError('The number of constraints are different from that of objective functions');
    end
    
end