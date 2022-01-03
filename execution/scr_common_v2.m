% setting.plot.savefile = true;
setting.plot.savefile = false;

% setting.nRepeats = 100;
setting.nRepeats = 2;


setting.show_default = false;
setting.show_exact_fmin = true;


%% Constrained Settings
setting.type_cell_settings = {...
        {TypeInfillOptimizer.DirectConstE2CBO}...
        , ...
        {TypeInfillOptimizer.DirectConstPatternSearch}...
        , ...
        {TypeInfillOptimizer.DirectConstSobol}...
        , ...
        {TypeInfillOptimizer.DirectConstHWEI}...
        , ...
        {TypeInfillOptimizer.DirectConstNEI}...
    };




% setting.type_gp_fit = TypeEstimatorGp.GpmlLogSampleVar;
% setting.type_gp_fit = TypeEstimatorGp.DirectLogSampleVar;
% setting.type_gp_fit = TypeEstimatorGp.DirectGivenNoiseVariance_SeparateLengthscale;
setting.type_gp_fit = TypeEstimatorGp.DirectLogSampleVar_SeparateLengthscale;

% setting.type_infill_opt = TypeInfillOptimizer.GridUnconst;
% setting.type_infill_opt = TypeInfillOptimizer.DirectUnconst;
% setting.type_infill_opt = TypeInfillOptimizer.GridConstSampling;
setting.type_infill_opt = TypeInfillOptimizer.DirectConstSampling;
% setting.type_infill_opt = TypeInfillOptimizer.DirectConstWEI;

% setting.type_acqusition = TypeAcquisition.EI;
setting.type_acquisition = TypeAcquisition.LCB;

% setting.type_obj_func = TypeObjectiveFuncion.Branin2;
% setting.type_obj_func = TypeObjectiveFuncion.Thc2;
% setting.type_obj_func = TypeObjectiveFuncion.Shc2;
% setting.type_obj_func = TypeObjectiveFuncion.Levy6;

setting.type_initial_sampler = TypeInitializer.Lhd;
% setting.type_initial_sampler = TypeInitializer.UnifRnd;




% setting.type_gp_kernel = TypeGpKernel.Matern52;
setting.type_gp_kernel = TypeGpKernel.Se;

% =========================================================================
setting.is_conservative = true;
% setting.is_conservative = false;

% setting.meanFuncConservative = @(meanV, minV, maxV, rangeV) minV - 2* rangeV;
%     setting.meanFuncConservative = @(meanV, minV, maxV, rangeV) minV - 1* rangeV;

setting.wantDetailedOptimalValues = true;
% setting.max_eval_diRect = 1000;
%% ----- Visualization ------------------------------------------------------

% setting.show_default = true;
% % setting.show_default = false;

% if show_default is false, follow below...
setting.show_all_grid_visualization = true;
% setting.show_all_grid_visualization = false;

setting.show_all_search_trails = true;
% setting.show_all_search_trails = false;


% -----------------------------------------------------------
% -----------------------------------------------------------

% -------------------------------------------------------------
% ----- For algorithm. Not user pecified ------------------------------------------------------------

if(setting.typeProblem ~= TypeProblem.RobustDesignOptimization)
    setting.samplesize_for_variance_init = 1;
end

cell_mean_func = cell(numel(cell_type_mean_func),1);
for idx_prob = 1:size(cell_type_mean_func,1)
    setting_mean_func = cell_type_mean_func{idx_prob};
    nCons = numel(setting_mean_func);
    cell_sel_prob_mean_func = cell(numel(setting_mean_func),1);
    for idx_con = 1:nCons
        sel_con = setting_mean_func{idx_con};
%         setting_mean_func = cell_type_mean_func{idx_prob};
        if(numel(sel_con)>1)
            type_mean_func = sel_con{1};
            mod_mean_func = sel_con{2};
            cell_sel_prob_mean_func{idx_con} = FunctionFactory.getFunction(type_mean_func,[],[], mod_mean_func);
        elseif(numel(sel_con)==1)
            type_mean_func = getElementFromArrayOrCell(sel_con,1);
            mod_mean_func = [];
            if( isa(type_mean_func ,'TypeFunction') )
                cell_sel_prob_mean_func{idx_con} = FunctionFactory.getFunction(type_mean_func,[],[], mod_mean_func);
            elseif( isa(cell_obj_func{idx_prob}, 'GrapheneModelSolver'))
                % Special case. No allocation
            else
                error('Unknown type of mean function')
            end
        else
            type_mean_func = [];
            mod_mean_func = [];
            cell_sel_prob_mean_func{idx_con} = [];
        end
    end
    cell_mean_func{idx_prob} = cell_sel_prob_mean_func;
end




if(exist('type_cell_const_func','var') && numel(cell_const_func)>0)
    if(numel(cell_const_func) ~= numel(cell_obj_func))
        error('The number of constraints are different from that of objective functions');
    end
    
end


nRepeats = setting.nRepeats;

type_cell_samplesize = setting.type_cell_samplesize;
type_cell_settings = setting.type_cell_settings;
% - For algorithm. Not user pecified------------------------------------------------------------
if(~exist('type_cell_settings','var') || numel(setting.type_cell_settings) ==0)
    nTypeSettings = 1;
    setting.type_cell_settings = [];
else
    nTypeSettings = numel(setting.type_cell_settings);
end




if(setting.typeProblem ~= TypeProblem.RobustDesignOptimization)
    setting.type_cell_samplesize = {1};
end

results_all = cell(numel(setting.type_cell_samplesize), setting.nRepeats, numel(cell_obj_func), nTypeSettings);
fprintf('Result Dimension:');
fprintf('%d ',size(results_all));
fprintf('\n');
