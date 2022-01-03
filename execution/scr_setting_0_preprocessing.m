clear all;close all;setting.isMaximizeObj = false; setting.power = 1;
A_RUN_PURPOSE = ...
    'RDOs';
%     '';
%     'ariafar19';
%     'gelbart14';

% =========================================================================
% setting.nRepeats = 50;
% setting.budget = 100;
% setting.nRepeats = 1;
% setting.budget = 100;
% 
% setting.isDebug = true;
% setting.isInfo = true;
% setting.isDebug = false;
% setting.isInfo = false;
% 
% % setting.max_opt_iter= 10;
% % setting.max_opt_iter= 4000;
% setting.max_opt_iter= 200;
%------------------------------------------------------
%% Debug
setting.nRepeats = 1;
setting.budget = 50;
% setting.isDebug = true;
setting.isInfo = true;
setting.isDebug = false;
% setting.isInfo = false;
setting.max_opt_eval= 1000;
setting.max_opt_iter= 1e5;
setting.show_plot_online = true;
% setting.show_plot_online = false;
setting.grid_bin = 25;

% %------------------------------------------------------
% %% Test at Server
% setting.nRepeats = 20;
% setting.budget = 100;
% setting.max_opt_iter= 500;
% setting.isDebug = false;setting.isInfo = false;setting.show_plot_online = false;

%------------------------------------------------------
% %% FAST MODE FOR DEBUG
% setting.nRepeats = 2;
% setting.budget = 20;
% setting.max_opt_iter= 10;
% setting.isInfo = true; setting.isDebug = false; setting.show_plot_online = false;

% %% Fast Plot
% setting.show_plot_online = true;
% setting.grid_bin = 5;

% =========================================================================
%% PLOT SETTING
% setting.show_plot_online = true;
% setting.show_plot_online = false;
% setting.show_default = true;
setting.show_default = false;
setting.show_exact_fmin = true;
% % setting.show_exact_fmin = false; % show difference between eval fmin and true fmin
% setting.grid_bin = 25;

setting.typeProblem = TypeProblem.ConventionalBoDeterministic;
% setting.typeProblem = TypeProblem.ConventionalBoStochastic;
% setting.typeProblem = TypeProblem.RobustDesignOptimization;
% setting.sigmaObj = .01;
% setting.sigmaCons = {.01};
setting.sigmaObj = .5;
setting.sigmaCons = {1};
% setting.typeProblem = TypeProblem.RobustDesignOptimization;

%% Sample sizes
setting.type_cell_samplesize = {...
%     {TypeSampleSize.Adjust,[]}
%     {TypeSampleSize.Adjust,[]} ...
%     ,...
%     {TypeSampleSize.Fixed,3} ...
    {TypeSampleSize.Fixed,5} ...
% %     , ...
%     {TypeSampleSize.Fixed,10} ...
% %     , ...
% %     {TypeSampleSize.Fixed,50}...
    };

%% Constrained Settings
setting.type_cell_settings = {...
% %         {TypeInfillOptimizer.DirectConstTesting2, TypeAcquisition.EI, TypeGpMean.Liberal}...
%         {TypeInfillOptimizer.DirectConstAdmm, TypeAcquisition.EI, TypeGpMean.Liberal}...
        {TypeInfillOptimizer.DirectConstWEI, TypeAcquisition.EI, TypeGpMean.Liberal}...
        , ...
        {TypeInfillOptimizer.DirectConstTesting8, TypeAcquisition.EI, TypeGpMean.Liberal}...
%         , ...
%         {TypeInfillOptimizer.DirectConstAdmmNew, TypeAcquisition.EI, TypeGpMean.Liberal}...
%         , ...
%         {TypeInfillOptimizer.DirectConstSampling, TypeAcquisition.LCB, TypeGpMean.Liberal}...
%         , ...
%         {TypeInfillOptimizer.DirectConstTesting, TypeAcquisition.EI, TypeGpMean.Liberal}...
%         , ...
%         {TypeInfillOptimizer.DirectConstTesting2, TypeAcquisition.EI, TypeGpMean.Liberal}...
%         {TypeInfillOptimizer.DirectConstTesting7, TypeAcquisition.EI, TypeGpMean.Liberal}...
%         , ...
%         {TypeInfillOptimizer.DirectConstTesting5, TypeAcquisition.EI, TypeGpMean.Liberal}...
%         , ...
%         {TypeInfillOptimizer.DirectConstTesting3, TypeAcquisition.EI, TypeGpMean.Liberal}...
%         , ...
%         {TypeInfillOptimizer.DirectConstTesting4, TypeAcquisition.EI, TypeGpMean.Liberal}...
%         , ...
%         {TypeInfillOptimizer.DirectConstTesting6, TypeAcquisition.EI, TypeGpMean.Liberal}...
%         , ...
%         {TypeInfillOptimizer.DirectConstTesting, TypeAcquisition.EI, TypeGpMean.Liberal}...
%         , ...
%         , ...
%         {TypeInfillOptimizer.DirectConstSampling, TypeAcquisition.LCB, TypeGpMean.Liberal}...
%         , ...
%         {TypeInfillOptimizer.DirectConstSampling, TypeAcquisition.EI, TypeGpMean.Liberal}...
%         , ...
%         {TypeInfillOptimizer.DirectConstSampling, TypeAcquisition.LCB, TypeGpMean.Conservative}...
%         , ...
%         {TypeInfillOptimizer.DirectConstWEI, TypeAcquisition.EI, TypeGpMean.Liberal}...
    };

% % Unconstrained
% setting.type_cell_settings = {...
%         {TypeInfillOptimizer.DirectUnconst, TypeAcquisition.LCB, TypeGpMean.Liberal}...
%         , ...
%         {TypeInfillOptimizer.DirectUnconst, TypeAcquisition.EI, TypeGpMean.Conservative}...
%         , ...
%         {TypeInfillOptimizer.DirectUnconst, TypeAcquisition.LCB, TypeGpMean.Conservative}...
% %         , ...
% %         {TypeInfillOptimizer.DirectConstWEI, TypeAcquisition.EI, TypeGpMean.Liberal}...
%     };


%% Problems
clearProblems;
import_problem_const_ariafar19;     addProblems;
% import_problem_const_gramacy;       addProblems;
% import_problem_const_gardener14_1;  addProblems;
% import_problem_const_gardener14_2;  addProblems;
% import_problem_const_gardener14_2_medium;  addProblems;

% import_problem_const_gelbart14;     addProblems;

% -----------------------------------------------------------------------
%% Temporarily do not use gardener's. They have too small variations over x

% % import_problem_RDO_const_gardener14_1;addProblems;
% import_problem_const_gardener14_2;addProblems;
% import_problem_RDO_const_gramacy;addProblems;
% import_problem_RDO_const_gelbart14;addProblems;
% import_problem_const_gardener14_1_easy;  addProblems;
% % import_problem_const_gardener14_1_medium;  addProblems;
% % import_problem_const_gardener14_1_hard;  addProblems;
% import_problem_const_gardener14_2;  addProblems;
% 
% import_problem_const_gramacy_easy;       addProblems;
% % import_problem_const_gramacy_medium;       addProblems;
% % import_problem_const_gramacy_hard;       addProblems;
% % 
% % 
import_problem_const_ariafar19_easy;     addProblems;
% % import_problem_const_ariafar19_medium;     addProblems;
% % import_problem_const_ariafar19_hard;     addProblems;
% % 
% import_problem_const_gelbart14_easy;     addProblems;
% % import_problem_const_gelbart14_medium;     addProblems;
% % import_problem_const_gelbart14_hard;     addProblems;
% 
% % import_problem_unconst_simple4func;     addProblems;
% % import_problem_custom;

% =========================================================================

% setting.type_gp_fit = TypeEstimatorGp.GpmlLogSampleVar;
setting.type_gp_fit = TypeEstimatorGp.DirectLogSampleVar;

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

setting.type_initial_sample_size = TypeInitialSampleSize.Fixed;
setting.n_init_sample_fixed = 2; %% if # inital sample points are fixed, this is used.
% setting.type_initial_sample_size = TypeInitialSampleSize.MultiplyDim10;

setting.type_simulator = TypeSimulator.LogSampleVar;

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

setting.samplesize_for_variance_init = 3;
% -----------------------------------------------------------

setting.casestudy.validation.samplesize = 100;
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

