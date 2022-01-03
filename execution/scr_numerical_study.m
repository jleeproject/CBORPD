% clear all;
close all;setting.isMaximizeObj = false; setting.power = 1;
A_RUN_PURPOSE = ...
    'NumericalStudy';

%% Debug
setting.type_budget = TypeBudget.MultiplyDimBy;
setting.budget = 100;
setting.max_opt_eval= 500;

setting.samplesize_for_variance_init = 2;

setting.typeProblem = TypeProblem.RobustDesignOptimization;

%% Sample sizes
setting.type_cell_samplesize = {...
    {TypeSampleSize.Fixed,setting.samplesize_for_variance_init} ...
    };


%% Problems
clearProblems;

% -----------------------------------------------------------------------
import_problem_gardner;    addProblems;
import_problem_gardner2;        addProblems;
import_problem_gelbart;       addProblems;
import_problem_hartmann6;     addProblems;


setting.type_initial_sample_size = TypeInitialSampleSize.Fixed;
setting.n_init_sample_fixed = 2; %% if # inital sample points are fixed, this is used.

