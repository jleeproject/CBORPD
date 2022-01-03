close all;setting.isMaximizeObj = false; setting.power = 1;
A_RUN_PURPOSE = ...
    'Casestudy_Graphene_ub_05_b3000_ev200';

%% Debug
setting.budget = 3000;
setting.max_opt_eval= 1000;
% setting.casestudy.validation.samplesize = 1000;
setting.casestudy.validation.samplesize = 200;
% setting.casestudy.validation.samplesize = 10;
setting.n_init_sample_fixed = 100; %% if # inital sample points are fixed, this is used.
% setting.n_init_sample_fixed = 5; %% if # inital sample points are fixed, this is used.
% setting.type_initial_sample_size = TypeInitialSampleSize.MultiplyDim10;
% setting.samplesize_for_variance_init = 10 ;
setting.samplesize_for_variance_init = 2 ;

% setting.typeProblem = TypeProblem.ConventionalBoDeterministic;
% setting.typeProblem = TypeProblem.ConventionalBoStochastic;
setting.typeProblem = TypeProblem.RobustDesignOptimization;
% setting.sigmaObj = .01;
% setting.sigmaCons = {.01};
setting.sigmaObj = .5;
setting.sigmaCons = {1};
% setting.typeProblem = TypeProblem.RobustDesignOptimization;


%% Sample sizes
setting.type_cell_samplesize = {...
    {TypeSampleSize.Adjust,[]} ...
%     {TypeSampleSize.Adjust,[]}
%     ,...
%     {TypeSampleSize.Fixed,1} ...
%     {TypeSampleSize.Fixed,3} ...
% %     , ...
%     {TypeSampleSize.Fixed,10} ...
% %     , ...
%     {TypeSampleSize.Fixed,50}...
%     {TypeSampleSize.Fixed,100}...
    };

% setting.type_simulator = TypeSimulator.LogSampleVar;
% setting.type_simulator = TypeSimulator.SimulatorGrapheneSameSensorOnOffRatio;
setting.type_simulator = TypeSimulator.SimulatorGrapheneSameSensorResponseRatio;


%% Problems
clearProblems;
problem_name = ...
    'Graphene Sensor Model';

% multiplier2Mean = 10;

isMaximizeObj = false; 
% power = 1;
cell_obj_func = {...
    GrapheneModelSolver(  GrapheneSimulFixedParameters(), GrapheneSimulChangingParameters(TypeGrapheneSensingSolution.WATER_C0) , GrapheneSimulChangingParameters(TypeGrapheneSensingSolution.WATER_C20)	)  ...
    };

cell_type_mean_func = {...
    {setting.type_simulator}...
%     {TypeSimulator.SimulatorGrapheneSensor}...
    };
cell_const_func = {...
    MultipleConstraints(...
        {Constraint(TypeConstraint.UnknownMeanConstraint, .05, [] ) }...
    )...
};
addProblems;
% import_problem_const_ariafar19;     addProblems;
% import_problem_const_gramacy;       addProblems;
% import_problem_const_gardener14_1;  addProblems;
% import_problem_const_gardener14_2;  addProblems;
% import_problem_const_gardener14_2_medium;  addProblems;

% import_problem_const_gelbart14;     addProblems;

% -----------------------------------------------------------------------
%% Temporarily do not use gardener's. They have too small variations over x


% % import_problem_RDO_const_gardener14_1;addProblems;
% % import_problem_const_gardener14_2;addProblems;
% % import_problem_RDO_const_gramacy;addProblems;
% % import_problem_RDO_const_gelbart14;addProblems;
% % import_problem_const_gardener14_1_easy;  addProblems;
% import_problem_const_gardener14_1_medium;  addProblems;
% % % import_problem_const_gardener14_1_hard;  addProblems;
% 
% import_problem_const_gelbart14_easy;     addProblems;
% % % import_problem_const_gelbart14_medium;     addProblems;
% % % import_problem_const_gelbart14_hard;     addProblems;
% 
% import_problem_const_gramacy_easy;       addProblems;
% % % import_problem_const_gramacy_medium;       addProblems;
% % % import_problem_const_gramacy_hard;       addProblems;
% 
% import_problem_const_gardener14_2;  addProblems;
% % 
% % % 
% % % 
% import_problem_const_ariafar19_easy;     addProblems;
% % % import_problem_const_ariafar19_medium;     addProblems;
% % % import_problem_const_ariafar19_hard;     addProblems;
% % % 
% % 
% % % import_problem_unconst_simple4func;     addProblems;
% % % import_problem_custom;

% =========================================================================


setting.type_initial_sample_size = TypeInitialSampleSize.Fixed;
