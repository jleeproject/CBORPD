% clear all;
close all;setting.isMaximizeObj = false; setting.power = 1;
A_RUN_PURPOSE = ...
    'Sto_NEI_Challenge_sigx10';

%% Debug
setting.budget = 300;
% setting.budget = 10;
setting.max_opt_eval= 700;


% setting.typeProblem = TypeProblem.ConventionalBoDeterministic;
setting.typeProblem = TypeProblem.ConventionalBoStochastic;
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
    {TypeSampleSize.Fixed,1} ...
%     {TypeSampleSize.Fixed,3} ...
% %     , ...
%     {TypeSampleSize.Fixed,10} ...
% %     , ...
% %     {TypeSampleSize.Fixed,50}...
    };



%% Problems
clearProblems;
% import_problem_const_ariafar19;     addProblems;
% import_problem_const_gramacy;       addProblems;
% import_problem_const_gardener14_1;  addProblems;
% import_problem_const_gardener14_2;  addProblems;
% import_problem_const_gardener14_2_medium;  addProblems;

% import_problem_const_gelbart14;     addProblems;

% -----------------------------------------------------------------------
%% Temporarily do not use gardener's. They have too small variations over x
import_problem_const_gramacy_letham_medium;       addProblems;
import_problem_const_ariafar19_letham_medium;     addProblems;
import_problem_const_gelbart14_letham_medium;     addProblems;
import_problem_const_gardener14_1_letham_medium;  addProblems;
import_problem_const_ariafar19_letham_variant_medium;     addProblems;

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
% setting.n_init_sample_fixed = 2; %% if # inital sample points are fixed, this is used.
setting.type_initial_sample_size = TypeInitialSampleSize.MultiplyDim10;

setting.samplesize_for_variance_init = 1;
