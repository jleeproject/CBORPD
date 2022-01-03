close all;setting.isMaximizeObj = false; setting.power = 1;
A_RUN_PURPOSE = ...
    'Casestudy';

%% Debug
setting.budget = 500;
setting.max_opt_eval= 1000;
setting.casestudy.validation.samplesize = 100;
setting.n_init_sample_fixed = 2; %% if # inital sample points are fixed, this is used.
% setting.type_initial_sample_size = TypeInitialSampleSize.MultiplyDim10;
setting.samplesize_for_variance_init = 2 ;

setting.typeProblem = TypeProblem.RobustDesignOptimization;
setting.sigmaObj = .5;
setting.sigmaCons = {1};


%% Sample sizes
setting.type_cell_samplesize = {...
    {TypeSampleSize.Fixed,setting.samplesize_for_variance_init} ...
    };

setting.type_simulator = TypeSimulator.SimulatorGrapheneSameSensorResponseRatio;


%% Problems
clearProblems;
problem_name = ...
    'Graphene Sensor Model';


isMaximizeObj = false; 
% power = 1;
cell_obj_func = {...
    GrapheneModelSolver(  GrapheneSimulFixedParameters(), GrapheneSimulChangingParameters(TypeGrapheneSensingSolution.WATER_C0) , GrapheneSimulChangingParameters(TypeGrapheneSensingSolution.WATER_C20)	)  ...
    };

cell_type_mean_func = {...
    {setting.type_simulator}...
    };
cell_const_func = {...
    MultipleConstraints(...
        {Constraint(TypeConstraint.UnknownMeanConstraint, .10, [] ) }...
    )...
};
addProblems;


setting.type_initial_sample_size = TypeInitialSampleSize.Fixed;
