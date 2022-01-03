optProb.name = cell_prob_names{idx_function};
sigmaObj = cell_sigmaObj{idx_function};
sigmaCons = cell_sigmaCons{idx_function};
setting.sigmaObj = sigmaObj;
setting.sigmaCons = sigmaCons;

% objFunc = cell_obj_func{idx_function};
% sigmaFunc = cell_obj_func{idx_function};
% meanFunc = cell_mean_func{idx_function};
optProb.objFunc = cell_obj_func{idx_function};
if(numel(cell_const_func{idx_function})>0)
    optProb.con = cell_const_func{idx_function};
    optProb.nCon =  optProb.con.nConstraints;
else
    optProb.nCon =  0;
end

if( isa(cell_obj_func{idx_function}, 'GrapheneModelSolver'))
%     disp('');
    grapheneVar = GrapheneDecisionVariables;
    optProb.xDomain = grapheneVar.domains.all;
    optProb.nVariable = size(optProb.xDomain,1);
    
    
    if (  setting.typeProblem ==  TypeProblem.RobustDesignOptimization)
            optProb.maximize = false;

%             simul = struct();
        simul = cell_obj_func{idx_function};
%             simul.obj.sigmaFunc = cell_obj_func{idx_function};
%             %% Temporary
%             optProb.conFunc = cell_mean_func{idx_function};
%             simul.obj.meanFunc = cell_mean_func{idx_function}{1};
% 
%             simul.nCon = 0;
            if(~contains(setting.type_gp_fit.char , 'LogSampleVar'))
                error(sprintf('[WARNING] Wrong type of GP fit : %s', setting.type_gp_fit.char));
            end
    %         if(~contains(setting.type_simulator.char , 'LogSampleVar'))
    %             error('[WARNING] Wrong type of simulator. Use LogSampleVar.');
    %         else
    %         end
        if(~contains(setting.type_simulator.char, 'Graphene'))
            error('Wrong Simulator');
        end
%             setting.type_simulator = TypeSimulator.SimulatorGrapheneSensor;
    else
        error('[scr_setting1_in_functions_sameObj] Graphene sensor model is only compatible with RDO');
    end
    
else
    optProb.xDomain = optProb.objFunc.getXDomain;
    optProb.nVariable = optProb.objFunc.getDim;

    switch setting.typeProblem
        case TypeProblem.RobustDesignOptimization
            optProb.maximize = false;

            simul.obj.sigmaFunc = cell_obj_func{idx_function};
            %% Temporary
            optProb.conFunc = cell_mean_func{idx_function};
            simul.obj.meanFunc = cell_mean_func{idx_function}{1};

            simul.nCon = 0;
            if(~contains(setting.type_gp_fit.char , 'LogSampleVar'))
                error('[WARNING] Wrong type of GP fit');
            end
    %         if(~contains(setting.type_simulator.char , 'LogSampleVar'))
    %             error('[WARNING] Wrong type of simulator. Use LogSampleVar.');
    %         else
    %         end
            setting.type_simulator = TypeSimulator.SimulatorStochasticMeanSigma;

        case TypeProblem.ConventionalBoDeterministic
            optProb.maximize = setting.isMaximizeObj;

            simul.obj.sigmaFunc = 0;
            simul.obj.meanFunc = cell_obj_func{idx_function};

    %         cell_mean_func{idx_function};

            simul.nCon = cell_const_func{idx_function}.nConstraints;
    %         if(simul.nCon ~= numel(cell_mean_func{idx_function}))
    %             error('Dimension mismatch');
    %         end
            simul.con = cell(simul.nCon,1);
            for idx_con = 1: simul.nCon
                simul.con{idx_con}.meanFunc = cell_mean_func{idx_function}{idx_con};
                simul.con{idx_con}.sigmaFunc = 0;
            end
            optProb.conFunc = cell_mean_func{idx_function};

            setting.type_simulator = TypeSimulator.SimulatorDeterministic;
            if(contains(setting.type_gp_fit.char , 'LogSampleVar'))
                disp('[WARNING] Wrong type of GP fit: change into stochastic one.');
                if(contains(setting.type_gp_fit.char , 'Direct'))
                    setting.type_gp_fit = TypeEstimatorGp.DirectDeterministic;
                elseif(contains(setting.type_gp_fit.char , 'Gpml'))
                    setting.type_gp_fit = TypeEstimatorGp.GpmlDeterministic;
                else
                    throwUndefinedTypeError(setting.type_gp_fit)
                end
    %         elseif(contains(setting.type_simulator.char , 'LogSampleVar'))
    %             setting.type_simulator = TypeSimulator.SimulatorDeterministic;
    %             fprintf('[WARNING] Wrong type of simulator. Use Conventional BO simulator.');
    %         else
            end

            %% SAMPLE SIZE = 1;
            setting.type_cell_samplesize = {{TypeSampleSize.Fixed,1}};
        case TypeProblem.ConventionalBoStochastic
            optProb.maximize = setting.isMaximizeObj;

            simul.obj.sigmaFunc = setting.sigmaObj;
            simul.obj.meanFunc = cell_obj_func{idx_function};

    %         cell_mean_func{idx_function};

            simul.nCon = cell_const_func{idx_function}.nConstraints;
            if(simul.nCon ~= numel(cell_mean_func{idx_function}))
                error('Dimension mismatch');
            end
            simul.con = cell(simul.nCon,1);
            for idx_con = 1: simul.nCon
                simul.con{idx_con}.meanFunc = cell_mean_func{idx_function}{idx_con};
                if(numel(setting.sigmaCons)>1)
                    simul.con{idx_con}.sigmaFunc = setting.sigmaCons{idx_con};
                else
                    simul.con{idx_con}.sigmaFunc = setting.sigmaCons;
                end
            end
            %% Temporary
            optProb.conFunc = cell_mean_func{idx_function};

            if(contains(setting.type_gp_fit.char , 'LogSampleVar'))
                disp('[WARNING] Wrong type of GP fit: change into stochastic one.');
                if(contains(setting.type_gp_fit.char , 'Direct'))
                    setting.type_gp_fit = TypeEstimatorGp.DirectStochastic;
                elseif(contains(setting.type_gp_fit.char , 'Gpml'))
                    setting.type_gp_fit = TypeEstimatorGp.GpmlStochastic;
                else
                    throwUndefinedTypeError(setting.type_gp_fit)
                end
            elseif(contains(setting.type_simulator.char , 'LogSampleVar'))
                setting.type_simulator = TypeSimulator.SimulatorStochasticMeanSigma;
                fprintf('[WARNING] Wrong type of simulator. Use Conventional BO simulator.');
            else
            end

            %% SAMPLE SIZE = 1;
            setting.type_cell_samplesize = {{TypeSampleSize.Fixed,1}};
            setting.type_simulator = TypeSimulator.SimulatorStochasticMeanSigma;

        otherwise
            throwUndefinedTypeError(setting.typeProblem);
    end
end







setting.optProb = optProb;
setting.simul = simul;

switch setting.type_initial_sample_size
    case TypeInitialSampleSize.MultiplyDim10
        n_init_sample = optProb.nVariable * 10; % Will be used on.
    case TypeInitialSampleSize.Fixed
        n_init_sample = setting.n_init_sample_fixed;
    otherwise
        throwError('Undefined type_initial_sample_size');
end
setting.n_init_sample = n_init_sample;


xDomain = optProb.xDomain;
simulator = BoFactory.getSimulator(setting.type_simulator, simul, optProb, setting.typeProblem);
initializer = BoFactory.getInitializer(setting.type_initial_sampler, simulator, xDomain, n_init_sample, setting.optProb.nCon, setting.samplesize_for_variance_init);



setting_default = ObjectCopier.copy(setting);



