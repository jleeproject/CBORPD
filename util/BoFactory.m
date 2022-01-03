classdef BoFactory
    %BOOBJECTFACTORY Summary of this class goes here
    %   Detailed explanation goes here
    
    methods (Static)
        %% AcquisitionFunction Objects
        function obj = getAcquisitionFunction(type, varargin)
            switch type
                case TypeAcquisition.EI
                    obj = AcqFuncEI(varargin{:});
                case TypeAcquisition.UCB
                    obj = AcqFuncUCB(varargin{:});
                case TypeAcquisition.LCB
                    obj = AcqFuncLCB(varargin{:});
                otherwise
                    error('Undefined Class')
            end
        end
        

        %% EstimatorGP Objects
        function obj = getEstimatorGP(type, varargin)
            switch type
                case TypeEstimatorGp.DirectLogSampleVar
                    obj = EstimatorGpDirectLogSampleVar(varargin{:});
                case TypeEstimatorGp.GpmlLogSampleVar
                    obj = EstimatorGpGpmlLogSampleVar(varargin{:});
                case TypeEstimatorGp.GpmlLogSampleVarRandomLS
%                     obj =  (varargin{:});
                case TypeEstimatorGp.DirectLogSampleVarRandomLS
                    obj = EstimatorGpDirectLogSampleVar_randomLS(varargin{:});
                case TypeEstimatorGp.DirectStochastic
                    obj = EstimatorGpDirectStochastic(varargin{:});
                case TypeEstimatorGp.DirectDeterministic
                    obj = EstimatorGpDirectDeterministic(varargin{:});
                case TypeEstimatorGp.GpmlStochastic
                    obj = EstimatorGpGpmlStochastic(varargin{:});
                case TypeEstimatorGp.GpmlDeterministic
                    obj = EstimatorGpGpmlDeterministic(varargin{:});
                case TypeEstimatorGp.DirectStochasticGivenNoiseVar
                    obj = EstimatorGpDirectStochasticNoiseVar(varargin{:});
                case TypeEstimatorGp.DirectLogSampleVar_SeparateLengthscale
                    obj = EstimatorGpDirectLogSampleVar_SeparateLengthscale(varargin{:});
                case TypeEstimatorGp.DirectGivenNoiseVariance_SeparateLengthscale
                    obj = EstimatorGpDirectStochasticNoiseVar_SeparateLengthscale(varargin{:});
%                 case TypeEstimatorGp.DirectGivenNoiseVariance_SeparateLengthscale
%                     obj = EstimatorGpDirectStochasticNoiseVar_SeparateLengthscale(varargin{:});
%                 case TypeEstimatorGp.GpmlStochasticGivenNoiseVar
%                     obj = EstimatorGpGpmlStochastic(varargin{:});
                otherwise
                    error(sprintf('Undefined Class:%s', type.char))
            end
        end
        
        %% Infill Optimizer Objects
        function obj = getInfillOptimizer(type, varargin)
            switch type
                case TypeInfillOptimizer.DirectUnconst
                    obj = InfillOptimizerDirectUnconst(varargin{:});
                case TypeInfillOptimizer.GridUnconst
                    obj = InfillOptimizerGridUnconst(varargin{:});
                case TypeInfillOptimizer.DirectConstSampling
                    obj = InfillOptimizerDirectConstSampling(varargin{:});
                case TypeInfillOptimizer.GridConstSampling
                    obj = InfillOptimizerGridConstSampling(varargin{:});
                case TypeInfillOptimizer.DirectConstWEI
                    obj = InfillOptimizerDirectConstWEI(varargin{:});
                case TypeInfillOptimizer.DirectConstWEI_f_opt_within_xtest
                    obj = InfillOptimizerDirectConstWEI_f_opt_within_xtest(varargin{:});
                case TypeInfillOptimizer.DirectConstWEI_f_opt_within_xtest_using_mean
                    obj = InfillOptimizerDirectConstWEI_f_opt_within_xtest_using_mean(varargin{:});
                case TypeInfillOptimizer.DirectConstNEI
                    obj = InfillOptimizerDirectConstNEI(varargin{:});
                case TypeInfillOptimizer.DirectConstTesting
                    obj = InfillOptimizerDirectConst_Testing(varargin{:});
                case TypeInfillOptimizer.DirectConstTesting2
                    obj = InfillOptimizerDirectConst_Testing2(varargin{:});
                case TypeInfillOptimizer.DirectConstTesting3
                    obj = InfillOptimizerDirectConst_Testing3(varargin{:});
                case TypeInfillOptimizer.DirectConstTesting4
                    obj = InfillOptimizerDirectConst_Testing4(varargin{:});
                case TypeInfillOptimizer.DirectConstTesting5
                    obj = InfillOptimizerDirectConst_Testing5(varargin{:});
                case TypeInfillOptimizer.DirectConstTesting6
                    obj = InfillOptimizerDirectConst_Testing6(varargin{:});
                case TypeInfillOptimizer.DirectConstTesting7
                    obj = InfillOptimizerDirectConst_Testing7(varargin{:});
                case TypeInfillOptimizer.DirectConstTesting8
                    obj = InfillOptimizerDirectConst_Testing8(varargin{:});
                case TypeInfillOptimizer.DirectConstTesting9
                    obj = InfillOptimizerDirectConst_Testing9(varargin{:});
                case TypeInfillOptimizer.DirectConstE2CBO
                    obj = InfillOptimizerDirectConstE2CBO(varargin{:});
                case TypeInfillOptimizer.DirectConstHWEI
                    obj = InfillOptimizerDirectConstHWEI(varargin{:});
                case TypeInfillOptimizer.DirectConstAdmm
                    %
                    obj = InfillOptimizerDirectConstSampling(varargin{:});
                case TypeInfillOptimizer.DirectConstAdmmNew
                    %
                    obj = InfillOptimizerDirectConstSampling(varargin{:});
                case TypeInfillOptimizer.DirectConstSobol
                    obj = InfillOptimizerDirectConstSobol(varargin{:});
                case TypeInfillOptimizer.DirectConstPatternSearch
                    obj = [];
                otherwise
                    error('Undefined Class')
            end
        end
        
        %% Initilizer Objects
        function obj = getInitializer(type, varargin)
            switch type
                case TypeInitializer.Lhd
                    obj = InitializerLHD(varargin{:});
                case TypeInitializer.UnifRnd
                    obj = InitializerUnifRnd(varargin{:});
                case TypeInitializer.Sobol
                    obj = InitializerSobol(varargin{:});
                otherwise
                    error('Undefined Class')
            end
        end
        
        %% Objective Function Objects
        function obj = getObjectiveFunction(type, varargin)
            switch type
                case TypeObjectiveFuncion.Branin2
                    obj = ObjFuncBranin2(varargin{:});
                case TypeObjectiveFuncion.Levy6
                    obj = ObjFuncLevy6(varargin{:});
                case TypeObjectiveFuncion.Shc2
                    obj = ObjFuncShc2(varargin{:});
                case TypeObjectiveFuncion.Thc2
                    obj = ObjFuncT(varargin{:});
                otherwise
                    error('Undefined Class')
            end
        end
        
        %% Simulator Objects
        function obj = getSimulator(type, varargin)
            switch type
                case TypeSimulator.LogSampleVar
                    obj = SimulatorLogSampleVar(varargin{:});
                case TypeSimulator.SimulatorStochasticMeanSigma
                    obj = SimulatorStochasticMeanSigma(varargin{:});
                case TypeSimulator.SimulatorDeterministic
                    obj = SimulatorDeterministic(varargin{:});
                case TypeSimulator.SimulatorGrapheneSensor
                    obj = SimulatorGrapheneSensor(varargin{:});
                case TypeSimulator.SimulatorGrapheneSameSensorOnOffRatio
                    obj = SimulatorGrapheneSameSensorOnOffRatio(varargin{:});
                case TypeSimulator.SimulatorGrapheneSameSensorResponseRatio
                    obj = SimulatorGrapheneSameSensorResponseRatio(varargin{:});
                otherwise
                    error('Undefined Class')
            end
        end
        
    end

end

