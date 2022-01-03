classdef Param_BO_v2<handle & matlab.mixin.Copyable
    
    properties
        typeGpFit
        typeInfillOptimization
        typeAcquisition;
        typeInitialSampler
        typeSimulator
        typeGpKernel
        %
        nInitSample = -1;
        is_conservative
        nRepeats
        budget
        %
        %
        isTrueFunctionSet = false;
        
        objFunc
        %
        debug
        info
        plotOnline = false;
        %
        initSmallSamplesize = true;
        samplesize_for_init_stage = 3;
        %
        gridBin = 100;
        % diRect
        maxOptIter = 1e4;
        maxOptEval = 200
        %
        hasConstraints = false;
        constraints = [];
        nConstraints = 0;
        
        typeProblem
        setting
        
        gpSettingMaker 
    end
    
    methods
        function this = Param_BO_v2(setting)            %PARAMETER Construct an instance of this class
%             this.typeGpFit = typeGpFit;
%             this.typeInfillOptimization = typeInfillOptimization;
%             this.nInitSample = nInitSample;
            this.budget = setting.budget;
            this.nRepeats = setting.nRepeats;
            

            this.setting = setting;            
%             this.isMaximizeObj = setting.isMaximizeObj;
%             this.power = setting.power;
            this.nRepeats = setting.nRepeats;
            this.budget = setting.budget;
            this.debug = setting.isDebug;
            this.info = setting.isInfo;
            this.plotOnline = setting.show_plot_online;
%             this.showDefault = setting.show_default;
%             this.showExactFmin = setting.show_exact_fmin;
            this.maxOptEval = setting.max_opt_eval;
            this.maxOptIter = setting.max_opt_iter;
%             setting.max_opt_eval= 500;
%               setting.max_opt_iter= 100000;
            this.samplesize_for_init_stage = setting.samplesize_for_variance_init;

            this.typeGpFit = setting.type_gp_fit;
            this.typeInfillOptimization = setting.type_infill_opt;
            this.typeAcquisition = setting.type_acquisition;
            this.typeInitialSampler = setting.type_initial_sampler;
%             this.typeInitialSampleSize = setting.type_initial_sample_size;
            this.typeSimulator = setting.type_simulator;
            this.typeGpKernel = setting.type_gp_kernel;
            this.is_conservative = setting.is_conservative;
%             this.wantDetailedOptimalValues = setting.wantDetailedOptimalValues;
%             this.showAllGridVisualization = setting.show_all_grid_visualization;
%             this.showAllSearchTrails = setting.show_all_search_trails;
            this.gridBin = setting.grid_bin;
            this.nInitSample = setting.n_init_sample;
            this.typeProblem = setting.typeProblem;
            
            this.objFunc = setting.optProb.objFunc;
            if(setting.optProb.nCon>0)
                this.setConstraints(setting.optProb.con);
            else
                this.setConstraints([]);
            end
      
            init(this);
        end
        
        function init(this)
        end
        
        function setConstraints(this, constraints);
            this.constraints = constraints;
            if(numel(constraints)>0)
                this.hasConstraints = true;
                this.nConstraints = constraints.nConstraints;
            else
                this.hasConstraints = false;
                this.nConstraints = 0;
            end
        end
        
        function out = getConstraints(this)
            out =  this.constraints;
        end

        function logic = isGpFitGpml(this)
%             logic = (this.typeGpFit==TypeEstimatorGp.GpmlLogSampleVar)...
%                 || (this.typeGpFit==TypeEstimatorGp.GpmlLogSampleVarRandomLS);
            logic = contains(this.typeGpFit.char, 'Gpml');
        end
        
        function logic = isGpFitDirect(this)
%             logic = (this.typeGpFit==TypeEstimatorGp.DirectLogSampleVar)...
%                 || (this.typeGpFit==TypeEstimatorGp.DirectLogSampleVarRandomLS);
            logic = contains(this.typeGpFit.char, 'Direct');
        end
        
        function logic = isInfillSearchGrid(this)
%             logic = (this.typeInfillOptimization==TypeInfillOptimizer.GridConstSampling)...
%                 ||(this.typeInfillOptimization==TypeInfillOptimizer.GridUnconst);
            logic = contains(this.typeInfillOptimization.char, 'Grid');
        end
        
        function logic = isInfillSearchDirect(this)
%             logic = (this.typeInfillOptimization==TypeInfillOptimizer.DirectConstSampling)...
%                 ||(this.typeInfillOptimization==TypeInfillOptimizer.DirectUnconst);
            logic = contains(this.typeInfillOptimization.char, 'Direct');
        end
        
        function setConservative(this,inputArg)
            if(nargin>1)
                this.is_conservative = inputArg;
            else
                this.is_conservative = true;
            end
        end
        function logic = isConservative(this)
            logic = (this.is_conservative);
        end
    
      
        function output = getEstUpperLimitNumIter(this, typeSamplesize, givenSamplesize)
            if(this.typeProblem == TypeProblem.RobustDesignOptimization)
                switch typeSamplesize
                    case TypeSampleSize.Adjust
                        output = floor(this.budget / 2) + this.nInitSample;
                    case TypeSampleSize.Fixed
                        output = floor(this.budget / givenSamplesize) + this.nInitSample;
                    otherwise
                end
            else
                output = this.budget + 1;
            end
        end
        
        function logic = isReady(this)
            logic = true;
            % If not satisfied one of below, not ready
            % Is true function set?
            if(~this.isTrueFunctionSet)
                logic = false;
            end
            if(this.nInitSample<=0)
                logic = false;
            end
        end
        
        function logic = isGridSearch(this)
            logic = (this.typeInfillOptimization == TypeInfillOptimizer.GridConstSampling )...
                ||(this.typeInfillOptimization == TypeInfillOptimizer.GridUnconst );
        end
%% =========================================================================
        %% Related to True Functions            
        function setDebug(this,debug)
            this.debug = debug;
        end

        function output = isDebug(this)
            output = this.debug;
        end


        function setInfo(this,info)
            this.info = info;
        end

        function output = isInfo(this)
            output = this.info;
        end
        
        
%% =========================================================================
%% ===== GETTERS AND SETTERS ===============================================
%% =========================================================================
%--------------------- [samplesizeForInitStage] ---------------------
function setSamplesizeForInitStage(this,samplesizeForInitStage)
    this.samplesizeForInitStage = samplesizeForInitStage;
    if(samplesizeForInitStage>3)
        this.initSmallSamplesize = false;
    end
end

function output = getSamplesizeForInitStage(this)
    output = this.samplesize_for_init_stage;
end

%--------------------- [initSmallSamplesize] ---------------------
function setInitSmallSamplesize(this,initSmallSamplesize)
    if(nargin>1)
        this.initSmallSamplesize = initSmallSamplesize;
    else
        this.initSmallSamplesize = true;
    end
    if(this.initSmallSamplesize)
        this.samplesize_for_init_stage = 3;
    end
end

function output = isInitSmallSamplesize(this)
    output = this.initSmallSamplesize;
end


function setObjFunc(this,objFunc)
    this.objFunc = objFunc;
end

function out = getObjFunc(this)
    out = this.objFunc;
end
%--------------------- [budget] ---------------------
function setBudget(this,budget)
    this.budget = budget;
end

function output = getBudget(this)
    output = this.budget;
end

%--------------------- [strFnName] ---------------------
% function setStrFnName(this,strFnName)
%     this.strFnName = strFnName;
% end
% 
% function output = getStrFnName(this)
%     output = this.str_fn_name;
% end

%--------------------- [nInitSample] ---------------------
function setNInitSample(this,nInitSample)
    this.nInitSample = nInitSample;
end

function output = getNInitSample(this)
    output = this.nInitSample;
end

%--------------------- [plotOnline] ---------------------
function setPlotOnline(this,plotOnline)
    this.plotOnline = plotOnline;
end

function output = isPlotOnline(this)
    output = this.plotOnline;
end

%--------------------- [gridBin] ---------------------
function setGridBin(this,gridBin)
    this.gridBin = gridBin;
end

function output = getGridBin(this)
    output = this.gridBin;
end

%--------------------- [typeAcquisition] ---------------------
function setTypeAcquisition(this,typeAcquisition)
    this.typeAcquisition = typeAcquisition;
end

function setTypeAcquisitionEI(this)
    this.typeAcquisition = TypeAcquisition.EI;
end

function setTypeAcquisitionUCB(this)
    this.typeAcquisition = TypeAcquisition.UCB;
end

function setTypeAcquisitionLCB(this)
    this.typeAcquisition = TypeAcquisition.LCB;
end

function output = isTypeAcquisitionEI(this)
    output = (this.typeAcquisition == TypeAcquisition.EI);
end

function output = isTypeAcquisitionUCB(this)
    output = (this.typeAcquisition == TypeAcquisition.UCB);
end

function output = isTypeAcquisitionLCB(this)
    output = (this.typeAcquisition == TypeAcquisition.LCB);
end

function output = getStrTypeAcquisition(this)
    switch this.typeAcquisition
        case TypeAcquisition.EI
            output = 'EI';
        case TypeAcquisition.UCB
            output = 'UCB';
        case TypeAcquisition.LCB
            output = 'LCB';
        otherwise
            output = 'ERROR:unknown';
    end
end
%--------------------- [maxOptIter] ---------------------
function setMaxOptIter(this,maxOptIter)
    this.maxOptIter = maxOptIter;
end

function output = getMaxOptIter(this)
    output = this.maxOptIter;
end

%--------------------- [typeGpFit] ---------------------
function setTypeGpFit(this,typeGpFit)
    this.typeGpFit = typeGpFit;
end

function output = getTypeGpFit(this)
    output = this.typeGpFit;
end


%--------------------- [typeInfillOptimization] ---------------------
function setMethodSearch(this,typeInfillOptimization)
    this.typeInfillOptimization = typeInfillOptimization;
end

function output = getMethodSearch(this)
    output = this.typeInfillOptimization;
end


%--------------------- [typeAcquisition;] ---------------------

function output = getTypeAcquisition(this)
    output = this.typeAcquisition;;
end

    end
end

