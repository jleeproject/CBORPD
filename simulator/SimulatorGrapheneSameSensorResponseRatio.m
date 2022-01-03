classdef SimulatorGrapheneSameSensorResponseRatio  < matlab.mixin.Copyable
% classdef SimulatorStochasticMeanSigma < AbsSimulator

%% Designed for Conventiona BO problem.
    properties
        type;
%         meanFunc
%         sigmaFunc

        grapheneSolver
        
        
        nCon
        
        obj
        cell_con_meanFunc
        cell_con_sigmaFunc
        hasMean
        multiConsObj
        typeProblem
        
        hasOnOffChangingParameters = false;
        changingParameterSensorOn % pure water
        changingParameterSensorOff % contaminated water
    end

    methods
        function this = SimulatorGrapheneSameSensorResponseRatio(grapheneSolver, optProb, typeProblem, changingParameterSensorOn, changingParameterSensorOff)

            if(optProb.nCon>0)
                this.multiConsObj = optProb.con;
            end

            this.grapheneSolver = grapheneSolver;
            this.typeProblem = typeProblem;
            
            if(nargin>3)
                this.hasOnOffChangingParameters = true;
                this.changingParameterSensorOn = changingParameterSensorOn;
                this.changingParameterSensorOff = changingParameterSensorOff;
            else
                if( ~grapheneSolver.hasSecondChangingParam )
                    error('[SimulatorGrapheneSameSensorOnOffRatio] Changing parameter of contaminated sensor is needed either in solver or in this object.');
                end
            end
        end
        
        function out = evaluate(this, xx, samplesize)
            out = [];
            if(contains(this.typeProblem.char,'Conventional'))
                samplesize = 1;
            end
            if(samplesize>0)

                dim_eval = size(xx,2);
                nSamples = size(xx, 1);

                if(this.typeProblem==TypeProblem.RobustDesignOptimization && samplesize == 1)
                    error('[ERROR : SimulatorStochasticMeanSigma] For Robust Design Optimization, need >1 samplesize');
                end
                
                if(this.typeProblem==TypeProblem.ConventionalBoStochastic)
                    error('[SimulatorStochasticGrapheneSensor] Not defined in Conventional BO model');
                elseif(this.typeProblem==TypeProblem.ConventionalBoDeterministic)
                    error('[SimulatorStochasticGrapheneSensor] Not defined in Conventional BO model');
                else
                    if(nSamples > 1)
                        error('[SimulatorStochasticGrapheneSensor] No rules are defined when given x has multiple values');
                    end
                    samples = zeros(samplesize, 1);
                    for i = 1:samplesize
%                         samples(i) = this.grapheneSolver.evalWithVecX(xx);
%                         samples(i) = this.grapheneSolver.evalOnOffRatioWithVecX(xx, this.typeSensorOn, this.typeSensorOff);

%                         samples(i) = this.grapheneSolver.evalResponseRatioWithVecX(xx, this.typeSensorOn, this.typeSensorOff);
                        if(this.hasOnOffChangingParameters)
                            samples(i) = this.grapheneSolver.evalSameSensorResponseRatioWithVecX(xx, this.changingParameterSensorOn, this.changingParameterSensorOff);
                        else
                            samples(i) = this.grapheneSolver.evalSameSensorResponseRatioWithVecX(xx);
                        end
                    end
                    out = DataTransformer(this.multiConsObj, samples, samplesize, this.typeProblem);
                end
%                 res = sample;
%                 logs = log(std(sample, 0, 'all'));
            end
        end
        
%         function sample = getNormalSample(this, meanFunc, sigmaFunc, xx_args, samplesize)
%             if(contains(class(this.obj.sigmaFunc),'TypeFunction'))
%                 sFunc = FunctionFactory.getFunction(this.obj.sigmaFunc,[]).getFnEval;
%                 eval_true_sigma = sFunc(xx_args{:});
%             elseif( isa(sigmaFunc,'AbsFunction') )
%                 sFunc = sigmaFunc.getFnEval;
%                 eval_true_sigma = sFunc(xx_args{:});
%             elseif( isnumeric(sigmaFunc) )
%                 eval_true_sigma = sigmaFunc;
%             else
%                 error('[ERROR] Undefined type');
%             end
% %             if(this.hasMean)
%                 %%
%                 if(contains(class(meanFunc),'TypeFunction'))
%                     mFunc = FunctionFactory.getFunction(meanFunc,[]).getFnEval;
%                     eval_true_mean = mFunc(xx_args{:});    
%                 elseif( isa(meanFunc,'AbsFunction') )
%                     mFunc = meanFunc.getFnEval;
%                     eval_true_mean = mFunc(xx_args{:});    
%                 elseif( isnumeric(meanFunc) )
%                     eval_true_mean = meanFunc;
%                 else
%                     error('[ERROR] Undefined type');
%                 end
%                 if(samplesize==1)
%                     sample = normrnd(eval_true_mean, eval_true_sigma);
%                 elseif(numel(eval_true_mean)==1)
%                     sample = normrnd(eval_true_mean, eval_true_sigma, samplesize,1);
%                 else
%                     if(numel(eval_true_mean) ~= samplesize)
%                         error('Dimension of evaluated mean and samplesize does not match.');
%                     else
%                         sample = normrnd(eval_true_mean, eval_true_sigma, samplesize,1);
%                     end
%                 end
% %             else
% %                 sample = normrnd(0, eval_true_sigma, samplesize,1);
% %             end
%         end
        
%         function sample = getDeterministicSample(this, meanFunc, xx_args)
%                 if(contains(class(meanFunc),'TypeFunction'))
%                     mFunc = FunctionFactory.getFunction(meanFunc,[]).getFnEval;
%                     eval_true_mean = mFunc(xx_args{:});    
%                 elseif( isa(meanFunc,'AbsFunction') )
%                     mFunc = meanFunc.getFnEval;
%                     eval_true_mean = mFunc(xx_args{:});    
%                 elseif( isnumeric(meanFunc) )
%                     eval_true_mean = meanFunc;
%                 else
%                     error('[ERROR] Undefined type');
%                 end
%                 sample = eval_true_mean;
%         end
        
    end
end

