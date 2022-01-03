classdef SimulatorStochasticMeanSigma  < matlab.mixin.Copyable
% classdef SimulatorStochasticMeanSigma < AbsSimulator

%% Designed for Conventiona BO problem.
    properties
        type;
%         meanFunc
%         sigmaFunc
        nCon
        
        obj
        cell_con_meanFunc
        cell_con_sigmaFunc
        hasMean
        multiConsObj
        typeProblem
    end

    methods
        function this = SimulatorStochasticMeanSigma(simul, optProb, typeProblem)

            if(optProb.nCon>0)
                this.multiConsObj = optProb.con;
            end
            this.obj.meanFunc = simul.obj.meanFunc;
            this.obj.sigmaFunc = simul.obj.sigmaFunc;
            
            this.nCon = simul.nCon;
            this.cell_con_meanFunc = cell(simul.nCon,1);
            this.cell_con_sigmaFunc = cell(simul.nCon,1);
            for i=1:simul.nCon
                this.cell_con_meanFunc{i} = simul.con{1}.meanFunc;
                this.cell_con_sigmaFunc{i} = simul.con{1}.sigmaFunc;
            end
            
            if(numel(this.obj.meanFunc)>0)
                this.hasMean = true;
            else
                this.hasMean = false;
                this.obj.meanFunc = @(x) 0;
            end
            this.typeProblem = typeProblem;
        end
        
        function out = evaluate(this, xx, samplesize)
            out = [];
            if(contains(this.typeProblem.char,'Conventional'))
                samplesize = 1;
            end
            if(samplesize>0)
    %             if(nargin<4)
                dim_eval = size(xx,2);
                nSamples = size(xx, 1);
    %             end
%                 args = mat2cell(xx,[1] ,ones(1,dim_eval));
                args = mat2arg(xx);
                %
%                 %%
%                 if(contains(class(this.obj.sigmaFunc),'TypeFunction'))
%                     sFunc = FunctionFactory.getFunction(this.obj.sigmaFunc,[]).getFnEval;
%                     eval_true_sigma = sFunc(args{:});
%                 elseif( isa(this.obj.sigmaFunc,'AbsFunction') )
%                     sFunc = this.obj.sigmaFunc.getFnEval;
%                     eval_true_sigma = sFunc(args{:});
%                 elseif( isnumeric(this.obj.sigmaFunc) )
%                     eval_true_sigma = this.obj.sigmaFunc;
%                 else
%                     error('[ERROR] Undefined type');
%                 end
%                 if(this.hasMean)
%                     %%
%                     if(contains(class(this.obj.meanFunc),'TypeFunction'))
%                         mFunc = FunctionFactory.getFunction(this.obj.meanFunc,[]).getFnEval;
%                         eval_true_mean = mFunc(args{:});    
%                     elseif( isa(this.obj.meanFunc,'AbsFunction') )
%                         mFunc = this.obj.meanFunc.getFnEval;
%                         eval_true_mean = mFunc(args{:});    
%                     elseif( isnumeric(this.obj.meanFunc) )
%                         eval_true_mean = this.obj.meanFunc;
%                     else
%                         error('[ERROR] Undefined type');
%                     end
%                     sampleObj = normrnd(eval_true_mean, eval_true_sigma, samplesize,1);
%                 else
%                     sampleObj = normrnd(0, eval_true_sigma, samplesize,1);
%                 end
                if(this.typeProblem==TypeProblem.RobustDesignOptimization && samplesize == 1)
                    error('[ERROR : SimulatorStochasticMeanSigma] For Robust Design Optimization, need >1 samplesize');
                end
                
                if(this.typeProblem==TypeProblem.ConventionalBoStochastic)
                    if(this.hasMean)
                        sampleObj = getNormalSample(this, this.obj.meanFunc, this.obj.sigmaFunc, args, samplesize);
                    else
                        sampleObj = getNormalSample(this, 0, this.obj.sigmaFunc, args, samplesize);
                    end
                    sampleCon = zeros(nSamples, this.nCon);
                    for i = 1:this.nCon
                        sampleCon(:,i) = getNormalSample(this, this.cell_con_meanFunc{i}, this.cell_con_sigmaFunc{i}, args, samplesize);
                    end
                    out = DataTransformer(this.multiConsObj, sampleObj, samplesize, this.typeProblem, sampleCon);
                elseif(this.typeProblem==TypeProblem.ConventionalBoDeterministic)
                    sampleObj = getDeterministicSample(this, this.obj.meanFunc, args);
                    sampleCon = zeros(nSamples, this.nCon);
                    for i = 1:this.nCon
                        sampleCon(:,i) = getDeterministicSample(this, this.cell_con_meanFunc{i}, args);
                    end
                    out = DataTransformer(this.multiConsObj, sampleObj, samplesize, this.typeProblem, sampleCon);
                else
                    if(this.hasMean)
                        sampleObj = getNormalSample(this, this.obj.meanFunc, this.obj.sigmaFunc, args, samplesize);
                    else
                        sampleObj = getNormalSample(this, 0, this.obj.sigmaFunc, args, samplesize);
                    end
                    out = DataTransformer(this.multiConsObj, sampleObj, samplesize, this.typeProblem);
                end
%                 res = sample;
%                 logs = log(std(sample, 0, 'all'));
            end
        end
        
        function sample = getNormalSample(this, meanFunc, sigmaFunc, xx_args, samplesize)
            if(contains(class(this.obj.sigmaFunc),'TypeFunction'))
                sFunc = FunctionFactory.getFunction(this.obj.sigmaFunc,[]).getFnEval;
                eval_true_sigma = sFunc(xx_args{:});
            elseif( isa(sigmaFunc,'AbsFunction') )
                sFunc = sigmaFunc.getFnEval;
                eval_true_sigma = sFunc(xx_args{:});
            elseif( isnumeric(sigmaFunc) )
                    if(numel(sigmaFunc)>0)
                        eval_true_sigma = sigmaFunc;
                    else
                        error('[ERROR] Empty Sigma Function');
                    end
            else
                error('[ERROR] Undefined type');
            end
%             if(this.hasMean)
                %%
                if(contains(class(meanFunc),'TypeFunction'))
                    mFunc = FunctionFactory.getFunction(meanFunc,[]).getFnEval;
                    eval_true_mean = mFunc(xx_args{:});    
                elseif( isa(meanFunc,'AbsFunction') )
                    mFunc = meanFunc.getFnEval;
                    eval_true_mean = mFunc(xx_args{:});    
                elseif( isnumeric(meanFunc) )
                    if(numel(meanFunc)>0)
                        eval_true_mean = meanFunc;
                    else
                        error('[ERROR] Empty Mean Function');
                    end
                else
                    error('[ERROR] Undefined type');
                end
                if(samplesize==1)
                    sample = normrnd(eval_true_mean, eval_true_sigma);
                elseif(numel(eval_true_mean)==1)
                    sample = normrnd(eval_true_mean, eval_true_sigma, samplesize,1);
                else
                    if(numel(eval_true_mean) ~= samplesize)
                        error('Dimension of evaluated mean and samplesize does not match.');
                    else
                        sample = normrnd(eval_true_mean, eval_true_sigma, samplesize,1);
                    end
                end
%             else
%                 sample = normrnd(0, eval_true_sigma, samplesize,1);
%             end
        end
        
        function sample = getDeterministicSample(this, meanFunc, xx_args)
                if(contains(class(meanFunc),'TypeFunction'))
                    mFunc = FunctionFactory.getFunction(meanFunc,[]).getFnEval;
                    eval_true_mean = mFunc(xx_args{:});    
                elseif( isa(meanFunc,'AbsFunction') )
                    mFunc = meanFunc.getFnEval;
                    eval_true_mean = mFunc(xx_args{:});    
                elseif( isnumeric(meanFunc) )
                    eval_true_mean = meanFunc;
                else
                    error('[ERROR] Undefined type');
                end
                sample = eval_true_mean;
        end
        
    end
end

