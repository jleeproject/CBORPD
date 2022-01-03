classdef SimulatorDeterministic   < matlab.mixin.Copyable

    properties
        type = TypeSimulator.SimulatorDeterministic;
        nCon
        
        obj
        cell_con_meanFunc
        cell_con_sigmaFunc
        multiConsObj
        typeProblem
    end

    methods
        function this = SimulatorDeterministic(simul, optProb, typeProblem)

            if(optProb.nCon>0)
                this.multiConsObj = optProb.con;
            end
            this.obj.meanFunc = simul.obj.meanFunc;
            this.obj.sigmaFunc = simul.obj.sigmaFunc;
            
            this.nCon = simul.nCon;
            this.cell_con_meanFunc = cell(simul.nCon,1);
            this.cell_con_sigmaFunc = cell(simul.nCon,1);
            for i=1:simul.nCon
                this.cell_con_meanFunc{i} = simul.con{i}.meanFunc;
                this.cell_con_sigmaFunc{i} = simul.con{i}.sigmaFunc;
            end
            
            this.typeProblem = typeProblem;
        end
        
        function out = evaluate(this, xx, samplesize)
            dim_eval = size(xx,2);
            args = mat2cell(xx,[1] ,ones(1,dim_eval));
            sampleObj = getDeterministicSample(this, this.obj.meanFunc, args);
            sampleCon = zeros(1, this.nCon);
            for i = 1:this.nCon
                sampleCon(1,i) = getDeterministicSample(this, this.cell_con_meanFunc{i}, args);
            end
            out = DataTransformer(this.multiConsObj, sampleObj, samplesize, this.typeProblem, sampleCon);
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

