classdef SimulatorLogSampleVar < AbsSimulator

    properties
        type = TypeSimulator.LogSampleVar;
    end

    methods
        
        function res = evaluate(this, xx, samplesize)
            if(samplesize>1)
    %             if(nargin<4)
                dim_eval = size(xx,2);
    %             end
                args = mat2cell(xx,[1] ,ones(1,dim_eval));
                fnEval = this.objFunc.getFnEval;
                eval_true_sigma = fnEval(args{:});       
            %         eval_true_sigma = fn_eval_true(xx{:});
                sample = normrnd(0, eval_true_sigma, samplesize,1);
                logs = log(std(sample, 0, 'all'));
                res = logs   - get_bias(samplesize);
            %     fprintf('[debug]  [sampled vs true]\tlog(s)=\t%f \t %f\n', logs, log(eval_true_sigma));
                disp('[WARNING] deprecated. need updates. SimulatorStochasticMeanSigma is up to date.');
            else
                res = [];
            end
        end
    end
end

