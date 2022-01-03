classdef SimulatorSimpleMeanError < AbsSimulator

    properties
        type = TypeSimulator.LogSampleVar;
        sigerr = 1e-3;
    end

    methods
        
        function res = evaluate(this, xx)
            if(samplesize>1)
    %             if(nargin<4)
                dim_eval = size(xx,2);
    %             end
                args = mat2cell(xx,[1] ,ones(1,dim_eval));
                fnEval = this.objFunc.getFnEval;
                eval_true_mu = fnEval(args{:});       
            %         eval_true_sigma = fn_eval_true(xx{:});
                sample = normrnd(eval_true_mu, sigerr);
                res = sample;
            %     fprintf('[debug]  [sampled vs true]\tlog(s)=\t%f \t %f\n', logs, log(eval_true_sigma));
            else
                res = [];
            end
        end
    end
end

