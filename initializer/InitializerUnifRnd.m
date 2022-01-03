classdef InitializerUnifRnd < AbsInitializer
    properties
        type = TypeInitializer.UnifRnd;
    end
    
    methods
        function [yy_arr, xx_arr, samplesize_arr] = getSamples(this)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
%             [yy_arr, xx_arr, samplesize_arr] = fn_init_get_unifrnd_values(param_BO);
            fn_eval_true = this.objFunc.getFnEval();
            fn_x_domain = this.objFunc.getFnXDomain();
            dim_eval = this.objFunc.getFnDimEval();
            samplesize = this.param_BO.getSamplesizeForInitStage();
            n_init_sample = this.param_BO.getNInitSample();



            xx_arr = zeros(n_init_sample, dim_eval);
            yy_arr = zeros(n_init_sample, 1);
%             true_sig = zeros(n_init_sample, 1);
            % points to evaluate
            for i=1:dim_eval
                xx_arr(:,i) = unifrnd(fn_x_domain(i,1), fn_x_domain(i,2), n_init_sample,1);
            end

            samplesize_arr = ones(n_init_sample, 1)*samplesize;
            for i = 1:n_init_sample
                logs = eval_simul(fn_eval_true, xx_arr(i,:), samplesize, dim_eval);        
                yy_arr(i) = logs;
            end
        end
    end
end

