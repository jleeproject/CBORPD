classdef InitializerSobol < AbsInitializer
    properties
        type = TypeInitializer.Lhd;
    end
            
    methods
        function [yy_arr, xx_arr, samplesize_arr, yy_cons, y_orig_cell] = getSamples(this)
%             vals_init_standardized = (lhsdesign( this.nInitSample,size(this.xDomain,1)  ));
%             vals_init_standardized = (lhsdesign( this.nInitSample,size(this.xDomain,1)  ));
            vals_init_standardized = net(sobolset(size(this.xDomain,1),'skip',0),this.nInitSample);

            xx_arr = resizeRange(vals_init_standardized, [0, 1], this.xDomain  );
            samplesize = this.samplesizeForVariance;
%             dim_eval = this.objFunc.getDim();
            n_init_sample = this.nInitSample;
            samplesize_arr = ones(n_init_sample, 1)*samplesize;
            
            yy_arr = zeros(n_init_sample,1);
            yy_cons = zeros(n_init_sample,this.nConstraints);
            y_orig_cell = cell(n_init_sample,1);
            for i = 1:n_init_sample
                [y_obj, y_cons, y_orig_cell{i}] = MultipleSamplerObjCon.evaluate(this.simulator, xx_arr(i,:), samplesize);
                yy_arr(i) = y_obj;
                yy_cons(i,:) = y_cons;
            end
        end
    end
end

