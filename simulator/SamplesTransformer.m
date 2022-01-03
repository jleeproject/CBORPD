classdef SamplesTransformer
    
    properties
        samples
        samplesize
    end
    
    methods
        function this = SamplesTransformer(samples, samplesize)
            this.samples = samples;
            this.samplesize = samplesize;
        end
        
        function out = getUnbiasSampleLogVar(this)
            if(this.samplesize>0)
                stdval = std(this.samples, 0, 'all');
%                 if(stdval ==0)
%                     logs = log(mean(this.samples,'all').*1e-100);
%                 else
                    logs = log(stdval);
%                 end
                out = logs - get_bias_logv(this.samplesize);
            else
                out = [];
            end
        end
        
        function out = getSampleMean(this)
            if(this.samplesize>0)
                out = mean(this.samples);
            else
                out = [];
            end
        end
        
        function out = getSamples(this)
            out = this.samples;
        end
    end
end

