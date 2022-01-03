classdef MultipleSamplerObjCon
    
    properties
    end
    
    methods (Static)
        function [yy_arr, yy_cons, yy_obs] = evaluate(simulator, xx, samplesize)
            data = simulator.evaluate(xx, samplesize);
            if(data.hasConstraints)                
                yy_arr = data.getDataForObjFunc();
                yy_cons = data.getDataForConstraints();
            else                
                yy_arr = data.getDataForObjFunc();
                yy_cons = [];
            end
            if(sum(isinf(yy_arr))>0)
                % if sample variance is zero and log(s) = -inf, try one more
                data = simulator.evaluate(xx, samplesize);
                if(data.hasConstraints)                
                    yy_arr = data.getDataForObjFunc();
                    yy_cons = data.getDataForConstraints();
                else                
                    yy_arr = data.getDataForObjFunc();
                    yy_cons = [];
                end
                % if still infeasible,
                if(sum(isinf(yy_arr))>0)
                    yy_arr = log(abs(mean(data.getOriginalData,'all')).*1e-5)-100;
                end
                yy_obs = data.getOriginalData();
            else
                yy_obs = data.getOriginalData();
            end
        end
    end
end

