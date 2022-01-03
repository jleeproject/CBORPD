classdef AbsSimulator
    properties (Abstract)
        type
    end

    properties
        param_BO
        objFunc
    end
    methods (Abstract)
        res = evaluate(this, xx, samplesize)
    end
    methods
        function this = AbsSimulator(param_BO, objFunc)
            this.param_BO = param_BO;
            this.objFunc = objFunc;
        end
    end
end

