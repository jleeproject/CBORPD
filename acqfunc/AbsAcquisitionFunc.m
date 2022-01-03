classdef AbsAcquisitionFunc < handle  & matlab.mixin.Copyable
    
    properties (Abstract)
        name
        type
        isMaximizeAcq
    end
    
    properties
        isMaximizeObj
        param_BO
        objFunc
    end
    
    methods (Abstract)
        validate(this)
        acquire(this, predictor, x_pred, iter)
        init(this)
    end
    
    methods
        function this = AbsAcquisitionFunc(param_BO, objFunc, isMaximizeObj)
            this.param_BO = param_BO;
            this.isMaximizeObj = isMaximizeObj;
%             this.isMaximizeAcq = isMaximizeAcq;
            this.objFunc = objFunc;
            validate(this);
            init(this)
        end
        
        function updateFOpt(this, f_opt)
        end
        
    end
end

