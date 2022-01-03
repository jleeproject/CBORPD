classdef FuncShc2 < AbsFunction
    properties (Constant)
    end
    
    properties
        fnEval;
        xDomain = [-3 3;-2,2];
        dim = 2;
        name = 'SHC';
        optSol = [.0898, -0.7126;-0.0898, 0.7126];
        optVal;
        %
        % rng_branin_all = [0,1];
        % rng_branin = repmat(rng_branin_all ,dim_branin,1);
        %
        power = 1;
        type = TypeFunction.Shc2;
    end
    
    methods
        function init(this,power)
            this.fnEval = @(x1,x2)  (   (4-2.1 .*x1.^2 + x1.^4/3).*x1.^2 +x1.*x2 + (-4+4.*x2.^2).*x2.^2 + 10 ).^power;
            this.optVal = this.fnEval(this.optSol(1,1) ,this.optSol(1,2));
        end
    end
end

