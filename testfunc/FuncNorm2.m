classdef FuncNorm2 < AbsFunction
    properties (Constant)
    end
    
    properties
        fnEval = @(x1,x2) sqrt(x1.^2 + x2.^2);
        xDomain = repmat([-5, 5] ,2,1);
        dim = 2;
        name = 'Norm2';
        optSol = [ 0 0];
        optVal = 0;
        %
        % rng_branin_all = [0,1];
        % rng_branin = repmat(rng_branin_all ,dim_branin,1);
        %
        power = 1;
        type = TypeFunction.Norm2;
    end
    
    methods
        function init(this,power)
%             this.fnEval = @(x1,x2) (   2.*x1.^2 - 1.05 .* x1.^4 + x1.^6 / 6 + x1.*x2 + x2.^2  + 10).^power;
            this.optVal = this.fnEval(this.optSol(1,1) ,this.optSol(1,2));
        end
    end
end

