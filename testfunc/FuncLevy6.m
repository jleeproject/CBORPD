classdef FuncLevy6 < AbsFunction
    properties (Constant)
    end
    
    properties
        fnEval;
%         rng_levy_all= [-10,10];
        xDomain =  repmat([-10,10] ,6,1);
        dim = 6;
        name = 'levy';
        optSol = [ 1, 1, 1, 1, 1, 1];
        optVal;
        %
        % rng_branin_all = [0,1];
        % rng_branin = repmat(rng_branin_all ,dim_branin,1);
        %
        power = 1;
        type = TypeFunction.Levy6;
    end
    
    methods
        function init(this,power)
            this.fnEval = @(x1, x2, x3, x4, x5, x6) ( sin(pi.* (1+(x1-1)/4) ).^2 + (  1+(x1-1)/4  -1).^2 .* (1+ 10 .* sin(pi.* (1+(x1-1)/4) + 1).^2 ) + (  1+(x2-1)/4  -1).^2 .* (1+ 10 .* sin(pi.* (1+(x2-1)/4) + 1).^2 )+ (  1+(x3-1)/4  -1).^2 .* (1+ 10 .* sin(pi.* (1+(x3-1)/4) + 1).^2 )+ (  1+(x4-1)/4  -1).^2 .* (1+ 10 .* sin(pi.* (1+(x4-1)/4) + 1).^2 )+ (  1+(x5-1)/4  -1).^2 .* (1+ 10 .* sin(pi.* (1+(x5-1)/4) + 1).^2 )+ (  1+(x6-1)/4  -1).^2 .* (1+ sin(2.* pi.* (1+(x6-1)/4) ).^2 ) + 10  ).^power;
            this.optVal = this.fnEval(this.optSol(1,1) ,this.optSol(1,2) ,this.optSol(1,3) ,this.optSol(1,4) ,this.optSol(1,5) ,this.optSol(1,6));
        end
    end
end

