classdef FuncGardner14_1C1d2 < AbsFunction
  
    properties
        fnEval;
        dim = 2;
        name = 'Gardener2014 Prob 1. Const.1';
        optVal;
        %
        % -------- Original -----------
        xDomain = [ 0 0; 0 0];
        optSol = [0 0];
        % -------- Adjusted -----------
%         xDomain = repmat([0,1] ,2 ,1);
%         optSol = [ 1, 1];
        % ----------------------------------------------
        power = 1;
        type = TypeFunction.Gardener14_1C1d2;
    end
    
    methods
        function init(this,power)
            this.fnEval = @(x1, x2) cos(x1).*cos(x2) - sin(x1).*sin(x2) - .5;
            this.optVal = this.fnEval(this.optSol(1,1) ,this.optSol(1,2));
        end

    end
end

