classdef FuncSum2 < AbsFunction
  
    properties
        fnEval;
        dim = 2;
        name = 'Sum';
        optVal;
        %
        % -------- Original -----------
        xDomain = [ 0 1; 0 1];
        optSol = [ 0.1927    0.4074 ];
        % -------- Adjusted -----------
%         xDomain = repmat([0,1] ,2 ,1);
%         optSol = [ 1, 1];
        % ----------------------------------------------
        power = 1;
        type = TypeFunction.Branin2;
    end
    
    methods
        function init(this,power)
            this.fnEval = @(x1, x2) ( x1 + x2 );
            this.optVal = this.fnEval(this.optSol(1,1) ,this.optSol(1,2));
        end

    end
end

