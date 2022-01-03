classdef FuncSum4 < AbsFunction
  
    properties
        fnEval;
        dim = 4;
        name = 'Sum';
        optVal;
        %
        % -------- Original -----------
        xDomain = [ 0 1; 0 1; 0 1; 0 1];
        optSol = [ 0 0  0 0 ];
        % -------- Adjusted -----------
%         xDomain = repmat([0,1] ,2 ,1);
%         optSol = [ 1, 1];
        % ----------------------------------------------
        power = 1;
        type = TypeFunction.Branin2;
    end
    
    methods
        function init(this,power)
            this.fnEval = @(x1, x2, x3, x4) ( x1 + x2 + x3 + x4 );
            this.optVal = this.fnEval(this.optSol(1,1) ,this.optSol(1,2), this.optSol(1,3) ,this.optSol(1,4));
        end

    end
end

