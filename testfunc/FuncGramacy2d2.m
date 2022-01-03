classdef FuncGramacy2d2 < AbsFunction
  
    properties
        fnEval;
        dim = 2;
        name = ' x1^2 + x2^2 - 1.5 <=0';
        optVal;
        %
        % -------- Original -----------
        xDomain = [ 0 1; 0 1];
        optSol = [ 0 0  ; 0 0];
        % -------- Adjusted -----------
%         xDomain = repmat([0,1] ,2 ,1);
%         optSol = [ 1, 1];
        % ----------------------------------------------
        power = 1;
        type = TypeFunction.Gramacy2d2;
    end
    
    methods
        function init(this,power)
            this.fnEval = @(x1, x2) ( x1.^2 + x2.^2 - 1.5 );
            this.optVal = this.fnEval(this.optSol(1,1) ,this.optSol(1,2));
        end

    end
end

