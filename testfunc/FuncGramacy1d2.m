classdef FuncGramacy1d2 < AbsFunction
  
    properties
        fnEval;
        dim = 2;
        name = 'Gramacy C1 <=0';
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
        type = TypeFunction.Gramacy1d2;
    end
    
    methods
        function init(this,power)
            this.fnEval = @(x1, x2) ( 1.5 - x1 - 2.*x2 - .5 .* sin( 2.*pi .*(x1.^2 - 2.* x2)) );
            this.optVal = this.fnEval(this.optSol(1,1) ,this.optSol(1,2));
        end

    end
end

