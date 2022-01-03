classdef ConsFuncBranin2 < AbsObjectiveFunction
  
    properties
        fnEval = @(x1, x2) (   (x2-5.1/(4.*pi.^2) .* x1.^2 + 5/pi.*x1 -6 ).^2 + 10.*(1-1/(8.*pi)).*cos(x1)+10   );
        dim = 2;
        name = 'branin';
%         optVal;
        %
        % -------- Original -----------
%         xDomain = [ -5 10; 0 15];
%         optSol = [ -pi, 12.275; pi, 2.275; 9.42478, 2.475];
        % -------- Adjusted -----------
        xDomain = repmat([0,1] ,2 ,1);
%         optSol = [ 1, 1];
        % ----------------------------------------------
%         power = 1;
        type = TypeConstraintFuncion.Branin2;
    end
    
    methods
%         function init(this)
% %             this.fnEval = @(x1, x2) (   (x2-5.1/(4.*pi.^2) .* x1.^2 + 5/pi.*x1 -6 ).^2 + 10.*(1-1/(8.*pi)).*cos(x1)+10   );
% %             this.optVal = this.fnEval(this.optSol(1,1) ,this.optSol(1,2));
%         end

    end
end

