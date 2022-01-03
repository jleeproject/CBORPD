classdef FuncHartmann4 < AbsFunction
  
    properties
        fnEval;
        dim = 4;
        name = 'Hartmann (d=4)';
        optVal;
        %
        % -------- Original -----------
        xDomain = [ 0 1; 0 1;0 1 ; 0 1];
%         optSol = [ 0.1873, 0.1906, 0.5566, 0.2647];
        optSol = [0.3897    0.5000    0.6111    0.8328]; % unconstrained
        % -------- Adjusted -----------
%         xDomain = repmat([0,1] ,2 ,1);
%         optSol = [ 1, 1];
        % ----------------------------------------------
        power = 1;
        type = TypeFunction.Hartmann4;
    end
    
    methods
        function init(this, power)
%             this.fnEval = @(x1, x2, x3, x4) ( x1 + x2 + x3 + x4 );
            this.fnEval = @(x1, x2, x3, x4) eval(this,x1,x2,x3,x4);
            this.optVal = this.fnEval(this.optSol(1,1) ,this.optSol(1,2), this.optSol(1,3) ,this.optSol(1,4));
%             this.optVal = this.fnEval(mat2arg(this.optSol));
        end
        
        function y = eval(this, x1, x2, x3, x4)
            xx = [x1, x2, x3, x4];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % HARTMANN 4-DIMENSIONAL FUNCTION
            %
            % Authors: Sonja Surjanovic, Simon Fraser University
            %          Derek Bingham, Simon Fraser University
            % Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
            %
            % Copyright 2013. Derek Bingham, Simon Fraser University.
            %
            % THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
            % FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
            % derivative works, such modified software should be clearly marked.
            % Additionally, this program is free software; you can redistribute it 
            % and/or modify it under the terms of the GNU General Public License as 
            % published by the Free Software Foundation; version 2.0 of the License. 
            % Accordingly, this program is distributed in the hope that it will be 
            % useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
            % of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
            % General Public License for more details.
            %
            % For function details and reference information, see:
            % http://www.sfu.ca/~ssurjano/
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % INPUT:
            %
            % xx = [x1, x2, x3, x4]
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            alpha = [1.0, 1.2, 3.0, 3.2]';
            A = [10, 3, 17, 3.5, 1.7, 8;
                 0.05, 10, 17, 0.1, 8, 14;
                 3, 3.5, 1.7, 10, 17, 8;
                 17, 8, 0.05, 10, 0.1, 14];
            P = 10^(-4) * [1312, 1696, 5569, 124, 8283, 5886;
                           2329, 4135, 8307, 3736, 1004, 9991;
                           2348, 1451, 3522, 2883, 3047, 6650;
                           4047, 8828, 8732, 5743, 1091, 381];

            outer = 0;
            for ii = 1:4
                inner = 0;
                for jj = 1:4
                    xj = xx(:,jj);
                    Aij = A(ii, jj);
                    Pij = P(ii, jj);
                    inner = inner + Aij*(xj-Pij).^2;
                end
                new = alpha(ii) * exp(-inner);
                outer = outer + new;
            end

            y = (1.1 - outer) / 0.839;

               

        end

    end
end

