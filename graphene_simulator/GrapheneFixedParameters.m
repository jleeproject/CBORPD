classdef GrapheneFixedParameters < handle
    %FIXEDPARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Rd=3450;% Parasitic drain (Fixed parameter 1)
        Rs=3450;% Parasitic source (Fixed parameter 2, this parameter is typically set the same as the parameter 1)
%         Rs=Rd;% Parasitic source (Fixed parameter 2, this parameter is typically set the same as the parameter 1)
        Eg=0.1;  %(Fixed parameter 3)
        alpha=0.7; % (Fixed parameter 4)
        alpha1=0.079; %(Fixed parameter 5)
    end
    
    methods
        function this = GrapheneFixedParameters(Rd, Eg, alpha, alpha1)
            if(nargin>0)
                this.Rd = Rd;
                this.Rs = Rd;
                this.Eg = Eg;
                this.alpha = alpha;
                this.alpha1 = alpha1;
            end
        end
    end
end

