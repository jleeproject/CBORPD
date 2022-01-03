classdef CombinedHSPPrior <handle
    %CombinedHSPPrior Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        label='cHS+';
    end
    
    methods
        function this = CombinedHSPPrior()
        end
                
        function res = paramLabels(this)
            res = sprintf('(%s)', this.method);
        end
    end
end