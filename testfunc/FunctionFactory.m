classdef FunctionFactory
    %BOOBJECTFACTORY Summary of this class goes here
    %   Detailed explanation goes here
    
    methods (Static)
        %% AcquisitionFunction Objects
        function obj = getFunction(type, varargin)
            switch type
                case TypeFunction.Branin2
                    obj = FuncBranin2(varargin{:});
                case TypeFunction.Levy6
                    obj = FuncLevy6(varargin{:});
                case TypeFunction.Shc2
                    obj = FuncShc2(varargin{:});
                case TypeFunction.Thc2	
                    obj = FuncThc2(varargin{:});
                case TypeFunction.Norm2
                    obj = FuncNorm2(varargin{:});
                case TypeFunction.Norm4
                    obj = FuncNorm4(varargin{:});
                case TypeFunction.Norm6
                    obj = FuncNorm6(varargin{:});
                case TypeFunction.SumSquare2
                    obj = FuncSumSquare2(varargin{:});
                case TypeFunction.Gramacy1d2
                    obj = FuncGramacy1d2(varargin{:});
                case TypeFunction.Gramacy2d2
                    obj = FuncGramacy2d2(varargin{:});
                case TypeFunction.Sum2
                    obj = FuncSum2(varargin{:});
                case TypeFunction.Sum4
                    obj = FuncSum4(varargin{:});
                case TypeFunction.Gardener14_1C1d2
                    obj = FuncGardner14_1C1d2(varargin{:});
                case TypeFunction.Gardener14_1Objd2
                    obj = FuncGardner14_1Objd2(varargin{:});
                case TypeFunction.Gardener14_2C1d2
                    obj = FuncGardner14_2C1d2(varargin{:});
                case TypeFunction.Gardener14_2Objd2
                    obj = FuncGardner14_2Objd2(varargin{:});
                case TypeFunction.Hartmann4
                    obj = FuncHartmann4(varargin{:});
                case TypeFunction.Hartmann6
                    obj = FuncHartmann6(varargin{:});
%                 case TypeFunction.
%                     obj = Func(varargin{:});
                otherwise
                    error(sprintf('Undefined Class :%s', type.char))
            end
        end
        
    end

end

