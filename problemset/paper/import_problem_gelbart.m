problem_name = ...
    'Gelbart';


isMaximizeObj = false; power = 1;
cell_obj_func = {...
    FunctionFactory.getFunction(TypeFunction.Branin2, isMaximizeObj,power) ...
    };

cell_type_mean_func = {...
    {{TypeFunction.SumSquare2, FuncModifier(@(y) 15*(y-50), {@(x) x-2.5, @(x) x-7.5} )}}...
    };
cell_const_func = {...
    MultipleConstraints(...
        {Constraint(TypeConstraint.UnknownMeanConstraint, [], 0 ) }...
    )...
};
