problem_name = ...
    'Hartmann6';

isMaximizeObj = false; power = 1;
cell_obj_func = {...
    FunctionFactory.getFunction(TypeFunction.Hartmann6, isMaximizeObj,power, FuncModifier(@(y) y+3,[],[],[0.0000    0.0000    0.0000    0.5556 0 0] )) ...
    };

cell_type_mean_func = {...
    {{TypeFunction.Norm6, FuncModifier(@(y) (y+40)*80, {@(x) x, @(x) x} )}}...
    };
cell_const_func = {...
    MultipleConstraints(...
        {Constraint(TypeConstraint.UnknownMeanConstraint, [], 3280 ) }...
    )...
};

sigmaObj = 0.2;
sigmaCons = 0.2;
