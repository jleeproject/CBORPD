problem_name = ...
    'unconstrained 4 func';

isMaximizeObj = false; power = 1;
cell_obj_func = {...
    FunctionFactory.getFunction(TypeFunction.Branin2, isMaximizeObj,power) ...
    ,...
    FunctionFactory.getFunction(TypeFunction.Thc2, isMaximizeObj,power) ...
    , ...
    FunctionFactory.getFunction(TypeFunction.Shc2, isMaximizeObj,power) ...
    , ...
    FunctionFactory.getFunction(TypeFunction.Levy6, isMaximizeObj,power)...
    };

%% Make sure unconstrained.
cell_type_mean_func = {...
    {}...
};
cell_const_func = {...
    {}...
};