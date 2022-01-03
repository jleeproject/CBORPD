if(~exist('multipleProblems','var'))
    cell_obj_func_comb        = cell(1,1);
    cell_type_mean_func_comb  = cell(1,1);
    cell_const_func_comb      = cell(1,1);
    n_problems = 0;
    cell_prob_names = cell(1,1);
    cell_sigmaObj = cell(1,1);
    cell_sigmaCons = cell(1,1);
end

if(~exist('problem_name','var'))
    problem_name = '';
end

if(setting.typeProblem ~= TypeProblem.ConventionalBoDeterministic && setting.typeProblem ~= TypeProblem.RobustDesignOptimization)
    if(exist('sigmaCons','var'))
        prev_sigmaCons = sigmaCons;
    else
        if(exist('prev_sigmaCons','var'))
            fprintf('[WARNING] sigmaCons is not defined. Define previous one:%g\n', prev_sigmaCons);
            sigmaCons = prev_sigmaCons;
        else
            if(setting.typeProblem == TypeProblem.ConventionalBoStochastic)
                error('sigmaCons has never been defined.');
            else
                sigmaCons = 0;
            end
        end
    end
    if(exist('sigmaObj','var'))
        prev_sigmaObj = sigmaObj;
    else
        if(exist('prev_sigmaObj','var'))
            fprintf('[WARNING] sigmaObj is not defined. Define previous one:%g\n', prev_sigmaObj);
            sigmaObj = prev_sigmaObj;
        else
            if(setting.typeProblem == TypeProblem.ConventionalBoStochastic)
                error('sigmaObj has never been defined.');
            else
                sigmaObj = 0;
            end
        end
    end
else
    sigmaObj = 0;
    sigmaCons = 0;
end

multipleProblems = true;
n_problems = n_problems + 1;

if(numel(cell_obj_func)>1)
    for i = 1:numel(cell_obj_func)
        cell_obj_func_comb(n_problems,1) = cell_obj_func(i);
        if(numel(cell_type_mean_func)>1)
            cell_type_mean_func_comb(n_problems,1) = cell_type_mean_func(i);
        else
            cell_type_mean_func_comb(n_problems,1) = cell_type_mean_func(1);
        end
        if(numel(cell_type_mean_func)>1)
            cell_const_func_comb(n_problems,1) = cell_const_func(i);
        else
            cell_const_func_comb(n_problems,1) = cell_const_func(1);
        end
        cell_prob_names{n_problems,1} = problem_name;
        cell_sigmaObj{n_problems,1} = sigmaObj;
        cell_sigmaCons{n_problems,1} = sigmaCons;
%         n_problems = n_problems + 1;
    end
else
    cell_obj_func_comb(n_problems,1) = cell_obj_func;
    cell_type_mean_func_comb(n_problems,1) = cell_type_mean_func;
    cell_const_func_comb(n_problems,1) = cell_const_func;
    cell_prob_names{n_problems,1} = problem_name; 
    cell_sigmaObj{n_problems,1} = sigmaObj;
    cell_sigmaCons{n_problems,1} = sigmaCons;
%     n_problems = n_problems + 1;
end
clear problem_name sigmaCons;

cell_obj_func = cell_obj_func_comb;
cell_type_mean_func = cell_type_mean_func_comb;
cell_const_func = cell_const_func_comb;
