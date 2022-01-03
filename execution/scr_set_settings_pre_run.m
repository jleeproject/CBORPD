isMaximizeObj = setting.isMaximizeObj;
power = setting.power;
nRepeats = setting.nRepeats;

if isfield(setting,'type_budget')
    switch setting.type_budget
        case TypeBudget.FixedBudget
            budget = setting.budget;
        case TypeBudget.MultiplyDimBy
            budget = setting.budget *size(setting.optProb.xDomain,1);
        otherwise
            throwError('Unknown Type of type_budget');
    end
else
    budget = setting.budget;
end
param_BO.budget = budget;

isDebug = setting.isDebug;
isInfo = setting.isInfo;
show_plot_online = setting.show_plot_online;
show_default = setting.show_default;
show_exact_fmin = setting.show_exact_fmin;
max_opt_iter = setting.max_opt_iter;
type_gp_fit = setting.type_gp_fit;
type_infill_opt = setting.type_infill_opt;
type_acquisition = setting.type_acquisition;
type_initial_sampler = setting.type_initial_sampler;
type_initial_sample_size = setting.type_initial_sample_size;
type_simulator = setting.type_simulator;
type_gp_kernel = setting.type_gp_kernel;
is_conservative = setting.is_conservative;
wantDetailedOptimalValues = setting.wantDetailedOptimalValues;
show_all_grid_visualization = setting.show_all_grid_visualization;
show_all_search_trails = setting.show_all_search_trails;
grid_bin = setting.grid_bin;




