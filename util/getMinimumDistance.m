function diff_x_f_min = fn_get_minimum_distance(xmin_history_sel, opt_sol_fn)
% opt_sol_fn [row]: multiple solutions : 
% xmin_history [row]: history : 
% column: dimensions
dim = size(opt_sol_fn,1);
diff_x_f_mins = zeros(size(xmin_history_sel,1),dim);
    for i = 1:dim
        diff_x_f_mins(:,i) = vecnorm( xmin_history_sel - opt_sol_fn(i,:), 2, 2);
    end
    diff_x_f_min = min(diff_x_f_mins,[],2);
end