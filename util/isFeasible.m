function [feas, strCon] = isFeasible(constraint, fn_eval_true, xx)
    if(constraint.hasUb && constraint.hasLb)
        ub = constraint.ub;
        lb = constraint.lb;
    elseif(constraint.hasUb)
        ub = constraint.ub;
    elseif(constraint.hasLb)
        lb = constraint.lb;
    else
        throwError('scr_visualize_BO_online... constraints bounds are not specified.')
    end
%     args = mat2arg(xx);
    vals = fn_eval_true(xx);
    if(constraint.hasUb && constraint.hasLb)
        strCon = sprintf('ub=%.2g, lb=%.2g', ub, lb);
        feas = ...
            (  vals >lb  )  .*...
            (  vals <ub  )  ;
    elseif(constraint.hasUb)
        strCon = sprintf('ub=%.2g', ub);
        feas = (  vals <ub  )  ;
    elseif(constraint.hasLb)
        strCon = sprintf('lb=%.2g', lb);
        feas = (  vals >lb   ) ;
    else
        throwError('scr_visualize_BO_online... constraints bounds are not specified.')
    end
end