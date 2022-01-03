function gap = feasibilityGap(constraint, fn_eval_true, xx)
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
%     if(size(xx,2)==2 )
%         if(constraint.hasUb && constraint.hasLb)
%             gap = ...
%                 max([(  lb - fn_eval_true(xx(:,1), xx(:,2))  )  ,...
%                 (  fn_eval_true(xx(:,1), xx(:,2)) - ub  ),0])  ;
%         elseif(constraint.hasUb)
%             gap= max(  [fn_eval_true(xx(:,1), xx(:,2)) -ub,0]  )  ;
%         elseif(constraint.hasLb)
%             gap = max( [ lb - fn_eval_true(xx(:,1), xx(:,2)) ,0]    ) ;
%         else
%             throwError('scr_visualize_BO_online... constraints bounds are not specified.')
%         end
%     else
%         arg = mat2arg(xx);
        if(constraint.hasUb && constraint.hasLb)
            val1 = (  lb - fn_eval_true(xx)  ) ;
            val2 = (  fn_eval_true(xx) - ub  );
            gap= max(  [val1, val2 ,zeros(size(val1)) ], [], 2  )  ;
        elseif(constraint.hasUb)
            val = fn_eval_true(xx) -ub;
            gap= max(  [val,zeros(size(val)) ], [], 2  )  ;
        elseif(constraint.hasLb)
            val = lb - fn_eval_true(xx);
            gap= max(  [val,zeros(size(val)) ], [], 2  )  ;
        else
            throwError('scr_visualize_BO_online... constraints bounds are not specified.')
        end
%     else
%         throwError('Not implemented: more than two dimensions ');
%     end
end