function [handler, lb, ub] = showCIPlot(x_range, mu_range, sig, coeff, strLegend, color, alpha)
    if(nargin<6)
        color = [128 193 219]./255;
    end
    if(nargin<7)
        alpha = 0.3;
    end
    range_ci = [x_range' fliplr(x_range')]';    
    lb = mu_range'-coeff*sqrt(sig)';
    ub = mu_range'+coeff*sqrt(sig)';
    ci_pred = [lb fliplr(ub)]';
    handler = fill(range_ci, ci_pred, color , 'LineWidth',1, 'FaceAlpha',alpha, 'LineStyle','none');
    if(nargin>3)
        handler.DisplayName = strLegend;
    end;
end        