function res = resizeRange(value ,from_range, to_range)
    if(size(value,2)~= size(to_range,1))
        disp('ERROR: Inconsistant dimensions');
        return;
    end
    res = zeros(size(value));
    for irow = 1:size(to_range,1)
   
        from_lower = from_range(1);
        from_upper = from_range(2);
        to_lower = to_range(irow,1);
        to_upper = to_range(irow,2);

    %     slope = U-L;
    %     intercept =  L;
        slope_intercept = ([from_lower 1; from_upper 1])^(-1)* [to_lower to_upper]';
    %     slope = slope_intercept(1);
    %     intercept = slope_intercept(2);
        res(:,irow) = slope_intercept(1) * value(:,irow)  + slope_intercept(2);
    end
end