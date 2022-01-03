function [index] = binarySearch_function(func, num, min_range, max_range, tolerance, ascending, maxIter, debug)

    if(nargin ==7)
        debug = false;
    end
    %--------------------------------------------------------------------------
    % Syntax:       [index] = binarySearch(A, n, num);
    %               
    % Inputs:       A: Array (sorted) that you want to search
    %               n: Length of array A
    %               num: Number you want to search in array A
    %               
    % Outputs:      index: Return position in A that A(index) == num
    %                      or -1 if num does not exist in A
    %               
    % Description:  This function find number in array (sorted) using binary
    %               search
    %               
    % Complexity:   O(1)    best-case performance
    %               O(log_2 (n))    worst-case performance
    %               O(1)      auxiliary space

    %--------------------------------------------------------------------------
    % if(  min(size(func)) ~= 1  )
    %     disp('[ ERROR ] input is not a vector')
    %     return;
    % else
    %     dim = max(size(func));
    % end

    if(nargin>6)
        hasMaxIter = true;
        iter = 1;
    else
        hasMaxIter = false;
    end
    left = min_range;
    right = max_range;
    flag = 0;
    % while left <= right
    mid = ((left + right) / 2);
    % mid =@(left,right) ((left + right) / 2);
    yy = func(mid);
    while abs( yy - num ) > tolerance
        mid = ((left + right) / 2);
        yy = func(mid);
    %     if abs( func(mid) - num ) < tolerance
    %         index = mid;
    %         flag = 1;
    %         break;
    %     else if func(mid) > num
    %         right = mid;
    %         else
    %             left = mid;
    %         end
    %     end
        if(ascending)
            if yy < num 
                left = mid;
            else
                right = mid;
            end
        else
            if yy < num 
                right = mid;
            else
                left = mid;
            end
        end
        if(debug)
            fprintf('#=%3d: left = %g\t mid = %g\t right = %g\t Value = %g\t:\tTarget = %g \t (Diff = %g)\n',iter, left,mid,right,yy, num, yy - num);
        end
        iter = iter+1;
        if(hasMaxIter)
            if(iter == maxIter)
    %             return;
                error('reached MAX ITERATION');
            end
        end
    end
    if flag == 0;
    %     index = -1;
        index = mid;
    end

        
end
