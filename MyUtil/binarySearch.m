function [index] = binarySearch(vector, num)
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
if(  min(size(vector)) ~= 1  )
    disp('[ ERROR ] input is not a vector')
    return;
else
    dim = max(size(vector));
end
left = 1;
right = dim;
flag = 0;
while left <= right
    mid = ceil((left + right) / 2);
    
    if vector(mid) == num
        index = mid;
        flag = 1;
        break;
    else if vector(mid) > num
        right = mid - 1;
        else
            left = mid + 1;
        end
    end
end
if flag == 0;
%     index = -1;
    index = mid;
end
end
