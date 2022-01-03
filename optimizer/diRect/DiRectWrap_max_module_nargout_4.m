function [maxVal, maxPt, history, queries, queryVals, queryVals2, queryVals3, queryVals4] = ...
  DiRectWrap_max_module_nargout_4(func,bounds,opts,varargin)
% A wrapper function for diRect.m so that I could feed functions in my format.
% It now treats this as a maximization problem and the inputs to func should be
% row vectors.

    Problem.f = @(t, varargin) transFunctionReturnFirst2Negative_d4(func,t');
    
  [ret_minval, final_xatmin, history, queries, queryVals, queryVals2, queryVals3, queryVals4] = ...
    diRect_nargout_func_4(Problem, bounds, opts, varargin);
  queryVals = -queryVals;
  history(:,3) = -history(:,3);

%   [maxVal, maxIdx] = max(queryVals);
%   maxPt = queries(maxIdx, :);
    maxVal = -ret_minval;
    maxPt = final_xatmin';
    fprintf('%.4f, %.4f, %.4f, %.4f\n',final_xatmin, ret_minval, maxPt, maxVal);
end

