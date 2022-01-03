function [maxVal, maxPt, history, queries, queryVals, queryVals2, queryVals3] = ...
  DiRectWrap_max_constrained_nargout_3(Problem, func,bounds,opts)
% A wrapper function for diRect.m so that I could feed functions in my format.
% It now treats this as a maximization problem and the inputs to func should be
% row vectors.

    Problem.f = @(t, varargin) transFunctionReturnFirst2Negative_d3_module(func,t');
    
  [ret_minval, final_xatmin, history, queries, queryVals, queryVals2, queryVals3] = ...
    diRect_nargout_func_3_module(Problem, bounds, opts);
  queryVals = -queryVals;
  history(:,3) = -history(:,3);

%   [maxVal, maxIdx] = max(queryVals);
%   maxPt = queries(maxIdx, :);
    maxVal = -ret_minval;
    maxPt = final_xatmin';
%     fprintf('%.4f, %.4f, %.4f, %.4f\n',final_xatmin, ret_minval, maxPt, maxVal);
end

