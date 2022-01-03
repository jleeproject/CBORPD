function [minVal, minPt, history, queries, queryVals, queryVals2, queryVals3, queryVals4] = ...
  DiRectWrap_min_constrained_nargout_3(Problem, func,bounds,opts)
% A wrapper function for diRect.m so that I could feed functions in my format.
% It now treats this as a maximization problem and the inputs to func should be
% row vectors.

  Problem.f = @(t, varargin) func(t');
  [ret_minval, final_xatmin, history, queries, queryVals, queryVals2, queryVals3] = ...
    diRect_nargout_func_3_module(Problem, bounds, opts);

%   [minVal, minIdx] = min(queryVals);
%   minPt = queries(minIdx, :);
    minVal = ret_minval;
    minPt = final_xatmin';

end

