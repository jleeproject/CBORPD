function [minVal, minPt, history, queries, queryVals, queryVals2, queryVals3, queryVals4] = ...
  DiRectWrap_min_module_nargout_4(func,bounds,opts,varargin)
% A wrapper function for diRect.m so that I could feed functions in my format.
% It now treats this as a maximization problem and the inputs to func should be
% row vectors.

  Problem.f = @(t, varargin) func(t');
  [ret_minval, final_xatmin, history, queries, queryVals, queryVals2, queryVals3, queryVals4] = ...
    diRect_nargout_func_4(Problem, bounds, opts, varargin);
%   queryVals = queryVals;
%   history(:,3) = history(:,3);

%   [minVal, minIdx] = min(queryVals);
%   minPt = queries(minIdx, :);
    minVal = ret_minval;
    minPt = final_xatmin';

end

