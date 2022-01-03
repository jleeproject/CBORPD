function [maxVal, maxPt, history, queries, queryVals] = ...
  diRectWrapMax_module(func,bounds,opts,varargin)
% A wrapper function for diRect.m so that I could feed functions in my format.
% It now treats this as a maximization problem and the inputs to func should be
% row vectors.

  Problem.f = @(t, varargin) -func(t');
  [ret_minval, final_xatmin, history, queries, queryVals] = ...
    diRectOpt(Problem, bounds, opts, varargin);
  queryVals = -queryVals;
  history(:,3) = -history(:,3);
 
%   [maxVal, maxIdx] = max(queryVals);
%   maxPt = queries(maxIdx, :);
    maxVal = - ret_minval;
    maxPt = final_xatmin';

end

