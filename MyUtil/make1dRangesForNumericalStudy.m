function [ranges] = make1dRangesForNumericalStudy(sttNum,endNum, numBin)
% makeRangesForNumericalStudy(sttNum,endNum)
%   sttNum: Starting number
%   endNum: Last number
%   numBin: Number of Bins
% 
% Example 1)
% ranges = make1dRangesForNumericalStudy(1,4)
% 
% ranges =
% 
%      1     2     3     4
% 
% Example 2)
% make1dRangesForNumericalStudy(1,2,4)
% 
% ans =
% 
%     1.0000    1.3333    1.6667    2.0000
    numBin = round(numBin);
    if(nargin==2)
        numBin = endNum - sttNum + 1;
    end

    ranges = [0:(numBin-1)]/(numBin-1)*(endNum-sttNum)+sttNum;
    
end

