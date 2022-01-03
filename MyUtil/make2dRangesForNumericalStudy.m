function [xx,yy, ranges, ranges2] = make2dRangesForNumericalStudy(sttNum,endNum, numBin, sttNum2,endNum2, numBin2)
    scr_makeLoggerIfNoExist;
% make2dRangesForNumericalStudy(sttNum,endNum)
%   sttNum: Starting number
%   endNum: Last number
%   numBin: Number of Bins
% 
% Example:
% [xrange,yrange] = make2dRangesForNumericalStudy(1,3)
% 
% xrange =
% 
%      1     2     3
%      1     2     3
%      1     2     3
% 
% yrange =
% 
%      1     1     1
%      2     2     2
%      3     3     3
%     xx=ones(endNum-sttNum+1,1)*(sttNum:endNum);
%     yy=xx';
% 
% Example 2)
% [xrange,yrange] = make2dRangesForNumericalStudy(1,2,4)
% 
% xrange =
% 
%     1.0000    1.3333    1.6667    2.0000
%     1.0000    1.3333    1.6667    2.0000
%     1.0000    1.3333    1.6667    2.0000
%     1.0000    1.3333    1.6667    2.0000
% 
% 
% yrange =
% 
%     1.0000    1.0000    1.0000    1.0000
%     1.3333    1.3333    1.3333    1.3333
%     1.6667    1.6667    1.6667    1.6667
%     2.0000    2.0000    2.0000    2.0000

    if(nargin==2)
        numBin = endNum - sttNum + 1;
    end

    if(nargin == 3)
        numBin= round(numBin);
        ranges = [0:(numBin-1)]/(numBin-1)*(endNum-sttNum)+sttNum;

        xx=ones(numBin,1)*ranges;
        yy=xx';
    %     r_row =[1:nsamples]./nsamples*max_num_r;
    elseif(nargin == 6)
        numBin= round(numBin);
        numBin2 = round(numBin2);
        % xdim = numBin;
        % ydim = numBin2;
        ranges = [0:(numBin-1)]/(numBin-1)*(endNum-sttNum)+sttNum;
        ranges2 = [0:(numBin2-1)]/(numBin2-1)*(endNum2-sttNum2)+sttNum2;

        xx=ones(numBin2,1)*ranges;
        yy=ranges2'*ones(1,numBin);
%         yy=xx';
    %     r_row =[1:nsamples]./nsamples*max_num_r;
    else
        logger.error(  sprintf('WRONG NUMBER OF ARGUMENTS : %d', nargin)  );
    end

end

