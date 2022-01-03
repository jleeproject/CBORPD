function [mm,ss] = fnGetParam4Dist(mu, variance, priorDist, verbose)
%FNGETPARAM4DIST Summary of this function goes here
%   Detailed explanation goes here
% outputArg1 = inputArg1;
% outputArg2 = inputArg2;

    if(priorDist == Prior.LogNormal)
        if(exist('verbose','var')); if(verbose); disp('');end; end;
        mm = log(mu) - 1/2*log( variance./mu.^2 + 1 );
%         mm = log(mu^2/sqrt(variance + mu^2));
        vv = log(variance./mu.^2 + 1);
        ss = sqrt(vv);
        % TEST
%         samples = lognrnd(mm,ss, [1000000, 1]); figure(1);clf;
%         histogram(samples,'Normalization','pdf','BinLimits',[0,30],'BinWidth',0.5);
%         disp(mean(samples));disp(var(samples))
    else
        disp('[ERROR] No distribution')
    end
        

end

