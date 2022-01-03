function out = lognrndFromMeanStdDev(M, sig, varargin);
    V = sig^2;
    MU = log(M^2 / sqrt(V+M^2));
    SIGMA = sqrt(log(V/M^2 + 1));
    out = lognrnd(MU,SIGMA,varargin{:});
end
