function out = nthOutputReturnWrapper( fn, idx, varargin)
    nout = nargout(fn);
    if(nout<idx)
        error('Output Index is larger than # output arguments in function.');
    end
    outcell = cell(nout,1);
	[outcell{:}] = fn(varargin{:});
    out = outcell{idx};
end