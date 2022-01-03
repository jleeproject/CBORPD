function [out1, out2, out3] = transFunctionReturnFirst2Negative_d3_module(func,xx)
    [out1, out2, out3] = func(xx);
    out1 = -out1;
end
