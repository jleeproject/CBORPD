function args = mat2arg(xs)
    args = mat2cell(xs,[size(xs,1)],ones(1,size(xs,2)));
end