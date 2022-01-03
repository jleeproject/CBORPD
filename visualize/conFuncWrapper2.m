function y_cons = conFuncWrapper2(simulator, xnew, samplesize)
fprintf('A');
    out = zero(size(xnew,1),1);
    for i=1:size(xnew,1)
        xnew(i,:)
        [~, y_cons, ~] = MultipleSamplerObjCon.evaluate(simulator, xnew(i,:), samplesize);
        out(i) = y_cons;
    end
end