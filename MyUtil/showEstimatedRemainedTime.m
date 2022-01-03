function showEstimatedRemainedTime(elapsedT, iter, nIter)
    remainedT = elapsedT/iter* (nIter-iter) ;
    fprintf('Expected Remained Time:%s\n', showPrettyElapsedTime(remainedT));
%     fprintf('Expected End Time:%s\n', showPrettyElapsedTime(now + seconds(remainedT)) );
end