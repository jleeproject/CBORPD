function showEstimatedEndTime(elapsedT, iter, nIter)
    remainedT = elapsedT/iter* (nIter-iter) ;
%     fprintf('Expected Remained Time:%s\n', showPrettyDateTime(remainedT));
    fprintf('Expected End Time:%s\n', showPrettyElapsedTime(now + seconds(remainedT)) );
end