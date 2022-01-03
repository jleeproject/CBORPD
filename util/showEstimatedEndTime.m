function showEstimatedEndTime(elapsedT, iter, nIter)
    remainedT = elapsedT/iter* (nIter-iter) ;
%     fprintf('Expected Remained Time:%s\n', showPrettyDateTime(remainedT));
    fprintf('[%s] Expected End Time:%s\n', showPrettyDateTime(now), showPrettyElapsedTime(now + seconds(remainedT)) );
end