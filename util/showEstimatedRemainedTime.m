function showEstimatedRemainedTime(elapsedT, iter, nIter)
    remainedT = elapsedT/iter* (nIter-iter) ;
    fprintf('[%s] Expected Remained Time:%s\n', showPrettyDateTime(now), showPrettyElapsedTime(remainedT));
%     fprintf('Expected End Time:%s\n', showPrettyElapsedTime(now + seconds(remainedT)) );
end