function showEstimatedEndRemainedTime(elapsedT, iter, nIter)
%     remainedT = elapsedT/iter* (nIter-iter) ;
    remainedT = getEstimatedRemainedTime(elapsedT, iter, nIter);
    finT = getEstimatedEndTime(elapsedT, iter, nIter);
    fprintf('[%s] Expected End Time: %s (Remained:%s)\n', showPrettyDateTime(now), showPrettyDateTime(finT), showPrettyElapsedTime(remainedT));
%     fprintf('Expected End Time:%s\n', showPrettyElapsedTime(now + seconds(remainedT)) );
end