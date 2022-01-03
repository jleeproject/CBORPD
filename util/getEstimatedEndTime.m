function finT = getEstimatedEndTime(elapsedT, iter, nIter)
    remainedT = getEstimatedRemainedTime(elapsedT, iter, nIter);
    finT = now + seconds(remainedT);
end