function remainedT = getEstimatedRemainedTime(elapsedT, iter, nIter)
    remainedT = elapsedT/iter* (nIter-iter) ;
end