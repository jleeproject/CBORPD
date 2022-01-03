function [pos,val]=argmax(aa)
    [pos,val]=argRank(aa,1);
%     if(size(aa,1)==1)
%         idx = [1:size(aa,2)]
%         ret=idx(  abs(aa)==max(abs(aa))  )
%     elseif(size(aa,2)==1)
%         idx = [1:size(aa,2)]'
%         ret=idx(  abs(aa)==max(abs(aa))  )
%     else
%         idx1 = [1:size(aa,2)]'
%         ret=idx(  abs(aa)==max(max(abs(aa)))  )
%     end
%     
