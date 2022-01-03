function [pos,val]=argRank(x,rr)
nRows=size(x,1);
nCols=size(x,2);
nVal = nRows*nCols;
zz=reshape(x,nVal,1);
% crt=prctile(zz,percent);
% idx = [1:size(aa,2)]'
% aa=x;
ss=sort(zz,'descend');
crt=ss(rr,1);
    
    if(size(x,1)==1)
        idx = [1:nCols];
        pos2=idx(  x>=crt  )';
        val2=x(pos2)';
        orig= [val2,pos2];
        res=sortrows(orig,-1);
        val=res(:,1);
        pos=res(:,2:end);
        
    elseif(size(x,2)==1)
        idx = [1:nRows]';
        pos2=idx(  x>=crt  );
%         val=x(pos);
        val2=x(pos2)';
        orig= [val2,pos2];
        res=sortrows(orig,-1);
        val=res(:,1);
        pos=res(:,2:end);
    else
        idx2 = ones(nRows,1)*[1:nCols];
        idx1 = [1:nRows]'*ones(1,nCols);

        r1=idx1(  x>=crt  );
        r2=idx2(  x>=crt  );
        pos2=[r1,r2];
        val3=x(r1,r2);
        val2=diag(val3);
        
        orig= [val2,pos2];
        res=sortrows(orig,-1);
        val=res(:,1);
        pos=res(:,2:end);
    end
    
