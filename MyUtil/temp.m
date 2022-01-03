close all;
xxl = repmat(cc,1,nIter-nBurnin+1);
idx = 1;
zz=zeros(M,nIter-nBurnin+1);zz(:,:)=theta(:,idx,nBurnin:end);
figure(1);hold on;
    lb = mean(zz,2)-2*std(zz,0,2);
    ub = mean(zz,2)+2*std(zz,0,2);
   
    plot(xxl,zz,'-');
    plot(xxl,(true_theta(:,idx)),'ok'); 
    plot(xxl,mean(zz,2),'-r','LineWidth',2);
    plot(xxl,ub,'-b','LineWidth',2);
    plot(xxl,lb,'-b','LineWidth',2);

    title('\theta(c) (o:data, -:est)');  ylabel('\theta(.)'); xlabel('c');

hold off;