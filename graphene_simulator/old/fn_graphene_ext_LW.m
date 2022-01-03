function [F] = fn_graphene_ext_LW(x, Vb, fixedParams, chngParams, constants ,L,W)
%     global q hbar epso T L W k Vf t_top t_back ka_top ka_back Vg0 m Ct Cb Vg;
%     global q hbar  T k Vf Vg0 m Ct Cb Vg;

    Rd      = fixedParams(1);
    Rs      = fixedParams(2);
    Eg      = fixedParams(3);
    alpha   = fixedParams(4);
    alpha1  = fixedParams(5);
    
    Vb0      = chngParams(1);
    mu_n     = chngParams(2);
    mu_p     = chngParams(3);
    delta    = chngParams(4);

    q     =  constants.q   ;
    hbar  =  constants.hbar;

    T     =  constants.T   ;
%     L     =  constants.L   ;
%     W     =  constants.W   ;
    k     =  constants.k   ;
    Vf    =  constants.Vf  ;
    Vg0   =  constants.Vg0 ;
    m     =  constants.m   ;
    Ct    =  constants.Ct  ;
    Cb    =  constants.Cb  ;
    Vg    =  constants.Vg  ;

    Vds   =  constants.Vds ;
    
    n_pud=delta^2*q^2/(pi*(hbar*Vf)^2);
    a=pi*(k*T)^2/(3*(hbar*Vf)^2);
    b=q^2/(pi*(hbar*Vf)^2);
    c=q^2/(k*T*log(4))^2;
    h=(mu_n+mu_p)/2;
    z=mu_p-mu_n;
    e=0.285;g=0.058;
    f=(q/(5*k*T))^2;
    d=2*q^2*k*T*log(4)/((Ct+Cb)*pi*(hbar*Vf)^2);
%     global Vds


    
    F = zeros(3,1);
    F(1)=200*(x(1)*(Ct+Cb)+1/2*Ttrans(x(1),Eg)*2*q^2*k*T*log(4)/pi/(hbar*Vf)^2*sqrt(1+(q*Ttrans(x(1),Eg)/(k*T*log(4)))^2)+Ct*(Vg-Vg0-x(3)*Rs)+Cb*(Vb-Vb0-x(3)*Rs));

    F(2)=200*(x(2)*(Ct+Cb)+1/2*Ttrans(x(2),Eg)*2*q^2*k*T*log(4)/pi/(hbar*Vf)^2*sqrt(1+(q*Ttrans(x(2),Eg)/(k*T*log(4)))^2)+Ct*(Vg-Vg0-(Vds-x(3)*Rd))+Cb*(Vb-Vb0-(Vds-x(3)*Rd)));

    func2=@(t)q*(a+b*Ttrans(t,Eg).^2).*(h+14*z*Ttrans(t,Eg)./(sqrt(1+c*Ttrans(t,Eg).^2))).*(m./(m+Ttrans(t,Eg).^2)).*(1+(d*alpha*(2+c*Ttrans(t,Eg).^2.*(3+2*c*Ttrans(t,Eg).^2)))./(1+c*Ttrans(t,Eg).^2).^(1.5));

    func3=@(tt) 1/Vf*(1+f*Ttrans(tt,Eg).^2)./(e+(1+f*Ttrans(tt,Eg).^2)*g).*(m./(m+Ttrans(tt,Eg).^2)).*(h+(14*z*(a+b*Ttrans(tt,Eg).^2).*Ttrans(tt,Eg))./((a+b*Ttrans(tt,Eg).^2+n_pud).*sqrt(1+c*Ttrans(tt,Eg).^2))+(h*d*alpha*(2+c*Ttrans(tt,Eg).^2.*(3+2*c*Ttrans(tt,Eg).^2)))./((1+c*Ttrans(tt,Eg).^2).^(1.5))+(d*alpha*(2+c*Ttrans(tt,Eg).^2.*(3+2*c*Ttrans(tt,Eg).^2))*14*z.*Ttrans(tt,Eg).*(a+b*Ttrans(tt,Eg).^2))./((1+c*Ttrans(tt,Eg).^2).^2.*(a+b*Ttrans(tt,Eg).^2+n_pud)));

    F(3)=80*(x(3)-alpha1*(W*(integral(func2,x(1),x(2))+Vds*q*n_pud*h*(m/(m+(x(1)+x(2))^2/4))))/(L+abs(integral(func3,x(1),x(2)))));
    
    clear  Rd Rs Eg alpha alpha1 Vb0 mu_n mu_p delta n_pud a b c h z e f d 
end

