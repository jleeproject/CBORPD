% 5 Fixed variables
% 4 Changing Variables
%% Constants
% global Vds q hbar epso T k Vf t_top t_back ka_top ka_back Vg0 m Ct Cb Vg L W;
Vds=1;
q=1.60218*10^(-19);%elementary charge
hbar=1.05458*10^(-34);%Reduced Plank constant
epso=8.85418*10^(-12); %Vacuum permittivity
T=300; %Temperature
L=2*10^-6;
W=3*10^-6;

L_orig = L;
W_orig = W;

k=1.38064852*10^(-23);%Boltzmann constant
Vf=10^6; %Fermi velocity
t_top=15*10^(-9);% Top thickness
t_back=300*10^(-9); % Back thickness
ka_top=0;% top relative permittivity
ka_back=3.9;%back relative permittivity
Vg0=0; %flatband voltage top
m=2.511;%square of refence channel potential
Ct=epso*ka_top/t_top;
Cb=epso*ka_back/t_back;
Vg=0;

constants.q     =  q   ;
constants.hbar  =  hbar;

constants.T     =  T   ;
constants.L     =  L   ;
constants.W     =  W   ;
constants.k     =  k   ;
constants.Vf    =  Vf  ;
constants.Vg0   =  Vg0 ;
constants.m     =  m   ;
constants.Ct    =  Ct  ;
constants.Cb    =  Cb  ;
constants.Vg    =  Vg  ;

constants.Vds   =  Vds;
%% AIR
%%% Fixed Parameters
% 1)
Rd=3450;% Parasitic drain (Fixed parameter 1)
Rs=Rd;% Parasitic source (Fixed parameter 2, this parameter is typically set the same as the parameter 1)
Eg=0.1;  %(Fixed parameter 3)
alpha=0.7; % (Fixed parameter 4)
alpha1=0.079; %(Fixed parameter 5)

%%% Changing Parameters
% 1) c=0 : Water
Vb0=46; %flatband voltage back (changing parameter 1)
mu_n=0.001143;%low field mobiliy (changing parameter 2) : Curvature /
mu_p=0.00229;%low field mobiliy (changing parameter 3) : Curvature \
delta=1.06050515523;% Spatial inhomogeneity of potential (changing parameter 4)

fixedParams = getSettingFixedParam(Rd, Rs, Eg, alpha, alpha1);
chngParams_4 = getSettingChangingParam(Vb0, mu_n, mu_p, delta);

weight_mu = 1.05;
weight_delta = 0.975;
%% c = 0 : Pure Water
%%% Fixed Parameters
% % 1)
% Rd=3450;% Parasitic drain (Fixed parameter 1)
% Rs=Rd;% Parasitic source (Fixed parameter 2, this parameter is typically set the same as the parameter 1)
% Eg=0.1;  %(Fixed parameter 3)
% alpha=0.7; % (Fixed parameter 4)
% alpha1=0.079; %(Fixed parameter 5)

%%% Changing Parameters
% 1) c=0 : Water
Vb0=25; %flatband voltage back
mu_n=0.0011583;%low field mobiliy
mu_p=0.0031;%low field mobiliy
delta=0.8930515523;% Spatial inhomogeneity of potential
delta = delta * weight_delta*.97;
mu_n =  mu_n * weight_mu * 1.5;
mu_p = mu_p * weight_mu*0.9;

fixedParams = getSettingFixedParam(Rd, Rs, Eg, alpha, alpha1);
chngParams_1 = getSettingChangingParam(Vb0, mu_n, mu_p, delta);


%% c=0 : Contaminated
%%% Fixed Parameters
% % 1)
% Rd=3450;% Parasitic drain (Fixed parameter 1)
% Rs=Rd;% Parasitic source (Fixed parameter 2, this parameter is typically set the same as the parameter 1)
% Eg=0.1;  %(Fixed parameter 3)
% alpha=0.7; % (Fixed parameter 4)
% alpha1=0.079; %(Fixed parameter 5)

%%% Changing Parameters
% 1) c=10 : Water
Vb0=19; %flatband voltage back
mu_n=0.00233;%low field mobiliy
mu_p=0.0025;%low field mobiliy
delta=0.827515523;% Spatial inhomogeneity of potential
% delta = delta * weight_delta;
mu_n =  mu_n * weight_mu;
mu_p = mu_p * weight_mu;

% delta =     0.8374;
delta =     0.8344;
mu_n =  0.00189;
mu_p = 0.00285;
% fixedParams_1 = getSettingFixedParam(Rd, Rs, Eg, alpha, alpha1);
chngParams_2 = getSettingChangingParam(Vb0, mu_n, mu_p, delta);

%% 3) c=10 : Contaminated
%%% Fixed Parameters
% % 1)
% Rd=3450;% Parasitic drain (Fixed parameter 1)
% Rs=Rd;% Parasitic source (Fixed parameter 2, this parameter is typically set the same as the parameter 1)
% Eg=0.1;  %(Fixed parameter 3)
% alpha=0.7; % (Fixed parameter 4)
% alpha1=0.079; %(Fixed parameter 5)

%%% Changing Parameters
% 1) c=20 : Water
Vb0=9; %flatband voltage back
mu_n=0.00181;%low field mobiliy
mu_p=0.0027;%low field mobiliy
delta=0.851523;% Spatial inhomogeneity of potential
delta = delta * weight_delta;
mu_n =  mu_n * weight_mu;
mu_p = mu_p * weight_mu;

% fixedParams_1 = getSettingFixedParam(Rd, Rs, Eg, alpha, alpha1);
chngParams_3 = getSettingChangingParam(Vb0, mu_n, mu_p, delta);

% %% Set Variables
% global Vb Vds;
% n_pud=delta^2*q^2/(pi*(hbar*Vf)^2);
% a=pi*(k*T)^2/(3*(hbar*Vf)^2);
% b=q^2/(pi*(hbar*Vf)^2);
% c=q^2/(k*T*log(4))^2;
% h=(mu_n+mu_p)/2;
% z=mu_p-mu_n;
% e=0.285;g=0.058;
% f=(q/(5*k*T))^2;
% d=2*q^2*k*T*log(4)/((Ct+Cb)*pi*(hbar*Vf)^2);

zz = zeros(4,3);
    zz(1,:) = [   17	10	3	];
    zz(2,:) = [   0.0024583	0.00203	0.00162	]; % /
    zz(3,:) = [   0.0036	0.0031	0.0026	]; % \
    zz(4,:) = [   0.747052	0.803516	0.8808	];
    

    chngParams_1 = zz(:,1)';
    chngParams_2 = zz(:,2)';
    chngParams_3 = zz(:,3)';
    chngParams_all = {chngParams_1, chngParams_2, chngParams_3};
