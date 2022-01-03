classdef GrapheneModelConstants
    
    properties (Constant)
        q=1.60218*10^(-19);%elementary charge
        hbar=1.05458*10^(-34);%Reduced Plank constant
        epso=8.85418*10^(-12); %Vacuum permittivity
        k=1.38064852*10^(-23);%Boltzmann constant
        Vf=10^6; %Fermi velocity
        ka_top=0;% top relative permittivity
        ka_back=3.9;%back relative permittivity
        Vg0=0; %flatband voltage top
        m=2.511;%square of refence channel potential
        Vg = 0;
        Vds = 1;
        t_top=15*10^(-9);% Top thickness

    end
    
end

