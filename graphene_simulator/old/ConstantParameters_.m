classdef ConstantParameters
    % Constants used in the graphene model
    
    properties
        q=1.60218*10^(-19);%elementary charge
        hbar=1.05458*10^(-34);%Reduced Plank constant
        epso=8.85418*10^(-12); %Vacuum permittivity
        T=300; %Temperature
        L=2*10^-6;
        W=3*10^-6;
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
        Vds=1;
    end
    
    methods
        function this = ConstantParameters(inputArg1,inputArg2)
            %FIXEDPARAMETERS Construct an instance of this class
            %   Detailed explanation goes here
            this.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(this,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = this.Property1 + inputArg;
        end
    end
end

