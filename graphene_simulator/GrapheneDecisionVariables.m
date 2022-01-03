classdef GrapheneDecisionVariables < handle
    % Constants used in the graphene model
    
%     properties (Access = private)
%         %% Variables
%         T=300; %Temperature
%         L=2*10^-6; 
%         W=3*10^-6; % Micrometer 
%         t_top=15*10^(-9);% Top thickness
%         t_back=300*10^(-9); % Back thickness
%         %% 
% %         Ct=GrapheneModelConstants.epso*GrapheneModelConstants.ka_top/t_top;
% %         Cb=GrapheneModelConstants.epso*GrapheneModelConstants.ka_back/t_back;
% %         Ct = [];
% %         Cb = [];
%         Vg=0; %% Voltage at Gate:
%         Vds=1; %% Voltage at Drain-source
%         Vb = [];
%     end
    
    properties 
        T=300; %Temperature
        L=2*10^-6; 
        W=3*10^-6; % Micrometer 
        t_top=15*10^(-9);% Top thickness
        t_back=300*10^(-9); % Back thickness
        %% 
%         Ct=GrapheneModelConstants.epso*GrapheneModelConstants.ka_top/t_top;
%         Cb=GrapheneModelConstants.epso*GrapheneModelConstants.ka_back/t_back;
%         Ct = [];
%         Cb = [];
        Vg=0; %% Voltage at Gate:
        Vds=1; %% Voltage at Drain-source
        Vb = [];
        domains
        nVar  = 5;
        names = {'T', 'L', 'W', 't_back', 'Vb'};
    end
    
    methods
        function this = GrapheneDecisionVariables(T, L, W,  t_back, Vb)
            if(nargin>7)
                setParameters(this, T, L, W, t_back, Vb);
            elseif(nargin>5)
                setParameters(this, T, L, W, t_back);
%                 this = setVoltages(this, Vg, Vds);
            elseif(nargin>0)
                setParameters(this, T, L, W, t_back);
            end
            init_domains(this);
        end
        
        function init_domains(this)
            this.domains.T = [223.15, 373.15]; % -50 C ~ 100 C
%             this.domains.L = [0.02 2].*1e-6; % : Impact of Channel Length and Width for Charge Transportation of Graphene Field Effect Transistor
            this.domains.L = [0.2 3].*1e-6; % : Impact of Channel Length and Width for Charge Transportation of Graphene Field Effect Transistor
            this.domains.W = [1 6].*1e-6;  %  2~5 : Impact of Channel Length and Width for Charge Transportation of Graphene Field Effect Transistor
%             this.domains.t_top = [10 20].* 1e-9;
            this.domains.t_back = [10, 400].* 1e-9; %% https://www.sciencedirect.com/science/article/pii/S0375960121000037
%             this.domains.Vg = [-3 3];
%             this.domains.Vg = [-80 80];
%             this.domains.Vg = [-120 120]; % https://www.researchgate.net/publication/277332975_Functionalized_Graphene_and_Other_Two-Dimensional_Materials_for_Photovoltaic_Devices_Device_Design_and_Processing/figures?lo=1&utm_source=google&utm_medium=organic
%             this.domains.Vds = [-80 80];
%             this.domains.Vds = [-1 1];
            this.domains.Vb = [-80 80];
            this.domains.all = [this.domains.T;this.domains.L;this.domains.W; this.domains.t_back; this.domains.Vb];
%             this.domains.all = [this.domains.T;this.domains.L;this.domains.W;this.domains.t_top;this.domains.t_back;this.domains.Vg;this.domains.Vds;this.domains.Vb];
        end
        
        function out = getValues(this)
            out = [this.T; this.L; this.W; this.t_back; this.Vb ];
        end
        
%         function this = setVoltages(this, Vg, Vds)
%             %% Input Arguments: Vg, Vds
%             this.Vg  = Vg;
%             this.Vds = Vds;
%         end

        function this = setParameters(this, T, L, W, t_back, Vb)
            this.T = T;
            this.L = L;
            this.W = W;
%             this.t_top = t_top;
            this.t_back = t_back;
%             calcCtCb(this, t_top, t_back);
            
%             if(nargin>6)
%                 setVoltages(this, Vg, Vds);
%             end
            if (nargin>6)
                this.Vb =Vb;
            end
        end
        
%         function this = calcCtCb(this, t_top, t_back)
%             this.Ct=GrapheneModelConstants.epso*GrapheneModelConstants.ka_top/t_top;
%             this.Cb=GrapheneModelConstants.epso*GrapheneModelConstants.ka_back/t_back;
%         end
        
        function out = hasVb(this)
            out = (numel(this.Vb)>0);
        end
    end
end

