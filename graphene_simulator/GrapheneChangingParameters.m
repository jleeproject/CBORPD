classdef GrapheneChangingParameters < handle
    properties
        Vb0; %flatband voltage back (changing parameter 1)
        mu_n;%low field mobiliy (changing parameter 2) : Curvature /
        mu_p;%low field mobiliy (changing parameter 3) : Curvature \
        delta;% Spatial inhomogeneity of potential (changing parameter 4)
% Vb0, mu_n, mu_p, delta
    end
    
    methods
        function this = GrapheneChangingParameters(typeSensingSolution)
            switch typeSensingSolution
                case TypeGrapheneSensingSolution.AIR
                    setParametersInAir(this);
                case TypeGrapheneSensingSolution.WATER_C0
                    setParametersInWaterC0(this);
                case TypeGrapheneSensingSolution.WATER_C10
                    setParametersInWaterC10(this);
                case TypeGrapheneSensingSolution.WATER_C20
                    setParametersInWaterC20(this);
                otherwise
                    error('[ChangingParameters] Undefined TypeSensingSolution');
            end
        end
        
        function setParametersInAir(this)
            this.Vb0=46; %flatband voltage back (changing parameter 1)
            this.mu_n=0.001143;%low field mobiliy (changing parameter 2) : Curvature /
            this.mu_p=0.00229;%low field mobiliy (changing parameter 3) : Curvature \
            this.delta=1.06050515523;% Spatial inhomogeneity of potential (changing parameter 4)
        end
        
        function setParametersInWaterC0(this)
%             this.Vb0=21; %flatband voltage back
%             this.mu_n=0.0011583;%low field mobiliy
%             this.mu_p=0.0031;%low field mobiliy
%             this.delta=0.8930515523;% Spatial inhomogeneity of potential
            weight_mu = 1.05;
            weight_delta = 0.975;
            
            this.Vb0=25; %flatband voltage back
            this.mu_n=0.0011583  * weight_mu * 1.5;%low field mobiliy
            this.mu_p=0.0031  * weight_mu*0.9;%low field mobiliy
            this.delta=0.8930515523  * weight_delta*.97;% Spatial inhomogeneity of potential
%             delta = delta * weight_delta*.97;
%             mu_n =  mu_n * weight_mu * 1.5;
%             mu_p = mu_p * weight_mu*0.9;
            
        end
        function setParametersInWaterC10(this)
%             this.Vb0=20; %flatband voltage back
%             this.mu_n=0.00233;%low field mobiliy
%             this.mu_p=0.0025;%low field mobiliy
%             this.delta=0.827515523;% Spatial inhomogeneity of potential
            this.Vb0=19; %flatband voltage back
%             mu_n=0.00233;%low field mobiliy
%             mu_p=0.0025;%low field mobiliy
%             delta=0.827515523;% Spatial inhomogeneity of potential
            % delta = delta * weight_delta;
%             mu_n =  mu_n * weight_mu;
%             mu_p = mu_p * weight_mu;

            % delta =     0.8374;
            this.delta =     0.8344;
            this.mu_n =  0.00189;
            this.mu_p = 0.00285;
        end
        function setParametersInWaterC20(this)
%             this.Vb0=46; %flatband voltage back (changing parameter 1)
%             this.mu_n=0.00181;%low field mobiliy
%             this.mu_p=0.0027;%low field mobiliy
%             this.delta=0.851523;% Spatial inhomogeneity of potential
            weight_mu = 1.05;
            weight_delta = 0.975;

            this.Vb0=9; %flatband voltage back
            this.mu_n=0.00181 * weight_mu;%low field mobiliy
            this.mu_p=0.0027 * weight_mu;%low field mobiliy
            this.delta=0.851523  * weight_delta;% Spatial inhomogeneity of potential
%             delta = delta * weight_delta;
%             mu_n =  mu_n * weight_mu;
%             mu_p = mu_p * weight_mu;

        end
    end
end

