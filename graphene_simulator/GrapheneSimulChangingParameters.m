classdef GrapheneSimulChangingParameters < handle
    properties
        Vb0_; %flatband voltage back (changing parameter 1)
        mu_n_;%low field mobiliy (changing parameter 2) : Curvature /
        mu_p_;%low field mobiliy (changing parameter 3) : Curvature \
        delta_;% Spatial inhomogeneity of potential (changing parameter 4)
% Vb0_, mu_n_, mu_p_, delta_
        sig_Vb0;
        sig_mun;
        sig_mup;
        sig_delta;
    end
    
    methods
        function this = GrapheneSimulChangingParameters(typeSensingSolution)
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
            init(this)
        end
        
        function init(this)
            this.sig_Vb0 = 10;
            this.sig_mun = this.mu_n_/200;
            this.sig_mup = this.mu_p_/200;
            this.sig_delta = this.delta_/500;
%             this.sig_mun = 0;
%             this.sig_mup = 0;
%             this.sig_delta = 0;
        end
        
        function setStdDevs(this, sig_Vb0, sig_mun, sig_mup, sig_delta)
            this.sig_Vb0 = sig_Vb0;
            this.sig_mun = sig_mun;
            this.sig_mup = sig_mup;
            this.sig_delta = sig_delta;
        end
        
        function setParametersInAir(this)
            this.Vb0_=46; %flatband voltage back (changing parameter 1)
            this.mu_n_=0.001143;%low field mobiliy (changing parameter 2) : Curvature /
            this.mu_p_=0.00229;%low field mobiliy (changing parameter 3) : Curvature \
            this.delta_=1.06050515523;% Spatial inhomogeneity of potential (changing parameter 4)
        end
        
        function setParametersInWaterC0(this)
%             this.Vb0_=21; %flatband voltage back
%             this.mu_n_=0.0011583;%low field mobiliy
%             this.mu_p_=0.0031;%low field mobiliy
%             this.delta_=0.8930515523;% Spatial inhomogeneity of potential
            weight_mu = 1.05;
            weight_delta_ = 0.975;
            
            this.Vb0_=25; %flatband voltage back
            this.mu_n_=0.0011583  * weight_mu * 1.5;%low field mobiliy
            this.mu_p_=0.0031  * weight_mu*0.9;%low field mobiliy
            this.delta_=0.8930515523  * weight_delta_*.97;% Spatial inhomogeneity of potential
%             delta_ = delta_ * weight_delta_*.97;
%             mu_n_ =  mu_n_ * weight_mu * 1.5;
%             mu_p_ = mu_p_ * weight_mu*0.9;
            
        end
        function setParametersInWaterC10(this)
%             this.Vb0_=20; %flatband voltage back
%             this.mu_n_=0.00233;%low field mobiliy
%             this.mu_p_=0.0025;%low field mobiliy
%             this.delta_=0.827515523;% Spatial inhomogeneity of potential
            this.Vb0_=19; %flatband voltage back
%             mu_n_=0.00233;%low field mobiliy
%             mu_p_=0.0025;%low field mobiliy
%             delta_=0.827515523;% Spatial inhomogeneity of potential
            % delta_ = delta_ * weight_delta_;
%             mu_n_ =  mu_n_ * weight_mu;
%             mu_p_ = mu_p_ * weight_mu;

            % delta_ =     0.8374;
            this.delta_ =     0.8344;
            this.mu_n_ =  0.00189;
            this.mu_p_ = 0.00285;
        end
        function setParametersInWaterC20(this)
%             this.Vb0_=46; %flatband voltage back (changing parameter 1)
%             this.mu_n_=0.00181;%low field mobiliy
%             this.mu_p_=0.0027;%low field mobiliy
%             this.delta_=0.851523;% Spatial inhomogeneity of potential
            weight_mu = 1.05;
            weight_delta_ = 0.975;

            this.Vb0_=9; %flatband voltage back
            this.mu_n_=0.00181 * weight_mu;%low field mobiliy
            this.mu_p_=0.0027 * weight_mu;%low field mobiliy
            this.delta_=0.851523  * weight_delta_;% Spatial inhomogeneity of potential
%             delta_ = delta_ * weight_delta_;
%             mu_n_ =  mu_n_ * weight_mu;
%             mu_p_ = mu_p_ * weight_mu;

        end
        
        function out = Vb0(this)
            out = this.Vb0_ + normrnd(0, this.sig_Vb0);
        end
        
        function out = mu_n(this)
            out = this.mu_n_ + normrnd(0, this.sig_mun);
        end
        
        function out = mu_p(this)
            out = this.mu_p_ + normrnd(0, this.sig_mup);
        end
        
        function out = delta(this)
            out = this.delta_ + normrnd(0, this.sig_delta);
        end
        
    end
end

