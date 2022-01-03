classdef GrapheneSimulFixedParameters < handle
    %FIXEDPARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Rd_ =3450;% Parasitic drain (Fixed parameter 1)
        Rs_ =3450;% Parasitic source (Fixed parameter 2, this parameter is typically set the same as the parameter 1)
%         Rs=Rd;% Parasitic source (Fixed parameter 2, this parameter is typically set the same as the parameter 1)
        Eg_ =0.1;  %(Fixed parameter 3)
        alpha_ =0.7; % (Fixed parameter 4)
        alpha1_ =0.079; %(Fixed parameter 5)
        
        sig_Rd
        sig_Rs
        sig_Eg
        sig_alpha
        sig_alpha1
    end
    
    methods
        function this = GrapheneSimulFixedParameters(Rd, Eg, alpha, alpha1)
            if(nargin>0)
                this.Rd_ = Rd;
                this.Rs_ = Rd;
                this.Eg_ = Eg;
                this.alpha_ = alpha;
                this.alpha1_ = alpha1;
            end
            init(this)
        end
        function init(this)
%             this.sig_Rd = this.Rd_/500;
%             this.sig_Rs = this.Rs_/500;
%             this.sig_Eg = this.Eg_/500;
%             this.sig_alpha = this.alpha_/500;
%             this.sig_alpha1 = this.alpha1_/500;
            this.sig_Rd = 0;
            this.sig_Rs = 0;
            this.sig_Eg = 0;
            this.sig_alpha = 0;
            this.sig_alpha1 = 0;
        end
        
        function setStdDevs(this, sig_Rd, sig_Rs, sig_Eg, sig_alpha, sig_alpha1)
            this.sig_Rd = sig_Rd;
            this.sig_Rs = sig_Rs;
            this.sig_Eg = sig_Eg;
            this.sig_alpha = sig_alpha;
            this.sig_alpha1 = sig_alpha1;
        end
        
        function out = Rd(this)
            out = lognrndFromMeanStdDev(this.Rd_, this.sig_Rd);
%             out = this.Rd_ + normrnd(0, this.sig_Rd);
        end
        
        function out = Rs(this)
            out = lognrndFromMeanStdDev(this.Rs_, this.sig_Rs);
%             out = this.Rs_ + normrnd(0, this.sig_Rs);
        end
        
        function out = Eg(this)
            out = lognrndFromMeanStdDev(this.Eg_, this.sig_Eg);
%             out = this.Eg_ + normrnd(0, this.sig_Eg);
        end
                
        function out = alpha(this)
            out = lognrndFromMeanStdDev(this.alpha_, this.sig_alpha);
%             out = this.alpha_ + normrnd(0, this.sig_alpha);
        end
        
        function out = alpha1(this)
            out = lognrndFromMeanStdDev(this.alpha1_, this.sig_alpha1);
%             out = this.alpha1_ + normrnd(0, this.sig_alpha1);
        end

        
    end
end

