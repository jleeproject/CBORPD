classdef GrapheneModelSolver
    % When using SimulParameters, every evaluation, different random parameters will be used.
    
    properties
        name = 'Graphene Model Solver'
        fixedParam = [] % Depends on sensor structure
        changingParam = []; % Depends on sensing target
        changingParam2 = []; % Sensor on Contaminated Water.
        hasSecondChangingParam = false;
        x0=[0.25 0.3 10^-5];
        decisionVariables = GrapheneDecisionVariables;
        sig_err = .7e-7;
        verbose = false;
    end
    
    methods
        function this = GrapheneModelSolver(fixedParameters, changingParam, changingParam2)
            this.fixedParam = fixedParameters;
            if(nargin>1)
                this.changingParam = changingParam;
            end
            if(nargin>2)
                this.changingParam2 = changingParam2; 
                this.hasSecondChangingParam =  true;
            end
        end
        
        function this = setChangingParameters(this, changingParam)
            this.changingParam = changingParam;
        end
        
        function yys_new = solve(this, decisionVariables, Vbs)
            % if # arg = 1, arg1 = decisionVariables
            % if # arg = 2, arg1 = decisionVariables, arg2 = Vbs;
            if(nargin==2)
                Vbs = decisionVariables.Vb;
            end
            nVbs = numel(Vbs);
            
            T = decisionVariables.T;
            L = decisionVariables.L;
            W = decisionVariables.W;
%             t_top = decisionVariables.t_top;
            t_back = decisionVariables.t_back;
%             Vg = decisionVariables.Vg;
%             Vds = decisionVariables.Vds;
%             Vb = decisionVariables.Vb;
            
%             if(nVbs ==1)
%             elseif(nVbs ==0)
%             else
            yys_new = zeros(nVbs,1);
            for i=1:nVbs
                yys_new(i) = solve2(this, T, L, W, t_back, Vbs(i));
%                     yys_new(i) = fn_solve_eq_graphene(Vbs(i),fixedParams,chngParams_Sel, verbose);
            end
%             end
        end
        
        function yys_new = solve2(this, T, L, W, t_back, Vb)            
            yys_new = fn_solve_eq_graphene( this, T, L, W, t_back, Vb );
        end
        
        function [yys_new, yys_on, yys_off ] = evalSameSensorResponseRatioWithVecX(this, X, changingParameterSensorOn, changingParameterSensorOff)
            if(nargin>2)
                [yys_new, yys_on, yys_off ] = evalSameSensorOnOffRatioWithVecX(this, X, changingParameterSensorOn, changingParameterSensorOff, true);
            else
                [yys_new, yys_on, yys_off ] = evalSameSensorOnOffRatioWithVecX(this, X, [], [], true);
            end
        end
        
        function [yys_new, yys_on, yys_off ] = evalSameSensorOnOffRatioWithVecX(this, X, changingParameterSensorOn, changingParameterSensorOff, responseRatio)
%             chngParamOn = GrapheneChangingParametersExporter.getChangingParameter(typeSensorOn);
%             chngParamOff = GrapheneChangingParametersExporter.getChangingParameter(typeSensorOff);
            yys_on = zeros(size(X,1), 1);
            yys_off = zeros(size(X,1), 1);
            yys_new = zeros(size(X,1), 1);
            for i=1:size(X,1);
                T = X(i,1);
                L = X(i,2);
                W = X(i,3);
%                 t_top = X(i,4);
                t_back = X(i,4);
%                 Vg = X(i,6);
%                 Vds = X(i,7);
                Vb = X(i,5);
%                         names = {'T', 'L', 'W', 't_back', 'Vb'};

%                 yys_on = fn_solve_eq_graphene_samesensor( this, T, L, W, t_top, t_back, Vg, Vds, Vb , chngParamOn);
%                 yys_off = fn_solve_eq_graphene_samesensor( this, T, L, W, t_top, t_back, Vg, Vds, Vb , chngParamOff);

                Rd      = this.fixedParam.Rd;
                Rs      = this.fixedParam.Rs;
                Eg      = this.fixedParam.Eg;
                alpha   = this.fixedParam.alpha;
                alpha1  = this.fixedParam.alpha1;
                
                if(nargin>2 && numel(changingParameterSensorOn)>0 )
                    Vb0      = changingParameterSensorOn.Vb0;
                    mu_n     = changingParameterSensorOn.mu_n;
                    mu_p     = changingParameterSensorOn.mu_p;
                    delta    = changingParameterSensorOn.delta;
                else
                    Vb0      = this.changingParam.Vb0;
                    mu_n     = this.changingParam.mu_n;
                    mu_p     = this.changingParam.mu_p;
                    delta    = this.changingParam.delta;
                end

                yys_on = fn_solve_eq_graphene_given_all_parameters( this, T, L, W, t_back, Vb, Rd, Rs, Eg, alpha, alpha1, Vb0, mu_n, mu_p, delta)   ;
                
                if(nargin>2 && numel(changingParameterSensorOff)>0 )
                    Vb0      = changingParameterSensorOff.Vb0;
                    mu_n     = changingParameterSensorOff.mu_n;
                    mu_p     = changingParameterSensorOff.mu_p;
                    delta    = changingParameterSensorOff.delta;
                else
                    Vb0      = this.changingParam2.Vb0;
                    mu_n     = this.changingParam2.mu_n;
                    mu_p     = this.changingParam2.mu_p;
                    delta    = this.changingParam2.delta;
                end

                yys_off = fn_solve_eq_graphene_given_all_parameters( this, T, L, W, t_back, Vb, Rd, Rs, Eg, alpha, alpha1, Vb0, mu_n, mu_p, delta)   ;
                
                if(exist('responseRatio','var') && responseRatio)
                    yys_new(i,1) = (yys_on - yys_off)./yys_on;
                else
                    yys_new(i,1) = yys_on - yys_off;
                end
%                 yys_new(i,1) = yys_on - yys_off;
                yys_on(i,1) = yys_on;
                yys_off(i,1) = yys_off;
            end
        end
        function [yys_diff, yys_on, yys_off ] = evalSameSensorResponseRatioWithVecXVaryingVb(this, X, changingParameterSensorOn, changingParameterSensorOff, Vbs)
            if(nargin>2)
                [yys_diff, yys_on, yys_off ] = evalSameSensorOnOffRatioWithVecXVaryingVb(this, X, changingParameterSensorOn, changingParameterSensorOff, Vbs, true);
            else
                [yys_diff, yys_on, yys_off ] = evalSameSensorOnOffRatioWithVecXVaryingVb(this, X, [], [], Vbs, true);
            end
        end
                function [yys_diff, yys_on, yys_off ] = evalSameSensorOnOffRatioWithVecXVaryingVb(this, X, changingParameterSensorOn, changingParameterSensorOff, Vbs, responseRatio)
%             chngParamOn = GrapheneChangingParametersExporter.getChangingParameter(typeSensorOn);
%             chngParamOff = GrapheneChangingParametersExporter.getChangingParameter(typeSensorOff);
            yys_on = zeros(size(X,1), numel(Vbs));
            yys_off = zeros(size(X,1), numel(Vbs));
            yys_diff = zeros(size(X,1), numel(Vbs));
            for i=1:size(X,1);
                T = X(i,1);
                L = X(i,2);
                W = X(i,3);
                t_back = X(i,4);
%                 Vb = X(i,5);
%                 t_top = X(i,4);
%                 Vg = X(i,6);
%                 Vds = X(i,7);
%                         names = {'T', 'L', 'W', 't_back', 'Vb'};

%                 yys_on = fn_solve_eq_graphene_samesensor( this, T, L, W, t_top, t_back, Vg, Vds, Vb , chngParamOn);
%                 yys_off = fn_solve_eq_graphene_samesensor( this, T, L, W, t_top, t_back, Vg, Vds, Vb , chngParamOff);

                Rd      = this.fixedParam.Rd;
                Rs      = this.fixedParam.Rs;
                Eg      = this.fixedParam.Eg;
                alpha   = this.fixedParam.alpha;
                alpha1  = this.fixedParam.alpha1;
                
                if(nargin>2 && numel(changingParameterSensorOn)>0 )
                    Vb0      = changingParameterSensorOn.Vb0;
                    mu_n     = changingParameterSensorOn.mu_n;
                    mu_p     = changingParameterSensorOn.mu_p;
                    delta    = changingParameterSensorOn.delta;
                else
                    Vb0      = this.changingParam.Vb0;
                    mu_n     = this.changingParam.mu_n;
                    mu_p     = this.changingParam.mu_p;
                    delta    = this.changingParam.delta;
                end
                
                yy_on = fn_solve_eq_graphene_given_all_parameters_varyingVb( this, T, L, W, t_back, Vbs, Rd, Rs, Eg, alpha, alpha1, Vb0, mu_n, mu_p, delta)   ;

                if(nargin>2 && numel(changingParameterSensorOff)>0 )
                    Vb0      = changingParameterSensorOff.Vb0;
                    mu_n     = changingParameterSensorOff.mu_n;
                    mu_p     = changingParameterSensorOff.mu_p;
                    delta    = changingParameterSensorOff.delta;
                else
                    Vb0      = this.changingParam2.Vb0;
                    mu_n     = this.changingParam2.mu_n;
                    mu_p     = this.changingParam2.mu_p;
                    delta    = this.changingParam2.delta;
                end

                yy_off = fn_solve_eq_graphene_given_all_parameters_varyingVb( this, T, L, W, t_back, Vbs, Rd, Rs, Eg, alpha, alpha1, Vb0, mu_n, mu_p, delta)   ;
                yys_on(i,:) = reshape(yy_on , 1,[]);
                yys_off(i,:)  = reshape( yy_off, 1,[]);
                if(exist('responseRatio','var') && responseRatio)
                    yys_diff(i,:) = (yys_on - yys_off)./yys_on;
                else
                    yys_diff(i,:) = yys_on - yys_off;
                end
%                 yys_on(i,:) = yys_on;
%                 yys_off(i,:) = yys_off;
            end
        end
        
        function yys_new = evalDiffSensorResponseRatioWithVecX(this, X, chngParamOn, chngParamOff)
            if(nargin>2)
            	yys_new = evalDiffSensorOnOffRatioWithVecX(this, X, chngParamOn, chngParamOff, true);
            else
            	yys_new = evalDiffSensorOnOffRatioWithVecX(this, X, [], [], true);
            end
        end
        function yys_new = evalDiffSensorOnOffRatioWithVecX(this, X, chngParamOn, chngParamOff, responseRatio)
            yys_new = zeros(size(X,1), 1);
            for i=1:size(X,1);
                T = X(i,1);
                L = X(i,2);
                W = X(i,3);
                t_back = X(i,4);
                Vb = X(i,5);
%                 t_top = X(i,4);
%                 Vg = X(i,6);
%                 Vds = X(i,7);
                if(nargin>2 && numel(chngParamOn)>0 )
                    yys_on = fn_solve_eq_graphene_chngParam( this, T, L, W, t_back, Vb , chngParamOn);
                    yys_off = fn_solve_eq_graphene_chngParam( this, T, L, W, t_back, Vb , chngParamOff);
                else
                    yys_on = fn_solve_eq_graphene_chngParam( this, T, L, W, t_back, Vb , this.changingParam);
                    yys_off = fn_solve_eq_graphene_chngParam( this, T, L, W, t_back, Vb , this.changingParam2);
                end
%                 yys_new(i,1) = yys_on - yys_off;
                if(exist('responseRatio','var') && responseRatio)
                    yys_new(i,1) = (yys_on - yys_off)./yys_on;
                else
                    yys_new(i,1) = yys_on - yys_off;
                end
            end
        end
% 
%         function yys_new = evalResponseRatioWithVecX(this, X, typeSensorPure, typeSensorContam)
%             chngParamPure = GrapheneChangingParametersExporter.getChangingParameter(typeSensorPure);
%             chngParamContam = GrapheneChangingParametersExporter.getChangingParameter(typeSensorContam);
%             yys_new = zeros(size(X,1), 1);
%             for i=1:size(X,1);
%                 T = X(i,1);
%                 L = X(i,2);
%                 W = X(i,3);
%                 t_top = X(i,4);
%                 t_back = X(i,5);
%                 Vg = X(i,6);
%                 Vds = X(i,7);
%                 Vb = X(i,8);
%                 yys_pure = fn_solve_eq_graphene_samesensor( this, T, L, W, t_top, t_back, Vg, Vds, Vb , chngParamPure);
%                 yys_contam = fn_solve_eq_graphene_samesensor( this, T, L, W, t_top, t_back, Vg, Vds, Vb , chngParamContam);
%                 yys_new(i,1) = (yys_pure - yys_contam)./yys_pure;
%             end
%         end

        function yys_new = evalWithVecX(this, X)            
%             T, L, W, t_top, t_back, Vg, Vds, Vb
            yys_new = zeros(size(X,1), 1);
            for i=1:size(X,1);
                T = X(i,1);
                L = X(i,2);
                W = X(i,3);
                t_back = X(i,4);
                Vb = X(i,5);
%                 t_top = X(i,4);
%                 Vg = X(i,6);
%                 Vds = X(i,7);
                yys_new(i,1) = fn_solve_eq_graphene( this, T, L, W, t_back, Vb );
            end
        end
        
        function out = getDomainX(this)
            out = this.decisionVariables.domains.all;
        end

        function yys_new = evalSameSensorVaryingVb(this, decisionVariables, Vbs)            
%             T, L, W, t_top, t_back, Vg, Vds, Vb
            yys_new = zeros(numel(Vbs), 1);
            
            Rd      = this.fixedParam.Rd;
            Rs      = this.fixedParam.Rs;
            Eg      = this.fixedParam.Eg;
            alpha   = this.fixedParam.alpha;
            alpha1  = this.fixedParam.alpha1;
            Vb0      = this.changingParam.Vb0;
            mu_n     = this.changingParam.mu_n;
            mu_p     = this.changingParam.mu_p;
            delta    = this.changingParam.delta;

            
            
            for i=1:numel(Vbs);
                T = decisionVariables.T;
                L = decisionVariables.L;
                W = decisionVariables.W;
%                 t_top = decisionVariables.t_top;
                t_back = decisionVariables.t_back;
%                 Vg = decisionVariables.Vg;
%                 Vds = decisionVariables.Vds;
                Vb = Vbs(i);
%                 yys_new(i,1) = fn_solve_eq_graphene_samesensor( this, T, L, W, t_top, t_back, Vg, Vds, Vb );
                yys_new(i,1) = fn_solve_eq_graphene_given_all_parameters( this, T, L, W, t_back, Vb, Rd, Rs, Eg, alpha, alpha1, Vb0, mu_n, mu_p, delta)   ;
            end
        end
        function [out] = fn_solve_eq_graphene( this, T, L, W, t_back, Vb );

%             try
                if(this.verbose)
                    yyss=fsolve(@(y)fn_graphene(this, y, T, L, W, t_back, Vb ),this.x0);
%                     yyss=fsolve(@(y)fn_graphene(y, Vb, fixedParams, chngParams),x_init);
                else
                    yyss=fsolve(@(y)fn_graphene(this, y, T, L, W, t_back, Vb ),this.x0, optimset('Display','off'));
                end
%                 out = yyss(3);
                out = yyss(3) + normrnd(0,this.sig_err);

%             catch err
%                 disp(err)
%                 out=[];
%             end
        end
        

        
        function [out] = fn_solve_eq_graphene_chngParam( this, T, L, W, t_back, Vb , changingParam);
            this.changingParam = changingParam;
%             try
                if(this.verbose)
                    yyss=fsolve(@(y)fn_graphene(this, y, T, L, W, t_back, Vb ),this.x0);
%                     yyss=fsolve(@(y)fn_graphene(y, Vb, fixedParams, chngParams),x_init);
                else
                    yyss=fsolve(@(y)fn_graphene(this, y, T, L, W, t_back, Vb ),this.x0, optimset('Display','off'));
                end
%                 out = yyss(3);
                out = yyss(3) + normrnd(0,this.sig_err);

%             catch err
%                 disp(err)
%                 out=[];
%             end
        end
        
%         function [out] = fn_solve_eq_graphene_samesensor( this, T, L, W, t_top, t_back, Vg, Vds, Vb, Rd, Rs, Eg, alpha, alpha1, Vb0, mu_n, mu_p, delta, changingParam );
        function [out] = fn_solve_eq_graphene_given_all_parameters( this, T, L, W, t_back, Vb, Rd, Rs, Eg, alpha, alpha1, Vb0, mu_n, mu_p, delta);

            try
                if(this.verbose)
                    yyss=fsolve(@(y)fn_graphene_with_param(this, y, T, L, W, t_back, Vb, Rd, Rs, Eg, alpha, alpha1, Vb0, mu_n, mu_p, delta),this.x0);
%                     yyss=fsolve(@(y)fn_graphene(y, Vb, fixedParams, chngParams),x_init);
                else
                    yyss=fsolve(@(y)fn_graphene_with_param(this, y, T, L, W, t_back, Vb, Rd, Rs, Eg, alpha, alpha1, Vb0, mu_n, mu_p, delta),this.x0, optimset('Display','off'));
                end
%                 out = yyss(3);
                out = yyss(3) + normrnd(0,this.sig_err);

            catch err
                showErrors(err);
%                 disp(err)
                error('ERROR: cannot continue;')
                out=[];
            end
        end
        
        function [out] = fn_solve_eq_graphene_given_all_parameters_varyingVb( this, T, L, W, t_back, Vbs, Rd, Rs, Eg, alpha, alpha1, Vb0, mu_n, mu_p, delta);
            out = zeros(numel(Vbs),1);
            for i=1:numel(Vbs);
                Vb =Vbs(i);
                try
                    if(this.verbose)
                        yyss=fsolve(@(y)fn_graphene_with_param(this, y, T, L, W, t_back, Vb, Rd, Rs, Eg, alpha, alpha1, Vb0, mu_n, mu_p, delta),this.x0);
    %                     yyss=fsolve(@(y)fn_graphene(y, Vb, fixedParams, chngParams),x_init);
                    else
                        yyss=fsolve(@(y)fn_graphene_with_param(this, y, T, L, W, t_back, Vb, Rd, Rs, Eg, alpha, alpha1, Vb0, mu_n, mu_p, delta),this.x0, optimset('Display','off'));
                    end
%                     out(i) = yyss(3);
                    out(i) = yyss(3) + normrnd(0,this.sig_err);

                catch err
                    showErrors(err);
    %                 disp(err)
                    error('ERROR: cannot continue;')
%                     out=[];
                end
            end
        end
        
        function [F] = fn_graphene_with_param(this, x, T, L, W, t_back, Vb, Rd, Rs, Eg, alpha, alpha1, Vb0, mu_n, mu_p, delta)
            %% Constants;
            q     =  GrapheneModelConstants.q;
            hbar  =  GrapheneModelConstants.hbar;

            k     =  GrapheneModelConstants.k;
            Vf    =  GrapheneModelConstants.Vf;
            Vg0   =  GrapheneModelConstants.Vg0;
            m     =  GrapheneModelConstants.m;
            
            t_top =  GrapheneModelConstants.t_top;
            Vg    =  GrapheneModelConstants.Vg;
            Vds   =  GrapheneModelConstants.Vds;

            %% Calculated;
            Ct=GrapheneModelConstants.epso*GrapheneModelConstants.ka_top/t_top;
            Cb=GrapheneModelConstants.epso*GrapheneModelConstants.ka_back/t_back;


            n_pud=delta^2*q^2/(pi*(hbar*Vf)^2);
            a=pi*(k*T)^2/(3*(hbar*Vf)^2);
            b=q^2/(pi*(hbar*Vf)^2);
            c=q^2/(k*T*log(4))^2;
            h=(mu_n+mu_p)/2;
            z=mu_p-mu_n;
            e=0.285;g=0.058;
            f=(q/(5*k*T))^2;
            d=2*q^2*k*T*log(4)/((Ct+Cb)*pi*(hbar*Vf)^2);

            F = zeros(1,3);
            F(1,1)=200*(x(1)*(Ct+Cb)+1/2*Ttrans(x(1),Eg)*2*q^2*k*T*log(4)/pi/(hbar*Vf)^2*sqrt(1+(q*Ttrans(x(1),Eg)/(k*T*log(4)))^2)+Ct*(Vg-Vg0-x(3)*Rs)+Cb*(Vb-Vb0-x(3)*Rs));

            F(1,2)=200*(x(2)*(Ct+Cb)+1/2*Ttrans(x(2),Eg)*2*q^2*k*T*log(4)/pi/(hbar*Vf)^2*sqrt(1+(q*Ttrans(x(2),Eg)/(k*T*log(4)))^2)+Ct*(Vg-Vg0-(Vds-x(3)*Rd))+Cb*(Vb-Vb0-(Vds-x(3)*Rd)));

            func2=@(t)q*(a+b*Ttrans(t,Eg).^2).*(h+14*z*Ttrans(t,Eg)./(sqrt(1+c*Ttrans(t,Eg).^2))).*(m./(m+Ttrans(t,Eg).^2)).*(1+(d*alpha*(2+c*Ttrans(t,Eg).^2.*(3+2*c*Ttrans(t,Eg).^2)))./(1+c*Ttrans(t,Eg).^2).^(1.5));

            func3=@(tt) 1/Vf*(1+f*Ttrans(tt,Eg).^2)./(e+(1+f*Ttrans(tt,Eg).^2)*g).*(m./(m+Ttrans(tt,Eg).^2)).*(h+(14*z*(a+b*Ttrans(tt,Eg).^2).*Ttrans(tt,Eg))./((a+b*Ttrans(tt,Eg).^2+n_pud).*sqrt(1+c*Ttrans(tt,Eg).^2))+(h*d*alpha*(2+c*Ttrans(tt,Eg).^2.*(3+2*c*Ttrans(tt,Eg).^2)))./((1+c*Ttrans(tt,Eg).^2).^(1.5))+(d*alpha*(2+c*Ttrans(tt,Eg).^2.*(3+2*c*Ttrans(tt,Eg).^2))*14*z.*Ttrans(tt,Eg).*(a+b*Ttrans(tt,Eg).^2))./((1+c*Ttrans(tt,Eg).^2).^2.*(a+b*Ttrans(tt,Eg).^2+n_pud)));

            F(1,3)=80*(x(3)-alpha1*(W*(integral(func2,x(1),x(2))+Vds*q*n_pud*h*(m/(m+(x(1)+x(2))^2/4))))/(L+abs(integral(func3,x(1),x(2)))));

%             clear  Rd Rs Eg alpha alpha1       Vb0 mu_n mu_p delta q hbar    T L W k Vf Vg0 m Ct Cb Vg Vds     n_pud a b c h z e f d     func2  func3
        %                 showMemoryUsage(' > @ 1 - 2:') % 
        end
        
        function [F] = fn_graphene(this, x, T, L, W, t_back,Vb)
        %     global q hbar epso T L W k Vf t_top t_back ka_top ka_back Vg0 m Ct Cb Vg;
        %     global q hbar  T L W k Vf Vg0 m Ct Cb Vg;
        %                 showMemoryUsage(' > @ 1 - 1:') % 

            %% Fixed Parameters;
            Rd      = this.fixedParam.Rd;
            Rs      = this.fixedParam.Rs;
            Eg      = this.fixedParam.Eg;
            alpha   = this.fixedParam.alpha;
            alpha1  = this.fixedParam.alpha1;

            %% Changing Parameters;
            Vb0      = this.changingParam.Vb0;
            mu_n     = this.changingParam.mu_n;
            mu_p     = this.changingParam.mu_p;
            delta    = this.changingParam.delta;

            F = fn_graphene_with_param(this, x, T, L, W,  t_back, Vb, Rd, Rs, Eg, alpha, alpha1, Vb0, mu_n, mu_p, delta);
        end




    end
end

