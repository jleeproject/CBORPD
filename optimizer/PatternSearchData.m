classdef PatternSearchData < handle
    properties
        x_history = [];
        y_obj_history = [];
        y_con_history = [];
        y_obs_history = [];
        feasgap_history = [];
        x_iter_history = [];
        y_iter_history = [];
        iter_history = [];
%         x_prev = [];
        cnt_eval_obj = 0;
        cnt_eval_con = 0;
        cnt_eval = 0;
        y_obj_map = containers.Map();
        y_cons_map = containers.Map();
        y_obs_map = containers.Map();
        
        verbose = false;
        
        simulator
        samplesize
        feasgap
    end
    methods
        function this = PatternSearchData(simulator, samplesize, feasgap)
            this.simulator = simulator;
            this.samplesize = samplesize;
%             this.x_prev = x_train(:,end);
            this.feasgap = feasgap;
            
            this.x_history = [];
            this.y_obj_history = [];
            this.y_con_history = [];
            this.y_obs_history = [];
            this.feasgap_history = [];
            this.x_iter_history = [];
            this.y_iter_history = [];
            this.iter_history = [];
            this.cnt_eval_obj = 0;
            this.cnt_eval_con = 0;
            this.cnt_eval = 0;
            
            this.y_obj_map = containers.Map();
            this.y_cons_map = containers.Map();
            this.y_obs_map = containers.Map();
        end
        
        function out = objfcn(this, x)
            if this.verbose
                fprintf('[obj] x=(');
                fprintf('%g\t',x);
                fprintf(')\n');
            end
            if(~this.y_obj_map.isKey(num2str(x)))
                [y_obj, y_cons, y_obs] = MultipleSamplerObjCon.evaluate(this.simulator, x, this.samplesize);
                this.x_history = [this.x_history; x];
                this.y_obj_history = [this.y_obj_history; y_obj];
                this.y_con_history = [this.y_con_history; y_cons];
                this.y_obs_history = [this.y_obs_history; y_obs];
                this.feasgap_history = [this.feasgap_history; this.feasgap(y_cons)];
                this.cnt_eval = this.cnt_eval +1;
%                 this.y_obj_map(num2str(x)) = y_obj;
                this.y_cons_map(num2str(x)) = y_cons;
%                 this.y_obs_map(num2str(x)) = y_obs;
            else
                y_obj = this.y_obj_map(num2str(x));
                this.y_obj_map.remove(num2str(x));
%                 y_cons = this.y_cons_map(num2str(x));
%                 y_obs = this.y_obs_map(num2str(x));
                if this.verbose
                    disp('USED Existing y_obs');
                end
            end
            out = y_obj;
%             (out);
            if this.verbose
                fprintf('[obj] y=%g\n',out);
            end
            if(numel(out)>1)
                disp('WARNING');
            end
            this.cnt_eval_obj = this.cnt_eval_obj + 1;
        end

        function [out, ceq] = confcn(this, x)
            if this.verbose
                fprintf('[con] x=(');
                fprintf('%g\t',x);
                fprintf(')\n');
            end
            if(~this.y_cons_map.isKey(num2str(x)))
                this.x_history = [this.x_history; x];
                [y_obj, y_cons, y_obs] = MultipleSamplerObjCon.evaluate(this.simulator, x, this.samplesize);
                this.y_obj_history = [this.y_obj_history; y_obj];
                this.y_con_history = [this.y_con_history; y_cons];
                this.y_obs_history = [this.y_obs_history; y_obs];
                this.feasgap_history = [this.feasgap_history; this.feasgap(y_cons)];
                this.cnt_eval = this.cnt_eval +1;
                this.y_obj_map(num2str(x)) = y_obj;
%                 this.y_cons_map(num2str(x)) = y_cons;
%                 this.y_obs_map(num2str(x)) = y_obs;
            else
%                 y_obj = this.y_obj_map(num2str(x));
                y_cons = this.y_cons_map(num2str(x));
                this.y_cons_map.remove(num2str(x));
%                 y_obs = this.y_obs_map(num2str(x));
                if this.verbose
                    disp('USED Existing y_cons');
                end
            end
%             out = this.y_cons;
            out = this.feasgap(y_cons);
%             this.feasgap_history = [this.feasgap_history; out];
            ceq = [];
            if this.verbose
                fprintf('[con] y=%g\n',out);
            end
            this.cnt_eval_con = this.cnt_eval_con + 1;
        end 
        
        function [stop, nothing, change] = myoutput(this, x,optimvalues,state);
            nothing = optimvalues;
            change = false;
            stop = false;
            if isequal(state,'iter')
              this.x_iter_history = [this.x_iter_history; x.x];
              this.y_iter_history = [this.y_iter_history; x.fval];
              this.iter_history = [this.iter_history; x];
            elseif isequal(state,'init')
                if this.verbose
                    disp('init');
                end
            else
                if this.verbose
                    fprintf('[warning] exception state=%s',state);
                end
%                 disp(x);
%                 disp(optimvalues);
%                 disp(state);
            end
        end

    end
end