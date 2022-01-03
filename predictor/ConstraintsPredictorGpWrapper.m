classdef ConstraintsPredictorGpWrapper
    
    properties
        cell_predictors
        nPredictor
    end
    
    methods
        
        function this = ConstraintsPredictorGpWrapper(cell_predictors)
            this.cell_predictors = cell_predictors;
            this.nPredictor = numel(cell_predictors);
        end
        
       
        function [mus,sigs,covMats] = predict(this, xx)
            nPredictor = this.nPredictor;
%             mus = cell(nPredictor, 1);
%             covMats = cell(nPredictor, 1);
%             sigs = cell(nPredictor, 1);
            nEvals = size(xx,1);
            mus = zeros(nEvals, nPredictor);
%             covMats = zeros(nEvals, nPredictor);
            covMats = cell(nPredictor,1);
            sigs = zeros(nEvals, nPredictor);
            for i=1:this.nPredictor
%                 [mus{i}, covMats{i}, sigs{i}] = this.cell_predictors{i}.predict(xx);
%                 [mus(:,i), covMats(:,i), sigs(:,i)] = this.cell_predictors{i}.predict(xx);
%                 [mu, sig, covMat] = this.cell_predictors{i}.predict(xx);
%                 mus(:,i) = mu;
%                 sigs(:,i) = sig;
%                 covMats{i} = covMat;
                [mus(:,i), sigs(:,i), covMats{i}] = this.cell_predictors{i}.predict(xx);
            end
        end
    end
end

