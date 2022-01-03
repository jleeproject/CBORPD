classdef SampleSizeSelector
    %SAMPLESIZESELECTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type
        isFixed
        fixedSamplesize = [];
    end
    
    methods
        function this = SampleSizeSelector(type, samplesize)
            this.type = type;
            switch type
                case TypeSampleSize.Fixed
                    this.isFixed = true;
                    this.fixedSamplesize = samplesize;
                case TypeSampleSize.Adjust
                    this.isFixed = false;
                    this.fixedSamplesize = [];
                otherwise
                    throwError('[SampleSizeSelector] Not defined type');
            end
        end
        
        function samplesize = selectSamplesize(this, xnew, predictor)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
%             samplesize = obj.Property1 + inputArg;
            switch this.type
                case TypeSampleSize.Fixed
                    samplesize = this.fixedSamplesize;
                case TypeSampleSize.Adjust
                    [~, sig_xnew] = predictor.predict(xnew);
                    samplesize = getAdjustedSamplesize(sig_xnew);
                otherwise
                    throwError('[SampleSizeSelector] Not defined type');
            end
        end
    end
end

