classdef Visualizer_BO <handle & matlab.mixin.Copyable
    %VISUALIZER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        isGridSet = false;
        gridCoord
        gridSquareD1
        gridSquareD2
        gridDomainD1
        gridDomainD2
        gridBin
        %
        trueGlobalMinF
        trueGlobalMinX
        fnEval
        %
%         options
        showAllGrid = false;
        showSearchTrails = true;
        %
%         isReady = true;
    end
    
    methods
%         function obj = Visualizer(inputArg1,inputArg2)
        function this = Visualizer_BO()
%             this.options = struct();
%             this.options.alpha = 0.5;
%             this.options.color = [125 125 125]./255;
%             this.options.str_title='Bayesian Optimization';
        end
        
        function setOtherGridFieldsByCoord(this)
            this.gridSquareD1=reshape(this.gridCoord(:,1),this.gridBin,this.gridBin);
            this.gridSquareD2=reshape(this.gridCoord(:,2),this.gridBin,this.gridBin);
            this.gridDomainD1=this.gridSquareD1(:,1);
            this.gridDomainD2= this.gridSquareD2(1,:)';
%             zz1sq=reshape(xs(:,1),dim,dim);
%             zz2sq=reshape(xs(:,2),dim,dim);
%             zz1_1d=zz1sq(:,1);
%             zz2_1d= zz2sq(1,:)';
%             if(numel(this.fnEval)>0)
%                 args = mat2arg(this.gridCoord);
%                 [true_fmin,true_min_idx ] = min(log(this.fnEval(args{:})));
%                 true_minx = this.gridCoord(true_min_idx,:);
%                 setTrueGlobalMinValueNX(this, true_fmin, true_minx)
%             else
%                 disp('ERROR: fn_eval not set');
%             end
        end
        
        function setTrueGlobalMinValueNX(this, true_fmin, true_minx)
            this.trueGlobalMinF = true_fmin;
            this.trueGlobalMinX = true_minx;
        end
        
        
        function ready = isReady(this)
            ready = true;
            if (this.showAllGrid &&~ this.isGridSet )
                ready = false;
                disp('Please allocate grids.')
            end
        end
        %--------------------- [gridCoord] ---------------------
        function domain = getGridDomains(this,counts, xDomain)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            dmn = [1:counts];
            a = [  ( xDomain(:,2) - xDomain(:,1) )/(counts-1)   ];
            b = [  ( xDomain(:,1).*counts - xDomain(:,2) )/(counts-1)   ];
            domain = a.*dmn+b;
        end
        
        function output = getGridInCoordCombinationArr(this,counts, objFunc)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            % COMBINATIONS
%             if(objFunc.dim>2)
%                 xdim = objFunc.dim ;
%                 mainCnt = ceil(sqrt((counts.^2)./2.^(xdim-2)));
%                 anxCnt = 2;
%                 
%                 xDomain = objFunc.getXDomain();
%                 
%                 mainDomain = this.getGridDomains(mainCnt, xDomain(1:2,:) );
%                 anxDomain = this.getGridDomains(anxCnt, xDomain(3:end,:)  );
%                 
%                 xs = combvec(mainDomain(1,:), mainDomain(2,:));
%                 for i = 3:objFunc.getDim()
%                     xs = combvec(xs, anxDomain(i-2,:));
%                 end
%             else
                if( isa(objFunc,'AbsFunction'))
                    domain = this.getGridDomains(counts, objFunc.getXDomain());
                elseif( isa(objFunc,'GrapheneModelSolver'))
                    domain = this.getGridDomains(counts, objFunc.decisionVariables.domains.all);
                else
                    error('[AbsInfillOptimizer] Undefined Type.');
                end
%                 domain = this.getGridDomains(counts, objFunc.getXDomain());
                xs = combvec(domain(1,:), domain(2,:));
%                 for i = 2:objFunc.getDim()
%                     xs = combvec(xs, domain(i,:));
%                 end                
%             end
%             xs = combvec(domain(1,:), domain(1,:));
            output = xs';
        end
        
        
        function setGrid(this,counts, objFunc)
            this.gridCoord = getGridInCoordCombinationArr(this,counts, objFunc );
%             this.gridBin = sqrt(size(this.gridCoord,1));
            this.gridBin = counts;
            setOtherGridFieldsByCoord(this)
            this.isGridSet = true;
        end

        function output = getGridCoord(this,current_x, plotIdx)
            if(size(current_x,2)==2)
                output = this.gridCoord;
            else
                if(nargin>1)
                    if(nargin<3)
                        plotIdx = [1 2];
                    else
                        if(sum(size(current_x,2)<plotIdx))
                            plotIdx = [1 2];
                            fprintf('[Warning] Plot index is larger than the problem dimension\n plotIdx  is replaced by [1 2].');
                        end
                    end
                    %This is for multiple Grid.
                    if(size(current_x,1)==1)
                        if( sum(plotIdx == [1,2])==2)
                            output = [this.gridCoord (current_x(3:end).*ones(size(this.gridCoord,1),1))];
                        else
                            output = (current_x.*ones(size(this.gridCoord,1),1));
                            output(:,plotIdx) = this.gridCoord ;
    %                         output = [this.gridCoord (current_x(3:end).*ones(size(this.gridCoord,1),1))];
                        end
                    else
                        error('cannot have X with multiple rows');
                    end
                else
                    output = this.gridCoord;
                end
            end
        end


        %--------------------- [gridSquareD1] ---------------------
        function setGridSquareD1(this,gridSquareD1)
            this.gridSquareD1 = gridSquareD1;
        end

        function output = getGridSquareD1(this)
            output = this.gridSquareD1;
        end


        %--------------------- [gridSquareD2] ---------------------
        function setGridSquareD2(this,gridSquareD2)
            this.gridSquareD2 = gridSquareD2;
        end

        function output = getGridSquareD2(this)
            output = this.gridSquareD2;
        end


        %--------------------- [gridDomainD1] ---------------------
        function setGridDomainD1(this,gridDomainD1)
            this.gridDomainD1 = gridDomainD1;
        end

        function output = getGridDomainD1(this)
            output = this.gridDomainD1;
        end


        %--------------------- [gridDomainD2] ---------------------
        function setGridDomainD2(this,gridDomainD2)
            this.gridDomainD2 = gridDomainD2;
        end

        function output = getGridDomainD2(this)
            output = this.gridDomainD2;
        end
            
        %--------------------- [trueGlobalMinF] ---------------------
        function setTrueGlobalMinF(this,trueGlobalMinF)
            this.trueGlobalMinF = trueGlobalMinF;
        end

        function output = getTrueGlobalMinF(this)
            output = this.trueGlobalMinF;
        end


        %--------------------- [trueGlobalMinX] ---------------------
        function setTrueGlobalMinX(this,trueGlobalMinX)
            this.trueGlobalMinX = trueGlobalMinX;
        end

        function output = getTrueGlobalMinX(this)
            output = this.trueGlobalMinX;
        end

        %--------------------- [fnEval] ---------------------
        function setFnEval(this,fnEval)
            this.fnEval = fnEval;
        end

        function output = getFnEval(this)
            output = this.fnEval;
        end
        
        %--------------------- [showAllGrid] ---------------------
        function setShowAllGrid(this,showAllGrid)
            if(nargin>1)
                this.showAllGrid = showAllGrid;
            else
                this.showAllGrid  = true;
            end
        end

        function output = isShowAllGrid(this)
            output = this.showAllGrid;
        end

        %--------------------- [dim] ---------------------
        function setDim(this,dim)
            this.gridBin = dim;
        end

        function output = getDim(this)
            output = this.gridBin;
        end

        %--------------------- [showSearchTrails] ---------------------
        function setShowSearchTrails(this,showSearchTrails)
            this.showSearchTrails = showSearchTrails;
        end

        function output = isShowSearchTrails(this)
            output = this.showSearchTrails;
        end

%         function outputArg = setGridCoord(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

