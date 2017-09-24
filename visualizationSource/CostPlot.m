classdef CostPlot
    % Contains creation and updating of the cost plot of the main figure
    % includes differnt colors for the different optimizers
    
    properties
        costPlotGlobal
        costPlotLocal
        periodLines
    end
    
    methods
        function this = CostPlot(color1,color2,size1)
        
            costSub = subplot(3,3,6);
            costSub.Position = costSub.Position + [-.01, -.03, 0, 0];
            hold on
            this.periodLines = loglog(NaN,NaN,'color',[.8 .6 .6]);
             
            this.costPlotGlobal = loglog(NaN,NaN,'o','MarkerSize',3,'color',color1,'MarkerFaceColor',color1);
            this.costPlotLocal = loglog(NaN,NaN,'o','MarkerSize',3,'color',color2,'MarkerFaceColor',color2);
            
            
            
            costXlabel = xlabel('iteration \#','interpreter','latex','fontsize',size1);
            costYlabel = ylabel('Cost','interpreter','latex','fontsize',size1);
            costYlabel.Units = 'normalized';
            costYlabel.Position = [-.09, .5, 0];
                    
            title('Cost function output','interpreter','latex','fontsize',size1)
            cA = gca;
            cA.XScale = 'linear';
            cA.YScale = 'log';
            costSub.YLim = [1e-5 1];
            costSub.YTick = [1e-5 1e-4 1e-3 1e-2 1e-1 1];
            grid on
            
        end
        
        function UpdateCostPlot(this,costTracker,costMode,periodLinesXdata)
            nperiod = numel(periodLinesXdata);
            periodX = zeros(3*nperiod,1);
            periodX(1:3:3*nperiod) = periodLinesXdata;
            periodX(2:3:3*nperiod) = periodLinesXdata;
            periodX(3:3:3*nperiod) = NaN;
            periodY = repmat([1e-5;1;NaN],nperiod,1);
            this.periodLines.XData = periodX;
            this.periodLines.YData = periodY;
            this.costPlotGlobal.YData = costTracker.*costMode;
            this.costPlotGlobal.XData = 1:numel(costTracker);
            this.costPlotLocal.YData = costTracker.*(1-costMode);
            this.costPlotLocal.XData = 1:numel(costTracker);
        end
    end
    
end

