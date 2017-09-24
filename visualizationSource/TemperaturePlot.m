classdef TemperaturePlot < handle
    % Contains creation and updating of the cost plot of the main figure
    % includes differnt colors for the different optimizers
    
    properties
        parent
        SubPlot
        plotTx
        plotTr
        plotTt
        plotTm
    end
    
    methods
        function this = TemperaturePlot(parent,color1,color2,color3,size1)
            
            this.parent = parent;
            this.SubPlot = subplot(3,4,11:12);
            this.SubPlot.Position = this.SubPlot.Position + [0, -.03, 0, 0];
            hold on
            this.plotTx = semilogy(NaN,NaN,'s-','color',color1,'MarkerFaceColor',color1);
            this.plotTr = semilogy(NaN,NaN,'s-','color',color2,'MarkerFaceColor',color2);
            this.plotTt = semilogy(NaN,NaN,'s-','color',color3,'MarkerFaceColor',color3);
            this.plotTm = semilogy(NaN,NaN,'s-','color','k','MarkerFaceColor','k');
            xlabel('$t\,[\mu\textup{s}]$','interpreter','latex','fontsize',size1)
            yLabelPower = ylabel('$T/E_0$','interpreter','latex','fontsize',size1);
            yLabelPower.Units = 'normalized';
            
            tempTitle = title('Temperatures - \color[rgb]{0.85  0  0}T_x, \color[rgb]{0 0.4470 0.7410}T_r, \color[rgb]{.3  .7  .15}T_\theta, \color[rgb]{0 0 0}T_{total}','interpreter','tex','fontsize',12);
            tempTitle.Units = 'normalized';
            tempTitle.Position = [.5 .98 0];
            
            eA = gca;
            set(eA,'YScale','log')
            grid on
            
            this.SubPlot.YLim = [1e-7 1];
            this.SubPlot.YTick = [1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1];
        end
        
        function UpdateTemperaturePlot(this,timeTemp,tempX,tempR,tempT)
%             keyboard
            this.plotTx.XData = timeTemp*1e6;
            this.plotTx.YData = tempX;
            this.plotTr.XData = timeTemp*1e6;
            this.plotTr.YData = tempR;
            this.plotTt.XData = timeTemp*1e6;
            this.plotTt.YData = tempT;
            this.plotTm.XData = timeTemp*1e6;
            this.plotTm.YData = tempX+tempR+tempT;
            try
                this.SubPlot.XLim = [0 max(timeTemp)*1e6];
            catch
                this.SubPlot.XLim = [0 1]; % if we are at the first timestep
            end

        end
    end
    
end

