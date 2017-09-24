classdef ElectrodePlot
    % Contains creation and updating of the cost plot of the main figure
    % includes differnt colors for the different optimizers
    
    properties
        parent
        voltagePlot
        bestVoltagePlot
        fusionVoltagePlot
    end
    
    methods
        function this = ElectrodePlot(parent,color1,color2,size1)
            this.parent = parent;
            evSub = subplot(3,3,3);
            evSub.Position = evSub.Position + [-.01, -.03, 0, 0];
            hold on
            this.fusionVoltagePlot = plot([...
                parent.ElectrodeRadialPositions(1),...
                parent.ElectrodeRadialPositions(end),...
                NaN,...
                parent.ElectrodeRadialPositions(1),...
                parent.ElectrodeRadialPositions(end)],NaN(1,5),'-.','color',[.4 .6 .4]);
            plot(parent.ElectrodeRadialPositions,parent.LB,'--','color',[.5 .5 .5]);
            plot(parent.ElectrodeRadialPositions,parent.UB,'--','color',[.5 .5 .5]);
            this.bestVoltagePlot = plot(parent.ElectrodeRadialPositions,NaN*parent.ElectrodeRadialPositions,'-o','color',color1,'MarkerFaceColor',color1);
            this.voltagePlot = plot(parent.ElectrodeRadialPositions,NaN*parent.ElectrodeRadialPositions,'-o','color',color2,'MarkerFaceColor',color2);
            xlabel('ElectrodePosition','interpreter','latex','fontsize',size1)
            ylabel('ElectrodeVoltage','interpreter','latex','fontsize',size1)
            ylim([parent.AllowableVoltageRange(1) parent.AllowableVoltageRange(2)])
            title('Wall electrode voltage','interpreter','latex','fontsize',size1)
            grid on
        end
        
        function UpdateElectrodePlot(this,bestElectrodeVoltagesToPlot,ElectrodeVoltagesToPlot,fusionVoltage)
            ElectrodeVoltagesRescaled = ElectrodeVoltagesToPlot*(this.parent.AllowableVoltageRange(2)-this.parent.AllowableVoltageRange(1)) + this.parent.AllowableVoltageRange(1);
            bestElectrodeVoltagesRescaled = bestElectrodeVoltagesToPlot*(this.parent.AllowableVoltageRange(2)-this.parent.AllowableVoltageRange(1)) + this.parent.AllowableVoltageRange(1);
            this.bestVoltagePlot.YData = bestElectrodeVoltagesRescaled;
            this.voltagePlot.YData = ElectrodeVoltagesRescaled;
            this.fusionVoltagePlot.YData = [...
                ElectrodeVoltagesRescaled(1),...
                ElectrodeVoltagesRescaled(1),...
                NaN,...
                ElectrodeVoltagesRescaled(1)+fusionVoltage,...
                ElectrodeVoltagesRescaled(1)+fusionVoltage];
                
        end
    end
    
end

