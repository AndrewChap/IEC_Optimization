classdef EnergyPlot
    % Contains creation and updating of the conservation of energy plot of the main figure
    
    properties
        plotKE
        plotPE
        plotTE
    end
    
    methods
        function this = EnergyPlot(color1,color2,size1)
        
            costSub = subplot(3,3,6);
            costSub.Position = costSub.Position + [-.01, -.03, 0, 0];
            hold on
            this.plotKE = plot(NaN,NaN,'color',color1);
            this.plotPE = plot(NaN,NaN,'color',color2);
            this.plotTE = plot(NaN,NaN,'color','k');
            
            tempEnergy = title('Energies - \color[rgb]{0.85  0  0}KE, \color[rgb]{0 0.4470 0.7410}PE, \color[rgb]{0 0 0}KE+PE','interpreter','tex','fontsize',12);
%             title('Cons. of Energy','interpreter','latex','fontsize',size1)
            xlabel('$t\,[\mu\textup{s}]$','interpreter','latex','fontsize',size1)
            ylabel('$E/E_0$','interpreter','latex','fontsize',size1);
            grid on
            
        end
        
        function UpdateEnergyPlot(this,time,KE,PE,Eloss)
            this.plotKE.XData = time*1e6;
            this.plotKE.YData = KE;
            this.plotPE.XData = time*1e6;
            this.plotPE.YData = PE;
            this.plotTE.XData = time*1e6;
            this.plotTE.YData = KE+PE+0*Eloss;  %stopped using Eloss because it only has the potential energy of the lost particles, not the kinetic.
        end
    end 
end