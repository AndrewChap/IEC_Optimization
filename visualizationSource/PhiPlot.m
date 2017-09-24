classdef PhiPlot < handle
    % Contains creation and updating of the electric potential plot
    % Also plots the circle of neutralization if applicable
    % Does not scatter the particle positions atop the colorplot, that is
    % done by UpdateParticlePositions
    
    properties
        parent
        PlotHandle
        SubPlot
    end
    
    properties (Constant)
       zoffcircle = linspace(0,pi/2,50); 
    end
    
    methods
        function this = PhiPlot(parent,xGrid,rGrid,ElectrodeCoordinates,AllowableVoltageRange,zOffRadius,size1,size2,color1)
            this.parent = parent;
            circleX = zOffRadius*cos(this.zoffcircle);
            circleY = zOffRadius*sin(this.zoffcircle);
            
            
            this.SubPlot = subplot(3,3,1:2);
            this.SubPlot.Position = this.SubPlot.Position + [-.06, -.03, 0.09, 0];
            
            ElectrodeScatter = DrewScatter(ElectrodeCoordinates);
            ElectrodeScatter.MarkerSize = 30;
            ElectrodeScatter.Color = color1;
            ElectrodeScatter.MarkerFaceColor = color1;
            ElectrodeScatter.MarkerEdgeColor = color1;
            hold on
            this.PlotHandle = surf(xGrid,rGrid,0*xGrid);
            this.PlotHandle.LineStyle = 'none';
            this.PlotHandle.FaceColor = 'interp';
            view(0,90);
            caxis([AllowableVoltageRange(1) AllowableVoltageRange(2)])
            phiTitle = title('Electric Potential [V]','interpreter','latex','fontsize',size1);
            phiTitle.Units = 'normalized';
            phiTitle.Position = [.55  .9   0];
            xlabel('$x\,[\textup{m}]$','interpreter','latex','fontsize',size2)
            ylabel('$r\,[\textup{m}]$','interpreter','latex','fontsize',size2)
            grid off
            
            phiBar = colorbar;
            phiBar.Position = [.605 phiBar.Position(2:4)];
            axis equal
            axis([0 1 0 .3]*xGrid(end,1))
            drawnow
            
            if (any(circleX~=0))
                patchline(circleX,circleY,zeros(size(this.zoffcircle)),'edgecolor','w','linewidth',3,'edgealpha',0.3);
            end
        end
    end
    
    methods
        function UpdatePhiPlot(this,phiVec)
            phi = single(zeros(size(this.parent.xGrid)));
            phi(this.parent.vec2grid) = phiVec;
            this.PlotHandle.CData = phi;
        end
    end
    
end

