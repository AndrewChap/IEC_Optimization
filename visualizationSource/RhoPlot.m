classdef RhoPlot < handle
    % Contains creation and updating of the electric potential plot
    % Also plots the circle of neutralization if applicable
    % Does not scatter the particle positions atop the colorplot, that is
    % done by UpdateParticlePositions
    
    properties
        parent
        PlotHandle
        SubPlot
        AxesHandle
        minDens
        maxDens = -Inf;
    end
    
    properties (Constant)

    end
    
    methods
        function this = RhoPlot(parent,xGrid,rGrid,size1,size2)
            this.parent = parent;
            
            this.SubPlot = subplot(3,3,4:5);
            this.SubPlot.Position = this.SubPlot.Position + [-.06, -.03, 0.09, 0];
            this.PlotHandle = surf(xGrid,rGrid,0*xGrid);
            this.PlotHandle.LineStyle = 'none';
            this.PlotHandle.FaceColor = 'interp';
            view(0,90);
            
            rhoTitle = title('$\log_{10}(\textup{Density}\,[\textup{m}^{-3}])$','interpreter','latex','fontsize',size1);
            rhoTitle.Units = 'normalized';
            rhoTitle.Position = [.55  .9   0];
            xlabel('$x\,[\textup{m}]$','interpreter','latex','fontsize',size2)
            ylabel('$r\,[\textup{m}]$','interpreter','latex','fontsize',size2)
            this.AxesHandle = gca;
            grid off
            rhoBar = colorbar;
            rhoBar.Position = [.605 rhoBar.Position(2:4)];
            axis equal
            axis([0 1 0 .3]*parent.L)
        end
    end
    
    methods
        function UpdateRhoPlot(this,rhoVecPM)
            
            rhoVec = single(zeros(size(this.parent.vec2grid)));
            rhoVec(this.parent.vec2vecPM) = double(rhoVecPM)./double(this.parent.cellVolumesPM);
            rho = single(zeros(size(this.parent.xGrid)));
            rho(this.parent.vec2grid) = rhoVec;
            log10rho = log10(rho);
            
            this.maxDens = max(max(max(log10rho)),this.maxDens);
            this.minDens = this.maxDens - 5;
            log10rho(log10rho<this.minDens) = this.minDens;
            
            this.PlotHandle.CData = log10rho;

            if this.maxDens > this.minDens
                caxis(this.AxesHandle,[this.minDens  this.maxDens]);
            end
            
        end
    end
    
end

