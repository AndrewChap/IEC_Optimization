classdef ScatterParticles
    % Contains creation and updating of the particle scatters for all plots
    
    properties
        parent
        XXRX
        XXVX_R
        XXVX_L
        VXVRVT_R
        VXVRVT_L
        VRVT_S
        VXVT_S
        VXVR_S
        particleSizeMaster = 3;
    end
    
    properties (Constant)
       limVXplot = 0.15;
       limVRplot = 0.1;
       limVTplot = 0.1;
    end
    
    methods
        function this = ScatterParticles(parent,partColors,size3,size4,size8)
            this.parent = parent;
            
            subplot(this.parent.PhiHandle.SubPlot);
            this.XXRX = DrewScatter([NaN,NaN],2,[0 0 0]);
            
            xxvxSub = subplot(3,4,9);
            xxvxSub.Position = xxvxSub.Position + [0, -.03, 0, 0];
            hold on
            for i = 1:2
                legDummies(i) = plot(NaN,NaN,'LineStyle','none','Marker','.','MarkerSize',10,'Color',partColors(i+1,:));
            end
            
            this.XXVX_R = DrewScatter([NaN,NaN],this.particleSizeMaster,partColors(2,:));
            this.XXVX_L = DrewScatter([NaN,NaN],this.particleSizeMaster,partColors(3,:));
            leg = legend(legDummies,'$v_x\geq0$','$v_x<0$');
            leg.Interpreter = 'latex';
            leg.FontSize = size8;
            leg.Location = 'southeast';
            xlabel('$x\,[\textup{m}]$','interpreter','latex','fontsize',size4)
            ylabel('$|v_x|/v_\textup{fusion}$','interpreter','latex','fontsize',size4)
            title('$(x,v_x)$ Phase Space','interpreter','latex','fontsize',size3)
            
            xlim([0 .2*parent.L])
            ylim([.9  1.1])
            grid on
            
            vxvrvtSub = subplot(3,4,10);
            limVXplot = .15;
            limVRplot = .1;
            limVTplot = .1;
            vxvrvtSub.Position = vxvrvtSub.Position + [0, -.03, 0, 0];
            hold on
%             for i = 1:2
%                 legDummies(i) = plot(NaN,NaN,'LineStyle','none','Marker','.','MarkerSize',dummySize,'Color',partColors(i+1,:));
%             end
            this.VXVRVT_R = DrewScatter([NaN,NaN,NaN],this.particleSizeMaster,partColors(2,:));
            this.VXVRVT_L = DrewScatter([NaN,NaN,NaN],this.particleSizeMaster,partColors(3,:));
            this.VRVT_S = DrewScatter([limVRplot*ones(parent.np,1),limVTplot*ones(parent.np,1),-limVXplot*ones(parent.np,1)],this.particleSizeMaster,partColors(1,:));
            this.VXVT_S = DrewScatter([limVRplot*ones(parent.np,1),limVTplot*ones(parent.np,1),-limVXplot*ones(parent.np,1)],this.particleSizeMaster,partColors(1,:));
            this.VXVR_S = DrewScatter([limVRplot*ones(parent.np,1),limVTplot*ones(parent.np,1),-limVXplot*ones(parent.np,1)],this.particleSizeMaster,partColors(1,:));
            leg.Interpreter = 'latex';
            leg.Location = 'southeast';
            leg.FontSize = size8;
            xlabel('$v_r$','interpreter','latex','fontsize',size4)
            ylabel('$v_\theta$','interpreter','latex','fontsize',size4)
            zlabel('$v_x$','interpreter','latex','fontsize',size4)
            title('Local $(v_x,v_r,v_\theta)$ Space','interpreter','latex','fontsize',size3)
            view(-24,36)
            axis equal
            
            zlim([-limVXplot limVXplot])
            ylim([-limVTplot limVTplot])
            xlim([-limVRplot limVRplot])
            grid on
        end
    end
    
    methods
        function UpdateParticles(this,np,px,pr,vx,vr,vt)
            
            vx = vx/this.parent.fusionVelocity;
            vr = vr/this.parent.fusionVelocity;
            vt = vt/this.parent.fusionVelocity;
            
            % velocities must be normalized to the fusion velocity!
            range = 1:np;
            vxMean = mean(abs(vx(range)));
            this.XXRX.XData = px(range);
            this.XXRX.YData = pr(range);
               
            pRight = vx(range)>0;
            pLeft = vx(range)<=0;
    
            this.XXVX_R.XData = px(pRight);
            this.XXVX_R.YData = vx(pRight);
            this.XXVX_L.XData = px(pLeft);
            this.XXVX_L.YData = abs(vx(pLeft));
            
            this.VXVRVT_R.XData = vr(pRight);
            this.VXVRVT_R.YData = vt(pRight);
            this.VXVRVT_R.ZData = abs(vx(pRight))-vxMean;
            this.VXVRVT_L.XData = vr(pLeft);
            this.VXVRVT_L.YData = vt(pLeft);
            this.VXVRVT_L.ZData = abs(vx(pLeft))-vxMean;
            this.VRVT_S.XData = vr(range);
            this.VRVT_S.YData = vt(range);
            this.VXVT_S.YData = vt(range);
            this.VXVT_S.ZData = abs(vx(range))-vxMean;
            this.VXVR_S.XData = vr(range);
            this.VXVR_S.ZData = abs(vx(range))-vxMean;
        end
    end
    
end

