function vols = FindCellVolumesSkewMissingChunk(x,r)
    
    % Finds cell volumes by splitting up the cell grid into a fine cell
    % grid with nodes defining the half boundaries
    nx = size(x,1);
    nr = size(x,2);
    
    nxFine = 2*nx-1;
    nrFine = 2*nr-1;
    
    xCenters = x(:,1);
    xBoundaries = 0.5*(xCenters(2:end)+xCenters(1:end-1));
    xSpace = sort([xCenters;xBoundaries]);
    xFine = repmat(xSpace,1,nrFine);
    rFine = zeros(size(xFine));
    for i = 1:nx-1
        ii = 2*i-1;
        rCentersL = r(i,:);
        rBoundariesL = 0.5*(rCentersL(2:end)+rCentersL(1:end-1));
        rSpaceL = sort([rCentersL,rBoundariesL]);
        rCentersR = r(i+1,:);
        rBoundariesR = 0.5*(rCentersR(2:end)+rCentersR(1:end-1));
        rSpaceR = sort([rCentersR,rBoundariesR]);
        rFine(ii:ii+2,:) = [rSpaceL;0.5*(rSpaceL+rSpaceR);rSpaceR];
    end
%     %%
%     xFineVals = 0*xFine;
%     xFineVals(1:2:end,1:2:end) = NaN;
%     xFineVals2 = 0*xFine;
%     xFineVals2(2:end-1,2:end-1) = NaN;
%     nrFineInside = find(isnan(rFine(end,:)),1)-1
%     nxFineInside = find(isnan(rFine(:,end)),1)-1
%     xFineVals2(nxFineInside:end,nrFineInside) = 0;
%     xFineVals2(nxFineInside,nrFineInside:end) = 0;
% %     Plot cell centers and boundaries
%     figure(324)
%     cla
% %     ma= mesh(xFine([1,end],[1,end]),rFine([1,end],[1,end]),0*xFine([1,end],[1,end]));
% %     ma.EdgeColor = 'k';
%     hold on
%     mb1 = mesh(xFine,rFine,xFineVals2);
%     mb1.FaceAlpha = 0;
%     mb1.EdgeColor = [.3 .3 .3];
%     mb2 = mesh(xFine,rFine,xFineVals);
%     mb2.EdgeColor = [.6 .6 .6 ];
% %     mb = mesh(xFine([1,2:2:end-1,end],[1,2:2:end-1,end]),rFine([1,2:2:end-1,end],[1,2:2:end-1,end]),0*xFine([1,2:2:end-1,end],[1,2:2:end-1,end]));
%     
%     dots = scatter(x(:),r(:),8,'filled');
% %     dots2 = scatter(xFine(:),rFine(:),3,[.2 .2 .7],'filled');
% %     mb.EdgeColor = [.7 .3 .3];
% %     mb.FaceAlpha = 0;
%     view(0,90)
%     axis equal
%     axis off
%     set(gcf,'color','w')
%     matlab2tikz('DomainGrid.tex')
%     export_fig('DomainGrid','-pdf')
%     hold on
    %%
%     keyboard
    
    xFineBuffer = nan(nxFine+2,nrFine+2);
    rFineBuffer = nan(nxFine+2,nrFine+2);
    xFineBuffer(2:end-1,2:end-1) = xFine;
    rFineBuffer(2:end-1,2:end-1) = rFine;
    volsFine = zeros(nxFine+1,nrFine+1);
    for j = 1:nrFine+1
        for i = 1:nxFine+1
            volsFine(i,j) = FindVolumeExtrude(xFineBuffer(i,j),xFineBuffer(i+1,j),rFineBuffer(i,j),rFineBuffer(i,j+1),rFineBuffer(i+1,j),rFineBuffer(i+1,j+1));
        end
    end
    volsFine(isnan(volsFine))=0;
    vols = nan(nx,nr);
    for j = 1:nr
        for i = 1:nx
            ii = 2*i-1;
            jj = 2*j-1;
            vols(i,j) = sum(sum(volsFine(ii:ii+1,jj:jj+1)));
        end
    end
end

function volume = FindVolumeExtrude(x1,x2,r1i,r1o,r2i,r2o)
    dx = x2-x1;
    a = (r2o-r1o)/dx;
    b = r1o;
    c = (r2i-r1i)/dx;
    d = r1i;
    volume = pi*(dx*(b^2 - d^2) + dx*dx*(a*b - c*d) + dx*dx*dx*(a^2-c^2)/3);
end