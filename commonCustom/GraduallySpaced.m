function [x,c,C,IndexFinder] = GraduallySpaced( xStart, xEnd, dxStart, dxEnd)
    L = xEnd - xStart;
    dxAvg = (dxStart+dxEnd)*.5;
    nx = ceil(L/dxAvg)+1;
    if dxStart == dxEnd
        x = linspace(xStart,xEnd,nx)';
        c = []; C = []; IndexFinder=[];
    else
        c = lsqnonlin(@cFunction,1,0,Inf,optimset('display','off'));
        C = ((c*(nx-1)+1)^2 - 1)/L;
        x = zeros(nx,1);
        for i = 1:nx-1
            x(i+1) = (c*i + 1)^2 - 1;
        end
        x = x/C+xStart;
        OneOverC = 1/c;
        IndexFinder = @(x) OneOverC*(sqrt((x-xStart)*C+1)-1)+1;
    end
    function Cost = cFunction(c)
        CC = ((c*(nx-1)+1)^2 - 1);
        dxStartHere = 1/CC*c*(2+c);
        dxEndHere = 1/CC*c*(2+c*(2*nx-3));
        Cost = abs(dxStartHere/dxStart-dxEndHere/dxEnd);
    end
end