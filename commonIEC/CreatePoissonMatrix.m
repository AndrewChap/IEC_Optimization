function [PM,PM_DIA,PM_OFF,EX,EX_DIA,EX_OFF,ER,ER_DIA,ER_OFF,bDirichlet] = CreatePoissonMatrix2DCylindricalNonUniformSkewMissingChunkCUDA4(x,r,dirichlet)

    % Drew Chap 11/28/2016
    % updated 4/19/2017 for the DIA format for CUDA (https://www.bu.edu/pasi/files/2011/01/NathanBell1-10-1000.pdf)
    %  Coefficients taken from "Finite Difference schemes on non-Uniform meshes
    %  for Hyperbolic Conservation Laws" Nikolaos Sfakianakis PhD thesis 2009
    %  (http://www2.mathematik.uni-mainz.de/~sfaknikj/files/SfakianakisThesis.pdf)
    % update 2/6/17 now can be done with domains that omit parts.  really
    % only designed to work with my particular IEC simulation
    % update 4/28/17 Setting outer R boundary (core region) to the innermost potential
    % update 5/1/17 CUDA vectors for electric field added
    % update 5/15/17 there should be no electric field on the neumann boundaries
    
    nx = size(x,1);
    nr = size(x,2);
    scale = (x(2,1)-x(1,1))^-2;

    nxBulk = sum(isinf(dirichlet(:,1)));
    nrBulk = sum(isinf(dirichlet(1,:)));
    nxShort = sum(isinf(dirichlet(:,end-1)));
    nrShort = sum(isinf(dirichlet(end-1,:)));
    nBulk = nxBulk*nrShort;
    
    
    % number of matrix rows.  includes all interior nodes as well as 
    N = sum(sum(isinf(dirichlet)));
    
    bDirichlet = zeros(N,1);
    PM = sparse(N,N);
    PM_DIA = single(zeros(size(PM,1),7));   %This code uses a seven point stencil
    PM_OFF = int32(zeros(size(PM,1),7));   %We might not need all 7 offset rows.  The constant ones can be hard-coded in
%     
%     PM_DIA = zeros(size(PM,1),7);   %This code uses a seven point stencil
%     PM_OFF = zeros(size(PM,1),7);   %We might not need all 7 offset rows.  The constant ones can be hard-coded in

    
    for a = 1:N
        if a <= nBulk
            i = mod(a-1,nxBulk)+1;
            j = ceil(a/(nxBulk));
        else
            i = mod(a-nBulk-1,nxShort)+1;
            j = ceil((a-nBulk)/nxShort)+nrShort;
        end
         
        if i == 1 && j == 1  % reflective boundary condition at origin
            hr = r(1,2) - r(1,1);
            hx = x(2,1) - x(1,1);
            PM(a,a) = -2/(hx*hx)-2/(hr*hr);
            PM(a,a+1) = 2/(hx*hx);
            PM(a,a+nxBulk) = 2/(hr*hr);
            PM_OFF(a,5) = 1;
            PM_OFF(a,7) = nxBulk;
            PM_DIA(a,4) = -2/(hx*hx)-2/(hr*hr);
            PM_DIA(a,5) = 2/(hx*hx);
            PM_DIA(a,7) = 2/(hr*hr);
        elseif i == 1 && j == nrBulk  % reflective in x direction, dirichlet in r direction
            hr = r(1,nrBulk) - r(1,nrBulk-1);
            hx = x(2,nrBulk) - x(1,nrBulk);
            PM(a,a) = -2/(hx*hx)-2/(hr*hr);
            PM(a,a-nxShort) = -1/r(i,j)*1/(2*hr) + 1/(hr*hr);
            PM(a,a+1) = 2/(hx*hx);
            PM_OFF(a,1) = -nxShort;
            PM_OFF(a,5) = 1;
            PM_DIA(a,1) = -1/r(i,j)*1/(2*hr) + 1/(hr*hr);
            PM_DIA(a,4) = -2/(hx*hx)-2/(hr*hr);
            PM_DIA(a,5) = 2/(hx*hx);
            bDirichlet(a) = bDirichlet(a) - (1/r(i,j)*1/(2*hr) + 1/(hr*hr));
        elseif i == 1 % reflective boundary condition along x=0 
            hr1 = r(i,j) - r(i,j-1);
            hr2 = r(i,j+1) - r(i,j);
            hx = x(2,j) - x(1,j);  
            
            alphaR1 = -hr2/(hr1*(hr1+hr2));
            gammaR1 = hr1/(hr2*(hr1+hr2));
            betaR1 = -alphaR1-gammaR1;
            %coefficients for second derivative of r(page 10 of Sfakianakis thesis)
            alphaR2 = 2/(hr1*(hr1+hr2));
            gammaR2 = 2/(hr2*(hr1+hr2));
            betaR2 = -alphaR2-gammaR2;
            
            PM(a,a) = 1/r(i,j)*betaR1 + betaR2 - 2/(hx*hx);
            PM_DIA(a,4) = 1/r(i,j)*betaR1 + betaR2 - 2/(hx*hx);
            
            PM(a,a+1) = 2/(hx*hx);
            PM_OFF(a,5) = 1;
            PM_DIA(a,5) = 2/(hx*hx);
            
            if j <= nrShort+1
                PM(a,a-nxBulk) = 1/r(1,j)*alphaR1 + alphaR2;
                PM_OFF(a,1) = -nxBulk;
                PM_DIA(a,1) = 1/r(1,j)*alphaR1 + alphaR2;
            else
                PM(a,a-nxShort) = 1/r(1,j)*alphaR1 + alphaR2;
                PM_OFF(a,1) = -nxShort;
                PM_DIA(a,1) = 1/r(1,j)*alphaR1 + alphaR2;
            end
            if j <= nrShort
                PM(a,a+nxBulk) = 1/r(i,j)*gammaR1 + gammaR2;
                PM_OFF(a,7) = nxBulk;
                PM_DIA(a,7) = 1/r(i,j)*gammaR1 + gammaR2;
            else
                PM(a,a+nxShort) = 1/r(i,j)*gammaR1 + gammaR2;
                PM_OFF(a,7) = nxShort;
                PM_DIA(a,7) = 1/r(i,j)*gammaR1 + gammaR2;
            end
            
        %elseif i == nx  % dirichlet here so don't need to cover
        elseif j == 1   % at r = 0, use values from chao 1997 paper
            hr = r(i,j+1) - r(i,j);
            hx1 = x(i,j) - x(i-1,j);
            hx2 = x(i+1,j) - x(i,j);
            alphaX2 = 2/(hx1*(hx1+hx2));
            gammaX2 = 2/(hx2*(hx1+hx2));
            betaX2 = -alphaX2-gammaX2;
            
            PM(a,a) = -4/(hr*hr) + betaX2;
            PM(a,a+nxBulk)  = 4/(hr*hr);
            PM(a,a-1) = alphaX2;
            PM_DIA(a,4) = -4/(hr*hr) + betaX2;
            PM_DIA(a,7) = 4/(hr*hr);
            PM_DIA(a,3) = alphaX2;
            PM_OFF(a,7) = nxBulk;
            PM_OFF(a,3) = -1;
            if isinf(dirichlet(i+1,j))
                PM(a,a+1) = gammaX2;
                PM_DIA(a,5) = gammaX2;
                PM_OFF(a,5) = 1;
            else
                bDirichlet(a) = gammaX2*(dirichlet(i+1,j)==1);
            end
            
        else
            hr1 = r(i,j) - r(i,j-1);
            hr2 = r(i,j+1) - r(i,j);
            hx1 = x(i,j) - x(i-1,j);
            hx2 = x(i+1,j) - x(i,j);
            
            %coefficients for first derivative of r (page 9 of Sfakianakis thesis)
            alphaR1 = -hr2/(hr1*(hr1+hr2));
            gammaR1 = hr1/(hr2*(hr1+hr2));
            betaR1 = -alphaR1-gammaR1;
            %coefficients for second derivative of r(page 10 of Sfakianakis thesis)
            alphaR2 = 2/(hr1*(hr1+hr2));
            gammaR2 = 2/(hr2*(hr1+hr2));
            betaR2 = -alphaR2-gammaR2;
            %coefficients for first derivative of x (page 9 of Sfakianakis thesis)
            alphaX1 = -hx2/(hx1*(hx1+hx2));
            gammaX1 = hx1/(hx2*(hx1+hx2));
            betaX1 = -alphaX1-gammaX1;
            %coefficients for second derivative of x(page 10 of Sfakianakis thesis)
            alphaX2 = 2/(hx1*(hx1+hx2));
            gammaX2 = 2/(hx2*(hx1+hx2));
            betaX2 = -alphaX2-gammaX2;
            %linear interpolation stuff needed for skewed grid (skewed only in the r-direction)
            rPoint = r(i,j);
            Al = min(max(1 - (rPoint-r(i-1,j-1))/(r(i-1,j)-r(i-1,j-1)),0),1);
            Cl = min(max((rPoint-r(i-1,j))/(r(i-1,j+1)-r(i-1,j)),0),1);
            Bl = 1 - Al - Cl;
            Ar = min(max(1 - (rPoint-r(i+1,j-1))/(r(i+1,j)-r(i+1,j-1)),0),1);
            Cr = min(max((rPoint-r(i+1,j))/(r(i+1,j+1)-r(i+1,j)),0),1);
            Br = 1 - Ar - Cr;
            
            if j <= nrShort+1
                nxLeft = nxBulk;
            else
                nxLeft = nxShort;
            end
            if j <= nrShort
                nxRight = nxBulk;
            else
                nxRight = nxShort;
            end
            
            PM(a,a) = 1/r(i,j)*betaR1 + betaR2 + betaX2;
            PM_DIA(a,4) = 1/r(i,j)*betaR1 + betaR2 + betaX2;

            if size(PM,2) == 3546
                keyboard
            end
            if isinf(dirichlet(i-1,j))
                PM(a,a-1) = Bl*alphaX2;
                PM_DIA(a,3) = Bl*alphaX2;
                PM_OFF(a,3) = -1;
            else
                bDirichlet(a) = bDirichlet(a) - Bl*alphaX2*(dirichlet(i-1,j)==1);
            end
            if isinf(dirichlet(i-1,j+1))
                PM(a,a-1+nxRight) = Cl*alphaX2;
                PM_DIA(a,6) = Cl*alphaX2;
                PM_OFF(a,6) = nxRight-1;
            else
                bDirichlet(a) = bDirichlet(a) - Cl*alphaX2*(dirichlet(i-1,j+1)==1);
            end
            if isinf(dirichlet(i+1,j-1))
                PM(a,a+1-nxLeft) = Ar*gammaX2;
                PM_DIA(a,2) = Ar*gammaX2;
                PM_OFF(a,2) = -nxLeft+1;
            else
                bDirichlet(a) = bDirichlet(a) - Ar*gammaX2*(dirichlet(i+1,j-1)==1);
            end
            if isinf(dirichlet(i+1,j))
                PM(a,a+1) = Br*gammaX2;
                PM_DIA(a,5) = Br*gammaX2;
                PM_OFF(a,5) = 1;
            else
                bDirichlet(a) = bDirichlet(a) - Br*gammaX2*(dirichlet(i+1,j)==1);
            end

            if isinf(dirichlet(i,j-1))
                PM(a,a-nxLeft) = 1/r(i,j)*alphaR1 + alphaR2;
                PM_DIA(a,1) = 1/r(i,j)*alphaR1 + alphaR2;
                PM_OFF(a,1) = -nxLeft;
            else
                bDirichlet(a) = bDirichlet(a) - (1/r(i,j)*alphaR1 + alphaR2)*(dirichlet(i,j-1)==1);
            end
            if isinf(dirichlet(i,j+1))
                PM(a,a+nxRight) = 1/r(i,j)*gammaR1 + gammaR2;
                PM_DIA(a,7) = 1/r(i,j)*gammaR1 + gammaR2;
                PM_OFF(a,7) = nxRight;
            else
                bDirichlet(a) = bDirichlet(a) - (1/r(i,j)*gammaR1 + gammaR2)*(dirichlet(i,j+1)==1);
            end

        end
        if size(PM,2) == 3546
            keyboard
        end
    end
    
    % number of matrix rows in EX and EM.  Includes all interior nodes and
    % boundaries
    N = sum(sum(~isnan(dirichlet)));
    [EX,ER] =  deal(sparse(N,N));
    EX_DIA = zeros(N,5);
    EX_OFF = int32(zeros(size(PM,1),5));
    ER_DIA = zeros(N,3);
    ER_OFF = int32(zeros(size(PM,1),3));
    
    nxBulk = sum(~isnan(dirichlet(:,1)));
    nrBulk = sum(~isnan(dirichlet(1,:)));
    nxShort = sum(~isnan(dirichlet(:,end-1)));
    nrShort = sum(~isnan(dirichlet(end-1,:)));
    nBulk = nxBulk*nrShort;
    
    for a = 1:N
        if a <= nBulk
            i = mod(a-1,nxBulk)+1;
            j = ceil(a/(nxBulk));
        else
            i = mod(a-nBulk-1,nxShort)+1;
            j = ceil((a-nBulk)/nxShort)+nrShort;
        end
        
        if i == 1 && j == 1
%             hr = r(1,2) - r(1,1);
%             hx = x(2,1) - x(1,1);
%             EX(a,a) = -1/hx;
%             EX(a,a+1) = 1/hx;
%             EX_DIA(a,3) = -1/hx;
%             EX_DIA(a,4) = 1/hx;
%             EX_OFF(a,4) = 1;
%             ER(a,a) = -1/hr;
%             ER(a,a+nxBulk) = 1/hr;
%             ER_DIA(a,2) = -1/hr;
%             ER_DIA(a,3) = 1/hr;
%             ER_OFF(a,3) = nxBulk;
        elseif i == 1 && j == nr
            hr = r(1,nr) - r(1,nr-1);
%             hx = x(2,nr) - x(1,nr);
%             EX(a,a) = -1/hx;
%             EX(a,a+1) = 1/hx;
%             EX_DIA(a,3) = -1/hx;
%             EX_DIA(a,4) = 1/hx;
%             EX_OFF(a,4) = 1;
            ER(a,a) = 1/hr;
            ER(a,a-nxShort) = -1/hr;
            ER_DIA(a,2) = 1/hr;
            ER_DIA(a,1) = -1/hr;
            ER_OFF(a,1) = -nxShort;
        elseif i == nx && j == 1 
%             hr = r(nx,2) - r(nx,1);
            hx = x(nx,1) - x(nx-1,1);
            EX(a,a) = 1/hx;
            EX(a,a-1) = -1/hx;
            EX_DIA(a,3) = 1/hx;
            EX_DIA(a,2) = -1/hx;
            EX_OFF(a,2) = -1;
%             ER(a,a) = -1/hr;
%             ER(a,a+nxBulk) = 1/hr;
%             ER_DIA(a,2) = -1/hr;
%             ER_DIA(a,3) = 1/hr;
%             ER_OFF(a,3) = nxBulk;
        elseif i == nx && j == nrShort 
            hr = r(nx,nrShort) - r(nx,nrShort-1);
            hx = x(nx,nrShort) - x(nx-1,nrShort);
            EX(a,a) = 1/hx;
            EX(a,a-1) = -1/hx;
            EX_DIA(a,3) = 1/hx;
            EX_DIA(a,2) = -1/hx;
            EX_OFF(a,2) = -1;
            ER(a,a) = 1/hr;
            ER(a,a-nxBulk) = -1/hr;
            ER_DIA(a,2) = 1/hr;
            ER_DIA(a,1) = -1/hr;
            ER_OFF(a,1) = -nxBulk;
        elseif i == nxShort && j == nr 
            hr = r(nxShort,nr) - r(nxShort,nr-1);
            hx = x(nxShort,nr) - x(nxShort-1,nr);
            EX(a,a) = 1/hx;
            EX(a,a-1) = -1/hx;
            EX_DIA(a,3) = 1/hx;
            EX_DIA(a,2) = -1/hx;
            EX_OFF(a,2) = -1;
            ER(a,a) = 1/hr;
            ER(a,a-nxShort) = -1/hr;
            ER_DIA(a,2) = 1/hr;
            ER_DIA(a,1) = -1/hr;
            ER_OFF(a,1) = -nxShort;
        elseif i == 1
            hr1 = r(1,j) - r(1,j-1);
            hr2 = r(1,j+1) - r(1,j);
%             hx = x(2,j) - x(1,j);  
            alphaR1 = -hr2/(hr1*(hr1+hr2));
            gammaR1 = hr1/(hr2*(hr1+hr2));
            betaR1 = -alphaR1-gammaR1;
            
%             EX(a,a) = -1/hx;
%             EX(a,a+1) = 1/hx;
%             
%             EX_DIA(a,3) = -1/hx;
%             EX_DIA(a,4) = 1/hx;
%             EX_OFF(a,4) = 1;
            
            ER(a,a) = betaR1;
            ER_DIA(a,2) = betaR1;
            
            if j <= nrShort+1
                ER(a,a-nxBulk) = alphaR1;
                ER_DIA(a,1) = alphaR1;
                ER_OFF(a,1) = -nxBulk;
            else
                ER(a,a-nxShort) = alphaR1;
                ER_DIA(a,1) = alphaR1;
                ER_OFF(a,1) = -nxBulk;
            end
            if j <= nrShort
                ER(a,a+nxBulk) = gammaR1;
                ER_DIA(a,3) = gammaR1;
                ER_OFF(a,3) = nxBulk;
            else
                ER(a,a+nxShort) = gammaR1;
                ER_DIA(a,3) = gammaR1;
                ER_OFF(a,3) = nxShort;
            end
            
            
            
        elseif i == nx
            hr1 = r(nx,j) - r(nx,j-1);
            hr2 = r(nx,j+1) - r(nx,j);
            hx = x(nx,j) - x(nx-1,j);
            alphaR1 = -hr2/(hr1*(hr1+hr2));
            gammaR1 = hr1/(hr2*(hr1+hr2));
            betaR1 = -alphaR1-gammaR1;
            
            EX(a,a) = 1/hx;
            EX(a,a-1) = -1/hx;
            EX_DIA(a,3) = 1/hx;
            EX_DIA(a,2) = -1/hx;
            EX_OFF(a,2) = -1;
            
            ER(a,a) = betaR1;
            ER(a,a-nxBulk) = alphaR1;
            ER(a,a+nxBulk) = gammaR1;
            ER_DIA(a,2) = betaR1;
            ER_DIA(a,1) = alphaR1;
            ER_DIA(a,3) = gammaR1;
            ER_OFF(a,1) = -nxBulk;
            ER_OFF(a,3) = nxBulk;
        elseif i == nxShort && j > nrShort   % need to adjust this for skewed r
            hr1 = r(nxShort,j) - r(nxShort,j-1);
            hr2 = r(nxShort,j+1) - r(nxShort,j);
            hx = x(nxShort,j) - x(nxShort-1,j);
            alphaR1 = -hr2/(hr1*(hr1+hr2));
            gammaR1 = hr1/(hr2*(hr1+hr2));
            betaR1 = -alphaR1-gammaR1;
            
            EX(a,a) = 1/hx;
            EX(a,a-1) = -1/hx;
            EX_DIA(a,3) = 1/hx;
            EX_DIA(a,2) = -1/hx;
            EX_OFF(a,2) = -1;
            
            ER(a,a) = betaR1;
            ER_DIA(a,2) = betaR1;
            if j == nrShort+1
                ER(a,a-nxBulk) = alphaR1;
                ER_DIA(a,1) = alphaR1;
                ER_OFF(a,1) = -nxBulk;
            else
                ER(a,a-nxShort) = alphaR1;
                ER_DIA(a,1) = alphaR1;
                ER_OFF(a,1) = -nxShort;
            end
            ER(a,a+nxShort) = gammaR1;
            ER_DIA(a,3) = gammaR1;
            ER_OFF(a,3) = nxShort;
        elseif j == 1
%             hr = r(i,2) - r(i,1);
            hx1 = x(i,j) - x(i-1,j);
            hx2 = x(i+1,j) - x(i,j);
            alphaX1 = -hx2/(hx1*(hx1+hx2));
            gammaX1 = hx1/(hx2*(hx1+hx2));
            betaX1 = -alphaX1-gammaX1;
            
            EX(a,a) = betaX1;
            EX(a,a-1) = alphaX1;
            EX(a,a+1) = gammaX1;
            EX_DIA(a,3) = betaX1;
            EX_DIA(a,2) = alphaX1;
            EX_DIA(a,4) = gammaX1;
            EX_OFF(a,2) = -1;
            EX_OFF(a,4) = 1;
            
%             ER(a,a) = -1/hr;
%             ER(a,a+nxBulk) = 1/hr;
%             ER_DIA(a,2) = -1/hr;
%             ER_DIA(a,3) = 1/hr;
%             ER_OFF(a,3) = nxBulk;
            
        elseif j == nrShort && i > nxShort  % need to adjust this for skewed r
            hr = r(i,nrShort) - r(i,nrShort-1);
            hx1 = x(i,j) - x(i-1,j);
            hx2 = x(i+1,j) - x(i,j);
            alphaX1 = -hx2/(hx1*(hx1+hx2));
            gammaX1 = hx1/(hx2*(hx1+hx2));
            betaX1 = -alphaX1-gammaX1;
            
            EX(a,a) = betaX1;
            EX(a,a-1) = alphaX1;
            EX(a,a+1) = gammaX1;
            EX_DIA(a,3) = betaX1;
            EX_DIA(a,2) = alphaX1;
            EX_DIA(a,4) = gammaX1;
            EX_OFF(a,2) = -1;
            EX_OFF(a,4) = 1;
            
            ER(a,a) = 1/hr;
            ER(a,a-nxBulk) = -1/hr;
            ER_DIA(a,2) = 1/hr;
            ER_DIA(a,1) = -1/hr;
            ER_OFF(a,1) = -nxBulk;
        elseif j == nr
            hr = r(i,nr) - r(i,nr-1);
            hx1 = x(i,nr) - x(i-1,nr);
            hx2 = x(i+1,nr) - x(i,nr);
            alphaX1 = -hx2/(hx1*(hx1+hx2));
            gammaX1 = hx1/(hx2*(hx1+hx2));
            betaX1 = -alphaX1-gammaX1;
            
            EX(a,a) = betaX1;
            EX(a,a-1) = alphaX1;
            EX(a,a+1) = gammaX1;
            EX_DIA(a,3) = betaX1;
            EX_DIA(a,2) = alphaX1;
            EX_DIA(a,4) = gammaX1;
            EX_OFF(a,2) = -1;
            EX_OFF(a,4) = 1;
            
            ER(a,a) = 1/hr;
            ER(a,a-nxShort) = -1/hr;
            ER_DIA(a,2) = 1/hr;
            ER_DIA(a,1) = -1/hr;
            ER_OFF(a,1) = -nxShort;
        else
            hr1 = r(i,j) - r(i,j-1);
            hr2 = r(i,j+1) - r(i,j);
            hx1 = x(i,j) - x(i-1,j);
            hx2 = x(i+1,j) - x(i,j);
            
            %coefficients for first derivative of r (page 9 of Sfakianakis thesis)
            alphaR1 = -hr2/(hr1*(hr1+hr2));
            gammaR1 = hr1/(hr2*(hr1+hr2));
            betaR1 = -alphaR1-gammaR1;
            %coefficients for first derivative of x (page 9 of Sfakianakis thesis)
            alphaX1 = -hx2/(hx1*(hx1+hx2));
            gammaX1 = hx1/(hx2*(hx1+hx2));
            betaX1 = -alphaX1-gammaX1;
            %linear interpolation stuff needed for skewed grid (skewed only in the r-direction)
            rPoint = r(i,j);
            Al = min(max(1 - (rPoint-r(i-1,j-1))/(r(i-1,j)-r(i-1,j-1)),0),1);   %weight of the point (i-1,j-1)
            Cl = min(max((rPoint-r(i-1,j))/(r(i-1,j+1)-r(i-1,j)),0),1);         %weight of the point (i-1,j+1)
            Bl = 1 - Al - Cl;                                                   %weight of the point (i-1,j)
            Ar = min(max(1 - (rPoint-r(i+1,j-1))/(r(i+1,j)-r(i+1,j-1)),0),1);   %weight of the point (i+1,j-1)
            Cr = min(max((rPoint-r(i+1,j))/(r(i+1,j+1)-r(i+1,j)),0),1);         %weight of the point (i+1,j+1)
            Br = 1 - Ar - Cr;                                                   %weight of the point (i+1,j)
            
            if j <= nrShort+1
                nxLeft = nxBulk;
            else
                nxLeft = nxShort;
            end
            if j <= nrShort
                nxRight = nxBulk;
            else
                nxRight = nxShort;
            end
            
            EX(a,a) = betaX1;
            EX(a,a-1) = Bl*alphaX1;
%             EX(a,a-1-nxLeft) = Al*alphaX1;
            EX(a,a-1+nxRight) = Cl*alphaX1;
            EX(a,a+1) = Br*gammaX1;
            EX(a,a+1-nxLeft) = Ar*gammaX1;
%             EX(a,a+1+nxRight) = Cr*gammaX1;
            
            EX_DIA(a,3) = betaX1;
            EX_DIA(a,2) = Bl*alphaX1;
            EX_DIA(a,5) = Cl*alphaX1;
            EX_DIA(a,4) = Br*gammaX1;
            EX_DIA(a,1) = Ar*gammaX1;
            EX_OFF(a,2) = -1;
            EX_OFF(a,5) = -1+nxRight;
            EX_OFF(a,4) = 1;
            EX_OFF(a,1) = +1-nxLeft;
            
            ER(a,a) = betaR1;
            ER(a,a-nxLeft) = alphaR1;
            ER(a,a+nxRight) = gammaR1;
            ER_DIA(a,2) = betaR1;
            ER_DIA(a,1) = alphaR1;
            ER_DIA(a,3) = gammaR1;
            ER_OFF(a,1) = -nxLeft;
            ER_OFF(a,3) = nxRight;
        end
        
    end
    
end