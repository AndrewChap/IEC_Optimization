function [bDirichletControl, phiInitialVecPM, phiInitialVec,phiInitial,efxInitialVec,efrInitialVec] = NormalizedElectrodeVoltages2PhiVecEMGeneral(...
    normalizedElectrodeVoltages,...
    normalizedElectrodePositions,...
    AllowableVoltageRange,...
    normalizedControlRangePositions,...
    PM,bDirichlet,vec2vecPM,vec2grid,...
    EX,ER,...
    nPM,nEM,nx,nr,nxBound,nrInside)
    
    normalizedControlWallVoltage = interp1(normalizedElectrodePositions,normalizedElectrodeVoltages,normalizedControlRangePositions);
        
    controlVec = normalizedControlWallVoltage*(AllowableVoltageRange(2) - AllowableVoltageRange(1)) + AllowableVoltageRange(1);
    controlVec = controlVec(1:end-1);
    
    % Set up control voltage vector of size nPM that will be part of "b" in "Ax=b"
    bControl = zeros(nPM,1);
    rangeControl = (nx-1)*(nrInside-1)-nx+1+(nxBound:nx-1);  %these are the rows of PM that are affected by the control voltage
    rangeInner = (nx-1)*(nrInside-1)+nxBound-2+(1:(nxBound-1):(nxBound-1)*(nr-2-nrInside)+1);   %these are the inner radius ones that have just the first voltage of control
    rangeInner = [rangeInner (nPM-(nxBound-2)):nPM];  % add on the fusion core wall at r = LrZero
    bControl(rangeControl) = controlVec;
    bControl(rangeInner) = controlVec(1);
    
    % Insert known potential values into vector of indices
    phiVecDirichlet = zeros(nEM,1);
    phiVecDirichlet(nx*(nrInside-1)+(nxBound:nx-1)) = controlVec;
    phiVecDirichlet(nx*nrInside + (nxBound:nxBound:nxBound*(nr-nrInside))) = controlVec(1);
    phiVecDirichlet(nEM-(nxBound-1):nEM) = controlVec(1);
    
    % Find initial potential
    phiInitialVecPM = PM\(bDirichlet.*bControl);
    phiInitialVec = phiVecDirichlet;
    phiInitialVec(vec2vecPM) = phiInitialVecPM;
    
    % Find initial electric field
    efxInitialVec = single(-EX*phiInitialVec);
    efrInitialVec = single(-ER*phiInitialVec);
    
    phiInitialVecPM = single(phiInitialVecPM);
    phiInitialVec = single(phiInitialVec);
    phiInitial = single(zeros(nx,nr));
    phiInitial(vec2grid) = phiInitialVec;
    
    % Define control vector used for simulation
    bDirichletControl = single(bDirichlet.*bControl);
    


end