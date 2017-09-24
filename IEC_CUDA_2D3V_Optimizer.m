clear all
close all
clc
% Main file, run this one!
%
% Choose "OptimizeMode = true" to run the optimizer and plot the simulation
% state each time the cost function is evaluated
%
% Choose "OptimizeMode = false" to run the simulation with the optimized
% electrodes as given by "preOptimizedElectrodes.mat"

%%%%%%%%%%%%%%%%%%%%%
OptimizeMode = false;
%%%%%%%%%%%%%%%%%%%%%

% type any useful notes here.  It will appear in notes.txt in your output files folder
userNotes = 'testing opt results WITH mag ON but raising CFL, also fixed energy x-axis (time) problem';

addpath('commonIEC')        % functions specific to this simulation
addpath('commonCustom')     % functions I wrote
addpath('commonFEX')        % functions downloaded from MATLAB file exchange
addpath('commonFiles')      % non-code files
addpath('magneticField')    % functions and data for magnetic field

gitInfo = getGitInfo;

if ~OptimizeMode
%     load('preOptimizedElectrodes.mat')
%     load('F:\Users\TeslaSSD\Documents\MATLAB\IEC_CUDA_2D3V_Optimizer\outputFiles\Opt2017-09-20_15-58-54\optimizedVoltages\bestVoltgagesPeriod8.mat')
    load('F:\Users\TeslaSSD\Documents\MATLAB\IEC_CUDA_2D3V_Optimizer\outputFiles\Opt2017-09-20_16-02-31\optimizedVoltages\bestVoltgagesPeriod8.mat')
end

if OptimizeMode
    oString = 'Opt';
else
    oString = 'Run';
end
filename = [mfilename,oString];
timename = clockstring;

loadPrevMode = true;  % when optimization of the new period starts, load the particle positions and velocities from the previous optimization results (see thesis)
energyBump = true;    % at period completion, modify energy so that our energy doesn't continuoually increase over time
if ~OptimizeMode
    loadPrevMode = false;  % won't do anything and will mess up initial temperature value if you try loadPrevMode when not in OptimizeMode
    energyBump = false;
end

mexFile = 'IEC_ParticleInCell';       % MEX file name (don't include extension)
rng('default')  % seed random number generator for reproducability
cudaRandSeed = int32(0);  % cuda random number seed for reproducability

% define physical constants
e = 1.6e-19;
eV = 1.6e-19;
AMU = 1.67e-27;
eps0 = 8.85e-12;

%define fusion fuel
m = (11/6)*AMU; % ion mass
mu = m/2;       % reduce mass of m-m pair
Z = sqrt(5);    % ionization level
q = Z*e;        % ion charge
fusionEnergyReleased = 8.7e6*eV;  % fusion energy released per fusion event (J)
fusionPeakRelativeVelocity = 1.1126e7;
fusionVelocity = fusionPeakRelativeVelocity/2;
fusionKE = .5*m*fusionVelocity.^2;
fusionVoltage = fusionKE/q;

%define fusor parameters
L = 1;  % length scale (m)
a_coeff = q^2/(4*pi*eps0*mu);
Remanence = 1;
MagnetThickness = 0.02;             % in radians, changing this will change the magnetization but won't change the wall angle (ideally it should)
numReal = 2e9*L;                    % number of real particles in the simulation (scales nicely with L, see thesis for explanation)
Lx = 1*L;                           % non-dimensionalized length units
innerRadius = .25*L;                % non-dimensionalized length units
outerRadius = 1*L;                  % non-dimensionalized length units
angleDegrees = 17.6885;             % found in IEC diagram (with WallThickness = 0.04)
angle = deg2rad(angleDegrees);

% initial distribution variations in position and velocity
sigX = .02*L;
sigR = .01*L;
sigVX = .001;
sigVR = .001;

% initial number of particles
if OptimizeMode
    np = 5000;
    CFL = .5;               % Courant-Friedrichs-Lewy condition
    FixEnergy = true;
else
    np = 100000;
    CFL = .06;
    FixEnergy = false;
end

maxParticlesPerCell = int32(ceil(np/6)); 

itPM = int32(30);           % # of Poisson matrix iterations

% define allowable voltages at inner and outer radii
LowerVoltageInside = -3*fusionVoltage;
HigherVoltageInside = -.8*fusionVoltage;
LowerVoltageOutside = -.6*fusionVoltage;
HigherVoltageOutside = 0;

%define electrode positions
ElectrodeRadialPositions = GraduallySpaced(innerRadius, outerRadius, .05*L, .16*L)';
%define intial guess of voltages a function of position
ElectrodeVoltages = (scale01(ElectrodeRadialPositions).^2-1)*fusionVoltage*1.3;
% LB and UB for each voltage to be used in optimizer
LB = scale01(ElectrodeRadialPositions)*(LowerVoltageOutside-LowerVoltageInside) + LowerVoltageInside;
UB = scale01(ElectrodeRadialPositions)*(HigherVoltageOutside-HigherVoltageInside) + HigherVoltageInside;
AllowableVoltageRange = [min(LB) max(UB)];

if ~OptimizeMode && exist('bestElectrodeVoltagesToOptimize')
    ElectrodeVoltages = [bestElectrodeVoltagesToOptimize 1]*(AllowableVoltageRange(2) - AllowableVoltageRange(1)) + AllowableVoltageRange(1);
end

% don't need last voltage for the optimizer (it will just be ground) so we toss it out
ElectrodeRadialPositions = ElectrodeRadialPositions(1:end-1);
ElectrodeVoltages = ElectrodeVoltages(1:end-1);
LB = LB(1:end-1);
UB = UB(1:end-1);

% Electrode X-R coordinates
ElectrodeCoordinates = ElectrodeRadialPositions'*[cos(angle) sin(angle)];


%%%%%%%%%%%%%%%%%%%%%%%
%    Set up domain    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    Domain is a skewed grid with non-constant spacing in x, so the
%    following parameters are highly custom to this specific domain
%
%
disp('Setting up domain')
nrInside = 30;            % number of nodes perpendicular to the entire beamline
nrExtra = 8;              % number of extra nodes that extend out in the perpendicular direction in the core
dx0 = 0.005*L;            % value of dx from x=0 to x=xBound
dxR = 0.02*L;             % value of dx at x=Lx
numElectrodes = numel(ElectrodeRadialPositions);
nr = nrInside+nrExtra;
xBound = cos(angle)*innerRadius;                                                % x-value of the inner radius of the device
rBound = sin(angle)*innerRadius;                                                % r-value of the inner radius of the device
zOffRadius = (Remanence>0)*innerRadius;                                         % radius at which charge of ions is "turned off" due to electron neutralization
LrOuter = Lx*tan(angle);                                                        % the maximum r-value at x=Lx
[xSpaceInner,~,~,~] = GraduallySpaced(0,xBound,dx0,dx0);                        % Equivalent to linspace since the 3rd and 4th arguments are the same
dx0 = mean(xSpaceInner(2:end)-xSpaceInner(1:end-1));                            % update our dx0 to the actual value
[xSpaceOuter,c,C,IndexFinderX] = GraduallySpaced(xBound,Lx,dx0,dxR);            % Spacing increases as you get farther away from fusion core
xSpace = unique([xSpaceInner;xSpaceOuter]);                                     % Complete list of x-coordinates of values
dx0real = xSpace(2)-xSpace(1);                                                  % Minimum spacing between x-points
nxBound = numel(xSpaceInner);                                                   % Number of nodes between x=0 anc x=xBound
nx = numel(xSpace);                                                             % Number of x-coordinates
nBulk = nrInside*nx;                                                            % Total number of nodes in the beamline (excluding the extra fusion core extension)
nTotal = nBulk + nxBound*nrExtra;                                               % Total number of nodes in the simulation (I think equivalent to nPM later)

rSpaceInner = linspace(0,rBound,nrInside); %line from [x=0,r=0] to [x=0,r=rBound] - this is rSpace in the beamline region
drInner = rSpaceInner(2)-rSpaceInner(1);    % dr at x<=xBound

rSpaceExtra = linspace(rBound+drInner,rBound+nrExtra*drInner,nrExtra); %line from [x=0,r=rBound] to [x=0,r=rBound] - this is rSpace in the extra fusion core extension
LrZero = rSpaceExtra(end);  % maximum r-value at x=0
Lr = max(LrOuter,LrZero); % maximum r-value overall

rSpaceOuter = linspace(0,LrOuter,nrInside); % List of r-coordinates at x=Lx
drOuter = rSpaceOuter(2)-rSpaceOuter(1);    % dr at x=Lx
drSlope = (drOuter-drInner)/(Lx-xBound);    % since dr varies linearly between x=xBound and x=Lx, this is its slope

% Set up simulation grid
xGrid = repmat(xSpace,1,nr);
rGridInner = zeros(nx,nrInside);
%rGridInner is rGrid along the beamline
for ii = 1:nx
    rGridInner(ii,:) = max(rSpaceInner,rSpaceInner+(rSpaceOuter-rSpaceInner)/(Lx-xBound)*(xSpace(ii)-xBound));
end
%rGridExtra is rGrid in the fusion region
rGridExtra = repmat(rSpaceExtra,nxBound,1);
rGrid = [rGridInner,[rGridExtra;NaN(nx-nxBound,nrExtra)]];  %use nans for the part of the grid not inluded in the simulation
xGrid(isnan(rGrid)) = NaN;  %nan out parts of the xGrid not included in the simulation

% Calculate time-step from smallest dx value
dx = min(min(xGrid(2:end,:)-xGrid(1:end-1,:)));
dt = CFL*dx/fusionVelocity;

% Calculate "normalized"/"scaled" Electrode positions and voltages for
% input into the optimizer, so that optimizer is dealing with things on the
% order of 1
normalizedElectrodeRadialPositions = (ElectrodeRadialPositions - ElectrodeRadialPositions(1))/(ElectrodeRadialPositions(end) - ElectrodeRadialPositions(1));
normalizedControlRangePositions = (xGrid(nxBound:end,nrInside).^2 + rGrid(nxBound:end,nrInside).^2).^.5;
normalizedControlRangePositions = (normalizedControlRangePositions - normalizedControlRangePositions(1))/(normalizedControlRangePositions(end)-normalizedControlRangePositions(1));
normalizedElectrodeVoltages = (ElectrodeVoltages - AllowableVoltageRange(1))/(AllowableVoltageRange(2) - AllowableVoltageRange(1));
normalizedLowerBound = (LB - AllowableVoltageRange(1))/(AllowableVoltageRange(2) - AllowableVoltageRange(1));
normalizedUpperBound = (UB - AllowableVoltageRange(1))/(AllowableVoltageRange(2) - AllowableVoltageRange(1));

% matrix used to hand over some domain information to the poisson matrix generator
% Inf - unknown voltage location
% 0   - ground
% 1   - Voltage control point
% NaN - Outside of domain
PoissonIndicator = inf(nx,nr);              
PoissonIndicator(end,1:nrInside) = 0;
PoissonIndicator(nxBound+1:nx,nrInside+1:end) = NaN;
PoissonIndicator(nxBound,nrInside:nr) = 1;
PoissonIndicator(1:nxBound,nr) = 1;
PoissonIndicator(nxBound+1:nx-1,nrInside) = 1;

% Calculate poisson matrix
[PM,PM_DIA,PM_OFF,EX,EX_DIA,EX_OFF,ER,ER_DIA,ER_OFF,bDirichlet] = CreatePoissonMatrix(xGrid,rGrid,PoissonIndicator);

nPM = size(PM,1);   % Number of nodes of unknown potential
nEM = size(EX,1);   % Total number of simulation nodes (nodes where we need to know the electric field
nBulkPM = (nx-1)*(nrInside-1);  % Number of nodes in the beam region (

% Create vector of indices to convert between the full grid indices and the
% vector indices of all simulation cells
vec2grid = zeros(nEM,1);
vec2grid(1:nBulk) = (1:nBulk)';
for j = 1:(nr-nrInside)
    vec2grid(nBulk+nxBound*(j-1)+(1:nxBound)) = nBulk+nx*(j-1)+(1:nxBound);
end 

% Create vector of indices to convert between the vector indices of all
% simulation cells and the vector indices of unknown potential cells
% - Did this section somewhat by trial and error.  in contains indices to map
% a vec value to a vecPM value, i.e. myVecPM = myVec(vec2vecPM).
% Also it works in the reverse, i.e. myVec(vec2vecPM) = myVecPM.
vec2vecPM = zeros(nPM,1);
for jj = 1:nrInside-1
    PMrange = (nx-1)*(jj-1) + (1:nx-1)';
    vecRange = nx*(jj-1) + (1:nx-1)';
    vec2vecPM(PMrange) = vecRange;
end
PMrange = nBulkPM + (1:nxBound-1)';
vecRange = nx*(nrInside-1) + (1:nxBound-1)';
vec2vecPM(PMrange) = vecRange;
for jj = (nrInside+1:nr-1)-nrInside+1
    PMrange = nBulkPM + (nxBound-1)*(jj-1) + (1:nxBound-1)';
    vecRange = nx*nrInside + nxBound*(jj-2) + (1:nxBound-1)';
    vec2vecPM(PMrange) = vecRange;
end

% Define specific function to pass to optimizer for finding initial
% voltages in MATLAB
Ev2PhiSpecific = @(EV) Ev2Phi(...
    EV,...
    normalizedElectrodeRadialPositions,...
    AllowableVoltageRange,...
    normalizedControlRangePositions,...
    PM,bDirichlet,vec2vecPM,vec2grid,...
    EX,ER,...
    nPM,nEM,nx,nr,nxBound,nrInside);

% Find initial voltages
[bDirichletControl, phiInitialVecPM, phiInitialVec,phiInitial,efxInitialVec,efrInitialVec] = Ev2PhiSpecific(normalizedElectrodeVoltages);

efxInitial = single(zeros(nx,nr));
efrInitial = single(zeros(nx,nr));
efxInitial(vec2grid) = efxInitialVec;
efrInitial(vec2grid) = efrInitialVec;

% Find the ion starting point so that the velocity in the core is the fusion velocity
disp('Estimating the period and turnaround point')
[bounceTime, beampoint, birthAccel] = PeriodEstimator(xGrid(:,1),fusionVelocity,q/m*efxInitialVec(1:nx),dt);


beampoint = single(beampoint);
minPhi = min(phiInitialVec);
maxPhi = max(phiInitialVec);
disp('Setting up Poisson Matrix')
% Calculate cell volumes (needed for calculating the density)
cellVolumes = FindCellVolumes(xGrid,rGrid);
cellVolumesVec = cellVolumes(vec2grid);
cellVolumesPM = cellVolumesVec(vec2vecPM);
bmaxVec = cellVolumesVec.^(1/3);
bmaxPM = single(cellVolumesPM.^(1/3));

%grid points where we don't calculate the ion density
zOff = (xGrid.^2 + rGrid.^2) > zOffRadius^2;
zOffVec = zOff(vec2grid);
zOffVecPM = single(zOffVec(vec2vecPM));

% Initialize density
rho = zeros(nx,nr);
rhoVec = zeros(nEM,1);

% Decompose PM_DIA etc. into matrices useful for CUDA
U1 = unique(PM_OFF(:,1));
U2 = unique(PM_OFF(:,2));
U3 = unique(PM_OFF(:,3));
U4 = unique(PM_OFF(:,4));
U5 = unique(PM_OFF(:,5));
U6 = unique(PM_OFF(:,6));
U7 = unique(PM_OFF(:,7));

X1 = unique(EX_OFF(:,1));
X2 = unique(EX_OFF(:,2));
X3 = unique(EX_OFF(:,3));
X4 = unique(EX_OFF(:,4));
X5 = unique(EX_OFF(:,5));
R1 = unique(ER_OFF(:,1));
R2 = unique(ER_OFF(:,2));
R3 = unique(ER_OFF(:,3));

U1 = U1(U1~=0);
U2 = U2(U2~=0);
U3 = U3(U3~=0);
U4 = U4(U4~=0);
U5 = U5(U5~=0);
U6 = U6(U6~=0);
U7 = U7(U7~=0);

X1 = X1(X1~=0);
X2 = X2(X2~=0);
X3 = X3(X3~=0);
X4 = X4(X4~=0);
X5 = X5(X5~=0);
R1 = R1(R1~=0);
R2 = R2(R2~=0);
R3 = R3(R3~=0);

cellResidents = int32(zeros(nPM*maxParticlesPerCell,1));   % array that holds the indices of all particles that reside in that cell
partNs = single(zeros(nPM*maxParticlesPerCell,1));   % array that holds the indices of all particles that reside in that cell
partAs = single(zeros(nPM*maxParticlesPerCell,1));   % array that holds the indices of all particles that reside in that cell

cellOccupants = int32(zeros(nPM,1));                       % array that holds the number of occupants in each cell

fusionRate = single(zeros(nPM,1));                          % array for stores results of local fusion rate calculations

w = numReal/np;

numIt = ceil(bounceTime/dt)+1;
dt = single(bounceTime/numIt);


maxDensity = 2*numReal/(2*pi*sigX*sigR*sigR);

px = abs(normrnd(0,sigX,np,1));
py = abs(normrnd(0,sigR,np,1));
pz = abs(normrnd(0,sigR,np,1));

pr = sqrt(py.^2 + pz.^2);

randLeftRight = [-ones(np/2,1);ones(np/2,1)];  % this is done to ensure equal amounts of positive and negative moving particles
randLeftRight = randLeftRight(randperm(np));

vx = fusionVelocity*(randLeftRight + normrnd(0,sigVX,np,1));
vy = fusionVelocity*normrnd(0,sigVR,np,1);
vz = fusionVelocity*normrnd(0,sigVR,np,1);

%calculate r and theta velocities (from http://www.maths.usyd.edu.au/u/MOW/vectors/vectors-10/v-10-2.html)
vrVec = bsxfun(@times,[py pz],(py.*vy + pz.*vz)./(py.^2 + pz.^2));
vtVec = [vy vz] - vrVec;
vr = (2*round(rand(np,1))-1).*sum(vrVec.^2,2).^.5;
vt = (2*round(rand(np,1))-1).*sum(vtVec.^2,2).^.5;
vr(isnan(vr)) = 0;
vt(isnan(vt)) = 0;

px = single(px);
pr = single(pr);
vx = single(vx);
vr = single(vr);
vt = single(vt);
KE = single(zeros(np,1));
PE = single(zeros(np,1));

creationRate = single(0);
nt = int32(numIt);

minDensity = .1*w/max(max(cellVolumes));

%%
MexIEC_Name = CompileIfNecessary(mexFile); 
MexIEC_Handle = str2func(MexIEC_Name);

efxVec = single(efxInitialVec);
efrVec = single(efrInitialVec);
efxGPU = single(efxInitial);
efrGPU = single(efrInitial);
% rhoVec = single(zeros(nEM,1));

creationRate = single(creationRate);

dt = single(dt);

nxBound = int32(nxBound);
xBound = single(xBound);
OneOverC = single(1/c);
c = single(c);
C = single(C);
drSlope = single(drSlope);
dx = single(dx0);
dr = single(drInner);
nx = int32(nx);
nr = int32(nr);
nrBound = int32(nrInside);
nInside = int32(nBulk);
phicpu = 1*phiInitial;
NumTimers = 16;
CudaFunctionTimers = single(zeros(NumTimers,1));

PM_OFF_1A = int32(U1(1));
PM_OFF_1B = int32(U1(2));
PM_OFF_1T = int32(find(PM_OFF(:,1)==PM_OFF_1B,1))-1;
PM_OFF_2A = int32(U2(1));
PM_OFF_2B = int32(U2(2));
PM_OFF_2T = int32(find(PM_OFF(:,2)==PM_OFF_2B,1))-1;
PM_OFF_6A = int32(U6(2));
PM_OFF_6B = int32(U6(1));
PM_OFF_6T = int32(find(PM_OFF(:,6)==PM_OFF_6B,1))-1;
PM_OFF_7A = int32(U7(2));
PM_OFF_7B = int32(U7(1));
PM_OFF_7T = int32(find(PM_OFF(:,7)==PM_OFF_7B,1))-1;

EX_OFF_1A = int32(X1(1));
EX_OFF_1B = int32(X1(2));
EX_OFF_1T = int32(find(EX_OFF(:,1)==EX_OFF_1B,1))-1;
EX_OFF_5A = int32(X5(2));
EX_OFF_5B = int32(X5(1));
EX_OFF_5T = int32(find(EX_OFF(:,5)==EX_OFF_5B,1))-1;

ER_OFF_1A = int32(R1(1));
ER_OFF_1B = int32(R1(2));
ER_OFF_1T = int32(find(ER_OFF(:,1)==ER_OFF_1B,1))-1;
ER_OFF_3A = int32(R3(2));
ER_OFF_3B = int32(R3(1));
ER_OFF_3T = int32(find(ER_OFF(:,3)==ER_OFF_3B,1))-1;

PM_OFF_1A = int32(abs(PM_OFF_1A));
PM_OFF_1B = int32(abs(PM_OFF_1B));
PM_OFF_1T = int32(abs(PM_OFF_1T));
PM_OFF_2A = int32(abs(PM_OFF_2A));
PM_OFF_2B = int32(abs(PM_OFF_2B));
PM_OFF_2T = int32(abs(PM_OFF_2T));
PM_OFF_6A = int32(abs(PM_OFF_6A));
PM_OFF_6B = int32(abs(PM_OFF_6B));
PM_OFF_6T = int32(abs(PM_OFF_6T));
PM_OFF_7A = int32(abs(PM_OFF_7A));
PM_OFF_7B = int32(abs(PM_OFF_7B));
PM_OFF_7T = int32(abs(PM_OFF_7T));

EX_OFF_1A = int32(abs(EX_OFF_1A));
EX_OFF_1B = int32(abs(EX_OFF_1B));
EX_OFF_1T = int32(abs(EX_OFF_1T));
EX_OFF_5A = int32(abs(EX_OFF_5A));
EX_OFF_5B = int32(abs(EX_OFF_5B));
EX_OFF_5T = int32(abs(EX_OFF_5T));

ER_OFF_1A = int32(abs(ER_OFF_1A));
ER_OFF_1B = int32(abs(ER_OFF_1B));
ER_OFF_1T = int32(abs(ER_OFF_1T));
ER_OFF_3A = int32(abs(ER_OFF_3A));
ER_OFF_3B = int32(abs(ER_OFF_3B));
ER_OFF_3T = int32(abs(ER_OFF_3T));


PM1 = single(-PM_DIA(:,1));
PM2 = single(-PM_DIA(:,2));
PM3 = single(-PM_DIA(:,3));
PM4 = single(1./PM_DIA(:,4));
PM5 = single(-PM_DIA(:,5));
PM6 = single(-PM_DIA(:,6));
PM7 = single(-PM_DIA(:,7));

EX1 = single(-EX_DIA(:,1));
EX2 = single(-EX_DIA(:,2));
EX3 = single(-EX_DIA(:,3));
EX4 = single(-EX_DIA(:,4));
EX5 = single(-EX_DIA(:,5));
ER1 = single(-ER_DIA(:,1));
ER2 = single(-ER_DIA(:,2));
ER3 = single(-ER_DIA(:,3));

fusionVelocity = single(fusionVelocity);
cellVolumesPM = single(cellVolumesPM);
phiVecPM = single(phiInitialVecPM);
phiVecGPU = single(phiInitialVec);
phiVec4GPU = single(phiInitialVec);
phiGPU = single(phiInitial);
rhoVecPM = single(zeros(nPM,1));
cellVolumesPM = single(cellVolumesPM);

macroWeight = single(w);
W = single(1);
LxZero = single(Lx);
LrZero = single(LrZero);
tanangle = single(tan(angle));
qom = single(q/m);
creationRadius = single(rBound/2);
np = int32([np 0]);
npR = single([0 0]);
L = single(L);
a_coeff = single(a_coeff);
numPeriods = int32(1);

%% Magnetic field

if Remanence > 0
    load('TrunciData.mat')
    dpmMag = 0.01;
    inner2outer = innerRadius/outerRadius;
    [ptsDip, dip, dipVols, streamLineStartPoints, fileStringD] = CreateMagneticDipoles(dpmMag,MagnetThickness,inner2outer,Fpent,Fhex,V,PH,PV,PF,distBetweenPentAndVert);
    pts = [xGrid(:), rGrid(:), zeros(size(xGrid(:)))];
    pts(isnan(xGrid(:)),:) = [];
    stringID = ['n=',num2str(size(pts,1)),'i=',num2str(inner2outer),'o=',num2str(1),'X=',num2str(nx),'R=',num2str(nr),'x=',num2str(nxBound),'r=',num2str(nrInside),'t=',num2str(angle),'d=',num2str(dpmMag)];
        
    numTheta = 40;   % average the magnetic field around the axis
    dtheta = 2*pi/numTheta;
    thetaSpace = 0:dtheta:(2*pi-dtheta);
    pts = zeros(numel(xGrid)*numTheta,3);
    for i = 1:numTheta
        range = (1:numel(xGrid)) + (i-1)*numel(xGrid);
        pts(range,:) = [-rGrid(:)*cos(thetaSpace(i)), -xGrid(:), -rGrid(:)*sin(thetaSpace(i))];
    end
    
    pts = [-rGrid(:), -xGrid(:), zeros(size(xGrid(:)))];
    
    pts(isnan(xGrid(:)),:) = [];
    
    % Calclulate (or load) magnetic field
    bfxr = single((1/(4*pi))*Remanence*MagFromDipoles(ptsDip,dip,(3*dipVols/(4*pi)).^(1/3),pts,stringID));
    
    bfx = zeros(nx,nr);
    bfr = zeros(nx,nr);
    bfm = zeros(nx,nr);
    bfxVec = bfxr(:,2);
    bfrVec = bfxr(:,1);
    bfmVec = sqrt(bfxVec.^2 + bfrVec.^2);
    bfx(vec2grid) = bfxVec;
    bfr(vec2grid) = bfrVec;
    bfm(vec2grid) = bfmVec;
    bfmPlot = bfm;
    bfmPlot(bfm>.25) = .25;
else
    bfxVec = single(zeros(size(efxVec)));
    bfrVec = bfxVec;
end

%%
disp('Packing structs for optimizer...')
periodEst = bounceTime;
periods = [];
LL = single(1);
SimulationVariables = v2struct(np,npR,periods,px,pr,vx,vr,vt,KE,PE,efxVec,efrVec,bfxVec,bfrVec,rhoVecPM,phiVecPM,phiVec4GPU);
SimulationConstants = v2struct(dx,dr,nx,nr,nxBound,nrBound,nInside,xBound,zOffRadius,c,C,drSlope,dt,nt,numPeriods,...
                               creationRate,beampoint,creationRadius,fusionVelocity,LxZero,LrZero,tanangle,q,m,qom,...
                               PM1,PM2,PM3,PM4,PM5,PM6,PM7,...
                               PM_OFF_1A,PM_OFF_1B,PM_OFF_1T,...
                               PM_OFF_2A,PM_OFF_2B,PM_OFF_2T,...
                               PM_OFF_6A,PM_OFF_6B,PM_OFF_6T,...
                               PM_OFF_7A,PM_OFF_7B,PM_OFF_7T,...
                               EX1,EX2,EX3,EX4,EX5,ER1,ER2,ER3,...
                               EX_OFF_1A,EX_OFF_1B,EX_OFF_1T,...
                               EX_OFF_5A,EX_OFF_5B,EX_OFF_5T,...
                               ER_OFF_1A,ER_OFF_1B,ER_OFF_1T,...
                               ER_OFF_3A,ER_OFF_3B,ER_OFF_3T,...
                               cellVolumesPM,bmaxPM,a_coeff,macroWeight,W,itPM,L,bDirichletControl,zOffVecPM,...
                               maxParticlesPerCell,cellResidents,partNs,partAs,cellOccupants,cudaRandSeed,FixEnergy,loadPrevMode,energyBump);

SimulationPerformance = v2struct(CudaFunctionTimers,fusionRate);

OptimizerNecessities = v2struct(...
    fusionVelocity,...
    fusionEnergyReleased,...
    normalizedElectrodeVoltages,...
    normalizedLowerBound,...
    normalizedUpperBound,...
    ElectrodeRadialPositions,...
    Ev2PhiSpecific,...
    xGrid,rGrid,L,AllowableVoltageRange,vec2vecPM,vec2grid,cellVolumesPM,...
    periodEst,zOffRadius,macroWeight,...
    OptimizeMode,filename,timename,fusionKE,loadPrevMode,...
    MexIEC_Handle);

np = double(np(1)); %make a MATLAB version of scalar np  for plotter (previous one had np and npKill as int32 array)

PlotterNecessities = v2struct(xGrid,rGrid,L,LB,UB,AllowableVoltageRange,ElectrodeCoordinates,ElectrodeRadialPositions,vec2vecPM,vec2grid,cellVolumesPM,...
    zOffRadius,macroWeight,np,fusionVelocity,fusionVoltage);

% write input variables to notes file
notes = []; 
notes = userNotes;
notes = [notes,evalc('gitInfo')];
notes = [notes,evalc('SimulationVariables')];
notes = [notes,evalc('SimulationConstants')];
notes = [notes,evalc('SimulationPerformance')];
notes = [notes,evalc('OptimizerNecessities')];

disp('launching optimizer function...')
OptimizerWrapper(...
    SimulationVariables,...
    SimulationConstants,...
    SimulationPerformance,...
    OptimizerNecessities,...
    PlotterNecessities,...
    notes);