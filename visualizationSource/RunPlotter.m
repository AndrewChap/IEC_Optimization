
% IEC_CUDA_Optimizer writes the value of 'filename' to this file

addpath('..\..\commonIEC')        % functions specific to this simulation
addpath('..\..\commonCustom')     % functions I wrote
addpath('..\..\commonFEX')        % functions downloaded from MATLAB file exchange
addpath('..\..\commonFiles')      % non-code files
addpath('..\..\visualizationSource')     % plotter class

periodStateFolder = 'periodStateFolder\';
mkdir(periodStateFolder)
videoRecord = false;
visibility = true;
deleteUsed = false;
disp('loading figure...')
dataPath = 'data\';
load([dataPath,fileName,'_initialData.mat']);
disp('figure loading complete...')
fileExists = true;
frameNumber = 1;
oldFrameName = '';
oldNumPeriods = 0;
if videoRecord
    videoRecorder = VideoRecorder(fileName,pwd,'-fullTitle');
    videoRecorder.zoom = '-m1.4';
    videoRecorder.frameRate = 30;
end

PlotterObject = PlotterClass(PN);
SetVisibility(PlotterObject,visibility);
CreateTopRow(PlotterObject,false)
% CreateCostPlot(PlotterObject)
CreatePhiPlot(PlotterObject)
CreateRhoPlot(PlotterObject)
CreateParticlePlots(PlotterObject)
CreateElectrodePlot(PlotterObject)
CreateTemperaturePlot(PlotterObject)
CreateEnergyPlot(PlotterObject)
CreateSimText(PlotterObject)

%%

while fileExists
    frameName = ['data\',fileName,num2str(frameNumber,'_%06i'),'.mat'];
    try
        load(frameName);
    catch
        disp('file doesn''t exist yet. waiting...')
        pause(2)
        continue
    end
    disp(['Frame # ',num2str(frameNumber)])
    v2struct(saveStruct);
    numSteps = sum(energyLoss>=0);
    numPeriods = sum(periods>0);
    if numPeriods > oldNumPeriods
        save([periodStateFolder,'PeriodState',num2str(numPeriods,'%04i')],'px','pr','vx','vr','vt')
    end
    UpdatePhiPlot(PlotterObject.PhiHandle,phiVec);
    UpdateRhoPlot(PlotterObject.RhoHandle,rhoVecPM);
    UpdateParticles(PlotterObject.ParticleHandle,np,px,pr,vx,vr,vt);
    UpdateElectrodePlot(PlotterObject.ElectrodeHandle,bestElectrodeVoltagesToPlot,ElectrodeVoltagesToPlot,PN.fusionVoltage)
    UpdateTemperaturePlot(PlotterObject.TemperatureHandle,timeTemp,tempX,tempR,tempT);
    UpdateEnergyPlot(PlotterObject.EnergyHandle,EnergyTime(1:frameNumber),KE(1:frameNumber),PE(1:frameNumber),Eloss(1:frameNumber));
    UpdateTopRow(PlotterObject.TopRowHandle,numPeriods,currentCompTime);
    UpdateSimText(PlotterObject.SimTextHandle,timeVector(numSteps),npKill,'',np);

    drawnow

    if videoRecord
        videoWrite(videoRecorder,frameNumber);
    end
    if deleteUsed && ~isempty(oldFrameName)
        delete(oldFrameName);  % only delete the previous one so that we never delete our last frame!
    end
    oldFrameName = frameName;
    frameNumber = frameNumber+1;
end