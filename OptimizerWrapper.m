% Gateway function for both optimizing the IEC and running a long-timescale
% IEC simulation

function OptimizerWrapper(...
        SV,...
        SC,...
        SP,...
        ON,...
        PN,... // used as output information
        notes)
    
    init = int32(1);
    ticTotal = tic;
    numPeriodsCompleted=1;
    
    pxState = 1*SV.px;
    prState = 1*SV.pr;
    vxState = 1*SV.vx;
    vrState = 1*SV.vr;
    vtState = 1*SV.vt;
    
    npInitial = double(SV.np(1));
    
    outputPath = ['outputFiles'];
    if ~exist(outputPath,'dir')
        mkdir(outputPath)
        disp(['creating a folder ''\outputFiles'' in ',pwd])
    end

    
    filenameGlobal = [ON.filename(end-2:end),ON.timename];
    masterDirectory = [outputPath,'\',filenameGlobal,'\'];
%     fileNameAndPathGlobal = [masterDirectory,'\',filenameGlobal];
    mkdir(masterDirectory)
    frameDataDirectory = [masterDirectory,'data\'];
    mkdir(frameDataDirectory)
    
    plotterFileName = ['RUN_ME_',ON.filename(end-2:end),'.m'];

    if ON.OptimizeMode
        plotterSourceName = 'visualizationSource\OptPlotter';
    else
        plotterSourceName = 'visualizationSource\RunPlotter';
    end
    
    % got the following from https://stackoverflow.com/questions/13846491/how-to-textscan-to-read-all-the-lines-in-a-file
    dlmwrite([masterDirectory,'\',plotterFileName], ['clear all; close all; clc; fileName = ''',filenameGlobal,''';'], 'delimiter', '');
    dlmwrite([masterDirectory,'\',plotterFileName], fileread([plotterSourceName,'.m']), '-append', 'delimiter', '');
%     copyfile('PlotterClassMaster.m',[masterDirectory,'\PlotterClass.m']);
    % write notes to text file
    fid = fopen([masterDirectory,'\notes',filenameGlobal,'.txt'],'w');
    fprintf(fid, notes);
    fclose(fid);
    ntOrig = int32(SC.numPeriods)*SC.nt;   % we change nt to a lower value when recording video, so we need to save it here

    dt = SC.dt;

    save([frameDataDirectory,filenameGlobal,'_initialData'],'PN')
    drawnow

    frameNumber = 1;
    costTracker = [];
    costMode = [];
    numIterations = 0;
    frameSkip = 1;
  
    bestElectrodeVoltagesToOptimize = ON.normalizedElectrodeVoltages;
    periodLinesXdata = [];
if ON.OptimizeMode
    
    while true  % infinite loop, just keep optimizing until stopped by user
        
        numOptFrames = 0;
        addToSkip = 0;
        
        if ON.loadPrevMode
            iterGlobal = 85;
            iterLocal = 70;
            iterLocalCatch = 70;
        else
            iterGlobal = 50;
            iterLocal = 40;
            iterLocalCatch = 70;
        end
            
        optimizerExists = false;
        
        if ON.loadPrevMode
            disp('Loading best particle state')
            SV.px = 1*pxState;
            SV.pr = 1*prState;
            SV.vx = 1*vxState;
            SV.vr = 1*vrState;
            SV.vt = 1*vtState;
        end
        
        
        
        try
            
            globalOpts = saoptimset('simulannealbnd');
            globalOpts.InitialTemperature = .5*ones(size(ON.normalizedElectrodeVoltages));
            globalOpts.ReannealInterval = 43;
            globalOpts.MaxFunEvals = iterGlobal;
            globalOpts.Display = 'iter';
            globalOpts.TolFun = 1e-4;
%             keyboard
            globalOpts.AnnealingFcn = @annealingFunction;
            costThisPeriod = [];
            optimizerText = '{\color[rgb]{0.60 0.30 0.90}Simulated Annealing (Global)}';
            LocalIt = 1;
            Global = true;
            ElectrodeVoltagesOptimizerInput = bestElectrodeVoltagesToOptimize;
            disp(['SIMULATED ANNEALING Optimizer for ',num2str(SC.numPeriods),' periods'])
            % Start Simulated Annealing optimizer
            inputScalingOffset = 0*ElectrodeVoltagesOptimizerInput;
            [OptimizedElectrodeVoltages,fval,exitFlag,output] = simulannealbnd(...
                @CostFunction,...
                ElectrodeVoltagesOptimizerInput,...
                ON.normalizedLowerBound,...
                ON.normalizedUpperBound,...
                globalOpts);
            optimizerExists = true;
        catch
            warning('Either the cost function failed, or you do not have the global optimization toolbox, only the local optimizer will run')
            iterLocal = iterLocalCatch;
        end
        
        try
            localStepScale = 1./(.5+sqrt(double(SC.numPeriods)));
            localOpts = optimset('fminsearch');
            localOpts.MaxFunEvals = iterLocal;
            localOpts.Display = 'iter';
            optimizerText = '{\color[rgb]{0.30 0.70 0.15}Nelder-Mead Simplex (Local)}';
            ElectrodeVoltagesOptimizerInput = bestElectrodeVoltagesToOptimize;
            LocalIt = 1;
            Global = false;
            disp(['Nelder-Mead Simplex Optimizer for ',num2str(SC.numPeriods),' periods'])
            
            inputScalingOffset = 0*ElectrodeVoltagesOptimizerInput;%(1-localStepScale)*ElectrodeVoltagesOptimizerInput;
            % Start Nelder-Mead
            [OptimizedElectrodeVoltages,fval,exitFlag,output] = fminsearchbndDrew(...
                @CostFunction,...
                ElectrodeVoltagesOptimizerInput-inputScalingOffset,...
                ON.normalizedLowerBound,...
                ON.normalizedUpperBound,...
                localStepScale,...
                localOpts);
        catch
            if optimizerExists
                warning('Either the cost function failed, or you do not have the optimization toolbox, only the global optimizer will run')
            else
                error('You do not have any of the required optimization functions')
            end
        end
        ElectrodeVoltagesRescaled = OptimizedElectrodeVoltages*(ON.AllowableVoltageRange(2)-ON.AllowableVoltageRange(1)) + ON.AllowableVoltageRange(1);
        periodLinesXdata(SC.numPeriods,:) = numIterations;
        filename = [frameDataDirectory,'FullDump_',num2str(SC.numPeriods,'%03i')];
        save(filename)
        SC.numPeriods = SC.numPeriods+round(SC.numPeriods^2*1.3e-3)+1;
        numPeriodsCompleted = SC.numPeriods;
    end
else
    
    disp('Long timescale simulation')
    optimizerText = [];
    LocalIt = [];
    SC.numPeriods = 1000;
    [SC.bDirichletControl, SVcopy.phiVecPM, SVcopy.phiVec4GPU,~,SVcopy.efxVec,SVcopy.efrVec]...
        = ON.Ev2PhiSpecific(ON.normalizedElectrodeVoltages);
     
    NumFrames = SC.numPeriods*4*30; % constant number of frames per period
    ntTotalEst = ceil(ON.periodEst*double(SC.numPeriods)/SC.dt);  % upper estimate total number of timesteps needed for the video
    ntTotalSafe = ceil(1.2*ntTotalEst);
    SC.nt = int32(ceil(ntTotalEst/double(NumFrames)));  % number of time-steps per video frame
    
    SVcopy = CopyStruct(SV);
    SVcopy.periods = single(zeros(SC.numPeriods+1,1));
    MexInputs = [struct2cell(SVcopy);struct2cell(SC);struct2cell(SP)];
    
    timeVector = double(1:ntTotalSafe)'*double(SC.dt);
    fusionOutput = single(-ones(ntTotalSafe,1));
    energyLoss = single(-ones(ntTotalSafe,1));
    
    IT = 1;
    
    init = int32(1);
    time = 0;
    numPeriodsCompleted = sum(SVcopy.periods>0);
    while numPeriodsCompleted < SC.numPeriods % sum(SV.periods>0) is how many periods have been completed by the CUDA simulation
         disp(['Starting frame = ',num2str(frameNumber),', nt = ',num2str(IT*double(SC.nt)),])
        ntRange = (1:SC.nt) + (IT-1)*SC.nt;

        tic
        
        [energyLoss(ntRange), fusionOutput(ntRange),timeTemp,tempX,tempR,tempT] = ON.MexIEC_Handle(init,MexInputs{:});
        
        np = 1*SVcopy.np(1);
        npKill = 1*SVcopy.np(2);
        [costX, costV, costLoss] = calculateCostFunction(SVcopy.px(1:np),SVcopy.pr(1:np),SVcopy.vx(1:np),SVcopy.vr(1:np),SVcopy.vt(1:np),npKill);
        % For checking conservation of Energy
%         SV.np(1)
        KE(IT) = sum(SVcopy.KE)/(ON.macroWeight*ON.fusionKE*double(SV.np(1)));
        PE(IT) = sum(SVcopy.PE)/(ON.macroWeight*ON.fusionKE*double(SV.np(1)));
%         if any(energyLoss(ntRange))
%             keyboard
%         end
           
        Eloss(IT) = sum(energyLoss(ntRange))/(ON.macroWeight*ON.fusionKE*double(SV.np(1)));
        EnergyTime(IT) = double(SC.nt)*IT*double(SC.dt);
%         Eloss(isnan(Eloss)) = 0;
        time = time + double(SC.nt)*double(SC.dt);
        
        disp(['nt = ',num2str(IT*double(SC.nt)),', t = ',num2str(time)]);
        UpdatePlot(SVcopy.px, SVcopy.pr, SVcopy.vx, SVcopy.vr, SVcopy.vt, SVcopy.phiVec4GPU, SVcopy.rhoVecPM,SVcopy.np(1),...
            ON.normalizedElectrodeVoltages,bestElectrodeVoltagesToOptimize,SC.numPeriods,SVcopy.periods,timeVector,EnergyTime,energyLoss,fusionOutput,timeTemp,tempX,tempR,tempT,npKill,KE,PE,Eloss);
        init = int32(0);
        IT  = IT+1;
        
    end

end

    function cost = CostFunction(ElectrodeVoltagesToOptimize)
        
        numIterations = numIterations+1;
        
        disp(['Optimization Iteration: ',optimizerText(29:end)])
        tic

        SV.periods = single(zeros(SC.numPeriods,1));
        SVcopy = CopyStruct(SV);
       
        [SC.bDirichletControl, SVcopy.phiVecPM, SVcopy.phiVec4GPU,~,SVcopy.efxVec,SVcopy.efrVec]...
            = ON.Ev2PhiSpecific(ElectrodeVoltagesToOptimize);
        
        SC.nt = round (1.5*ntOrig*SC.numPeriods);
        timeVector = double(1:SC.nt)'*double(SC.dt);
        MexInputs = [struct2cell(SVcopy);struct2cell(SC);struct2cell(SP)];
        
        [energyLoss, fusionOutput, timeTemp, tempX, tempR, tempT] = ON.MexIEC_Handle(int32(1),MexInputs{:});
        
        np = SVcopy.np(1);
        disp(['moved in ',num2str(toc)])
         
        np = SVcopy.np(1);
        npKill = SVcopy.np(2);
        [costX, costV, costLoss] = calculateCostFunction(SVcopy.px(1:np),SVcopy.pr(1:np),SVcopy.vx(1:np),SVcopy.vr(1:np),SVcopy.vt(1:np),npKill);
        cost = min(1,double(costX + costV + costLoss));
        costTracker(numIterations) = cost;
        costMode(numIterations) = Global;
        
        costThisPeriod(numOptFrames+1) = cost;
        [~,best] = min(costThisPeriod);
        
        if best == numOptFrames+1
            disp(['Writing BEST COST and Particle State for period ',numPeriodsCompleted])
            bestElectrodeVoltagesToOptimize = ElectrodeVoltagesToOptimize;
            pxState(1:np) = 1*SVcopy.px(1:np);  % the 1:np probably isnt necessary, but it doesn't hurt either
            prState(1:np) = 1*SVcopy.pr(1:np);
            vxState(1:np) = 1*SVcopy.vx(1:np);
            vrState(1:np) = 1*SVcopy.vr(1:np);
            vtState(1:np) = 1*SVcopy.vt(1:np);
            if ~exist([masterDirectory,'optimizedVoltages'],'dir')
                mkdir([masterDirectory,'optimizedVoltages']);
            end
            save([masterDirectory,'optimizedVoltages\bestVoltgagesPeriod',num2str(numPeriodsCompleted)],'bestElectrodeVoltagesToOptimize','numPeriodsCompleted')
        end
        
        if mod(numOptFrames,frameSkip) == 0
            UpdatePlot(SVcopy.px, SVcopy.pr, SVcopy.vx, SVcopy.vr, SVcopy.vt, SVcopy.phiVec4GPU, SVcopy.rhoVecPM,SVcopy.np(1),...
                ElectrodeVoltagesToOptimize,bestElectrodeVoltagesToOptimize,SC.numPeriods,SVcopy.periods,timeVector,[],energyLoss,fusionOutput,timeTemp,tempX,tempR,tempT,npKill,[],[],[]);
        end
        numOptFrames = numOptFrames+1;
        LocalIt = LocalIt + 1;
    end
    
    function UpdatePlot(px, pr, vx, vr, vt, phiVec,rhoVecPM,np,ElectrodeVoltagesToPlot,bestElectrodeVoltagesToPlot,numPeriods,periods,timeVector,EnergyTime,energyLoss,fusionOutput,timeTemp,tempX,tempR,tempT,npKill,KE,PE,Eloss)
        disp(['Writing frame ',num2str(frameNumber)]);
        filename = [filenameGlobal,'_',num2str(frameNumber,'%06i')];
        currentCompTime = toc(ticTotal);
        saveStruct = v2struct(px,pr,vx,vr,vt,phiVec,rhoVecPM,np,ElectrodeVoltagesToPlot,bestElectrodeVoltagesToPlot,periodLinesXdata,numPeriods,numPeriodsCompleted,periods,timeVector,EnergyTime,energyLoss,fusionOutput,costX,costV,costLoss,costTracker,costMode,currentCompTime,optimizerText,LocalIt,timeTemp,tempX,tempR,tempT,npKill,KE,PE,Eloss);
        save([frameDataDirectory,filename],'saveStruct');
        frameNumber = frameNumber + 1;
    end

    function [CX, CV, CL] = calculateCostFunction(px,pr,vx,vr,vt,npKill)
        CX = single(10)*norm(px.^2 + pr.^2)/single(npInitial);
        CV = single(10)*norm(((abs(vx)-ON.fusionVelocity).^2 + vr.^2 + vt.^2)/ON.fusionVelocity^2)/single(npInitial);
        CL = single(npKill)/single(npInitial);
        CX(isnan(CX)) = 0;
        CV(isnan(CV)) = 0;
        CL(isnan(CL)) = 1;
        disp(['Cost = ',num2str(CX+CV+CL)])
    end
    
end