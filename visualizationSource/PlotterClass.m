classdef PlotterClass < handle
   
    properties (Access=public)
        Fig
        
        TopRowHandle
        PhiHandle
        RhoHandle
        CostHandle
        ParticleHandle
        ElectrodeHandle
        TemperatureHandle
        SimTextHandle
        EnergyHandle
        
        np
        vec2grid
        vec2vecPM
        cellVolumesPM
        xGrid
        rGrid
        L
        zOffRadius
        partColors
        AllowableVoltageRange
        fusionVelocity
        UB
        LB
        ElectrodeRadialPositions
        ElectrodeCoordinates
        macroWeight
        
    end
    
    properties (Constant)
        colorRed = [.85 0 0];
        colorGray = [.4 .4 .4];
        colorPurple = [.6 .3 .9];
        colorGreen = [.3 .7 .15];
        colorBlue = [0 0.4470 0.7410];
        colorOrange = [0.85 0.325 0.098];
        colorLightGray = [.6 .6 .6];
        colorLightRed = [.85 .7 .7];
        colorFaintRed = [.75 .3 .2];        
        size1 = 13;
        size2 = 12;
        size3 = 14;
        size4 = 15;
        size5 = 11;
        size6 = 18;
        size7 = 12;
        size8 = 11;
        size9 = 9;
        size10 = 10;
        size11 = 8;
    end
    
    methods
        function this = PlotterClass(PN)
            
            this.Fig = figure(1);

            this.np = PN.np;
            set(this.Fig,'outerposition', [ 575  140      1239    906],'color','w')
            
            this.vec2grid = PN.vec2grid;
            this.vec2vecPM = PN.vec2vecPM;
            this.cellVolumesPM = PN.cellVolumesPM;
            this.xGrid = PN.xGrid;
            this.rGrid = PN.rGrid;
            this.AllowableVoltageRange = PN.AllowableVoltageRange;
            this.zOffRadius = PN.zOffRadius;
            this.partColors = [.65*[1 1 1]; this.colorBlue; this.colorOrange];
            this.L = PN.L;
            this.fusionVelocity = PN.fusionVelocity;
            this.UB = PN.UB;
            this.LB = PN.LB;
            this.ElectrodeRadialPositions = PN.ElectrodeRadialPositions;
            this.ElectrodeCoordinates = PN.ElectrodeCoordinates;
            this.macroWeight = PN.macroWeight;

        end
        
        function SetVisibility(this,visibility)
            if visibility
                this.Fig.Visible = 'on';
            else
                this.Fig.Visible = 'off';
            end
        end
        
        function Delete(this)
            close(this.Fig);
            delete(this);
        end
        
        function CreateTopRow(this,Optimizer)
            if Optimizer
                this.TopRowHandle = TopRowExtend(this.size6,this.size7,this.colorRed,this.colorBlue);
            else
                this.TopRowHandle = TopRow(this.size6,this.size7,this.colorRed,this.colorBlue);
            end
        end
        
        function CreateCostPlot(this)
            this.CostHandle = CostPlot(this.colorPurple,this.colorGreen,this.size5);
        end
        
        function CreatePhiPlot(this)
            this.PhiHandle = PhiPlot(this,this.xGrid,this.rGrid,this.ElectrodeCoordinates,this.AllowableVoltageRange,this.zOffRadius,this.size1,this.size2,this.colorRed);
        end
        
        function CreateRhoPlot(this)
            this.RhoHandle = RhoPlot(this,this.xGrid,this.rGrid,this.size1,this.size2);
        end
        
        function CreateParticlePlots(this)
            this.ParticleHandle = ScatterParticles(this,this.partColors,this.size3,this.size4,this.size8);
        end
        
        function CreateElectrodePlot(this)
            this.ElectrodeHandle = ElectrodePlot(this,this.colorLightRed,this.colorRed,this.size5);
        end
        
        function CreateTemperaturePlot(this)
            this.TemperatureHandle = TemperaturePlot(this,this.colorRed,this.colorBlue,this.colorGreen,this.size5);
        end
        
        function CreateEnergyPlot(this)
            this.EnergyHandle = EnergyPlot(this.colorRed,this.colorBlue,this.size5);
        end
        
        function CreateSimText(this)
            this.SimTextHandle = SimText(this,this.size10,this.colorGray);
        end

    end
            
    
end
    
