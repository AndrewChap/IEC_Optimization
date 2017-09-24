classdef TopRowExtend
    % Writes the Period and time information to top row of output window

    properties
        textPeriod
        textTime
        textCX
        textCV
        textCL
        
        maxPeriods = 0;
    end
    
    methods
        function this = TopRowExtend(size1,size2,color1,color2)
            
            % construct superclass PlotterTopRow
%             this@PlotterTopRow(size1,size2,color1,color2);
            textPeriod = annotation('textbox','interpreter','latex');
            textPeriod.Position = [.11 .98 0 0];
            textPeriod.String = ['\# Periods = 1'];
            textPeriod.EdgeColor = 'none';
            textPeriod.FitBoxToText = 'on';
            textPeriod.FontSize = size1;
            textPeriod.Color = color1;
            this.textPeriod = textPeriod;
            drawnow
            textTime = annotation('textbox','Interpreter','latex');
            textTime.Position = [textPeriod.Position(1)+1.1*textPeriod.Position(3) .98 0 0];
            textTime.String = ['Computation time = ',timestring(1)];
            textTime.EdgeColor = 'none';
            textTime.FitBoxToText = 'on';
            textTime.FontSize = size2;
            textTime.Color = color2;
            this.textTime = textTime;
            drawnow            

            textCX = annotation('textbox','Interpreter','latex');
            textCX.Position = [this.textTime.Position(1)+1.4*this.textTime.Position(3)  .98 0 0];
            textCX.String = ['$C_X$ = ',num2str(1,'%1.2e')];
            textCX.EdgeColor = 'none';
            textCX.FitBoxToText = 'on';
            textCX.FontSize = size2;
            textCX.Color = color2;
            this.textCX = textCX;
            drawnow
            textCV = annotation('textbox','Interpreter','latex');
            textCV.Position = [textCX.Position(1)+1.1*textCX.Position(3)  .98 0 0];
            textCV.String = ['$C_V$ = ',num2str(1,'%1.2e')];
            textCV.EdgeColor = 'none';
            textCV.FitBoxToText = 'on';
            textCV.FontSize = size2;
            textCV.Color = color2;
            this.textCV = textCV;
            drawnow
            textCL = annotation('textbox','Interpreter','latex');
            textCL.Position = [textCV.Position(1)+1.1*textCV.Position(3)  .98 0 0];
            textCL.String = ['$C_\textup{loss}$ = ',num2str(1,'%1.2e')];
            textCL.EdgeColor = 'none';
            textCL.FitBoxToText = 'on';
            textCL.FontSize = size2;
            textCL.Color = color2;
            this.textCL = textCL;
        end
        
        function UpdateTopRowExtend(this,numPeriods,compTime,costX,costV,costLoss)
%             this@PlotterTopRow(size1,size2,color1,color2);
            this.maxPeriods = max(numPeriods,this.maxPeriods);
            this.textPeriod.String = ['\# Periods = ',num2str(this.maxPeriods)];
            this.textTime.String = ['Computation time = ',timestring(compTime)];
            this.textCX.String = ['$C_X$ = ',num2str(costX,'%1.2e')];
            this.textCV.String = ['$C_V$ = ',num2str(costV,'%1.2e')];
            this.textCL.String = ['$C_\textup{loss}$ = ',num2str(costLoss,'%1.2e')];
        end
        
    end    
end

