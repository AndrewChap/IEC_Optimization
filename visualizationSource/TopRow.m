classdef TopRow
    % Writes the Period and time information to top row of output window
    
    properties
        textPeriod
        textTime
    end
    
    methods
        
        function this = TopRow(size1,size2,color1,color2)
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
            
        end
        
        function UpdateTopRow(this,numPeriods,compTime)
            this.textPeriod.String = ['\# Periods = ',num2str(numPeriods)];
            this.textTime.String = ['Computation time = ',timestring(compTime)];
        end
    end
    
end

