classdef SimText < handle
    % Writes the Period and time information to top row of output window
    
    properties
        parent
        
        textSimTime
        textOpt
        textMacroparticles
        textMacroweight
        textLoss
    end
    
    methods
        
        function this = SimText(parent,size1,color1)
            this.parent = parent;
            textSimTime = annotation('textbox');
            textSimTime.Position = [parent.PhiHandle.SubPlot.Position(1)+.06, parent.PhiHandle.SubPlot.Position(2)+parent.PhiHandle.SubPlot.Position(4)-.005, .2,0];
            textSimTime.String = 't';
            textSimTime.EdgeColor = 'none';
            textSimTime.FitBoxToText = 'on';
            textSimTime.FontSize = size1;
            textSimTime.Color = color1;
            this.textSimTime = textSimTime;
            drawnow
            textLoss = annotation('textbox');
            textLoss.Position = [textSimTime.Position(1), textSimTime.Position(2)+.01, .2,0];
            textLoss.String = ['Macroparticles lost ='];
            textLoss.EdgeColor = 'none';
            textLoss.FitBoxToText = 'on';
            textLoss.FontSize = size1;
            textLoss.Color = color1;
            this.textLoss = textLoss;
            drawnow
            textMacroparticles = annotation('textbox');
            textMacroparticles.Position = [parent.RhoHandle.SubPlot.Position(1)+.06, parent.RhoHandle.SubPlot.Position(2)+parent.RhoHandle.SubPlot.Position(4)-.005, .2,0];
            textMacroparticles.String = ['M'];
            textMacroparticles.EdgeColor = 'none';
            textMacroparticles.FitBoxToText = 'on';
            textMacroparticles.FontSize = size1;
            textMacroparticles.Color = color1;
            this.textMacroparticles = textMacroparticles;
            drawnow
            textMacroweight = annotation('textbox');
            textMacroweight.Position = [textMacroparticles.Position(1), textMacroparticles.Position(2)+.01, .2,0];
            textMacroweight.String = ['macroweight = ',num2str(parent.macroWeight,'%1.1e')];
            textMacroweight.EdgeColor = 'none';
            textMacroweight.FitBoxToText = 'on';
            textMacroweight.FontSize = size1;
            textMacroweight.Color = color1;
            this.textMacroweight = textMacroweight;
            drawnow
            textOpt = annotation('textbox');
            textOpt.Interpreter = 'tex';
            textOpt.Position = [textMacroweight.Position(1), textMacroweight.Position(2)+.01, .2,0];
            textOpt.String = 'O';
            textOpt.EdgeColor = 'none';
            textOpt.FitBoxToText = 'on';
            textOpt.FontSize = size1;
            textOpt.Color = color1;
            this.textOpt = textOpt;
        end
        
        function UpdateSimText(this,time,npKill,OptimizerText,np,it)
            if isempty(OptimizerText)
                this.textOpt.String = '';
            else
                this.textOpt.String = [OptimizerText,', ',num2str(it)];
            end
            this.textSimTime.String = ['t [s] = ',num2str(time,'%1.2e')];
            this.textLoss.String = ['Macroparticles lost = ',num2str(npKill)];
            this.textMacroparticles.String =  ['# macroparticles =',num2str(np)];
        end
        
    end
    
end

