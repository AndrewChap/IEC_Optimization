function [ output_string ] = timestring(varargin)
% takes the time and displays it as a string like a normal person would
% Drew Chap 9-16-14 (last edited 3/19/17
% 
% 1st input: time in seconds
% 2nd input (optional): number of decimals to include (defaults to 2 if not
% set)
% does not handle negative times

    time = varargin{1};
    
    latex = false;

    numDecimals = 2;
        
    if nargin > 1
        if isnumeric(varargin{2})
            numDecimals = varargin{2};
            if nargin == 3
                if strcmp(varargin{3},'latex')
                    latex = true;
                end
            end
        else
            if strcmp(varargin{2},'latex')
                latex = true;
            end
            if nargin == 3
                if isnumeric(varargin{3})
                    numDecimals = varargin{3};
                end
            end
        end
    end

    stringSpecifier = ['%1.',num2str(numDecimals,'%i'),'f'];

    if latex
        lB = '\textup{';
        lE = '}';
    else
        lB = '';
        lE = '';
    end
    
    if time == 0
        output_string = [num2str(0,stringSpecifier),lB,' s',lE];
    elseif time < 1e-21
        output_string = [num2str(time*1e24,stringSpecifier),lB,' ys',lE];
    elseif time >= 1e-21 && time < 1e-18
        output_string = [num2str(time*1e21,stringSpecifier),lB,' zs',lE];
    elseif time >= 1e-18 && time < 1e-15
        output_string = [num2str(time*1e18,stringSpecifier),lB,' as',lE];
    elseif time >= 1e-15 && time < 1e-12
        output_string = [num2str(time*1e15,stringSpecifier),lB,' fs',lE];
    elseif time >= 1e-12 && time < 1e-9
        output_string = [num2str(time*1e12,stringSpecifier),lB,' ps',lE];
    elseif time >= 1e-9 && time < 1e-6
        output_string = [num2str(time*1e9,stringSpecifier),lB,' ns',lE];
    elseif time >= 1e-6 && time < 1e-3
        output_string = [num2str(time*1e6,stringSpecifier),' ',char(181),lB,'s',lE];
    elseif time >= 1e-3 && time < 1
        output_string = [num2str(time*1e3,stringSpecifier),lB,' ms',lE];
    elseif time >= 1 && time < 60
        output_string = [num2str(time,stringSpecifier),lB,' s',lE];
    elseif time >= 60 && time < 3600
        output_string = [num2str(time/60,stringSpecifier),lB,' minutes',lE];
    elseif time >=3600 && time < 86400
        output_string = [num2str(time/3600,stringSpecifier),lB,' hours',lE];
    elseif time >= 86400 && time < 604800
        output_string = [num2str(time/86400,stringSpecifier),lB,' days',lE];
    elseif time >= 604800 && time < 31536000
        output_string = [num2str(time/604800,stringSpecifier),lB,' weeks',lE];
    else
        output_string = [num2str(time/31536000,stringSpecifier),lB,' years',lE];
    end


end

