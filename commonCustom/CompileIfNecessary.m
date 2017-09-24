function mexFunction = CompileIfNecessary(cudaname)
    % Drew Chap 5/4/2017
    % compiles a mexFile (specified by input cudaname) if it needs to be compiled and 
    % returns a handle to the function
    
    % get information about our source and exe files
    fileSource = dir([cudaname,'.cu']);
    fileExe = dir([cudaname,'.mexw64']);
    
    %initialize booleans
    compileit = false;
    compileUnsuccessful = true;
    compileTries = 0;
    if isempty(fileExe)  % if executable doesn't exist, we must compile the source
        compileit = true;
    else
        sourceTime = GetFileTime(fileSource.name);
        exeTime = GetFileTime(fileExe.name);
        % check if the file has been compiled since we last modified it
        gt = sourceTime.Write > exeTime.Write;
        eq = sourceTime.Write == exeTime.Write;

        if gt(find(~eq,1,'first'))
            compileit = true;
        end
    end

    if compileit
        disp([cudaname, ' was created or modified since last compilation, compiling mex file...'])
%         eval(['mexcuda ',cudaname,'.cu -compatibleArrayDims -largeArrayDims -lcurand'])   % include '-lcurand' because of http://wlandau.github.io/gpu/lectures/curand/curand.pdf
%         eval(['mexcuda ',cudaname])
        outputName = cudaname;
        while compileUnsuccessful && compileTries < 10
            try
                mexcuda('-output',outputName,[cudaname,'.cu']); % compile and specify output name
                compileUnsuccessful = false;
            catch ME
                
                if strcmp(ME.message(1:4), 'LINK')
                    compileTries = compileTries + 1;
                    message1 = ['Cannot compile ',outputName,'. '];
                    outputName = [cudaname,num2str(compileTries,'_%02i')];
                    message2 = ['Going to try to compile as ',outputName];
                    warning([message1,message2]);
                else
                    error(ME.message);
                end
            end
        end
        if compileUnsuccessful
            error('Cannot compile after multiple tries')
        else
            disp('Compiling complete')
        end
        mexFunction = outputName;
    else
        disp([cudaname, ' has not been modified since last compilation'])
        mexFunction = cudaname;
    end
    
end