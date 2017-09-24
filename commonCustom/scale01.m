function [ xScaled, minX, maxX ] = scale01( x )

    % Drew Chap 5/26/17
    % scales an array so it ranges from zero to one.  outputs the min and
    % max from before this scaling was done
    minX =  min(x(:));
    maxX = max(x(:));
    xScaled = (x - minX)/(maxX-minX);

end