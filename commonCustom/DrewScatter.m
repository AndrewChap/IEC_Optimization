function handle = DrewScatter( varargin )

    % Drew Chap 12/21/16
    % Scatter points easily without having to break it up into x,y vectors
    % or x,y,z vectors.  Just DrewScatter(points) where points is size NxD
    % where N is your number of points and D is either 2 or 3 (2D or 3D)
    % 
    % Will use "line" if a single color is specified, will use "scatter" if
    % an array of colors is specified
    
    pts = varargin{1};
    marker = [];
    color = [];
    pSize = [];
    
    for n = 2:nargin
        if ischar(varargin{n})
            marker = varargin{n};
        elseif numel(varargin{n})>1
            color = varargin{n};
        else
            pSize = varargin{n};
        end
    end
    
    if isempty(marker)
        marker = '.';
    end
    if isempty(color)
        color = [0 0 0];
    end
    if isempty(pSize)
        pSize=1;
    end
    
    if numel(color)==3
        if size(pts,2) == 2
            handle = line('XData',pts(:,1),'YData',pts(:,2),'linestyle','none','marker',marker,'markersize',pSize,'markeredgecolor',color);
        elseif size(pts,2) == 3
            handle = line('XData',pts(:,1),'YData',pts(:,2),'ZData',pts(:,3),'linestyle','none','marker',marker,'markersize',pSize,'markeredgecolor',color);
        end
    else
        if size(pts,2) == 2
            handle = scatter(pts(:,1),pts(:,2),pSize,color,'filled');
        elseif size(pts,2) == 3
            handle = scatter3(pts(:,1),pts(:,2),pts(:,3),pSize,color,'filled');
        end
    end

end