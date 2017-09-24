function newx = annealingFunction(optimValues,problem)
%ANNEALINGFAST Generates a point using Student's t distribution
%   NEWX = ANNEALINGFAST(optimValues,problem) generates a point based on
%   the current point and the current temperature using Student's t
%   distribution. 
%
%   OPTIMVALUES is a structure containing the following information:
%              x: current point 
%           fval: function value at x
%          bestx: best point found so far
%       bestfval: function value at bestx
%    temperature: current temperature
%      iteration: current iteration 
%             t0: start time
%              k: annealing parameter
%
%   PROBLEM is a structure containing the following information:
%      objective: function handle to the objective function
%             x0: the start point
%           nvar: number of decision variables
%             lb: lower bound on decision variables
%             ub: upper bound on decision variables
%
%   Example:
%    Create an options structure using ANNEALINGFAST as the annealing
%    function
%    options = optimoptions('simulannealbnd','AnnealingFcn',@annealingfast);


%   Copyright 2006-2015 The MathWorks, Inc.

currentx = optimValues.x;
nvar = numel(currentx);
newx = currentx;
y = randn(nvar,1);
y = y./norm(y);
newx(:) = currentx(:) + optimValues.temperature.*y;


newx = sahonorbounds(newx,optimValues,problem);
deltax = newx-currentx;
disp(['Step size = ',num2str(norm(deltax)),', all = ',num2str(deltax),', prev = ',num2str(currentx)]);
% if deltax < 1e-4
%     keyboard
% end

function newx = sahonorbounds(newx,optimValues,problem)
%SAHONORBOUNDS ensures that the points that SIMULANNEAL  move forward with
%   are always feasible.  It does so by checking to see if the given point
%   is outside of the bounds, and then if it is, creating a point called
%   which is on the bound that was being violated and then generating a new
%   point on the line between the previous point and the projnewx. It is
%   assumed that optimValues.x is within bounds.

%   Copyright 2006-2010 The MathWorks, Inc.

% Return if the problem is unbounded
if ~problem.bounded
    return
end

xin = newx; % Get the shape of input
newx = newx(:); % make a column vector
lb = problem.lb;
ub = problem.ub;
lbound = newx < lb;
ubound = newx > ub;
alpha = rand;
% Project newx to the feasible region; get a random point as a convex
% combination of proj(newx) and currentx (already feasible)
if any(lbound) || any(ubound)
    projnewx = newx;
    projnewx(lbound) = lb(lbound);
    projnewx(ubound) = ub(ubound);
    newx = alpha*projnewx + (1-alpha)*optimValues.x(:);
    % Reshape back to xin
    newx = reshapeinput(xin,newx);
else
    newx = xin;
end

function X = reshapeinput(Xin, X)
% RESHAPEINPUT reshape X to match the shape of Xin 

%   Copyright 2003-2012 The MathWorks, Inc.


[m,n] = size(Xin);
% Scalar input is row major
if m == n && m == 1 
    X = X(:); 
    return; 
end

% Single point evaluation
if isvector(X) % X is a vector so use shape information
    Xin(:) = X;
    X = Xin;
    return;
elseif isvector(Xin)
    if  m > n && n == 1  % col major input
        return;
    elseif n > m && m == 1 % row major input
        X = X';
    end
else % Xin is a matrix; not a documented feature
    p = size(Xin,2);
    if p > 1          % Matrix Xin with vectorized 'on' 
        X = reshape(X,m,n,p);
    else
        X = reshape(X,m,n);
    end
end
