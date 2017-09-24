function [period, extent, accel] = PeriodEstimator(x,vInitial,a,dt)
    
% Drew Chap 1/31/2017
% x - position points
% xInitial=0 - assumed particle starts out at x = 0
% vInitial - initial velocity of particle
% a - acceleration of particle as a function of x
% dt - time-step
%
%    - Estimates the bounce-time in an IEC given the electric field of the
% symmetric potential well and the fusion velocity
%    - Also returns "extent" the x-location where the particle turns around
%    - Also returns "accel" the acceleration of the particle when it is
% turning around

    DT = dt/5;
    
    px = 0;
    pv = vInitial;
    
    t = 0;
    while pv > 0
        px = px + 0.5*pv*DT;
        eHere = interp1(x,a,px);
        pv = pv + eHere*DT;
        px = px + 0.5*pv*DT;
        t = t + DT;
    end
    
    period = 2*t;
    extent = px;
    accel = eHere;
    
end