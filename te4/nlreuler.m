function [x]=nlreuler(plant,t,x0,u0,h,p1,p2,p3,p4,p5)

% function [x]=nlreuler(plant,t,x0,u0,h,p1,p2,p3,p4,p5)

% Euler integration of "plant" with fixed step size, h.
% Integration is from t(1) to t(2).
% p1 to p5 are optional parameters to be passed to "plant".
%
% x is FINAL state.


% ---------Extra parameters for S-function systems ------
outstr=[];
for i=1:nargin-5
	outstr=[outstr,',p',int2str(i)];
end

% Set up calls for state derivatives

call1=['dxdt=feval(plant,time,x,u0,1',outstr,');'];

% Figure out how many steps are needed.

dt=t(2)-t(1);
nstep=ceil(dt/h);
hused=dt/nstep;        % Actual step size
x=x0;                  % Initial condition
time=t(1);
for i=1:nstep
   eval(call1)
   x=x+dxdt*hused;
   time=time+hused;
end
