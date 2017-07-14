% Matlab script to test the TE_4_SS function.

% See comments in that function for more details on usage.

% copyright N. L. Ricker
% University of Washington
% Chemical Engineering
% Box 351750
% Seattle, WA 98195-1750
% ricker@cheme.washington.edu

% Here are the specified conditions for the steady-state
% calculation:

ya3=.47;          % mol fraction of A in the purge
F4=100;           % product rate (kmol/h)
P=2700;           % pressure (kPa)

Tdelay=0.1;       % Sampling delay for gas composition
ya1=0.485;        % Feed 1 mole fraction of A
yb1=0.005;        % Feed 1 mole fraction of B
u2max=100;        % Maximium possible valve position, Feed 2
u4bar=47.0813;    % Nominal steady-state for product valve
VLsp=43.4225;     % Setpoint for liquid level (%)
KcVL =-1.4;       % Level controller gain (%/%)
kpar=0.00117;     % Pre-exponential
nCpar=0.4;        % Exponent on Pc

p1=[Tdelay ya1 yb1 u2max KcVL u4bar 0 kpar nCpar zeros(1,41)];

% Now do the calculation:

[xss,VL]=te_4_ss(ya3,F4,P,p1);

% Use the TE_4 function to verify that it's really at steady-state.

ss_states=xss'       % this just displays the states on the screen
VLpct=VL*100/30;
u0=[];x0=xss;
u0(1:3,1)=x0(5:7,:);
u0(4)=VLpct;
u4bar=x0(8);
p1_new=[Tdelay ya1 yb1 u2max KcVL u4bar 0 kpar nCpar zeros(1,41)];
te_4(-1,[],[],0,p1_new);            % initialize model
uss=u0;
y0=te_4(0,xss,uss,3,p1_new);        % get outputs at this condition
dxdt0=te_4(0,xss,uss,1,p1_new)'     % get state derivatives at this condition

% NOTE:  the dxdt0 values should be close to zero, e.g. of order e-13.
