% Matlab script to simulate multiloop PI control of 
% 2-phase reactor process.

% See: 'MPC of a continuous, nonlinear, two-phase reactor',
% N. L. Ricker, J. Process Control, vol. 3, 109-123(1993).

% NOTE:  set up to integrate ODEs using a crude Euler method
%        with fixed step size.  See NLREULER function for
%        more details.  Use better integration routine if
%        you have one.

% By default, runs Scenario I (disturbance in yA1 feed composition).
% Other scenarios are easy to do, and some can be switched on by
% changing comment lines in the code, below.

% A 30-hour simulation takes about 3 minutes to run on a Mac II-fx.
% You can reduce this by using a better integration routine.

% The results plotted at the end of the simulation should reproduce
% those shown in Fig. 6 of the paper cited above.


% copyright N. L. Ricker
% University of Washington
% Chemical Engineering
% Box 351750
% Seattle, WA 98195-1750
% ricker@cheme.washington.edu


% Initialize model

Tdelay=0.1;               % Sampling delay for gas composition
ya1=0.485;                % Feed 1 mole fraction of A
yb1=0.005;                % Feed 1 mole fraction of B
u2max=100;                % Maximium possible valve position, Feed 2
u4bar=47.03024823457651;  % Nominal steady-state for product valve
VLsp=44.17670682730923;   % Setpoint for liquid level (%)
KcVL =-1.4;               % Level controller gain (%/%)
kpar=0.00117;             % Pre-exponential
nCpar=0.4;                % Exponent on Pc

% Get steady-state condition

ya3=0.47;
F4=100;
P=2700;
p1=[Tdelay ya1 yb1 u2max KcVL u4bar 0 kpar nCpar zeros(1,41)];
[x0,VL]=te_4_ss(ya3,F4,P,p1);

VLpct=VL*100/30;
u0=[];
u0(1:3,1)=x0(5:7,:);
u0(4)=VLpct;
u4bar=x0(8);

% Define possible disturbances and drift...

del_kpar=0; tdrift=48;
del_nCpar=0;
ya1=0.45       % Causes a disturbance in feed composition
%u2max=0       % Causes loss of F2
%kpar=0.00117,  del_kpar=-0.0002, tdrift=48,  % Causes drift
%nCpar=0.4;     del_nCpar=-0.1,               % in kinetics

% Get outputs at the steady-state initial condition.

p1=[Tdelay ya1 yb1 u2max KcVL u4bar 0 kpar nCpar zeros(1,41)];
te_4(-1,[],[],0,p1);
y0=te_4(0,x0,u0,3,p1);

% Initialize plant states and inputs for simulation.

up=u0;
xp=x0;

% Other parameters

dt=0.1;                  % controller sampling period (hours)
stepsz=0.002;            % step size for nlreuler (hours)

t=0;                     % initial time
tend=30;                 % final time, hours
thour=0;                 % used to decide whether to display to screen

% Define parameters for the basic control loops:

iydisp=[4:9];            % indeces of displayed outputs
icv=[4,5,7];             % indeces of controlled outputs
imv=[1,3,2];             % indeces of manipulated inputs
setpts=y0(icv)';         % setpoint are initialized to measurements

tauis=[1.0, 1.5, 3.0];       % integral times
Kcs=[.1, -.25, 2.0];         % gains
errn1=zeros(1,length(Kcs));  % initialize integral errors to zero

% similarly for pressure override loop

errn1PC=0;
KcPC=0.7;
TauiPC=3;
F4sp_adj=0;
F4sp=[];

% Set up for simulation

nstep=round(tend/dt)+1;    % Number of sampling periods in simulation.
y=[]; u=[];  yB3sp=[];
dxdt0=te_4(0,x0,u0,1,p1)'  % displays initial state derivatives -- should be 
                           % near zero, e.g., of order e-13.

for ii=1:nstep          % MAIN LOOP

% The following simulates drift in the kinetic parameters:

   if t(ii) < tdrift
      p1(8)=kpar+del_kpar*(t(ii)/tdrift);
      p1(9)=nCpar+del_nCpar*(t(ii)/tdrift);
   else
      p1(8)=kpar+del_kpar;
      p1(9)=nCpar+del_nCpar;
   end

% Get current plant outputs

   yp=te_4(t(ii),xp,up,3,p1)';
   y(ii,:)=yp;

% Pressure control override loop

   errnPC=2900-yp(5);
   F4sp_adj=min([F4sp_adj+KcPC*(errnPC-errn1PC+(dt*errnPC)/TauiPC),0]);
   errn1PC=errnPC;
   setpts_save=setpts;
   setpts(1)=setpts(1)+F4sp_adj;
   F4sp(ii,1)=setpts(1);

% PI control

   errn=setpts-yp(icv);
   delu=Kcs.*(errn-errn1+(dt*errn)./tauis);
   up(imv,1)=up(imv,1)+delu(1:3)';
   errn1=errn;
   up(imv,1)=min([max([up(imv,1)';zeros(1,3)]);100*ones(1,3)])';
   u(ii,:)=up';

% Display to screen during simulation

   if t(ii) >= thour
      disp(['Time = ',num2str(t(ii))])
      disp('Pressure & setpoint:'),disp([y(ii,5), setpts(2)]/100)
      disp('Product rate & setpoint:'),disp([y(ii,4) setpts(1)])
      disp('YA3 & setpoint:'),disp([y(ii,7) setpts(3)])
      disp('Compositions:'),disp(y(ii,7:9))
      disp('u:'),disp([up'])
      thour=thour+1;
   end

% Restore setpoints to original values

   setpts=setpts_save;

% Break out of loop here if ii == nstep.  Prevents extra integration of
% plant equations.

   if ii == nstep
      break
   end

% Simulate plant to next sampling time.

   t(ii+1,1)=t(ii,1)+dt;
   tvec=t(ii:ii+1,1);          % Starting and ending time.
   xp=nlreuler('te_4',tvec,xp,up,stepsz,p1);

% Loop back for next step

end

% Plot results at end of simulation

clg
subplot(221)
plot(t,[y(:,4) F4sp]),title('Product Rate and Setpoint [kmol/h]')
xlabel('Time [h]')
subplot(222)
plot(t,y(:,5)),title('Pressure [kPa]')
xlabel('Time [h]')
subplot(223)
plot(t,y(:,7)),title('A in Purge [mole %]')
xlabel('Time [h]')
subplot(224)
plot(t,u(:,1:3)),title('Manipulated Variables [%]')
xlabel('Time [h]')
