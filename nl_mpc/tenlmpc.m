% CONTROL THE TEC PROCESS USING NL MPC.
% Kalman estimator 

% Initialize the plant model

idv=zeros(21,1);
p2=zeros(153,1);
p3=zeros(586,1);
p4=zeros(139,1);

% The following call to te_mex initializes p2, p3, and p4, and
% sets x0 to the default initial state vector.

[sys,x0]=te_mex_TC(-1,[],[],0,idv,p2,p3,p4); 

% Override default state with another...  Do this by changing
% the following definition of the Mode variable.  Mode 0 is
% the initial condition of Downs and Vogel.  Modes 1-6 are
% the optimal steady-state conditions determined by Ricker
% (Computers & Chemical Engineering, vol. 19, pp. 949-959(1995).)

Mode='0';
if exist('Mode') == 1 & Mode ~= '0'
    fprintf(['\n*** Loading Mode ',Mode,' ***\n'])
    eval(['load mode',Mode,'.dat -ascii;  x0=mode',Mode,';'])
end

% Set default inputs and initial outputs for plant:

u0=x0(39:50);
y0=TE_mex_TC(0,x0,u0,3,idv,p2,p3,p4);
u0(10:11)=y0([9,11]);
x0(51)=(x0(48)-41.10581288)/(-64);	% set integrator states in feedback
x0(52)=(x0(49)-18.11349055)/(-16);	% loops.
dxdt0=TE_mex_TC(0,x0,u0,1,idv,p2,p3,p4);
[max_dxdt0,i]=max(dxdt0)

xp=x0;
up=u0;
yp=y0;
nyp=length(yp);

% Initialize the estimator

Fconv=[11.2/.25052; 114.5/3664; 98./4509.3; 417.5/9.3477; ...
       15.1/.33712; 259.5/38.1; 211.3/46.53; 1; 1];
 
uest =[y0(1:4,1).*Fconv(1:4,1)
   y0(5)*44.79
   y0(10)*Fconv(5)
   u0(7:8,1).*Fconv(6:7,1)
   y0(9)
   y0(11)
   48.5
   .5
   100*ones(13,1)];

iymeas=[7,8,13,12,15,16,6,23:36,40,41];
yest=zeros(32,1);
yest(1:23)=y0(iymeas);

% calculate steady-state of estimator for specified inputs & outputs

[x,u1,y_ss]=TEest3_ss(yest,uest);

dt=0.01;
ef=[];  eg=[]; ep=[];
F=[];
ygas=-1;yprod=-1;
xest=x;
dest=u1(11:25);
yest=TEest3(0,xest,u1,3);
dxdt_est_max=max(abs(TEest3(0,xest,u1,1)))

% The following statements load the estimator gain matrix and
% pointers used by the estimator.  These are ASCII files.
% The load may not work properly on all systems.

load Kest1.dat
load idest.dat
load iyest.dat
load iymeas.dat

dt_est=0.1;

% Set up slave control loops (PI control)

icv=[1,2,3,4,10, 8,12,15];
imv=[3,1,2,4, 6,11, 7, 8];
setpts=yp(icv)';
yscPI=[1 .01 .01 1 1 1 1 1];

Kcs  =[  0.03, 3e-6, 2e-6,1.2e-3,0.015, 0.8, -0.2 -0.2];
tauis=[ 0.001,0.001,0.001, 0.001,0.001,  60,   60   60]/60;
nloops=length(Kcs);
errn1=zeros(1,nloops);

% Set up override loop(s)

Prod_adj=0;
ePa=0;

% Load saved run here!

%load lastrun

% Reset commons in case of a re-load

sys=TE_mex_TC(-1,[],[],0,idv,p2,p3,p4); 

% Set up NL MPC

% outputs:
%
% 1:  Production               [kmol/h]
% 2:  Pressure                 [kPa]
% 3:  A in feed                [mole %]
% 4:  E in feed                [mole %]
% 5:  B in purge               [mole %]
% 6:  G in product             [mole %]
% 7:  Liq in Separator         [%]
% 8:  Liq in Stripper          [%]


% Pointers and conversion factors

iy=[17,7,23,27,30,40,12,15];      % points to controlled variables
                               % NOTE:  17 is a dummy ... replaced by u(8)
iy_est=[24,1,8,12,15,22,4,5];    % points to corresponding estimator outputs
ny=length(iy);
iu=[1:4,6,9,7,8];                % points to manip. variables in estimator
nu=length(iu);
iv=[];
id=[];
yscMPC=[1 .01 1 1 1 1 1 1];

% MPC parameters

dt_MPC=0.1;
Nhor=10;
M=[2 3 5];
ywt=  [.6 .034 .2 .2 .2 2 .05 .05];
uwt=2*[0.1 0.1 0.05 0.1 0.15 0.2 .01 .01];
recalc_incr=5;
recalc_cntr=recalc_incr;

% Constraints

umin=[.1   1   1   1  .1  115 10 10];
umax=[45 181 181 681  32  129 90 90];  % Normal settings
dumx=[10  10  10  20  10    2 10 10];
ulim=[umin  umax dumx ];
ylim=[];
ylim_wt=[];

Pmax=2920;

% Set up to lag the outputs that otherwise have instantaneous
% responses to input variations.

af=-300*eye(2);
bf=zeros(2,24);  bf(1,1)=1;    bf(2,24)=1;
cf=zeros(24,2);  cf(1,1)=300;  cf(24,2)=300;
df=eye(24);      df(1,1)=0;    df(24,24)=0;

% Linearization setup

upert=1e-5+(1e-8*abs(u1));
upert(11:25)=zeros(15,1);

% Simulation time setup

t0=clock;
t=0;
tend=72;  % This defines desired elapsed time
format short

nstep=round(tend/dt)+1;        % Number of sampling periods in simulation.
inc_MPC=round(dt_MPC/dt);      % ratio of MPC to PI control interval
ii_MPC=[1:inc_MPC:nstep];
N_MPC=length(ii_MPC);          % Number of MPC moves in simulation
ii_next=1;
iiMPC=1;
isave=1;
save_inc=10;        
save_cnt=save_inc;
t_save=0;

% Setpoints

up(5)=0;      % maximize recycle rate
up(9)=0;      % turn off steam
up(12)=100;   % set agitator speed at max
setpts(6)=65; % set reactor liquid level

Nset=N_MPC+Nhor+1;
set_MPC=ones(Nset,1)*[339 2800 36.4  17  15.9  53.7 50 50];
delset= ones(Nset,1)*[  0    0   0    0    0      0  0  0];

% The following are example alternative setpoint definitions to
% handle other cases:
%delset= ones(Nset,1)*[    0    0   0      0     0     0   0  0];
%delset=ones(Nset,1)*[-32.3  145 -2   -14  10    36.4    0  0];
%NN=Nset-250;
%delset(251:Nset,:)= ones(NN,1)*[ -15.3  145  2.6    7    -2   -42  0  0];
%delset= ones(Nset,1)*[  150  145   4    -1   0     0   0  0];

% Calculate setpoint trajectory

phir=diag([.98 .95 .98 .98 .98 .98 .96 .96]); 
gamr=eye(ny)-phir; 
cr=eye(ny); dr=zeros(ny,ny);
set_MPC=set_MPC(:,1:ny)+dlsimm(phir,gamr,cr,dr,delset);
clear delset

% Constant Disturbances -- currently set for simultaneous
% disturbances 8 and 13.  Change as you wish.

%idv(1)=1;            % Step in A/C ratio
%idv(6)=1;            % Loss of A feed
idv(8)=1;            % random variation of A,B,C in stream 4
idv(13)=1;           % drift in kinetics
%idv(16:20)=ones(5,1);        % Mystery disturbances

idv_0=idv'

%++++++++++++++++  MAIN LOOP +++++++++++++++++++++

% Allocate storage

y=zeros(N_MPC,nyp);
u=zeros(N_MPC,12);
d_s=zeros(N_MPC,15);
yest=zeros(N_MPC,32); 
slave_set=zeros(N_MPC,length(setpts)); 
t_save=zeros(N_MPC,1); 
MPC_set=zeros(N_MPC,ny);
yP1=[];Prod_adj_s=zeros(N_MPC,1);
u1_s=zeros(nstep,10);

for ii=1:nstep  

%  Time-varying disturbances & setpoints -- can be turned on
%  by removing comment characters.

%   if t(ii) < 2
%      up(5)=15;      % recycle rate
%      up(9)=40;      % steam
%      up(12)=60;     % agitator speed
%      setpts(6)=70;  % reactor liquid level
%   elseif t(ii) < 4
%      up(5)=10;      % recycle rate
%      up(9)=30;      % steam
%      up(12)=70;     % agitator speed
%      setpts(6)=65;  % reactor liquid level
%   elseif t(ii) < 6
%      up(5)=5;      % recycle rate
%      up(9)=30;      % steam
%      up(12)=80;     % agitator speed
%      setpts(6)=60;  % reactor liquid level
%   elseif t(ii) < 8
%      up(5)=0;      % recycle rate
%      up(9)=20;      % steam
%      up(12)=90;     % agitator speed
%      setpts(6)=60;  % reactor liquid level
%   elseif t(ii) < 10
%      up(5)=0;      % recycle rate
%      up(9)=0;      % steam
%      up(12)=100;     % agitator speed
%      setpts(6)=60;  % reactor liquid level
%   else
%      up(5)=0;      % recycle rate
%      up(9)=0;      % steam
%      up(12)=100;     % agitator speed
%      setpts(6)=60;  % reactor liquid level
%   end
   
% following is for maximum production case.  Reduces offset in other
% controlled variables when production setpoint is infeasible.

   if max(up([1:2,4,6])) > 98
      ywt(1)=0.06;	% reduce weight on production setpoint
   end

%  Get current plant outputs

   yp=TE_mex_TC(t(ii),xp,up,3,idv,p2,p3,p4);
   ymeas=yp(iy);
   y(isave,:)=yp';
   MPC_set(isave,:)=set_MPC(iiMPC,:);

%  Run the estimator(s) to get the "measurement update"

   y1=TEest3(t(ii),xest,u1,3);  % current estimator outputs (before
                                % state correction).
   if ygas ~= yp(23)
      xkf=[xest;u1(idest)];
      xkf=xkf+Kest1*(yp(iymeas)-y1(iyest));
      xest=xkf(1:26);
      dest(idest-10)=xkf(27:length(xkf));
      u1(11:25)=dest;
      y1=TEest3(t(ii),xest,u1,3);  % current estimator outputs (after
                                   % state correction).
      ygas=yp(23);
   end

   yest(isave,:)=y1';
   d_s(isave,:)=dest';

%  If this is time for an MPC move, do it!

   if ii == ii_next

%  Linearize NL model.  Account for inventory loops and non-steady
%  initial condition.

      [a1,b1,c1,d1]=linmod_nlr('TEest3',xest,u1,[],[],upert);
      c1=c1(1:24,:);
      d1=d1(1:24,:);
      icv1=[2,4,5];
      imv1=[10,7,8];
      K1=[ 0.8 -1.36 -0.908];
      Tau1=[60, 60    60 ]/60;
      [a1,b1,c1,d1]=PI_feedback(a1,b1,c1,d1,K1,Tau1,icv1,imv1);
      [a1,b1,c1,d1]=series(a1,b1,c1,d1,af,bf,cf,df);
      c=c1(iy_est,:);
      d=d1(iy_est,iu);
      if max(max(abs(d))) < 1e-12
         d=zeros(ny,nu);
      else
         error('Bad D in linearization')
      end
      f0=TEest3(t(ii),xest,u1,1);  % current derivatives of states
      e0=setpts(1,6:8)'-y1(icv1,1);
      [phi,gam]=c2d(a1,[b1 [f0;e0;0;0]],dt_MPC);
      [rows,cols]=size(gam);
      del0=gam(:,cols);  % save effect of non-steady initial condition
      gam=gam(:,iu);     % pick out other desired columns

% Linear projection of plant outputs at current state and inputs.
      
      yP1=dstep(phi,[del0],c,zeros(ny,1),1,Nhor+1)+...
                       ones(Nhor+1,1)*y1(iy_est,1)';
      yP=reshape(yP1(2:Nhor+1,:)',ny*Nhor,1);

% Update uncertain flow limits

      u_now=[setpts(1,1:5)'.*Fconv(1:5,1); up(10,1);...
                     setpts(1,7:8)'];
      if up(3) > 50 & yp(1) < 0.2 
         umax(1)=0.2;
         u_now(1)=umax(1);
         setpts(1,1)=u_now(1)/Fconv(1);
      else
         umax(1)=45;
      end
      if up(6) > 99
         umax(5)=u_now(5);
         setpts(1,5)=u_now(5)/Fconv(5);
      else
         umax(5)=32;
      end

      ulim=[umin  umax  dumx ];

% Control action
      
      [nn,nn]=size(phi);
      imod=ss2mod(phi,gam,c,d,[dt_MPC,nn,nu,0,0,ny,0]);
      [nc,TAB,a,rhscon,iumin,iumax,iymin,iymax,ymin_wt,ymax_wt,...
      dumax,SuTQ]=smpcqp0s(imod,ywt,uwt,M,Nhor,ulim,ylim,ylim_wt);
      r=reshape(set_MPC(iiMPC:iiMPC+Nhor-1,:)',ny*Nhor,1);

%         Override to decrease production setpoint when P is
%         getting too high

      ePan=Pmax-yp(7);
      Prod_adj=Prod_adj-0.3*(ePan-ePa+(dt_MPC*(ePan+ePa))/(20/60));
      ePa=ePan;
      Prod_adj=max([0 Prod_adj]);
      Prod_adj_s(iiMPC,1)=Prod_adj;
      r([1:ny:ny*Nhor])=r([1:ny:ny*Nhor])-Prod_adj;
      deltau=mpcqps(nc,TAB,a,rhscon,iumin,iumax,iymin,...
                       iymax,ymin_wt,ymax_wt,dumax,SuTQ,r,yP,u_now);
      setpts(1,1:5)=setpts(1,1:5)+(deltau(1:5,1)./Fconv(1:5,1))';  
      up(10,1)=up(10,1)+deltau(6,1);
      setpts(1,7:8)=setpts(1,7:8)+deltau(7:8,1)';

%  Display 

fprintf('\n*** TIME = %4.1f/%4.1f hours ***\n\n',t(ii),tend)
fprintf(['        kmol/h    Pr      pctA    pctE',...
         '    pctB   pctG    SepLev  PrdLev \n'])
fprintf('setpt:  %6.2f  %6.1f  %6.2f',set_MPC(iiMPC,1)-Prod_adj,...
         set_MPC(iiMPC,2),set_MPC(iiMPC,3))
fprintf(      '  %6.2f  %6.2f  %6.2f',set_MPC(iiMPC,4),...
         set_MPC(iiMPC,5),set_MPC(iiMPC,6))
fprintf(      '  %6.2f  %6.2f\n',set_MPC(iiMPC,7),set_MPC(iiMPC,8))
fprintf('  PVs:  %6.2f  %6.1f  %6.2f',xp(46,1)*4.541,y(isave,7),y(isave,23))
fprintf(      '  %6.2f  %6.2f  %6.2f',y(isave,27),y(isave,30),y(isave,50))
fprintf(      '  %6.2f  %6.2f\n',y(isave,12),y(isave,15))
fprintf('Trend:  %6.2f  %6.1f  %6.2f',yP1(Nhor,1)-yP1(1,1),...
           yP1(Nhor,2)-yP1(1,2),yP1(Nhor,3)-yP1(1,3))
fprintf(      '  %6.2f  %6.2f  %6.2f',yP1(Nhor,4)-yP1(1,4),...
           yP1(Nhor,5)-yP1(1,5),yP1(Nhor,6)-yP1(1,6))
fprintf(      '  %6.2f  %6.2f\n',yP1(Nhor,7)-yP1(1,7),...
           yP1(Nhor,8)-yP1(1,8))
fprintf(' Bias:  %6.2f  %6.1f  %6.2f',xp(46,1)*4.541-y1(24),...
         y(isave,7)-y1(1),y(isave,23)-y1(8))
fprintf(      '  %6.2f  %6.2f  %6.2f',y(isave,27)-y1(12),...
         y(isave,30)-y1(15),y(isave,40)-y1(22))
fprintf('  %6.2f  %6.2f\n\n',y(isave,12)-y1(4),y(isave,15)-y1(5))

fprintf('           A       D       E      AC     Purge   ReacT  SepLev   PrdLev\n') 
fprintf('    U:  %6.2f  %6.2f  %6.2f',u_now(1),u_now(2),u_now(3))
fprintf(      '  %6.2f  %6.2f  %6.2f',u_now(4),u_now(5),u_now(6))
fprintf(      '  %6.2f  %6.2f\n',       u_now(7),u_now(8))
fprintf(' delU:  %6.2f  %6.2f  %6.2f',deltau(1),deltau(2),deltau(3))
fprintf(      '  %6.2f  %6.2f  %6.2f',deltau(4),deltau(5),deltau(6))
fprintf(      '  %6.2f  %6.2f\n\n',       deltau(7),deltau(8))

fprintf(' Loop:  React  React  Feed   Feed   Feed   Feed   Purge  Separ  Strip\n')
fprintf('        Level  Temp.    1      2      3      4           Level  Level\n')
fprintf('setpt:  %5.1f  %5.1f  %5.3f',setpts(6),up(10),setpts(1))
fprintf('  %5.0f  %5.0f  %5.2f',setpts(2),setpts(3),setpts(4))
fprintf('  %5.3f  %5.1f  %5.1f\n',setpts(5),setpts(7),setpts(8))
fprintf('   PV:  %5.1f  %5.1f  %5.3f',yp(8),yp(9),yp(1))
fprintf('  %5.0f  %5.0f  %5.2f',yp(2),yp(3),yp(4))
fprintf('  %5.3f  %5.1f  %5.1f\n',yp(10),yp(12),yp(15))
fprintf('pctMV:  %5.1f  %5.1f  %5.1f',xp(49),xp(48),xp(41))
fprintf('  %5.1f  %5.1f  %5.1f',xp(39),xp(40),xp(42))
fprintf('  %5.1f  %5.1f  %5.1f\n\n',xp(44),xp(45),xp(46))
fprintf(' Feed:    A      B      C      D      E      F   PROD:    E      F\n')
fprintf(' Meas:  %5.2f  %5.2f  %5.2f',yp(23),yp(24),yp(25))
fprintf('  %5.2f  %5.2f  %5.2f',yp(26),yp(27),yp(28))
fprintf('         %5.2f  %5.2f\n',yp(38),yp(39))
fprintf(' Bias:  %5.2f  %5.2f  %5.2f',yp(23)-y1(8),yp(24)-y1(9),yp(25)-y1(10))
fprintf('  %5.2f  %5.2f  %5.2f\n\n',yp(26)-y1(11),yp(27)-y1(12),yp(28)-y1(13))
fprintf('Purge:    A      B      C      D      E      F      G      H\n')
fprintf(' Meas:  %5.2f  %5.2f  %5.2f',yp(29),yp(30),yp(31))
fprintf('  %5.2f  %5.2f  %5.2f',yp(32),yp(33),yp(34))
fprintf('  %5.2f  %5.2f\n',yp(35),yp(36))
fprintf(' Bias:  %5.2f  %5.2f  %5.2f',yp(29)-y1(14),yp(30)-y1(15),yp(31)-y1(16))
fprintf('  %5.2f  %5.2f  %5.2f',yp(32)-y1(17),yp(33)-y1(18),yp(34)-y1(19))
fprintf('  %5.2f  %5.2f\n',yp(35)-y1(20),yp(36)-y1(21))

      ii_next=ii_next+inc_MPC;        % point to next time
      iiMPC=iiMPC+1;
   end


%  Run PI loops

   errn=setpts-[yp(icv(1:5));y1(icv1)]';
   delu=Kcs.*(errn-errn1+(dt*(errn+errn1))./tauis);
   up(imv,1)=up(imv,1)+delu';
   errn1=errn;
   up(imv,1)=min([max([up(imv,1)';zeros(1,nloops)]);100*ones(1,nloops)])';
   slave_set(isave,:)=setpts;
   u(isave,:)=[up(1:9,1)' xp(48:50,1)'];
   t_save(isave,1)=t(ii);

% Break out of loop here if ii == nstep.  Prevents extra integration of
% plant equations.

   if ii == nstep
      break
   end

% Filter manipulated flowrates to remove noise.  
% Then integration to update estimator states.

   [F,u1]=flowfilt(up,yp,F);
   u1=[u1;dest];
   tvec=[0;dt];
   u1_s(ii,:)=u1(1:10)';
   xest=nlr_euler('TEest3',tvec,xest,u1,1/3600);

% Simulate plant to next sampling time.

   t(ii+1,1)=t(ii,1)+dt;
   tvec=t(ii:ii+1,1);          % Starting and ending time.
   xp=nlr_euler('TE_mex_TC',tvec,xp,up,1/3600,idv,p2,p3,p4);

% Update save counter

   if save_cnt == save_inc
      isave=isave+1;
      save_cnt=1;
   else
      save_cnt=save_cnt+1;
   end

% Loop back for next step
end

% Save the data to disk.  This allows the run to be restored for a
% re-start at the last time point.  To do this, remove the comment
% character from the "load lastrun" command near the start of the
% script.

save lastrun
e_time=etime(clock,t0)/60

% The following commands produce a series of "custom" plots.
% PLOTPAIR plots the MPC setpoints (as defined by the setpoint
% trajectory) vs the measured outputs.  The PLOTEACH command
% plots the manipulated variables.  The ESTPLT1 function produces
% plots of other key variables.

[ii,dum]=size(u);
set_MPC(1:isave,1)=set_MPC(1:isave,1)-Prod_adj_s(1:isave,1);
plotpair(t_save(1:isave),MPC_set(1:isave,:),[u(1:isave,8)*4.541 y(1:isave,iy(2:ny))])
pause
ploteach([],u(1:isave,:),t_save(1:isave,:))
est_plt1(0.1,y(1:isave,:),yest(1:isave,:),d_s(1:isave,:))

% Optional saving of results in binary format.

%save IDV6_MPC.mat y u d_s yest slave_set t_save MPC_set ny iy isave u1_s
