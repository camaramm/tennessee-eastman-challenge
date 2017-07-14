function [x,u,y,x7]=TEest5_ss(y0,u0)

% Calculates steady state of estimator to match certain measured outputs and 
% specified inputs.  Outputs, states, and inputs are arranged as in TEest5.

%	STATE VARIABLES are all molar holdups [kmol] in a certain
%	location.

%	 1	A in reactor
%	 2	B in reactor
%	 3	C in reactor
%	 4	D in reactor
%	 5	E in reactor
%	 6	F in reactor
%	 7	G in reactor
%	 8	H in reactor
%	 9	A in separator
%	10	B in separator
%	11	C in separator
%	12	D in separator
%	13	E in separator
%	14	F in separator
%	15	G in separator
%	16	H in separator
%	17	A in feed zone
%	18	B in feed zone
%	19	C in feed zone
%	20	D in feed zone
%	21	E in feed zone
%	22	F in feed zone
%	23	G in feed zone
%	24	H in feed zone
%	25  D in product
%	26  E in product
%	27  F in product
%	28	G in product
%	29	H in product
%	30  internal energy in separator [MJ/kmol]
%	31  internal energy in stripper [MJ/kmol]

%	INPUTS are flows [kmol/h] unless noted otherwise:

%	 1	Feed 1 (pure A)
%	 2	Feed 2 (pure D)
%	 3	Feed 3 (pure E)
%	 4	Feed 4 (A & C)
%	 5	Compressor valve [Percent]
%	 6	Purge (stream 9)
%	 7	Separator underflow (stream 10)
%	 8	Product rate (stream 11)
%	 9	Steam Flow [kg/h]
%	10	Reactor temperature [deg C]
%	11	Condenser Valve  [%]
%	12	A in stream 4 [Mole Percent]
%	13	B in stream 4 [Mole Percent]
%	14	Reaction 1 activity factor [Percent].
%	15	Reaction 2 activity factor [Percent].
%	16	Stream 10 bias flow [kmol/h].
%	17	Reactor/Sep flow parameter [Percent].
%	18	Feed/Reactor flow parameter [Percent].
%	19	Not used in this version.
%	20	Adjustment to VLE of D to G in separator [Percent]
%	21	Adjustment to H VLE in separator [Percent]
%	22	C bias flow to feed zone [kmol/h]
%	23	D bias flow to feed zone [kmol/h]
%	24	E bias flow to feed zone [kmol/h]
%	25	Reaction 3 activity factor [Percent]
%	26	Adjustment to VLE in reactor [Percent]
%	27	Bias correction to compressor performance curve
%	28	Bias correction to compressor power prediction
%	29	Product density correction [Percent]
%	30	# equil. stages for D in stripper
%	31	# equil. stages for E in stripper
%	32	# equil. stages for F in stripper
%	33	Condenser coolant inlet temperature [C]
%	34	UA in condenser [kW/K]
%	35	Heat loss in stripper [kW]

%	OUTPUTS are mole % in stream unless noted otherwise:

%	 1	Reactor pressure [kPa]
%	 2	Reactor liq. holdup [Percent]
%	 3	Separator pressure [kPa]
%	 4	Separator liq. holdup [Percent]
%	 5	Product liq. holdup [Percent]
%	 6	Feed zone pressure [kPa]
%	 7	Total feed entering reactor (stream 6) [kscmh]
%	 8	A in reactor feed (stream 6)
%	 9	B in reactor feed (stream 6)
%	10	C in reactor feed (stream 6)
%	11	D in reactor feed (stream 6)
%	12	E in reactor feed (stream 6)
%	13	F in reactor feed (stream 6)
%	14	A in purge (stream 9)
%	15	B in purge (stream 9)
%	16	C in purge (stream 9)
%	17	D in purge (stream 9)
%	18	E in purge (stream 9)
%	19	F in purge (stream 9)
%	20	G in purge (stream 9)
%	21	H in purge (stream 9)
%	22	D in product (stream 11)
%	23	E in product (stream 11)
%	24	F in product (stream 11)
%	25	G in product (stream 11)
%	26	H in product (stream 11)
%	27	Volumetric production rate [m^3/h].
%	28	Recycle rate (stream 8) [kscmh]
%	29	Compressor power [kW]
%	30	Stripper temperature [C]
%	31	Separator temperature [C]
%	32	Condenser coolant outlet temp [C]
%	33	Cost [cents/kmol product]
%	34	Rate of reaction 1 [kmol G produced/h]
%	35	Rate of reaction 2 [kmol H produced/h]
%	36	Rate of reaction 3 [kmol F produced/h]
%	37	Partial pressure of A in reactor [kPa]
%	38	Partial pressure of C in reactor [kPa]
%	39	Partial pressure of D in reactor [kPa]
%	40	Partial pressure of E in reactor [kPa]

% Constant parameters:

AVP=[0;0;0;20.81;21.24;21.24;21.32;22.10];  		% A in Antoine eqn.
BVP=[0;0;0;-1444;-2114;-2114;-2748;-3318];  		% B in Antoine eqn.
CVP=[0;0;0;259;266;266;233;250];            		% C in Antoine eqn.
mwts=[2;25.4;28;32;46;48;62;76];            		% molecular wts.
molvol=[0;0;0;0.1070;0.1260;0.1463;0.1013;0.1231];  % Mol. volumes [m3/kmol]
CpG=[14.6;2.04;1.05;1.85;1.87;2.02;.712;.628];    	% vapor ht cap [kJ/kg-K]
CpL=[0;0;0;7.66;4.17;4.45;2.55;2.45];             	% liquid ht cap [kJ/kg-K]
dhv=[0;0;0;202;372;372;523;486];                  	% latent heats [kJ/kg]

% Initialize result vectors

u=u0;
y=y0;
x=zeros(31,1);

% Known temperatures and pressures

Pr=y(1)+101;       % Reactor pressure [kPa]
Ps=y(3)+101;       % Separator pressure [kPa]
Pf=y(6)+101;       % Stripper/feed pressure [kPa]

Tcr=u(10);         % Reactor temp [C]
Tkr=Tcr+273.2;     %              [K]
Tcs=y(31);         % Separator temp [C]
Tks=Tcs+273.2;     %                [K]
Tkf=86.1+273.2;    % Feed zone temp [K]
Tcc=y(30);         % Stripper temp [C]

% Vapor pressures in separator, reactor, and stripper

Pvs=[zeros(3,1);
    0.001*exp(AVP(4:8)+(BVP(4:8)./(CVP(4:8)+Tcs)))];	
Pvr=[zeros(3,1);
    0.001*exp(AVP(4:8)+(BVP(4:8)./(CVP(4:8)+Tcr)))];	
Pvc=[zeros(3,1);
    0.001*exp(AVP(4:8)+(BVP(4:8)./(CVP(4:8)+Tcc)))];

% Set specified flows [kmol/h].  Also mole fractions
% of feeds and their molar flows [kmol/h].

F1=u(1);  x1=[1;zeros(7,1)];     Fi1=x1*F1;   % Pure A
F2=u(2);  x2=[0;0;0;1;0;0;0;0];  Fi2=x2*F2;   % Pure D
F3=u(3);  x3=[0;0;0;0;1;0;0;0];  Fi3=x3*F3;   % Pure E
F4=u(4);                                      % A, B, C
x4=[u(12:13)/100;1-sum(u(12:13)/100);zeros(5,1)];
Fi4=x4*F4;
F8=y(28)*44.79;
F9=u(6);
F10=u(7);
F11=u(8);

% Product mole fractions (purge and liquid product)
% Normalize them to guarantee sum of mole fractions is 1.0.
% Note that product composition is based on assumption of
% zero A, B, and C.

x9=0.01*y(14:21);
x9=x9/sum(x9);   
x8=x9;            					% Recycle = purge
PurgMwt=sum(x8.*mwts);				% Mol Wt of purge
x11=[0;0;0;y(22:26)]/sum(y(22:26));
Fi11=x11*F11;

% Compute molar rate of each component leaving process.
% Also get molar rates of each component in purge and recycle

Fi8=x8*F8;
Fi9=x9*F9;
out_i=Fi9+Fi11;

% Use known stream 11 to calculate G and H in streams 5 and 10.

Fi5=zeros(8,1);
Fi10=zeros(8,1);
S=(Pvc/Pf)*F4/F10;                     % Stripping factors, Ki*V/L.
Fi10(7:8)=Fi11(7:8)./(1-(S(7:8).*[0.809;0.720]));
Fi5(7:8)=Fi10(7:8)-Fi11(7:8);

% Partial pressures in separator

Pis=x9*Ps;

% Use the known ratio of xH/xG in stream 10 and the separator
% partial and vapor pressures to get the VLE adjustment factors

xH2xG=Fi10(8)/Fi10(7);   % ratio of H to G in stream 10
Prat=Pis(4:8)./Pvs(4:8); % ratio of partial pressure to vapor pressure
gamG=sum(Prat(1:3))+(1+xH2xG)*Prat(4);
gamH=gamG*Prat(5)/Prat(4)/xH2xG;
u(20:21)=100*[gamG;gamH]; % VLE adjustment factors for separator

% Now we can get the mole fractions, etc., in stream 10:

x10=zeros(8,1);
x10(4:7)=Prat(1:4)/gamG;
x10(8)=Prat(5)/gamH;
if abs(sum(x10)-1) > 10*eps
   error('x10 does not sum to unity')
end
F10=Fi10(8)/x10(8);

molvols=sum(x10.*molvol);         % molar volume of liquid in sep. [m^3]
Fi10=x10*F10;
u(16)=u(7)-F10;                   % bias on stream 10 flow [kmol/h]

% Now that stream 10 is known, we can complete the stripper
% balance.

Fi5(1:3)=Fi4(1:3)-Fi11(1:3);
Fi5(4:6)=Fi10(4:6)-Fi11(4:6);

% We can also calculate the reactor effluent (stream 7),
% vapor mole fractions in reactor, and partial pressures in
% reactor.

Fi7=Fi8+Fi9+Fi10;     % Molar flows [kmol/h]
F7=sum(Fi7);          % Total flow [kmol/h]
x7=Fi7/F7;            % Mole fraction in stream 7
mwt7=sum(x7.*mwts);   % average molecular weight
Pir=Pr*x7;            % partial pressures [kPa]

% Get liquid composition in reactor.  Force VLE at the specified
% temperature and pressure by adjusting activity coefficients.

xr=zeros(8,1);
Prat=Pir(4:8)./Pvr(4:8);
gamDH=sum(Prat);
xr(4:8)=Prat/gamDH;            
u(26)=100*gamDH;       % Nonideal VLE adjustment.  Same for D to H.

% Calculate the flow coefficient that gives the correct F7 for
% the known pressure drop between the reactor and separator

if Pr <= Ps
   error('Reactor pressure <= Separator pressure')
end
u(17)=100*(F7/(5722/mwt7)/sqrt(Pr-Ps));

% Use known liquid holdup in reactor to get vapor volume.
% Need to calculate liquid molar volume.

molvolr=sum(xr.*molvol);    % Molar volume [m^3]
VLR=(y(2)+12.105)/5.263;    % Liq vol conversion from % to m^3
VVR=36.8-VLR;               % Vapor vol in reactor [m^3]

% Now use kinetic laws to get reaction rates [kmol/h]

RR(1)=VVR*exp(44.06-42600.0/1.987/Tkr)*(Pir(1)^1.08) ...
            *(Pir(3)^0.311)*(Pir(4)^0.874);
RR(2)=VVR*exp(10.27-19500.0/1.987/Tkr)*(Pir(1)^1.15) ...
            *(Pir(3)^0.370)*(Pir(5)^1.00);
RR(3)=VVR*exp(59.50-59500.0/1.987/Tkr)*Pir(1) ...
            *(0.77*Pir(4)+Pir(5));
		
% Calculate adjustment factors for reactions 1, 2, and 3 to make
% correct production of F, G and H

u(14)=100*out_i(7)/RR(1);
u(15)=100*out_i(8)/RR(2);
u(25)=100*out_i(6)/RR(3);
RR(1)=out_i(7);
RR(2)=out_i(8);
RR(3)=out_i(6);

% Use stoichiometry and known stream 7 to calculate reactor
% influent (stream 6)

nu1=[-1;0;-1;-1;0;0;1;0];             			% Stoich. coeffs. for reaction 1
nu2=[-1;0;-1;0;-1;0;0;1];              			% Stoich. coeffs. for reaction 2
nu3=[-0.333;0;0;-1;-0.333;1;0;0];      			% Stoich. coeffs. for reaction 2
genratei=nu1*RR(1) + nu2*RR(2) + nu3*RR(3);   	% generation rates [kmol/h]
Fi6=Fi7 - genratei ;  							% molar flows [kmol/h]
F6=sum(Fi6);                                  	% total flow [kmol/h]
x6=Fi6/F6;                                    	% mole fractions
mwt6=sum(x6.*mwts);                           	% average mol. weight

% Calculate flow coefficient between feed/stripper and reactor

if Pf <= Pr
   error('Feed pressure <= Reactor pressure')
end
u(18)=100*(F6/(2413.7/mwt6)/sqrt(Pf-Pr));

% Now calculate bias flows to satisfy the material balance
% at the feed point.  F, G and H should already be correct.

Fibias = Fi6 - Fi1 - Fi2 - Fi3 - Fi5 - Fi8;

if any(abs(Fibias([6:8])) > 1e-13)
   disp('WARNING:  Bad mass balance for F, G or H')
   fprintf(' F bias = %10.4e,  G bias = %10.4e,  H bias = %10.4e\n',...
   		Fibias(6),Fibias(7),Fibias(8))
end

% For the bias in A, recalculate the mole fraction of A in stream 4

Fi4(1)=Fi4(1)+Fibias(1);        % increase A by amount of bias
Fi4(3)=Fi4(3)-Fibias(1);        % reduce C to keep same total flow
Fibias(3)=Fibias(3)+Fibias(1);  % increase C bias accordingly
Fibias(1)=0;                    % this bias is now accounted for...

% Same idea for B

Fi4(2)=Fi4(2)+Fibias(2);        % increase A by amount of bias
Fi4(3)=Fi4(3)-Fibias(2);        % reduce C to keep same total flow
Fibias(3)=Fibias(3)+Fibias(2);  % increase C bias accordingly
Fibias(2)=0;                    % this bias is now accounted for...

Fi5(1:3)=Fi4(1:3)-Fi11(1:3);	% reset stream 5 A,B,C flows
F5=sum(Fi5);
x5=Fi5/F5;
x4=Fi4/F4;
u(12:13)=100*x4(1:2);           % corrected mole % A, B in stream 4

% Now set the bias parameters in the U vector

u(22:24)=Fibias(3:5);           % for C, D, and E

% Calculate the states.  First the reactor.  Holdup of A, B, and C
% sets the vapor moles.  That of D, E, F, G and H sets the liquid.
% Use ideal gas to calculate vapor moles, and liquid molar volume
% to get liquid.

Rg=8.314;                         % gas constant [kPa*m3/kmol/K]
x(1:3)=Pir(1:3)*VVR/(Rg*Tkr);     % A, B, and C in reactor [kmol]
NLr=VLR/molvolr;                  % kmol liquid in reactor
x(4:8)=NLr*xr(4:8);               % D to H in reactor [kmol]

% The separator uses the same approach

VLS=(y(4)+10.53)/12.28;           % Liquid in separator [m^3]
VVS=99.1-VLS;                     % Vapor in separator [m^3]
x(9:11)=Pis(1:3)*VVS/(Rg*Tks);    % A, B, and C in separator [kmol]
NLs=VLS/molvols;                  % kmol liquid in separator
x(12:16)=NLs*x10(4:8);            % D to H in separator [kmol]

% The feed zone is just an ideal gas

VVf=150;                          % Volume of feed zone
NVf=Pf*VVf/(Rg*Tkf);              % kmol gas in feed zone
x(17:24)=NVf*x6;                  % kmol A to H in feed zone

% Calculate the volume and holdup of G and H
% in the stripper bottoms

molvolp=sum(x11.*molvol);     		% molar volume [m^3/kmol]
VLp=(y(5)+49.03)/22.58;          	% liquid volume [m^3]
NLp=VLp/molvolp;                  	% kmol liquid
x(25:29)=NLp*x11(4:8);				% kmol of D-H 

% Set product density adjustment factor

Qcalc=F11*molvolp;
u(29)=100*y(27)/Qcalc;

% Reset the reactor feed composition and flow

y(8:13)=100*x6(1:6);
y(7)=F6/44.79;

% Reset reaction rates and reactant partial pressures

y(34:36)=RR(1:3)';
y(37:40)=Pir([1,3,4,5]);

% Set the parameters in the compressor model.

F8mass=F8*PurgMwt;					% mass flow of recycle [kg/h]
dP=max([0, y(6)-y(3)]);				% pressure rise across compressor
Mb=60.83*u(5)*sqrt(dP);				% back flow [kg/h]
Cflow=F8mass+Mb;					% total compressor throughput [kg/h]
ucCflow=-0.05152*dP*dP - 92.52*dP + 1.17e5;	% uncorrected throughput
u(27)=Cflow-ucCflow;						% Flow bias correction
molCflow=Cflow/PurgMwt;
CPnorm=y(29)*(3600/molCflow)/(Tks*Rg);
CPbias=CPnorm + 9.365e-8*dP*dP - 4.300e-4*dP + 8.667e-3;
u(28)=CPbias;

fprintf('Compressor flow  correction = %8.2f percent\n',100*u(27)/Cflow)
fprintf('Compressor power correction = %8.2f percent\n',100*u(28)/CPnorm)

% Calculate the condenser duty and related parameters.
% NOTE:  all enthalpy calculations assume a reference temperature
%        of 100 degrees C.

T4=45;                          % Assumption from base case [C]
H4=F4*sum(((T4-100)*CpG+dhv).*mwts.*x4);	 	% [kJ/h]
H5=F5*sum(((Tcc-100)*CpG+dhv).*mwts.*x5);	 	% [kJ/h]
H7=F7*sum(((Tcr-100)*CpG+dhv).*mwts.*x7);		% [kJ/h]
H89=(F8+F9)*sum(((Tcs-100)*CpG+dhv).*mwts.*x9);	% [kJ/h]
H10=F10*((Tcs-100)*sum(CpL.*mwts.*x10));		% [kJ/h]
H11=F11*((Tcc-100)*sum(CpL.*mwts.*x11));		% [kJ/h]

Qcon=(H89 + H10 - H7)/3600;        	% [kW]
Fw=2725*u(11);						% [kg/h]
Tw=y(32);
Tw0=Tw + Qcon*3600/(Fw*4.18);      	% condenser coolant inlet Temp [C].
UA=Qcon/(0.5*(Tw + Tw0) - Tcs);		% UA [kW/K]					 
cfact=Fw*4.18/(UA*3600);

u(33)=Tw0;
u(34)=UA;

% Stripper model parameters

A=ones(3,1)./S(4:6);
u(30:32)=log((1-A).*((F10*x10(4:6))./(x11(4:6)*F11)) + A)./log(S(4:6));

Qstrip=H11+H5-H4-H10;			% [kJ/h]
Qsteam=u(9)*2260;				% [kJ/h]
u(35)=(Qsteam-Qstrip)/3600;		% [kW]

% Energy states

x(30)=0.001*H10/F10;
x(31)=0.001*H11/F11;

% Balances

%overall_bal=Fi1+Fi2+Fi3+Fi4+Fibias-Fi9-Fi11+genratei
%strip_bal=Fi4+Fi10-Fi5-Fi11
%sep_bal=Fi7-Fi8-Fi9-Fi10
%feed_bal=Fi1+Fi2+Fi3+Fi5+Fi8-Fi6+Fibias
%reac_bal=Fi6-Fi7+genratei


tab=zeros(10,11);
tab(1,:)=[1:11];
tab(2:9,:)=[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11];
tab(10,:)=[F1 F2 F3 F4 F5 F6 F7 F8 F9 F10 F11];
comps=['A','B','C','D','E','F','G','H'];

fprintf('\n          1        2        3        4        5        6') 
for i=2:9
  fprintf(['\n',comps(i-1),':    %7.5f  %7.5f  %7.5f'],...
              tab(i,1),tab(i,2),tab(i,3))
  fprintf('  %7.5f  %7.5f  %7.5f',tab(i,4),tab(i,5),tab(i,6))
end
fprintf('\nkmol: %7.1f  %7.1f  %7.1f',tab(10,1),tab(10,2),tab(10,3))
fprintf('  %7.1f  %7.1f  %7.1f',tab(10,4),tab(10,5),tab(10,6))
fprintf('\n\n          7        8        9       10       11')
for i=2:9
  fprintf(['\n',comps(i-1),':    %7.5f  %7.5f  %7.5f'],...
            tab(i,7),tab(i,8),tab(i,9))
  fprintf('  %7.5f  %7.5f',tab(i,10),tab(i,11))
end
fprintf('\nkmol: %7.1f  %7.1f %7.1f',tab(10,7),tab(10,8),tab(10,9))
fprintf(' %7.1f  %7.1f\n\n',tab(10,10),tab(10,11))

% Get cost

purgls=0.01*F9*sum(y(14:21).*[221;0;618;2210;1460;1790;3040;2290]);	% Purge cost [$/h]
prodls=0.01*F11*sum(y(22:24).*[2210;1460;1790]);					% Product cost [$/h]
     	
y(33)=(5.36*y(29)+3.18*u(9)+purgls+prodls)/F11;  % [cents/kmol]
