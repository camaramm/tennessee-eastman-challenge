function [x,u,y]=TEest3_ss(y0,u0)

% Calculates steady state of estimator to match certain measured outputs and 
% specified inputs.  Outputs, states, and inputs are arranged as in TEest3.

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
%	25	G in product reservoir (stripper bottoms)
%	26	H in product reservoir

%	INPUTS are flows [kmol/h] unless noted otherwise:

%	 1	Feed 1 (pure A)
%	 2	Feed 2 (pure D)
%	 3	Feed 3 (pure E)
%	 4	Feed 4 (A & C)
%	 5	Recycle (stream 8)
%	 6	Purge (stream 9)
%	 7	Separator underflow (stream 10)
%	 8	Product rate (stream 11)
%	 9	Reactor temperature [deg C]
%	10	Separator temperature [deg C]
%	11	A in stream 4 [Mole %]
%	12	B in stream 4 [Mole %]
%	13	Reaction 1 activity factor [%].
%	14	Reaction 2 activity factor [%].
%	15	Stream 10 bias flow [kmol/h].
%	16	Reactor/Sep flow parameter [%].
%	17	Feed/Reactor flow parameter [%].
%	18	Product G+H purity parameter [%].
%	19	Adjustment to VLE of D to G in separator [%]
%	20	Adjustment to H VLE in separator [%]
%	21	C bias flow to feed zone [kmol/h]
%	22	D bias flow to feed zone [kmol/h]
%	23	E bias flow to feed zone [kmol/h]
%	24	F bias flow to feed zone [kmol/h]
%	25	Adjustment to VLE in reactor [%]

%	OUTPUTS are mole % in stream unless noted otherwise:

%	 1	Reactor pressure [kPa]
%	 2	Reactor liq. holdup [%]
%	 3	Separator pressure [kPa]
%	 4	Separator liq. holdup [%]
%	 5	Product liq. holdup [%]
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
%	22	G in product (stream 11)
%	23	H in product (stream 11)
%	24	Molar production rate (same as input 8) [kmol/h].
%	25	Cost [cents/kmol product]
%	26	Rate of reaction 1 [kmol G produced/h]
%	27	Rate of reaction 2 [kmol H produced/h]
%	28	Rate of reaction 3 [kmol F produced/h]
%	29	Partial pressure of A in reactor [kPa]
%	30	Partial pressure of C in reactor [kPa]
%	31	Partial pressure of D in reactor [kPa]
%	32	Partial pressure of E in reactor [kPa]

% Constant parameters:

AVP=[0;0;0;20.81;21.24;21.24;21.32;22.10];  % A in Antoine eqn.
BVP=[0;0;0;-1444;-2114;-2114;-2748;-3318];  % B in Antoine eqn.
CVP=[0;0;0;259;266;266;233;250];            % C in Antoine eqn.
mwts=[2;25.4;28;32;46;48;62;76];            % molecular wts.
molvol=[0;0;0;0.1070;0.1260;0.1463;0.1013;0.1231];  % Mol. volumes [m3/kmol]
SF=[ones(6,1);0.07;0.04];                   % Stripping factors

% Initialize result vectors

u=u0;
y=y0;
x=zeros(26,1);

% Set specified flows [kmol/h].  Also mole fractions
% of feeds and their molar flows [kmol/h].

F1=u(1);  x1=[1;zeros(7,1)];     Fi1=x1*F1;   % Pure A
F2=u(2);  x2=[0;0;0;1;0;0;0;0];  Fi2=x2*F2;   % Pure D
F3=u(3);  x3=[0;0;0;0;1;0;0;0];  Fi3=x3*F3;   % Pure E
F4=u(4);                                      % A, B, C
x4=[u(11:12)/100;1-sum(u(11:12)/100);zeros(5,1)];
Fi4=x4*F4;
F8=u(5);
F9=u(6);
F10=u(7);
F11=u(8);

% Product mole fractions (purge and main product)
% Normalize x9 to guarantee sum of mole fractions is 1.0

x9=0.01*y(14:21);
x9=x9/sum(x9);   
x8=x9;            % Recycle = purge
x11=0.01*[zeros(6,1);y(22:23)];
Fi11=x11*F11;

% Product purity parameter

u(18)=100*sum(x11);

% Compute molar rate of each component leaving process.
% Also get molar rates of each component in purge and recycle

out_i=F9*x9+F11*x11;
Fi9=x9*F9;
Fi8=x8*F8;

% Use known stream 11 to calculate G and H in streams 5 and 10

Fi5=zeros(8,1);
Fi10=zeros(8,1);
Fi10(7:8)=Fi11(7:8)./(1-SF(7:8));
Fi5(7:8)=Fi10(7:8)-Fi11(7:8);

% Known temperatures and pressures

Pr=y(1)+101;       % Reactor pressure [kPa]
Ps=y(3)+101;       % Separator pressure [kPa]
Pf=y(6)+101;       % Stripper/feed pressure [kPa]

Tcr=u(9);          % Reactor temp [C]
Tkr=Tcr+273.2;     %              [K]
Tcs=u(10);         % Separator temp [C]
Tks=Tcs+273.2;     %                [K]
Tkf=86.1+273.2;    % Feed zone temp [K]

% Vapor pressures in separator and reactor

Pvs=[zeros(3,1);
    0.001*exp(AVP(4:8)+(BVP(4:8)./(CVP(4:8)+Tcs)))];	
Pvr=[zeros(3,1);
    0.001*exp(AVP(4:8)+(BVP(4:8)./(CVP(4:8)+Tcr)))];	

% Partial pressures in separator

Pis=x9*Ps;

% Use the desired ratio of xH/xG in stream 10 and the separator
% partial and vapor pressures to get the VLE adjustment factors

xH2xG=Fi10(8)/Fi10(7);   % ratio of H to G in stream 10
Prat=Pis(4:8)./Pvs(4:8); % ratio of partial pressure to vapor pressure
gamG=sum(Prat(1:3))+(1+xH2xG)*Prat(4);
gamH=gamG*Prat(5)/Prat(4)/xH2xG;
u(19:20)=100*[gamG;gamH]; % VLE adjustment factors for separator

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
u(15)=u(7)-F10;                   % bias on stream 10 flow [kmol/h]

% Now that stream 10 is known, we can complete the stripper
% balance, assuming that all the D, E, and F
% in stream 10 go into stream 5.

Fi5(1:3)=Fi4(1:3);
Fi5(4:6)=Fi10(4:6);

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
u(25)=100*gamDH;       % Nonideal VLE adjustment.  Same for D to H.

% Calculate the flow coefficient that gives the correct F7 for
% the known pressure drop between the reactor and separator

if Pr <= Ps
   error('Reactor pressure <= Separator pressure')
end
u(16)=100*(F7/(5722/mwt7)/sqrt(Pr-Ps));

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
		
% Calculate adjustment factors for reactions 1 and 2 to make
% correct production of G and H

u(13)=100*out_i(7)/RR(1);
u(14)=100*out_i(8)/RR(2);
RR(1)=out_i(7);
RR(2)=out_i(8);

% Use stoichiometry and known stream 7 to calculate reactor
% influent (stream 6)

nu1=[-1;0;-1;-1;0;0;1;0];              % Stoich. coeffs. for reaction 1
nu2=[-1;0;-1;0;-1;0;0;1];              % Stoich. coeffs. for reaction 2
nu3=[-0.333;0;0;-1;-0.333;1;0;0];      % Stoich. coeffs. for reaction 2
genratei=nu1*RR(1) + nu2*RR(2) + nu3*RR(3);   % generation rates [kmol/h]
Fi6=Fi7 - genratei ;  % molar flows [kmol/h]
F6=sum(Fi6);                                  % total flow [kmol/h]
x6=Fi6/F6;                                    % mole fractions
mwt6=sum(x6.*mwts);                           % average mol. weight

% Calculate flow coefficient between feed/stripper and reactor

if Pf <= Pr
   error('Feed pressure <= Reactor pressure')
end
u(17)=100*(F6/(2413.7/mwt6)/sqrt(Pf-Pr));

% Now calculate bias flows to satisfy the material balance
% at the feed point.  G and H should already be correct.

Fibias = Fi6 - Fi1 - Fi2 - Fi3 - Fi5 - Fi8;

if any(abs(Fibias([7,8])) > 100*eps)
   disp('WARNING:  Bad mass balance for G or H')
   disp(['Biases are:  G = ',num2str(Fibias(7)),'   H = ',...
         num2str(Fibias(8))])
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

Fi5(1:3)=Fi4(1:3);              % reset stream 5 A,B,C flows
F5=sum(Fi5);
x5=Fi5/F5;
x4=Fi4/F4;
u(11:12)=100*x4(1:2);           % corrected mole % A, B in stream 4

% Now set the bias parameters in the U vector

u(21:24)=Fibias(3:6);           % for C, D, E, and F

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

% Finally, calculate the volume and holdup of G and H
% in the stripper bottoms

molvolp=sum(x11.*molvol);         % molar volume [m^3]
VLp=(y(5)+49.03)/22.58;           % liquid volume [m^3]
NLp=VLp/molvolp;                  % kmol liquid
x(25:26)=NLp*x11(7:8);            % kmol of G and H 

% Reset the reactor feed composition and flow

y(8:13)=100*x6(1:6);
y(7)=F6/44.79;

% Reset reaction rates and reactant partial pressures

y(26:28)=RR(1:3)';
y(29:32)=Pir([1,3,4,5]);

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
fprintf('  %7.1f  %7.1f  %7.1f',tab(10,4),tab(10,5),tab(10,7))
fprintf('\n\n          7        8        9       10       11')
for i=2:9
  fprintf(['\n',comps(i-1),':    %7.5f  %7.5f  %7.5f'],...
            tab(i,7),tab(i,8),tab(i,9))
  fprintf('  %7.5f  %7.5f',tab(i,10),tab(i,11))
end
fprintf('\nkmol: %7.1f  %7.1f %7.1f',tab(10,7),tab(10,8),tab(10,9))
fprintf(' %7.1f  %7.1f\n\n',tab(10,10),tab(10,11))



