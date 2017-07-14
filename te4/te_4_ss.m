function [x,VL]=te_4_ss(ya3,F4,P,p1)

% [x,VL]=te_4_ss(ya3,F4,P,p1)
%
% Finds a steady-state of the 2-phase reactor for
% the following specified conditions:
%
% ya3   mol fraction of A in purge
% F4    product rate (kmol/h)
% P     pressure (kPa)
% p1    parameter vector.  See TE_4.f for details.
%
% where x and VL are the calculated steady-state
% values of the states and the liquid volume (m^3),
% respectively.

% copyright N. L. Ricker
% University of Washington
% Chemical Engineering
% Box 351750
% Seattle, WA 98195-1750
% ricker@cheme.washington.edu

% Get yc3 from production rate

Pa=P*ya3;
Pc=(F4/(p1(8)*(Pa^1.2)))^(1/p1(9));
yc3=Pc/P;
Pb=P-Pa-Pc;
yb3=Pb/P;

ya1=p1(2);
yb1=p1(3);
yc1=1-ya1-yb1;

% Solve mat. bal. eqns. for flowrates

A=[ya1 1 -ya3
   yb1 0 -yb3
   yc1 0 -yc3];
b=[F4
    0
   F4];

F=A\b;

Cv=[3.3046, 0.2246, 0.00352, 0.0417];

x(5,1)=F(1)/Cv(1);
x(6,1)=F(2)/Cv(2);
x(7,1)=F(3)/(Cv(3)*sqrt(P-100));
x(8,1)=F4/(Cv(4)*sqrt(P-100));

Nd=110;
VL=Nd/8.3;
x(4,1)=Nd;
VV=122-VL;
NV=P*VV/(8.314*373);
x(1,1)=NV*ya3;
x(2,1)=NV*yb3;
x(3,1)=NV*yc3;
