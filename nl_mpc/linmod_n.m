function [A,B,C,D]=linmod_nlr(fun,x,u,para,xpert,upert,...
                              coeff1,coeff2,coeff3,coeff4,coeff5)   
% Modification of LINMOD to check whether any upert elements are zero.
% if so, skips those inputs.
%
%LINMOD Obtains linear models from systems of ord. diff. equations (ODEs).
%	[A,B,C,D]=LINMOD('SFUNC') obtains the state-space linear model of the 
%	system of ordinary differential equations described in the 
%       S-function 'SFUNC' when  the state variables and inputs are set
%       to zero. SFUNC may be a SIMULINK or other model. See, SFUNC.
%
%	[A,B,C,D]=LINMOD('SFUNC',X,U) allows the state vector, X, and 
%	input, U, to be specified. A linear model will then be obtained
%	at this operating point. 
%
%	[A,B,C,D]=LINMOD('SFUNC',X,U,PARA) allows a vector of parameters
%	to be set.  PARA(1) sets the perturbation level for obtaining the
%	linear model (default PARA(1)=1e-5).  For systems that are functions
%	of time PARA(2) may be set with the value of t at which the linear 
%	model is to be obtained (default PARA(2)=0).

%	Copyright (c) 1990-1992 by the MathWorks, Inc.
%	Andrew Grace 11-12-90.
%	Revised ACWG 3-9-91 

%	[A,B,C,D]=LINMOD('SFUNC',X,U,PARA,XPERT,UPERT) allows the perturbation
%	levels for all of the elements of X and U to be set. 
%	The default is otherwise  XPERT=PARA(1)+1e-3*PARA(1)*abs(X),
%	UPERT=PARA(1)+1e-3*PARA(1)*abs(U).
%
%	[A,B,C,D]=LINMOD('SFUNC',X,U,PARA,XPERT,UPERT,COEFF) allows 
%	coefficients, COEFF, to be passed directly to SFUNC at each call from 
%	LINMOD: [X,T]=SFUNC(X,T,U,FLAG,COEFF).
%
%	Any or all of  PARA, XPERT, UPERT may be empty matrices in which case
%	these parameters will be assumed to be undefined and the default 
%	option will be used.
 

% ---------Extra paramaters for S-function systems ------
outstr=[];
for i=1:nargin-6
	outstr=[outstr,',coeff',int2str(i)];
end
% ---------------Options--------------------
ifun = [fun,'([],[],[],0', outstr,')'];
sizes = eval(ifun);
sizes=[sizes(:); zeros(6-length(sizes),1)];
nxz=sizes(1)+sizes(2); nu=sizes(4); ny=sizes(3);
nx=sizes(1);

if nargin<2, x=[]; end
if nargin<3, u=[]; end
if nargin<4, para=[]; end
if nargin<5, xpert=[]; end
if nargin<6, upert=[]; end

if isempty(u), u=zeros(nu,1); end
if isempty(x), x=zeros(nxz,1); end
if isempty(para) , para=[0;0]; end
if para(1)==0, para(1)=1e-5; end
if isempty(upert), upert=para(1)+1e-3*para(1)*abs(u); end
if isempty(xpert), xpert=para(1)+1e-3*para(1)*abs(x); end
if length(para)>1, t=para(2); else t=0; end
if length(x)<nxz 
	disp('Warning: Extra states being set to zero.')
	x=[x(:); zeros(nxz-length(x),1)];
end

if nxz > nx 
	disp('Warning: Extra discrete states ignored (see dlinmod for conversion).');
end

fun1=[fun,'(t, x,u,1', outstr, ')'];
fun2=[fun, '(t, x,u,3', outstr, ')'];

% Check for zero upert elements

iu=find(upert);
nu=length(iu);
A=zeros(nx,nx); B=zeros(nx,nu); C=zeros(ny,nx); D=zeros(ny,nu);

% Initialization
y=0;
dx = eval(fun1);
if ny > 0, y = eval(fun2); end
olddx=dx; oldu=u; oldy=y; oldx=x;

% A amd C matrices
for i=1:nx;
	x(i)=x(i)+xpert(i);
	dx = eval(fun1);
 	A(:,i)=(dx-olddx)./xpert(i);
	if ny > 0
		y = eval(fun2);
		C(:,i)=(y-oldy)./xpert(i);
	end
	x=oldx;
end

% B and D matrices
for i=1:nu
	  u(iu(i))=u(iu(i))+upert(iu(i));
	  dx = eval(fun1);
	  B(:,i)=(dx-olddx)/upert(iu(i));
	  if ny > 0
		    y = eval(fun2); 
		    D(:,i)=(y-oldy)/upert(iu(i));
	  end
	  u=oldu;
end

