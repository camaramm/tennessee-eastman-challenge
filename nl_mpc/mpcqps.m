function deltau=mpcqps(nc,TAB,a,rhscon,iumin,iumax,iymin,...
                 iymax,ymin_wt,ymax_wt,dumax,SuTQ,r,y0,u0)

%deltau=mpcqps(nc,TAB,a,rhscon,iumin,iumax,iymin,...
%                 iymax,ymin_wt,ymax_wt,dumax,SuTQ,r,y0,u0)
%
% nc, TAB, a, ..., SuTQ should have been calculated by the smpcqp0s function.
% r  is a setpoint projection vector.  Length must be number of outputs times
%    length of the prediction horizon.  Must be a column vector.
% y0 is an output projection vector (i.e., what the outputs would do if there
%    were no change in the manipulated variables.  Same size as r.
% u0 is column vector current manipulated variables.  
%
% deltau will contain the "optimal" changes in the manipulated variables.

[nui,dum]=size(u0);
[mnu,pny]=size(SuTQ);
[tabrow,tabcol]=size(TAB);
nvars=tabrow-nc;
M=mnu/nui;
del=reshape(u0(:,ones(1,M)),mnu,1); % creates "stacked" u0 repeated M times

% Calculate starting basis vector for the QP

rhsa=a+SuTQ*(r-y0);
if nvars > mnu
   rhsa=[rhsa;0];   % adjusts for soft output constraints
end

% Update the RHS of the inequality constraints

rhsc=zeros(mnu,1);           

if ~ isempty(iumin)    % Equations for lower bound on u
   rhsc=[rhsc;del(iumin,:)];
end
if ~ isempty(iumax)    % Equations for upper bound on u
   rhsc=[rhsc;-del(iumax,:)];
end
if ~ isempty(iymin)    % Equations for lower bound on y
   rhsc=[rhsc;ymin_wt*y0(iymin,:)];
end
if ~ isempty(iymax)    % Equations for upper bound on y
   rhsc=[rhsc;-ymax_wt*y0(iymax,:)];
end

rhsc=rhsc+rhscon;      % Add on the constant part computed earlier.

% Set up and solve the QP;

basisi=[    -TAB(1:nvars,1:nvars)*rhsa
           rhsc-TAB(nvars+1:nvars+nc,1:nvars)*rhsa];
ibi=-[1:nvars+nc]';
ili=-ibi;
[basis,ib,il,iter]=dantzgmp(TAB,basisi,ibi,ili);
if iter < 0
   error('Infeasible QP.  Check constraints.');
end
deltau=-dumax(1:nui,1);
j=find(il([1:nui],1) > 0);
for j=1:nui
   if il(j) > 0
     deltau(j,1)=basis(il(j))-dumax(j,1);
   end
end
