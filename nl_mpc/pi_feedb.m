function [a,b,c,d]=PI_feedback(a1,b1,c1,d1,K,Tau,icv,imv)

% Connect diagonal PI feedback controller to specified inputs 
% & outputs of plant in LTI (continuous) form.
%
%[a,b,c,d]=PI_feedback(a1,b1,c1,d1,K,Tau,icv,imv)
%
%a1, b1, c1, d1:  state-space model of plant
%icv:   index of plant outputs to be controlled
%imv:   index of plant inputs to use as manipulated variables
%K, Tau:  vectors of gains and integral times for each loop in icv/imv
%         pairing.
%a, b, c, d:  resulting closed-loop system.  Inputs imv are replaced
%             by setpoint inputs.

% error checks

nl=length(icv);
if length(K) ~= nl
   error('K and icv must be equal length')
elseif length(Tau) ~= nl
   error('Tau and icv must be equal length')
elseif length(imv) ~= nl
   error('imv and icv must be equal length')
end

error(abcdchkm(a1,b1,c1,d1))

[ny1,n1]=size(c1);
[n1,nu1]=size(b1);

if min(icv) < 1 | max(icv) > ny1
   error('icv has element out of range')
elseif min(imv) < 1 | max(imv) > nu1
   error('imv has element out of range')
end

% form the controller state-space matrices (diagonal PI control)

a2=zeros(nl,nl);
b2=eye(nl);
c2=diag(K./Tau);
d2=diag(K);

% set up pointers to unconnected inputs and outputs

iuu=[1:nu1];
iyu=[1:ny1];

iuu(imv)=[];        % unconnected inputs (not used for feedback)
nuu=length(iuu);
iyu(icv)=[];        % unconnected outputs (not used for feedback)
nyu=length(iyu);

% re-order state-space matrices

b1=b1(:,[imv,iuu]);
c1=c1([icv,iyu],:);
d1=d1([icv,iyu],[imv,iuu]);

% pointers to states, etc. in combined, re-ordered system

ix1=[1:n1];        % plant states
ix2=[n1+1:n1+nl];  % controller states
iu1=[1:nl];        % feedback inputs
iu2=[nl+1:nl+nuu]; % unconnected inputs
iy1=[1:nl];        % feedback outputs
iy2=[nl+1:nl+nyu]; % unconnected outputs

% Pick out desired submatrices

B1c=b1(:,iu1);
D1c=d1(:,iu1);
B1u=b1(:,iu2);
D1u=d1(:,iu2);

% Algebra to form combined system

t=eye(ny1)+[D1c*d2  zeros(ny1,nyu)];

c=t\[c1 D1c*c2];
d=t\[D1c*d2 D1u];

alph=[-d2*c(iy1,ix1)   c2-d2*c(iy1,ix2)];
bet=[(d2*(eye(nl)-d(iy1,iu1)))  -d2*d(iy1,iu2)];

a=[ (a1+B1c*alph(:,ix1))   B1c*alph(:,ix2)
       -b2*c(iy1,ix1)       (a2-b2*c(iy1,ix2))];

b=[  B1c*bet(:,iu1)            (B1c*bet(:,iu2)+B1u)
     b2*(eye(nl)-d(iy1,iu1))   -b2*d(iy1,iu2)];

% Put matrices in original input/output order

b(1:n1+nl,[imv,iuu])=b;
c([icv,iyu],1:n1+nl)=c;
d([icv,iyu],[imv,iuu])=d;





