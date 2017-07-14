function [nc,TAB,a,rhscon,iumin,iumax,iymin,iymax,dumax,SuTQ]=...
         smpcqp0(imod,ywt,uwt,blocks,p,ulim,ylim)

%[nc,TAB,a,rhscon,iumin,iumax,iymin,iymax,dumax,SuTQ]=...
%         smpcqp0(imod,ywt,uwt,blocks,p,ulim,ylim);
%
%Designs an MPC-type controller for constrained problems.
%Uses a state-space model of the process (MOD format).
%
%Output is intended as input for the MPCQP function.
%
%imod:     proces model in MOD format.
%ywt:      Penalty weighting for setpoint tracking.
%uwt:      Penalty weighting for changes in manipulated variables.
%M:        Number of moves OR sequence of blocking factors.
%P:        Length of prediction horizon.
%ulim:     [ Ulow Uhigh delU].  NOTE:  delU must be finite.
%ylim:     [ Ylow Yhigh].

% Copyright by N. L. Ricker, 12/91.

% +++ Check input arguments for errors and set default values. +++

if nargin == 0
   disp('USAGE: [nc,TAB,a,rhscon,iumin,iumax,iymin,iymax,dumax,SuTQ]=...')
   disp('                   smpcqp0(imod,ywt,uwt,blocks,p,ulim,ylim)')
   return
elseif nargin ~= 7
   error('Incorrect number of input arguments')
end

if isempty (imod)
   error('Process model must be supplied.')
end

% Get internal model specifications.

[phii,gami,ci,di,minfoi]=mod2ss(imod);
tsamp=minfoi(1);
ni=minfoi(2);
nui=minfoi(3);
nvi=minfoi(4);
mi=nui+nvi;
nwi=minfoi(5);
if nwi > 0             % Strip off unmeasured disturbance inputs if they exist.
   gami=gami(:,1:mi);
   di=di(:,1:mi);
end
nymi=minfoi(6);
nyui=minfoi(7);
nyi=nymi+nyui;

% Check for errors and inconsistencies in the models.

if any(any(di(:,1:nui)))
   error(['IMOD:  first nu=',int2str(nui),' columns of D must be zero.'])
end

if isempty(p)
   p=1;
elseif p < 1
   error('Prediction horizon is less than 1')
end

if isempty(ywt)
   ywt=ones(1,nyi);
   nywt=1;
else
   [nywt,ncol]=size(ywt);
   if ncol ~= nyi | nywt <= 0
      error('YWT is wrong size')
   end
   if any(any(ywt < 0))
      error('One or more elements of YWT are negative')
   end
end

if isempty(uwt),
   uwt=zeros(1,nui);
   nuwt=1;
else
   [nuwt,ncol]=size(uwt);
   if ncol ~= nui | nuwt <= 0
      error('UWT is wrong size')
   end
   if any(any(uwt < 0))
      error('UWT is negative')
   end
end

if isempty(blocks)
   blocks=ones(1,p);
   nb=p;
else
   [nrow,nb]=size(blocks);
   if nrow ~= 1 | nb < 1 | nb > p
      error('M vector is wrong size')
   end
   if any(blocks < 1)
      error('M contains an element that is < 1')
   end

   if nb == 1

%  This section interprets "blocks" as a number of moves, each
%  of one sampling period duration.

      if blocks > p
         disp('WARNING: M > P.  Truncated.')
         nb=p;
      elseif blocks <= 0
         disp('WARNING: M <= 0.  Set = 1.')
         nb=1;
      else
         nb=blocks;
      end
      blocks=[ones(1,nb-1) p-nb+1];

   else

% This section interprets "blocks" as a vector of blocking factors.

      sumblocks=sum(blocks);
      if sumblocks > p
               disp('WARNING:  sum(M) > P.')
               disp('          Moves will be truncated at P.')
               nb=find(cumsum(blocks) > p);
		   nb=nb(1);
               blocks=blocks(1,1:nb);
      elseif sumblocks < p
         nb=nb+1;
         blocks(nb)=p-sumblocks;
         disp('WARNING:  sum(M) < P.  Will extend to P.')
      end
   end
end

% Check the constraint specifications.  First set up some indices to pick out
% certain columns of the ulim and ylim matrices.

iumin=[1:nui];     % Points to columns of ulim containing umin.
iumax=iumin+nui;   % Points to columns of ulim containing umax.
idumax=iumax+nui;  % Points to columns of ulim containing delta u max.
iymin=[1:nyi];     % Points to columns of ylim containing ymin.
iymax=iymin+nyi;   % Points to columns of ylim containing ymax.

% Now check the values supplied by the user for consistency.

[nulim,ncol]=size(ulim);
if ncol ~= 3*nui | nulim <= 0
   error('ULIM matrix is empty or wrong size.')
elseif any(any(ulim(:,idumax) < 0))
   error('A constraint on DELTA U was < 0')
elseif any(any(ulim(:,iumax) < ulim(:,iumin)))
   error('A lower bound on U was greater than its upper bound')
end

% When using the DANTZGMP routine for the QP problem, we must have all
% bounds on delta u finite.  A bound that is finite but large can cause
% numerical problems.  The following loop checks for this.

ichk=0;
for i=idumax
   ifound=find(ulim(:,i) > 1e6);
   if ~ isempty(ifound)
      ichk=1;
      ulim(ifound,i)=1e6*ones(length(ifound),1);
   end
end
if ichk      
   disp('One or more constraints on delta_u were > 1e6.')
   disp('Reduced to 1e6 to prevent numerical problems in QP.')
end

if isempty(ylim)
   ylim=[-inf*ones(1,nyi) inf*ones(1,nyi)];
else
   [nylim,ncol]=size(ylim);
   if ncol ~= 2*nyi | nylim <= 0
      error('YLIM matrix is wrong size')
   elseif any(any(ylim(:,iymax) < ylim(:,iymin) ))
      error('A lower bound on y was greater than its upper bound')
   end 
end

% ++++ Beginning of controller design calculations. ++++

% The following index vectors are used to pick out certain columns
% or rows in the state-space matrices.

iu=[1:nui];         % columns of gami, gamp, di, dp related to delta u.
iv=[nui+1:nui+nvi];  % points to columns for meas. dist. in gamma.
iym=[1:nymi];       % index of the measured outputs.

% +++ Augment the internal model state with the outputs.

[PHI,GAM,C,D,N]=mpcaugss(phii,gami,ci,di);

% +++ Calculate the basic projection matrices +++

pny=nyi*p;          % Total # of rows in the final projection matrices.
mnu=nb*nui;         % Total # of columns in final Su matrix.

Cphi=C*PHI;
Su=[    C*GAM(:,iu)
    zeros(pny-nyi,nui)];

r1=nyi+1;
r2=2*nyi;
for i=2:p
   Su(r1:r2,:)=Cphi*GAM(:,iu);
   Cphi=Cphi*PHI;
   r1=r1+nyi;
   r2=r2+nyi;
end

%disp('Step Response'),disp(Su)  % Debugging output

Sdel=eye(nui); % Sdel is to be a block-lower-triangular matrix in which each
               % block is an identity matrix.  Used in constraint definition.
for i=2:nb
   Sdel=[Sdel;eye(nui)];
end

% If number of moves > 1, fill the remaining columns of Su and Sdel, 
% doing "blocking" at the same time.

if nb > 1
   k = nui;
   blocks=cumsum(blocks);
   for i = 2:nb
      row0=blocks(i-1)*nyi;
      row1=(i-1)*nui;
      Su(row0+1:pny,k+1:k+nui)=Su(1:pny-row0,1:nui);
      Sdel(row1+1:mnu,k+1:k+nui)=Sdel(1:mnu-row1,1:nui);
      k=k+nui;
   end
end

% Set up weighting matrix on outputs.  Q is a column vector
% containing the diagonal elements of the weighting matrix, SQUARED.

irow=0;
for i=1:p
   Q(irow+1:irow+nyi,1)=ywt(min(i,nywt),:)';
   irow=irow+nyi;
end
Q=Q.*Q;

% Set up weighting matrix on manipulated variables.  R
% is a column vector containing the diagonal elements, SQUARED.

uwt=uwt+10*sqrt(eps);  %for numerical stability
irow=0;
for i=1:nb
   R(irow+1:irow+nui,1)=uwt(min(i,nuwt),:)';
   irow=irow+nui;
end
R=R.*R;

%A=[diag(Q)*Su; diag(R)];disp('A:'),disp(A)  % Debugging output

% Usually, some of the general inequality constraints are not used.
% This section sets up index vectors for each type of constraint to
% pick out the ones that are actually needed for the problem.  This
% helps to minimize the size of the QP.

% First set up column vectors containing the bounds for each type of
% constraint over the entire prediction horizon.  For the inputs, the
% resulting vectors must be length mnu.  For outputs, length is pny.

umin=ulim(:,iumin)';
umin=umin(:);         % Stetches the matrix out into one long column
umax=ulim(:,iumax)';
umax=umax(:);
dumax=ulim(:,idumax)';
dumax=dumax(:);
ymin=ylim(:,iymin)';
ymin=ymin(:);
ymax=ylim(:,iymax)';
ymax=ymax(:);
clear ulim ylim       % Releases memory no longer needed.

lenu=length(umin);
if lenu > mnu         % Has user specified more bounds than necessary?
   disp('WARNING:  too many rows in ULIM matrix.')
   disp('          Extra rows deleted.')
   umin=umin(1:mnu);
   umax=umax(1:mnu);
   dumax=dumax(1:mnu);
elseif lenu < mnu     % If fewer rows than needed, must copy last one.
   r2=[lenu-nui+1:lenu];
   for i=1:round((mnu-lenu)/nui)
      umin=[umin;umin(r2,:)];
      umax=[umax;umax(r2,:)];
      dumax=[dumax;dumax(r2,:)];
   end
end

leny=length(ymin);
if leny > pny         % Has user specified more bounds than necessary?
   disp('WARNING:  too many rows in YLIM matrix.')
   disp('          Extra rows deleted.')
   ymin=ymin(1:pny);
   ymax=ymax(1:pny);
elseif leny < pny     % If fewer rows than needed, must copy last one.
   r2=[leny-nyi+1:leny];
   for i=1:round((pny-leny)/nyi)
      ymin=[ymin;ymin(r2,:)];
      ymax=[ymax;ymax(r2,:)];
   end
end

% The bounds on delta u must always be included in the problem.  The
% other bounds should only be included as constraints if they're finite.
% Generate vectors that contain a list of the finite constraints.

iumin=find(umin ~= -inf);
iumax=find(umax ~=  inf);
iymin=find(ymin ~= -inf);
iymax=find(ymax ~=  inf);

% Delete the infinite values.  At the same time, form the coefficient
% matrix for the inequality constraints.  Do this by picking out only
% the equations actually needed according to the lists established above.
% Finally, calculate the constant part of the RHS of the inequality
% constraints for these equations.

A=eye(mnu);        % These are the equations that are always present.
rhscon=2*dumax;    % They are the bounds on delta u.  A is the coefficient
                   % matrix and rhscon is the constant part of the RHS.

if ~ isempty(iumin)    % Add equations for lower bound on u
   umin=umin(iumin);
   A=[A;-Sdel(iumin,:)];
   rhscon=[rhscon;-Sdel(iumin,:)*dumax-umin];
else
   umin=[];
end
if ~ isempty(iumax)    % Add equations for upper bound on u
   umax=umax(iumax);
   A=[A;Sdel(iumax,:)];
   rhscon=[rhscon;Sdel(iumax,:)*dumax+umax];
else
   umax=[];
end
if ~ isempty(iymin)    % Add equations for lower bound on y
   ymin=ymin(iymin);
   A=[A;-Su(iymin,:)];
   rhscon=[rhscon;-Su(iymin,:)*dumax-ymin];
else
   ymin=[];
end
if ~ isempty(iymax)    % Add equations for upper bound on y
   ymax=ymax(iymax);
   A=[A;Su(iymax,:)];
   rhscon=[rhscon;Su(iymax,:)*dumax+ymax];
else
   ymax=[];
end

[nc,dumdum]=size(A);   % Save total number of inequality constraints.

% +++ Define the matrices needed for the QP +++


SuTQ=Su'*diag(Q);
B=SuTQ*Su+diag(R);
clear Su
a=B'*dumax;   % This is a constant term that adds to the initial basis
              % in each QP.
B=B\eye(mnu);
TAB=[-B   B*A' ;A*B   -A*B*A'];

% disp('C:'),disp(A)  % Debugging output
