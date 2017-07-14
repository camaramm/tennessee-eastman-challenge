c	Process model for the system described in "MPC of a continuous,
c	nonlinear, 2-phase reactor" by N. L. Ricker, J. of Process Control,
c	Vol. 3, pp. 109-123(1993).  Please see this paper for modeling details.

c	Copyright by N. L. Ricker, Professor
c	University of Washington
c	Department of Chemical Engineering
c	Box 351750
c	Seattle, WA 98195-1750

c	E-mail:  ricker@cheme.washington.edu

C	This is a standard Matlab system function in the format required by
c	SIMULINK and the associated analysis and integration routines.

c	written in:  Language Systems Fortran Version 3.0.1
c	Should work with most Fortran compilers.

c	NOTE:	 	For version 3.5 of Matlab.  WILL NOT work with version 4.x 
c			or later.  If you are running a later version, see
c			te_4.fv4

c	Typical MATLAB calls are:

c		dxdt=te_4(t,x,u,1,p);
c		   y=te_4(t,x,u,3,p);

c	where "t" is the current time, dxdt is a calculated vector of time 
c	derivatives of the states, and the other variables are as described
c	below:

c	STATE VARIABLES

c	x(1)	Molar holdup of A in the reactor (kmol)
c	 (2)	Molar holdup of B in the reactor (kmol)
c	 (3)	Molar holdup of C in the reactor (kmol)
c	 (4)	Molar holdup of D in the reactor (kmol)
c	 (5)  Current Feed 1 valve position (%)
c	 (6) 	Current Feed 2 valve position (%)
c	 (7)	Current Purge valve position (%)
c	 (8)	Current Product valve position (%)

c	INPUT VARIABLES

c	u(1)	Desired Feed 1 valve position (%)
c	 (2)	Desired Feed 2 valve position (%)
c	 (3)	Desired Purge valve position (%)
c	 (4)	Setpoint for liquid level controller (%)

c	OUTPUT VARIABLES

c	y(1)	Feed 1 flow measurement (kmol/h)
c	 (2)	Feed 2 flow measurement (kmol/h)
c	 (3)	Purge flow measurement (kmol/h)
c	 (4)	Product flow measurement (kmol/h)
c	 (5)	Pressure (kPa)
c	 (6)	Liquid level (%)
c	 (7)	Mol % A in purge
c	 (8)	Mol % B in purge
c	 (9)	Mol % C in purge
c	(10)	Instantaneous operating cost ($/kmol)

c	PARAMETERS

c	p(1)	Sampling delay in gas composition measurement (h)
c	 (2)	Mol fraction A in Feed 1.
c	 (3)	Mol fraction B in Feed 1.
c	 (4)	Maximum allowed position for product valve (%)
c	 (5)	Gain for the level controller (%/%)
c	 (6)	Nominal steady-state for product valve (%)
c	 (7)	NOT USED -- user should set == 0
c	 (8)	Kinetic parameter -- pre-exponential k0 of eq. 5.
c	 (9)	Kinetic parameter -- power on Pc in eq. 5.
c	(10)-(50)  User should initialize to zero, then
c		leave unchanged until end of simulation.  These
c		variables are used to store intermediate results.

c	NOTE:	See provided Matlab scripts for example uses of this code.

C---------------------------------------------------------------------

      SUBROUTINE USRFCN(NLHS, PLHS, NRHS, PRHS)
      INTEGER PLHS(*), PRHS(*)
      INTEGER NLHS, NRHS
C
      INTEGER CRTMAT, REALP, IMAGP, GETGLO, ALREAL, ALINT, GETSTR
C
      DOUBLE PRECISION GETSCA
C
C THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS ARE STANDARD
C IN ALL FORTRAN MEX FILES.
C---------------------------------------------------------------------
C

c	The following will hold pointers:

	integer t,x,u,flag,p1,p2

c		The following are local variables.  In this case we know
c		that there will be 11 continuous states, 5 inputs, and 17 outputs.
c		These values are specified by the parameters nx, nu, and ny,
c		respectively.  Also define corresponding character constants
c		lenx, lenu, to use in error messages, and double precision
c		constants for the SIZES vector.

	integer iflag,m,n,nx,nu,ny
	character*(*) lenu,lenx
	double precision dnx, dnu, dny
	parameter (nx=8, nu=4, ny=10, lenx='8', lenu='4', dnx=8.,
     +		dnu=4.,dny=10.)
	
	
	double precision 	rx(nx),rx0(nx),dxdt(nx),sizes(6),rt,ry(ny),ru(nu)


	character*80 message	! For debug messages to screen
	
c		The following is the appropriate definition of the SIZES
c		vector for SIMULAB system functions when there are no
c		discrete-time states and no direct feed-through from input
c		to output.

	data sizes/dnx, 0., dny, dnu, 0., 0./
	
c		Define default initial state vector.  The first 4 are the molar 
c		holdups of the 4 components.  The next 4 are the
c		valve positions.

	data rx0/	4.449999958429348e+01,1.353296996509594e+01,
     +		3.664788062995841e+01,1.100000000000000e+02,
     +		6.095327313484253e+01,2.502232231706676e+01,
     +		3.925777017606444e+01,4.703024823457651e+01/

	INCLUDE 'te4com.inc'
	
c		Internal storage

	double precision 	F(4)
     
     	double precision RR1,Ntot,NL,NV,VL,VV,Pa,Pc,PT,y3(3),u4,VLpct
	
	integer iter,iratio


c		Check the number of input arguments

	if (nrhs .ne. 5) then
		call mexerr('TE_4 needs 5 input variables.')
	else

c			There are 5 input arguments.  Set pointers and copy the
c			inputs to local storage.

		call getsiz(prhs(1),m,n)		! Get size of "t" array.
		if (m.gt.0 .and. n.gt.0) then		! If it's non-empty...
			rt=getsca(prhs(1))		! gets time value (a scalar).
		else
			rt=0.0				! Otherwise set to zero.
		end if
		
		call getsiz(prhs(4),m,n)		! "flag"
		if (m.gt.0 .and. n.gt.0) then
			if (m.eq.1 .and. n.eq.1) then
				iflag=abs(nint(getsca(prhs(4))))
			else
				call mexerr('FLAG must be a scalar.')
			end if
		else
			iflag=0
		end if
			
		call getsiz(prhs(2),m,n)		! "x" vector.
		if (m.gt.0 .and. n.gt.0) then		
			if (min(m,n).ne.1 .or. max(m,n) .ne. nx) then	
				call mexerr('X must be a vector, length '//lenx)
			end if
			x=realp(prhs(2))
			call cpout(x,rx,nx)		! copy x into local storage
		else if (iflag.eq.1 .or. iflag.eq.3) then
			call mexerr('X must be supplied.')
		end if
	
			
		call getsiz(prhs(3),m,n)		! "u" vector.
		if (m.gt.0 .and. n.gt.0) then		
			if (min(m,n).ne.1 .or. max(m,n) .ne. nu) then	
				call mexerr('U must be a vector, length '//lenu)
			end if
			u=realp(prhs(3))
			call cpout(u,ru,nu)		! copy u into local storage
		else if (iflag.eq.1 .or. iflag.eq.3) then
			call mexerr('U must be supplied.')
		end if
		
		call getsiz(prhs(5),m,n)		! p1
		if (m.eq.1 .and. n.eq.np1 .or. m.eq.np1 .and. n.eq.1) then
			p1=realp(prhs(5))
			call cpout(p1,dp1,np1)		! copy p1 to local REAL*8 array 
 		else
			call mexerr('P1 must be a vector, length 50')
		end if
	
	end if
	
	if (rt.lt.0.0.and.iflag.eq.0) then

c		Initialize certain parameters that change from one call to the next.

		tpurge=tdgas-1.0d-4				! time of next purge sample
				
	end if
	
	
c		Now return outputs according to the value of "flag"

	if (iflag.eq.0) then
	
		plhs(1)=crtmat(6,1,0)				! "sizes"
		call cpin(sizes,realp(plhs(1)),6)
		if (nlhs.gt.1) then
			plhs(2)=crtmat(nx,1,0)			! "x0"
			call cpin(rx0,realp(plhs(2)),nx)
		end if
		if (nlhs.gt.2) then
			plhs(3)=crtmat(0,0,0)			! "xstr" (empty)
		end if
		
	else if (iflag.eq.1 .or. iflag.eq.3) then
	
c		First calculate the output variables (since they are needed in
c		either case).

		NL=rx(4)			! Total liquid moles [kmol]
		VL=NL/Lden			! Liquid volume [m^3]
		VV=VT-VL			! Vapor volume [m^3]
		yc1=1.0d0-ya1-yb1		! Mole fraction C in feed 1
		NV=rx(1)+rx(2)+rx(3)	! Total vapor moles [kmol]
		PT=NV*Rgas*Tgas/VV	! Total pressure [kPa]
		VLpct=VL*100.0d0/30.0d0	! Liquid volume as a percentage of capacity
		
c		Flowrates
		
		F(1)=Cv(1)*rx(5)
		F(2)=Cv(2)*rx(6)
		F(3)=Cv(3)*rx(7)*sqrt(PT-100.0d0)
		F(4)=Cv(4)*rx(8)*sqrt(PT-100.0d0)
		
		if (rx(4) .le. 0.0) then
			F(4)=0.0d0
		end if

c		Vapor composition

		y3(1)=rx(1)/NV		! ya3
		y3(2)=rx(2)/NV		! yb3
		y3(3)=rx(3)/NV		! yc3
		
c		Special check for t=0 to set up delayed sampling

		if (rt .le. 0.0) then
			do i=1,3
				ylast(i)=y3(i)
				ymeas(i)=ylast(i)
			end do
		
c		Check whether purge gas should be "sampled"
		
		else if (rt.gt.tpurge .and. iflag.eq.3) then
			tpurge=tpurge+tdgas
			do i=1,3
				ymeas(i)=ylast(i)
				ylast(i)=y3(i)
			end do
		end if
		
c		Reaction rate.  Two parameters can be user specified (the constant
c		of proportionality and the exponent on the C partial pressure.
c		If these are <= 0, defaults are used.
		
		if (kpar .le. 0.0d0) then
			kpar=0.00117
		end if
		if (nCpar .le. 0.0d0) then
			nCpar=0.4
		end if
		Pa=y3(1)*PT
		Pc=y3(3)*PT
		RR1=kpar*(Pa**1.2)*(Pc**nCpar)
		
     
c		If flag=1, calculate dx/dt

		if (iflag.eq.1) then

c			Composition dynamics

			dxdt(1)=ya1*F(1)+F(2)-y3(1)*F(3)-RR1
			dxdt(2)=yb1*F(1)-y3(2)*F(3)
			dxdt(3)=yc1*F(1)-y3(3)*F(3)-RR1
			dxdt(4)=RR1-F(4)
		
c			Valve dynamics.

c			Enforce saturation limits on manipulated variables.  Normally
c			these are 0 and 100%, but input 2 is a special case, for which the
c			upper limit is user-specified (in the p1 parameter list).  Note that
c			exponential functions are used to prevent discontinuities.

			
c			The u(4) value is the setpoint for the liquid level, going to
c			an analog P controller.

			u4=u4bar+KcVL*(ru(4)-VLpct)
			ru(4)=u4
			
			do i=1,4
				if (ru(i) .le. 1.0d0) then
					ru(i)=exp(ru(i)-1.0d0)
				else if (ru(i) .ge. 99.0d0) then
					ru(i)=100.0d0-exp(99.0d0-ru(i))
				end if
				if (i.eq.2) then
					ru(i)=min(ru(i),u2max)
				end if
				dxdt(i+4)=(ru(i)-rx(i+4))/tauvlv
			end do
						
			plhs(1)=crtmat(nx,1,0)
			call cpin(dxdt,realp(plhs(1)),nx)  		! returns dxdt to Matlab
		
		else
		
c		If flag=3, calculate outputs

			ry(1)=F(1)
			ry(2)=F(2)
			ry(3)=F(3)
			ry(4)=F(4)
			ry(5)=PT
			ry(6)=VLpct
			ry(7)=ymeas(1)*100.0d0
			ry(8)=ymeas(2)*100.0d0
			ry(9)=ymeas(3)*100.0d0
			ry(10)=F(3)*(ymeas(1)*2.206+ymeas(3)*6.177)/F(4)  ! cost/mole product
			
			plhs(1)=crtmat(ny,1,0)
			call cpin(ry,realp(plhs(1)),ny)  	! returns y to Matlab
		end if		
	else

c			For all other flags, return an empty matrix.	

		plhs(1)=crtmat(0,0,0)
	
	end if
	
c		Save certain variables for next call.

	call cpin(dp1,p1,np1)
	
	return
	end
