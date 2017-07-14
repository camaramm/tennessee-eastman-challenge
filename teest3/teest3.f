c	Simplified TE model in .MEX format
c	This one is designed for use as state estimator
c	The Matlab MEX interface is for version 3.5 of Matlab

c	Copyright N. L. Ricker, November, 1993.

c	STATE VARIABLES are all molar holdups [kmol] in a certain
c	location.

c	 1	A in reactor
c	 2	B in reactor
c	 3	C in reactor
c	 4	D in reactor
c	 5	E in reactor
c	 6	F in reactor
c	 7	G in reactor
c	 8	H in reactor
c	 9	A in separator
c	10	B in separator
c	11	C in separator
c	12	D in separator
c	13	E in separator
c	14	F in separator
c	15	G in separator
c	16	H in separator
c	17	A in feed zone
c	18	B in feed zone
c	19	C in feed zone
c	20	D in feed zone
c	21	E in feed zone
c	22	F in feed zone
c	23	G in feed zone
c	24	H in feed zone
c	25	G in product reservoir (stripper bottoms)
c	26	H in product reservoir

c	INPUTS are flows [kmol/h] unless noted otherwise:

c	 1	Feed 1 (pure A)
c	 2	Feed 2 (pure D)
c	 3	Feed 3 (pure E)
c	 4	Feed 4 (A & C)
c	 5	Recycle (stream 8)
c	 6	Purge (stream 9)
c	 7	Separator underflow (stream 10)
c	 8	Product rate (stream 11)
c	 9	Reactor temperature [deg C]
c	10	Separator temperature [deg C]
c	11	A in stream 4 [Mole %]
c	12	B in stream 4 [Mole %]
c	13	Reaction 1 activity factor [%].
c	14	Reaction 2 activity factor [%].
c	15	Bias subtracted from u(7) to get stream 10 flow [kmol/h].
c	16	Reactor/Sep flow parameter [%].
c	17	Feed/Reactor flow parameter [%].
c	18	Product G+H purity parameter [%].
c	19	Adjustment to VLE of D to G in separator [%]
c	20	Adjustment to H VLE in separator [%]
c	21	C bias flow to feed zone [kmol/h]
c	22	D bias flow to feed zone [kmol/h]
c	23	E bias flow to feed zone [kmol/h]
c	24	F bias flow to feed zone [kmol/h]
c	25	VLE correction to D-H in reactor [%]

c	OUTPUTS are mole % in stream unless noted otherwise:

c	 1	Reactor pressure [kPa]
c	 2	Reactor liq. holdup [%]
c	 3	Separator pressure [kPa]
c	 4	Separator liq. holdup [%]
c	 5	Product liq. holdup [%]
c	 6	Feed zone pressure [kPa]
c	 7	Total feed entering reactor (stream 6) [kscmh]
c	 8	A in reactor feed (stream 6)
c	 9	B in reactor feed (stream 6)
c	10	C in reactor feed (stream 6)
c	11	D in reactor feed (stream 6)
c	12	E in reactor feed (stream 6)
c	13	F in reactor feed (stream 6)
c	14	A in purge (stream 9)
c	15	B in purge (stream 9)
c	16	C in purge (stream 9)
c	17	D in purge (stream 9)
c	18	E in purge (stream 9)
c	19	F in purge (stream 9)
c	20	G in purge (stream 9)
c	21	H in purge (stream 9)
c	22	G in product (stream 11)
c	23	H in product (stream 11)
c	24	Molar production rate (same as input 8) [kmol/h].
c	25	Cost [cents/kmol product]
c	26	Rate of reaction 1 [kmol G produced/h]
c	27	Rate of reaction 2 [kmol H produced/h]
c	28	Rate of reaction 3 [kmol F produced/h]
c	29	Partial pressure of A in reactor [kPa]
c	30	Partial pressure of C in reactor [kPa]
c	31	Partial pressure of D in reactor [kPa]
c	32	Partial pressure of E in reactor [kPa]
c
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

	integer t,x,u,flag,p1,p2,p3,p4

c		The number of states, inputs, and outputs are specified
c		by the parameters nx, nu, and ny, respectively.  
c		Also define corresponding character constants
c		lenx, lenu, to use in error messages, and double precision
c		constants for the SIZES vector.

	integer iflag,m,n,nx,nu,ny
	character*(*) lenu,lenx
	double precision dnx, dnu, dny
	parameter (nx=26, nu=25, ny=32, lenx='26', lenu='25', dnx=26.,
     +		dnu=25.,dny=32.)
	
	
	double precision 	rx(nx),rx0(nx),dxdt(nx),sizes(6),rt,ru(nu),
     +			ry(ny)

	character*80 message	! For debug messages to screen
	
c		The following is the appropriate definition of the SIZES
c		vector for SIMULAB system functions when there are no
c		discrete-time states and no direct feed-through from input
c		to output.

	data sizes/dnx, 0., dny, dnu, 0., 0./
	
	data rx0/26*0.d0/

c		Check the number of input arguments

	if (nrhs .gt. 0 .and. nrhs .ne. 4) then
		call mexerr('TEest3 needs 4 input variables.')
	else if (nrhs .le. 0) then
		iflag=0
	else

c			There are 4 input arguments.  Set pointers and copy the
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
		
	end if	
	
c		Now return outputs according to the value of "flag"

	if (iflag.eq.0) then
	
c			Return default initial states and size info
	
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
	
		call TEsimf(nx,nu,ny,rx,ru,ry,dxdt)    	! Main simulation call

		if (iflag.eq.1) then
			plhs(1)=crtmat(nx,1,0)
			call cpin(dxdt,realp(plhs(1)),nx)  	! returns dxdt to Matlab
		else
			plhs(1)=crtmat(ny,1,0)
			call cpin(ry,realp(plhs(1)),ny)  	! returns y to Matlab
		end if		
	else

c			For all other flags, return an empty matrix.	

		plhs(1)=crtmat(0,0,0)
	
	end if
	
	return
	end
	
c	*******************************************************************************

	SUBROUTINE TEsimf(nx,nu,ny,x,u,y,dxdt)
	
c	Calculates outputs and rates of change of states.

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)

	integer nx,nu,ny,nc
	double precision cf2cm,dzero,pctGH
	
	parameter (nc=8,cf2cm=0.028317,TauTr=0.08333,TauTs=0.08333,
     +		TauF=0.005,TauC=0.08333,dzero=0.0d0)
	
	double precision x(nx),u(nu),y(ny),dxdt(nx)
	
	double precision 	rlmol(nc),slmol(nc),rvmol(nc),svmol(nc),
     +			fvmol(nc),xlr(nc),xls(nc),RR(4),
     +			ppr(nc),pps(nc),xvr(nc),xvs(nc),xvf(nc),
     +			Fcmol(nc)

	double precision 	mwts(nc),molvol(nc),avp(nc),bvp(nc),
     +			cvp(nc),sfr(nc),gami(nc),gamr(nc)
	integer code
	
	data AVP/3*0.0,20.81,2*21.24,21.32,22.10/			! A in Antoine eqn.
      data BVP/3*0.0,-1444.0,2*-2114.0,-2748.0,-3318.0/	! B in Antoine eqn.
	data CVP/3*0.0,259.0,2*266.0,233.0,250.0/			! C in Antoine eqn.
	data mwts/2.0,25.4,28.0,32.0,46.0,48.0,62.0,76.0/	! molec. wts.
	data molvol/3*0.0,0.1070,0.1260,0.1463,0.1013,0.1231/	! [m3/kmol]
	data sfr/6*1.0,0.07,0.04/					! stripping factors
	data gami/8*1.0/							! sep. activ. coeffs.
	data gamr/8*1.0/							! react. activ. coeffs.



c		Get disturbances from input vector. Flow disturbances are scaled so
c		they're comparable to manipulated variables (0-100% basis)

	xA4=u(11)/100.			! mol frac B in stream 4 (A+C feed)
	xB4=u(12)/100.			! mol frac C in stream 4 (A+C feed)
	xC4=1.0-xA4-xB4			! get A in stream 4 by difference
	R1F=u(13)/100.			! Reactivity factor, reaction 1
	R2F=u(14)/100.			! Reactivity factor, reaction 2
	F10b=u(15)				! Bias adjustment for stream 10 flow [kmol/h]
	RSflow=u(16)/100.			! Adjusts flow between reactor and separator
	FRflow=u(17)/100.			! Adjusts flow between feed zone and reactor
	fracGH=u(18)/100.			! Adjusts for impurity in product (components other
						! than G and H)
	do i=4,7
		gami(i)=u(19)/100.	! Adjusts VLE of D to G in separator
	end do
	gami(8)=u(20)/100.		! Adjusts VLE of H in separator
	Cbias=u(21)				! C flow bias at feed point
	Dbias=u(22)				! D flow bias at feed point
	Ebias=u(23)				! E flow bias at feed point
	Fbias=u(24)				! F flow bias at feed point
	do i=4,8
		gamr(i)=u(25)/100.	! Adjusts VLE in reactor
	end do
	
c		Temperatures in reactor & separator

	TCR=u(9)				! reactor temp [C]
	TKR=TCR+273.2			! reactor temp [K]
	TCS=u(10)				! separator temp [C]
	TKS=TCS+273.2			! separator temp [K]
	
c		Product reservoir states (stripper bottoms)

	plG=dmax1(x(25),dzero)			! G holdup [kmol]
	plH=dmax1(x(26),dzero)			! H holdup [kmol]
	VLP=plG*molvol(7)+plH*molvol(8)	! Volume of liquid [m3]
	xGp=fracGH*(plG/(plG+plH))		! mol frac G in product.
							! adjusted for assumed impurity level (fracGH)
	xHp=fracGH*(plH/(plG+plH))		! mol frac H in product

c		Get molar holdups in reactor, separator, and feed vapor.  
c		Assumes A, B, C insoluble in liquid.

	rlsum=0.0
	slsum=0.0
	rvsum=0.0
	svsum=0.0
	VLR=0.0
	VLS=0.0
	
	do i=1,3
		rlmol(i)=0.0
		slmol(i)=0.0
	end do
	
	do i=4,nc	
		rlmol(i)=dmax1(x(i),dzero)	! gets liq. mol from state vector
		slmol(i)=dmax1(x(i+8),dzero)
		rlsum=rlsum+rlmol(i)		! accumulate total moles
		slsum=slsum+slmol(i)
		VLR=VLR+rlmol(i)*molvol(i)	! accumulate liquid volume [m3]
		VLS=VLS+slmol(i)*molvol(i)
	end do
	Pvol=VLS/slsum				! molar volume of product [m3/kmol]
	
	fvsum=0.0
	do i=1,nc
		fvmol(i)=dmax1(x(i+16),dzero)	! gets feed vapor from state vector
		fvsum=fvsum+fvmol(i)		! total feed vapor moles
	end do
		
c		Volumes in reactor, separator, & feed zone

	VTR=36.8				! total reactor volume [m3]
      VTS=99.1				! total separator volume [m3]
      VTV=150.0				! total feed vapor holdup volume [m3]
	VVR=VTR-VLR				! vapor in reactor [m3]
	VVS=VTS-VLS				! vapor in separator [m3]
	
c		Mole fractions, partial pressures, etc.

	PTR=0.0
	PTS=0.0
	Rg=8.314					! gas constant [kPa*m3/kmol/K]
	PTV=fvsum*Rg*(273.2+86.1)/VTV		! Total pressure in feed zone [kPa]
							! (assumes constant temperature of 86.1 C)

	do i=1,3							! for A, B, C
		ppr(i)=dmax1(x(i),dzero)*Rg*TKR/VVR		! pp in reactor [kPa]
		pps(i)=dmax1(x(i+8),dzero)*Rg*TKS/VVS	! pp in separator [kPa]
		PTR=PTR+ppr(i)					! total P in reactor [kPa]
		PTS=PTS+pps(i)					! total P in separator [kPa]
		xlr(i)=0.0
		xls(i)=0.0
	end do

	do i=4,nc								! For D,E,F,G,H
		xlr(i)=rlmol(i)/rlsum					! reactor liquid mol frac
		pvap=1.e-3*dexp(avp(i)+bvp(i)/(cvp(i)+tcr))	! vap. pres. [kPa]
		ppr(i)=pvap*xlr(i)*gamr(i)				! VLE
		xls(i)=slmol(i)/slsum					! Repeat for separator
		pvap=1.e-3*dexp(avp(i)+bvp(i)/(cvp(i)+tcs))	
		pps(i)=pvap*xls(i)*gami(i)				! Note use of VLE adjustment
		PTR=PTR+ppr(i)						! Reactor total press [kPa]
		PTS=PTS+pps(i)						! Separator total press [kPa]
	end do
	
	wtmolr=0.0
	wtmolf=0.0
	do i=1,8
		xvr(i)=ppr(i)/PTR			! reactor vapor mol frac
		xvs(i)=pps(i)/PTS			! sep vapor mol frac
		xvf(i)=fvmol(i)/fvsum		! feed vapor mol frac
		wtmolr=wtmolr+xvr(i)*mwts(i)	! average mol wt of reactor vapor
		wtmolf=wtmolf+xvf(i)*mwts(i)	! average mol wt of feed vapor
	end do

c		Rate Laws.  All rates in [kmol/h].

      RR(1)=VVR*DEXP(44.06-42600.0/1.987/TKR)*R1F*(PPR(1)**1.08)
     +		*(PPR(3)**0.311)*(PPR(4)**0.874)
      RR(2)=VVR*DEXP(10.27-19500.0/1.987/TKR)*R2F*(PPR(1)**1.15)
     +		*(PPR(3)**0.370)*(PPR(5)**1.00)
      RR(3)=VVR*DEXP(59.50-59500.0/1.987/TKR)*PPR(1)
     +		*(0.77*PPR(4)+PPR(5))
	

c		Define total flows [kmol/h]

	F1=u(1)
	F2=u(2)
	F3=u(3)
	F4=u(4)
	F8=u(5)
	F9=u(6)
	F10=u(7) - F10b
	F11=u(8)

c		Stripper balance to get combined feeds (streams 1+2+3+5)
c		Also account for flow biases.

	Fcmol(1)=F1+xA4*F4
	Fcmol(2)=xB4*F4
	Fcmol(3)=xC4*F4+Cbias
	Fcmol(4)=F2+sfr(4)*F10*xls(4)+Dbias
	Fcmol(5)=F3+sfr(5)*F10*xls(5)+Ebias
	Fcmol(6)=sfr(6)*F10*xls(6)+Fbias
	Fcmol(7)=sfr(7)*F10*xls(7)
	Fcmol(8)=sfr(8)*F10*xls(8)
	Fcomb=0.0
	do i=1,8
		Fcomb=Fcomb+Fcmol(i)
	end do
		
c		Pressure drop vs. flow to get flows of streams 6 & 7.

	F6=(2413.7/wtmolf)*sqrt(abs(PTV-PTR))*FRflow
	if (ptr .gt. ptv) then
		F6=-F6
	end if
	F7=(5722.0/wtmolr)*sqrt(abs(PTR-PTS))*RSflow
	if (pts .gt. ptr) then
		F7=-F7
	end if

c		Balances to get rates of change of states.  First do  
c		convection terms (in - out).

	do i=1,8
		dxdt(i)=F6*xvf(i) - F7*xvr(i)					! Reactor
		dxdt(i+8)=F7*xvr(i) - (F8+F9)*xvs(i)
     +		 		- F10*xls(i)				! Separator
     		dxdt(i+16)=Fcmol(i)+F8*xvs(i)-F6*xvf(i)			! Feed zone
	end do
	
	dxdt(25)=(1.0-sfr(7))*F10*xls(7) - F11*xGp		! G in product
	dxdt(26)=(1.0-sfr(8))*F10*xls(8) - F11*xHp		! H in product

c		Now account for reaction.

	dxdt(1)=dxdt(1) - rr(1) - rr(2) - 0.333*rr(3)
	dxdt(2)=dxdt(2)
	dxdt(3)=dxdt(3) - rr(1) - rr(2)
	dxdt(4)=dxdt(4) - rr(1) - rr(3)
	dxdt(5)=dxdt(5) - rr(2) - 0.333*rr(3)
	dxdt(6)=dxdt(6) + rr(3)
	dxdt(7)=dxdt(7) + rr(1)
	dxdt(8)=dxdt(8) + rr(2)
	
c		Calculate the outputs.

	y(1)=PTR-101.					! Reactor pressure 	[kPa guage]
	y(2)=5.263*VLR - 12.105				! Reactor liq. 		[%]
	y(3)=PTS-101.					! Separator pressure	[kPa guage]
	y(4)=12.28*VLS - 10.53				! Separator liquid	[%]
	y(5)=22.58*VLP - 49.03				! Product liq. holdup 	[%]
	y(6)=PTV-101.					! Pressure in feed zone [kPa guage]
	y(7)=F6/44.79					! Stream 6 flowrate 	[kscmh]
	do i=1,6
		y(i+7)=xvf(i)*100.0			! Reactor feed 		[mol %]
		y(i+13)=xvs(i)*100.0			! Purge 			[mol %]
	end do

	y(20)=xvs(7)*100.0				! G in purge		[mol %]
	y(21)=xvs(8)*100.0				! H in purge		[mol %]
	y(22)=xGp*100.0					! G in product		[mol %]
	y(23)=xHp*100.0					! H in product		[mol %]
	y(24)=F11						! Molar production    	[kmol/h]
	y(25)=(F9/(F11+0.01))*(221.*xvs(1) 
     +				+ 618.*xvs(3) + 2210.*xvs(4)
     +				+ 1460.*xvs(5) + 1790.*xvs(6) 
     +				+ 3040.*xvs(7) + 2290.*xvs(8))          ! [cents/kmol]

	y(26)=RR(1)						! G production		[kmol/h]
	y(27)=RR(2)						! H production		[kmol/h]
	y(28)=RR(3)						! F production		[kmol/h]
	y(29)=PPR(1)					! A PP in reactor		[kPa]
	y(30)=PPR(3)					! C PP in reactor		[kPa]
	y(31)=PPR(4)					! D PP in reactor		[kPa]
	y(32)=PPR(5)					! F PP in reactor		[kPa]

	return
	end
	
