	integer ncomp,np1

	parameter (ncomp=4,np1=50)

	double precision 	tpurge,tdgas,VT,VLmax,Rgas,Tgas,Lden,
     +			dp1(np1), ylast(ncomp),ymeas(ncomp),pmeas(4),
     +			tauvlv,ya1,yb1,yc1,Cv(4),u2max,KcVL,u4bar,
     +			kpar,nCpar
     
     	data cv/3.3046, 0.2246, 0.00352, 0.0417/
     	
	parameter 	(tauvlv=2.77d-3, VT=122.0d0, VLmax=30.0d0)
	parameter	(Rgas=8.314d0, Tgas=373.0d0, Lden=8.3d0)

     	equivalence (dp1(1),tdgas),(dp1(2),ya1),(dp1(3),yb1)
	equivalence (dp1(4),u2max),(dp1(5),KcVL),(dp1(6),u4bar)
	equivalence (dp1(8),kpar),(dp1(9),nCpar)
	equivalence	(dp1(10), tpurge)
	equivalence (dp1(15),ylast(1)),(dp1(25),ymeas(1))

	common /TE_pars/dp1
