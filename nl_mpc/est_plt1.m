function est_plt1(dt,y,yest,d_s)
 
[M,N]=size(y);
tplt=dt*[0:M-1]';

clg
subplot(211)
plot(tplt,yest(1:M,1),'-',tplt,y(1:M,7),'--')
title('React Pressure')
subplot(212)
plot(tplt,yest(1:M,2),'-',tplt,y(1:M,8),'--')
title('Reactor level')
pause

clg
subplot(221)
plot(tplt,yest(1:M,3),'-',tplt,y(1:M,13),'--')
title('Sep Pressure')
subplot(223)
plot(tplt,d_s(1:M,6))
title('R/S Flow factor')
subplot(222)
plot(tplt,yest(1:M,4),'-',tplt,y(1:M,12),'--')
title('Separator level')
subplot(224)
plot(tplt,yest(1:M,5),'-',tplt,y(1:M,15),'--')
title('Stripper level')
pause

clg
subplot(221)
plot(tplt,yest(1:M,22),'-',tplt,y(1:M,40),'--')
title('Product (G %)')
subplot(222)
plot(tplt,yest(1:M,23),'-',tplt,y(1:M,41),'--')
title('Product (H %)')
subplot(223)
plot(tplt,yest(1:M,14),'-',tplt,y(1:M,29),'--')
title('A in Purge')
subplot(224)
plot(tplt,yest(1:M,15),'-',tplt,y(1:M,30),'--')
title('B in Purge')
pause

clg
subplot(221)
plot(tplt,yest(1:M,16),'-',tplt,y(1:M,31),'--')
title('C in Purge')
subplot(222)
plot(tplt,yest(1:M,17),'-',tplt,y(1:M,32),'--')
title('D in Purge')
subplot(223)
plot(tplt,yest(1:M,18),'-',tplt,y(1:M,33),'--')
title('E in Purge')
subplot(224)
plot(tplt,yest(1:M,19),'-',tplt,y(1:M,34),'--')
title('F in Purge')
pause

clg
subplot(221)
plot(tplt,yest(1:M,20),'-',tplt,y(1:M,35),'--')
title('G in Purge')
subplot(223)
plot(tplt,yest(1:M,21),'-',tplt,y(1:M,36),'--')
title('H in Purge')
subplot(222)
plot(tplt,yest(1:M,8),'-',tplt,y(1:M,23),'--')
title('A in Feed')
subplot(224)
plot(tplt,d_s(1:M,1))
title('A in Stream 4')
pause

clg
subplot(221)
plot(tplt,yest(1:M,9),'-',tplt,y(1:M,24),'--')
title('B in Feed')
subplot(222)
plot(tplt,yest(1:M,10),'-',tplt,y(1:M,25),'--')
title('C in Feed')
subplot(223)
plot(tplt,d_s(1:M,2))
title('B in Stream 4')
pause

clg
subplot(221)
plot(tplt,yest(1:M,11),'-',tplt,y(1:M,26),'--')
title('D in Feed')
subplot(222)
plot(tplt,yest(1:M,12),'-',tplt,y(1:M,27),'--')
title('E in Feed')
subplot(223)
plot(tplt,yest(1:M,13),'-',tplt,y(1:M,28),'--')
title('F in Feed')
pause

clg
subplot(221)
plot(tplt,yest(1:M,6),'-',tplt,y(1:M,16),'--')
title('Stripper/Feed zone Pressure [kPa]')
subplot(222)
plot(tplt,yest(1:M,7),'-',tplt,y(1:M,6),'--')
title('Feed to reactor [kscmh]')
subplot(223)
plot(tplt,d_s(1:M,7))
title('F/R flow factor')
subplot(224)
plot(tplt,yest(1:M,25),'-',tplt,y(1:M,42),'--')
title('Operating cost')
pause

clg
subplot(221)
plot(tplt,yest(1:M,26),'-',tplt,y(1:M,43),'--')
title('Reaction 1 [kmol/h]')
subplot(222)
plot(tplt,yest(1:M,27),'-',tplt,y(1:M,44),'--')
title('Reaction 2 [kmol/h]')
subplot(223)
plot(tplt,d_s(1:M,3))
title('Reaction 1 factor')
subplot(224)
plot(tplt,d_s(1:M,4))
title('Reaction 2 factor')
pause


clg
subplot(221)
plot(tplt,yest(1:M,28),'-',tplt,y(1:M,45),'--')
title('Reaction 3 [kmol/h]')
subplot(223)
plot(tplt,d_s(1:M,5))
title('F10 bias [kmol/h]')
subplot(222)
plot(tplt,100-yest(1:M,22)-yest(1:M,23),'-',...
        tplt,100-y(1:M,40)-y(1:M,41),'--')
title('Product impurities')
subplot(224)
plot(tplt,d_s(1:M,8))
title('Impurity factor')
pause

clg
subplot(221)
plot(tplt,yest(1:M,29),'-',tplt,y(1:M,46),'--')
title('PA in Reactor')
subplot(222)
plot(tplt,yest(1:M,30),'-',tplt,y(1:M,47),'--')
title('PC in Reactor')
subplot(223)
plot(tplt,yest(1:M,31),'-',tplt,y(1:M,48),'--')
title('PD in Reactor')
subplot(224)
plot(tplt,yest(1:M,32),'-',tplt,y(1:M,49),'--')
title('PE in Reactor')
pause

clg
subplot(221)
plot(tplt,d_s(1:M,9))
title('E VLE in separator')
subplot(222)
plot(tplt,d_s(1:M,10))
title('G VLE in separator')
subplot(223)
plot(tplt,d_s(1:M,15))
title('H VLE in reactor')
pause

clg
subplot(221)
plot(tplt,d_s(1:M,11))
title('C feed bias [kmol/h]')
subplot(222)
plot(tplt,d_s(1:M,12))
title('D feed bias [kmol/h]')
subplot(223)
plot(tplt,d_s(1:M,13))
title('E feed bias [kmol/h]')
subplot(224)
plot(tplt,d_s(1:M,14))
title('F feed bias [kmol/h]')
