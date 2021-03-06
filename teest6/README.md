teest6
=====

There are at least 2 ways you can use the TEest6 model.

1. As a MEX file in Matlab.

2. Called directly in Fortran.

Both methods are illustrated here. The following is a Matlab script.  See comments in the script to understand what it does. You will need the following special Matlab functions:

- `TEest6.mex`: The MEX version of the `TEest6` model. To generate this, you must compile and link TEest6.f in the proper way, according to the instructions for MEX files in your Matlab manual. Note that `TEest6.f` is written for Matlab version 4.x, and will not work for earlier versions. In that case, try `TEest6.f.V35`.
    
- `TEest6_ss.m`: available in this directory. Note that the file name is longer than 8 characters, which may cause trouble on some machines ... sorry about that. You can get around that by changing the file name *AND* the function statement name (first line in the file) *AND* the name used in the script below.
    
```matlab
% START OF TEST SCRIPT

% Set y0 equal to the steady-state outputs from the Downs
% and Vogel base case.

y0=[0.2505213768
 3664.0388848630
 4509.2683901682
    9.3477423425
   26.9024403996
   42.3388488563
 2705.0000000464
   74.9999999809
  120.4000000017
    0.3371170513
   80.1094022593
   49.9999999957
 2633.7278831266
   25.1601733849
   50.0000000045
 3102.2131372370
   22.9487300898
   65.7310290341
  230.3078204048
  341.4336508454
   94.5992754900
   77.2969835300
   32.1876130945
    8.8932726456
   26.3832894332
    6.8820117573
   18.7763059474
    1.6567273566
   32.9575372908
   13.8229458674
   23.9781442398
    1.2565424847
   18.5792697025
    2.2633413851
    4.8436409545
    2.2985780751
    0.0178655250
    0.8357022623
    0.0985769135
   53.7239787893
   43.8280218936
   80.7495835086
  114.2295345460
   92.9394845628
    0.5152329090
  762.3066026486
  554.6135778697
   30.1595684440
  497.3423475977
   53.7239787893
   43.8280218936];

% iymeas is an index vector that gives the correspondence
% between the outputs of TEest6 and TEPROB (the Downs and
% Vogel model).  For example, output 1 of TEest6 is equivalent
% to output 7 of TEPROB ... the reactor pressure.

iymeas=[7,8,13,12,15,16,6,23:41,17,5,20,18,11,22];

% The goal is to determine vectors of states and
% parameters for TEest6 such that we have a steady state, 
% and match as many of the outputs of TEPROB as possible.

% yest will contain the DESIRED outputs from TEest6.
% The first 32 elements have equivalents in the TEPROB model --
% see definition of iymeas, above.

yest=zeros(40,1); 
yest(1:32)=y0(iymeas);  % Sets the first 32 outputs of TEest6
                        % equal to equivalent outputs of TEPROB

% The measured product concentration from TEest6 is normalized
% to sum to unity.  Hence, we need the following statement
% for consistency.

yest(22:26)=yest(22:26)*100/sum(yest(22:26));

% Define the known inputs.  These are from the base-case
% values of Downs and Vogel.  Note that only the first
% 11 of the 35 inputs to TEest6 are known.  The rest will be
% calculated so as to match outputs at steady state.

uest=[11.2208524662
  114.5012151520
   98.0275736993
  418.6853795209
   22.2100000000
   15.0994727278
  259.5017882903
  211.3007075421
  230.3078204048
  120.4000000017
   18.1134905500
   zeros(24,1)];

% The following function calculates the states
% and outputs that allow TEest6 to match TEPROB
% outputs at steady-state.  It will print a 
% mass balance on the screen.

[x_ss,u_ss,y_ss,xvr_ss] = TEest6_ss(yest,uest);  

% The following are diagnostic printouts.
% First display the calculated parameters.

fprintf('\nCalculated Inputs (TEest6 parameters):\n\n')
for i=12:35
  fprintf('%3.0f  %15.10f\n',i,u_ss(i))
end

% Now check whether the outputs of TEest6 match the desired
% values.  NOTE:  With the given parameters, it is impossible
% to match the flow and composition of the stream entering
% the reactor (outputs 7-13 of TEest6).  All others (1-6,14-32) 
% should be matched.  The errors in outputs 7-13 are usually small,
% the order of 1% or less.

fprintf('\n\nComparison of outputs:\n\n')
for i=1:32
  fprintf('%3.0f TEest6 = %15.10f, TEPROB = %15.10f, Diff = %12.8f\n',...
      i,y_ss(i),yest(i),y_ss(i)-yest(i))
end

% Verify that we actually get the correct outputs when we call
% TEest6 with the states and inputs calculated above.

y=TEest6(0,x_ss,u_ss,3);
maxyerr=max(abs(y-y_ss));
if maxyerr > 1e-10
  fprintf('\n Warning:  Inconsistent outputs.  Values are:\n\n')
  for i=1:32
      fprintf('%3.0f  %15.10f  %15.10f  %15.10f\n',i,...
        y(i),y_ss(i),y(i)-y_ss(i))
  end 
else
  fprintf('\n max(abs(y-y_ss)) = %10.2e.  Consistent.\n',maxyerr)
end


% Also verify that it's a steady-state -- all derivatives should
% be very small, of the order of e-14.

dxdt=TEest6(0,x_ss,u_ss,1);
maxdxdt=max(abs(dxdt));
if maxdxdt > 1e-10
  fprintf('\n Warning:  NOT a steady state.  Derivatives are:\n\n')
  for i=1:31
      fprintf('%3.0f  %15.10f\n',i,dxdt(i))
  end
else
  fprintf('\n max(abs(dxdt)) = %10.2e.  Steady state verified.\n',maxdxdt)
end

% END OF TEST SCRIPT
```

The following is a copy of the output written to the screen when running the above test script on Matlab 4.1 for Macintosh.

```matlab
Compressor flow  correction =    -3.62 percent
Compressor power correction =     3.52 percent

          1        2        3        4        5        6
A:    1.00000  0.00000  0.00000  0.48790  0.43402  0.32169
B:    0.00000  0.00000  0.00000  0.00499  0.00443  0.08856
C:    0.00000  0.00000  0.00000  0.50712  0.45111  0.26407
D:    0.00000  1.00000  0.00000  0.00000  0.00117  0.06955
E:    0.00000  0.00000  1.00000  0.00000  0.07246  0.18753
F:    0.00000  0.00000  0.00000  0.00000  0.00884  0.01651
G:    0.00000  0.00000  0.00000  0.00000  0.01981  0.03554
H:    0.00000  0.00000  0.00000  0.00000  0.00815  0.01656
kmol:    11.2    114.5     98.0    418.7    470.7   1904.4

          7        8        9       10       11
A:    0.27108  0.32958  0.32958  0.00000  0.00000
B:    0.11370  0.13823  0.13823  0.00000  0.00000
C:    0.19722  0.23978  0.23978  0.00000  0.00000
D:    0.01073  0.01257  0.01257  0.00223  0.00018
E:    0.17702  0.18579  0.18579  0.13635  0.00848
F:    0.02156  0.02263  0.02263  0.01661  0.00100
G:    0.12382  0.04844  0.04844  0.47314  0.54540
H:    0.08487  0.02299  0.02299  0.37167  0.44494
kmol:  1483.3   1205.0    15.1   263.3    211.3
```

Calculated Inputs (`TEest6` parameters):

```matlab
 12    48.7899733652
 13     0.4985108254
 14   105.1335109913
 15   102.8564839615
 16    -3.7747408443
 17   100.1053234816
 18   100.5239821232
 19     0.0000000000
 20    99.8879615715
 21    98.9675736481
 22     1.6355754505
 23     2.2544283060
 24     1.1170022020
 25   109.5177344776
 26   100.3001455996
 27  -2175.6029934786
 28     0.0062788380
 29    97.6196636041
 30     1.3813877433
 31     5.2493305510
 32     5.3285702384
 33    20.0913058891
 34   104.3616523046
 35  -105.3038384479
```

Comparison of outputs:

```matlab
  1 TEest6 = 2705.0000000464, TEPROB = 2705.0000000464, Diff =   0.00000000
  2 TEest6 =   74.9999999809, TEPROB =   74.9999999809, Diff =   0.00000000
  3 TEest6 = 2633.7278831266, TEPROB = 2633.7278831266, Diff =   0.00000000
  4 TEest6 =   49.9999999957, TEPROB =   49.9999999957, Diff =   0.00000000
  5 TEest6 =   50.0000000045, TEPROB =   50.0000000045, Diff =   0.00000000
  6 TEest6 = 3102.2131372370, TEPROB = 3102.2131372370, Diff =   0.00000000
  7 TEest6 =   42.5179315447, TEPROB =   42.3388488563, Diff =   0.17908269
  8 TEest6 =   32.1691665278, TEPROB =   32.1876130945, Diff =  -0.01844657
  9 TEest6 =    8.8558148155, TEPROB =    8.8932726456, Diff =  -0.03745783
 10 TEest6 =   26.4067510403, TEPROB =   26.3832894332, Diff =   0.02346161
 11 TEest6 =    6.9547721538, TEPROB =    6.8820117573, Diff =   0.07276040
 12 TEest6 =   18.7526435763, TEPROB =   18.7763059474, Diff =  -0.02366237
 13 TEest6 =    1.6506096434, TEPROB =    1.6567273566, Diff =  -0.00611771
 14 TEest6 =   32.9575372908, TEPROB =   32.9575372908, Diff =   0.00000000
 15 TEest6 =   13.8229458674, TEPROB =   13.8229458674, Diff =   0.00000000
 16 TEest6 =   23.9781442398, TEPROB =   23.9781442398, Diff =   0.00000000
 17 TEest6 =    1.2565424847, TEPROB =    1.2565424847, Diff =   0.00000000
 18 TEest6 =   18.5792697025, TEPROB =   18.5792697025, Diff =   0.00000000
 19 TEest6 =    2.2633413851, TEPROB =    2.2633413851, Diff =   0.00000000
 20 TEest6 =    4.8436409545, TEPROB =    4.8436409545, Diff =   0.00000000
 21 TEest6 =    2.2985780751, TEPROB =    2.2985780751, Diff =   0.00000000
 22 TEest6 =    0.0181368255, TEPROB =    0.0181368255, Diff =   0.00000000
 23 TEest6 =    0.8483929880, TEPROB =    0.8483929880, Diff =   0.00000000
 24 TEest6 =    0.1000738732, TEPROB =    0.1000738732, Diff =   0.00000000
 25 TEest6 =   54.5398151317, TEPROB =   54.5398151317, Diff =   0.00000000
 26 TEest6 =   44.4935811817, TEPROB =   44.4935811817, Diff =   0.00000000
 27 TEest6 =   22.9487300898, TEPROB =   22.9487300898, Diff =   0.00000000
 28 TEest6 =   26.9024403996, TEPROB =   26.9024403996, Diff =   0.00000000
 29 TEest6 =  341.4336508454, TEPROB =  341.4336508454, Diff =   0.00000000
 30 TEest6 =   65.7310290341, TEPROB =   65.7310290341, Diff =   0.00000000
 31 TEest6 =   80.1094022593, TEPROB =   80.1094022593, Diff =   0.00000000
 32 TEest6 =   77.2969835300, TEPROB =   77.2969835300, Diff =   0.00000000

 max(abs(y-y_ss)) =   5.43e-11.  Consistent.

 max(abs(dxdt)) =   5.85e-12.  Steady state verified.
```

---

Here's a trivial Fortran program that calls the model to verify a steady-state for the Downs and Vogel base case conditions.
To run it, you will need to compile and link with all subroutines in the `TEest6.f` file EXCEPT the `mexFunction` subroutine.

```fortran
c *******************************************************************

c Example call to the TEsimf subroutine (simplified TE problem model).

c N. L. Ricker, Department of Chemical Engineering, University of
c Washington, ricker@cheme.washington.edu

  implicit double precision (a-h, o-z)
  
  integer nx, ny, nu
  
  parameter (nx=31, ny=40, nu=35)
  
  double precision x(nx), y(ny), u(nu), dxdt(nx)
  
  
c   Define vectors of states and inputs.  The values given
c   below are chosen such that the model outputs match (where possible)
c   those of the Downs and Vogel base case.

  u(1)=  11.220852466199998
  u(2)= 114.501215152000015
  u(3)=  98.027573699300007
  u(4)= 418.685379520899971
  u(5)=  22.210000000000001
  u(6)=  15.099472727800000
  u(7)= 259.501788290299999
  u(8)= 211.300707542099985
  u(9)= 230.307820404799998
  u(10)= 120.400000001700008
  u(11)=  18.113490550000002
  u(12)=  48.789973365214209
  u(13)=   0.498510825435812
  u(14)= 105.133510991278769
  u(15)= 102.856483961462800
  u(16)=  -3.774740844302983
  u(17)= 100.105323481555942
  u(18)= 100.523982123155832
  u(19)=   0.000000000000000
  u(20)=  99.887961571546484
  u(21)=  98.967573648069276
  u(22)=   1.635575450523589
  u(23)=   2.254428306027727
  u(24)=   1.117002201967097
  u(25)= 109.517734477616472
  u(26)= 100.300145599611085
  u(27)=-2175.602993478598364
  u(28)=   0.006278838046563
  u(29)=  97.619663604063504
  u(30)=   1.381387743286999
  u(31)=   5.249330550998074
  u(32)=   5.328570238368596
  u(33)=  20.091305889102607
  u(34)= 104.361652304621074
  u(35)=-105.303838447935192

  x(1)=   4.706892104820388
  x(2)=   1.974149773223099
  x(3)=   3.424483353070932
  x(4)=   0.180387437729033
  x(5)=  10.231473705387032
  x(6)=   1.246406249479704
  x(7)=  66.215392688690144
  x(8)=  67.847366802496580
  x(9)=  28.894829211633251
  x(10)=  12.118977717175163
  x(11)=  21.022334785139851
  x(12)=   0.096802998265292
  x(13)=   5.919711670229184
  x(14)=   0.721143975281563
  x(15)=  20.542538419416164
  x(16)=  16.136872440118314
  x(17)=  51.742769165001889
  x(18)=  14.244210566299735
  x(19)=  42.474163024869526
  x(20)=  11.186462348546154
  x(21)=  30.162848856015202
  x(22)=   2.654937102127430
  x(23)=   5.717063005097626
  x(24)=   2.663395389524711
  x(25)=   0.007149632465752
  x(26)=   0.334440998864720
  x(27)=   0.039449649597940
  x(28)=  21.499883319486074
  x(29)=  17.539604810951840
  x(30)=  -3.466083444426376
  x(31)=  -5.858626637868833

  call TEsimf(nx,nu,ny,x,u,y,dxdt)
  
  dxer=0.0
  do i=1,nx
    dxer=max(dxer,abs(dxdt(i)))
  end do
  write (*,*) 'Maximum abs(dxdt) = ',dxer
  write (*,*) 'Calculated Outputs:'
  do i=1,ny
    write(*,1000) i,y(i)
  end do
  stop
 1000 format(1x,i3,f20.10)
  end
```

---

The above program was compiled and linked with the code in the file `TEest6.f` (Macintosh MPW environment, Language Systems fortran compiler).  When run, it writes the following results to the screen:

```fortran
Maximum abs(dxdt) =    5.722450291401060D-12
Calculated Outputs:
  1     2705.0000000464
  2       74.9999999809
  3     2633.7278831266
  4       49.9999999957
  5       50.0000000045
  6     3102.2131372370
  7       42.5179315447
  8       32.1691665278
  9        8.8558148155
 10       26.4067510403
 11        6.9547721538
 12       18.7526435763
 13        1.6506096434
 14       32.9575372908
 15       13.8229458674
 16       23.9781442398
 17        1.2565424847
 18       18.5792697025
 19        2.2633413851
 20        4.8436409545
 21        2.2985780751
 22        0.0181368255
 23        0.8483929880
 24        0.1000738732
 25       54.5398151317
 26       44.4935811817
 27       22.9487300898
 28       26.9024403996
 29      341.4336508454
 30       65.7310290341
 31       80.1094022593
 32       77.2969835300
 33       81.0469709631
 34      115.9743795103
 35       94.3623250173
 36        0.5532094172
 37      760.6483049192
 38      553.4070889515
 39       30.1109792530
 40      496.7077359268

STOP 
```

---


There are two intended uses for the model:

1. State/parameter estimation.  Determine parameters and states such that the model matches the plant in a defined sense.  Can be done on line in a variety of ways, e.g., an extended Kalman filter.
The resulting model can be used for on-line optimization and/or control.

2. As a simplified version of the TE plant for dynamic simulations. Add your own integration package...
