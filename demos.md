# USING THE MEX FILE IN MATLAB

> For general information on S-functions, see the SIMULINK manual.

In addition to the standard input arguments (time, states, inputs, flag), `TE_MEX` requires 4 parameters, as follows:

- `p1`: `idv(20)` -- a vector of length 20 containing the desired values of the disturbance flags (see paper by Downs and Vogel).

The following are all used by the TE routines to save data between Matlab calls. They are initialized by a call to this function with `time<0` and `flag=0`.

You must create Matlab variables of the appropriate size (a vector of zeros is convenient), and pass them to `te_mex` as input arguments.
Once you have initialized these variables, DO NOT modify them during the simulated plant operation.
- `p2`: Vector of length 153
- `p3`: Vector of length 586
- `p4`: Vector of length 139

Here are some example calls from Matlab:

```matlab
t=-1;
x=[];
u=[];
flag=0;
idv=zeros(20,1);
p2=zeros(153,1);
p3=zeros(586,1);
p4=zeros(139,1);
[sys,x0]=te_mex(t,x,u,flag,idv,p2,p3,p4);
      
u=[specify a vector, length 12];
flag=1;
t=0;
dxdt=te_mex(t,x0,u,flag,idv,p2,p3,p4);
flag=3;
y=te_mex(t,x0,u,flag,idv,p2,p3,p4);
```

The *first call* to `te_mex` (with `t=-1`) initializes the storage vectors, `p2` to `p4`, and returns the (default) initial state vector in `x0`. This is a good test to verify that the code is working (especially for those of you who have compiled and installed it yourselves).
You can compare `x0` to the values listed for `YY(1)` to `YY(50)` in the `TEINIT` subroutine supplied by Downs and Vogel.  Doing this, for example, helped me to figure out how to port the code from a Mac, where it was developed, to a VAX. Our installation requires the `G_Float` option (see Matlab/VAX manual), and the `x0` vector was incorrect unless this was set properly for the compilation (and the `fmexg.com` file was used for linking).

The *second call* to `te_mex` (with `flag=1`) calculates the time-derivatives of the states at the initial conditions (`x0`) and for a set of inputs specified by the "`u`" vector.  The derivatives are returned in the "`dxdt`" variable.

The *third call* (with `flag=3`) returns the values of the outputs (all 41 of them - see paper by Downs and Vogel) for the given `x0` and `u`.  As explained in the
paper, when `t>0`, measurement noise is added to the outputs automatically.  
Also, the composition analyzers are only sampled at certain times.  Thus,
proper specification of the time variable is important.

The typical use of `TE_MEX.MEX` is to simulate operation of the plant for
a given set of conditions.  

For example, the following sequence of Matlab statements shows how to reproduce the simulations of Figure 3 of Downs and Vogel.
The SIMULINK function, `RK45`, is used to integrate the differential equations for a period of 0.3 hours.  See the SIMULINK manual for details on this function.

```matlab
% Initialize variables:

idv=zeros(20,1);
p2=zeros(153,1);
p3=zeros(586,1);
p4=zeros(139,1);
[sys,x0]=te_mex(-1,[],[],0,idv,p2,p3,p4); 
u0=[63.53,  53.98, 24.644, 61.302, 22.21, ...
    40.064, 38.1,  46.534, 47.446, 41.106, ...
    18.114, 50];       % Base case inputs (XMV variables).
u0(10)=38;             % Specifies a step in XMV(10)

% Set up for integration:

tvec=[0; 0.3];         % Starting and ending time.
ut=[tvec [u0;u0] ];    % Specifies constant inputs.
options=[.001, .00001, .001, 0, 0, 0, 0];

% Integrate over entire 0.3 hour time period:

[tt,xt,yt]=rk45('te_mex',tvec,x0,options,ut,idv,p2,p3,p4);

% Plot results: 

clf
subplot(221)
plot(tt,yt(:,6)),title('Reactor Feed'),xlabel('Time (hours)')
subplot(222)
plot(tt,yt(:,7)),title('Reac Pressure'),xlabel('Time (hours)')
subplot(223)
plot(tt,yt(:,8)),title('Reac Level'),xlabel('Time (hours)')
subplot(224)
plot(tt,yt(:,9)),title('Reac Temp'),xlabel('Time (hours)')
pause

clf
subplot(221)
plot(tt,yt(:,11)),title('Prod Sep Temp'),xlabel('Time (hours)')
subplot(222)
plot(tt,yt(:,12)),title('Prod Sep Level'),xlabel('Time (hours)')
subplot(223)
plot(tt,yt(:,10)),title('Purge Rate'),xlabel('Time (hours)')
subplot(224)
plot(tt,yt(:,18)),title('Stripper Temp'),xlabel('Time (hours)')
```

I found that a specification of a maximum step size of about 0.001 was necessary to get reliable accuracy from the `RK45` function (see use of "options" variable, above).  Other settings may give better performance.  I didn't experiment with it much. Elapsed time for the integration was about 1 minute on a Mac IIfx.

The following additional plot commands show how the concentrations from the analyzers vary with time.  You can observe the discrete sampling behavior, e.g., the product concentrations vary every 0.25 hours, and the others vary every 0.1 hour.

```matlab
clf
subplot(221)
plot(tt,yt(:,23)),title('A in Feed'),xlabel('Time (hours)')
subplot(222)
plot(tt,yt(:,30)),title('B in Purge'),xlabel('Time (hours)')
subplot(223)
plot(tt,yt(:,40)),title('G in Product'),xlabel('Time (hours)')
subplot(224)
plot(tt,yt(:,41)),title('H in Product'),xlabel('Time (hours)')
```

Many other uses of `TE_MEX.MEX` should be obvious to those familiar with the analysis tools in SIMULINK and the various Matlab toolboxes.

# USING THE SOURCE CODE TO COMPILE/LINK YOUR OWN `.MEX` FILE:

The linking procedure is machine-dependent. See your Matlab manual.  I have verified that the code works on a Macintosh and a VAX. I would expect it to work on other machines, but one can never be sure. You're on your own if you attempt it (or modify the source code).

The following assumes Matlab 3.x, so is now outdated.

Some hints follow:

On the Macintosh, I used MPW 3.2 with Language Systems Fortran 3.0. I compiled `TEPROB.F` and `TE_MEX.F` using the following options:

```bash
-i '{MexDir}' -mc68020 -mc68881
```

The symbol `{MexDir}` was set as follows:

```bash
Set MexDir      "{Boot}MATLAB:MEX:"
```

You may need to modify this, depending on your directory structure. The Link step is tricky.  I used the following command:

```bash
  Link -o "te_mex.mex" -d -f -ad 4 -rt 'MEX0'=0 -sg MEX -m mexmain
        -t MEX0 -c MATL -srt -br on
       "{MexDir}MPW:cmex.c.o"
       "{MexDir}MPW:LS:fmex.c.o"
       'te_mex.f.o'
       'teprob.f.o'
       "{FLibs}FORTRANLib.o"
       "{FLibs}IntrinsicLibFPU.o"
       "{FLibs}Interface.o"
       "{FLibs}Runtime.o"
```

with the symbol `{FLibs}` set to "`{Boot}MPW:Libraries:Libraries:FLibraries:`".
This assumes that you have the compiled object files `te_mex.f.o`, and 
`teprob.f.o`, stored in the default directory, and it creates `te_mex.mex` in 
the same directory. Modify as appropriate. The key to making it work is the "`-br on`" option, which allows linking of large code segments. It is apparently unavailable in MPW prior to version 3.2.
