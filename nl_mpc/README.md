nl_mpc
=====

This directory contains code related to the nonlinear MPC of the Tennessee Eastman industrial challenge process as described in [*Nonlinear Model Predictive Control of the Tennessee Eastman Challenge Process*, **Computers & Chemical Engineering**, Vol. 19, No. 9, pp. 961-981 (1995)](https://doi.org/10.1016/0098-1354(94)00105-W).
The following is a brief description of each file.

- `TENLMPC.m`: Matlab script that coordinates the calculations. To run it, you will need the files in this directory PLUS the following resources:
  1. The MPC Toolbox from Mathworks
  2. The Control Toolbox from Mathworks
  3. `TEest3_ss.m` from the `TEest3` directory in this directory
  4. `TEest3` MEX file, which you must create by compiling and linking the Fortran code in the `TEest3` directory.
  5. `TE_mex_TC` MEX file, which you must create by compiling and linking the Fortran code in the *Basic TEC Code* or *MATLAB 7.x Code* area.

> **NOTE:** Those running on Windows, MATLAB 6.0 or later, can use the files `TEEST3.DLL` and `te_mex_tc.dll` as items (4) and (5) respectively. These were created by Weilu Lin (mailto:wlin@che.utexas.edu) of the University of Texas with Digital Visual Fortran V 6.0. Weilu was able to compile the Fortran code with no problems.

> **WARNING:**  By default, the simulation will be for 72 hours. This will take many hours to run on most computers. You may wish to reduce the value of `tend` in the script, at least for your first trials.

See comments in the script to understand how to vary setpoints, turn disturbances on and off, etc.  

- `TEEST3.DLL`:  See NOTE above.
- `te_mex_tc.dll`: See NOTE above.
- `est_plt1.m`: does custom plotting for the TE problem. It is called at the end  the simulation.
- `flowfilt.m`:  filters measured flowrates using an exponential filter.
- `linmod_nlr.m` is a variant of the SIMULINK `linmod` function for linearization of a nonlinear state-space model.
- `nlr_euler.m`: Crude Euler (fixed-step) integration of the ODEs in the plant and the `TEest3` model.  You might want to try a better integrator if you have one.
- `PI_feedback.m`: Adds linear PI feedback controllers to a linear plant.
- `plotpair.m`: Plots pairs of variables.

The following 4 ASCII data (`.dat`) files are loaded by the script. They define the estimator gain, and pointers used by the estimator.

- `idest.dat`
- `iyest.dat`
- `iymeas.dat`
- `Kest1.dat`

The following 6 ASCII data files define optional steady-state initial conditions for the plant.  You can load the desired one by setting the "Mode" variable in the script.  See additional comments there.

- `mode1.dat`
- `mode2.dat`
- `mode3.dat`
- `mode4.dat`
- `mode5.dat`
- `mode6.dat`

The following are special versions of code in the MPC toolbox, intended for nonlinear MPC and other custom applications. They are  similar to the `SCMPC` function in the toolbox. A brief explanation follows. For more details, see the comments in the code, and calling sequence in the `TENLMPC` script.

- `smpcqp0.m`: Defines an MPC controller.  Inputs are the linear state-space model of the plant (in MOD format), and the MPC tuning parameters.  Outputs are constant matrices defining the controller, which are needed as inputs to MPCQP.  Note that if all the inputs to this function are constant throughout your simulation, you only need to call this function once.
- `mpcqp.m`: Solves the quadratic programming problem in MPC.  Requires previous call to the SMPCQP0 function.  Call MPCQP at eac sampling period to get the latest control action to be applied to the plant.
- ` SMPCQP0S.m`:  As for `SMPCQP0`, but allows SOFT constraints on outputs.
- `mpcqps.m`: As for `MPCQP`, but for problems with soft output constraints.
