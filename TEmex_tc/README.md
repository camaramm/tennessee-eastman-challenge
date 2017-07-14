TEmex_tc
=====

> **NOTE:** see updated information at end of file.

#### Contents

- `temex_tc.c`: C source code for compiling TEMEX mex file.
- `temex_tc.h`: C header file needed when compiling.
- `temex_tc.dll`: Mex file ready for use on a Windows machine. It was compiled and built using MS Visual C++ and tested in MATLAB release 11.
- `temex_tc.mex`: Mex file ready for use on a PowerMac. It was built using the MPW MrC compiler and tested in MATLAB 5.2.1.
- `te_test.mdl`: Simulink model illustrating use of the TE challenge process simulation.
- `teplot.m`: Script that plots the results at the end of a simulation.

#### System Requirements:

Matlab Version 5.2 or higher, with Simulink version 2 or higher.

#### Installation

1. Unzip the archive into a new directory. If you are on a Windows computer, you may delete the `.mex` file. If you are on a PowerMac you may delete the `.dll` file. Otherwise you may delete both the `.mex` and `.dll` files.

2. Start MATLAB and make the new directory the default (use the CD command).

3. If you are using a system other than Windows or PowerMac, you must now compile and build the code. You will need a C compiler (probably already installed if you're using a Unix machine). If you have one, type the following in the MATLAB command window:
      ```matlab
      mex temex_tc.c
      ```
      
  If all goes well, after a short delay the compiled-and-built mex file will be in the default directory (verify this). You may see a couple of warning messages during the compile step, but you can ignore them (unless they are errors rather than warnings).

  If you have trouble with the compile & build, make sure you have a C compiler approved for use with MATLAB c-mex files, and that it is set up properly. See MATLAB's "Application Program Interface Guide" (available on-line in PDF format) for more details.

#### Testing the code

4. In the MATLAB command window type `te_test` to bring up a Simulink window containing the example.
  > This assumes that the directory created in step 1 is still the MATLAB default.

5. Start the simulation.

  The initial condition is the base case defined in the Downs and Vogel paper.
  
  In this case, however, the reactor and separator temperatures are under feedback control. This stabilizes the plant with respect to _minor_ disturbances.
  
  See the `TEMEX` archive for the original open-loop unstable case. The variables will be plotted (the format of this  depends on whether or not you have the MPC Tools).
  
  You will note some transients because disturbance \#8 is active (see the definition of the second parameter in the `TEMEX_TC` block).
  
  If you re-run with an increased simulation time you will find that the plant shuts down after about 6 hours, i.e., disturbance 8 is a serious one, and the two temperature control loops are insufficient to keep the plant operating.

---

Please see additional discussion in the `README` file of the `TEMEX` archive.

If you have trouble installing the code or it seems to give incorrect results, please contact N. L. Ricker giving full details of your problem: (ricker@u.washington.edu)

---

#### UPDATE 25 FEBRUARY 2002 -- changes to following files

- `temex_tc.dll`: Renamed to `temex_tc.dll.v5`. A new version of `temex_tc.dll` was created for Matlab version 6.1 (Release 12). To use the original, delete the new one, and remove the appended .v5 from the original's name.
- `temex_tc.c`: Revised source code used to create the new dll file. Major difference is that the `TEMEX_TC` S-function now considers the input to include a direct feedthrough to the output. This was necessary to prevent segmentation errors.