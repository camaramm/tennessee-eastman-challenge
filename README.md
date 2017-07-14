# tennessee-eastman-challenge
Tennessee Eastman Challenge Problem Archive

> **Disclaimer:**
> 
> This public repository was created based on the author's webpage (http://depts.washington.edu/control/LARRY/TE/download.html), which is freely available to anyone with internet access. Similarly, the original zip files and source codes are equally freely available on the author's webpage.
>
> That said, this repository is intended to serve as a host to the Tennessee Eastman Challenge studies publicly available, providing the development platform needed to spread, improve and share the source code efficiently.
>
> All links to download files in this *readme* file refer to the original urls in the author's webpage.
>
> Some *readme* files have been edited a little with the purpose of better displaying its content.
>
> The author, Prof. N. L. Ricker, has given permission to the creation of this repository.

## Table of contents

1. [Author](#author)
2. [Links to topics](#table-of-contents)
3. [General comments](#general-comments)
4. [Updated TE Code](#updated-te-code)
5. [New Simulink models of two decentralized control strategies](#new-simulink-models-of-two-decentralized-control-strategies)
6. [MATLAB 7.x Code](#matlab-7x-code)
7. [Basic TEC Code](#basic-tec-code)
8. [Rivera Group’s MATLAB Simulation](#rivera-groups-matlab-simulation)
9. [Optimal steady states](#optimal-steady-states)
10. [Nonlinear MPC](#nonlinear-mpc)
11. [Simplified TE process](#simplified-te-process)
12. [Decentralized control](#decentralized-control)
  1. [Control strategy simulation](#control-strategy-simulation)
  2. [Zip archives](#zip-archives)
13. [Using the TE code in Matlab](#using-the-te-code-in-matlab)

Author
---------

> N. Lawrence Ricker
>
> Professor Emeritus, Chemical Engineering
>
> University of Washington
>
> Box 351750 Seattle, WA 98195-1750 USA
>
> ricker@u.washington.edu
>
> Source: http://depts.washington.edu/control/LARRY/TE/download.html
 
General comments
-----

1. Many of the codes were written for Matlab 3.x and Matlab 4.x.  Only a few have been updated to work with newer versions.  I will consider updating on a case-by-case basis if there is sufficient demand.
2. If you would like me to link to your own TE archive, send the URL and a description of your site.

*E-mail  your comments and questions to  [N. Lawrence Ricker](mailto:ricker@u.washington.edu)*
 
Updated TE Code
-----
Original url link: [temexd_mod.zip](http://depts.washington.edu/control/LARRY/TE/temexd_mod.zip)
(version:  January 23, 2015)


A new version of the C code has been developed to allow the use of variable-step integration methods in MATLAB/Simulink and other simulation codes. 
The use of such methods with Downs and Vogel’s original Fortran code (or the C translation provided in [Basic TEC Code](#basic-tec-code)) leads to inconsistent results.
The updated code has been tested with MATLAB/Simulink 2014b.

The updated code also has the following new features:

 - Documentation of most internal variables.  This violates Downs and Vogel’s intent to keep model details secret, but a careful reading of their original code allows one to decipher it.  As more than 20 years have elapsed since the problem’s release, it no longer seems necessary to keep the details obscure.
 - Additional disturbances and measurements.  These may be of use but are optional.

Note that disturbance switches (IDV values) must be supplied as model input signals, not parameters.  

For more details, download the following archive, and see its `README.TXT` file.
The archive contains a draft paper submitted to ADCHEM 2015 that explains the rationale for the updated model and other technical details.  The archive also contains several demonstrations of open-loop and closed-loop operation with variable-step-size integration.
These may be compared to similar demonstrations using the original code.

[Back to Topics](#table-of-contents)

New Simulink models of two decentralized control strategies
-----

As of 2 December 2002, two new Simulink models have been added to the `temex` archive ([see below](#matlab-7.x-code)).
For more details, download the [archive](http://depts.washington.edu/control/LARRY/TE/temex.zip), and see its `README.TXT` file.
Additional changes made on 24 February 2005.

[Back to Topics](#table-of-contents)
 
MATLAB 7.x Code
------

Two of the key interface routines were revised in December, 1998, to operate under MATLAB 5.2 or greater in Simulink 2 or greater.
Click the links below to download the corresponding Zip archive.
See the README file therein for details on installation and use.
The following comments supersede the instructions in the README file.

Original url link:
- [temex](http://depts.washington.edu/control/LARRY/TE/temex.zip)
- [temex_tc](http://depts.washington.edu/control/LARRY/TE/temex_tc.zip)


> **NOTE:**  these files were updated 23 February 2002 and tested successfully in MATLAB 6.1 (Release 12).

> **NOTE:**  the temex  archive was updated 2 December 2002, and is compatible with MATLAB 6.5 (Release 13). It was tested on 24 Feb 2005 and seemed to work with MATLAB 7.04 (Release 14).

> **NOTE:**  the temex archive was updated on 2 February 2005. The Simulink models for the Multi-loop strategies now require Simulink version 6.2 or higher.

[Back to Topics](#table-of-contents)

Basic TEC Code
------

Original url link: [`tecode.zip`](http://depts.washington.edu/control/LARRY/TE/tecode.zip)

Contents of the `tecode` directory, all in ASCII format:

- `teprob.f`: Fortran code provided by Tennessee Eastman

- `tecommon.inc`: An "include" file needed to compile `TEPROB.F`

- `te_mex.f`: Fortran code needed to generate the .MEX interface to Matlab.*FOR MATLAB 3.x ONLY!!*

- `te_mextc.fv4`: This modification of `te_mex.f` was written for Matlab 4.2c on a MacIntosh. It has not been tested on other machines. It includes PI control of the reactor temperature and separator temperature. There are also modifications to the output vector. See comments in the code for more information. The state vector for this version has 52 elements. The last 2 are the integrated errors for the 2 PI controllers. You must initialize these properly for the system to be at steady state.

[Back to Topics](#table-of-contents)

Rivera Group's MATLAB Simulation
-----

Marty Braun, working with Prof. Daniel Rivera (ASU, Chem Engr) converted Down and Vogel's code from Fortran to MATLAB.  This allows you to work in the convenient Simulink environment without having to compile a MEX file.  You can also modify the MATLAB code so as to make the simulation more amenable to your specific research goals.  

Here is a link to their site, where you can download the code: (http://www.eas.asu.edu/~csel/Software-TennEast.htm)

The disadvantages of this approach are:

- It runs much slower than the equivalent MEX code (see "temex" above). In one test involving a PC running MATLAB Release 12, the MEX version was about 100 times faster.
- Numerical results will not match those of Downs and Vogel.  Braun was unable to duplicate the random noise and disturbance generation of the original code. The results are very similar, however, and the simulation's overall accuracy has been endorsed by Downs and Vogel.

[Back to Topics](#table-of-contents)

Optimal steady states
-----

Original url link: [`tables.zip`](http://depts.washington.edu/control/LARRY/TE/tables.zip)

Extract the following 3 files, which are in tab-delimited text (ASCII) format. It should be possible to read them into a spreadsheet program or text editor.  They contain both text (row labels and units) and numerical data, so you may need to re-format before reading into Matlab, etc. The files give operating  conditions at 6 different steady state operating conditions described in the paper [*Optimal Steady state Operation of the Tennessee Eastman Challenge Process*, N. L. Ricker, **Computers & Chemical Engineering**, Vol. 19, No. 9, pp.  949-959 (1995)](https://doi.org/10.1016/0098-1354(94)00043-N).

- `table2.txt`: First 38 state variables. To generate the complete state vector (50) variables, you must append the 12 manipulated variables (from `table4.txt`) to the end of the states in `table2.txt`.

- `table3.txt`: Output variables (41)

- `table4.txt`: Manipulated variables (12)

To verify a steady state I suggest that you call subroutine tefunc (part of [`teprob.f`](#basic-tec-code)) to calculate the derivatives for the given states and manipulated variables.  All derivatives should  be less than `1.e-3`. You can also use the `te_mex` function in Matlab.

[Back to Topics](#table-of-contents)

Nonlinear MPC
-----

The following zip archives contain Matlab and Fortran code described in the papers:
> - [*Nonlinear Model Predictive Control of the Tennessee Eastman Challenge Process*, **Computers & Chemical Engineering**, Vol. 19, No. 9, pp. 961-981 (1995)](https://doi.org/10.1016/0098-1354(94)00105-W), and
> - [*Nonlinear Modeling and State Estimation for the Tennessee Eastman Challenge Process*, **Computers & Chemical Engineering**, pp. 983-1005 (1995)](https://doi.org/10.1016/0098-1354(94)00113-3).

Click on the link to download the archive.

> **NOTE:**  Since 2002, some researchers have had trouble getting this software to work.  This appears to be caused by changes in MATLAB. I have been unable to correct this problem. If you encounter problems with instability, etc., and are able to correct them, please let me know.

- [`teest3.zip`](http://depts.washington.edu/control/LARRY/TE/teest3.zip): Directory containing code related to the above papers. After expanding the archive, see the file `teest3.doc` for more information on these files.

- [`teest6.zip`](http://depts.washington.edu/control/LARRY/TE/teest6.zip): Directory containing code related to a newer version of the `teest3` model.  After expanding the archive, see the file example.doc for more information.

- [`nl_mpc.zip`](http://depts.washington.edu/control/LARRY/TE/nl_mpc.zip): Code to implement nonlinear MPC. After expanding the archive, see `contents.txt` for more information.

[Back to Topics](#table-of-contents)

Simplified TE process
-----

The following directory contains Matlab and Fortran code described in the paper [*MPC of a continuous, nonlinear, two-phase reactor*, N. L. Ricker, **J. Process Control**, vol. 3, 109-123 (1993)](https://doi.org/10.1016/0959-1524(93)80006-W).  Click on the link to download the archive.  Use Win Zip or equivalent to expand it.

- [`te4.zip`](http://depts.washington.edu/control/LARRY/TE/te4.zip): Expand the archive, then see the file contents.txt for more information on the contents of the TE4 directory.

[Back to Topics](#table-of-contents)

Decentralized control
-----

#### Control strategy simulation

See Simulink code available in the [`temex`](#matlab-7x-code) archive.

#### Zip archives

The following zip archives document the performance of the decentralized control strategy described in [*Decentralized control of the Tennessee Eastman Challenge Process*, N. L. Ricker, **J. Proc. Cont.**, Vol. 6, pp. 205-221 (1996)](https://doi.org/10.1016/0959-1524(96)00031-5).
Each archive is a recording of the response to a particular disturbance (numbered 1-15 in the original problem description of Downs and Vogel).  Click on an archive to download it.

[`idv1.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv1.zip): Step in A/C feed ratio in stream 4

[`idv2.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv2.zip): Step in B composition in stream 4

[`idv3.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv3.zip): Step in D feed temperature (stream 3)

[`idv4.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv4.zip): Step in reactor cooling water inlet temperature

[`idv5.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv5.zip): Step in condenser cooling water inlet temperature

[`idv6.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv6.zip): Sudden loss of A feed (stream 1). This is a tough one!

[`idv7.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv7.zip): Stream 4 header pressure loss. (step change)

[`idv8.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv8.zip): Random variations in A,B,C compositions in stream 4 (another tough one).

[`idv9.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv9.zip): Random variations in D feed temperature

[`idv10.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv10.zip): Random variations in C feed temperature

[`idv11.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv11.zip): Random variations in reactor cooling water inlet temperature

[`idv12.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv12.zip): Random variations in condenser cooling water inlet temperature

[`idv13.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv13.zip): Slow drift in reaction kinetics (also difficult)

[`idv14.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv14.zip): Reactor cooling water valve sticking

[`idv15.zip`](http://depts.washington.edu/control/LARRY/TE/IDVs/idv15.zip): Condenser cooling water valve sticking

Click [HERE](http://depts.washington.edu/control/LARRY/TE/IDVs/format.txt) to see a text file describing the data format of the above archives. You can save this for future reference.

[Back to Topics](#table-of-contents)

Using the TE code in Matlab
-----

Click [HERE](http://depts.washington.edu/control/LARRY/TE/demos.txt) to see a text file containing additional suggestions on using the Basic (original) TE code with MATLAB.

> **NOTE:**  The information has not been updated since the release of MATLAB 5.  You may need to change some aspects to make it work in the latest version of MATLAB.

[Back to Topics](#table-of-contents)
