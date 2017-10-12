# `IDVx`
Directory containing data files from a simulation of the decentralized strategy with disturbance `x` of Downs and Vogel.

The disturbance occurs at time = 1 hour and lasts for the remainder of the 50-hour simulation.
Variables are sampled every 10 minutes.
There are 4 data files, organized as follows:

- `y.dat`

	ASCII (TEXT) file containing the measured outputs.
    Each row represents a time instant, starting at time = 0.
    Each column is a measured variable.
    The first 41 are the `XMEAS` variables from the TE model (measured outputs).
    Columns are delimited by TAB characters to make it easier to load the files into plotting and analysis packages (such as Matlab).
    The remaining columns are as follows:

        (42) is for cost [cents/kmol product].
        (43) is production rate of G [kmol G generated/h]
        (44) is production rate of H [kmol H generated/h]
        (45) is production rate of F [kmol F generated/h]
        (46) is partial pressure of A in reactor [kPa]
        (47) is partial pressure of C in reactor [kPa]
        (48) is partial pressure of D in reactor [kPa]
        (49) is partial pressure of E in reactor [kPa]
        (50) is true (delay free) mole % G in product
        (51) is true (delay free) mole % H in product

	These "extra" outputs are not used in the control strategy.
	They may be useful in interpreting the results, however.

- `u.dat`

	As for `y`, but containing the 12 manipulated variable signals (`XMV`) of the TE model.
	All values are in the range 0-100 %.

- `r.dat`

	As for `y`, but containing 36 variables used in the decentralized control strategy.
	These include setpoints for the various loops as well as other signals.
	They are probably not of general interest, but just in case, here's a description.

        (1)   Setpoint for A feed flow (stream 1), kscmh
        (2)   Setpoint for D feed flow (stream 2), kg/h
        (3)   Setpoint for E feed flow (stream 3), kg/h
        (4)   Setpoint for C+A feed flow (stream 4), kscmh
        (5)   Setpoint for purge rate (stream 9), kscmh
        (6)   Setpoint for sep. underflow (stream 10), m^3/h
        (7)   Setpoint for product rate (stream 11), m^3/h
        (8)   Setpoint for reactor pressure, kPa
        (9)   Setpoint for reactor level, %
        (10)  Setpoint for separator level, %
        (11)  Setpoint for stripper level, %
        (12)  No longer used
        (13)  Production rate target (stream 11), m^3/h, supplied by the operator.  Note that this may be overridden under some conditions.  Actual setpoint for stream 11 is in (7), above.
        (14)  Target for mole % G in product (stream 11), as supplied by the operator.  May be modified (rate-of-change constraint). The setpoint used in the feedback loop is (34), below.
        (15)  No longer used
        (16)  Setpoint for %A/(%A + %C) in reactor feed (stream 6), %
        (17)  Setpoint for %A + %C in reactor feed (stream 6), %
        (18)  No longer used
        (19)  Maximum reactor pressure -- setpoint for pressure override, Loop 18 as described in the paper (kPa).
        (20)  Minimum value for separator coolant valve (reactor level override control), %.  See discussion of Loop 19 in the paper.
        (21)  Maximum value for separator coolant valve (reactor level override control), %.  See discussion of Loop 19 in the paper.
        (22)  Ratio setpoint defining variable (1)  = (22) * (32).
        (23)  Ratio setpoint defining variable (2)  = (23) * (32).
        (24)  Ratio setpoint defining variable (3)  = (24) * (32).
        (25)  Ratio setpoint defining variable (4)  = (25) * (32).
        (26)  Ratio setpoint defining variable (5)  = (26) * [(32) - (33)].
        (27)  Ratio setpoint defining variable (6)  = (27) * [(32) - (33)].
        (28)  Ratio setpoint defining variable (7)  = (28) * [(32) - (33)].
        (29)  No longer used
        (30)  Output of reactor pressure override loop.  Used to decrease production rate when reactor pressure is too high.
        (31)  Production rate index.  Nominal value = 100*(var 13)/23.  May be modified to limit rate-of-change.
        (32)  Production rate index after adjusting for overrides, contraint on rate of change, and feedback from production rate loop. Used to determine setpoints 1-7. See formulas for variables (22)-(28), above.
        (33)  Feedback adjustment to production rate index (Loop 8 of paper).
        (34)  Current setpoint used in feedback control of mol % G in stream 11.
        (35)  Eadj value used in eqs. 5 and 6 of the paper.
        (36)  Reactor temperature setpoint, degrees C.


- `t.dat`

	The first column gives the time instants corresponding to each row in `y`, `u` and `r`.
	The other 2 columns in `t.dat` are not of general interest and can be ignored.
