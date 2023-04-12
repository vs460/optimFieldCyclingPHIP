# optimFieldCyclingPHIP

Optimization starts by running main.m

The length of the B0 field increase can be set by the variable T in main.m

By default the simulation performs the following steps prior to optimization based on the paper doi.org/10.1021/acs.jpcb.5b06222:
- Step 1: density matrix averaging at Earth's magnetic field for 1s
- Step 2: diabatic transport to zero field
- Step 3: density matrix averaging at zero field for 1s

These are followed by the optimization of the adiabatic passage that returns B0 from 50 nT to 2 uT.
Step 1-3 can be swithced ON and OFF by setting pars.prepPhase = false in main.m
