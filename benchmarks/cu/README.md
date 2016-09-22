# Benchmark cases for MFIX-Exa


## Overview of benchmark cases

### DEM01: _Homogeneous Cooling System_ [1]

User files create the initial particle layout for all processors and collect/report the granular temperature of the system over time. The run script (runtests) is mainly for local testing but could be modified to generate queue submission scripts.

### DEM02: _Settling_ [1]

Initially, particles are uniformly, randomly distributed throughout the domain with near-zero initial velocity. The gas velocity (if present) is zero everywhere. Particles are allowed to settle under gravity for 0.05 seconds.

### DEM03: _Fluidized Bed_ [1]


### FLD01: _Lid-driven cavity flow_ [2]


## References
1. Peiyuan Liu, Timothy Brown, William D. Fullmer, Thomas Hauser, and Christine Hrenya, "A Comprehensive Benchmark Suite for Simulations of Particle Laden Flows Using the Discrete Element Method with Performance Profiles from the Multiphase Flow with Interface eXchanges (MFiX) Code," Technical Report NREL/TP-2C00-65637, January 2016.

2. J. Musser and A. Choudhary, "MFIX Documentation Volume 3: Verification and Validation Manual," from URL http://mfix.netl.doe.gov
