# Electro-Magnetic Field Solver

This section configures the Electro-Magnetic field solver settings and is
optional. If not present the code will use the default values. It
accepts the following data if one uses the "stencil" solver in the Electro Magnetic field
settings (otherwise, it accepts no data):

- **k1**, **k2**, real, default = 0.0

**k1**, **k2**, specify the parameters for the "stencil" field solver.
Interesting values for these parameters are k1 \< 0, and k2 = 0 or k2 =
2 k1. Setting k1 = k2 = 0 recovers the Yee solver, and setting k1 =
-1/8, k2 = 0, implements a 4th order accurate spatial derivative
approximation.
