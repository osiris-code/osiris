# OSIRIS User guide

## Code Units and Normalization

OSIRIS simulations are done in normalized units. For details on the units and normalization used in OSIRIS please read the [units section](Units.md) of the user guide.

## Boundary conditions

OSIRIS implements a number of boundary conditions for particles and fields, including periodic, open-space, reflecting and thermal bath boundaries. For a detailed discussion of the boundary conditions available please read the [boundary conditions section](Boundary_Conditions.md) of the user guide.

## Shared Memory (OpenMP) parallelism

In addition to the standard spatial decomposition, OSIRIS can use a hybrid distributed memory (MPI) / shared memory (OpenMP) algorithm that relies on the fact that usually the computer nodes used have a number of cores that share the node memory, and may improve performance in some simulations. To read more about shared memory parallelism in Osiris please consult the [shared memory](Shared_Memory.md) section of the User Guide.

## Diagnostics

### Particle Tracking

Particle tracking diagnostics allow you to follow the detailed trajectories (tracks) of individual particles in your simulation. To read more about particle tracking diagnostics in Osiris please consult the [particle tracking](Particle_Tracking.md) section of the User Guide.

### Average/Envelope Grids

Average/Envelope grid diagnostics allow you to significantly reduce grid diagnostics file size by spatially averaging/taking the envelope of the original data. To read more about particle tracking diagnostics in Osiris please consult the [average](Average_Envelope_Grids.md) section of the User Guide.
