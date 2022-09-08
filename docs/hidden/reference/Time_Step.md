# Time step

This section configures the time_step and dump frequency settings and
must be present in the input file. It accepts the following data:

- **dt**, float, default = 0.0
- **ndump**, integer, default = 0

**dt** specifies the time step in simulation units. It must be different
from 0.0.

**ndump** specifies the number of iterations between any diagnostic or
restart file dumps. If set to 0 all dumps will be turned off. This value
will be multiplied by each of the ndump_fac_\* parameters in the
restart and diagnostic sections.

Here's an example of a time_step section, using a 0.0032 time step and
doing dumps every 200 time steps.

```text
time_step 
{
 dt = 0.0032,
 ndump = 200,
}
```
