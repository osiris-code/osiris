---
layout: single
classes: wide
title: Time step
permalink: /reference/time_step
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the time_step and dump frequency settings and
must be present in the input file. It accepts the following data:

- **dt**, float, default = 0.0
- **ndump**, integer, default = 0
- **dump_start**, integer, default = 0

**dt** specifies the time step in simulation units. It must be different
from 0.0.

**ndump** specifies the number of iterations between any diagnostic or
restart file dumps. If set to 0 all dumps will be turned off. This value
will be multiplied by each of the ndump_fac\_\* parameters in the
restart and diagnostic sections.

**dump_start** specifies the number of iterations after which to begin any diagnostic
file dumps. Note this does not affect restart file dumps.

Here's an example of a time_step section, using a 0.0032 time step and
doing dumps every 200 time steps, where the first two dump events will be ignored
due to the way the `dump_start` parameter is set.

```text
time_step 
{
  dt = 0.0032,
  ndump = 200,
  dump_start = 500,
}
```
