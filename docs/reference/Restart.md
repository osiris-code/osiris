---
layout: single
classes: wide
title: Restart
permalink: /reference/restart
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the restart settings and is optional. If not
present the code will use the default values. It accepts the following
data:

- **ndump_fac**, integer, default = 0
- **debug_iter**, integer, default = -1
- **if_restart**, bool, default = .false.
- **if_remold**, bool, default = .false.
- **ndump_time**, real, default = 0.0

**ndump_fac** - specifies the frequency at which to write restart
information. This value is multiplied by the ndump value specified in
the time_step section to determine the number of iterations between each
restart dump. If set to 0 the writing of restart information is
disabled.

**debug_iter**, sets a specific iteration to write restart information.
This is the exact value of the iteration, and is meant to be used for
debugging purposes.

**if_restart** - specifies whether the code should attempt to read
information from restart files previously saved in order to restart the
run exactly as it was at the time the restart information was saved.

**if_remold** - specifies whether the code should remove older restart
files after it has successfuly saved restart information on all the
nodes. This is extremely usefull in saving disk space.

**ndump_time** - Instead of specifying the frequency at which to write restart information
by using with `ndump_fac` as listed above, this parameter allows one to specify the number of seconds between each restart dump. Note, this requires setting
`ndump_fac=-1`.

Here's an example of a restart section that will write restart
information every 5 \* ndump iterations, will not attempt to use restart
information from previously saved restart information, and will delete
old restart files as it writes new ones.

```text
restart
{
  ndump_fac = 5,
  if_restart=.false.,
  if_remold=.true.,
}
```
