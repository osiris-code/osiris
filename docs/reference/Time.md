---
layout: single
classes: wide
title: Time
permalink: /reference/time
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the spatial information and moving window
settings of the simulation and must be present in the input file. It
accepts the following data:

- **tmin**, float, default = 0.0
- **tmax**, float, default = 0.0

**tmin**, **tmax** specify the initial and final time of the simulation
in normalized units.

Here's an example of a time section that sets the initial and final time
of the simulation to 0.0 and 70.0 respectively.

```text
time
{
  tmin = 0.0d0,
  tmax  = 70.d0,
}
```
