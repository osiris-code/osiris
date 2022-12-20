---
layout: single
classes: wide
title: Antenna Array
permalink: /reference/antenna_array
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures the number of antennas to use in the simulation
is optional. If not present the code will not use any antennas. It
accepts the following data:

- **n_antenna**, integer, default = 0

**n_antenna** - specifies the number of antennas to use in the simulation.

This section must be followed by n_antenna antenna sections describing
each individual antenna.

Here's an example of an antenna_array section setting the number of
antennas to use to 1.

```text
antenna_array{
  n_antenna = 1,
}
```
