---
layout: single
classes: wide
title: Electric Current
permalink: /reference/current
usemathjax: true

sidebar:
  nav: "ref"
---

This section is a placeholder for current-specific initializations that may be necessary in the future. As of the current release, all parameters related to current belong in other namelists, namely [`smooth`](Smooth.md) and [`diag_current`](Electric_Current_Diagnostics.md). The `current` namelist should be empty. This section of the input deck should be structured as follows:

```text
current 
{
    ! no parameters (section still needed to use smooth or diag_current)
}

smooth
{
    ! optional current smoothing parameters
}

diag_current
{
    ! optional current diagnostics
}

```
