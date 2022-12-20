---
layout: single
classes: wide
title: Other topics
permalink: /other/
usemathjax: true

sidebar:
  nav: "other"
---

## Analytical Function Parser

OSIRIS includes an analytical function parser so that the user can specify parameters for the simulation in the form of an analytical expression that is evaluated at run time. For details on the capabilities and syntax of the function parser please check the [function parser](function_parser) section.

## Boundary conditions

OSIRIS implements a number of boundary conditions for particles and fields, including periodic, open-space, reflecting and thermal bath boundaries. For a detailed discussion of the boundary conditions available please read the [boundary conditions section](boundary_conditions).

## Shared Memory (OpenMP) parallelism

In addition to the standard spatial decomposition, OSIRIS can use a hybrid distributed memory (MPI) / shared memory (OpenMP) algorithm that relies on the fact that usually the computer nodes used have a number of cores that share the node memory, and may improve performance in some simulations. To read more about shared memory parallelism in Osiris please consult the [shared memory](shared_memory) section.

## Particle Tracking

Particle tracking diagnostics allow you to follow the detailed trajectories (tracks) of individual particles in your simulation. To read more about particle tracking diagnostics in Osiris please consult the [particle tracking](particle_tracking) section.

## Average/Envelope Grids

Average/Envelope grid diagnostics allow you significantly reduce grid diagnostics file size by spatially averaging/taking the envelope of the original data. To read more about particle tracking diagnostics in Osiris please consult the [average](average_envelope) section.

## Grid diagnostics

Diagnostics of grid quantities in OSIRIS share a common interface that allows the user to do several different types of grid diagnostics besides the simple dumping of the grid quantity to disk. The user may choose to do data reduction operations (such as spatial averaging and lineouts) and also to perform time averaging operations to filter out higher (time) frequency components. For a full description of the syntax and available diagnostic types see the [grid diagnostics](grid_diagnostics) section.
