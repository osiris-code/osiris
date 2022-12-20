---
layout: single
title: OSIRIS Reference Guide

toc: true

sidebar:
  nav: "ref"
---

## Input file structure

An Osiris input file is made out of the following sections, placed in the order defined below. Follow the links for a detailed description of each section. Sections marked with (optional) are not required for an osiris run and can be left out of the input file. In this situation the code will issue a warning message stating that it didn't find the section and is using the defaults. The order of the sections CAN NOT be changed and doing so will result in a run-time error.

### General Simulation Parameters

* [simulation](ref/Simulation.md) - (optional) -  global simulation parameters
* [node_conf](ref/Node_Configuration.md) - node configuration, and periodic boundary settings
* [grid](ref/Grid.md) - grid and coordinate system
* [time_step](ref/Time_Step.md) - time step and dump frequency
* [restart](ref/Restart.md) - (optional) - restart settings
* [space](ref/Space.md) - spacial dimensions and moving window
* [time](ref/Time.md) - time limits

### Electro-Magnetic Fields

* [el_mag_fld](ref/Electro-Magnetic_Field.md) - (optional) - solver type, initial/external fields, subcycling and
  smoothing
* [emf_bound](ref/Electro-Magnetic_Field_Boundaries.md) - electro-magnetic field boundaries
* [pgc](ref/Ponderomotive_Guiding_Center_Laser.md) - (optional) - ponderomotive guiding center laser
* [smooth](ref/Smooth.md) - (optional) - electro-magnetic field smoothing
* [diag_emf](ref/Electro-Magnetic_Field_Diagnostics.md) - (optional) - electo-magnetic field diagnostics

### Particles

* [particles](ref/Particles.md) - number of species, cathodes, neutrals and neutral_mov_ions in the simulation and global particle diagnostics.

The following items are repeated for every species

* [species](ref/Species.md) - information on the
  particle species
* [udist](ref/Momentum_Distribution.md) -
  (optional) - momentum distribution for the species
* [profile](ref/Profile.md) - (optional) - density
  profile information for the species
* [spe_bound](ref/Species_Boundary.md) - boundary
  conditions for the species
* [diag_species](ref/Species_Diagnostics.md) -
  (optional) - species diagnostics

The following items are repeated for every cathode

* [cathode](ref/Cathode.md) - information on the
  cathode
* [species](ref/Species.md) - information on the
  particle species to inject electrons from the cathode
* [udist](ref/Momentum_Distribution.md) -
  (optional) - momentum distribution for the species
* [spe_bound](ref/Species_Boundary.md) - boundary
  conditions for the species holding injected electrons
* [diag_species](ref/Species_Diagnostics.md) -
  (optional) - species diagnostics for the species holding injected
  electrons

The following items are repeated for every neutral

* [neutral](ref/Neutrals.md) - information on
  ionization gas and ion model
* [profile](ref/Profile.md) - (optional) - density
  profile information for the neutral gas
* [diag_neutral](ref/Neutrals_Diagnostics.md) -
  (optional) - neutral gas diagnostics
* [species](ref/Species.md) - information on the
  particle species to inject electrons from ionization
* [spe_bound](ref/Species_Boundary.md) - boundary
  conditions for the species holding injected electrons
* [diag_species](ref/Species_Diagnostics.md) -
  (optional) - species diagnostics for the species holding injected
  electrons

The following items are repeated for every neutral_mov_ions

* [neutral_mov_ions](ref/Neutrals_with_Moving_Ions.md) -
  information on ionization gas and ion model
* [profile](ref/Profile.md) - (optional) - density
  profile information for the neutral gas
* [diag_neutral]ref/(Neutrals_Diagnostics.md) -
  (optional) - neutral gas diagnostics
* [species](ref/Species.md) - information on the
  particle species to inject electrons from ionization
* [spe_bound](ref/Species_Boundary.md) - boundary
  conditions for the species holding injected electrons
* [diag_species](ref/Species_Diagnostics.md) -
  (optional) - species diagnostics for the species holding injected
  electrons
* [species](ref/Species.md) - information on the
  particle species to inject ions from ionization
* [spe_bound](ref/Species_Boundary.md) - boundary
  conditions for the species holding injected ions
* [diag_species](ref/Species_Diagnostics.md) -
  (optional) - species diagnostics for the species holding injected ions

The following item is specified after all the species, cathode, neutral,
and neutral_mov_ions information.

* [collisions](ref/Collisions.md) - (optional)
  information on the binary collisions between particles

### Laser Pulses

The following item is repeated for every laser pulse to be launched.
There is no section specifying the number of laser pulses; this is
determined automatically from the number of zpulse sections in the input
deck.

* [zpulse](ref/Laser_Pulse.md) - (optional) -
  laser pulse parameters

### Electrical Current

* [current](ref/Electric_Current.md) -
  (optional) - dummy section for separating the electrical current
  pfield section
* [smooth](ref/Smooth.md) - (optional) - current
  smoothing information
* [diag_current](ref/Electric_Current_Diagnostics.md)-
  (optional) - electrical current diagnostics

### Antennas

* [antenna_array](ref/Antenna_Array.md)-
  (optional) - number of antennas to use

The following item is repeated for every antenna

* [antenna](ref/Antenna.md) - information on the
  antenna to use
