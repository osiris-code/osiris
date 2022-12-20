---
layout: single
classes: wide
title: Input file structure
permalink: /input/
usemathjax: true

sidebar:
  nav: "input"
---

An OSIRIS input file is made out of the following sections, placed in the order defined below. Follow the links for a detailed description of each section. Sections marked with (optional) are not required for an osiris run and can be left out of the input file. In this situation the code will issue a warning message stating that it didn't find the section and is using the defaults. The order of the sections CANNOT be changed and doing so will result in a run-time error.

### General Simulation Parameters

* [simulation](../reference/simulation) - (optional) -  global simulation parameters
* [node_conf](../reference/node_conf) - node configuration, and periodic boundary settings
* [grid](../reference/grid) - grid and coordinate system
* [time_step](../reference/time_step) - time step and dump frequency
* [restart](../reference/restart) - (optional) - restart settings
* [space](../reference/space) - spacial dimensions and moving window
* [time](../reference/time) - time limits

### Electro-Magnetic Fields

* [el_mag_fld](../reference/el_mag_fld) - (optional) - solver type, initial/external fields,
smoothing, etc
* [emf_bound](../reference/emf_bound) - electro-magnetic field boundaries
* [emf_solver](../reference/emf_solver) - (optional) - electro-magnetic field solver parameters
* [smooth](../reference/smooth) - (optional) - electro-magnetic field smoothing
* [diag_emf](../reference/diag_emf) - (optional) - electo-magnetic field diagnostics

### Particles

* [particles](../reference/particles) - number of species, cathodes, neutrals and neutral_mov_ions in the simulation and global particle diagnostics.

The following items are repeated for every species

* [species](../reference/species) - information on the
  particle species
* [udist](../reference/udist) -
  (optional) - momentum distribution for the species
* [profile](../reference/profile) - (optional) - density
  profile information for the species
* [spe_bound](../reference/spe_bound) - boundary
  conditions for the species
* [diag_species](../reference/diag_species) -
  (optional) - species diagnostics

The following items are repeated for every cathode

* [cathode](../reference/cathode) - information on the
  cathode
* [species](../reference/species) - information on the
  particle species to inject electrons from the cathode
* [udist](../reference/udist) -
  (optional) - momentum distribution for the species
* [spe_bound](../reference/spe_bound) - boundary
  conditions for the species holding injected electrons
* [diag_species](../reference/diag_species) -
  (optional) - species diagnostics for the species holding injected
  electrons

The following items are repeated for every neutral

* [neutral](../reference/neutral) - information on
  ionization gas and ion model
* [profile](../reference/profile) - (optional) - density
  profile information for the neutral gas
* [diag_neutral](../reference/diag_neutral) -
  (optional) - neutral gas diagnostics
* [species](../reference/species) - information on the
  particle species to inject electrons from ionization
* [spe_bound](../reference/spe_bound) - boundary
  conditions for the species holding injected electrons
* [diag_species](../reference/diag_species) -
  (optional) - species diagnostics for the species holding injected
  electrons

The following items are repeated for every neutral_mov_ions

* [neutral_mov_ions](../reference/neutral_mov_ions) -
  information on ionization gas and ion model
* [profile](../reference/profile) - (optional) - density
  profile information for the neutral gas
* [diag_neutral](../reference/diag_neutral) -
  (optional) - neutral gas diagnostics
* [species](../reference/species) - information on the
  particle species to inject electrons from ionization
* [spe_bound](../reference/spe_bound) - boundary
  conditions for the species holding injected electrons
* [diag_species](../reference/diag_species) -
  (optional) - species diagnostics for the species holding injected
  electrons
* [species](../reference/species) - information on the
  particle species to inject ions from ionization
* [spe_bound](../reference/spe_bound) - boundary
  conditions for the species holding injected ions
* [diag_species](../reference/diag_species) -
  (optional) - species diagnostics for the species holding injected ions

The following item is specified after all the species, cathode, neutral,
and neutral_mov_ions information.

* [collisions](../reference/collisions) - (optional)
  information on the binary collisions between particles

### Laser Pulses

The following item is repeated for every laser pulse to be launched.
There is no section specifying the number of laser pulses; this is
determined automatically from the number of zpulse sections in the input
deck.

* [zpulse](../reference/zpulse) - (optional) -
  laser pulse parameters

### Electrical Current

* [current](../reference/current) -
  (optional) - dummy section for separating the electrical current section
* [smooth](../reference/smooth) - (optional) - current
  smoothing information
* [diag_current](../reference/diag_current)-
  (optional) - electrical current diagnostics

### Antennas

* [antenna_array](../reference/antenna_array)-
  (optional) - number of antennas to use

The following item is repeated for every antenna

* [antenna](../reference/antenna) - information on the
  antenna to use
