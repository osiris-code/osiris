---
layout: single
classes: wide
title: Antenna
permalink: /reference/antenna
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures antenna settings and must be present in the
input file if any atnennas were specified in the antenna_array section.
One of these sections must exist for every atenna we intend to use. It
accepts the following data:

- **a0**, float, default = 0.0
- **pol**, float, default =0.0
- **side**, integer, default =1
- **delay**, float, default =0.0
- **tilt**, float, default = 0.0
- **phase**, float, default = 0.0
- **focus**, float, default =0.0
- **ant_type**, integer, default = 1
- **t_rise**, float, default = 0.0
- **t_fall**, float, default = 0.0
- **t_flat**, float, default = 0.0
- **omega0**, float, default = 0.0
- **x0**, float, default = 0.0
- **y0**, float, default = 0.0
- **rad_x**, float, defaullt = 0.0
- **rad_y**, float, default = 0.0
- **t_offset**, float, default = 0.0
- **direction**, integer, default = 1
- **chirp**, float, default = 0.0
- **jitterw**, float, default = 0.0
- **jitteromg**, float, default = 1.0

Presently antennas can only be used along the x1 direction. Also note
that using antennas imposes limitation on the boundary conditions you
can use for the electro-magnetic fields. The boundary conditions for the
electro-magnetic fields MUST be Lindmman open-space boundaries (type in
the emf_bound section must be set to 30 for this direction) on both the
lower and upper boundaries along the direction where the EM wave is
emitted. This also means that the remaining boundary conditions CANNOT
be of this type because of the corner limitation mentioned in the
emf_bound section.

**a0** specifies the normalized vector potential for the peak intensity
of the EM wave in units of c (same as vosc in the pulse section)

**pol** specifies the EM wave polarization in degrees at t = 0. For pol
= 0.0, the EM wave is polarized in the x2 direction and for 90.0 the
laser pulse is polarized in the x3 direction.

**side** specifies the side from which to launch the EM wave. The wave
can be launched from the left (1) or from the right (2), corresponding
the the lower or upper boundary in the propagation direction.

**delay** specifies the time in the simulation in normalized units at
which the EM wave will begin to be launched.

**tilt** specifies the angle at which to launch the EM wave in degrees.

**phase** specifies the initial phase of the EM wave in degrees.

**focus** specifies the distance of the focal point of the EM wave to
the antenna in normalized units.

**ant_type** specifies the type of EM Wave to launch. Valid values are

- 1 - launch a pulse with a gaussian envelope
- 2 - launch a pulse with a linear envelope
- 3 - launch a longitudinal (Alfven) wave with a gaussian envelope

**t_rise**, **t_fall**, **t_flat** specify the time parameters of the
laser pulse. t_rise and t_fall specify the laser pulse rise and fall
times. The shapes of the rise and fall are approximately gaussian if
ant_type equals 1 or 3 and linear if ant_type equals 2. t_flat specifies
the time the EM wave holds its peak value. All values are in normalized
simulation units.

**omega0** specifies the ratio between the EM wave frequency and plasma
frequency (i.e. the EM wave frequency in normalized units.

**x0**, **y0** specify the position of the center of the antenna in
normalized units. y0 only affects 3D runs.

**rad_x**, **rad_y** specify the full width of the launched EM wave
(which always has a gaussian transverse profile) in normalized units.
rad_y only affects 3D runs.

**t_offset** t_offset provides a time offset for the envelope function.  In situations where the rise time is very very large but the system is below the LPI threshold for a long time, you can turn on the t_offset so the system does not spend a long period of time under threshold.  In short, the envelope function is advanced by t_offset when this parameter is set to a non-zero value.  t_offset is given in simulation unit.

**direction** (integer) Direction of the Alfven wave propagation.  If direction==1 then the Alfven wave will propagate along the x1 direction, if direction==2 then the Alfven wave will propagate along the x2 direction.  This parameter applies to ant_type==3 only.

**chirp** amplitude of the frequency shift, the time dependent frequency omega(T) = (omega0 + chirp \* T)

**jitterw** amplitude of the jitter frequency

**jitteromg** the frequency of the jitter.  Using both jitterw and jitteromg, you can control the frequency jitter by the formula omega(t) = omega0 + jitterw \* cos(Time \* jitteromg)

Here's an example of an antenna section for a 2D run.

```text
antenna
{
  a0 = 16.0d0,
  t_rise = 1.8831d+02,
  t_fall = 1.8831d+02,
  omega0 = 1.0d0,
  x0 = 1.8850d-01,
  ant_type = 1,
  rad_x = 1000.0d0,
  pol = 90.0
}
````
