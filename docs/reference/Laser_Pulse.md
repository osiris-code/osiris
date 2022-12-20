---
layout: single
classes: wide
title: Laser Pulse
permalink: /reference/zpulse
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures "zpulse" laser pulses. One of these sections must
exist for all the required zpulses . Note that contrary to what happens
in the antenna sections it is not necessary to write a preceding section
specifying the number of zpulses to use. The user simply needs to write
a zpulse section for every zpulse required. This section accepts the
following data: (**note that the parameters specified are for the
electric and magnetic field, not the intensity**)

- **type**, character(\*), default = "box"
- **iflaunch**, bool, default = .true.
- **a0**, real, default = 1.0
- **omega0**, real, default = 10.0
- **phase**, real, default = 0.0
- **chirp_order**, integer, default = 0
- **chirp_coefs**, real(4), def ault = 0
- **per_chirp_order**, integer, default = 0
- **per_chirp_coefs**, integer, default = 0
- **pol_type**, integer, default = 0
- **pol**, real, default = 90.0
- **propagation**, character(\*), default = "forward"
- **direction**, integer, default = 1
- **lon_type**, character(\*), default = "polynomial"
- **lon_start**, float, default = 0.0
- **lon_rise**, float, default = 3.1415
- **lon_fall**, float, default = 3.1415
- **lon_flat**, float, default = 0.0
- **lon_duration**, float, default = 0.0
- **lon_x0**, float, default = -huge(1.0)
- **lon_range**, float, default = 0.0
- **lon_math_func**, character(\*), default = "NO_FUNCTION_SUPPLIED!"
- **lon_tilt**, float(2), default = 0.0
- **per_type**, character(\*), default = "plane"
- **per_center**, float(2), default = -huge(1.0)
- **per_w0**, float, default = 90.0
- **per_fwhm**, float, default = 90.0
- **per_focus**. float, default = 0.0
- **per_tem_mode**, integer(2), default = 0,0
- **per_w0_asym**, float(2,2), default = 90.0
- **per_fwhm_asym**, float(2,2), default = 90.0
- **per_focus_asym**. float(2,2), default = 0.0
- **per_asym_trans**, float(2), default = 0.0
- **per_n**, integer, default = 0
- **per_0clip**, integer, default = 1
- **per_kt**, float, default = 0.0
- **time**, float, default = 0.0
- '''no_div_corr, bool, default = .false.

## General Laser Parameters

**type** specifies the type of pulse initialization to use. Valid values
are:

- "box" - The laser pulse is initialized in a single timestep inside the
  simulation box (this is the default)
- "wall" - The laser pulse is launched gradually from the simulation box
  wall.

**time** specifies the time, in normalized units, at which the pulse
will be initialized for "box" type initialization, or start to be
launched for "wall" type initialization. The default is launch (box) or
start launching (wall) the laser pulse at t=0.0.

**a0** specifies the normalized vector potential for the peak intensity
of the laser in units of c.

**omega0** specifies the ratio between the laser frequency and plasma
frequency (i.e. the laser frequency normalized to plasma frequency
units).

**phase** specifies the initial phase of the laser in degrees.

**chirp_order** specifies the order of the spatial chirp to be used on
the pulse. The wavenumber for the injected laser pulse will be in the
form k = k0 + k1 z + k2 z^2 + ... + kN z^N where z is the longitudinal
distance from the center the laser pulse, k0 is the fundamental
wavenumber defined with the omega0 parameter, and k1, k2, etc. are the
chirp coeficients defined with the chirp_coefs parameter, and N is
defined through this parameter (chirp_order). A pulse with chirp_order =
1 will have a linear chirp. The position of the center of the pulse
depends on the longitudinal profile used. Currently the maximum
chirp_order that can be used is 4, but this can be easily changed in the
os-zpulse file. A value of 0 turns off the chirp. The center of the
laser pulse depends on the type of longitudinal profile chosen. Check
the lon_type parameter for details.

**chirp_coefs** specifies the coeficients for the spatial chirp to be
used on the pulse. See chirp_order for details.

**per_chirp_order** specifies the order of the perpendicular spatial
chirp to be used n the pulse. Similar to chirp_order parameter. Note
that in 3D the perpendicular spatial chirp will always be in the first
perpendicular direction (e.g. for a laser pulse launched along x1, the
spatial chirp will only be along x2).

**per_chirp_coefs** specifies the coeficients for the perpendicular
spatial chirp to be used on the pulse. See per_chirp_order for details.

**pol_type** specifies the laser polarization type. Valid values are:

- +1 - clockwize polarization
- 0 - linear polarization
- -1 - counter-clockwize polarization

**pol** specifies the laser polarization in degrees. For pol = 0.0, the
laser pulse is polarized in the x1x2 plane and for 90.0, the laser pulse
is polarized in the x3 plane.

**propagation** specifies the propagation direction for the pulse. Valid
values are: "forward" and "backward", specifying forward and backward
propagation respectively.

**direction** specifies the direction along which to launch the laser
pulse.

## Longitudinal Profile

**lon_type** specifies the type of longitudinal profile to use. Valid
values are:

- *"polynomial"* - gaussian like 5th order symmetric polynomial profile.
  The envelope is given by $10 t^3 - 15 t^4 + 6 t^5$ where $t = (x -
  lon_{start})/lon_{rise}$ for the laser pulse envelope rise, and $t
  = (lon_{start} + lon_{flat} - x)/lon_{fall}$ for the laser pulse
  envelope fall. See the parameters lon_rise, lon_flat, lon_fall and
  lon_start for details. The center of the pulse is automatically
  defined as the center of the flat region.
- *"sin2"* - $\sin^2$ profile. See the parameters lon_rise,
  lon_fall, lon_fall and lon_start for details. The center of the pulse
  is automatically defined as the center of the flat region.
- *"gaussian"* - gaussian profile. See the parameters lon_duration,
  lon_x0 and lon_range for details. The center of the pulse is defined
  as lon_x0.
- *"math"* - analytical math function. See the parameter lon_math_func,
  lon_x0 and lon_range for details. The center of the pulse is defined
  as lon_x0.

**lon_start** specifies the position of the front of the laser pulse,
for the a longitudinal profile of type polynomial, at the moment of its
initialization, in normalized units. This parameter is silently ignored
if the laser is being launched from the wall.

**lon_rise**, **lon_fall**, **lon_flat** specify the time parameters of
the laser pulse for a longitudinal profile of types polynomial and sin2.
lon_rise and lon_fall specify the laser pulse rise and fall times. The
shape of the rise and fall edges of the laser pulse follow either a 6th
order polynomial, with a gaussian like shape if type is set to
polynomial, or a sin^2 function if type is set to sin2. This function
will go from 0.0 to 1.0 in the specified lengths. lon_flat specifies the
time the laser pulse holds its peak value. All values are in normalized
simulation units.

**lon_duration**, **lon_x0**, **lon_range** the time parameters of the
laser pulse for a longitudinal profile of type gaussian. The
longitudinal envelope will be set to \$exp\left( - 2
\frac{(x-lon\\_x0)^2}{lon\\_duration^2}\right)\$ if \|x-lon_x0\| \<
lon_range/2, and 0.0 otherwise. This means the lon_duration parameter
specifies the FWHM for the envelope, lon_x0 the central position and
lon_range the full length where the laser pulse is initialized. When
launching laser pulses from the wall the lon_x0 parameter is silently
ignored.

**lon_math_func**, **lon_x0**, **lon_range** specifies the analytical
function to be used to define the longitudinal profile for the math
longitudinal profile type. This expression must be function of the
variable "x", defined as the distance to the pulse center, defined by
the parameter lon_x0. The longitudinal envelope will be set to zero if
\|x-lon_x0\| \>= lon_range/2. See the documentation on the analytical
function parser for details on the mathematical expression. Note that
this function should only depend on the variable along the axis the
laser is propagating; otherwise the laser pulse being initialized may be
a non-physical one.

**lon_tilt**, specifies a tilt in the longitudinal envelope, in degrees.
In 3D the user can specify a different tilt for each perpendicular
direction. In effect, this parameter results in adding
tan(lon_tilt)\*rho to the parameter of the longitudinal envelope
calculations, resulting in a tilted beam. This has no impact in the
phase information of the pulse.

## Perpendicular Profile

**per_type** specifies the type of perpendicular (transverse) profile to
use. This parameter has no effect in 1D simulations. Valid values are:

- *"plane"* - plane wave . No extra parameters required.
- *"gaussian"* or *"hermite"* - gaussian/ hermite gaussian beam
  (astigmatic beams are also supported). See the parameters per_center,
  per_w0, per_fwhm, per_focus and per_tem_mode for details.
- *"laguerre"* - Laguerre gaussian beam (3D only). See the parameters
  per_center, per_w0, per_fwhm, per_focus and per_tem_mode for details.
- *"bessel"* - Bessel beam (3D only). See the parameters per_center,
  per_n, per_0clip and per_kt for details.
- *"asymmetric"* - Asymmetric gaussian beam. Note that asymmetric
  gaussian beams are only allowed in 2D and 3D simulations. See the
  parameters per_center, per_w0_asym, per_fwhm_asym, per_focus_asym and
  per_asym_trans for details.

**per_center** specifies the position of the optical axis for the beam
in normalized units. If not specified the beam will be launched along
the corresponding central axis of the simulation box at t=0.

**per_w0**/**per_fwhm** specify the parameters for the perpendicular
profile of a hermite/laguerre gaussian beam. The per_w0/per_fwhm
parameters specify the FWHM of the laser spot size at the focal plane.
When using the per_w0 parameter the user is defining the pulse waist
*w0* as defined (for example) in Saleh and Teich, Fundamentals of
Photonics, chapter 3.1.B; when using per_fwhm you are defining the FWHM
of the field. In 3D, for gaussian beams only, the user may specify 2
values for these parameters (*e.g.* per_w0(1:2) = 3.2, 4.0) defining a
beam with astigmatic properties. For a laser propagating along z, the
first value reports to the x waist, *w0x*, and the second to the y
waist, *w0y*. Finally, the per_fwhm parameter has precedence over
per_w0.

**per_focus** specifies the position of the focal plane of the pulse in
normalized units. In 3D, for gaussian beams only, the user may specify 2
values for these parameters (*e.g.* per_focus(1:2) = 0.0, 1.0) defining
a beam with astigmatic properties. For a laser propagating along z, the
first value reports to the position of the x focal plane, and the second
to the y focal plane.

**per_tem_mode** parameter specifies the TEM mode of the
hermite/laguerre-gaussian beam. If not specified it defaults to TEM00,
the fundamental gaussian mode.

**per_w0_asym**/**per_fwhm_asym** have a similar meaning to
**per_w0**/**per_fwhm** for regions above/below the optical axis for
asymmetric gaussian beams. These parameters are in the form per_w0_asym(
zone, coord ), where *zone* can be 1 or 2, specifying the region
below/above the optical axis, and *coord* can be 1 or 2 specifiying the
transverse coordinate. In 2D *coord* must always be 1.

**per_focus_asym** has a similar meaning to **per_focus** for regions
above/below the optical axis for asymmetric gaussian beams. This
parameter is in the form per_focus_asym( zone, coord ), where *zone* can
be 1 or 2, specifying the region below/above the optical axis, and
*coord* can be 1 or 2 specifiying the transverse coordinate. In 2D
*coord* must always be 1.

**per_asym_trans** specifies the size of the transition region between
the different values of per_w0_asym and per_focus_asym. This insures
that there is no phase discontinuity in the transition between different
parameter regions and must be \> 0.

**per_n**, **per_0clip**, **per_kt** specify the parameters for the
perpendicular profile of a bessel beam. The per_n parameter specifies
the order of the bessel function to use. The per_0clip parameter
specifies at wich zero the bessel beam should be clipped. Setting this
value to 0 means that the beam is not to be clipped. The default is 1,
meaning that the bessel beam is clipped at the first zero. This
functionality is only implemented for bessel beams of up to 5th order,
and the beam can only be clipped in one of the first 5 zeros. The per_kt
parameter specifies the transverse wavenumber in simulation parameters
and should be smaller than the longitudinal wavenumber, k (which has the
same numerical value as omega0). When per_kt \> 0.5 k the code will
issue a warning. The definition of the bessel beam follows Saleh and
Teich, Fundamentals of Photonics, chapter 3.4, Wiley & Sons, 1991, and
Hafizi, Esarey and Sprangle, Physical Review E, vol. 44, no. 3, pp.
3539-3545, 1997.

## Other

**no_div_corr** specfies whether the divergence correction algorithm
should be turned off ("box" type initialization only). The algorithm for
injecting laser pulses used by osiris, injects only the transverse
component of the electro-magnetic field. The longitudinal components,
present for example in 0 order bessel beams or tightly focused gaussian
beams, are calculated by solving the equations \$ \nabla \cdot \bf{E} =
\nabla \cdot \bf{B} = 0\$. Setting this parameter to .true. will turn
off the calculation of the longitudinal field components.

## Examples

Here's an example of a zpulse section for a 2D run, initializing a laser
pulse directly on the box:

```text
zpulse 
{
  a0 = 1.0,                  ! use a normalized peak value of 1
  omega0 = 5.0,              ! use a normalized central frequency of 5.0
  pol = 45.0d0,              ! use a 45 deg polarization angle

  chirp_order = 1,           ! add a linear chirp to the pulse
  chirp_coefs(1) = 0.5,      ! the wavenumber will increase 0.5 per unit length

  propagation = "forward",   ! launch the pulse with forward propagation (default)
  
  lon_type = "polynomial",   ! use a polynomial longitudinal envelope
  lon_rise = 6.4,            ! use a 6.4 1/wp rise time
  lon_flat = 0.0,            ! hold the peak value for 0.0 length (i.e. don't hold it)
  lon_fall = 6.4,            ! use a 6.4 1/wp fall time. The pulse will have a total length of 12.8
  lon_start = 0.0,           ! start the pulse exactly at the rigth edge of the simulation 

  per_type = "gaussian",     ! use a gaussian transverse profile
  per_w0 = 1.0,              ! with a spot size (radius) of 1.0
  per_focus = 6.4,           ! place the focal plane at 6.4 from the rigth edge of the simulation
}
```

Here's another example of a zpulse section for a 3D run, launching a
laser pulse from the wall:

```text
 zpulse
 {
     type = "wall",
     direction = 2,
     a0 = 1.0,
     omega0 = 5.0,
     pol = 45.0d0,
     propagation = "backward",
     lon_type = "polynomial",
     lon_rise = 6.4,
     lon_flat = 0.0,
     lon_fall = 6.4,
     lon_start = 12.8,
     per_type = "gaussian",
     per_w0 = 1.0,
     per_focus = 6.4,
     per_center(1:2) = 2.0, 6.4,
 }
````
