---
layout: single
classes: wide
title: Laser Pulse
permalink: /reference/zpulse
usemathjax: true

sidebar:
  nav: "ref"
---

This section configures `zpulse` laser pulses. The standard `zpulse` algorithm initializes the electric and magnetic fields throughout the simulation domain during a single specified timestep -- typically at `t=0`. This is in contrast to the [`antenna`](Antenna.md), which injects the pulse from a simulation boundary over a range of timesteps. The `zpulse` and `antenna` functionality are in certain ways redundant and are maintained separately for legacy purposes; future versions of OSIRIS may have a single module for launching laser pulses. As of this release, `zpulse` injection is recommended for pulses that fit entirely within the box without significantly extending the computational domain.

In the `zpulse` injection scheme, the full transverse electric field is initialized at a fixed timestep $t_\text{launch}$ as $ \mathbf E_\perp(\xi, z, \mathbf x_\perp) \approx B(\xi) C(z,\mathbf x_\perp) \exp[ i( k_0 z - \omega_0 t + \phi_0  ) ] \mathbf { \hat{\epsilon}} + \text{c.c.}$. The `zpulse` module allows the user to separately specify the pulse's polarization $ \mathbf { \hat{\epsilon}}$, longitudinal profile $B(\xi)$, and transverse (or "perpendicular") profile $C(z,\mathbf x_\perp)$. For example, $B(\xi)$ may be a Gaussian with duration $\tau$ and $C(z,\mathbf x_\perp)$ may be a Hermite- or Laguerre-Gaussian beam mode. The functions $B$ and $C$ are controlled through `lon_type`, `per_type`, and related parameters. The magnetic field is initialized using a similar equation and can also be set using the above equation for $\mathbf E_\perp$ to exactly satisfy the discretized Faraday's law. The longitudinal components of the fields are initialized by numerically integrating to satisfy Gauss's law $\nabla \cdot E = \nabla \cdot B = 0$; in the code, this is referred to as "divergence correction".

The above separable initialization works because if $C(z, \mathbf x_\perp)$ is a solution of the paraxial beam equation $ik_0 \partial_z C + \frac 1 2 \nabla_\perp^2 C = 0$ (as are Hermite- and Laguerre-Gaussian modes), then the above equation for $\mathbf E_\perp(\xi, z, x_\perp) $ is an approximate solution of Maxwell's equations over a length scale $L_\text{stc} = \frac 1 2 \omega_0 \tau k_0 w_0^2$ which typically far exceeds the Rayleigh range $z_R \equiv \frac 1 2 k_0 w_0^2$. Thus, for pulses separably initialized within several Rayleigh ranges of their focal point, the above form approximately holds *even for $t>t_\text{launch}$*. Thus, an ultrashort pulse with a paraxial-beam transverse profile will maintain the properties of the beam -- such as spot size and focal point -- as it propagates forward. See reference \[1\] for more details on separable solutions to Maxwell's equations.

The `zpulse` module should only be used to initialize laser pulses within vacuum regions before they enter the plasma. Initializing the pulse within plasma will result in an unphysical plasma response and an incorrect longitudinal electric field.

The `zpulse` is not designed for Lorentz-boosted simulations, for which pulse injection over a range of timesteps is optimal (see \[2\] for an example). A separate module, the `zpulse_mov_wall`, is designed for injection of pulses into Lorentz-boosted simulations and will be added to the open-source version of Osiris in the future.

Note that is not necessary to write a preceding section specifying the number of `zpulses` to use. The user simply needs to write a separate zpulse section for every `zpulse` required.

Each `zpulse` section accepts the following parameters:

- __a0__, real, default = 1.0
- __omega0__, real, default = 10.0
- __phase__, real, default = 0.0
- __pol_type__, integer, default = 0
- __pol__, real, default = 90.0
- __propagation__, character(\*), default = "forward"
- __direction__, integer, default = 1

<!-- -->

- __if_launch__, bool, default = .true.
- __launch_time__, float, default = 0.0
- __b_type__, character(\*), default = "normal"
- __no_div_corr__, logical, default = .false.

<!-- -->

- __lon_type__, character(\*), default = "polynomial"
- __lon_start__, float, default = 0.0
- __lon_rise__, float, default = 3.1415
- __lon_fall__, float, default = 3.1415
- __lon_flat__, float, default = 0.0
- __lon_fwhm__, float, default = -1.0
- __lon_duration__, float, default = 0.0
- __lon_x0__, float, default = -huge(1.0)
- __lon_range__, float, default = 0.0
- __lon_math_func__, character(\*), default = "NO_FUNCTION_SUPPLIED!"
- __lon_tilt__, float(2), default = 0.0
- __chirp_order__, integer, default = 0
- __chirp_coefs__, real(4), def ault = 0

<!-- -->

- __per_type__, character(\*), default = "plane"
- __per_center__, float(2), default = -huge(1.0)
- __per_focus__. float, default = 0.0
- __per_w0__, float, default = 90.0
- __per_fwhm__, float, default = 90.0
- __per_tem_mode__, integer(2), default = 0,0
- __per_w0_asym__, float(2,2), default = 90.0
- __per_fwhm_asym__, float(2,2), default = 90.0
- __per_focus_asym__. float(2,2), default = 0.0
- __per_asym_trans__, float(2), default = 0.0
- __per_n__, integer, default = 0
- __per_0clip__, integer, default = 1
- __per_kt__, float, default = 0.0
- __per_chirp_order__, integer, default = 0
- __per_chirp_coefs__, integer, default = 0
- __per_tem_mode__, integer(2), default = 0

### Basic parameters

__a0__ specifies the normalized vector potential for the peak intensity
of the laser in units of c. The electric and magnetic fields will have amplitudes of $ E_0 = B_0 = a_0 \omega_0 $ in normalized units. In quasi-3D simulations, the amplitude of each real and imaginary component will be $ \frac 1 2 a_0 \omega_0 $, which combine to produce a total field amplitude of $ E_0 = B_0 = a_0 \omega_0 $.

__omega0__ specifies the ratio between the laser frequency and plasma frequency (i.e. the laser frequency normalized to plasma frequency units).

__phase__ specifies the initial phase of the laser in degrees.

__pol_type__ specifies the laser polarization type. Valid values are:

- +1 - clockwise polarization
- 0 - linear polarization
- -1 - counter-clockwise polarization

__pol__ specifies the laser polarization in degrees. For pol = 0.0, the laser pulse is polarized in the $x_2$ direction and for 90.0, the laser pulse is polarized in the $x_3$ direction.

__propagation__ specifies the propagation direction for the pulse. Valid values are: "forward" and "backward", specifying forward and backward propagation respectively.

__direction__ specifies the direction along which to launch the laser pulse; for example, 1 will launch the pulse in the x1 direction. Currently, only propagation in the x1 direction is supported.

### Launch Parameters

__if_launch__ specifies whether or not to launch the pulse in the simulation.

__launch_time__ specifies the time, in normalized units, at which the pulse will be initialized. The default is to launch the laser pulse at $t=0.0$.

__b_type__ specifies the algorithm used to initialize the magnetic field. Valid values are:

- *"normal"* - the magnetic field is set the same way as the electric field: perpendicular field components are initialized as the above product of the longitudinal and perpendicular profiles.
- *"int"* - The magnetic field is set to satisfy Faraday's law.

__no_div_corr__ specifies whether to skip the divergence correction step which sets the $x_1$ component of the electric and magnetic fields to ensure that Gauss's law holds numerically in the vacuum initialization region, i.e. $\nabla \cdot E = \nabla \cdot B = 0 $.

### Longitudinal Profile

__lon_type__ specifies the type of longitudinal profile to use. Valid values are:

- *"polynomial"* - gaussian-like 5th order polynomial profile. The longitudinal profile is given by $10 u^3 - 15 u^4 + 6 u^5$ where $u = (z - \text{lon_start})/\text{lon_rise}$ for the laser pulse envelope rise and $u = (\text{lon_start} + \text{lon_flat} - z)/\text{lon_fall}$ for the laser pulse envelope fall. Here, $z$ is the absolute position of the spatial coordinate along the propagation axis. See the parameters `lon_rise`, `lon_flat`, `lon_fall` and `lon_start` for details. The start of the pulse is defined as `lon_start`.
- *"sin2"* - $\sin^2$ profile. See the parameters `lon_rise`, `lon_fall`, `lon_fall` and `lon_start` for details. The start of the pulse is defined by `lon_start`.
- *"gaussian"* - gaussian profile. See the parameters `lon_duration`, `lon_x0` and `lon_range` for details. The center of the pulse is defined as `lon_x0`.
- *"math"* - analytical math function. See the parameters `lon_math_func`, `lon_x0` and `lon_range` for details. The center of the pulse is defined as `lon_x0`, meaning that the coordinates are defined using this value as the reference.

__lon_start__ specifies the position of the front of the laser pulse for a longitudinal profile of type "polynomial" or "sin2", at the moment of its initialization, in normalized units.

__lon_rise__, __lon_fall__, __lon_flat__ specify the time parameters of the laser pulse for a longitudinal profile of types "polynomial" and "sin2". `lon_rise` and `lon_fall` specify the laser pulse rise and fall times. The shape of the rising and falling edges of the laser pulse follow either a 6th order polynomial, with a Gaussian-like shape if type is set to polynomial, or a $sin^2$ function if type is set to "sin2". `lon_flat` specifies the time the laser pulse holds its peak value. All values are in normalized simulation units. This function will go from 0 to 1 and back to 0 within the specified initialization region.

__lon_duration__, __lon_x0__, __lon_range__ specify the time parameters of the laser pulse for a longitudinal profile of type Gaussian. The longitudinal envelope will be set to $\exp\left( - 2 (x-\text{lon_x0})^2 / \text{lon_duration}^2 \right)$ if $\|x-\text{lon_x0}\| < \text{lon_range}/2$, and 0.0 otherwise. This means the `lon_duration` parameter specifies the duration for the envelope, `lon_x0` the central position, and `lon_range` the full length where the laser pulse is initialized.

__lon_math_func__, __lon_x0__, __lon_range__ specifies the analytical function to be used to define the longitudinal profile for the math longitudinal profile type. This expression must be a function of the variable "x", defined as the distance to the pulse center, defined by the parameter `lon_x0`. The longitudinal envelope will be set to 0 if $\|x-\text{lon_x0}\| \ge \text{lon_range}/2$. See the documentation on the [analytical function parser](../other/Function_Parser.md) for details on the mathematical expressions. This function should only depend on the variable of the direction in which the laser is propagating.

__lon_fwhm__ allows the user to directly specify the longitudinal FWHM of the pulse for different pulse types. For "polynomial" and "sin2" longitudinal types, specifying `lon_fwhm` is equivalent to setting `lon_rise = lon_fwhm/2`, `lon_flat = 0`, and `lon_fall = lon_fwhm/2`. For a "gaussian" longitudinal type, specifying `lon_fwhm` is equivalent to setting `lon_duration  = lon_fwhm/sqrt( 2.0_p_double * log(2.0_p_double))`. In any case, `lon_fwhm` overrides the other duration parameters. `lon_fwhm` cannot be used when using a "math" longitudinal profile.

__lon_tilt__, specifies a tilt in the longitudinal envelope, in degrees. In 3D the user can specify a different tilt for each perpendicular direction. This parameter results in adding $\tan(\text{lon_tilt}) \|x_\perp\| $ to the $\xi$ parameter of the longitudinal envelope calculation, resulting in a tilted envelope. This has no impact in the phase information of the pulse.

__chirp_order__ specifies the order $n$ of the spatial chirp to be used on the pulse. The wavenumber for the injected laser pulse will be in the form $k = k_0 + k_1 z + \dots + k_n z^n$ where $z$ is the longitudinal distance from the center of the laser pulse, $k_0$ is the fundamental wavenumber defined with the `omega0` parameter, $k_j$ are the chirp coefficients defined with the `chirp_coefs` parameter, and $N$ is defined through the `chirp_order` parameter. For example, a pulse with `chirp_order = 1` will have a linear chirp. Currently, the maximum `chirp_order` that can be used is 4. Longitudinal chirp specified by `chirp_coefs` may be combined with a transverse chirp specified by `per_cherp_coefs`.

__chirp_coefs__ specifies the coefficients for the spatial chirp to be used on the pulse. See `chirp_order` for details.

### Perpendicular Profile

__per_type__ specifies the type of perpendicular (transverse) profile to use. This parameter has no effect in 1D simulations, where only plane waves can be used. Valid values are:

- *"plane"* - plane wave. No extra parameters are required.
- *"gaussian"* or *"hermite"* - Gaussian / Hermite-Gaussian beam (astigmatic beams are also supported). See the parameters `per_center`, `per_w0`, `per_fwhm`, `per_focus` and `per_tem_mode` for details.
- *"laguerre"* - Laguerre gaussian beam (3D only). See the parameters `per_center`, `per_w0`, `per_fwhm`, `per_focus` and `per_tem_mode` for details.
- *"bessel"* - Bessel beam (3D only). See the parameters `per_center`, `per_n`, `per_0clip` and `per_kt` for details.
- *"asymmetric"* - Asymmetric gaussian beam. See the parameters `per_center`, `per_w0_asym`, `per_fwhm_asym`, `per_focus_asym` and `per_asym_trans` for details.

__per_center__ specifies the transverse coordinates of the focal position of the pulse in normalized units. If not specified the pulse will be launched along the corresponding central axis of the simulation box (using initial the box boundaries, if `launch_time > 0`).

__per_focus__ specifies the longitudinal position of the focal plane of the pulse in normalized units. In 3D, for Gaussian beams only, the user may specify 2 values for these parameters (e.g. `per_focus(1:2) = 0.0, 1.0`) to define an astigmatic beam. For a laser propagating along $z$, the first value reports to the position of the $x$ focal plane, and the second to the $y$ focal plane.

__per_w0__, __per_fwhm__ specify the parameters for the perpendicular profile of a Hermite- or Laguerre-Gaussian beam. The `per_w0` and `per_fwhm` parameters specify the laser spot size at the focal plane. When using the `per_w0` parameter the user is defining the pulse waist $w_0$ [as defined here](https://en.wikipedia.org/wiki/Gaussian_beam). When using `per_fwhm` you are defining the FWHM of the field, related by $\text{FWHM} = w_0 \sqrt{4}\log(2)$. In 3D, for Gaussian beams only, the user may specify 2 values for these parameters (*e.g.* `per_w0(1:2) = 3.2, 4.0`) to define an astigmatic beam with separate spot sizes $w_{0x}$ and $w_{0y}$. The `per_fwhm` parameter has precedence over `per_w0`.

__per_tem_mode__ specifies the TEM mode of the Hermite- or Laguerre-Gaussian profile. For Hermite-Gaussian modes, these are the mode numbers of the $x$ and $y$ Hermite polynomial factors. For Laguerre-Gaussian modes, `per_tem_mode = (p,l)`; the second entry is the azimuthal quantum number. If not specified, `per_tem_mode` defaults to `0,0`, which is a Gaussian beam for either Hermite- or Laguerre-Gaussian modes.

__per_w0_asym__, __per_fwhm_asym__, __per_focus_asym__ allow the user to prescribe an asymmetric Gaussian transverse profile with different spot sizes or longitudinal focal points in regions above and below the optical axis. These parameters are of the form `per_w0_asym( zone, coord )` and `per_focus_asym( zone, coord )`, where *zone* can be 1 or 2, specifying the region below/above the optical axis, and *coord* can be 1 or 2 specifying the transverse coordinate. In 2D *coord* must always be 1.

__per_asym_trans__ specifies the size of the transition region between the different values of per_w0_asym and per_focus_asym. This ensures that there is no phase discontinuity in the transition between different parameter regions. `per_asym_trans >  0` is required when using an asymmetric Gaussian beam.

__per_n__, __per_0clip__, __per_kt__ specify the parameters for the perpendicular profile of a Bessel beam. The `per_n` parameter specifies the order of the bessel function to use. The `per_0clip` parameter specifies at which zero the bessel beam should be clipped. Setting this value to 0 means that the beam is not to be clipped. The default is 1, meaning that the Bessel beam is clipped at the first zero. This functionality is only implemented for Bessel beams of up to 5th order, and the beam can only be clipped in one of the first 5 zeros. The `per_kt` parameter specifies the transverse wavenumber in simulation parameters and should be smaller than the longitudinal wavenumber $k$ (which has the same numerical value as `omega0`). When $\text{per_kt} > 0.5 \text{omega0}$ the code will issue a warning. The definition of the Bessel beam follows \[3\] and \[4\].

__per_chirp_order__ specifies the order $n$ of the perpendicular spatial chirp to be used by the pulse. This parameter functions like the `chirp_order` parameter: the local wavevector $k$ is varied as $k = k_0 + k_1 r + \dots + k_n r^n$ where $r$ is the transverse distance from the propagation axis, $k_0$ is the fundamental wavenumber defined with the `omega0` parameter, $k_j$ are the chirp coefficients defined with the `per_chirp_coefs` parameter, and $N$ is defined through the `per_chirp_order` parameter. Note that in 3D the perpendicular spatial chirp will always be in the first perpendicular direction (e.g. for a laser pulse launched along $x_1$, the spatial chirp will only be along $x_2$). Longitudinal chirp specified by `chirp_coefs` may be combined with a transverse chirp specified by `per_cherp_coefs`.

__per_chirp_coefs__ specifies the coefficients for the perpendicular spatial chirp to be used on the pulse. See per_chirp_order for details.

## Examples

Here's an example of a `zpulse`` section for a 2D run:

```text
 zpulse
 {
     a0 = 1.0,
     omega0 = 5.0,
     pol = 45.0d0,
     propagation = "forward",
     lon_type = "polynomial",
     lon_rise = 6.0,
     lon_flat = 0.0,
     lon_fall = 6.0,
     lon_start = -4,
     per_type = "gaussian",
     per_w0 = 4.0,
     per_focus = 6.4,
     per_center(1:2) = 2.0, 6.4,
 }
```

## References

\[1\] Pierce et al., Arbitrarily Structured Laser Pulses, Arxiv, 2022.

\[2\] Yu et al., Enabling Lorentz boosted frame particle-in-cell simulations of laser wakefield acceleration in quasi-3D geometry, Journal of Computational Physics 316, pp. 747–759, 2016

\[3\] Saleh and Teich, Fundamentals of Photonics, chapter 3.4, Wiley & Sons, 1991,

\[4\] Hafizi, Esarey and Sprangle, Physical Review E, vol. 44, no. 3, pp. 3539-3545, 1997.
