# Ponderomotive Guiding Center Laser

This section specifies the laser parameters for the ponderomotive
guiding center (PGC) algorithm. For the algorithm details please see D.
Gordon, W.B. Mori and T. Antonsen, IEEE-TPS Vol 28 No 4 August 2000. It
accepts the following data:

- **w0**, real, default = 1.0
- **tau**, real, default = 1.0
- **omega**, real, default = 1.0
- **lon_center**, real, default = 0.0
- **per_focus**, real, default = 0.0
- **per_center**, real, default = 0.0
- **a0**, real, default = 0.0
- **free_stream** logical, default = .false.

**w0** specifies the laser spot size ($1/e$) in the fields.

**tau** specifies the laser duration in the fields. Only polynomial
longitudinal profiles are available (please see zpulse parameters for
additional details on the polynomial).

**omega** specifies the laser central frequency.

**lon_center** specifies the center of the laser pulse in the x1
propagation direction.

**per_focus** specifies the laser focal plane in the x1 propagation
direction.

**per_center** specifies the center of the laser pulse in the direction
perpendicular to the propagation direction (x2).

**a0** specifies peak laser vector potential. A linearly polarized laser
is assumed.

**free_stream** if set to .true. then uses free streaming ponderomotive
guiding center laser pulse.

Here's an example of a pgc section for a 2D run:

```text
pgc{
 w0         = 2.8,
 tau        = 0.25,
 omega      = 185.0,
 lon_center = 19.74,
 per_focus  = 20.,
 per_center = 0.0,
 a0         = 2.0,
 free_stream = .false.,
}
```
