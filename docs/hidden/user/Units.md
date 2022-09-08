# OSIRIS Units

## Code Units and Normalization

OSIRIS simulations are done in normalized units. This has 2 distinct
advantages: i) multiplication by several constants (like $m_e$, $e$ and
$c$, for example) is avoided, resulting in a significant performance
increase and ii) by expressing the simulation quantities in terms of the
fundamental plasma quantities the results are general and not bound to
some specific units we may choose. In our case we chose to normalize our
quantities to $\omega_p$, $m_e$, $c$ and $e$, the electron plasma
frequency, the electron rest mass, the speed of light and the electron
charge, respectively.

### Simulation Units

| Quantity       | Conversion                                                                                        |
|----------------|---------------------------------------------------------------------------------------------------|
| Time           | $t' = t \times \omega_p$                                                                          |
| Frequency      | $\omega' = \omega / \omega_p$                                                                     |
| Position       | $\mathbf{x}' = \frac{\omega_p}{c} \mathbf{x}$                                                     |
| Momenta        | $\mathbf{u'} = \frac{\mathbf{p}}{m_{sp} c} = \frac{\mathbf{\gamma v}}{c} = \frac{\mathbf{u}} {c}$ |
| Electric field | $\mathbf{E'} = e \frac{c / \omega_p}{m_e c^2} \mathbf{E}$                                         |
| Magnetic field | $\mathbf{B'} = e \frac{c / \omega_p}{m_e c^2} \mathbf{B}$                                         |

Where $m_{sp}$ is the mass of the species being considered. In this
situation, the relativistic Lorentz factor $\gamma$ can be calculated as
$\gamma = (1 + u'^2)^{1/2}$. Density will be normalized to the plasma
density that yields the reference plasma frequency i.e. a plasma with a
density of 1 in simulation units will have a plasma frequency
$\omega_p$.

### Laboratory Units

In practical units the physical quantities will relate to simulation
quantities as:

| Quantity       | Conversion                                                                                             |
|----------------|--------------------------------------------------------------------------------------------------------|
| Position       | $\mathbf{x} [\mathrm{cm}]=  2.998 \times 10^{10}\, \mathbf{x}' \, \omega_p ^{-1} [\mathrm{rad / s}]$   |
|                | $\mathbf{x} [\mathrm{cm}]=  0.532 \times 10^{6}\, \mathbf{x}' \, n_0^{-1/2} [\mathrm{cm} ^{-3}]$       |
| Momenta        | $\mathbf{p} [\mathrm{g \, cm / s}]=  2.731 \times 10^{-17}\, \frac{m_{sp}}{m_e} \mathbf{u}'$           |
| Electric field | $\mathbf{E} [\mathrm{GV/cm}] = 1.704 \times 10 ^{-14}\, \mathbf{E'} \, \omega_p [\mathrm{rad / s}]$    |
|                | $\mathbf{E} [\mathrm{GV/cm}] = 9.613 \times 10 ^{-10}\, \mathbf{E'} \, n_0 ^{1/2} [\mathrm{cm} ^{-3}]$ |
| Magnetic Field | $\mathbf{B} [\mathrm{gauss}] = 5.681 \times 10 ^{-8}\, \mathbf{B'} \, \omega_p [\mathrm{rad / s}]$     |
|                | $\mathbf{B} [\mathrm{gauss}] = 3.204 \times 10 ^{-3}\, \mathbf{B'} \, n_0 ^{1/2} [\mathrm{cm} ^{-3}]$  |

Also note that for high kinetic energies, $p \gg m_{sp} c$, the
relativistic energy is reduced to $W \simeq m_{sp} c^2 p'$, where
$m_{sp} c^2$ is the rest mass energy for the species being considered.

### Normalizing to another frequency

Alternatively, the user may chose another frequency such as the laser
frequency, $\omega_0$ as the norm. In this case just replace $\omega_p$
with $\omega_0$ in the above formulae. A simulation plasma density of
1.0 will correspond to the critical density.
