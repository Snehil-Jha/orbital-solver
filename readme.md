# Units
We are working in Hatree atomic units:
$$
    \hbar = m_e = e = 4\pi\epsilon_0 = 1
$$

## Consequences:

Hydrogenic Energy Levels:
$$
    E_n = - \frac{1}{2} \frac{Z^2}{n^2} \\
$$

Classical Atomic Radii:
$$
    r_n = \frac{n^2}{Z}
$$

# Coordinate System
We have employed a logarithmic grid
$$x = ln(r)$$

## Consequences
### Wavefunction
If the physical wavefunction is given by
$$\Psi_{nlm} = R_{nl} Y_{lm}$$

Then the appropriate wavefunction for this coordinate system is
$$\phi(x) = \sqrt{r} R(r)$$

### Schrodinger Equation
The Schrodinger Equation transforms into
$$
    \left[-\frac{1}{2} \partial_{xx} + \frac{1}{8} + \frac{l(l+1)}{2} + e^{2x}V(e^x)\right] \phi = E e^{2x} \phi
$$

### Scaled Electron Density
If the spherical charge density is such that
$$
    N_e = \int_0^{\infty} 4\pi r^2 n(r) dr
$$

Then the scaled density used in the algorithm is
$$
    \rho(x) = 4\pi e^{3x} n(e^x)
$$

Note that we also have
$$
    n(r) = \frac{||R(r)||^2}{4\pi}
$$

So, we ultimately get
$$
    \rho(x) = e^{2x} \cdot ||\phi(x)||^2
$$

### Exchange Potential
Under the Local (Spin) Density Approximation, we have
$$
    V_{x, \sigma}(r) = - \left(\frac{6}{\pi}\right)^{1/3} n_{\sigma}(r)^{1/3} = - \left(\frac{3}{2\pi^2}\right)^{1/3} \frac{\rho(x)^{1/3}}{r}
$$