# Wrapped phase

*(Extracts from the [codedlight
documentation](https://bitbucket.org/nicolasmartin3d/codedlight))*

The phase is now computed using [Chen] :

$$
    \phi(x') = \mathrm{atan}\left(\frac{\sum_{k=0}^s \sin(\dfrac{2 k \pi}{s}) I^c_k(x'))}{\sum_{k=0}^s \cos(\dfrac{2 k \pi}{s}) I^c_k(x'))}\right)
$$


# Phase unwrapping

One of the main issue with this method is that we can only compute
$\phi(x)$ which is the wrapped phase (ie $0 \le \phi(x) \le 2 \pi$).

There are several ways to **unwrap** the phase [Judge]. Spatial
approachs which use some sort of region growing algorithms are not
(yet) implemented in the library.  Instead, two temporal approachs are
implemented in the library.


## Unwrapping from coarse phase

The first one uses an external match map which has unwrapped
(unambiguous) phases but poor precision. This map can be computed with
any other coded light method like **Gray codes**.

It can also be computed with a very low frequency phase shifting
algorithm, even though it is not a really good idea. In this case, the
final phase computed is :

$$
    x =  \left( \left\lfloor \frac{M(x')}{f^{-1}} \right\rfloor + \frac{\phi(x')}{2\pi} \right) f^{-1}
$$

where $M(x')$ is the absolute correspondence of $x'$ given by the
external map.

Note that the performance of this algorithm are greatly dependant on
the level of noise, as this kind of unwrapping is not very robust.

\newpage

## Temporal phase-unwrapping [Huntley]

The other uses several wrapped phases at different frequencies to
compute a final unwrapped phase that is as precise as the phase
computed with the highest frequency.

With this method, the final phase is computed as [Huntley] :

\begin{eqnarray}
\Phi_{n-1}(x') &=& \phi_{n-1}(x') \\
\Phi_{t}(x') &=& \phi_{t}(x') - 2 \pi \, \mathrm{round} \left( \frac{\phi_{t}(x') - \frac{f_t}{f_{t+1}} \Phi_{t+1}(x')}{2 \pi} \right) \\
x' &=& \frac{\Phi_0(x')}{2\pi} f_0^{-1}
\end{eqnarray}

where $f_t$ is the $\mathsf{t}^{\mathsf{th}}$ frequency with $f_0$
being the highest and $f_{n-1}$ the lowest, $\phi_{t}(x)$ the wrapped
phase computed from the $\mathsf{t}^{\mathsf{th}}$ sine wave alone
using [Chen] and $\Phi_{t}(x)$ the absolute phase after $t$ cumulative
unwraps. This method is more robust, since the error can be bounded by
the period difference between two consecutive frequencies, but it
requires more patterns.

\newpage

# References

- [Huntley] Huntley, JM and Saldner, HO, **Error-reduction methods for
shape measurement by temporal phase unwrapping**, *Journal of the
Optical Society of America A*, 1997.

- [Chen] Chen, T and Seidel, HP and Lensch, H, **Modulated
phase-shifting for 3D scanning**, *Computer Vision and Pattern
Recognition*, 2008.

- [Judge] Judge T.R. and Bryanston-Cross P. J., **A review of phase
unwrapping techniques in fringe analysis**, *Optics and Lasers in
Engineering*, 1994

