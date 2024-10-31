# Mechanically Driven Polymers

## Model
We model the polymer as a chain of N connected bonds with fixed length lb. The tangent of bond i is ti ≡ (ri+1 − ri)/lb, where ri is the position of the joint connecting bonds i−1 and i. We fix one end of the polymer at the origin. The polymer energy is given by:

$$
E = \sum_{i=0}^{N-2} \frac{\kappa}{2}\frac{(\mathbf{t}_{i+1} - \mathbf{t}_i)^2}{l_b} - \sum_{i=0}^{N-1} (\gamma z_i + f)(l_b \mathbf{t}_i \cdot \mathbf{x})
$$

where $\kappa$ is the bending modulus, $f$ is the stretching force applied in the $x$-direction, $\gamma$ is the shear ratio along the $z$ direction,
$z_i = \mathbf{r}_i\cdot\mathbf{z}$ is the $z$ component of the position of joint $i$, and
$(\mathbf{t}_i \cdot \mathbf{x})$ is the $x$ component of the bond tangent $\mathbf{t}_i$. A hard sphere interaction between polymer joints, with a sphere radius $l_b/2$, was used to account for self-avoidance of the polymer.

## Mechanical response of the polymer
Ref: [Off-Lattice Markov Chain Monte Carlo Simulations of Mechanically Driven Polymers](https://arxiv.org/abs/2409.15223)

In this work, we present an off-lattice model for sampling the configuration space and calculating the conformational properties of semiflexible polymers under mechanical deformation using Markov Chain Monte Carlo simulations. Our model resolves the orientational bias rooted in on-lattice model and enables precise calculation of polymer’s response to external force. Three states of the system have been studied, in the absence of external field, un- der an uniaxial stretching and in a steady shear. We demonstrate that the deformation of the polymer is captured by the scattering functions with expected anisotropic patterns and the polymer’s conformation in response to external forces calculated using our model is in excellent agreement with theoretical prediction.


## Machine Learning inversion from scattering
Ref: [Machine Learning Inversion from Scattering for Mechanically Driven Polymers](https://arxiv.org/abs/2410.05574)

In this work, we apply a ML inversion method to extract feature parameters from the scattering data of mechanically driven polymers. The ML inversion framework was trained based on the theoretically calculated data set of polymer system that is determined by the energy parameters: bending modulus $\kappa$, stretching force $f$ and shear rate $\gamma$. The inversion targets included these energy parameters and conformation variables such as end-to-end distance $R^2$, radius of gyration $R_g^2$ and off-diagonal component of the gyration tensor $R_{xz}$. The scattering function $I_{xz}(\mathbf{Q})$ of the polymer under different energy parameters was calculated using a MC method we previously developed\cite{ding2024off}. We demonstrate the feasibility of the ML inversion by carrying out PCA of the data set $\mathbf{F} = \left\{ I_{xz}(\mathbf{Q})\right\}$ and investigate the distribution of feature parameters by projecting the data set $\mathbf{F}$ to a 3 dimensional singular vectors space. The GPR was trained and validated, showing that inversion of the feature parameters can be achieved with high-precision.
