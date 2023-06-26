.. TurboRVB_manual documentation master file, created by
   sphinx-quickstart on Thu Jan 24 00:11:17 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

TurboRVB in a nutshell
======================================================

`TurboRVB <https://aip.scitation.org/doi/10.1063/5.0005037>`__ is a computational package for **ab initio Quantum Monte Carlo (QMC) simulations** of both molecular and bulk electronic systems.
The code was initially launched by **Prof. Sandro Sorella** and **Prof. Michele Casula** and has been continuously developed by many contributors for over 20 years. The code implements two types of well established QMC algorithms: Variational Monte Carlo (VMC), and Diffusion Monte Carlo in its robust and efficient lattice regularized variant (LRDMC).

TurboRVB sets itself apart from other QMC codes through several unique features:

- The program utilizes a resonating valence bond (RVB)-type wave function, such as the Jastrow Geminal/Jastrow Pfaffian. This wave function offers the ability to capture correlation effects extending beyond the scope of the commonly employed Jastrow-Slater wave function seen in other QMC codes.
- Incorporating cutting-edge optimization algorithms, like stochastic reconfiguration and the linear method, TurboRVB ensures stable optimization of the amplitude and nodal surface of many-body wave functions at the variational quantum Monte Carlo level.
- The code implements the so-called lattice regularized diffusion Monte Carlo method, which assures a numerically stable diffusion quantum Monte Carlo calculation.
- The integration of an adjoint algorithmic differentiation offers us the capability to efficiently differentiate many-body wave functions, facilitating structural optimizations and the calculation of molecular dynamics.
- TurboGenius and TurboWorkflows allow us to realize high-throughput QMC calculations based on TurboRVB. Both TurboWorkflows and TurboGenius are implemented by Python 3 in the object-oriented fashion; thus, they are readily extended.

When you publish a paper using TurboRVB, please cite the following paper(s).

- `TurboRVB: a many-body toolkit for ab initio electronic simulations by quantum Monte Carlo <https://doi.org/10.1063/5.0005037>`_

   | K. Nakano, C. Attaccalite, M. Barborini, L. Capriotti, M. Casula, E. Coccia, M. Dagrada, Y. Luo, G. Mazzola, A. Zen, and S. Sorella, *J. Chem. Phys.* 152, 204121 (2020)

- TurboGenius: Python suite for high-throughput calculations of the ab-initio quantum Monte Carlo methods

   | K. Nakano, O. Kohulak, A. Raghav, S. Sorella, M. Casula, *in preparation* (2023)
