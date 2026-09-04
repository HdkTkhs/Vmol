Sep. 2nd, 2026  
Hideaki Takahashi,  
Tohoku University,   
Sendai, Japan

# Overview 
'Vmol' is a software written with Fortran for the electronic density-functional theory (DFT) based on the real-space grid formalism. 
'Vmol' was originally combined with codes (produced by third parties) for the classical molecular dynamics method to 
build a QM/MM simulator. However, in the current distribution, only the core of the 'Vmol' is being placed on GitHub. 
Thus, only the original Kohn Sham-DFT machinery has been extracted from an integrated code. As a consequence, a lot of 
functions equipped on the code has been disabled. In the following, we listed several features of the present distribution.
* Optimized-effective potential for Hartree-Fock method (HF-OEP)  
  'Vmol' at the present distribution is specific to a parallelized HF-OEP calculation using MPI libraries. 
* For the parallel execution of HF calculation, f90 module: poisson_solver.f90 ([by S. Sakuraba](https://github.com/shunsakuraba/poisson_solver)) is to be invoked.
*   

