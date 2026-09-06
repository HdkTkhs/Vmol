Sep. 2nd, 2026  
Hideaki Takahashi,
Tohoku University, 
Sendai, Japan

# MPI-parallelized Vmol
'Vmol' is a code[1,2] written in Fortran for the electronic density-functional theory (DFT) calculation based on the real-space grid formalism[3]. 
'Vmol' had been combined with classical molecular dynamics codes (produced by third parties) to 
build a QM/MM simulator[2] in the early stage of the development. However, the current release only includes the original MPI-parallelized Kohn Sham-DFT[4] module extracted from the integrated code. Thus, the core programs of the 'Vmol' have been placed on GitHub. 
As a consequence, a lot of functions equipped on the original code have been disabled. In the following, we listed several features of the present distribution.
* Optimized-effective potential for Hartree-Fock method (HF-OEP)[5]  
  'Vmol' at the present distribution is specific to a parallelized HF-OEP calculation[6] using MPI libraries. 
* For the parallel execution of HF calculation[7], f90 module: [poisson_solver.f90](https://github.com/shunsakuraba/poisson_solver) is involved. 
* The 'Vmol' on the GitHub is specifically edited to perform the HF calculation that is followed by HF-OEP. However, it can also be
  used for normal KS-DFT calculations with slight modifications of the code.  
* The external subroutines and functions made by others were also excluded from the original code due to the copyright issues in the distribution. The names of the routines
  and the functions are provided below.
  
**References**  
```
[1] H. Takahashi, T. Hori, T. Wakabayashi, and T. Nitta, “Real space ab initio molecular dynamics simulations  
    for the reactions of OH radical/OH anion with formaldehyde,”  J. Phys. Chem. A 105, 4351 (2001).  
[2] H. Takahashi, T. Hori, H. Hashimoto, and T. Nitta, “A hybrid QM/MM method employing real space grids   
    for QM water in the TIP4P water solvents,” J. Comp. Chem. 22, 1252–1261 (2001).  
[3] J. R. Chelikowsky, N. Troullier, and Y. Saad, “Finite-difference pseudopotential method: electronic structure  
    calculations without a basis,” Phys. Rev. Lett. 72, 1240–1243 (1994).  
[4] W. Kohn and L. J. Sham, “Self-consistent equations including exchange and correlation effects,”  
    Phys. Rev. 140, A1133–A1138 (1965).  
[5] W. Yang and Q. Wu, “Direct method for optimized effective potentials in density-functional theory,”  
    Phys. Rev. Lett. 89, 143002 (2002).    
[6] H. Takahashi, “Comparison of optimized effective potential with inverse Kohn–Sham method for Hartree–Fock  
    exchange energy,” J. Chem. Phys. 161, 104108(11) (2024).  
[7] H. Takahashi, S. Sakuraba, and A. Morita, “Large-scale parallel implementation of Hartree−Fock  
    exchange energy on real-space grids using 3D-parallel fast Fourier transform,”  
    J. Chem. Inf. Model 60, 1376–1389 (2020).
```

# Prerequisites 
* Compilation of 'Vmol' requires [FFTW](https://www.fftw.org) and [pFFT](https://github.com/mpip/pfft). 
* Compilation of 'Vmol' requires Intel Math Kernel Libraries (MKL).
* Execution of 'Vmol' requires pseudopotential database [NCPS](http://www.bandstructure.jp/readmee.html).
* As indicated by an include sentence placed at the end of the 'Vmol' source code, 'ext_routines.f' file (**not** provided in the distribution due to the copyright issues) is required to incorporates the miscellaneous external subroutines and functions into 'Vmol'. The subroutines are those
  provided in 'Numerical Recipes' by William H. Press, et al (Cambridge University Press). Explicitly, the following programs are required;
  1. LSFIT.f (Least-Square Fitting)
  2. GAUSSJ.f (matrix inversion using Gauss Jordan elimination)
  3. POLINT.f (Polynomial interpolation)
  4. DDPOLY.f (Evaluation of derivative of a given polynomial)
  5. POLCOE.f (Polynomial coefficients)
     
  Most of these routines can be replaced with equivalent routines provided by e.g. MKL, and the rest will be developed by users. Anyway, the file 'ext_routines.f', which contains these routines, must be prepared by users. 
* The file 'ext_routines.f' also includes functions: RANFQ.f and GAUSS.f. RANFQ.f generates random real numbers R(0<R<1), while GAUSS.f generates
  random real numbers which forms a Gaussian distribution. RANFQ.f can be equivalently replaced by the intrinsic function 'random_number(r)' in Intel Fortran. 
  GAUSS.f might be also replaced by a function in MKL. Note, however, that KS-DFT calculation does not use GAUSS.f. Thus, the processes related to GAUSS.f can
  be safely omitted for the sole purpose to perform the KS-DFT calculation.

# Input files
Running 'Vmol' requires the following input files;
  1. qm.dat      (computational settings, molecular specifications, etc)
  2. PS_DATA/    (pseudopotential database [NCPS](http://www.bandstructure.jp/readmee.html))    
  3. basis.dat   (LCAO basis set data for constructing initial guess)
  4. valence.dat (LCAO coefficients  )

# Compile and link
Provided that the prerequisites noted above are fulfilled, 'Vmol' can be compiled and linked by invoking 
the script Compile.exe in Vmol_src/ .  
```mpif90 -O -c -I/path/to/FFTW/include -I/path/to/pfft/include -std=f2003 poission_solver.F90```

# Output files

# Citation
Please cite the following articles when you use the part of 'Vmol'.  

```
[1] H. Takahashi, T. Hori, T. Wakabayashi, and T. Nitta,  
    “Real space ab initio molecular dynamics simulations for the reactions  
    of OH radical/OH anion with formaldehyde,” J. Phys. Chem. A 105, 4351 (2001).  
[2] H. Takahashi, T. Hori, H. Hashimoto, and T. Nitta,  
    “A hybrid QM/MM method employing real space grids for QM water in the TIP4P
    water solvents,” J. Comp. Chem. 22, 1252–1261 (2001).  
[3] H. Takahashi, S. Sakuraba, and A. Morita,  
    “Large-scale parallel implementation of Hartree−Fock exchange energy on real-space  
    grids using 3D-parallel fast Fourier transform,” J. Chem. Inf. Model 60, 1376–1389 (2020).  
[4] H. Takahashi,  
    “Comparison of optimized effective potential with inverse Kohn–Sham method  
    for Hartree–Fock exchange energy,” J. Chem. Phys. 161, 104108(11) (2024).  
```

