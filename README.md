Sep. 2nd, 2026  
Hideaki Takahashi,
Tohoku University, 
Sendai, Japan

# MPI-parallelized Vmol
'Vmol' is a code[1,2] written in Fortran for the electronic density-functional theory (DFT) calculation based on the real-space grid formalism[3]. 
In the early stage of the development, 'Vmol' had been combined with classical molecular dynamics codes (produced by third parties) to 
build a QM/MM simulator[2]. However, the current release only includes the original MPI-parallelized Kohn Sham-DFT[4] module extracted from the integrated code. Thus, the core programs of the 'Vmol' have been placed on GitHub. 
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
* Compilation of 'Vmol' requires [FFTW](https://www.fftw.org) and [PFFT](https://github.com/mpip/pfft). 
* Compilation of 'Vmol' requires Intel Math Kernel Libraries (MKL).
* Execution of 'Vmol' requires pseudopotential database [NCPS](http://www.bandstructure.jp/readmee.html).
* As indicated by an include sentence placed at the end of the 'Vmol' source code, 'ext_routines.f' file (**not** provided in the distribution due to the copyright issues) is required to incorporates the miscellaneous external subroutines and functions into 'Vmol'. The subroutines are those
  provided in 'Numerical Recipes' by William H. Press, et al (Cambridge University Press). Explicitly, the following programs are required;
  1. LSFIT.f (Least-Square Fitting)
  2. GAUSSJ.f (matrix inversion using Gauss Jordan elimination)
  3. POLINT.f (Polynomial interpolation)
  4. DDPOLY.f (Evaluation of derivative of a given polynomial)
  5. POLCOE.f (Polynomial coefficients)
     
  Most of these routines can be replaced with equivalent routines provided by e.g. MKL, and the rest can be readily developed by users. Anyway, the file 'ext_routines.f', which contains these routines, must be prepared by users. 
* The file 'ext_routines.f' also includes functions: RANFQ.f and GAUSS.f. RANFQ.f generates random real numbers R(0<R<1), while GAUSS.f generates
  random real numbers which forms a Gaussian distribution. RANFQ.f can be equivalently replaced by the intrinsic function 'random_number(r)' in Intel Fortran. 
  GAUSS.f might be also replaced by a function in MKL. Note, however, that KS-DFT calculation does not use GAUSS.f. Thus, the processes related to GAUSS.f can
  be safely omitted for the sole purpose to perform the KS-DFT calculation.

# Usage
The specification of the real-space cell containing the uniform grids is provided in the include file 'QMpara.i'. The form of the grid is assumed to be cubic in our implementation. 
The grid width is defined by the cut off energy:COE in the input file 'qm.dat'. Shown below is the 'QMpara.i' file for a HF-OEP calculation of a water molecule. 
```
!     "QMpara.i"
!     QM parameter set for Vmol package

!     NNUC:    number of atoms
!     NLINK1:  number of link ( hydrogen ) atoms for QM/MM simulations (disabled)
!     NMAXX:   number of grid points along the x-axis, and similarly for the y- and z- axes.
!     NDIM:    number of dimensions
!     NORA:    total number of orbitals for alpha spin
!     NORB:    total number of orbitals for beta spin
!     MORA:    total number of occupied orbitals for alpha spin
!     MORB:    total number of occupied orbitals for beta spin
!     ITMAX:   maximum number of SCF iterations

      integer*4 NNUC, NLINK1, NMAXX, NMAXY, NMAXZ, NDIM
      integer*4 NORA, NORB, MORA, MORB, NBASIS, NCORE
      integer*4 ITMAX
      integer*4 NFUZZY
      real*8 EPSOPT, EPSMD

      PARAMETER ( NNUC =     3 ) 
      PARAMETER ( NLINK1 =   0 ) 
      PARAMETER ( NFUZZY =   3 )    ! number of fuzzy cells, usually NFUZZY = NNUC
      PARAMETER ( NMAXX =  120 )
      PARAMETER ( NMAXY =  120 )
      PARAMETER ( NMAXZ =  120 )
      PARAMETER ( NDIM =     3 ) 
      PARAMETER ( NORA =    22 )
      PARAMETER ( NORB =     1 )    ! set at 1 for spin-restricted calculations
      PARAMETER ( MORA =     4 ) 
      PARAMETER ( MORB =     1 )    ! set at 1 for spin-restricted calculations
      PARAMETER ( NBASIS =  92 )    ! number of LCAO basis functions for constructing initial guesses 
      PARAMETER ( NCORE =    2 )    ! number of core electrons
      PARAMETER ( ITMAX =  300 )
      PARAMETER ( EPSOPT= 5.0D-5 )
      PARAMETER ( EPSMD = 5.0D-3 )

      PARAMETER (maxatm = 16000)

      character(*), PARAMETER :: PS_DIR = "/home3/takahasi/PS_DATA"     ! directory that contains pseudopotential database 'NCPS'
```

The include file 'mpi.i' specifies the numbers of divisions of the rectangular QM cell along the x,y, and z directions to define subdomains for MPI parallel calculation. For example, for the parallel calculation with 16(=4x2x2) CPUs, the 'mpi.i' file becomes,   
```
!     "mpi.i"
!     parameters for mpi 

!     NX,NY,NZ: numbers of divisions of the rectangular QM cell along the x,y, and z directions

      PARAMETER ( NX = 4 ) ! NX  = 2  >= NY >= NZ
      PARAMETER ( NY = 2 ) ! NY >= NZ
      PARAMETER ( NZ = 2 ) 
```
Note that the definition of the subdomains is directly related to the specification and requirement of the poisson_solver. See also the site of the [poisson_solver.f90](https://github.com/shunsakuraba/poisson_solver).  

# Input files
Running 'Vmol' requires the following input files;
  1. qm.dat      (computational settings, molecular specifications, etc)
  ```
  $INIDAT
 COE   = 100.0      ! cut off energy
 NDEN  = 7          ! number of dense grids
 DTMD  = 41.3411054611601
 MDMAX = 35000
 TEMP  = 300.0
 NRVLC = 50
 NCHK  = 5000
 CONV  = 1.0D-5
 NRST  = 'NEW'
 NOPT  = 'ONE'
 EXC   = 'RHF'      ! RHF calculation followed by HF-OEP
 NQMMM = 'QM'
 PRINT = 'LARGE'
 DGF   = 'DG4'      ! 4th-order Lagrange interpolation is used in the double grid method
 FREEZE= 'TRUE'
 NMRDF = 0 
 NLINK = 0 
 /
 MMID,ZA,ZVAL,SIG,EPSQM    
   1   0   0   8    6    5.788D0  2.8778D-4     0.000000    0.000000    0.223395 
   2   0   0   1    1    5.326D0  0.7869D-4     0.000000    1.427096   -0.893581
   3   0   0   1    1    5.326D0  0.7869D-4     0.000000   -1.427096   -0.893581
      
! index  dummy dummy  # of electrons  # of val. electrons   sigma   epsilon    x   y   z

! Units are in atomic units.

! The LJ parameters are somewhat arbitrary because this input file is for one-point calculations 
! at isolation. Note, however, the LJ-sigma parameters are used in the OPTFC.f or OPTFC1.f routines
! to optimize the fractional charge on each atom. The charges are used to determine the boundary 
! conditions of the Hartree potential. Thus, the choice of the atomic size is not critical for the
! KS-DFT calculations. 
  ```
  3. PS_DATA/    (directory that contains pseudopotential database [NCPS](http://www.bandstructure.jp/readmee.html))    
  4. basis.dat   (LCAO basis set data for constructing initial guesses of the wave functions )  
     The basis set data can be provided by conducting Gaussian16 program suite with an option 'gfinput'. The data will
     be extracted from the Gaussian output file to form the 'basis.dat' file. When one uses 'aug-cc-pVTZ' basis set
     (Ref: T.H. Dunning Jr., J. Chem. Phys. 90, 1007, (1989)), you will have
     the data starting with;    
     ```
           1 0    
      S   7 1.00       0.000000000000   
      0.1533000000D+05  0.5198089434D-03   
      0.2299000000D+04  0.4020256215D-02   
      0.5224000000D+03  0.2071282673D-01   
      0.1473000000D+03  0.8101055358D-01   
      0.4755000000D+02  0.2359629851D+00   
      0.1676000000D+02  0.4426534455D+00      
      ........   
     ```
  6. valence.dat (LCAO coefficients for constructing initial guesses of the wave functions)  
     The LCAO coefficients for the basis set noted above can be obtained by
     invoking the option 'punch=mo' in the SCF calculation using Gaussian 16. Note that the data for
     the core electrons must be excluded from the file.  
     ```
     (5D15.8)  
     2 Alpha MO OE=-0.13531360D+01   
     -0.20132123D+00 0.48727616D+00 0.13269231D+00 0.29345991D+00 0.34914213D-01  
      0.00000000D+00 0.00000000D+00-0.44761941D-01 0.00000000D+00 0.00000000D+00  
     -0.63536184D-01 0.00000000D+00 0.00000000D+00-0.17243387D-01 0.00000000D+00  
      0.00000000D+00-0.28189134D-02 0.97194357D-03 0.00000000D+00 0.00000000D+00  
     -0.36183378D-02 0.00000000D+00 0.36449177D-02 0.00000000D+00 0.00000000D+00  
      ........   
     ```
     The files 'basis.dat' and 'valence.dat' are used only to generate the initial guess wave functions
     in subroutine 'trwf3.f'. Thus, these files can be omitted if another approach is available. 

# Compile and link
Provided that the prerequisites noted above are fulfilled, 'Vmol' can be compiled and linked by invoking 
the csh script 'Compile_Vmol_oep.exe' in Vmol_src/. Edit the following script depending on your computational environment.   
```
#!/bin/csh -e

#Standard MPI + Intel compiler
set FC="mpif90"
set ftrn_prgm="Vmol01-qm-hf-oep_mpi"
set FFLAGS="-save -O3 -xHost -extend-source -mcmodel=medium -shared-intel"
set FFTW_INC="$HOME/opt/fftw/include"
set FFTW_DIR="$HOME/opt/fftw/lib"
set PFFT_INC="$HOME/opt/pfft/include"
set PFFT_DIR="$HOME/opt/pfft/lib"
set BLAS_LAPACK="-mkl=parallel"
 
 $FC -c $FFLAGS -I$PFFT_INC -I$FFTW_INC poisson_solver.f90
 $FC -c $FFLAGS $ftrn_prgm.f
 $FC -c $FFLAGS oep.f
#$FC -c $FFLAGS Vmol01-T-12-mm-mpi.f
$FC $FFLAGS $ftrn_prgm.o oep.o poisson_solver.o -L$PFFT_DIR -lpfft -L$FFTW_DIR -lfftw3 -lfftw3_mpi $BLAS_LAPACK -o $ftrn_prgm.exe
 rm -f *.o
```
# How to run Vmol
  Use the script 'mpirun.exe' in /examples/HF_OEP/.   
  ```
  #! /bin/csh
  mpirun -machinefile machines -np $argv[1]  $argv[2] 
  ```
  If you use 16 MPI processes for the parallel execution, submit the following job script
  ```
  ./mpirun.exe 16 ./Vmol01-qm-hf-oep_mpi.exe > Vmol01-qm-hf-oep_mpi.out & 
  ```
  with the machines file being 
  ```
  machine_name      cpu=16
  ```

# Output files
  The output file 'Vmol01-qm-hf-oep_mpi.out' in /examples/HF_OEP/ provides the computational settings and the energies of interest during the SCF procedure
  for the HF-OEP calculation of a water molecule. After the first SCF iteration for the RHF method, the total energy is given by 
  
  ```
    Final Energy =   -16.8838097849969
  ```

# Citation
If your work incorporates any part of 'Vmol', please cite the following references:  

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

