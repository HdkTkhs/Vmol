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

