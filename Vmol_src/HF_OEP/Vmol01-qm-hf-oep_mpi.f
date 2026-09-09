
! Copyright 2026 Hideaki Takahashi

! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C---------------------------------------------------
C
C   VMol Ver.1.0
C
C   FIRST PRINCIPLES MOLECULAR DYNAMICS
C
C   BASED ON
C
C   SD( Steepest Descent ) METHOD
C
C   EMPLOYING REAL-SPACE GRIDS
C
C  
C   REFERENCES: 
C
C      J. R. Chelikowsky, N. Troullier, and Y. Saad,
C         “Finite-difference pseudopotential method: electronic structure  
C          calculations without a basis,” 
C          Phys. Rev. Lett. 72, 1240–1243 (1994).  
C
C      H. Takahashi, T. Hori, T. Wakabayashi, and T. Nitta,
C         “Real space ab initio molecular dynamics simulations  
C          for the reactions of OH radical/OH anion with formaldehyde,” 
C          J. Phys. Chem. A 105, 4351 (2001).  
C
C      H. Takahashi, T. Hori, H. Hashimoto, and T. Nitta,
C         “A hybrid QM/MM method employing real space grids   
C          for QM water in the TIP4P water solvents,”
C          J. Comp. Chem. 22, 1252–1261 (2001).  
C
C   PROGRAMMED BY : H.TAKAHASHI ( Aug. ~1997 )      
C
C---------------------------------------------------

C---------------------------------------------------
C   EDIT LOG
C   Ver.1.0 (2001.11.27)
C   Ver.1.0 (2002.01.11) 
C     optfc:10746-10750
C   Ver.1.0 (2002.01.23) 
C     minor corrections
C   Ver.1.0 (2003.02.03) 
C     minor corrections
C     subroutine:diag
C---------------------------------------------------

      PROGRAM  VMOL

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NNUC = 36 )

      DIMENSION CRD( NNUC,3 )
      DIMENSION VLC( NNUC,3 )

      CALL MPI_INIT(IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      CALL DIMMPI

      IF(MYID .EQ. 0) THEN
        CALL OPF                             ! open files
        CALL DEF                             ! default
        CALL INIT( CRD,VLC )                 ! initiate FPMD
      ENDIF

      CALL MD( CRD,VLC )                     ! first-principles MD

      CALL CLF                               ! close files
 
      CALL MPI_FINALIZE(IERR)

      STOP 
      END

C-------------------------------
C
C     SUBROUTINE for MPI
C
C-------------------------------

      SUBROUTINE DIMMPI

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI2 / LX(0:NX*NY*NZ-1),LY(0:NX*NY*NZ-1),
     *                 LZ(0:NX*NY*NZ-1)
      COMMON / MMPI3 / NIDX(1:2),NIDY(1:2),NIDZ(1:2)
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / MMPIPOIS / ICOMM_CART

      INTEGER IDIMS(3), ICOORDS(3), ICOMM_CART
      LOGICAL PERIODS(3), REORDER

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      ! Here we use C-order because it is necessary with PFFT in Poisson solver
      IDIMS(1) = NZ
      IDIMS(2) = NY
      IDIMS(3) = NX

      PERIODS(:) = .FALSE.
      REORDER = .FALSE.         ! I really want to turn reorder, but that means the WHOLE communicator must be replaced
      CALL MPI_CART_CREATE(MPI_COMM_WORLD, 3, IDIMS,
     *     PERIODS, REORDER, ICOMM_CART, IERR)

      IF(IERR /= MPI_SUCCESS) stop "DIMMPI: MPI_CART_CREATE"
      IF(ICOMM_CART == MPI_COMM_NULL) THEN
         stop "Requested # of processes < the world size"
      ENDIF
      CALL MPI_CART_GET(ICOMM_CART, 3,
     *     IDIMS, PERIODS, ICOORDS, IERR)
      IF(IERR /= MPI_SUCCESS) stop "DIMMPI: MPI_CART_GET"
      MZ = ICOORDS(1)
      MY = ICOORDS(2)
      MX = ICOORDS(3)

      DO ID = 0, NUMPROCS - 1
         CALL MPI_CART_COORDS(ICOMM_CART, ID, 3, ICOORDS, IERR)
         LZ(ID) = ICOORDS(1)
         LY(ID) = ICOORDS(2)
         LX(ID) = ICOORDS(3)
      ENDDO
      
      CALL MPI_CART_SHIFT(ICOMM_CART, 0, 1, NIDZ(1), NIDZ(2), IERR)
      WHERE(NIDZ == MPI_PROC_NULL) NIDZ = -1
      CALL MPI_CART_SHIFT(ICOMM_CART, 1, 1, NIDY(1), NIDY(2), IERR)
      WHERE(NIDY == MPI_PROC_NULL) NIDY = -1
      CALL MPI_CART_SHIFT(ICOMM_CART, 2, 1, NIDX(1), NIDX(2), IERR)
      WHERE(NIDX == MPI_PROC_NULL) NIDX = -1

      JX = NMAXX/NX
      JY = NMAXY/NY
      JZ = NMAXZ/NZ

      RETURN
      END

C--------------------------------------------
      SUBROUTINE DEF
C--------------------------------------------
      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'                         ! Vmol

C     PARAMETER ( NNUC   = 36 )
C     PARAMETER ( NLINK1 = 10 )                 ! Vmol

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          FREEZE*5,PRINT*5,DGF*3

      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
!     COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1), ! Vmol
!    *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PRMT3 / blch(NLINK1),blcc(NLINK1),  ! Vmol
     *                 bkch(NLINK1),bkcc(NLINK1)   ! Vmol
      COMMON / NRDF1 / NMRDF,NRDF(NNUC)

      COE    = 60.D0
      NDEN   = 5
      DTMD   = 41.3411054611601
      MDMAX  = 50000
      TEMP   = 298.15D0
      NRVLC  = 50
      NCHK   = 50
      CONV   = 1.0D-5
      NRST   = 'NEW'
      NOPT   = 'ONE'
      EXC    = 'UBLYP'
      NQMMM  = 'QM'
      PRINT  = 'LARGE'
      FREEZE = 'FALSE'
      DGF    = 'DG4'
      NMRDF  = 0
      NMM2   = 0                ! Vmol
      NLINK  = 0                ! Vmol

      DO I = 1, NLINK1          ! Vmol
        blch(I) = 0.D0          ! Vmol
        blcc(I) = 0.D0          ! Vmol
        bkch(I) = 0.D0          ! Vmol
        bkcc(I) = 0.D0          ! Vmol
      ENDDO                     ! Vmol

      RETURN
      END

C--------------------------------------
C
C   SUBROUTINE INITIATE AB INITIO MD       
C
C--------------------------------------

      SUBROUTINE INIT( CRD,VLC ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"                          ! mpi
      include 'mpi.i'                           ! Vmol
      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'                         ! Vmol
      include 'nlocd.i'                         ! Vmol   non-local d
      include 'pc_crr.i'                        ! Vmol   pcc

C     PARAMETER ( NNUC   = 36 )
C     PARAMETER ( NLINK1 = 10 )                 ! Vmol
C     PARAMETER ( NSL    = 3000  )
   
      DIMENSION CRD( NNUC,3 )
      DIMENSION VLC( NNUC,3 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1), ! Vmol
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PRMT3 / blch(NLINK1),blcc(NLINK1),  ! Vmol
     *                 bkch(NLINK1),bkcc(NLINK1)   ! Vmol
      COMMON / NRDF1 / NMRDF,NRDF(NNUC)
      COMMON / LJQMMM / SQMMM(NNUC,maxatm),EPQMMM(NNUC,maxatm),
     *                  chgmm(maxatm),chgqm(NNUC)
      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,                 ! Vmol
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,                 ! Vmol
C    *                  WNUC(NFUZZY,NMAX,NMAX,NMAX),RSIZE(NNUC),           ! Vmol
C    *                  WNUC(NFUZZY,NMAX/NX,NMAX/NY,NMAX/NZ),
     *                  RSIZE(NNUC),
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF                       ! Vmol

      NAMELIST/INIDAT/COE,NDEN,DTMD,MDMAX,TEMP,
     *                NRST,NOPT,EXC,NQMMM,NRVLC,
     *                NCHK,CONV,FREEZE,PRINT,
C    *                DGF,NMRDF                  
     *                DGF,NMRDF,NMM2,NLINK                    ! Vmol

C---- read qm.dat file ----

      READ(20,INIDAT)
      READ(20,*)

      IF(NQMMM .EQ. 'LINK') THEN                                     ! Vmol
        nlqm = NNUC-NLINK                                            ! Vmol
        DO I =1,nlqm                                                 ! Vmol
          READ(20,*) MMID(I),NPD(I),NPC(I),ZA(I),ZVAL(I),SIG(I),EPSQM(I)    ! Vmol     non-local d
C         RSIZE(I) = SIG(I)                                          ! activated 2016.05.17 takahashi
C------- Temporary treatment for Miki LJ parameter -----------------------
C        SIG(I)   = SIG(I)*2.D0/0.529177D0
C        EPSQM(I) = -1.D0*EPSQM(I)/627.51D0
C------- 2013.0801 takahashi ---------------------------------------------
        ENDDO                                                 ! Vmol

        READ(20,*)                                            ! Vmol
        NCOUNT = 0                                            ! Vmol
        DO I =nlqm+1,NNUC                                     ! Vmol
          NCOUNT = NCOUNT + 1                                 ! Vmol
          IF( NCOUNT .GT. NLINK1 ) THEN                       ! Vmol
            WRITE(*,*)                                        ! Vmol
     *      ' init.f : size of dimension too small (nlink1)'  ! Vmol
            STOP                                              ! Vmol
          ENDIF                                               ! Vmol
C         READ(20,*) ZA(I),ZVAL(I),SIG(I),EPSQM(I),           ! Vmol
C    *               NLAQM(NCOUNT),NLAMM(NCOUNT),             ! Vmol
C    *               blch(NCOUNT),blcc(NCOUNT),               ! Vmol
C    *               bkch(NCOUNT),bkcc(NCOUNT)                ! Vmol
          ZA(I)    = 1.0D0                                    ! Vmol
          ZVAL(I)  = 1.0D0                                    ! Vmol
          SIG(I)   = 0.0D0                                    ! Vmol
          EPSQM(I) = 0.0D0                                    ! Vmol
C         RSIZE(I) = 5.6200D0                                 ! Vmol (AMBER95,vdw 34, see amber-b.prm in tinker)
C         RSIZE(I) = 3.35482D0                                ! Vmol  activated 2016.05.17 takahashi
          READ(20,*) NLAQM(NCOUNT),NLAMM(NCOUNT),             ! Vmol
     *               blch(NCOUNT),blcc(NCOUNT),               ! Vmol
     *               bkch(NCOUNT),bkcc(NCOUNT)                ! Vmol
        ENDDO                                                 ! Vmol
      ELSEIF(NQMMM .EQ. 'MM') THEN                            ! Vmol
        nlqm = NNUC-NLINK                                     ! Vmol
        DO I =1,nlqm                                          ! Vmol
          READ(20,*) MMID(I),NPD(I),NPC(I),ZA(I),ZVAL(I),SIG(I),EPSQM(I),     ! Vmol     non-local d
     *               chgqm(I)                                 ! Vmol
C         RSIZE(I) = SIG(I)                                   ! Vmol
        ENDDO                                                 ! Vmol
        READ(20,*)                                            ! Vmol
        NCOUNT = 0                                            ! Vmol
        DO I =nlqm+1,NNUC                                     ! Vmol
          NCOUNT = NCOUNT + 1                                 ! Vmol
          IF( NCOUNT .GT. NLINK1 ) THEN                       ! Vmol
            WRITE(*,*)                                        ! Vmol
     *      ' init.f : size of dimension too small (nlink1)'  ! Vmol
            STOP                                              ! Vmol
          ENDIF                                               ! Vmol
          ZA(I)    = 1.0D0                                    ! Vmol
          ZVAL(I)  = 1.0D0                                    ! Vmol
          SIG(I)   = 0.0D0                                    ! Vmol
          EPSQM(I) = 0.0D0                                    ! Vmol
          READ(20,*) NLAQM(NCOUNT),NLAMM(NCOUNT),             ! Vmol
     *               blch(NCOUNT),blcc(NCOUNT),               ! Vmol
     *               bkch(NCOUNT),bkcc(NCOUNT)                ! Vmol
        ENDDO                                                 ! Vmol
      ELSE                                                    ! Vmol
C     WRITE(*,*) ' NLINK = ', NLINK
        DO I =1,NNUC
C         READ(20,*) ZA(I),ZVAL(I),SIG(I),EPS(I),
          READ(20,*) MMID(I),NPD(I),NPC(I),ZA(I),ZVAL(I),SIG(I),EPSQM(I),    ! Vmol     non-local d
     *               CRD(I,1),CRD(I,2),CRD(I,3)
C         RSIZE(I) = SIG(I)                                   ! Vmol
        ENDDO
      ENDIF                                                   ! Vmol

      DO I =1,NNUC        ! converted to au for qmbuf 
      DO J =1,3   
!        CRD(I,J) = CRD(I,J)/0.529177D0        ! input data are given in au.
      ENDDO 
      ENDDO 

      DO I =1,NNUC                                            ! Vmol Atomic radii in angstroms
        NZA = INT(ZA(I))                                      ! Vmol (JCP, 41, 3199 (1964))
        IF(NZA    .eq.1) then                                 ! Vmol
C         RSIZE(I) = 0.25D0                                   ! Vmol
          RSIZE(I) = 0.35D0                                   ! Vmol
        ELSEIF(NZA.eq.6) then                                 ! Vmol
          RSIZE(I) = 0.7D0                                    ! Vmol
        ELSEIF(NZA.eq.7) then                                 ! Vmol
          RSIZE(I) = 0.65D0                                   ! Vmol
        ELSEIF(NZA.eq.8) then                                 ! Vmol
          RSIZE(I) = 0.6D0                                    ! Vmol
        ELSEIF(NZA.eq.9) then                                 ! Vmol F added 20230707
          RSIZE(I) = 0.5D0                                    ! Vmol
        ELSEIF(NZA.eq.17) then                                ! Vmol
          RSIZE(I) = 1.0D0                                    ! Vmol
        ELSEIF(NZA.eq.20) then                                ! Vmol
          RSIZE(I) = 1.8D0                                    ! Vmol
        ELSEIF(NZA.eq.25) then                                ! Vmol
          RSIZE(I) = 1.4D0                                    ! Vmol
        ENDIF                                                 ! Vmol
      ENDDO                                                   ! Vmol
                                                              ! Vmol
      IF(NMRDF.NE.0) THEN
        READ(20,*)
        READ(20,*) ( NRDF(I),I=1,NMRDF )
      ENDIF

C-----------------------------------------------------------------------------------------------------------
C       $INIDAT
C       COE       Cut off energy / au
C       NDEN      Number of dense grids
C       DTMD      MD time step / au     
C       MDMAX     Number of MD iterations 
C       TEMP      Temperature / K
C       NRVLC     Step number of rescale velocity
C       NCHK      Step number of check point file
C       CONV      Convergence of Density Matrix
C       NRST      Restart flag( NEW,RESTART )
C       NOPT      Option flag( ONE, OPT, MD )
C       EXC       SCF flag1(  UBLYP, RBLYP, UPZ, RPZ, UXalpha, RXalpha  )
C       PRINT     SCF write flag ( LARGE(detail) ,SMALL(rough) )
C       NQMMM     SCF flag2(  QM, QMMM, LINK, MM, QMVD, QMVC  )                     ! Vmol
C       FREEZE    QM molecules freeze flag 
C       DGF       Double grid flag( DG, DG4, DG6 )
C       NMRDF     Number of molecules to calculate radial distribution function  
C       NMM2      Number of boundary QM atoms                                       ! Vmol
C       NLINK     Number of LINK atoms                                              ! Vmol
C       $END
C------------------------------------------------------------------------------------------------------------

       DO I=1,NNUC                               ! gaussian distribution( not restart )                    
       VLC(I,1) = GAUSS(DUMMY)
       VLC(I,2) = GAUSS(DUMMY)
       VLC(I,3) = GAUSS(DUMMY)
       ENDDO
       CALL REVLC( VLC )

      IF(NQMMM .NE. 'MM') THEN                                ! Vmol
C     OPEN(21,FILE='orb.dat',STATUS='OLD')       ! open g94 output
      ENDIF                                                   ! Vmol

C     IF(NQMMM .EQ. 'QMMM') THEN
      IF(NQMMM .EQ. 'QMMM' .OR. NQMMM .EQ. 'LINK' .OR.        ! Vmol
     *   NQMMM .EQ. 'MM' ) THEN                               ! Vmol
        CALL MMOPF
      ELSEIF(NQMMM.EQ.'QMVD') THEN                            ! Vmol
        OPEN(UNIT=220,file='VPCE-dense.dat',status='old')     ! Vmol
      ELSEIF(NQMMM.EQ.'QMVC') THEN                            ! Vmol
        OPEN(UNIT=230,file='VPCE-coarse.dat',status='old')    ! Vmol
      ENDIF

      WRITE(*,*)                             
      IF( NRST .EQ. 'NEW' ) THEN
      WRITE(*,*) '--------- Start ab initio Molecular Dynamics 
     *Simulation ---------'
      ELSE
      WRITE(*,*) '------- Continue ab initio Molecular Dynamics 
     *Simulation --------'
      ENDIF

      WRITE(*,*)                             
      WRITE(*,*) '------------------------------ Method -------------
     *--------------' 
      WRITE(*,*)                             
      IF( NOPT .EQ. 'ONE' ) THEN
      WRITE(*,*) 'One Point Calculation'
      ELSEIF( NOPT .EQ. 'OPT' ) THEN
      WRITE(*,*) 'Geometry Optimization'
      ELSEIF( NOPT .EQ. 'MD' ) THEN
      WRITE(*,*) 'Molecular Dynamics'
      ENDIF
      IF( EXC .EQ. 'RPZ' .OR. EXC .EQ. 'RBLYP' .OR. EXC.EQ.'RHF' 
     *                   .OR. EXC .EQ. 'RXalpha' ) THEN
      WRITE(*,*) 'Restricted DFT'
      ELSEIF( EXC .EQ. 'UPZ' .OR. EXC .EQ. 'UBLYP' .OR. EXC.EQ.'UHF' 
     *                       .OR. EXC .EQ. 'UXalpha' ) THEN
      WRITE(*,*) 'Unrestricted DFT'
      ENDIF
      WRITE(*,*) 'Exchange Correlation Functional = ',EXC                             
      WRITE(*,*) 'System                          = ',NQMMM                             
      IF(FREEZE.EQ.'TRUE') THEN 
      WRITE(*,*) 'QM molecules are fixed'                             
      ENDIF
      IF(NMRDF.NE.0) THEN
      WRITE(*,*) ' RDFs are considered for the atoms; ',
     *           (NRDF(I),I=1,NMRDF)                             
      ENDIF

      IF( ( NQMMM .EQ. 'QM' ) .AND. ( FREEZE .EQ. 'TRUE' ) ) THEN
       WRITE(*,*)                             
       WRITE(*,*) 'Warning!:Request is identical to
     * "QM and One point" calculation.'
       WRITE(*,*) 'Switched to "QM and One point" calculation.'
       NOPT = 'ONE'
      ENDIF
      
      WRITE(*,*)                             
      WRITE(*,*) '------------------- Computational Parameters 
     *for QM -------------'
      WRITE(*,*)                             
      WRITE(*,93) ' Cut Off Energy [a.u.]           = ',COE
      WRITE(*,94) ' Number of Dense Grids           = ',NDEN
      WRITE(*,95) ' MD Time Step   [a.u.]           = ',DTMD
      WRITE(*,96) ' MD Time Step   [fs]             = ',DTMD*2.4189D-02
      WRITE(*,97) ' Temperature    [K]              = ',TEMP
      WRITE(*,98) ' Step number of vlc rescaling    = ',NRVLC
      WRITE(*,98) ' Step number of check point file = ',NCHK
      WRITE(*,99) ' Convergence of Density Matrix   = ',CONV
      WRITE(*,*)                             

      IF(NQMMM .NE. 'LINK') THEN                              ! Vmol

        WRITE(*,*) 'Initial Coordinates [a.u.]'
        DO I=1,NNUC                                ! write coordinates
        WRITE(*,100) ZA(I),ZVAL(I),
     *      CRD(I,1),CRD(I,2),CRD(I,3)
        ENDDO
       
        IF( NOPT .EQ. 'MD' ) THEN
          WRITE(*,*)                             
          WRITE(*,*) 'Initial Velocities'
          DO I=1,NNUC                                ! write velocities
          WRITE(*,100) ZA(I),ZVAL(I), 
     *                 VLC(I,1),VLC(I,2),VLC(I,3)
          ENDDO
        ENDIF

      ENDIF                                                   ! Vmol

      IF(NQMMM .NE. 'MM') THEN                                ! Vmol
        WRITE(*,*)
        WRITE(*,*)'Lennard Jones Parameters'
        WRITE(*,*)'         SIG[a.u.]','   EPS[a.u.]'
        DO I = 1,NNUC
        IF(ZA(I).EQ.8) THEN
        WRITE(*,200) I,SIG(I),EPSQM(I)                        ! Vmol
        ELSEIF(ZA(I).EQ.6) THEN
        WRITE(*,201) I,SIG(I),EPSQM(I)                        ! Vmol
        ELSEIF(ZA(I).EQ.1) THEN
        WRITE(*,202) I,SIG(I),EPSQM(I)                        ! Vmol
        ELSEIF(ZA(I).EQ.15) THEN
        WRITE(*,205) I,SIG(I),EPSQM(I)                        ! Vmol
        ELSEIF(ZA(I).EQ.17) THEN
        WRITE(*,203) I,SIG(I),EPSQM(I)                        ! Vmol
        ELSEIF(ZA(I).EQ.7) THEN                               ! Vmol
        WRITE(*,204) I,SIG(I),EPSQM(I)                        ! Vmol
        ELSEIF(ZA(I).EQ.20) THEN                              ! Vmol
        WRITE(*,206) I,SIG(I),EPSQM(I)                        ! Vmol
        ELSEIF(ZA(I).EQ.25) THEN                              ! Vmol
        WRITE(*,207) I,SIG(I),EPSQM(I)                        ! Vmol
        ELSEIF(ZA(I).EQ.9) THEN                               ! Vmol
        WRITE(*,208) I,SIG(I),EPSQM(I)                        ! Vmol
        ENDIF
        ENDDO
      ELSE                                                    ! Vmol
        WRITE(*,*)                                            ! Vmol
        WRITE(*,*)'Lennard Jones Parameters and Charges'      ! Vmol
        WRITE(*,*)'         SIG[a.u.]','   EPS[a.u.]',        ! Vmol
     *            '   chgqm[e]'                               ! Vmol
        tchgqm = 0.d0                                         ! Vmol
        DO I = 1,NNUC                                         ! Vmol
        tchgqm = tchgqm + chgqm(I)
        IF(ZA(I).EQ.8) THEN                                   ! Vmol
        WRITE(*,300) I,SIG(I),EPSQM(I),chgqm(I)               ! Vmol
        ELSEIF(ZA(I).EQ.6) THEN                               ! Vmol
        WRITE(*,301) I,SIG(I),EPSQM(I),chgqm(I)               ! Vmol
        ELSEIF(ZA(I).EQ.1) THEN                               ! Vmol
        WRITE(*,302) I,SIG(I),EPSQM(I),chgqm(I)               ! Vmol
        ELSEIF(ZA(I).EQ.15) THEN                              ! Vmol
        WRITE(*,305) I,SIG(I),EPSQM(I),chgqm(I)               ! Vmol
        ELSEIF(ZA(I).EQ.17) THEN                              ! Vmol
        WRITE(*,303) I,SIG(I),EPSQM(I),chgqm(I)               ! Vmol
        ELSEIF(ZA(I).EQ.7) THEN                               ! Vmol
        WRITE(*,304) I,SIG(I),EPSQM(I),chgqm(I)               ! Vmol
        ELSEIF(ZA(I).EQ.20) THEN                              ! Vmol
        WRITE(*,306) I,SIG(I),EPSQM(I),chgqm(I)               ! Vmol
        ELSEIF(ZA(I).EQ.25) THEN                              ! Vmol
        WRITE(*,307) I,SIG(I),EPSQM(I),chgqm(I)               ! Vmol
        ENDIF                                                 ! Vmol
        ENDDO                                                 ! Vmol
        WRITE(*,*)                                            ! Vmol
        WRITE(*,*)'Total charge of the QM subsystem = ',tchgqm! Vmol
      ENDIF                                                   ! Vmol

        WRITE(*,*)

93    FORMAT( A35,3X,F5.1)
94    FORMAT( A35,3X,I1)
95    FORMAT( A35,3X,F16.13)
96    FORMAT( A35,3X,F16.13)
97    FORMAT( A35,3X,F6.2)
98    FORMAT( A35,3X,I5)
99    FORMAT( A35,3X,E8.1)
100   FORMAT( F6.1,F6.1,3F12.6 )
200   FORMAT( '  O',I3,F12.6,F12.6 )
201   FORMAT( '  C',I3,F12.6,F12.6 )
202   FORMAT( '  H',I3,F12.6,F12.6 )
203   FORMAT( ' Cl',I3,F12.6,F12.6 )
204   FORMAT( '  N',I3,F12.6,F12.6 )
205   FORMAT( '  P',I3,F12.6,F12.6 )
206   FORMAT( ' Ca',I3,F12.6,F12.6 )
207   FORMAT( ' Mn',I3,F12.6,F12.6 )
208   FORMAT( '  F',I3,F12.6,F12.6 )
300   FORMAT( '  O',I3,3F12.6 )
301   FORMAT( '  C',I3,3F12.6 )
302   FORMAT( '  H',I3,3F12.6 )
303   FORMAT( ' Cl',I3,3F12.6 )
304   FORMAT( '  N',I3,3F12.6 )
305   FORMAT( '  P',I3,3F12.6 )
306   FORMAT( ' Ca',I3,3F12.6 )
307   FORMAT( ' Mn',I3,3F12.6 )

      RETURN
      END

C---------------------------------------------------
C
C     SUBROUTINE STEEPEST DESCENT
C
C---------------------------------------------------

      SUBROUTINE MD( CRD,VLC )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"                          ! mpi
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'                         ! Vmol

C     PARAMETER ( NMAX =  80 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NDIM =   3 )
C     PARAMETER ( NORA =  49 )
      PARAMETER ( PI   = 3.14159265358979323D0 )
C     PARAMETER ( EPSOPT= 5.0D-5 )
C     PARAMETER ( EPSMD = 5.0D-3 )
C     PARAMETER ( NLINK1 = 10 )                 ! Vmol

      PARAMETER ( NSP = 255 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / GRID2 / DX, DY, DZ
      COMMON / ACELL / XL, YL, ZL
      COMMON / ZMASS / ZCM(100)
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),     ! Vmol
     *                 NMMSW(maxatm),MMID(NNUC)
C     COMMON / MMDAT / SCRD( NDIM,4*NSP ),SCRD1( NDIM,4*NSP ),
C    *                 FRCS( 4*NSP,NDIM )
      COMMON / MMDAT / SCRD( NDIM,maxatm ),SCRD1( NDIM,maxatm ),   ! Vmol
     *                 FRCS( maxatm,NDIM )                         ! Vmol
C     COMMON / OPTDAT / EHOMO,ELUMO,DIPQ1(4)
      COMMON / NRDF1 / NMRDF,NRDF(NNUC)

      DIMENSION WFA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION WFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
C     DIMENSION WFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,1    )
      DIMENSION WFNA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION WFNB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION CRD(  NNUC,NDIM )
      DIMENSION CRDOLD(  3,NDIM )
      DIMENSION VLC(  NNUC,NDIM )
      DIMENSION FRC1( NNUC,NDIM )
      DIMENSION FRC2( NNUC,NDIM )

      DIMENSION NMID1(NNUC)                                    ! Vmol

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      NDAT = NNUC*  NDIM
      MDAT = maxatm*NDIM

C---- COMPUTE PARAMETERS ----

      CALL CONFIG

C---- INITIAL WAVEFUNC. AND PSEUDOPOTENTIAL DATA ----

      KMD= 0

      IF( NQMMM .NE. 'MM' ) THEN                               ! Vmol
C     CALL TRWF(  WFA,WFB,WFNA,WFNB )                                    
      CALL TRWF3( CRD,WFA,WFB,WFNA,WFNB )                      ! 2015.03.09 takahashi  
      IF(MYID.EQ.0) THEN
        CALL PSDAT                                             ! read pseudopotential data
      ENDIF
      ENDIF                                                    ! Vmol

      CALL PUTPS                                               ! takahashi (new subroutine)
C     CALL PCC                                                 ! read and distribute PCC data

      IF(NQMMM .EQ. 'QM' .OR. NQMMM .EQ. 'QMMM' ) THEN         ! Vmol
      IF(MYID.EQ.0) THEN
        CALL DISMAT( CRD )
      ENDIF
      ENDIF                                                    ! Vmol

C---- READ COORDINATES OF MM SITES ----

      IF(NQMMM .EQ. 'QMMM') THEN
        IF(MYID.EQ.0) THEN
!         CALL START( SCRD,KMD    )
!         CALL TRANS( SCRD,FRCS,1 )
        ENDIF

      ELSEIF(NQMMM .EQ. 'LINK') THEN                       ! Vmol

C     IF(MYID.EQ.0) THEN                                   ! Vmol
!         call dynamic(1,KMD)                              ! Vmol jflag=( 1,2:initialize, 3,4:dynamics )
C     ENDIF                                                ! Vmol
!         call spqmmm(NMID1)                               ! Vmol
!         call qmmmlj                                      ! Vmol
!         call trntin( scrd,crd )                          ! Vmol
                                                           ! Vmol
      IF(MYID.EQ.0) THEN                                   ! Vmol
          WRITE(*,*)                                       ! Vmol
          WRITE(*,*) 'Initial Coordinates [a.u.]'          ! Vmol
          DO I=1,NNUC-NLINK                                ! Vmol
            WRITE(*,98) NMID1(I),'=>',MMID(I),ZA(I),       ! Vmol
     *           ZVAL(I),CRD(I,1),CRD(I,2),CRD(I,3)        ! Vmol
          ENDDO                                            ! Vmol
      ENDIF                                                ! Vmol 
          CALL LINKA_FIX( CRD )                            ! Vmol adjust link atom position
C         CALL LINKA( CRD )                                ! Vmol adjust link atom position
          IF(MYID.EQ.0) CALL DISMAT( CRD )                 ! Vmol
                                                           ! Vmol
      ELSEIF(NQMMM .EQ. 'MM') THEN                         ! Vmol
                                                           ! Vmol
!         call dynamic(1,KMD)                              ! Vmol jflag=( 1,2:initialize, 3,4:dynamics )
!         call spqmmm(NMID1)                               ! Vmol
!         call qmmmlj                                      ! Vmol
!         call trntin( scrd,crd )                          ! Vmol
                                                           ! Vmol
          WRITE(*,*)                                       ! Vmol
          WRITE(*,*) 'Initial Coordinates [a.u.]'          ! Vmol
          DO I=1,NNUC-NLINK                                ! Vmol
            WRITE(*,98) NMID1(I),'=>',MMID(I),ZA(I),       ! Vmol
     *           ZVAL(I),CRD(I,1),CRD(I,2),CRD(I,3)        ! Vmol
          ENDDO                                            ! Vmol
                                                           ! Vmol
          CALL LINKA( CRD )                                ! Vmol adjust link atom position
          CALL DISMAT( CRD )                               ! Vmol
                                                           ! Vmol
      ENDIF

C----- construction of initial guess with LCAO coeff. -------

      IF( NQMMM .NE. 'MM' ) THEN                           ! Vmol
C       CALL TRWF3( CRD,WFA,WFB,WFNA,WFNB )                ! 2015.03.09 takahashi  
      ENDIF 

C------------------------------------------------------------

C       CALL MPI_BCAST( CRD(1,1), NDAT,MPI_DOUBLE_PRECISION,0,
C    *                                 MPI_COMM_WORLD,IERR )
C       CALL MPI_BCAST( SCRD(1,1),MDAT,MPI_DOUBLE_PRECISION,0,
C    *                                 MPI_COMM_WORLD,IERR )

C-----read initial Hartree potential --------              ! Vmol
C     CALL IVCOU( CRD )                                    ! Vmol
C--------------------------------------------              ! Vmol
                                            
C---- ONE POINT CALCULATION ----

      IF( NOPT .EQ. 'ONE' ) THEN
    
        CALL MPI_BCAST( CRD(1,1), NDAT,MPI_DOUBLE_PRECISION,0,    ! 2005.07.23 takahashi
     *                                 MPI_COMM_WORLD,IERR )
        CALL MPI_BCAST( SCRD(1,1),MDAT,MPI_DOUBLE_PRECISION,0,    ! 2005.07.23 takahashi
     *                                 MPI_COMM_WORLD,IERR )

        CALL GRND( WFA,WFB,WFNA,WFNB,CRD,FRC1,PENG,0 )     ! Hellmann-Feynman Force

      IF(MYID.EQ.0) THEN
!       CALL WCRD(CRD,SCRD,VLC,FRC1,0)
      ENDIF

        IF( NQMMM.EQ.'QMVD' .OR. NQMMM.EQ.'QMVC' ) THEN    ! Vmol
          CALL SSINT( WFA,WFB,IMD,SSE )                    ! Vmol
        ENDIF                                              ! Vmol


      RETURN  

C---- GEOMETRY OPTIMIZATION -----

      ELSEIF( NOPT .EQ. 'OPT' ) THEN
      
      IF(MYID.EQ.0) THEN
!       CALL WCRD(CRD,SCRD,VLC,FRC1,0)
      ENDIF

      CALL GRND( WFA,WFB,WFNA,WFNB,CRD,FRC1,PENG,0 )       ! Hellmann-Feynman Force

      KMD = KMD + 1

      DO 30 IMD=KMD,MDMAX                                  ! optimization loop

      IF(MYID.EQ.0) THEN

      PFE = PENG
      NEPS= 0
      
      DO NA=1,NNUC
      NZA = INT(ZA(NA))

        IF( DABS(FRC1(NA,1)) .GT. EPSMD ) THEN
         CRD(NA,1)= CRD(NA,1)
     *            +(FRC1(NA,1))*DTMD**2.D0/(2.D0*ZCM(NZA))
        ELSE
         NEPS = NEPS + 1
        ENDIF

        IF( DABS(FRC1(NA,2)) .GT. EPSMD ) THEN
         CRD(NA,2)= CRD(NA,2)
     *            +(FRC1(NA,2))*DTMD**2.D0/(2.D0*ZCM(NZA))
        ELSE
         NEPS = NEPS + 1
        ENDIF

        IF( DABS(FRC1(NA,3)) .GT. EPSMD ) THEN
         CRD(NA,3)= CRD(NA,3)
     *            +(FRC1(NA,3))*DTMD**2.D0/(2.D0*ZCM(NZA))
        ELSE
         NEPS = NEPS + 1
        ENDIF

      ENDDO

      WRITE(*,*)                             
      WRITE(*,*) '    QM Coordinates'
      DO J=1,NNUC                                ! write coordinates
        WRITE(*,99) ZA(J),ZVAL(J),
     *                CRD(J,1),CRD(J,2),CRD(J,3)
      ENDDO

      IF(NQMMM.EQ.'QMMM') THEN
       WRITE(*,*)
       WRITE(*,*) '!!!!!!!! START CLASSICAL MOLECULAR DYNAMICS !!!!!!!!'
       WRITE(*,*)
!      CALL TRANS(  SCRD,FRCS,2 )
!      CALL MDMPOL( SCRD,FRCS,0 )
!      CALL VALUES
!      CALL TRANS(  SCRD,FRCS,1 )
      ENDIF

      ENDIF

      CALL MPI_BCAST( CRD(1,1), NDAT,MPI_DOUBLE_PRECISION,0,
     *                               MPI_COMM_WORLD,IERR )
      CALL MPI_BCAST( SCRD(1,1),MDAT,MPI_DOUBLE_PRECISION,0,
     *                               MPI_COMM_WORLD,IERR )
      CALL MPI_BCAST( FRCS(1,1),MDAT,MPI_DOUBLE_PRECISION,0,
     *                               MPI_COMM_WORLD,IERR )

      CALL GRND( WFA,WFB,WFNA,WFNB,CRD,FRC1,PENG,IMD )       ! Hellmann-Feynman Force

      IF(MYID.EQ.0) THEN

      IF(NQMMM.EQ.'QMMM') THEN
      CALL CS( SCRD )
!     CALL WCRD(CRD,SCRD,VLC,FRC1,IMD)
       IF(MOD(IMD,NCHK).EQ.0) THEN
        CALL AGAIN(2,IMD)
       ENDIF
      ENDIF
      
      DEMD = PENG - PFE
      PFE  = PENG
        WRITE(*,*)                             
        WRITE(*,*) '  MD ITERATION    = ',IMD
        WRITE(*,*) '  DEMD            = ',DEMD

      IF( ( DABS(DEMD).LT.EPSOPT ) .AND.
     *    ( NEPS .EQ. 3*NNUC  ) .AND.
     *    ( DEMD .LT. 0.      ) ) THEN
        WRITE(*,*)
        WRITE(*,*) '!!!!!Optimization Converged!!!!!'
        RETURN
      ENDIF

      ENDIF

30    CONTINUE

C---- MOLECULAR DYNAMICS -----

      ELSEIF( NOPT .EQ. 'MD' ) THEN
      
      IF(MYID.EQ.0) THEN
C     CALL WCRD(CRD,SCRD,VLC,FRC1,0)
      ENDIF

      CALL GRND( WFA,WFB,WFNA,WFNB,CRD,FRC1,PENG,0 )       ! Hellmann-Feynman Force
     
      IF(NQMMM .EQ. 'LINK' .OR. NQMMM .EQ. 'MM') THEN      ! Vmol
!         IF(MYID.EQ.0) call trntin1( FRCS,FRC1 )          ! Vmol
!         call dynamic(2,KMD)                              ! Vmol jflag=( 1,2:initialize, 3,4:dynamics )
C         CALL ELECF( KMD )                                ! Vmol
          MDLOOP = 0                                       ! Vmol
      ELSE                                                 ! Vmol
          MDLOOP = NNUC                                    ! Vmol
      ENDIF                                                ! Vmol

      KMD = KMD+1

      CALL MPI_BCAST( MDLOOP,1,MPI_INTEGER,0,
     *                         MPI_COMM_WORLD,IERR )

      DO 10 IMD=KMD,MDMAX

      IF(MYID.EQ.0) THEN

      IF(FREEZE.EQ.'FALSE') THEN

C     DO NA=1,NNUC
      DO NA=1,MDLOOP                                       ! Vmol
      NZA = INT(ZA(NA))

      CRD(NA,1) = CRD(NA,1) + VLC(NA,1)*DTMD + FRC1(NA,1)            
     *          *DTMD**2.D0/(2.D0*ZCM(NZA))

      CRD(NA,2) = CRD(NA,2) + VLC(NA,2)*DTMD + FRC1(NA,2)            
     *          *DTMD**2.D0/(2.D0*ZCM(NZA))

      CRD(NA,3) = CRD(NA,3) + VLC(NA,3)*DTMD + FRC1(NA,3)            
     *          *DTMD**2.D0/(2.D0*ZCM(NZA))

      ENDDO

      ELSE

C     DO NA=1,NNUC
      DO NA=1,MDLOOP                                       ! Vmol
      NZA = INT(ZA(NA))

      CRD(NA,1) = CRD(NA,1)

      CRD(NA,2) = CRD(NA,2)

      CRD(NA,3) = CRD(NA,3)

      ENDDO

      ENDIF

      IF(NQMMM .NE. 'MM') THEN                             ! Vmol
      WRITE(*,*)
      WRITE(*,*) 'MD TIME( fs ) = ',DBLE(IMD)*DTMD*2.4189D-02
      ENDIF                                                ! Vmol

      ENDIF                                                ! IF(MYID.EQ.0) THEN                

      IF(NQMMM.EQ.'QMMM') THEN

        IF(MYID.EQ.0) THEN

        WRITE(*,*)
        WRITE(*,*) '!!!!!!! START CLASSICAL MOLECULAR DYNAMICS !!!!!!!'
        WRITE(*,*)

!       CALL TRANS( SCRD,FRCS,2 )
!       CALL MDMPOL( SCRD,FRCS,1 )
!       CALL VALUES
!       CALL TRANS( SCRD,FRCS,1 )

        ENDIF                                              ! IF(MYID.EQ.0) THEN                

      ELSEIF(NQMMM .EQ. 'LINK') THEN                       ! Vmol

!       call dynamic(3,IMD)                                ! Vmol jflag=( 1,2:initialize, 3,4:dynamics )
!       call trntin( scrd,crd )                            ! Vmol
!       CALL LINKA_FIX( CRD )                              ! Vmol adjust link atom position
C       CALL LINKA( CRD )                                  ! Vmol adjust link atom position

      ELSEIF(NQMMM .EQ. 'MM') THEN                         ! Vmol

!       call dynamic(3,IMD)                                ! Vmol jflag=( 1,2:initialize, 3,4:dynamics )
!       call trntin( scrd,crd )                            ! Vmol
        CALL LINKA( CRD )                                  ! Vmol adjust link atom position

      ENDIF

      CALL GRND( WFA,WFB,WFNA,WFNB,CRD,FRC2,PENG,IMD )     ! Hellmann-Feynman Force

C     CALL WCRD(CRD,SCRD,VLC,FRC1,IMD)

      IF(NQMMM .EQ. 'LINK' .OR. NQMMM .EQ. 'MM') THEN      ! Vmol
C       CALL SSINT( WFA,WFB,IMD,SSE )                      ! Vmol
C       CALL ELECF( IMD )                                  ! Vmol
!       IF(MYID.EQ.0) call trntin1( FRCS,FRC2 )            ! Vmol
!       call dynamic(4,IMD)                                ! Vmol jflag=( 1,2:initialize, 3,4:dynamics )
C       IF(NMRDF.NE.0) THEN                                ! Vmol
C       CALL RDFL(CRD,SCRD,IMD)                            ! Vmol radial distribution function
C       ENDIF                                              ! Vmol
      ENDIF                                                ! Vmol

      IF(MYID.EQ.0) THEN

!     CALL WCRD(CRD,SCRD,VLC,FRC1,IMD)

      IF(FREEZE.EQ.'FALSE') THEN
C       DO NA=1,NNUC
        DO NA=1,MDLOOP                                     ! Vmol
        NZA = INT(ZA(NA))
        DO ND=1,3
        VLC(NA,ND)=VLC(NA,ND)+(FRC2(NA,ND)+FRC1(NA,ND))    ! new velocities
     *                          * DTMD/(2.D0*ZCM(NZA))
        ENDDO
        ENDDO

        DO NA=1,NNUC
        DO ND=1,3
        FRC1(NA,ND) = FRC2(NA,ND)
        ENDDO
        ENDDO

        IF( MOD(IMD,NRVLC).EQ.0 ) THEN
          CALL REVLC( VLC )                                ! rescale velocities
        ENDIF

      ENDIF

      ENDIF                                                ! IF(MYID.EQ.0) THEN                

      IF(NQMMM.EQ.'QMMM') THEN

        CALL SSINT( WFA,WFB,IMD,SSE )

        IF(MYID.EQ.0) THEN                

        CALL WSS( SSE )

        CALL CS( SCRD )
        IF(NMRDF.NE.0) THEN
        CALL RDF(CRD,SCRD,IMD)                             ! radial distribution function
        ENDIF
      
        IF(MOD(IMD,NCHK).EQ.0) THEN
        CALL AGAIN(2,IMD)                                  ! 2 = WRITE RESTART FILE
        ENDIF

        ENDIF                                              ! IF(MYID.EQ.0) THEN                

      ENDIF

      IF(MOD(IMD,NCHK).EQ.0) THEN
      IF(MYID.EQ.0) THEN                
      CALL DATOUT(CRD)
      ENDIF
      ENDIF                                                ! IF(MYID.EQ.0) THEN                

10    CONTINUE

      ENDIF

9     FORMAT( 12X,3F12.6 )
98    FORMAT( I6,2X,A2,I6,F6.1,F6.1,3F12.6 )
99    FORMAT( F6.1,F6.1,3F12.6 )
999   FORMAT( F6.1,3F12.6 )
9999  FORMAT( 3X,3F12.6 )

      RETURN
      END

C-----------------------------------------------------------------
C     SUBROUTINE output input data to make orb.dat for gaussian98
C-----------------------------------------------------------------

      SUBROUTINE DATOUT( CRD )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"                          ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX =  80 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NORA =  49 )
C     PARAMETER ( NORB =  1  )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / GRID2 / DX, DY, DZ

      DIMENSION CRD(  NNUC,3 )

      REWIND(22)
      IF(EXC.EQ.'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                .OR. EXC.EQ.'RXalpha') THEN
      WRITE(22,*) '#p   RHF/3-21g'
      ELSEIF(EXC.EQ.'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                .OR. EXC.EQ.'UXalpha') THEN
      WRITE(22,*) '#p   UHF/3-21g'
      ENDIF
      WRITE(22,*) '#p  Cube=(cards,orbitals)'
      WRITE(22,*) '#p  Guess=Mix'
      WRITE(22,*) '#p  Nosymm'
      WRITE(22,*) '#p  Units=au'
      WRITE(22,*)
      WRITE(22,*) 'gaussian input data'
      WRITE(22,*)
 
      NZVAL = 0
      IF(EXC.EQ.'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                .OR. EXC.EQ.'RXalpha') THEN

      DO I=1,NNUC
      NZVAL = NZVAL+ZVAL(I)
      ENDDO
      N=NZVAL-2*NORA

      WRITE(22,85) N 

      ELSEIF(EXC.EQ.'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                .OR. EXC.EQ.'UXalpha') THEN
      
      DO I=1,NNUC
      NZVAL = NZVAL+ZVAL(I)
      ENDDO
      N=NZVAL-NORA-NORB

      L=NORA-NORB

      IF(L.EQ.0) THEN
      WRITE(22,85) N 
      ELSEIF(L.EQ.1) THEN
      WRITE(22,86) N
      ENDIF

      ENDIF
      
      DO I=1,NNUC
        IF(I.LT.10) THEN 
          IF(ZA(I).EQ.8) THEN
          WRITE(22,87)I,CRD(I,1),CRD(I,2),CRD(I,3)
          ELSEIF(ZA(I).EQ.1) THEN
          WRITE(22,88)I,CRD(I,1),CRD(I,2),CRD(I,3)
          ELSEIF(ZA(I).EQ.6) THEN        
          WRITE(22,89)I,CRD(I,1),CRD(I,2),CRD(I,3)
          ELSEIF(ZA(I).EQ.17) THEN
          WRITE(22,90)I,CRD(I,1),CRD(I,2),CRD(I,3)
          ENDIF
        ELSEIF(I.LT.100) THEN
          IF(ZA(I).EQ.8) THEN
          WRITE(22,91)I,CRD(I,1),CRD(I,2),CRD(I,3)
          ELSEIF(ZA(I).EQ.1) THEN
          WRITE(22,92)I,CRD(I,1),CRD(I,2),CRD(I,3)
          ELSEIF(ZA(I).EQ.6) THEN        
          WRITE(22,93)I,CRD(I,1),CRD(I,2),CRD(I,3)
          ELSEIF(ZA(I).EQ.17) THEN
          WRITE(22,94)I,CRD(I,1),CRD(I,2),CRD(I,3)
          ENDIF
        ENDIF

      ENDDO

      WRITE(22,*)
      WRITE(22,*)'orb.dat'
 
      AX=-1.D0*(DBLE(NMAXX)+1.D0)/2.D0*DX
      AY=-1.D0*(DBLE(NMAXY)+1.D0)/2.D0*DY
      AZ=-1.D0*(DBLE(NMAXZ)+1.D0)/2.D0*DZ
      
      WRITE(22,95) AX,AY,AZ
      WRITE(22,96) NMAXX,DX
      WRITE(22,97) NMAXY,DY
      WRITE(22,98) NMAXZ,DZ
      WRITE(22,*) 'VALENCE'
      WRITE(22,*) 

85    FORMAT(I2,' 1')
86    FORMAT(I2,' 2')
87    FORMAT('O' ,I1,2X,3F12.6)
88    FORMAT('H' ,I1,2X,3F12.6)
89    FORMAT('C' ,I1,2X,3F12.6)
90    FORMAT('Cl',I1,2X,3F12.6)
91    FORMAT('O' ,I2,2X,3F12.6)
92    FORMAT('H' ,I2,2X,3F12.6)
93    FORMAT('C' ,I2,2X,3F12.6)
94    FORMAT('Cl',I2,2X,3F12.6)
95    FORMAT('-1,',F10.5,',',F10.5,',',F10.5)
96    FORMAT(    '-',I2,',',F10.5,',0.0,0.0')
97    FORMAT(    ' ',I2,',0.0,',F10.5,',0.0')
98    FORMAT(    ' ',I2,',0.0,0.0,',F10.5   )

      RETURN
      END

C-------------------------------------------
C
C     SUBROUTINE RESCALE VELOCITIES
C
C-------------------------------------------

      SUBROUTINE REVLC( VLC )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
 
      include "mpif.h"                          ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NNUC =  36 )

      DIMENSION SUM(3)
      DIMENSION VLC(NNUC,3)

      COMMON / ZMASS / ZCM(100)
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV


      BC = 3.1668303D-06                 ! Boltzmann constant ( au/K )

C----- rescale momentum -----

      ZSUM = 0.D0

      DO NA=1,NNUC
      NZA  = INT(ZA(NA))
      ZSUM = ZSUM + ZCM(NZA) 
      ENDDO

      DO ND=1,3
      SUM(ND) = 0.D0
      ENDDO

      DO NA=1,NNUC
      NZA = INT(ZA(NA))
      DO ND=1,3
        SUM(ND) = SUM(ND) + ZCM(NZA)*VLC(NA,ND)
      ENDDO
      ENDDO

      DO ND=1,3
      DO NA=1,NNUC
        VLC(NA,ND) = VLC(NA,ND) - SUM(ND) / ZSUM
      ENDDO
      ENDDO

C----- rescale velocity -----

      SUMX = 0.D0
      SUMY = 0.D0
      SUMZ = 0.D0

      DO NA=1,NNUC
        NZA  = INT(ZA(NA))
        SUMX = SUMX + ZCM(NZA)*VLC(NA,1)**2
        SUMY = SUMY + ZCM(NZA)*VLC(NA,2)**2
        SUMZ = SUMZ + ZCM(NZA)*VLC(NA,3)**2
      ENDDO
      SUMA = SUMX + SUMY + SUMZ

C     ALX = DSQRT( DBLE(NNUC)*BC*TEMP / SUMX )
C     ALY = DSQRT( DBLE(NNUC)*BC*TEMP / SUMY )
C     ALZ = DSQRT( DBLE(NNUC)*BC*TEMP / SUMZ )
      AL  = DSQRT( 3.D0*DBLE(NNUC)*BC*TEMP / SUMA )

      DO NA=1,NNUC
C       VLC(NA,1) = ALX * VLC(NA,1) 
C       VLC(NA,2) = ALY * VLC(NA,2) 
C       VLC(NA,3) = ALZ * VLC(NA,3) 
        VLC(NA,1) = AL  * VLC(NA,1) 
        VLC(NA,2) = AL  * VLC(NA,2) 
        VLC(NA,3) = AL  * VLC(NA,3) 
      ENDDO

      RETURN
      END

C---------------------------------------------------
C     SUBROUTINE READ GAUSSIAN 94 OUTPUT
C---------------------------------------------------

      SUBROUTINE TRWF( WFA,WFB,WFNA,WFNB )
       
      IMPLICIT REAL*8 ( A-H,O-Z )
      IMPLICIT INTEGER*4 ( I-N )
      
      include "mpif.h"                          ! mpi
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX =  80 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NORA =  49 )
C     PARAMETER ( NORB =   1 )
C     PARAMETER ( MORA =  49 )
C     PARAMETER ( MORB =   1 )
      PARAMETER ( PI   =   3.14159265358979323D0 )
C     PARAMETER ( ALPHA = 3.D-1 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / AVRHO / RHOAV(NMAXX/NX,NMAXY/NY,NMAXZ/NZ)

      DIMENSION WFA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION WFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
C     DIMENSION WFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,1 )
      DIMENSION WFNA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION WFNB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION NRDA( NORA )
      DIMENSION NRDB( NORB )
      
      CHARACTER NAME1*11,NAME2*13

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C---- COMPUTE PARAMETERS ----

      DV = DX*DY*DZ

C----- open orbital files -----

      write (NAME1, '("orb.dat",   i4.4)') MYID
      write (NAME2, '("avrho.dat", i4.4)') MYID                 ! takahasi 2013.07.30

      OPEN(21,FILE=NAME1,STATUS='OLD',    FORM='UNFORMATTED')
      OPEN(24,FILE=NAME2,STATUS='UNKNOWN',FORM='UNFORMATTED')   ! takahasi 2013.07.30

C----- read -----

C     READ(21,*)
C     READ(21,*)
C     READ(21,*) N0,X0,Y0,Z0
C     READ(21,*) N1,X1,Y1,Z1
C     READ(21,*) N2,X2,Y2,Z2
C     READ(21,*) N3,X3,Y3,Z3

C     DO 10  L=1,NNUC
C     READ(21,*)
C  10 CONTINUE 

      IF(EXC.EQ.'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                .OR. EXC.EQ.'RXalpha') THEN

C       NSTEP = INT(NORA/10.D0) + 1

C       DO 20 L=1,NSTEP 
C         READ(21,*)
C  20   CONTINUE

C       DO 30 I=1,N1
C       DO 30 J=1,N2
C         READ(21,'(6E13.5)') ((WFA(I,J,K,N),N=1,NORA),K=1,N3)
C  30   CONTINUE

      DO N=1,NORA
      DO K=1,JZ
      DO J=1,JY
      DO I=1,JX
       READ(21) WFA(I,J,K,N)  
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ELSEIF(EXC.EQ.'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                    .OR. EXC.EQ.'UXalpha') THEN

C       NSTEP = INT((NORA+NORB)/10.D0) + 1

C       DO 25 L=1,NSTEP 
C         READ(21,*)
C  25   CONTINUE

C       DO 35 I=1,N1
C       DO 35 J=1,N2
C         READ(21,'(6E13.5)') ((WFA(I,J,K,N),N=1,NORA),
C    *                     (WFB(I,J,K,N),N=1,NORB),K=1,N3)
C  35   CONTINUE

      DO N=1,NORA
      DO K=1,JZ
      DO J=1,JY
      DO I=1,JX
       READ(21) WFA(I,J,K,N)  
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO N=1,NORB
      DO K=1,JZ
      DO J=1,JY
      DO I=1,JX
       READ(21) WFB(I,J,K,N)  
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C     SUM = 0.D0
C     DO I=1,JX
C     DO J=1,JY
C     DO K=1,JZ
C      READ(24) RHOAV(I,J,K)  
C      SUM = SUM + RHOAV(I,J,K)*DV
C     ENDDO
C     ENDDO
C     ENDDO
C     CALL MPI_REDUCE(SUM,SUMALL,1,MPI_DOUBLE_PRECISION,
C    *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
C     if(myid.eq.0) then
C     write(*,*) 'sum=',SUMALL
C     endif

      ENDIF

C     REWIND(24)
      DO K = 1,JZ  
      DO J = 1,JY
      DO I = 1,JX
C       WRITE(24) WFA(I,J,K,1)
      ENDDO
      ENDDO
      ENDDO

C----- Normalize for alpha spin -----

      DO 40 K = 1, NORA      

        ANORM = 0.D0

        DO 50 L = 1,JZ      
        DO 50 M = 1,JY      
        DO 50 N = 1,JX      
          ANORM = ANORM + WFA( N,M,L,K )**2*DV
50      CONTINUE

        CALL MPI_REDUCE(ANORM,BNORM,1,MPI_DOUBLE_PRECISION,
     *                  MPI_SUM,0,MPI_COMM_WORLD,IERR)

      if(myid.eq.0) then
      write(*,*) 'anorm=',BNORM
      endif
        ANORM = DSQRT(BNORM)

        CALL MPI_BCAST(ANORM,1,MPI_DOUBLE_PRECISION,
     *                 0,MPI_COMM_WORLD,IERR)

        DO 60 L = 1,JZ      
        DO 60 M = 1,JY      
        DO 60 N = 1,JX      
          WFA( N,M,L,K )  = WFA( N,M,L,K ) / ANORM 
          WFNA( N,M,L,K ) = WFA( N,M,L,K ) 
60      CONTINUE

40    CONTINUE

      DO K = 1, NORA
        NRDA(K) = K
      ENDDO

C     CALL  DIAG(      NRDA,NORA,WFNA,WFA )
      CALL MDIAG (     NRDA,NORA,WFNA,WFA )
C     CALL MDIAG_BLK ( NRDA,NORA,WFNA,WFA )

C----- normalize for beta spin -----

      IF(EXC.EQ.'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                .OR. EXC.EQ.'UXalpha') THEN

        DO 70 K = 1, NORB      

          ANORM = 0.D0

          DO 80 L = 1,JZ      
          DO 80 M = 1,JY      
          DO 80 N = 1,JX      
            ANORM = ANORM + WFB( N,M,L,K )**2*DV
80        CONTINUE

        CALL MPI_REDUCE(ANORM,BNORM,1,MPI_DOUBLE_PRECISION,
     *                  MPI_SUM,0,MPI_COMM_WORLD,IERR)

        ANORM = DSQRT(BNORM)

C       IF(MYID.EQ.0) THEN
C         write(*,*) 'anorm =',ANORM
C       ENDIF

        CALL MPI_BCAST(ANORM,1,MPI_DOUBLE_PRECISION,
     *                 0,MPI_COMM_WORLD,IERR)

          DO 90 L = 1,JZ           ! 2005.07.25 takahashi
          DO 90 M = 1,JY      
          DO 90 N = 1,JX      
            WFB( N,M,L,K ) = WFB( N,M,L,K ) / ANORM
            WFNB( N,M,L,K ) = WFB( N,M,L,K ) 
90        CONTINUE

70      CONTINUE

        DO K = 1, NORB
          NRDB(K) = K
        ENDDO

C       CALL  DIAG(      NRDB,NORB,WFNB,WFB )
        CALL MDIAG (     NRDB,NORB,WFNB,WFB )
C       CALL MDIAG_BLK ( NRDB,NORB,WFNB,WFB )

      ENDIF

      RETURN
      END

C---------------------------------------------------
C     SUBROUTINE OPEN FILES 
C---------------------------------------------------

      SUBROUTINE OPF
       
      IMPLICIT REAL*8 ( A-H,O-Z )
      IMPLICIT INTEGER*4 ( I-N )
      
      include "mpif.h"          ! mpi
      include "QMpara.i"

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF

C----- OPEN FILES -----

      OPEN(20,FILE='qm.dat',STATUS='OLD')                               ! open qm recipe file
C     orb.dat ( unit=21 ) will be opened in SUBROUTINE INIT
      OPEN(22,FILE='g98.orb.inp',STATUS='UNKNOWN')                      
C     OPEN(23,FILE='vcou.dat',STATUS='OLD')                             ! Vmol

      OPEN(30,FILE='ENE.DAT',STATUS='UNKNOWN')                          ! open output files 
      OPEN(31,FILE='DNS.DAT',STATUS='UNKNOWN')
      OPEN(32,FILE='CRD.DAT',STATUS='UNKNOWN')
      OPEN(33,FILE='VLC.DAT',STATUS='UNKNOWN')
      OPEN(34,FILE='FRC.DAT',STATUS='UNKNOWN')

      OPEN(40,FILE='/home3/takahasi/PS_DATA/H/HD.DAT',                  ! open pseudopotential files 
     *     STATUS='OLD')
      OPEN(41,FILE='/home3/takahasi/PS_DATA/H/HKB.DAT',
     *     STATUS='OLD')
 
      OPEN(42,FILE='/home3/takahasi/PS_DATA/C/CD.DAT',
     *     STATUS='OLD')
      OPEN(43,FILE='/home3/takahasi/PS_DATA/C/CKB.DAT',
     *     STATUS='OLD')

      OPEN(44,FILE='/home3/takahasi/PS_DATA/O/OD.DAT',
     *     STATUS='OLD')
      OPEN(45,FILE='/home3/takahasi/PS_DATA/O/OKB.DAT',
     *     STATUS='OLD')
     
      OPEN(66,FILE='/home3/takahasi/PS_DATA/F/FD.DAT'             !  2008.03.21
     *              ,STATUS='OLD')                                !  2003.03.21
      OPEN(67,FILE='/home3/takahasi/PS_DATA/F/FKB.DAT'            !  2008.03.21
     *              ,STATUS='OLD')                                !  2003.03.21

      OPEN(50,FILE='/home3/takahasi/PS_DATA/CL/CL_BHS/CL2/CLDBHS.DAT'
     *              ,STATUS='OLD')
      OPEN(51,FILE='/home3/takahasi/PS_DATA/CL/CL_BHS/CL2/CLKBHS.DAT'
     *              ,STATUS='OLD')
     
      OPEN(49,FILE='/home3/takahasi/PS_DATA/N/ND.DAT'
     *              ,STATUS='OLD')
      OPEN(53,FILE='/home3/takahasi/PS_DATA/N/NKB.DAT'
     *              ,STATUS='OLD')
 
      OPEN(68,FILE='/home3/takahasi/PS_DATA/P/P_BHS/PD.DAT'       ! 2008.03.21
     *              ,STATUS='OLD')                                ! 2003.03.21
      OPEN(69,FILE='/home3/takahasi/PS_DATA/P/P_BHS/PKB.DAT'      ! 2008.03.21
     *              ,STATUS='OLD')                                ! 2003.03.21

      OPEN(40,FILE=PS_DIR//'/H/HD.DAT',              ! open pseudopotential files 
     *     STATUS='OLD')
      OPEN(41,FILE=PS_DIR//'/H/HKB.DAT',
     *     STATUS='OLD')
 
      OPEN(42,FILE=PS_DIR//'/C/CD.DAT',
     *     STATUS='OLD')
      OPEN(43,FILE=PS_DIR//'/C/CKB.DAT',
     *     STATUS='OLD')

      OPEN(44,FILE=PS_DIR//'/O/OD.DAT',
     *     STATUS='OLD')
      OPEN(45,FILE=PS_DIR//'/O/OKB.DAT',
     *     STATUS='OLD')
     
      OPEN(66,FILE=PS_DIR//'/F/FD.DAT'             !  2008.03.21
     *              ,STATUS='OLD')                 !  2003.03.21
      OPEN(67,FILE=PS_DIR//'/F/FKB.DAT'            !  2008.03.21
     *              ,STATUS='OLD')                 !  2003.03.21

      OPEN(50,FILE=PS_DIR//'/CL/CL_BHS/CL2/CLDBHS.DAT'
     *              ,STATUS='OLD')
      OPEN(51,FILE=PS_DIR//'/CL/CL_BHS/CL2/CLKBHS.DAT'
     *              ,STATUS='OLD')
     
      OPEN(49,FILE=PS_DIR//'/N/ND.DAT'
     *              ,STATUS='OLD')
      OPEN(53,FILE=PS_DIR//'/N/NKB.DAT'
     *              ,STATUS='OLD')
 
      OPEN(68,FILE=PS_DIR//'/P/P_BHS/PD.DAT'       ! 2008.03.21
     *              ,STATUS='OLD')                 ! 2003.03.21
      OPEN(69,FILE=PS_DIR//'/P/P_BHS/PKB.DAT'      ! 2008.03.21
     *              ,STATUS='OLD')                 ! 2003.03.21

C     OPEN(75,FILE=PS_DIR//'/CA/CAD.DAT'       ! 2016.01.28
C    *              ,STATUS='OLD')             ! 2016.01.28
      OPEN(76,FILE=PS_DIR//'/CA/CAKB.DAT'      ! 2016.01.28
     *              ,STATUS='OLD')             ! 2016.01.28
C     OPEN(77,FILE=PS_DIR//'/MN/MNKBB.DAT'     ! 2016.01.28
C    *              ,STATUS='OLD')             ! 2016.01.28
      OPEN(77,FILE=PS_DIR//'/MN/MN_PCC/MNKB.DAT'     ! 2016.01.28
     *              ,STATUS='OLD')                   ! 2016.01.28

C     OPEN(40,FILE='/stfs1/uhome/y00372/PS_DATAA/H/HD.DAT', ! open pseudopotential files 
C    *     STATUS='OLD')
C     OPEN(41,FILE='/stfs1/uhome/y00372/PS_DATAA/H/HKB.DAT',
C    *     STATUS='OLD')
C
C     OPEN(42,FILE='/stfs1/uhome/y00372/PS_DATAA/C/CD.DAT',
C    *     STATUS='OLD')
C     OPEN(43,FILE='/stfs1/uhome/y00372/PS_DATAA/C/CKB.DAT',
C    *     STATUS='OLD')

C     OPEN(44,FILE='/stfs1/uhome/y00372/PS_DATAA/O/OD.DAT',
C    *     STATUS='OLD')
C     OPEN(45,FILE='/stfs1/uhome/y00372/PS_DATAA/O/OKB.DAT',
C    *     STATUS='OLD')
C    
C     OPEN(50,FILE=
C    *    '/stfs1/uhome/y00372/PS_DATAA/CL/CL_BHS/CL2/CLDBHS.DAT',
C    *     STATUS='OLD')
C     OPEN(51,FILE=
C    *    '/stfs1/uhome/y00372/PS_DATAA/CL/CL_BHS/CL2/CLKBHS.DAT',
C    *     STATUS='OLD')
C    
C     OPEN(49,FILE='/stfs1/uhome/y00372/PS_DATAA/N/ND.DAT'
C    *              ,STATUS='OLD')
C     OPEN(53,FILE='/stfs1/uhome/y00372/PS_DATAA/N/NKB.DAT'
C    *              ,STATUS='OLD')
C
C     OPEN(68,FILE='/stfs1/uhome/y00372/PS_DATAA/P/P_BHS/PD.DAT'   !  2008.03.21
C    *              ,STATUS='OLD')                                !  2003.03.21
C     OPEN(69,FILE='/stfs1/uhome/y00372/PS_DATAA/P/P_BHS/PKB.DAT'  !  2008.03.21
C    *              ,STATUS='OLD')                                !  2003.03.21

C     OPEN(75,FILE='/stfs1/uhome/y00372/PS_DATAA/CA/CAD.DAT'       !  2016.01.28
C    *              ,STATUS='OLD')                                !  2016.01.28
C     OPEN(76,FILE='/stfs1/uhome/y00372/PS_DATAA/CA/CAKB.DAT'      !  2016.01.28
C    *              ,STATUS='OLD')                                !  2016.01.28
C     OPEN(77,FILE='/stfs1/uhome/y00372/PS_DATAA/MN/MNKBB.DAT'     !  2016.01.28
C    *              ,STATUS='OLD')                                !  2016.01.28
C     OPEN(77,FILE='/stfs1/uhome/y00372/PS_DATAA/MN/MN_PCC/MNKB.DAT' ! 2016.01.28
C    *              ,STATUS='OLD')                                  ! 2016.01.28

C     OPEN(40,FILE='PS_DATAA/H/HD.DAT',STATUS='OLD')                 ! open pseudopotential files
C     OPEN(41,FILE='PS_DATAA/H/HKB.DAT',STATUS='OLD')
 
C     OPEN(42,FILE='PS_DATAA/C/CD.DAT',STATUS='OLD')
C     OPEN(43,FILE='PS_DATAA/C/CKB.DAT',STATUS='OLD')

C     OPEN(44,FILE='PS_DATAA/O/OD.DAT',STATUS='OLD')
C     OPEN(45,FILE='PS_DATAA/O/OKB.DAT',STATUS='OLD')
     
C     OPEN(50,FILE='PS_DATAA/CL/CL_BHS/CL2/CLDBHS.DAT',STATUS='OLD')
C     OPEN(51,FILE='PS_DATAA/CL/CL_BHS/CL2/CLKBHS.DAT',STATUS='OLD')
     
C     OPEN(49,FILE='PS_DATAA/N/ND.DAT',STATUS='OLD')
C     OPEN(53,FILE='PS_DATAA/N/NKB.DAT',STATUS='OLD')
     
      OPEN(46,FILE='hps.dat',STATUS='UNKNOWN')
      OPEN(47,FILE='cps.dat',STATUS='UNKNOWN')
      OPEN(48,FILE='ops.dat',STATUS='UNKNOWN')
      OPEN(52,FILE='clps.dat',STATUS='UNKNOWN')
      OPEN(54,FILE='nps.dat',STATUS='UNKNOWN')
      OPEN(57,FILE='caps.dat',STATUS='UNKNOWN')
      OPEN(58,FILE='mnps.dat',STATUS='UNKNOWN')

      OPEN(85,FILE='link.dat',STATUS='UNKNOWN')

      RETURN
      END

C----- OPEN FILES FOR MM SUBSYSTEM ----

      SUBROUTINE MMOPF
       
      IMPLICIT REAL*8 ( A-H,O-Z )
      IMPLICIT INTEGER*4 ( I-N )
 
      include "mpif.h"                          ! mpi
      include 'QMpara.i'                        ! Vmol

      CHARACTER NO1,NO2,KAZU*5

C     PARAMETER ( NNUC =  36 )
      
      COMMON / NRDF1 / NMRDF,NRDF(NNUC)
   
      OPEN(18,FILE='chk.out',STATUS='UNKNOWN')
      OPEN(55,FILE='PCH.DAT',STATUS='UNKNOWN')
      OPEN(60,FILE='QMFS.DAT',STATUS='UNKNOWN')
      OPEN(61,FILE='QMFN.DAT',STATUS='UNKNOWN')
      OPEN(63,FILE='SSE.DAT',STATUS='UNKNOWN')
      OPEN(64,FILE='SSE-EDF.DAT',STATUS='UNKNOWN')     ! 2005.11.17 takahashi
      OPEN(65,FILE='EDF.DAT',STATUS='UNKNOWN')
C     OPEN(66,FILE='EDF_QM.DAT',STATUS='UNKNOWN')
C     OPEN(67,FILE='EDF_MM.DAT',STATUS='UNKNOWN')
      OPEN(70,FILE='SCRD.DAT',STATUS='UNKNOWN')
      OPEN(80,FILE='ENE-MM.DAT',STATUS='UNKNOWN')
      OPEN(81,FILE='ENE-MM-1.DAT',STATUS='UNKNOWN')
      OPEN(99,FILE='CS.DAT',STATUS='UNKNOWN')

      IF(NMRDF.NE.0) THEN
      DO J = 1,NMRDF
        K=NRDF(J)
        N=J+100
        IF(J.LT.10) THEN
        NO1 = CHAR(J+48)
        KAZU = 'RDF'//NO1
        ELSEIF(J.LT.100) THEN
        IN = J/10
        NO2 = CHAR(IN+48)
        IN = MOD(J,10)
        NO1 = CHAR(IN+48)
        KAZU = 'RDF'//NO2//NO1
        ENDIF
        OPEN(N,FILE=KAZU,STATUS='UNKNOWN')
      ENDDO
      ENDIF

      RETURN
      END

C---------------------------------------------------
C     SUBROUTINE CLOSE FILES 
C---------------------------------------------------

      SUBROUTINE CLF
       
      IMPLICIT REAL*8 ( A-H,O-Z )
      IMPLICIT INTEGER*4 ( I-N )
      
      include "mpif.h"                          ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NNUC =  36 )
   
      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / NRDF1 / NMRDF,NRDF(NNUC)

      CLOSE(20)                                                 ! close recipe file
      IF(NQMMM .NE. 'MM') THEN                                  ! Vmol
      CLOSE(21)                                                 ! close g94 output 
      ENDIF                                                     ! Vmol
      IF( NOPT .EQ. 'MD' ) THEN
      CLOSE(22)                                                       
      ENDIF
      CLOSE(23)                                                 ! Vmol

      CLOSE(30)                                                 ! close output files 
      CLOSE(31)                                         
      CLOSE(32)                                        
      CLOSE(33)                                            
      CLOSE(34)                                            

      CLOSE(40)                                                 ! close pseudopotential files
      CLOSE(41)                                     
      CLOSE(42)                              
      CLOSE(43)
      CLOSE(44)
      CLOSE(45)
      CLOSE(50)
      CLOSE(51)
      CLOSE(49)
      CLOSE(53)

      CLOSE(46)                              
      CLOSE(47)
      CLOSE(48)
      CLOSE(52)
      CLOSE(54)

      IF(NQMMM.EQ.'QMVD') CLOSE(UNIT=220)
      IF(NQMMM.EQ.'QMVC') CLOSE(UNIT=230)

      IF(NQMMM .EQ. 'QMMM' .OR. NQMMM .EQ. 'LINK') THEN         ! Vmol

        CLOSE(18)
        CLOSE(55)
        CLOSE(60)
        CLOSE(61)
        CLOSE(63)
        CLOSE(64)                                               ! 2005.11.17 takahashi
        CLOSE(65)
        CLOSE(66)
        CLOSE(67)
        CLOSE(70)
        CLOSE(80)
        CLOSE(81)
        CLOSE(99)
        IF(NMRDF.NE.0) THEN
        DO I = 1,NMRDF
        N=I+100
        CLOSE(N)
        ENDDO
        ENDIF

      ENDIF

      RETURN
      END

C---------------------------------------------------------
C    BLOCK DATA FOR CORE MASS
C---------------------------------------------------------

      BLOCK DATA ZMSS
      IMPLICIT REAL*8 ( A-H,O-Z )      
      COMMON / ZMASS / ZCM(100)

      DATA ZCM(1),ZCM(6),ZCM(7),ZCM(8),ZCM(17),ZCM(20) 
C    &    / 1.83636D+03,2.1891D+04,2.91591D+04,6.46193D+4 /
     &    / 3.67272D+03,2.1891D+04,2.55307D+04,2.91591D+04,
     &      6.46193D+04,7.3598D+04 /

      END

C---------------------------------------------------
C
C   SUBROUTINE COMPUTING ELECTRONIC GROUND STATE      
C
C---------------------------------------------------

      SUBROUTINE GRND( TRFA,TRFB,TRFNA,TRFNB,CRD,TFRC,PENG,IMD ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"                          ! mpi
      include 'mpi.i'                           ! mpi

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'                         ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC = 36 )
C     PARAMETER ( NORA = 49 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / MPIFRC / FRCM( NNUC,NDIM )

      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / FRCE / FRC(NNUC,NDIM)
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),     ! Vmol
     *                 NMMSW(maxatm),MMID(NNUC)

      DIMENSION TRFA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )     ! trial wavefunctions for alpha spin
      DIMENSION TRFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )     ! trial wavefunctions for beta  spin
C     DIMENSION TRFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,1    )     ! trial wavefunctions for beta  spin
      DIMENSION TRFNA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )   
      DIMENSION TRFNB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )   

      DIMENSION CRD(   NNUC,3 )
      DIMENSION TFRC(  NNUC,3 )
      DIMENSION VNUC(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )         
      DIMENSION VNUCN( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )         

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      NDAT = NNUC*NDIM

C----- INPUT NUCLEAR CONFIGRATION -----

C     IF(IMD.EQ.0) THEN
      CALL INTCRD( CRD,IMD ) 
C     ENDIF
     
C     DO I=1,NNUC
C     DO J=1,NDIM
C       PNUC(I,J)=CRD(I,J)
C     ENDDO
C     ENDDO

C----- initialize -----

      DO J=1,NDIM
      DO I=1,NNUC
        FRC( I,J) = 0.D0
        FRCM(I,J) = 0.D0
      ENDDO
      ENDDO

      IF(NQMMM .EQ. 'MM') GOTO 10               ! Vmol

C     IF(IMD.NE.0) THEN
C       IF(MYID.EQ.0) THEN
C       CALL DISMAT( PNUC )                     ! comment 2005.11.22 takahashi
C       ENDIF      
C     ENDIF

C-------------------

      IF(IMD.EQ.0) THEN

C       WRITE(*,*)
        IF(MYID.EQ.0) THEN
          WRITE(*,*)
        ENDIF
          CALL NPOT                             ! nuclear-nuclear energy and force for non-periodic system
          CALL IVCOU(    PNUC )                 ! Vmol
C         CALL IVCOU_NW( PNUC )                 ! Vmol

        IF(DGF.EQ.'DG') THEN  
          CALL DG(VNUC)                         ! double grid method for non-periodic system
C         CALL FUZZYC                           ! Vmol fuzzy cell method for Hartree potential
        ELSEIF(DGF.EQ.'DG4') THEN
          IF(NQMMM.EQ.'LINK') THEN
C          CALL DG4QM                    
           CALL DG4_CUB( VNUCN )                 ! Warning: cannot be applied to systems with Link atoms 
           CALL DG4L_Y(VNUCN,VNUC)    
           CALL FUZZYC                           ! Vmol fuzzy cell method for estimation of oxidation number
          ELSE
           CALL FUZZYC                           ! Vmol fuzzy cell method for estimation of oxidation number
C          CALL DG (     VNUC)    
           CALL DG4(     VNUC)    
C          CALL DG4_CUB( VNUC)                   ! takahasi 2013.04.24 modified DG4
C          CALL DG4_MOD( VNUC)                   ! takahasi 2013.04.24 modified DG4
C          CALL DG4_MOD2(VNUC)                   ! takahasi 2013.04.30 modified DG4 with bilinear interpolation
          ENDIF
C         CALL FUZZYC                           ! Vmol fuzzy cell method for Hartree potential
        ELSEIF(DGF.EQ.'DG6') THEN
          CALL DG6(VNUC)          
C         CALL FUZZYC                           ! Vmol fuzzy cell method for Hartree potential
        ENDIF      

        CALL PCC                                ! read and distribute PCC data

          IF(MYID.EQ.0) THEN
          IF( NOPT .EQ. 'MD' ) THEN
          IF( NRST .EQ. 'NEW' ) THEN
            WRITE(*,*)                             
            WRITE(*,*) '  MD TIME( fs ) = ',0.D0
            WRITE(*,*)                             
          ENDIF
          ENDIF
          ENDIF

      ELSEIF(FREEZE.EQ.'FALSE') THEN

        IF(MYID.EQ.0) THEN
C       CALL DISMAT ( PNUC ) 
        WRITE(*,*)
        ENDIF
        CALL NPOT                               ! nuclear-nuclear energy and force for non-periodic system

        IF(DGF.EQ.'DG') THEN  
          CALL DG(VNUC)                         ! double grid method for non-periodic system
C         CALL FUZZYC                           ! Vmol fuzzy cell method for Hartree potential
        ELSEIF(DGF.EQ.'DG4') THEN
          CALL DG4(VNUC)    
C         CALL FUZZYC                           ! Vmol fuzzy cell method for Hartree potential
        ELSEIF(DGF.EQ.'DG6') THEN
          CALL DG6(VNUC)          
C         CALL FUZZYC                           ! Vmol fuzzy cell method for Hartree potential
        ENDIF

      ELSEIF(NQMMM.EQ.'LINK') THEN

        IF(MYID.EQ.0) THEN
C       CALL DISMAT ( PNUC ) 
        WRITE(*,*)
        ENDIF
        CALL NPOT                               ! nuclear-nuclear energy and force for non-periodic system

        IF(DGF.EQ.'DG') THEN  
C         CALL DG(VNUC)                         ! double grid method for non-periodic system
C         CALL FUZZYC                           ! Vmol fuzzy cell method for Hartree potential
        ELSEIF(DGF.EQ.'DG4') THEN
C         IF(NLINK.GT.0) CALL DG4L_Y(VNUCN,VNUC)  ! Double Grid for Link Atoms
C         CALL FUZZYC                           ! Vmol fuzzy cell method for Hartree potential
        ELSEIF(DGF.EQ.'DG6') THEN
C         CALL DG6(VNUC)          
C         CALL FUZZYC                           ! Vmol fuzzy cell method for Hartree potential
        ENDIF

      ENDIF

C----- electronic ground state -----

10    CONTINUE                                  ! Vmol

C     CALL PCC                                                 ! read and distribute PCC data

      CALL SCF( TRFA,TRFB,TRFNA,TRFNB,VNUC,PENG,IMD )

C     IF( MYID.EQ.0) THEN
C     WRITE(*,*) ' current PNUC:'
C     DO I = 1, NNUC
C       WRITE(*,99) ZA(I),ZVAL(I),PNUC(I,1),PNUC(I,2),PNUC(I,3)
C     ENDDO
C     ENDIF
99    format( 2f6.1,3f12.6 )

C----- replace -----

C     CALL MPI_REDUCE( FRC(1,1),TFRC(1,1),NDAT,MPI_DOUBLE_PRECISION,
C    *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

C     DO J=1,NDIM
C     DO I=1,NNUC
C       TFRC(I,J) = FRC(I,J)
C     ENDDO
C     ENDDO

      RETURN
      END

C--------------------------------------------
C   SUBROUTINE Internal Atomic Coordinates
C--------------------------------------------

      SUBROUTINE INTCRD( CRD,IMD ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"                          ! mpi
      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'                         ! Vmol
!     include 'atoms.i'                         ! Vmol

C     PARAMETER ( NDIM   =  3 )
C     PARAMETER ( NNUC   = 36 )
C     PARAMETER ( NLINK1 = 10 )                 ! Vmol

      PARAMETER ( NSP = 255 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),     ! Vmol
     *                 NMMSW(maxatm),MMID(NNUC)                    ! Vmol
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / GRID2 / DX, DY, DZ
C     COMMON / MMDAT / SCRD( NDIM,4*NSP ),SCRD1( NDIM,4*NSP ),
C    *                 FRCS( 4*NSP,NDIM )
      COMMON / MMDAT / SCRD( NDIM,maxatm ),SCRD1( NDIM,maxatm ),   ! Vmol
     *                 FRCS( maxatm,NDIM )                         ! Vmol

      DIMENSION CRD( NNUC,3 )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      IF(IMD.NE.0) GOTO 20

C----- X axis -----

      XMAX = CRD(1,1)
      XMIN = CRD(1,1)

      DO I = 2, NNUC
        XMAX = MAX( XMAX,CRD(I,1) )
        XMIN = MIN( XMIN,CRD(I,1) )
      ENDDO

C      Internal Coordinate      

      XADD = -0.5D0*(XMAX+XMIN)

      XADD = XADD/DX

      IF(XADD.LT.0.D0) THEN
        XADD = AINT(XADD) * DX - DX
      ELSE
        XADD = AINT(XADD) * DX
      ENDIF

C     DO I = 1, NNUC
C       PNUC(I,1) = CRD(I,1) + XADD
C     ENDDO

C----- Y axis -----

      YMAX = CRD(1,2)
      YMIN = CRD(1,2)

      DO I = 2, NNUC
        YMAX = MAX( YMAX,CRD(I,2) )
        YMIN = MIN( YMIN,CRD(I,2) )
      ENDDO

C      Internal Coordinate      

      YADD = -0.5D0*(YMAX+YMIN)

      YADD = YADD/DY

      IF(YADD.LT.0.D0) THEN
        YADD = AINT(YADD) * DY - DY
      ELSE
        YADD = AINT(YADD) * DY
      ENDIF

C     DO I = 1, NNUC
C       PNUC(I,2) = CRD(I,2) + YADD
C     ENDDO

C----- Z axis -----

      ZMAX = CRD(1,3)
      ZMIN = CRD(1,3)

      DO I = 2, NNUC
        ZMAX = MAX( ZMAX,CRD(I,3) )
        ZMIN = MIN( ZMIN,CRD(I,3) )
      ENDDO

C      Internal Coordinate      

      ZADD = -0.5D0*(ZMAX+ZMIN)

      ZADD = ZADD/DZ

      IF(ZADD.LT.0.D0) THEN
        ZADD = AINT(ZADD) * DZ - DZ
      ELSE
        ZADD = AINT(ZADD) * DZ
      ENDIF

C     DO I = 1, NNUC
C       PNUC(I,3) = CRD(I,3) + ZADD
C     ENDDO

C----- Internal Coordinates for QM molecule -----

20    CONTINUE

C     DO I = 1, NNUC
C       PNUC(I,1) = CRD(I,1) + XADD
C       PNUC(I,2) = CRD(I,2) + YADD
C       PNUC(I,3) = CRD(I,3) + ZADD
C     ENDDO

      DO I = 1, NNUC      ! no chnage in coordinates 2005.09.14 takahashi  
        PNUC(I,1) = CRD(I,1) 
        PNUC(I,2) = CRD(I,2) 
        PNUC(I,3) = CRD(I,3) 
      ENDDO

C----- Internal Coordinates for MM molecule -----

      IF(NQMMM.EQ.'QMMM') THEN
        NLOOP = 4*NSP
      ELSEIF(NQMMM.EQ.'LINK') THEN                  ! Vmol
        NLOOP = n-NNUC+NLINK                        ! Vmol
      ELSEIF(NQMMM.EQ.'MM') THEN                    ! Vmol
C       NLOOP = n-NNUC                              ! Vmol
        NLOOP = n-NNUC+NLINK                        ! Vmol
      ELSE                                          ! Vmol
        GOTO 10                                     ! Vmol
      ENDIF

C     DO I = 1, NLOOP
C       SCRD1(1,I) = SCRD(1,I) + XADD
C       SCRD1(2,I) = SCRD(2,I) + YADD
C       SCRD1(3,I) = SCRD(3,I) + ZADD
C     ENDDO

      DO I = 1, NLOOP       ! no chnage in coordinates 2005.09.14 takahashi
        SCRD1(1,I) = SCRD(1,I) 
        SCRD1(2,I) = SCRD(2,I) 
        SCRD1(3,I) = SCRD(3,I) 
      ENDDO

10    CONTINUE

C----- Write Internal Coordinates for QM molecule -----

C     WRITE(*,*)                             
C     WRITE(*,*) '    Internal Coordinates'
C     DO J=1,NNUC                                ! write coordinates
C     WRITE(*,99) ZA(J),ZVAL(J),
C    *            PNUC(J,1),PNUC(J,2),PNUC(J,3)
C     ENDDO

99    FORMAT( F6.1,F6.1,3F12.6 )

      RETURN
      END

C---------------------------------------------------
C   SUBROUTINE Distance Matrix For Atomic Nucleus       
C---------------------------------------------------

      SUBROUTINE DISMAT( PNUC ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"                          ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC =  36 )

      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV

      DIMENSION DISNUC(  NNUC,NNUC )
      DIMENSION PNUC(  NNUC,NDIM )

      CHARACTER ATOM*2

      WRITE(*,*) 
      WRITE(*,*) 'Distance Matrix [a.u.]   ' 

      DO I = 1, NNUC
      DO J = I, NNUC

      RX = PNUC( I,1 ) - PNUC( J,1 ) 
      RY = PNUC( I,2 ) - PNUC( J,2 ) 
      RZ = PNUC( I,3 ) - PNUC( J,3 ) 
      R2 = RX**2 + RY**2 + RZ**2
      R  = DSQRT( R2 )

      DISNUC(J,I) = R

      ENDDO
      ENDDO

      IMAT = INT( NNUC/5 ) + 1

      DO 10 I = 1,IMAT 

      J1MAT = 5*(I-1) +1
      IF(J1MAT .GT. NNUC) GOTO 10 

      J3MAT = J1MAT + 4
      IF(J3MAT .GT. NNUC) J3MAT=NNUC

      WRITE(*,98) ( N,N=J1MAT,J3MAT )

      DO 20 K = J1MAT, NNUC

      IF(K .LT. J3MAT) THEN
      J2MAT = K 
      ELSE
      J2MAT = J3MAT
      ENDIF

      IF(ZA(K).EQ.1.D0) THEN
      ATOM = 'H '
      ELSEIF(ZA(K).EQ.6.D0) THEN
      ATOM = 'C '
      ELSEIF(ZA(K).EQ.7.D0) THEN
      ATOM = 'N '
      ELSEIF(ZA(K).EQ.8.D0) THEN
      ATOM = 'O '
      ELSEIF(ZA(K).EQ.9.D0) THEN
      ATOM = 'F '
      ELSEIF(ZA(K).EQ.17.D0) THEN
      ATOM = 'Cl'
      ELSEIF(ZA(K).EQ.20.D0) THEN
      ATOM = 'Ca'
      ELSEIF(ZA(K).EQ.25.D0) THEN
      ATOM = 'Mn'
      ENDIF

      WRITE(*,99) I,K,ATOM,(DISNUC(K,N),N=J1MAT,J2MAT)

20    CONTINUE  
10    CONTINUE  

98    FORMAT( 18x,I2,3x,4(6x,I2,3x) )      
99    FORMAT( 2x,I2,2x,I2,2x,A2,5F11.6 )      

      RETURN
      END

C----------------------------------------------
C   SUBROUTINE Xalpha ( ELECTRONIC EXCHANGE ) 
C----------------------------------------------

      SUBROUTINE VXALP( RHO,RHOA,RHOB,VEXA,VEXB,PEXC,PEXT )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )

      PARAMETER ( PI   = 3.14159265358979323D0 )
      PARAMETER ( ALPH = 0.7D0 )                    ! empirical coefficient for Xalpha

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI3 / NIDX(1:2),NIDY(1:2),NIDZ(1:2)
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ            
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF

      DIMENSION RHO(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHOA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHOB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEXA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEXB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      COEFF1 = -(3.D0/2.D0)*ALPH*( (3.D0/PI)**(1.D0/3.D0) ) 
     *       *  (2.D0**(1.D0/3.D0))
      COEFF2 = -(9.D0/8.D0)*ALPH*( (3.D0/PI)**(1.D0/3.D0) ) 
     *       *  (2.D0**(1.D0/3.D0))

      DXYZ   = DX*DY*DZ

C----- local exchange -----

      DO 30 L = 1, JZ
      DO 20 M = 1, JY
      DO 10 N = 1, JX

        VEXA( N,M,L ) = COEFF1 * ( (RHOA( N,M,L ))**(1.D0/3.D0) ) 

10    CONTINUE
20    CONTINUE
30    CONTINUE

      PEXC = 0.D0                         ! total exchange
      PEXT = 0.D0                         ! external exchange

      IF(EXC.EQ.'UXalpha') THEN
   
      DO 35 L = 1, JZ
      DO 25 M = 1, JY
      DO 15 N = 1, JX

        VEXB( N,M,L ) = COEFF1 * ( (RHOB( N,M,L ))**(1.D0/3.D0) ) 

15    CONTINUE
25    CONTINUE
35    CONTINUE

      DO 60 L = 1, JZ
      DO 50 M = 1, JY
      DO 40 N = 1, JX

        PEXC = PEXC + ((RHOA( N,M,L ))**(4.D0/3.D0)  
     *              +  (RHOB( N,M,L ))**(4.D0/3.D0) )*DXYZ
        PEXT = PEXT + ((RHOA( N,M,L ))**(4.D0/3.D0)  
     *              +  (RHOB( N,M,L ))**(4.D0/3.D0) )*DXYZ

40    CONTINUE
50    CONTINUE
60    CONTINUE

      ELSEIF(EXC.EQ.'RXalpha') THEN
   
      DO 90 L = 1, JZ
      DO 80 M = 1, JY
      DO 70 N = 1, JX

        PEXC = PEXC + 2.D0* ((RHOA( N,M,L ))**(4.D0/3.D0))*DXYZ  
        PEXT = PEXT + 2.D0* ((RHOA( N,M,L ))**(4.D0/3.D0))*DXYZ  

70    CONTINUE
80    CONTINUE
90    CONTINUE

      ENDIF

      CALL MPI_REDUCE(PEXC,PPEXC,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(PEXT,PPEXT,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)

      PEXC = PPEXC
      PEXT = PPEXT
    
      CALL MPI_BCAST(PEXC,1,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(PEXT,1,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)

      PEXC = COEFF2 * PEXC
      PEXT = COEFF1 * PEXT

C     write(*,*) 'coeff1 =',COEFF1
C     write(*,*) 'coeff2 =',COEFF2

      RETURN
      END

C----------------------------------------------
C   SUBROUTINE PZ 
C   electronic exchange and correlation  
C   proposed by Perdew and Zunger
C   Phys. Rev. B23, (1981)5048
C----------------------------------------------

      SUBROUTINE PZ( RHOAB,VEXAB,PEXC ) ! Perdew & Zunger for spin density

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDOM =  3 )
      PARAMETER ( PI   = 3.14159265358979323D0 )

      COMMON / GRID2 / DX, DY, DZ

      DIMENSION RHOAB( NMAXX,NMAXY,NMAXZ )
      DIMENSION VEXAB( NMAXX,NMAXY,NMAXZ )
      
      DV = DX*DY*DZ

C----- coefficients -----

      A  = -0.4582D0            ! exchange
      A1 =  1.0D0               ! rs >= 1.
      A2 =  1.0529D0
      A3 =  0.3334D0
      A4 = -0.1423D0
      B1 = -0.048D0             ! rs =< 1.
      B2 =  0.0311D0
      B3 = -0.0116D0
      B4 =  0.002D0
      
C----- local/total exchange and correlation -----

      DO 10 L = 1, NMAXZ
      DO 10 M = 1, NMAXY
      DO 10 N = 1, NMAXX
          
        RS = (3.D0/(4.D0*PI))**(1.D0/3.D0)
     *     * (1.D0/(2.D0*RHOAB( N,M,L )))**(1.D0/3.D0)

        IF ( RS .GT. 1.D0 ) THEN
          DEN = ( A1+A2*RS**(0.5D0)+A3*RS ) 
          VEX  = A/RS + A4/DEN
          PEXC = PEXC + VEX*RHOAB( N,M,L )*DV
          VEXAB( N,M,L ) = VEX
     *                  - (RS/3.D0)*( -A/RS**2.D0  
     *                  - A4*( (A2/2.D0)*RS**(-0.5D0)+A3 )
     *                  / DEN**2.D0)
        ELSE 
          VEX = A/RS
     *        + B1 + B2*DLOG(RS) + B3*(RS) + B4*RS*DLOG(RS)
          PEXC = PEXC + VEX*RHOAB( N,M,L )*DV
          VEXAB( N,M,L ) = VEX
     *                  - (RS/3.D0)*( -A/RS**2 
     *                  + ( B2/RS + B3 + B4*(1.D0+DLOG(RS))) )
        ENDIF

10    CONTINUE

      RETURN
      END

C-----------------------------------------------------------
C            FUNCTIONS FOR BLYP                           
C-----------------------------------------------------------

C----- second deivative -----

      REAL*8 FUNCTION SDF(XYZ,X1,X2,X3,X4,X5,X6,X7,X8,
     *                        Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,
     *                        Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      COMMON / GRID2 / DX, DY, DZ

      A1  =  8.D0/ 5.D0
      A2  = -1.D0/ 5.D0
      A3  =  8.D0/ 315.D0
      A4  = -1.D0/ 560.D0

      AX1 = A1 / (DX**2.D0)
      AX2 = A2 / (DX**2.D0)
      AX3 = A3 / (DX**2.D0)
      AX4 = A4 / (DX**2.D0)

      AY1 = A1 / (DY**2.D0)
      AY2 = A2 / (DY**2.D0)
      AY3 = A3 / (DY**2.D0)
      AY4 = A4 / (DY**2.D0)

      AZ1 = A1 / (DZ**2.D0)
      AZ2 = A2 / (DZ**2.D0)
      AZ3 = A3 / (DZ**2.D0)
      AZ4 = A4 / (DZ**2.D0)

      AO  = -2.D0 * ( AX1+AX2+AX3+AX4 
     *          +   AY1+AY2+AY3+AY4
     *          +   AZ1+AZ2+AZ3+AZ4 )
 

          SDF = AX1 * ( X1+X2 )
     *        + AX2 * ( X3+X4 ) 
     *        + AX3 * ( X5+X6 ) 
     *        + AX4 * ( X7+X8 )
     *        + AY1 * ( Y1+Y2 )
     *        + AY2 * ( Y3+Y4 ) 
     *        + AY3 * ( Y5+Y6 ) 
     *        + AY4 * ( Y7+Y8 )
     *        + AZ1 * ( Z1+Z2 )
     *        + AZ2 * ( Z3+Z4 ) 
     *        + AZ3 * ( Z5+Z6 ) 
     *        + AZ4 * ( Z7+Z8 )
     *        + AO  *   XYZ

      RETURN
      END

C---------------- first deivative for x-axis ----------------

      REAL*8 FUNCTION DFX(X1,X2,X3,X4,X5,X6,X7,X8)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER ( I-N )      

      COMMON / GRID2 / DX, DY, DZ

      A1  = 449.D0 / 360.D0
      A2  = -11.D0 /  72.D0
      A3  =  11.D0 / 504.D0
      A4  =  -1.D0 / 560.D0

      AX1 = A1 / ( 2.D0*DX )
      AX2 = A2 / ( 2.D0*DX )
      AX3 = A3 / ( 2.D0*DX )
      AX4 = A4 / ( 2.D0*DX )

          DFX  = AX1 * ( X1-X2 )
     *         + AX2 * ( X3-X4 )
     *         + AX3 * ( X5-X6 )
     *         + AX4 * ( X7-X8 )

      RETURN
      END

C---------------- first deivative for y-axis ----------------

      REAL*8 FUNCTION DFY(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER ( I-N )      

      COMMON / GRID2 / DX, DY, DZ

      A1  = 449.D0 / 360.D0
      A2  = -11.D0 /  72.D0
      A3  =  11.D0 / 504.D0
      A4  =  -1.D0 / 560.D0

      AY1 = A1 / ( 2.D0*DY )
      AY2 = A2 / ( 2.D0*DY )
      AY3 = A3 / ( 2.D0*DY )
      AY4 = A4 / ( 2.D0*DY )

          DFY  = AY1 * ( Y1-Y2 )
     *         + AY2 * ( Y3-Y4 )
     *         + AY3 * ( Y5-Y6 )
     *         + AY4 * ( Y7-Y8 )

      RETURN
      END

C---------------- first deivative for z-axis ----------------

      REAL*8 FUNCTION DFZ(Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER ( I-N )      

      COMMON / GRID2 / DX, DY, DZ

      A1  = 449.D0 / 360.D0
      A2  = -11.D0 /  72.D0
      A3  =  11.D0 / 504.D0
      A4  =  -1.D0 / 560.D0

      AZ1 = A1 / ( 2.D0*DZ )
      AZ2 = A2 / ( 2.D0*DZ )
      AZ3 = A3 / ( 2.D0*DZ )
      AZ4 = A4 / ( 2.D0*DZ )

          DFZ  = AZ1 * ( Z1-Z2 )
     *         + AZ2 * ( Z3-Z4 )
     *         + AZ3 * ( Z5-Z6 )
     *         + AZ4 * ( Z7-Z8 )

      RETURN
      END

C--------------------- second deivative for x-axis ---------------------

      REAL*8 FUNCTION SDFX(XYZ,X1,X2,X3,X4,X5,X6,X7,X8)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER ( I-N )      

      COMMON / GRID2 / DX, DY, DZ

      A1  =  8.D0/ 5.D0
      A2  = -1.D0/ 5.D0
      A3  =  8.D0/ 315.D0
      A4  = -1.D0/ 560.D0

      AX1 = A1 / (DX**2.D0)
      AX2 = A2 / (DX**2.D0)
      AX3 = A3 / (DX**2.D0)
      AX4 = A4 / (DX**2.D0)

      AY1 = A1 / (DY**2.D0)
      AY2 = A2 / (DY**2.D0)
      AY3 = A3 / (DY**2.D0)
      AY4 = A4 / (DY**2.D0)

      AZ1 = A1 / (DZ**2.D0)
      AZ2 = A2 / (DZ**2.D0)
      AZ3 = A3 / (DZ**2.D0)
      AZ4 = A4 / (DZ**2.D0)

      AO  = -2.D0 * ( AX1+AX2+AX3+AX4 
     *          +   AY1+AY2+AY3+AY4
     *          +   AZ1+AZ2+AZ3+AZ4 )
 

          SDFX = AX1 * ( X1+X2 )
     *         + AX2 * ( X3+X4 ) 
     *         + AX3 * ( X5+X6 ) 
     *         + AX4 * ( X7+X8 )
     *         + AO  *   XYZ

      RETURN
      END

C--------------------- second deivative for y-axis ---------------------

      REAL*8 FUNCTION SDFY(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER ( I-N )      

      COMMON / GRID2 / DX, DY, DZ

      A1  =  8.D0/ 5.D0
      A2  = -1.D0/ 5.D0
      A3  =  8.D0/ 315.D0
      A4  = -1.D0/ 560.D0

      AY1 = A1 / (DY**2.D0)
      AY2 = A2 / (DY**2.D0)
      AY3 = A3 / (DY**2.D0)
      AY4 = A4 / (DY**2.D0)

          SDFY = AY1 * ( Y1+Y2 )
     *         + AY2 * ( Y3+Y4 ) 
     *         + AY3 * ( Y5+Y6 ) 
     *         + AY4 * ( Y7+Y8 )

      RETURN
      END

C--------------------- second deivative for z-axis ---------------------

      REAL*8 FUNCTION SDFZ(Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER ( I-N )      

      COMMON / GRID2 / DX, DY, DZ

      A1  =  8.D0/ 5.D0
      A2  = -1.D0/ 5.D0
      A3  =  8.D0/ 315.D0
      A4  = -1.D0/ 560.D0

      AZ1 = A1 / (DZ**2.D0)
      AZ2 = A2 / (DZ**2.D0)
      AZ3 = A3 / (DZ**2.D0)
      AZ4 = A4 / (DZ**2.D0)

          SDFZ = AZ1 * ( Z1+Z2 )
     *         + AZ2 * ( Z3+Z4 ) 
     *         + AZ3 * ( Z5+Z6 ) 
     *         + AZ4 * ( Z7+Z8 )

      RETURN
      END

C---------------------------------------------------------
C   SUBROUTINE COULOMB POTENTIAL FOR NON-PERIODIC SYSTEMS 
C   Phys.Rev.B 48 (1993) 2081 APPENDIX D
C---------------------------------------------------------

C      SUBROUTINE VCLMB( RHO )
C
C      IMPLICIT REAL*8 ( A-H,O-Z )      
C      IMPLICIT INTEGER*4 ( I-N )      
C      INTEGER*4 STATUS      
C
C      INCLUDE '/usr/include/DXMLDEF.FOR'
C
C      include 'QMpara.i'                        ! Vmol
C
CC     PARAMETER ( NMAX = 80 )
CC     PARAMETER ( NNUC =  36 )
CC     PARAMETER ( NDIM =  3 )
C
C      COMMON / GRID2 / DX, DY, DZ
C      COMMON / CLMB  / PCLMB, PEXC, PEXT
C      COMMON / RGRID / RR1(NMAX*2,NMAX*2,NMAX*2),
C     *                 RR2(NMAX*2,NMAX*2,NMAX*2)
C      COMMON / FFT / N_I, N_J, N_K, LDA_I, LDA_J
C      COMMON / OFC / VCOU( NMAX,NMAX,NMAX )
C
C      DIMENSION RHO(  NMAX,  NMAX,NMAX )
C      DIMENSION DRHO1( NMAX*2,NMAX*2,NMAX*2 )
C      DIMENSION DRHO2( NMAX*2,NMAX*2,NMAX*2 )
CC     DIMENSION TMP1( NMAX*2,NMAX*2,NMAX*2 )
CC     DIMENSION TMP2( NMAX*2,NMAX*2,NMAX*2 )
C
CC----- GRID -----
C
C      DO L = 1, NMAX*2
C      DO M = 1, NMAX*2
C      DO N = 1, NMAX*2
C
C          DRHO1( N,M,L ) =  0.D0 
C          DRHO2( N,M,L ) =  0.D0                     
C
C      ENDDO
C      ENDDO
C      ENDDO
C
C      DO L = NMAX/2+1, NMAX*3/2
C      DO M = NMAX/2+1, NMAX*3/2
C      DO N = NMAX/2+1, NMAX*3/2
C
C          DRHO1( N,M,L ) =  RHO( N-NMAX/2,M-NMAX/2,L-NMAX/2 ) 
C
C      ENDDO
C      ENDDO
C      ENDDO
C
CC----- FFT ----- 
C
C      STATUS = ZFFT_3D('R','R','F',DRHO1,DRHO2,DRHO1,DRHO2,
C     *                  N_I,N_J,N_K,LDA_I,LDA_J,1,1,1)
C
C      DO L = 1, NMAX*2
C      DO M = 1, NMAX*2
C      DO N = 1, NMAX*2
C
C          A1 = DRHO1( N,M,L )
C          A2 = DRHO2( N,M,L )
C
C          DRHO1( N,M,L ) =  A1 * RR1( N,M,L )    
C     *                   -  A2 * RR2( N,M,L )    
C          DRHO2( N,M,L ) =  A1 * RR2( N,M,L )    
C     *                   +  A2 * RR1( N,M,L )    
C
C      ENDDO
C      ENDDO
C      ENDDO
C
C      STATUS = ZFFT_3D('R','R','B',DRHO1,DRHO2,DRHO1,DRHO2,
C     *                  N_I,N_J,N_K,LDA_I,LDA_J,1,1,1)
C
CC----- replace -----
C
C      DO NZ = 1, NMAX*2
C        NZ0 = NZ + NMAX
C        IF( NZ .GT. NMAX ) NZ0 = NZ - NMAX
C        DO NY = 1, NMAX*2
C          NY0 = NY + NMAX
C          IF( NY .GT. NMAX ) NY0 = NY - NMAX
C          DO NX = 1, NMAX*2
C            NX0 = NX + NMAX
C            IF( NX .GT. NMAX ) NX0 = NX - NMAX
CC           TMP1( NX0,NY0,NZ0 ) = DRHO1( NX,NY,NZ )
CC           TMP2( NX0,NY0,NZ0 ) = DRHO2( NX,NY,NZ )
C            DRHO2( NX0,NY0,NZ0 ) = DRHO1( NX,NY,NZ )
C          ENDDO
C        ENDDO
C      ENDDO
C
CC--------------
C
C      DO L = NMAX/2+1, NMAX*3/2
C      DO M = NMAX/2+1, NMAX*3/2
C      DO N = NMAX/2+1, NMAX*3/2
C
CC         VCOU( N-NMAX/2,M-NMAX/2,L-NMAX/2 ) = TMP1( N,M,L )
C          VCOU( N-NMAX/2,M-NMAX/2,L-NMAX/2 ) = DRHO2( N,M,L )
C
C      ENDDO
C      ENDDO
C      ENDDO
C
CC----- classical coulomb potential -----   
C
C      PCLMB = 0.D0
C
C      DO L = 1, NMAX
C      DO M = 1, NMAX
C      DO N = 1, NMAX
C
C          PCLMB = PCLMB + VCOU( N,M,L )*RHO( N,M,L )*DX*DY*DZ
C
C      ENDDO
C      ENDDO
C      ENDDO
C
C      PCLMB = 0.5D0*PCLMB 
C
C      RETURN 
C      END

C-------------------------------------
C   SUBROUTINE COMPUTE PARAMETERS 
C-------------------------------------

      SUBROUTINE CONFIG
      USE POISSON_SOLVER,
     *     PSOLVER_INIT => INIT

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
      INTEGER*4 STATUS      

C     INCLUDE '/usr/include/DXMLDEF.FOR'

      include "mpif.h" 
      include 'mpi.i'       
      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
      include 'nlocd.i'
      include 'pc_crr.i'                        ! Vmol   pcc

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC =  36 )

      PARAMETER ( PI   =  3.14159265358979323D0  )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPIPOIS / ICOMM_CART

      COMMON / GRID2 / DX, DY, DZ
      COMMON / ACELL / XL, YL, ZL
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1), ! Vmol
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / ELCTRN / AX1, AX2, AY1, AY2, AZ1, AZ2, AO,
     *                  AX3, AX4, AY3, AY4, AZ3, AZ4
      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,                 ! Vmol
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,                 ! Vmol
C    *                  WNUC(NFUZZY,NMAX,NMAX,NMAX),RSIZE(NNUC),           ! Vmol
C    *                  WNUC(NFUZZY,NMAX/NX,NMAX/NY,NMAX/NZ),
     *                  RSIZE(NNUC),
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF                       ! Vmol
C     COMMON / RGRID / RR1(NMAX*2,NMAX*2,NMAX*2),
C    *                 RR2(NMAX*2,NMAX*2,NMAX*2)
C     COMMON / FFT / N_I, N_J, N_K, LDA_I, LDA_J


      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C----- Broadcast Parameters -----

      CALL MPI_BCAST(NRST,7,MPI_CHARACTER,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NOPT,3,MPI_CHARACTER,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(EXC,7,MPI_CHARACTER,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NQMMM,4,MPI_CHARACTER,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(FREEZE,5,MPI_CHARACTER,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(PRINT,5,MPI_CHARACTER,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DGF,3,MPI_CHARACTER,0,
     *               MPI_COMM_WORLD,IERR)

      CALL MPI_BCAST(COE,1,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(TEMP,1,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DTMD,1,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NDEN,1,MPI_INTEGER,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(MDMAX,1,MPI_INTEGER,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NRVLC,1,MPI_INTEGER,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NCHK,1,MPI_INTEGER,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(CONV,1,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NPD(1),NNUC,MPI_INTEGER,0,      ! non-local d
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NPC(1),NNUC,MPI_INTEGER,0,      ! partialc charge correction (PCC)
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ZA(1),NNUC,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ZVAL(1),NNUC,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(SIG(1),NNUC,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(EPSQM(1),NNUC,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(RSIZE(1),NNUC,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ZPOP(1),NFUZZY,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NF,1,MPI_INTEGER,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NLINK,1,MPI_INTEGER,0,
     *                MPI_COMM_WORLD,IERR )
      nlqm = NNUC-NLINK                                     ! Vmol
      CALL MPI_BCAST(MMID(1),nlqm,MPI_INTEGER,0,
     *               MPI_COMM_WORLD,IERR)

C----- PARAMETER -----

C     N_I = NMAX*2
C     N_J = NMAX*2
C     N_K = NMAX*2

C     LDA_I = NMAX*2
C     LDA_J = NMAX*2

C----- GRID2 -----

      GMAX = DSQRT( 2.D0*COE )

      DX = PI / ( GMAX )              
      DY = DX
      DZ = DY

      DV = DX*DY*DZ

      XL  = DX * DBLE(NMAXX)
      YL  = DY * DBLE(NMAXY)
      ZL  = DZ * DBLE(NMAXZ)
      
      VOL = XL*YL*ZL

      IF(MYID.EQ.0) THEN
      WRITE(*,*) "Grid Space and Volume [a.u.]" 
      WRITE(*,*) "  dx  = ",DX 
      WRITE(*,*) "  dy  = ",DY 
      WRITE(*,*) "  dz  = ",DZ 
      WRITE(*,*) "  dv  = ",DV 
      WRITE(*,*) 
      WRITE(*,*) "QM Cell Size [a.u.]" 
      WRITE(*,*) "  xl  = ",XL 
      WRITE(*,*) "  yl  = ",YL 
      WRITE(*,*) "  zl  = ",ZL 
      WRITE(*,*) "  vol = ",VOL 
      ENDIF

C----- ELCTRN -----
C
C  4-th Order Finite Difference Expression
C

      A1  =  8.D0/ 5.D0
      A2  = -1.D0/ 5.D0
      A3  =  8.D0/ 315.D0
      A4  = -1.D0/ 560.D0

      AX1 = -A1 / (2.D0*(DX**2))
      AX2 = -A2 / (2.D0*(DX**2))
      AX3 = -A3 / (2.D0*(DX**2))
      AX4 = -A4 / (2.D0*(DX**2))

      AY1 = -A1 / (2.D0*(DY**2))
      AY2 = -A2 / (2.D0*(DY**2))
      AY3 = -A3 / (2.D0*(DY**2))
      AY4 = -A4 / (2.D0*(DY**2))

      AZ1 = -A1 / (2.D0*(DZ**2))
      AZ2 = -A2 / (2.D0*(DZ**2))
      AZ3 = -A3 / (2.D0*(DZ**2))
      AZ4 = -A4 / (2.D0*(DZ**2))

      AO  = -2.D0 * ( AX1+AX2+AX3+AX4 
     *            +   AY1+AY2+AY3+AY4
     *            +   AZ1+AZ2+AZ3+AZ4 )
 
C
C  For Poisson Equation
C

      NF = 0
      DO NA = 1, NNUC
C       IF( ZA(NA) .GT. 1.D0 ) THEN
          NF = NF + 1
          NPOP(NF) = NA
C       ENDIF
      ENDDO

      TWO  = 2.d0
      TAO  = TWO*AO
      TAX1 = TWO*AX1
      TAX2 = TWO*AX2
      TAX3 = TWO*AX3
      TAX4 = TWO*AX4
      TAY1 = TWO*AY1
      TAY2 = TWO*AY2
      TAY3 = TWO*AY3
      TAY4 = TWO*AY4
      TAZ1 = TWO*AZ1
      TAZ2 = TWO*AZ2
      TAZ3 = TWO*AZ3
      TAZ4 = TWO*AZ4
      PI4  = 4.d0*PI

C----- CALCULATE INVERSE R OVER THE GRID -----

C     DO I = 1,NMAX*2 
C     DO J = 1,NMAX*2 
C     DO K = 1,NMAX*2 

C         FRX = DX*( I-NMAX-1 )  
C         FRY = DY*( J-NMAX-1 )  
C         FRZ = DZ*( K-NMAX-1 )  
C          FR  = FRX**2 + FRY**2 + FRZ**2

C     IF( FR .EQ. 0.D0 ) THEN

C         RR1( I,J,K ) = -(DX**2)*(PI/2.D0 
C    *                 +   3.D0*DLOG((3.D0**0.5D0-1.D0)
C    *                 /  (3.D0**0.5D0+1.D0)))
C         RR2( I,J,K ) = 0.D0 

C     ELSE

C         RR1( I,J,K ) = (DX*DY*DZ)/DSQRT(FR)
C         RR2( I,J,K ) = 0.D0 

C     ENDIF

C     ENDDO
C     ENDDO
C     ENDDO

C     STATUS = ZFFT_3D('R','R','F',RR1,RR2,RR1,RR2,
C    *                  N_I,N_J,N_K,LDA_I,LDA_J,1,1,1)

C Initialize Poisson solver
         CALL PSOLVER_INIT((/NMAXX, NMAXY, NMAXZ/),
     *     ICOMM_CART,
     *     (/DX, DY, DZ/),
     *     (/MX, MY, MZ/),
     *     poisson_bc_dirichlet_zero_staggered)
C     *     verbose = .TRUE.)

      RETURN
      END

C----------------------------------
C   SUBROUTINE COMPUTING DENSITY 
C----------------------------------

      SUBROUTINE DNST( RWFAB,RHOAB,MOR ) ! spin density for orbitals

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC = 36 )
C     PARAMETER ( NORA = 49 )

      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ
      COMMON / NCLR  / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV

      DIMENSION RWFAB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION RHOAB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      DV = DX*DY*DZ

C----- initialize -----

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

        RHOAB( I,J,K ) = 0.D0

      ENDDO
      ENDDO
      ENDDO

C----- spin electron density -----

      DO L = 1, MOR
      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

            RHOAB( I,J,K ) = RHOAB( I,J,K ) 
     *                     + RWFAB( I,J,K,L )**2

      ENDDO
      ENDDO
      ENDDO
      ENDDO

C----- normalize -----

      ANORM = 0.D0
      DMOR  = DBLE(MOR)

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

          ANORM = ANORM + RHOAB( I,J,K ) *DV

      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE(ANORM,BNORM,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST( BNORM,1,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

          RHOAB( I,J,K ) = DMOR * RHOAB( I,J,K ) / BNORM  

      ENDDO
      ENDDO
      ENDDO

      RETURN 
      END

      SUBROUTINE RDNST( RHO,RHOA,ZSUM )         ! total density for restricted orbitals

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC =  36 )

      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ
      COMMON / NCLR  / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV

      DIMENSION RHO ( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHOA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      DV = DX*DY*DZ

C----- total electron density -----

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

          RHO( I,J,K ) = 0.5D0*RHO(  I,J,K ) 
     *                 +       RHOA( I,J,K )            

      ENDDO
      ENDDO
      ENDDO

C----- normalize -----

      ANORM  = 0.D0

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

          ANORM = ANORM  + RHO( I,J,K ) *DV

      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE(ANORM,BNORM,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST( BNORM,1,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

          RHO(  I,J,K ) = ZSUM * RHO(  I,J,K ) / BNORM  

      ENDDO
      ENDDO
      ENDDO

      RETURN 
      END

      SUBROUTINE UDNST( RHO,RHOA,RHOB,ZSUM ) ! total density for unrestricted orbitals

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NORA =  49 )

      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV

      DIMENSION RHO ( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHOA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHOB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      DV = DX*DY*DZ

C----- total electron density -----

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

          RHO( I,J,K ) = 0.5D0*( RHO(  I,J,K )
     *                 +         RHOA( I,J,K )            
     *                 +         RHOB( I,J,K ) )          

      ENDDO
      ENDDO
      ENDDO

C----- normalize -----

      ANORM  = 0.D0

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

          ANORM  = ANORM  + RHO(  I,J,K ) *DV

      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE(ANORM,BNORM,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST( BNORM,1,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

          RHO( I,J,K ) = ZSUM * RHO( I,J,K ) / BNORM  

      ENDDO
      ENDDO
      ENDDO

      RETURN 
      END

C------------------------------------------
C   SUBROUTINE DIAGONALIZE WAVEFUNCTIONS
C   BY GRAM-SCHMIDT METHOD
C------------------------------------------

      SUBROUTINE DIAG( NRD,NOR,TWF,RWF ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                            ! mpi
      include 'QMpara.i'                         ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NORA =  49 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NOCC =  1 )

      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ

      DIMENSION NRD(    NOR )
      DIMENSION COEFF1( NOR-1 )
      DIMENSION COEFF2( NOR-1 )

      DIMENSION TWF( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NOR )
      DIMENSION RWF( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NOR )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      DV = DX*DY*DZ
C----- normalize -----

      ANORM = 0.D0

      DO 10 L = 1, JZ
        DO 20 M = 1, JY
          DO 30 N = 1, JX

          ANORM = ANORM + DABS( TWF( N,M,L,NRD(1) ) )**2*DV

30        CONTINUE
20      CONTINUE
10    CONTINUE

      CALL MPI_REDUCE(ANORM,BNORM,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)

      IF(MYID.EQ.0) THEN
      ANORM = DSQRT( BNORM )
      ENDIF

      CALL MPI_BCAST(ANORM,1,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

      DO 40 L = 1, JZ
        DO 50 M = 1, JY
          DO 60 N = 1, JX

          TWF( N,M,L,NRD(1) ) = TWF( N,M,L,NRD(1) ) / ANORM 
C         RWF( N,M,L,NOR )    = TWF( N,M,L,NRD(1) ) 
          RWF( N,M,L,1 )      = TWF( N,M,L,NRD(1) ) 

60        CONTINUE
50      CONTINUE
40    CONTINUE

C----- gram schmidt -----

      DO 300 NTMP = 2, NOR        

        DO 11 K = 1, NTMP-1
        COEFF1(K) = 0.D0
11      CONTINUE

        DO 70 K = 1,NTMP-1        
          DO 80 L = 1, JZ
            DO 90 M = 1, JY
              DO 100 N = 1, JX

            COEFF1(K) = COEFF1(K) + TWF(N,M,L,NRD(K))  
     *                * TWF(N,M,L,NRD(NTMP))*DV

100          CONTINUE
90         CONTINUE
80       CONTINUE
70     CONTINUE

      CALL MPI_REDUCE(COEFF1(1),COEFF2(1),NTMP-1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(COEFF2(1),NTMP-1,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

      DO 110 K = 1,NTMP-1        
        DO 120 L = 1, JZ
          DO 130 M = 1, JY
            DO 140 N = 1, JX

            TWF( N,M,L,NRD(NTMP) ) = TWF( N,M,L,NRD(NTMP) ) 
     *                             - COEFF2(K)*TWF( N,M,L,NRD(K) )

140         CONTINUE
130       CONTINUE
120     CONTINUE
110   CONTINUE

      PNORM = 0.D0

      DO L = 1, JZ
        DO M = 1, JY
          DO N = 1, JX

            PNORM = PNORM + TWF(N,M,L,NRD(NTMP))**2*DV
      
          ENDDO
        ENDDO
      ENDDO

      CALL MPI_REDUCE(PNORM,QNORM,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      IF(MYID.EQ.0) THEN
      PNORM = DSQRT( QNORM ) 
      ENDIF
      CALL MPI_BCAST(PNORM,1,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

      DO 180 L = 1, JZ
        DO 190 M = 1, JY
          DO 200 N = 1, JX

          TWF( N,M,L,NRD( NTMP ) )  = TWF( N,M,L,NRD( NTMP ) ) 
     *                              / PNORM
          RWF( N,M,L,NTMP )         = TWF( N,M,L,NRD( NTMP ) )  
C         RWF( N,M,L,(NOR-NTMP+1) ) = TWF( N,M,L,NRD( NTMP ) )  

200       CONTINUE
190     CONTINUE
180   CONTINUE

300   CONTINUE

      RETURN
      END

C--------------------------------------------------------
C   SUBROUTINE Double Grid
C   Phys. Rev. Lett. 48, 1425( 1982 )
C   Phys. Rev. B.    50, 12234(1994 )
C--------------------------------------------------------

      SUBROUTINE DG(VNUC)                  ! for non-periodic system

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol


C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NDIM =  3 )
      PARAMETER ( RCUT = 3.5D0 )
      PARAMETER ( NSL  = 9000  )

      PARAMETER ( PI    = 3.14159265358979323D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ
      COMMON / ACELL / XL, YL, ZL
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)
      COMMON / DBLG  / WNLOCS( NNUC,NSL ),WNLOCPX( NNUC,NSL ),
     *                 WNLOCPY(NNUC,NSL ),WNLOCPZ( NNUC,NSL )
      COMMON / PLOC / VLOC( NNUC,NSL )
      COMMON / BHS / C1(100),C2(100),AL1(100),AL2(100)

      DIMENSION VNUC( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) 'Double Grid(Linear Interpolation): Start'
      STIME = MPI_WTIME()
      ENDIF

C----- initialize -----

      DO J=1, NSL 
      DO I=1, NNUC  

        WNLOCS( I,J) = 0.D0
        WNLOCPX(I,J) = 0.D0
        WNLOCPY(I,J) = 0.D0
        WNLOCPZ(I,J) = 0.D0

        VLOC(I,J) = 0.D0

      ENDDO 
      ENDDO 

      DO K=1, JZ 
      DO J=1, JY
      DO I=1, JX

        VNUC( I,J,K ) = 0.D0                                       

      ENDDO 
      ENDDO 
      ENDDO 

C----- compute weight factor for each coarse grid -----  

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO 100 NA = 1, NNUC                  ! loop over atoms
        NUMZ = INT( ZA( NA ) )
        RLN  = DLOG( RAD( NUMZ,1 ) )
        NCOUNT = 0
      DO 10 K=1, JZ
        DO 20 J=1, JY
          DO 30 I=1, JX

          RX = DX*( I-NNX ) - PNUC(NA,1) 
          RY = DY*( J-NNY ) - PNUC(NA,2)
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)
          R2 = RX**2 + RY**2 + RZ**2
          R  = DSQRT(R2) 

          IF( R .LT. RCUT ) THEN
          NCOUNT = NCOUNT + 1
           IF( NCOUNT .GT. NSL ) THEN
           WRITE(*,*) ' dg.f : size of dimension too small '
           STOP
           ENDIF

          DO 40 KD = -NDEN,NDEN
          DO 50 JD = -NDEN,NDEN
          DO 60 ID = -NDEN,NDEN

          RXX = ( RX + DX*DBLE(ID)/DBLE(NDEN)) 
          RYY = ( RY + DY*DBLE(JD)/DBLE(NDEN)) 
          RZZ = ( RZ + DZ*DBLE(KD)/DBLE(NDEN)) 
          RD2 = RXX**2 + RYY**2 + RZZ**2
          RD  = DSQRT( RD2 )

          IF( RD .LT. RAD( NUMZ,1 ) ) THEN
          WGT1   =  RAD( NUMZ,1 ) - RD
          WGT2   =  RD
          VNLOCS = VNLS( NUMZ,1 )  
          VNLOCP = WGT2*VNLP( NUMZ,1 )  / ( WGT1 + WGT2 )     
          VLOCD2 = VLDC( NUMZ,1 )  
C         NR = NR + 1
          ELSE
          NRAD   =  INT( ( 0.5D0*DLOG( RD2 ) - RLN )  
     *           / DRLOG(NUMZ) ) + 1 
          WGT1   =  RAD( NUMZ,NRAD+1 ) - RD
          WGT2   = -RAD( NUMZ,NRAD )   + RD
          VNLOCS = ( WGT1*VNLS( NUMZ,NRAD ) 
     *                 +   WGT2*VNLS( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )     
          VNLOCP = ( WGT1*VNLP( NUMZ,NRAD ) 
     *                 +   WGT2*VNLP( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )     
          VLOCD2 = ( WGT1*VLDC( NUMZ,NRAD ) 
     *           +   WGT2*VLDC( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )                
           ENDIF

          WIJ    = ( 1.D0 - DABS(DBLE(ID)/DBLE(NDEN)))
     *           * ( 1.D0 - DABS(DBLE(JD)/DBLE(NDEN)))
     *           * ( 1.D0 - DABS(DBLE(KD)/DBLE(NDEN)))
     *           / DBLE(NDEN**3)

C--------------------------
C      non-local part      
C--------------------------

C----- for s-component -----

          WNLOCS( NA,NCOUNT )  = WNLOCS( NA,NCOUNT ) 
     *                         + WIJ*VNLOCS 

C----- for p-component -----

          IF( DSQRT(RD2) .LT. 1.0D-06 ) THEN

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 

C         NW = NW + 1

           ELSE

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RXX/RD 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RYY/RD
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RZZ/RD

           ENDIF

C----------------------------------------------
C      local part ( BHS + local-d )
C----------------------------------------------

          IF( RD .LT. 1.0D-06 ) THEN

            VT = -ZVAL(NA)*( C1(NUMZ)*2.D0*DSQRT(AL1(NUMZ))                        ! analytical function ( limit zero ) of BHS
     *         +             C2(NUMZ)*2.D0*DSQRT(AL2(NUMZ)) )  
     *         /  DSQRT(PI)           
            
          ELSE

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*RD2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*RD2) ) ) 
     *         /  RD           
            
          ENDIF

            VLOC( NA,NCOUNT ) = VLOC( NA,NCOUNT )
     *                        + WIJ*(VLOCD2+VT)    

 60        CONTINUE
 50        CONTINUE
 40        CONTINUE

            VNUC( I,J,K ) = VNUC( I,J,K ) 
     *                          + VLOC( NA,NCOUNT )

          ELSE

C----------------------------------------------
C      local part ( BHS )
C----------------------------------------------

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*R2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*R2) ) ) 
     *         /  R           
            
            VNUC( I,J,K ) = VNUC( I,J,K ) + VT                       

          ENDIF

 30        CONTINUE
 20      CONTINUE
 10    CONTINUE

 100   CONTINUE

      write(*,*) 'done! rank =',MYID

      CALL MPI_REDUCE(NCOUNT,NCOUNT1,1,MPI_INTEGER,MPI_SUM,0,
     *                MPI_COMM_WORLD,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) '  ncount = ', NCOUNT1
      ETIME = MPI_WTIME()
      write(*,*) 'Elapsed Time (dg4) = ',ETIME-STIME,' scnds'
      ENDIF


C     write(*,*) 'NR = ',NR 
C     write(*,*) 'NW = ',NW 

      RETURN
      END

C--------------------------------------------------------
C
C   SUBROUTINE Double Grid
C   employing fourth order Lagrange interpolation
C
C   Phys. Rev. Lett. 82, 5016( 1999 )
C   
C--------------------------------------------------------

      SUBROUTINE DG4(VNUC)                      ! for non-periodic system

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
      real*4 tim,ta(2)

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NNUC = 36 )
C     PARAMETER ( NDIM =  3 )
      PARAMETER ( RCUT = 3.5D0 )
      PARAMETER ( NSL  = 9000  )

      PARAMETER ( PI    = 3.14159265358979323D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ
      COMMON / ACELL / XL, YL, ZL
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)
      COMMON / DBLG  / WNLOCS( NNUC,NSL ),WNLOCPX( NNUC,NSL ),
     *                 WNLOCPY(NNUC,NSL ),WNLOCPZ( NNUC,NSL )
      COMMON / PLOC / VLOC( NNUC,NSL )
      COMMON / BHS / C1(100),C2(100),AL1(100),AL2(100)

      DIMENSION VNUC( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) 'Double Grid(4th-order Lagrange Interpolation): Start'
      STIME = MPI_WTIME()
      ENDIF

      RT3  = 1.D0/3.D0
      
C----- initialize -----

      DO J=1, NSL 
      DO I=1, NNUC  

        WNLOCS( I,J) = 0.D0
        WNLOCPX(I,J) = 0.D0
        WNLOCPY(I,J) = 0.D0
        WNLOCPZ(I,J) = 0.D0

        VLOC(I,J) = 0.D0

      ENDDO 
      ENDDO 

      DO K=1,JZ  
      DO J=1,JY 
      DO I=1,JX 

        VNUC( I,J,K ) = 0.D0                                       

      ENDDO 
      ENDDO 
      ENDDO 

C----- compute weight factor for each coarse grid -----  

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1
      DO 100 NA = 1, NNUC                  ! loop over atoms
        NUMZ = INT( ZA( NA ) )
        RLN  = DLOG( RAD( NUMZ,1 ) )
        NCOUNT = 0
      DO 10 K=1,JZ
        DO 20 J=1,JY
          DO 30 I=1,JX

          RX = DX*( I-NNX ) - PNUC(NA,1) 
          RY = DY*( J-NNY ) - PNUC(NA,2)
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)
          R2 = RX**2 + RY**2 + RZ**2
          R  = DSQRT(R2) 

          IF( R .LT. RCUT ) THEN
          NCOUNT = NCOUNT + 1
           IF( NCOUNT .GT. NSL ) THEN
           WRITE(*,*) ' dg4.f : size of dimension too small '
           STOP
           ENDIF

          DO 40 KD = -2*NDEN,2*NDEN
          DO 50 JD = -2*NDEN,2*NDEN
          DO 60 ID = -2*NDEN,2*NDEN

          DID = DABS(DBLE(ID)/DBLE(NDEN))         
          DJD = DABS(DBLE(JD)/DBLE(NDEN))         
          DKD = DABS(DBLE(KD)/DBLE(NDEN))         

          RXX = ( RX + DX*DBLE(ID)/DBLE(NDEN)) 
          RYY = ( RY + DY*DBLE(JD)/DBLE(NDEN)) 
          RZZ = ( RZ + DZ*DBLE(KD)/DBLE(NDEN)) 
          RD2 = RXX**2 + RYY**2 + RZZ**2
          RD  = DSQRT( RD2 )

          IF( RD .LT. RAD( NUMZ,1 ) ) THEN
C         WRITE(*,*) 'small rd!'
          WGT1   =  RAD( NUMZ,1 ) -  RD 
          WGT2   =  RD
          VNLOCS = VNLS( NUMZ,1 )  
          VNLOCP = WGT2*VNLP( NUMZ,1 )  / ( WGT1 + WGT2 )     
C         VNLOCP = VNLP( NUMZ,1 ) 
          VLOCD2 = VLDC( NUMZ,1 )  
C         NR = NR + 1
          ELSE
          NRAD   =  INT( ( 0.5D0*DLOG( RD2 ) - RLN )  
     *           / DRLOG(NUMZ) ) + 1 
          WGT1   =  RAD( NUMZ,NRAD+1 ) - RD
          WGT2   = -RAD( NUMZ,NRAD )   + RD
          VNLOCS = ( WGT1*VNLS( NUMZ,NRAD ) 
     *           +   WGT2*VNLS( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )     
          VNLOCP = ( WGT1*VNLP( NUMZ,NRAD ) 
     *           +   WGT2*VNLP( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )     
          VLOCD2 = ( WGT1*VLDC( NUMZ,NRAD ) 
     *           +   WGT2*VLDC( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )                
           ENDIF

          IF( (IABS(ID) .LE. NDEN)   .AND.                        ! 1
     *        (IABS(JD) .LE. NDEN)   .AND.
     *        (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 2
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 3
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 4
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 5
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 6
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 7
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DkD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 8
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ENDIF

C--------------------------
C      non-local part      
C--------------------------

C----- for s-component -----

          WNLOCS( NA,NCOUNT )  = WNLOCS( NA,NCOUNT ) 
     *                         + WIJ*VNLOCS 

C----- for p-component -----

          IF( RD .LT. 1.0D-06 ) THEN

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 

C         NW = NW + 1

          ELSE

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RXX/RD 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RYY/RD
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RZZ/RD

          ENDIF

C----------------------------------------------
C      local part ( BHS + local-d )
C----------------------------------------------

          IF( RD .LT. 1.0D-06 ) THEN

            VT = -ZVAL(NA)*( C1(NUMZ)*2.D0*DSQRT(AL1(NUMZ))                        ! analytical function ( limit zero ) of BHS
     *         +             C2(NUMZ)*2.D0*DSQRT(AL2(NUMZ)) )  
     *         /  DSQRT(PI)           
            
          ELSE

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*RD2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*RD2) ) ) 
     *         /  RD           
            
          ENDIF

            VLOC( NA,NCOUNT ) = VLOC( NA,NCOUNT )
     *                        + WIJ*(VLOCD2+VT)    

60        CONTINUE
50        CONTINUE
40        CONTINUE

            VNUC( I,J,K ) = VNUC( I,J,K ) 
     *                          + VLOC( NA,NCOUNT )

          ELSE

C----------------------------------------------
C      local part ( BHS )
C----------------------------------------------

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*R2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*R2) ) ) 
     *         /  R           
            
            VNUC( I,J,K ) = VNUC( I,J,K ) + VT                       

          ENDIF

30        CONTINUE
20      CONTINUE
10    CONTINUE

100   CONTINUE

      write(*,*) 'done! rank =',MYID

      CALL MPI_REDUCE(NCOUNT,NCOUNT1,1,MPI_INTEGER,MPI_SUM,0,
     *                MPI_COMM_WORLD,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) '  ncount = ', NCOUNT1
      ETIME = MPI_WTIME()
      write(*,*) 'Elapsed Time (dg4) = ',ETIME-STIME,' scnds'
      ENDIF

C     write(*,*) 'NR = ',NR 
C     write(*,*) 'NW = ',NW 

      RETURN
      END

      SUBROUTINE DG4QM                       ! for non-periodic system
                                             ! for QM atoms except for link atoms
      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NDIM =  3 )
      PARAMETER ( RCUT = 3.5D0 )
      PARAMETER ( NSL  = 9000  )

      PARAMETER ( PI    = 3.14159265358979323D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX, DY, DZ
      COMMON / ACELL / XL, YL, ZL
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)
      COMMON / DBLG  / WNLOCS( NNUC,NSL ),WNLOCPX( NNUC,NSL ),
     *                 WNLOCPY(NNUC,NSL ),WNLOCPZ( NNUC,NSL )
C     COMMON / DBLG1 / VNUCQM( NMAX,NMAX,NMAX )
      COMMON / PLOC / VLOC( NNUC,NSL )
      COMMON / BHS / C1(100),C2(100),AL1(100),AL2(100)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) 'Double Grid(4th-order Lagrange Interpolation): Start'
      WRITE(*,*) '--- for QM Atoms ( not include link atoms ) ---'
      ENDIF

      RT3  = 1.D0/3.D0
      ne = NNUC - NLINK
      
C----- initialize -----

      DO J=1, NSL 
      DO I=1, ne  

        WNLOCS( I,J) = 0.D0
        WNLOCPX(I,J) = 0.D0
        WNLOCPY(I,J) = 0.D0
        WNLOCPZ(I,J) = 0.D0

        VLOC(I,J) = 0.D0

      ENDDO 
      ENDDO 

C     DO K=1, NMAX  
C     DO J=1, NMAX 
C     DO I=1, NMAX 

C       VNUCQM( I,J,K ) = 0.D0                                       

C     ENDDO 
C     ENDDO 
C     ENDDO 

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

C----- compute weight factor for each coarse grid -----  

      DO 100 NA = 1, ne                  ! loop over atoms
        NUMZ = INT( ZA( NA ) )
        RLN  = DLOG( RAD( NUMZ,1 ) )
        NCOUNT = 0
      DO 10 K=1, JZ
        DO 20 J=1, JY
          DO 30 I=1, JX

          RX = DX*( I-NNX ) - PNUC(NA,1) 
          RY = DY*( J-NNY ) - PNUC(NA,2)
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)
          R2 = RX**2 + RY**2 + RZ**2
          R  = DSQRT(R2) 

          IF( R .LT. RCUT ) THEN
          NCOUNT = NCOUNT + 1
           IF( NCOUNT .GT. NSL ) THEN
           WRITE(*,*) ' dg4.f : size of dimension too small '
           STOP
           ENDIF

          DO 40 KD = -2*NDEN,2*NDEN
          DO 50 JD = -2*NDEN,2*NDEN
          DO 60 ID = -2*NDEN,2*NDEN

          DID = DABS(DBLE(ID)/DBLE(NDEN))         
          DJD = DABS(DBLE(JD)/DBLE(NDEN))         
          DKD = DABS(DBLE(KD)/DBLE(NDEN))         

          RXX = ( RX + DX*DBLE(ID)/DBLE(NDEN)) 
          RYY = ( RY + DY*DBLE(JD)/DBLE(NDEN)) 
          RZZ = ( RZ + DZ*DBLE(KD)/DBLE(NDEN)) 
          RD2 = RXX**2 + RYY**2 + RZZ**2
          RD  = DSQRT( RD2 )

          IF( RD .LT. RAD( NUMZ,1 ) ) THEN
          WGT1   =  RAD( NUMZ,1 ) -  RD 
          WGT2   =  RD
          VNLOCS = VNLS( NUMZ,1 )  
          VNLOCP = WGT2*VNLP( NUMZ,1 )  / ( WGT1 + WGT2 )     
          VLOCD2 = VLDC( NUMZ,1 )  
C         NR = NR + 1
          ELSE
          NRAD   =  INT( ( 0.5D0*DLOG( RD2 ) - RLN )  
     *           / DRLOG(NUMZ) ) + 1 
          WGT1   =  RAD( NUMZ,NRAD+1 ) - RD
          WGT2   = -RAD( NUMZ,NRAD )   + RD
          VNLOCS = ( WGT1*VNLS( NUMZ,NRAD ) 
     *                 +   WGT2*VNLS( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )     
          VNLOCP = ( WGT1*VNLP( NUMZ,NRAD ) 
     *                 +   WGT2*VNLP( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )     
          VLOCD2 = ( WGT1*VLDC( NUMZ,NRAD ) 
     *           +   WGT2*VLDC( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )                
           ENDIF

          IF( (IABS(ID) .LE. NDEN)   .AND.                        ! 1
     *        (IABS(JD) .LE. NDEN)   .AND.
     *        (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 2
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 3
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 4
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 5
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 6
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 7
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DkD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 8
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ENDIF

C--------------------------
C      non-local part      
C--------------------------

C----- for s-component -----

          WNLOCS( NA,NCOUNT )  = WNLOCS( NA,NCOUNT ) 
     *                         + WIJ*VNLOCS 

C----- for p-component -----

          IF( RD .LT. 1.0D-06 ) THEN

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 

C         NW = NW + 1

           ELSE

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RXX/RD 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RYY/RD
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RZZ/RD

           ENDIF

C----------------------------------------------
C      local part ( BHS + local-d )
C----------------------------------------------

          IF( RD .LT. 1.0D-06 ) THEN

            VT = -ZVAL(NA)*( C1(NUMZ)*2.D0*DSQRT(AL1(NUMZ))                        ! analytical function ( limit zero ) of BHS
     *         +             C2(NUMZ)*2.D0*DSQRT(AL2(NUMZ)) )  
     *         /  DSQRT(PI)           
            
          ELSE

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*RD2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*RD2) ) ) 
     *         /  RD           
            
          ENDIF

            VLOC( NA,NCOUNT ) = VLOC( NA,NCOUNT )
     *                        + WIJ*(VLOCD2+VT)    

60        CONTINUE
50        CONTINUE
40        CONTINUE

C           VNUCQM( I,J,K ) = VNUCQM( I,J,K ) 
C    *                          + VLOC( NA,NCOUNT )

C         ELSE

C----------------------------------------------
C      local part ( BHS )
C----------------------------------------------

C           VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*R2) )                  ! analytical function of BHS
C    *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*R2) ) ) 
C    *         /  R           
C           
C           VNUCQM( I,J,K ) = VNUCQM( I,J,K ) + VT                       

          ENDIF

30        CONTINUE
20      CONTINUE
10    CONTINUE

100   CONTINUE

      write(*,*) 'done! rank =',MYID

      CALL MPI_REDUCE(NCOUNT,NCOUNT1,1,MPI_INTEGER,MPI_SUM,0,
     *                MPI_COMM_WORLD,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) '  ncount = ', NCOUNT1
      ENDIF

C     write(*,*) 'NR = ',NR 
C     write(*,*) 'NW = ',NW 

      RETURN
      END

      SUBROUTINE DG4L(VNUCN,VNUC)            ! for non-periodic system
                                             ! for link atoms
      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NNUC = 36 )
C     PARAMETER ( NDIM =  3 )
      PARAMETER ( RCUT = 3.5D0 )
      PARAMETER ( NSL  = 9000  )

      PARAMETER ( PI    = 3.14159265358979323D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX, DY, DZ
      COMMON / ACELL / XL, YL, ZL
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)
      COMMON / DBLG  / WNLOCS( NNUC,NSL ),WNLOCPX( NNUC,NSL ),
     *                 WNLOCPY(NNUC,NSL ),WNLOCPZ( NNUC,NSL )
C     COMMON / DBLG1 / VNUCQM( NMAX,NMAX,NMAX )
      COMMON / PLOC / VLOC( NNUC,NSL )
      COMMON / BHS / C1(100),C2(100),AL1(100),AL2(100)

      DIMENSION VNUC(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VNUCL( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VNUCN( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) 'Double Grid(4th-order Lagrange Interpolation): Start'
      WRITE(*,*) '--- for Link Atoms ---'
      ENDIF

      RT3 = 1.D0/3.D0
      ns  = NNUC - NLINK + 1
      
C----- initialize -----

      DO J=1, NSL 
      DO I=ns, NNUC  

        WNLOCS( I,J) = 0.D0
        WNLOCPX(I,J) = 0.D0
        WNLOCPY(I,J) = 0.D0
        WNLOCPZ(I,J) = 0.D0

        VLOC(I,J) = 0.D0

      ENDDO 
      ENDDO 

      DO K=1, JZ  
      DO J=1, JY 
      DO I=1, JX 

        VNUCL( I,J,K ) = 0.0D0     ! VNUC for Link atoms 2016.03.18 takahashi
        VNUC(  I,J,K ) = 0.0D0   

      ENDDO 
      ENDDO 
      ENDDO 

      NCOUNT = 0
      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

C----- compute weight factor for each coarse grid -----  

      DO 100 NA = ns, NNUC                  ! loop over atoms
        NUMZ = INT(  ZA(  NA ) )
        RLN  = DLOG( RAD( NUMZ,1 ) )
        NCOUNT = 0
      DO 10 K=1, JZ
        DO 20 J=1, JY
          DO 30 I=1, JX

          RX = DX*( I-NNX ) - PNUC(NA,1) 
          RY = DY*( J-NNY ) - PNUC(NA,2)
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)
          R2 = RX**2 + RY**2 + RZ**2
          R  = DSQRT(R2) 

          IF( R .LT. RCUT ) THEN
          NCOUNT = NCOUNT + 1
           IF( NCOUNT .GT. NSL ) THEN
           WRITE(*,*) ' dg4.f : size of dimension too small '
           STOP
           ENDIF

          DO 40 KD = -2*NDEN,2*NDEN
          DO 50 JD = -2*NDEN,2*NDEN
          DO 60 ID = -2*NDEN,2*NDEN

          DID = DABS(DBLE(ID)/DBLE(NDEN))         
          DJD = DABS(DBLE(JD)/DBLE(NDEN))         
          DKD = DABS(DBLE(KD)/DBLE(NDEN))         

          RXX = ( RX + DX*DBLE(ID)/DBLE(NDEN)) 
          RYY = ( RY + DY*DBLE(JD)/DBLE(NDEN)) 
          RZZ = ( RZ + DZ*DBLE(KD)/DBLE(NDEN)) 
          RD2 = RXX**2 + RYY**2 + RZZ**2
          RD  = DSQRT( RD2 )

          IF( RD .LT. RAD( NUMZ,1 ) ) THEN
          WGT1   =  RAD( NUMZ,1 ) -  RD 
          WGT2   =  RD
          VNLOCS = VNLS( NUMZ,1 )  
          VNLOCP = WGT2*VNLP( NUMZ,1 )  / ( WGT1 + WGT2 )     
          VLOCD2 = VLDC( NUMZ,1 )  
C         NR = NR + 1
          ELSE
          NRAD   =  INT( ( 0.5D0*DLOG( RD2 ) - RLN )  
     *           / DRLOG(NUMZ) ) + 1 
          WGT1   =  RAD( NUMZ,NRAD+1 ) - RD
          WGT2   = -RAD( NUMZ,NRAD )   + RD
          VNLOCS = ( WGT1*VNLS( NUMZ,NRAD ) 
     *                 +   WGT2*VNLS( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )     
          VNLOCP = ( WGT1*VNLP( NUMZ,NRAD ) 
     *                 +   WGT2*VNLP( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )     
          VLOCD2 = ( WGT1*VLDC( NUMZ,NRAD ) 
     *           +   WGT2*VLDC( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )                
           ENDIF

          IF( (IABS(ID) .LE. NDEN)   .AND.                        ! 1
     *        (IABS(JD) .LE. NDEN)   .AND.
     *        (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 2
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 3
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 4
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 5
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 6
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 7
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DkD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 8
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ENDIF

C--------------------------
C      non-local part      
C--------------------------

C----- for s-component -----

          WNLOCS( NA,NCOUNT )  = WNLOCS( NA,NCOUNT ) 
     *                         + WIJ*VNLOCS 

C----- for p-component -----

          IF( RD .LT. 1.0D-06 ) THEN

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 

C         NW = NW + 1

           ELSE

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RXX/RD 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RYY/RD
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RZZ/RD

           ENDIF

C----------------------------------------------
C      local part ( BHS + local-d )
C----------------------------------------------

          IF( RD .LT. 1.0D-06 ) THEN

            VT = -ZVAL(NA)*( C1(NUMZ)*2.D0*DSQRT(AL1(NUMZ))                        ! analytical function ( limit zero ) of BHS
     *         +             C2(NUMZ)*2.D0*DSQRT(AL2(NUMZ)) )  
     *         /  DSQRT(PI)           
            
          ELSE

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*RD2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*RD2) ) ) 
     *         /  RD           
            
          ENDIF

            VLOC( NA,NCOUNT ) = VLOC( NA,NCOUNT )
     *                        + WIJ*(VLOCD2+VT)    

60        CONTINUE
50        CONTINUE
40        CONTINUE

            VNUCL( I,J,K ) = VNUCL( I,J,K ) 
     *                          + VLOC( NA,NCOUNT )

          ELSE

C----------------------------------------------
C      local part ( BHS )
C----------------------------------------------

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*R2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*R2) ) ) 
     *         /  R           
            
            VNUCL( I,J,K ) = VNUCL( I,J,K ) + VT                       

          ENDIF

30        CONTINUE
20      CONTINUE
10    CONTINUE

100   CONTINUE

      write(*,*) 'done! rank (Link) =',MYID

      CALL MPI_REDUCE(NCOUNT,NCOUNT1,1,MPI_INTEGER,MPI_SUM,0,
     *                MPI_COMM_WORLD,IERR)

      IF(MYID.EQ.0) THEN
        WRITE(*,*) '  ncount = ', NCOUNT1
      ENDIF

      VNUC(:,:,:) = VNUCN(:,:,:) + VNUCL(:,:,:)      ! add vnucl due to Link atoms 2016.03.18 takahashi

C     write(*,*) 'NR = ',NR 
C     write(*,*) 'NW = ',NW 

      ne = NNUC - NLINK
      DO 105 NA = 1, ne                  ! loop over atoms
        NUMZ = INT( ZA( NA ) )
        NCOUNT = 0
      DO 15 K=1, JZ
        DO 25 J=1, JY
          DO 35 I=1, JX

          RX = DX*( I-NNX ) - PNUC(NA,1) 
          RY = DY*( J-NNY ) - PNUC(NA,2)
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)
          R2 = RX**2 + RY**2 + RZ**2
          R  = DSQRT(R2) 

          IF( R .LT. RCUT ) THEN
            NCOUNT = NCOUNT + 1
C           VNUC( I,J,K ) = VNUC( I,J,K ) 
C    *                          + VLOC( NA,NCOUNT )

          ELSE

C----------------------------------------------
C      local part ( BHS )
C----------------------------------------------

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*R2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*R2) ) ) 
     *         /  R           
            
C           VNUC( I,J,K ) = VNUC( I,J,K ) + VT                       

          ENDIF

35        CONTINUE
25      CONTINUE
15    CONTINUE

105   CONTINUE

C     REWIND(24)
      DO 65 K = 1,JZ  
      DO 55 J = 1,JY
      DO 45 I = 1,JX
C       WRITE(24) VNUCL(I,J,K)
45    CONTINUE
55    CONTINUE
65    CONTINUE

      RETURN
      END

C--------------------------------------------------------
C
C   SUBROUTINE Double Grid
C   employing sixth order Lagrange interpolation
C
C   Phys. Rev. Lett. 82, 5016( 1999 )
C   
C--------------------------------------------------------

      SUBROUTINE DG6(VNUC)                  ! for non-periodic system

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NNUC = 36 )
C     PARAMETER ( NDIM =  3 )
      PARAMETER ( RCUT = 3.5D0 )
      PARAMETER ( NSL  = 9000  )

      PARAMETER ( PI    = 3.14159265358979323D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ
      COMMON / ACELL / XL, YL, ZL
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)
      COMMON / DBLG  / WNLOCS( NNUC,NSL ),WNLOCPX( NNUC,NSL ),
     *                 WNLOCPY(NNUC,NSL ),WNLOCPZ( NNUC,NSL )
      COMMON / PLOC / VLOC( NNUC,NSL )
      COMMON / BHS / C1(100),C2(100),AL1(100),AL2(100)

      DIMENSION VNUC( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) 'Double Grid(6th-order Lagrange Interpolation): Start'
      STIME = MPI_WTIME()
      ENDIF

C----- initialize -----

      DO J=1, NSL 
      DO I=1, NNUC  

        WNLOCS( I,J) = 0.D0
        WNLOCPX(I,J) = 0.D0
        WNLOCPY(I,J) = 0.D0
        WNLOCPZ(I,J) = 0.D0

        VLOC(I,J) = 0.D0

      ENDDO 
      ENDDO 

      DO K=1, NMAXZ 
      DO J=1, NMAXY
      DO I=1, NMAXX 

        VNUC( I,J,K ) = 0.D0                                       

      ENDDO 
      ENDDO 
      ENDDO 

C----- compute weight factor for each coarse grid -----  

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1
      DO 100 NA = 1, NNUC                  ! loop over atoms
        NUMZ = INT( ZA( NA ) )
        RLN  = DLOG( RAD( NUMZ,1 ) )
        NCOUNT = 0
      DO 10 K=1, JZ
        DO 20 J=1, JY
          DO 30 I=1, JX

          RX = DX*( I-NNX ) - PNUC(NA,1) 
          RY = DY*( J-NNY ) - PNUC(NA,2)
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)
          R2 = RX**2 + RY**2 + RZ**2
          R  = DSQRT(R2) 

          IF( R .LT. RCUT ) THEN
          NCOUNT = NCOUNT + 1
           IF( NCOUNT .GT. NSL ) THEN
           WRITE(*,*) ' dg6.f : size of dimension too small '
           STOP
           ENDIF

          DO 40 KD = -NDEN,NDEN-1
          DO 50 JD = -NDEN,NDEN-1
          DO 60 ID = -NDEN,NDEN-1

          DID = DBLE(ID)/DBLE(NDEN)         
          DJD = DBLE(JD)/DBLE(NDEN)         
          DKD = DBLE(KD)/DBLE(NDEN)         

          DO 45 KE = 1,3
          DO 55 JE = 1,3
          DO 65 IE = 1,3
          
          IF(ID.LT.0) THEN
          RXX = ( RX + DX*(DID-DBLE(IE)+1.D0) )
          ELSE
          RXX = ( RX + DX*(DID+DBLE(IE)-1.D0) )
          ENDIF

          IF(JD.LT.0) THEN
          RYY = ( RY + DY*(DJD-DBLE(JE)+1.D0) ) 
          ELSE
          RYY = ( RY + DY*(DJD+DBLE(JE)-1.D0) ) 
          ENDIF

          IF(KD.LT.0) THEN
          RZZ = ( RZ + DZ*(DKD-DBLE(KE)+1.D0) ) 
          ELSE
          RZZ = ( RZ + DZ*(DKD+DBLE(KE)-1.D0) ) 
          ENDIF

          RD2 = RXX**2 + RYY**2 + RZZ**2
          RD  = DSQRT( RD2 )

          IF( RD .LT. RAD( NUMZ,1 ) ) THEN
          WGT1   =  RAD( NUMZ,1 ) -  RD 
          WGT2   =  RD
          VNLOCS = VNLS( NUMZ,1 )  
          VNLOCP = WGT2*VNLP( NUMZ,1 )  / ( WGT1 + WGT2 )     
          VLOCD2 = VLDC( NUMZ,1 )  
C         NR = NR + 1
          ELSE
          NRAD   =  INT( ( 0.5D0*DLOG( RD2 ) - RLN )  
     *           / DRLOG(NUMZ) ) + 1 
          WGT1   =  RAD( NUMZ,NRAD+1 ) - RD
          WGT2   = -RAD( NUMZ,NRAD )   + RD
          VNLOCS = ( WGT1*VNLS( NUMZ,NRAD ) 
     *                 +   WGT2*VNLS( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )     
          VNLOCP = ( WGT1*VNLP( NUMZ,NRAD ) 
     *                 +   WGT2*VNLP( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )     
          VLOCD2 = ( WGT1*VLDC( NUMZ,NRAD ) 
     *           +   WGT2*VLDC( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )                
           ENDIF

          IF(ID.LT.0) THEN

          IF(IE.EQ.1) THEN
          WIJ = DBLE((3*NDEN-IABS(ID))*(2*NDEN-IABS(ID))
     *        *      (NDEN-IABS(ID))*(IABS(ID)+NDEN)
     *        *      (IABS(ID)+2*NDEN))
     *        / DBLE(12*(NDEN**6)) 
          ELSEIF(IE.EQ.2) THEN
          WIJ = -DBLE((3*NDEN-IABS(ID))*(2*NDEN-IABS(ID))
     *        *       (NDEN-IABS(ID))*IABS(ID)
     *        *       (IABS(ID)+2*NDEN))
     *        /  DBLE(24*(NDEN**6)) 
          ELSEIF(IE.EQ.3) THEN
          WIJ = DBLE((3*NDEN-IABS(ID))*(2*NDEN-IABS(ID))
     *        *      (NDEN-IABS(ID))*IABS(ID)
     *        *      (IABS(ID)+NDEN))
     *        / DBLE(120*(NDEN**6)) 
          ENDIF

          ELSE

          IF(IE.EQ.1) THEN
          WIJ = -DBLE((2*NDEN+IABS(ID))*(NDEN+IABS(ID))
     *        *       (IABS(ID)-NDEN)*(IABS(ID)-2*NDEN)
     *        *       (IABS(ID)-3*NDEN))
     *        /  DBLE(12*(NDEN**6)) 
          ELSEIF(IE.EQ.2) THEN
          WIJ = DBLE((2*NDEN+IABS(ID))*IABS(ID)*(IABS(ID)-NDEN)
     *        *      (IABS(ID)-2*NDEN)*(IABS(ID)-3*NDEN))
     *        / DBLE(24*(NDEN**6)) 
          ELSEIF(IE.EQ.3) THEN
          WIJ = -DBLE((NDEN+IABS(ID))*IABS(ID)
     *        *       (IABS(ID)-NDEN)*(IABS(ID)-2*NDEN)
     *        *       (IABS(ID)-3*NDEN))
     *        /  DBLE(120*(NDEN**6)) 
          ENDIF

          ENDIF

          IF(JD.LT.0) THEN

          IF(JE.EQ.1) THEN
          WIJ =  WIJ
     *        * ( DBLE((3*NDEN-IABS(JD))*(2*NDEN-IABS(JD))
     *        *        (NDEN-IABS(JD))*(IABS(JD)+NDEN)
     *        *        (IABS(JD)+2*NDEN))
     *        /   DBLE(12*(NDEN**6))) 
          ELSEIF(JE.EQ.2) THEN
          WIJ =  WIJ
     *        * (-DBLE((3*NDEN-IABS(JD))*(2*NDEN-IABS(JD))
     *        *        (NDEN-IABS(JD))*IABS(JD)
     *        *        (IABS(JD)+2*NDEN))
     *        /   DBLE(24*(NDEN**6))) 
          ELSEIF(JE.EQ.3) THEN
          WIJ =  WIJ
     *        * ( DBLE((3*NDEN-IABS(JD))*(2*NDEN-IABS(JD))
     *        *        (NDEN-IABS(JD))*IABS(JD)
     *        *        (IABS(JD)+NDEN))
     *        /   DBLE(120*(NDEN**6))) 
          ENDIF

          ELSE

          IF(JE.EQ.1) THEN
          WIJ =  WIJ
     *        * (-DBLE((2*NDEN+IABS(JD))*(NDEN+IABS(JD))
     *        *        (IABS(JD)-NDEN)*(IABS(JD)-2*NDEN)
     *        *        (IABS(JD)-3*NDEN))
     *        /   DBLE(12*(NDEN**6))) 
          ELSEIF(JE.EQ.2) THEN
          WIJ =  WIJ
     *        * ( DBLE((2*NDEN+IABS(JD))*IABS(JD)*(IABS(JD)-NDEN)
     *        *        (IABS(JD)-2*NDEN)*(IABS(JD)-3*NDEN))
     *        /   DBLE(24*(NDEN**6))) 
          ELSEIF(JE.EQ.3) THEN
          WIJ =  WIJ
     *        * (-DBLE((NDEN+IABS(JD))*IABS(JD)
     *        *        (IABS(JD)-NDEN)*(IABS(JD)-2*NDEN)
     *        *        (IABS(JD)-3*NDEN))
     *        /   DBLE(120*(NDEN**6))) 
          ENDIF

          ENDIF

          IF(KD.LT.0) THEN

          IF(KE.EQ.1) THEN
          WIJ =  WIJ
     *        * ( DBLE((3*NDEN-IABS(KD))*(2*NDEN-IABS(KD))
     *        *        (NDEN-IABS(KD))*(IABS(KD)+NDEN)
     *        *        (IABS(KD)+2*NDEN))
     *        /   DBLE(12*(NDEN**6))) 
          ELSEIF(KE.EQ.2) THEN
          WIJ =  WIJ
     *        * (-DBLE((3*NDEN-IABS(KD))*(2*NDEN-IABS(KD))
     *        *        (NDEN-IABS(KD))*IABS(KD)
     *        *        (IABS(KD)+2*NDEN))
     *        /   DBLE(24*(NDEN**6))) 
          ELSEIF(KE.EQ.3) THEN
          WIJ =  WIJ
     *        * ( DBLE((3*NDEN-IABS(KD))*(2*NDEN-IABS(KD))
     *        *        (NDEN-IABS(KD))*IABS(KD)
     *        *        (IABS(KD)+NDEN))
     *        /   DBLE(120*(NDEN**6))) 
          ENDIF

          ELSE

          IF(KE.EQ.1) THEN
          WIJ =  WIJ
     *        * (-DBLE((2*NDEN+IABS(KD))*(NDEN+IABS(KD))
     *        *        (IABS(KD)-NDEN)*(IABS(KD)-2*NDEN)
     *        *        (IABS(KD)-3*NDEN))
     *        /   DBLE(12*(NDEN**6))) 
          ELSEIF(KE.EQ.2) THEN
          WIJ =  WIJ
     *        * ( DBLE((2*NDEN+IABS(KD))*IABS(KD)*(IABS(KD)-NDEN)
     *        *        (IABS(KD)-2*NDEN)*(IABS(KD)-3*NDEN))
     *        /   DBLE(24*(NDEN**6))) 
          ELSEIF(KE.EQ.3) THEN
          WIJ =  WIJ
     *        * (-DBLE((NDEN+IABS(KD))*IABS(KD)
     *        *        (IABS(KD)-NDEN)*(IABS(KD)-2*NDEN)
     *        *        (IABS(KD)-3*NDEN))
     *        /   DBLE(120*(NDEN**6))) 
          ENDIF

          ENDIF

C--------------------------
C      non-local part      
C--------------------------

C----- for s-component -----

          WNLOCS( NA,NCOUNT )  = WNLOCS( NA,NCOUNT ) 
     *                         + WIJ*VNLOCS 

C----- for p-component -----

          IF( RD .LT. 1.0D-06 ) THEN

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 

C         NW = NW + 1

           ELSE

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RXX/RD 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RYY/RD
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RZZ/RD

           ENDIF

C----------------------------------------------
C      local part ( BHS + local-d )
C----------------------------------------------

          IF( RD .LT. 1.0D-06 ) THEN

            VT = -ZVAL(NA)*( C1(NUMZ)*2.D0*DSQRT(AL1(NUMZ))                        ! analytical function ( limit zero ) of BHS
     *         +             C2(NUMZ)*2.D0*DSQRT(AL2(NUMZ)) )  
     *         /  DSQRT(PI)           
            
          ELSE

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*RD2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*RD2) ) ) 
     *         /  RD           
            
          ENDIF

            VLOC( NA,NCOUNT ) = VLOC( NA,NCOUNT )
     *                        + WIJ*(VLOCD2+VT)    

65        CONTINUE
55        CONTINUE
45        CONTINUE

60        CONTINUE
50        CONTINUE
40        CONTINUE

            VNUC( I,J,K ) = VNUC( I,J,K ) 
     *                          + VLOC( NA,NCOUNT )

          ELSE

C----------------------------------------------
C      local part ( BHS )
C----------------------------------------------

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*R2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*R2) ) ) 
     *         /  R           
            
            VNUC( I,J,K ) = VNUC( I,J,K ) + VT                       

          ENDIF

30        CONTINUE
20      CONTINUE
10    CONTINUE

100   CONTINUE

      write(*,*) 'done! rank =',MYID

      CALL MPI_REDUCE(NCOUNT,NCOUNT1,1,MPI_INTEGER,MPI_SUM,0,
     *                MPI_COMM_WORLD,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) '  ncount = ', NCOUNT1
      ETIME = MPI_WTIME()
      write(*,*) 'Elapsed Time (dg6) = ',ETIME-STIME,' scnds'
      ENDIF

      RETURN
      END

C------------------------------------
C   SUBROUTINE SIMULATED ANNEALING 
C------------------------------------

      SUBROUTINE SCF( RWFA,RWFB,RWFNA,RWFNB,VNUC,PENG,IMD ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'
      include 'QMpara.i'                     ! Vmol
!     include 'sizes.i'                      ! Vmol
      include 'BHHpara.i'                    ! Vmol

      PARAMETER ( CONV_OEP = 1.0D-5 )        ! convergence criterion for OEP potential

C     PARAMETER ( NMAX =   80 )
C     PARAMETER ( NDIM =    3 )
C     PARAMETER ( NORA =   49 )
C     PARAMETER ( NORB =    1 )
C     PARAMETER ( MORA =   49 )
C     PARAMETER ( MORB =    1 )
C     PARAMETER ( NNUC =   36 )
C     PARAMETER ( RCUT = 2.5D0)
C     PARAMETER ( NSL  =  3000)
      PARAMETER ( NQM  =    1 )
C     PARAMETER ( NLINK1 = 10 )              ! Vmol
   
C     PARAMETER( ITMAX  = 2000  )
      PARAMETER( EPS1   =  1.D0 )

      PARAMETER ( PI   =  3.14159265358979323D0 )

      PARAMETER ( NSP = 255 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),                      ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / FRCE  / FRC( NNUC,NDIM )
      COMMON / GRID2 / DX, DY, DZ
      COMMON / CLMB  / PCLMB, PEXC, PEXT 
C     COMMON / LATTICE / POT
      COMMON / LTC     / POT                                       ! Vmol
      COMMON / SLT   / VPCE( NMAXX/NX,NMAXY/NY,NMAXZ/NZ),VPCZ,PLJ
C     COMMON / MMDAT / SCRD( NDIM,4*NSP ),SCRD1( NDIM,4*NSP ),
C    *                 FRCS( 4*NSP,NDIM )
      COMMON / MMDAT / SCRD( NDIM,maxatm ),SCRD1( NDIM,maxatm ),   ! Vmol
     *                 FRCS( maxatm,NDIM )                         ! Vmol
      COMMON / OFC   / VCOU( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      DIMENSION ORBEA(   NORA )
      DIMENSION ORBEB(   NORA )
      DIMENSION NRDA(    NORA )
      DIMENSION NRDB(    NORA )
      DIMENSION VNUC(    NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEXA(    NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEXB(    NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEFFA(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEFFB(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHO(     NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION TRHO(    NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHOA(    NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHOB(    NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
C     DIMENSION VHIJA(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA,NORA)
C     DIMENSION VHIJB(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA,NORA)
      DIMENSION VHIJA(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA)
      DIMENSION VHIJB(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA)
      DIMENSION VCOU_FA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )        ! Fermi-Amaldi potential as reference potential 20230623

      DIMENSION GNALMA(  NNUC,2,3,NORA )
      DIMENSION GNALMB(  NNUC,2,3,NORA )
      DIMENSION DGNALMA( NNUC,5,NORA )
      DIMENSION DGNALMB( NNUC,5,NORA )

      DIMENSION RWFA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION RWFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION RWFX(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
C     DIMENSION RWFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,1 )
      DIMENSION RWFNA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION RWFNB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )

      DIMENSION A1( NORA )
      DIMENSION A2( NORA )
      DIMENSION B1( NORA )
      DIMENSION B2( NORA )

      DIMENSION FRCN( NNUC,NDIM )

      CHARACTER NO1*1,NO2*1,KAZU*2,KAZU2*8

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      DV = DX*DY*DZ

      APOT = POT                                          ! nuclear - nuclear potential

      IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
     *                 .OR.EXC.EQ.'RXalpha') THEN         ! total charge of electrons
        ZSUM = DBLE(2*MORA)
      ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
     *                     .OR.EXC.EQ.'UXalpha') THEN   
        ZSUM = DBLE(MORA+MORB)
      ENDIF

      IF(NQMMM.EQ.'QMMM') THEN
        IF( MOD(IMD,NQM) .NE. 0 ) GOTO 777
!       CALL POCH ( SCRD )
        APOT = APOT + VPCZ + PLJ                          ! add nuclear - site and LJ potential
      ELSEIF(NQMMM.EQ.'LINK') THEN                        ! Vmol
C       CALL POCH1( SCRD  )                               ! Vmol
!       CALL POCH1_S( SCRD1 )                             ! Vmol
        APOT = APOT + VPCZ + PLJ                          ! Vmol
        IF( MOD(IMD,NQM) .NE. 0 ) GOTO 777
C       WRITE(*,*) 'rank',myid,'APOT=',APOT
      ELSEIF(NQMMM.EQ.'MM') THEN                          ! Vmol
C       CALL POCH2( SCRD )                                ! Vmol
!       CALL POCH2( SCRD1 )                               ! Vmol
C       CALL QMF2( SCRD, FRCS,FRCN )                      ! Vmol
!       CALL QMF2( SCRD1,FRCS,FRCN )                      ! Vmol
                                                          ! Vmol
        DO J=1,NDIM                                       ! Vmol
        DO I=1,NNUC                                       ! Vmol
          FRC(I,J) = FRC(I,J) + FRCN(I,J)                 ! Vmol
        ENDDO                                             ! Vmol
        ENDDO                                             ! Vmol
                                                          ! Vmol
        RETURN                                            ! Vmol

      ELSEIF( NQMMM.EQ.'QMVD' .OR. NQMMM.EQ.'QMVC' ) THEN ! Vmol
        CALL RVPCE( IMD )                                 ! Vmol
                                                          ! Vmol
      ENDIF

C----- initialize -----

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        RHO(  I,J,K ) = 0.D0
        TRHO( I,J,K ) = 0.D0
      ENDDO
      ENDDO
      ENDDO

      NJ   = 0

C----------------------------

      CALL OPTFC1(IMD)

      CALL LRCLMB0

C------ start SCF loop ------

      IF(MYID.EQ.0) THEN
      WRITE(*,*)                     
C     WRITE(*,994)                     
      ENDIF

      TOTE = 0.D0
      VHIJA( :,:,:,: ) = 0.D0       ! initialize the exchange pot. 20180309
      VHIJB( :,:,:,: ) = 0.D0       ! initialize the exchange pot. 20180309

      IOEP = 0    ! 20230623 initialize step number of OEP SCF

10    CONTINUE
      STIME = MPI_WTIME()

      NJ   = NJ + 1
      PREE = TOTE

C----- For Restricted orbitals -----

      IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
     *                 .OR.EXC.EQ.'RXalpha') THEN     

C      compute density       

        CALL DNST(  RWFA,RHOA,MORA )          ! spin density for orbitals
C       WRITE(*,*) 'rank ',MYID,'RHOA =',RHOA(1,1,1)
        CALL RDNST(  RHO,RHOA,ZSUM )          ! total density for restricted orbitals
C       WRITE(*,*) 'rank ',MYID,'RHO =',RHO(1,1,1)

C      coulomb potential

C       CALL VCLMB(  RHO )
        TIME1 = MPI_WTIME()
        CALL RVCLMB( RHO )
        TIME2 = MPI_WTIME()
C       IF(MYID.EQ.0) THEN
C         write(*,*) 'elapsed time:RVCLMB =',TIME2-TIME1
C       ENDIF
    
C      exchange and correlation      

        IF(EXC.EQ.'RXalpha') THEN
        
          CALL VXALP(RHO,RHOA,RHOB,VEXA,VEXB,PEXC,PEXT)
        
        ELSEIF(EXC.EQ.'RPZ') THEN
        
          PEXC = 0.D0
          CALL PZ( RHOA,VEXA,PEXC )           ! Perdew & Zunger for spin density
          PEXC = 2.D0*PEXC
      
          PEXT = 0.D0                         ! external exchange and correlation
          DO L = 1, JZ                
          DO M = 1, JY
          DO N = 1, JX
            PEXT = PEXT + RHOA(N,M,L)*VEXA(N,M,L)*DV
          ENDDO
          ENDDO
          ENDDO

          CALL MPI_REDUCE(PEXT,PEXT1,1,MPI_DOUBLE_PRECISION,
     *                  MPI_SUM,0,MPI_COMM_WORLD,IERR)
          PEXT = 2.D0*PEXT1

          CALL MPI_BCAST(PEXT,1,MPI_DOUBLE_PRECISION,
     *                 0,MPI_COMM_WORLD,IERR)

        ELSEIF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' ) THEN
      
C         TIME1 = MPI_WTIME()
!         CALL BLYP(RHO,RHOA,RHOB,VEXA,VEXB,PEXC,PEXT)    ! comment for pure HF calculation 20230619
C         TIME2 = MPI_WTIME()
C         IF(MYID.EQ.0) THEN
C           write(*,*) 'elapsed time:BLYP =',TIME2-TIME1
C         ENDIF

        ENDIF

C15    CONTINUE     ! loop for OEP-SCF   20230623

C      IF ( IOEP > 0 ) THEN                    ! increment step number for OEP 20230623
C        CALL DNST(  RWFA,RHOA,MORA )          ! spin density for orbitals
C        CALL RDNST(  RHO,RHOA,ZSUM )          ! total density for restricted orbitals
C      ENDIF

C      estimate effective potential (local potential)

        DO L = 1, JZ
        DO M = 1, JY
        DO N = 1, JX

          VEFFA(N,M,L) = VNUC(N,M,L)
     *                 + VCOU(N,M,L)
     *                 + VEXA(N,M,L) 

        ENDDO
        ENDDO
        ENDDO

C       WRITE(*,*) 'rank ',MYID,'VNUC =',VNUC(1,1,1)

C       IF(NQMMM.EQ.'QMMM') THEN
        IF(NQMMM.EQ.'QMMM' .OR. NQMMM.EQ.'LINK'  .OR.     ! Vmol
     *     NQMMM.EQ.'QMVD' .OR. NQMMM.EQ.'QMVC') THEN     ! Vmol

          DO L = 1, JZ
          DO M = 1, JY
          DO N = 1, JX

            VEFFA(N,M,L) = VEFFA(N,M,L)
     *                   + VPCE( N,M,L)                   ! add point charge contribution

          ENDDO
          ENDDO
          ENDDO
        ENDIF

C       WRITE(*,*) 'rank ',MYID,'VPCE =',VEFFA(1,1,1)

        DERHO = 0.D0

        DO L = 1, JZ
        DO M = 1, JY
        DO N = 1, JX
C         FRHO  = 2.D0*RHOA( N,M,L )
          GRHO  = RHO( N,M,L ) - TRHO( N,M,L ) 
          DERHO = DERHO + GRHO**2*DV
          TRHO( N,M,L ) = RHO( N,M,L )
        ENDDO
        ENDDO
        ENDDO

        CALL MPI_REDUCE(DERHO,DERHO1,1,MPI_DOUBLE_PRECISION,
     *                  MPI_SUM,0,MPI_COMM_WORLD,IERR)
        DERHO = DERHO1
        CALL MPI_BCAST(DERHO,1,MPI_DOUBLE_PRECISION,
     *                 0,MPI_COMM_WORLD,IERR)

C     IF(MYID.EQ.0) THEN
C     write(6,*) 'derho = ',derho 
C     ENDIF

C----- For Unrestricted orbitals -----

      ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
     *                 .OR.EXC.EQ.'UXalpha') THEN     

C      compute density       

        CALL DNST(  RWFA,RHOA,MORA )          ! spin density for orbitals
        CALL DNST(  RWFB,RHOB,MORB )          ! spin density for orbitals
        CALL UDNST( RHO,RHOA,RHOB,ZSUM )      ! total density for unrestricted orbitals

C      coulomb potential

C       CALL VCLMB( RHO )
        CALL RVCLMB( RHO )

C      exchange and correlation      
    
        IF(EXC.EQ.'UXalpha') THEN
        
          CALL VXALP(RHO,RHOA,RHOB,VEXA,VEXB,PEXC,PEXT)
        
        ELSEIF(EXC.EQ.'UPZ') THEN
        
        PEXC = 0.D0
        CALL PZ( RHOA,VEXA,PEXC ) ! Perdew & Zunger for spin density
        CALL PZ( RHOB,VEXB,PEXC ) ! Perdew & Zunger for spin density
      
        PEXT = 0.D0                                       ! external exchange and correlation
        DO L = 1, JZ                
        DO M = 1, JY
        DO N = 1, JX
          PEXT = PEXT + ( RHOA(N,M,L)*VEXA(N,M,L)
     *         +          RHOB(N,M,L)*VEXB(N,M,L) )*DV
        ENDDO
        ENDDO
        ENDDO
        CALL MPI_REDUCE(PEXT,PEXT1,1,MPI_DOUBLE_PRECISION,
     *                  MPI_SUM,0,MPI_COMM_WORLD,IERR)
        PEXT = PEXT1

        CALL MPI_BCAST(PEXT,1,MPI_DOUBLE_PRECISION,
     *                 0,MPI_COMM_WORLD,IERR)

        ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN
      
C         CALL BLYP(RHO,RHOA,RHOB,VEXA,VEXB,PEXC,PEXT)

        ENDIF

C      estimate effective potential (local potential)

        DO L = 1, JZ
        DO M = 1, JY
        DO N = 1, JX

          VEFFA(N,M,L) = VNUC(N,M,L)
     *                 + VCOU(N,M,L)
     *                 + VEXA(N,M,L) 
          VEFFB(N,M,L) = VNUC(N,M,L)
     *                 + VCOU(N,M,L)
     *                 + VEXB(N,M,L) 

        ENDDO
        ENDDO
        ENDDO

C       IF(NQMMM.EQ.'QMMM') THEN
        IF(NQMMM.EQ.'QMMM' .OR. NQMMM.EQ.'LINK') THEN     ! Vmol
          DO L = 1, JZ
          DO M = 1, JY
          DO N = 1, JX

            VEFFA(N,M,L) = VEFFA(N,M,L)
     *                   + VPCE(N,M,L)                    ! add point charge contribution
            VEFFB(N,M,L) = VEFFB(N,M,L)
     *                   + VPCE(N,M,L)                    ! add point charge contribution

          ENDDO
          ENDDO
          ENDDO
        ENDIF

        DERHO = 0.D0

        DO L = 1, JZ
        DO M = 1, JY
        DO N = 1, JX
C         FRHO  = RHOA( N,M,L ) + RHOB( N,M,L )
          GRHO  = RHO( N,M,L )  - TRHO( N,M,L ) 
          DERHO = DERHO + GRHO**2*DV
          TRHO( N,M,L ) = RHO( N,M,L )
        ENDDO
        ENDDO
        ENDDO

        CALL MPI_REDUCE(DERHO,DERHO1,1,MPI_DOUBLE_PRECISION,
     *                  MPI_SUM,0,MPI_COMM_WORLD,IERR)
        DERHO = DERHO1
        CALL MPI_BCAST(DERHO,1,MPI_DOUBLE_PRECISION,
     *                 0,MPI_COMM_WORLD,IERR)

      ENDIF

C----- convergence judge -----

      DERHO = DSQRT(DERHO)
C     DERHO = DERHO / DBLE( MORA + MORB )

      IF((DERHO.LT.CONV) .OR.
     *   (NJ.EQ.ITMAX)) THEN
        GOTO 100
      ENDIF

C----- FOR ALPHA SPIN ORBITALS -----

C----- estimate kinetic energy term -----
C----- estimate effective potential (non-local potential) -----
C----- operate hamiltonian --------
C----- orbital energy -----

C     WRITE(*,*) 'KNTCS:start'          ! Vmol
C     CALL KNTCS( RWFA,NORA,RWFNA )     ! Vmol
C     WRITE(*,*) 'KNTCS:end'            ! Vmol
C     WRITE(*,*) 'KNTC :start'          ! Vmol
      TIME7 = MPI_WTIME()
      CALL KNTC( RWFA,NORA,RWFNA )
      TIME8 = MPI_WTIME()
C     IF(MYID.EQ.0) THEN
C       write(*,*) 'elapsed time:KNTC =',TIME8-TIME7
C     ENDIF
C     WRITE(*,*) 'rank ',MYID,'RWFA  =',RWFA(1,1,1,1)
C     WRITE(*,*) 'rank ',MYID,'RWFNA =',RWFNA(1,1,1,1)

C     WRITE(*,*) 'KNTC :end'            ! Vmol
      
C     TIME1 = MPI_WTIME()
      CALL GNALM( RWFA,GNALMA,DGNALMA,NORA )

C----- Hartree-Fock exchange potential -----  20180309 Takahasi

      CALL HF_EXCH_PSN( RWFA,RWFX,NORA,MORA,VHIJA,HFEXA )        ! MPI parallel
      HFEXA = RHEX*HFEXA
      RWFNA(:,:,:,:) = RWFNA(:,:,:,:) + RHEX*RWFX(:,:,:,:)     ! kinetic + hf-exchange

C     TIME2 = MPI_WTIME()
C     write(*,*) 'elapsed time:GNALM,MYID =',TIME2-TIME1,MYID
C     IF(MYID.EQ.0) THEN
C       write(*,*) 'elapsed time:GNALM =',TIME2-TIME1
C     ENDIF
      CALL OPERATE(RWFA,RWFNA,NORA,VEFFA,GNALMA,DGNALMA)
 
C     TIME3 = MPI_WTIME()
C     write(*,*) 'elapsed time:OPERATE,MYID =',TIME3-TIME2,MYID
C     IF(MYID.EQ.0) THEN
C       write(*,*) 'elapsed time:OPERATE =',TIME3-TIME2
C     ENDIF
      CALL OENGY( RWFA,RWFNA,ORBEA,NORA ) 

20    CONTINUE

      DO K = 1, NORA
        A1(K) = 0.D0
        A2(K) = 0.D0
      ENDDO

      DO K = 1, NORA
      DO L = 1, JZ
      DO M = 1, JY
      DO N = 1, JX
        RWFNA( N,M,L,K ) = ORBEA( K ) * RWFA( N,M,L,K )
     *                   - RWFNA( N,M,L,K ) 
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL PC( RWFNA,RWFNB,NORA ) 

      CALL KNTC( RWFNA,NORA,RWFNB )
      CALL GNALM( RWFNA,GNALMA,DGNALMA,NORA )
      CALL OPERATE(RWFNA,RWFNB,NORA,VEFFA,GNALMA,DGNALMA)

C     WRITE(*,*) 'rank ',MYID,'NJ =',NJ
      DO K = 1, NORA
      DO L = 1, JZ
      DO M = 1, JY
      DO N = 1, JX
      A1(K) = A1(K) + RWFNA( N,M,L,K )**2*DV
      A2(K) = A2(K) 
     *            + RWFNA( N,M,L,K )
     *      * (RWFNB( N,M,L,K ) - ORBEA(K)*RWFNA( N,M,L,K ))*DV
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE( A1(1),B1(1),NORA,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )
      CALL MPI_REDUCE( A2(1),B2(1),NORA,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

      IF(MYID.EQ.0) THEN
      DO K = 1, NORA
        A2(K) = B1(K)/B2(K) 
      ENDDO
      ENDIF

      CALL MPI_BCAST( A2(1),NORA,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )

      DO K = 1, NORA
      DO L = 1, JZ
      DO M = 1, JY
      DO N = 1, JX
        RWFA( N,M,L,K ) = RWFA( N,M,L,K ) 
     *                  + A2(K) * RWFNA( N,M,L,K )                         
      ENDDO
      ENDDO
      ENDDO
      ENDDO

C     WRITE(*,*) 'rank ',MYID,'NJ =',NJ
      IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
     *                 .OR.EXC.EQ.'RXalpha') THEN     
        CALL TENGY0( ORBEA,HFEXA,TOTE )                    ! for restricted orbitals
        GOTO 30
      ENDIF

C----- FOR BETA SPIN ORBITALS -----

C----- estimate kinetic energy term -----
       CALL KNTC( RWFB,NORB,RWFNB )
C----- estimate effective potential (non-local potential) -----
C----- operate hamiltonian --------
C----- orbital energy -----
C----- total energy -----

      CALL GNALM( RWFB,GNALMB,DGNALMB,NORB )

      IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN
        CALL HF_EXCH_PSN( RWFB,RWFX,NORB,MORB,VHIJB,HFEXB )          ! MPI parallel
      ENDIF
      HFEXB = RHEX*HFEXB
      RWFNB(:,:,:,:) = RWFNB(:,:,:,:) + RHEX*RWFX(:,:,:,:)     ! kinetic + hf-exchange

      CALL OPERATE(RWFB,RWFNB,NORB,VEFFB,GNALMB,DGNALMB)
      CALL OENGY(  RWFB,RWFNB,ORBEB,NORB ) 
      HFEX = HFEXA + HFEXB
      CALL TENGY1( ORBEA,ORBEB,HFEX,TOTE )            ! for unrestricted orbitals

      DO K = 1, NORB
        B1(K) = 0.D0
        B2(K) = 0.D0
      ENDDO

      DO K = 1, NORB
      DO L = 1, JZ
      DO M = 1, JY
      DO N = 1, JX
        RWFNB( N,M,L,K ) = ORBEB( K ) * RWFB( N,M,L,K )
     *                   - RWFNB( N,M,L,K ) 
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL PC( RWFNB,RWFNA,NORB ) 

      CALL KNTC(   RWFNB,NORB,RWFNA )
      CALL GNALM(  RWFNB,GNALMB,DGNALMB,NORB )
      CALL OPERATE(RWFNB,RWFNA,NORB,VEFFB,GNALMB,DGNALMB)

      DO K = 1, NORB
      DO L = 1, JZ
      DO M = 1, JY
      DO N = 1, JX
        B1(K) = B1(K) + RWFNB( N,M,L,K )**2*DV
        B2(K) = B2(K) 
     *        + RWFNB( N,M,L,K )
     *        * (RWFNA( N,M,L,K ) - ORBEB(K)*RWFNB( N,M,L,K ))*DV
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE( B1(1),A1(1),NORA,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )
      CALL MPI_REDUCE( B2(1),A2(1),NORA,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

      IF(MYID.EQ.0) THEN
      DO K = 1, NORB
        B2(K) = A1(K)/A2(K) 
      ENDDO
      ENDIF

      CALL MPI_BCAST( B2(1),NORA,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )

      DO K = 1, NORB
      DO L = 1, JZ
      DO M = 1, JY
      DO N = 1, JX
        RWFB( N,M,L,K ) = RWFB( N,M,L,K ) 
     *                  + B2(K) * RWFNB( N,M,L,K )                         
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL  RORD( ORBEB,NRDB,NORB )
C     CALL  DIAG(      NRDB,NORB,RWFB,RWFNB )
      CALL MDIAG (     NRDB,NORB,RWFB,RWFNB )
C     CALL MDIAG_BLK ( NRDB,NORB,RWFB,RWFNB )

      IF(NORB.GT.1) THEN
C       CALL RITZ_S( RWFNB,RWFB,NORB,VEFFB,ORBEB ) 
C       CALL RITZ  ( RWFNB,RWFB,NORB,VEFFB,ORBEB ) 
      ELSE
        DO K = 1, NORB
        DO L = 1, JZ
        DO M = 1, JY
        DO N = 1, JX
          RWFB( N,M,L,K ) = RWFNB( N,M,L,K ) 
        ENDDO
        ENDDO
        ENDDO
        ENDDO
      ENDIF

30    CONTINUE

      CALL RORD( ORBEA,NRDA,NORA )

C     TIME4 = MPI_WTIME()

C     CALL  DIAG(      NRDA,NORA,RWFA,RWFNA )
      CALL MDIAG (     NRDA,NORA,RWFA,RWFNA )
C     CALL MDIAG_BLK ( NRDA,NORA,RWFA,RWFNA )

C     TIME5 = MPI_WTIME()
C     IF(MYID.EQ.0) THEN
C       write(*,*) 'elapsed time:DIAG =',TIME5-TIME4
C     ENDIF
C     WRITE(*,*) 'rank',myid,'is alive. 1'

      IF(NORA.GT.1) THEN
C       CALL RITZ_S( RWFNA,RWFA,NORA,VEFFA,ORBEA ) 
C       CALL RITZ  ( RWFNA,RWFA,NORA,VEFFA,ORBEA ) 
C       TIME6 = MPI_WTIME()
C       IF(MYID.EQ.0) THEN
C         write(*,*) 'elapsed time:RITZ =',TIME6-TIME5
C       ENDIF
C     WRITE(*,*) 'rank',myid,'is alive. 2'
      ELSE
        DO K = 1, NORA
        DO L = 1, JZ
        DO M = 1, JY
        DO N = 1, JX
          RWFA( N,M,L,K ) = RWFNA( N,M,L,K ) 
        ENDDO
        ENDDO
        ENDDO
        ENDDO
      ENDIF

      PENG = TOTE + APOT

      IF(MYID.EQ.0) THEN 

      IF(PRINT.EQ.'LARGE') THEN
C       WRITE(*,*)                     
        WRITE(*,*) ' SCF   Iteration = ', NJ   
        WRITE(*,*) ' Total Energy    = ', PENG
        IF(NJ.NE.1) THEN
          DE = TOTE - PREE
          WRITE(*,*) ' Delta E         = ', DE
          WRITE(*,*) ' Residual(rho)   = ', DERHO
        ENDIF
      ENDIF

      ENDIF

      ETIME = MPI_WTIME()
      IF(MYID.EQ.0) THEN
        write(*,*) ' elapsed time    =', ETIME-STIME
      ENDIF

C     IF( IOEP == 0 ) THEN
      GOTO 10    ! SCF for the reference density
C     ELSE
C       GOTO 15    ! SCF for OEP
C     ENDIF

100   CONTINUE

C     CALL DNSHOMO( RWFA )           ! 2005.12.02 takahashi

      CALL KNTC( RWFA,NORA,RWFNA )

C----- Hartree-Fock exchange potential -----

      CALL HF_EXCH_PSN( RWFA,RWFX,NORA,MORA,VHIJA,HFEXA )     ! MPI parallel
      HFEXA = RHEX*HFEXA
      RWFNA(:,:,:,:) = RWFNA(:,:,:,:) + RHEX*RWFX(:,:,:,:)     ! kinetic + hf-exchange

C----- Local exchange potential by J.C.Slater -----

!     CALL OEP_LOCX( RWFA,RWFX,MORA,VEXA )      ! this call sentence is moved below.

C-------------------------------------------

      CALL GNALM(  RWFA,GNALMA,DGNALMA,NORA )
      CALL OPERATE(RWFA,RWFNA,NORA,VEFFA,GNALMA,DGNALMA)
      CALL OENGY(  RWFA,RWFNA,ORBEA,NORA ) 

      DO K = 1, NORA
        A1(K) = 0.D0
      ENDDO

      DO K = 1, NORA
      DO L = 1, JZ
      DO M = 1, JY
      DO N = 1, JX
        RWFNB( N,M,L,K ) = ORBEA( K ) * RWFA( N,M,L,K )
     *                   - RWFNA( N,M,L,K ) 
        A1(K) = A1(K) + RWFNB( N,M,L,K )**2*DV
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE( A1(K),B1(K),1,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )
      CALL MPI_BCAST( B1(K),1,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )

        IF((B1(K) .GE. EPS1 ) .AND.
     *     (NJ.LT.ITMAX)) THEN
          IF(MYID.EQ.0) THEN
          WRITE(*,*)                     
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,*) ' SCF Convergence is not achieved'
          WRITE(*,*) '        Continue SCF loop'
          WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          ENDIF
          CONV = CONV*0.1D0
          GOTO 20
        ENDIF

      ENDDO

C----- Local exchange potential by J.C.Slater ----- 

      CALL HF_EXCH_PSN( RWFA,RWFX,NORA,MORA,VHIJA,HFEXA )     ! MPI parallel
      CALL OEP_LOCX( RWFA,RWFX,MORA,VEXA )      ! 20230803 

C-------------------------------------------

      IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
     *                 .OR.EXC.EQ.'RXalpha') THEN     

C       CALL TENGY0( ORBEA,TOTE )              ! for restricted orbitals
        CALL TENGY0( ORBEA,HFEXA,TOTE )        ! for restricted orbitals

      ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
     *                     .OR.EXC.EQ.'UXalpha') THEN     

        CALL KNTC( RWFB,NORB,RWFNB )

C----- Hartree-Fock exchange potential -----

        CALL HF_EXCH_PSN( RWFB,RWFX,NORB,MORB,VHIJB,HFEXB )     ! MPI parallel
        HFEXB = RHEX*HFEXB
        RWFNB(:,:,:,:) = RWFNB(:,:,:,:) + RHEX*RWFX(:,:,:,:)     ! kinetic + hf-exchange

C----- Local exchange potential by J.C.Slater -----

!       CALL OEP_LOCX( RWFB,RWFX,MORB,VEXB )      ! 20230702 

C-------------------------------------------

        CALL GNALM(  RWFB,GNALMB,DGNALMB,NORB )
        CALL OPERATE(RWFB,RWFNB,NORB,VEFFB,GNALMB,DGNALMB)
        CALL OENGY(  RWFB,RWFNB,ORBEB,NORB ) 

        DO K = 1, NORB
          B1(K) = 0.D0
        ENDDO

        DO K = 1, NORB
        DO L = 1, JZ
        DO M = 1, JY
        DO N = 1, JX
          RWFNB( N,M,L,K ) = ORBEB( K ) * RWFB( N,M,L,K )
     *                     - RWFNB( N,M,L,K ) 
          B1(K) = B1(K) + RWFNB( N,M,L,K )**2*DV
        ENDDO
        ENDDO
        ENDDO

        CALL MPI_REDUCE( B1(K),A1(K),1,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )
        CALL MPI_BCAST( A1(K),1,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )

          IF((A1(K) .GE. EPS1 ) .AND.
     *       (NJ.LT.ITMAX)) THEN
            IF(MYID.EQ.0) THEN
            WRITE(*,*)                     
            WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(*,*) ' SCF Convergence is not achieved'
            WRITE(*,*) '        Continue SCF loop'
            WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            ENDIF
            CONV = CONV*0.1D0
            GOTO 20
          ENDIF

        ENDDO

        HFEX = HFEXA + HFEXB
        CALL TENGY1( ORBEA,ORBEB,HFEX,TOTE ) ! for unrestricted orbitals

C----- Local exchange potential by J.C.Slater ----- 

        CALL HF_EXCH_PSN( RWFB,RWFX,NORB,MORB,VHIJB,HFEXB )     ! MPI parallel
        CALL OEP_LOCX( RWFB,RWFX,MORB,VEXB )      ! 20230803 

C-------------------------------------------

      ENDIF

      PENG = TOTE + APOT
      DE = TOTE - PREE

      IF(MYID.EQ.0) THEN

      WRITE(*,*) ' SCF   Iteration = ', NJ   
      WRITE(*,*) ' Total Energy    = ', PENG
      WRITE(*,*) ' Delta E         = ', DE
      WRITE(*,*) ' Residual(rho)   = ', DERHO
      WRITE(*,*)                     
C     WRITE(*,996)

      KAZU = ' A'
      KAZU2= ' Alpha  ' 
      
      WRITE(*,*)                     
C     WRITE(*,997)

      CALL OROUT1(ORBEA,NORA,MORA,KAZU,KAZU2)

      IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
     *                 .OR.EXC.EQ.'UXalpha') THEN     

C     WRITE(*,999)                     

      KAZU = ' B'
      KAZU2= ' Beta   ' 
      
      CALL OROUT1(ORBEB,NORB,MORB,KAZU,KAZU2)

C     WRITE(*,999)                     
      WRITE(*,*)                     
      
      ENDIF

C     WRITE(*,998)

      KAZU = ' A'
      KAZU2= ' Alpha  ' 
      
      CALL OROUT2(A1,NORA,MORA,KAZU,KAZU2)

      IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
     *                 .OR.EXC.EQ.'UXalpha') THEN     

C     WRITE(*,999)                     

      KAZU = ' B'
      KAZU2= ' Beta   ' 
      
C     ENDIF

      CALL OROUT2(B1,NORB,MORB,KAZU,KAZU2)

      ENDIF
      
C     WRITE(*,999)                     

C----- output final energy -----

      WRITE(*,*)                     
      WRITE(*,*) '  Final Energy = ', PENG
      WRITE(30,*)  PENG                    

      ENDIF                                         ! IF(MYID.EQ.0) THEN

!     CALL DPL( RHO,IMD )                           ! for QM subsystem

      IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
     *                 .OR.EXC.EQ.'RXalpha') THEN     
        NEL = 2*MORA
        VCOU_FA(:,:,:) = DBLE(NEL-1)/DBLE(NEL)*VCOU(:,:,:)        ! Fermi-Amaldi potential as reference potential 20230623
        VEXA(:,:,:) = VEXA(:,:,:) + VCOU(:,:,:)/DBLE(NEL)         ! subtract the long-range part of the local exchange by Slater
      ELSE
        NEL = MORA + MORB
        VCOU_FA(:,:,:) = DBLE(NEL-1)/DBLE(NEL)*VCOU(:,:,:)        ! Fermi-Amaldi potential as reference potential 20230623
        VEXA(:,:,:) = VEXA(:,:,:) + VCOU(:,:,:)/DBLE(NEL)         ! subtract the long-range part of the local exchange by Slater
        VEXB(:,:,:) = VEXB(:,:,:) + VCOU(:,:,:)/DBLE(NEL)         ! subtract the long-range part of the local exchange by Slater
      ENDIF

      IOEP = 0           ! initialize the OEP step number
      CALL OEP_SCF( RWFA,RWFB,RWFNA,RWFNB,ORBEA,ORBEB,VNUC,VCOU_FA,VEXA,VEXB,RHO )
!     CALL OEP_SCF( RWFA,RWFB,RWFNA,RWFNB,ORBEA,      VNUC,VCOU_FA,VEXA,VEXB,RHO )
!     CALL OEP_SCF( RWFA,RWFB,RWFNA,RWFNB,ORBEA,VNUC,VCOU_FA,VEXA,VEXB     )

15    CONTINUE

!     CALL OEP_YW( RWFA,ORBEA,VNUC,VHIJA,VEXA,DELV,IOEP,RHO )     ! for rsg-oep based on the Yang-Wu approach(PRL,2002)  20230623
!     CALL OEP_YW( RWFA,ORBEA,VNUC,VHIJA,VEXA,DELV,IOEP     )     ! for rsg-oep based on the Yang-Wu approach(PRL,2002)  20230623
      CALL OEP_YW_UHF( RWFA,RWFB,ORBEA,ORBEB,VNUC,
     *                 VHIJA,VHIJB,VEXA,VEXB,DELV,IOEP,RHO )      ! for rsg-oep that can be used both in rhf and uhf-oep calculations.

C     REWIND(24)
      DO K = 1,JZ  
      DO J = 1,JY
      DO I = 1,JX
C       WRITE(24) VEXA(I,J,K) - VCOU(I,J,K)/DBLE(NEL-1)
      ENDDO
      ENDDO
      ENDDO

      IF( IOEP > 1000 ) THEN
        IF( MYID == 0 ) THEN
          WRITE(*,*) 
          WRITE(*,*) ' # of OEP cycles exeeds 1000. '
          WRITE(*,*) ' Aborting OEP-SCF! '
          WRITE(*,*) 
        ENDIF

C       REWIND(24)
        DO K = 1,JZ  
        DO J = 1,JY
        DO I = 1,JX
C         WRITE(24) VEXA(I,J,K) - VCOU_FA(I,J,K)/DBLE(NEL-1)
C         WRITE(24) VEXA(I,J,K) 
        ENDDO
        ENDDO
        ENDDO

        RETURN
      ENDIF

!---- convergence judge for OEP-SCF ----

      IF( CONV_OEP < DELV ) THEN
        CALL OEP_SCF( RWFA,RWFB,RWFNA,RWFNB,ORBEA,ORBEB,VNUC,VCOU_FA,VEXA,VEXB,RHO )
!       CALL OEP_SCF( RWFA,RWFB,RWFNA,RWFNB,ORBEA,      VNUC,VCOU_FA,VEXA,VEXB,RHO )
!       CALL OEP_SCF( RWFA,RWFB,RWFNA,RWFNB,ORBEA,VNUC,VCOU_FA,VEXA,VEXB     )
        GOTO 15
      ELSE
        IF( MYID == 0 ) THEN
          WRITE(*,*) ' ---- SCF for OEP Converged! ----'
        ENDIF

        REWIND(24)
        DO K = 1,JZ  
        DO J = 1,JY
        DO I = 1,JX
          WRITE(24) VEXA(I,J,K) - VCOU_FA(I,J,K)/DBLE(NEL-1)
!         WRITE(24) VEXA(I,J,K) 
!         WRITE(24) RHO( I,J,K) 
        ENDDO
        ENDDO
        ENDDO

      ENDIF 

C----- output dipole moment -----

      CALL DPL( RHO,IMD )                           ! for QM subsystem

!     CALL FLUSH(06)                                ! for SX-ACE
!     CALL FLUSH(30)                                ! for SX-ACE

C----- compute energy distribution function -----

C     CALL EDF( RHO,SCRD1,PENG,IMD )                ! takahashi 2005.09.13

C----- output electronic density -----

C     IF(MYID.EQ.0) THEN
C     CALL WDNS( IMD,RHO )
C     ENDIF

C----------------------------------

      IF(FREEZE.EQ.'FALSE') THEN

      CALL LFRC( RHO ) ! Force (local part) for non-periodic system

      IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
     *                 .OR.EXC.EQ.'RXalpha') THEN     
        BF = 2.D0
        CALL NLFRC( RWFA,GNALMA,DGNALMA,MORA,BF ) ! Force (non-local part)
      ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
     *                 .OR.EXC.EQ.'UXalpha') THEN     
        BF = 1.D0
        CALL NLFRC( RWFA,GNALMA,DGNALMA,MORA,BF ) ! Force (non-local part)
        CALL NLFRC( RWFB,GNALMB,DGNALMB,MORB,BF ) ! Force (non-local part)
      ENDIF

      ELSEIF(NQMMM.EQ.'LINK') THEN

      CALL LFRC1( RHO ) ! Force (local part) for non-periodic system

      IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
     *                 .OR.EXC.EQ.'RXalpha') THEN     
        BF = 2.D0
C       CALL NLFRC1( RWFA,GNALMA,MORA,BF ) ! Force (non-local part)
        CALL NLFRC2( RWFA,GNALMA,MORA,BF,0 ) ! Force (non-local part)
      ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
     *                 .OR.EXC.EQ.'UXalpha') THEN     
        BF = 1.D0
C       CALL NLFRC1( RWFA,GNALMA,MORA,BF ) ! Force (non-local part)
C       CALL NLFRC1( RWFB,GNALMB,MORB,BF ) ! Force (non-local part)
        CALL NLFRC2( RWFA,GNALMA,MORA,BF,0 ) ! Force (non-local part)
        CALL NLFRC2( RWFB,GNALMB,MORB,BF,1 ) ! Force (non-local part)
      ENDIF

      ENDIF

C     CALL OPTFC
C     CALL OPTFC( IMD )
      CALL OPTFC3

      CALL AVDNST(RHO, IMD)
C     CALL AVDNST(RWFA,IMD)

      CALL FLUSH(24)                                ! for SX-ACE

      CALL REDDNS(RHO,IMD)

777   CONTINUE

C     CALL FUZZYC                             ! Vmol
      CALL EFTC(RHOA,RHOB)                          ! Vmol
C     WRITE(*,*) 'PCLMB =',PCLMB              ! Vmol
C----------------------------------------------------------------------
      IF(NQMMM.EQ.'QMMM') THEN

        CALL QMF( SCRD,FRCS,FRCN,RHO )

        DO J=1,NDIM
        DO I=1,NNUC
          FRC(I,J) = FRC(I,J) + FRCN(I,J)
        ENDDO
        ENDDO

      ELSEIF(NQMMM.EQ.'LINK') THEN            ! Vmol
                                              ! Vmol
        CALL SSINT1( RHO,IMD )                ! Vmol  2013.08.05 takahashi
C       CALL QMF1(   SCRD, FRCS,FRCN,RHO )    ! Vmol
C       CALL QMF1(   SCRD1,FRCS,FRCN,RHO )    ! Vmol
!       CALL RQMF1_S(SCRD1,FRCS,FRCN     )    ! Vmol
                                              ! Vmol
C        IF(MYID.EQ.0) THEN
CC         DO I=NNUC-NLINK+1,NNUC                           ! Vmol
C          DO I=1,NNUC                           ! Vmol
C            write(*,199) 'FRC=  ',FRC(I,1),FRC(I,2),FRC(I,3)     ! Vmol
C          ENDDO                                 ! Vmol
CC         DO I=NNUC-NLINK+1,NNUC                           ! Vmol
C          DO I=1,NNUC                           ! Vmol
C            write(*,199) 'FRCN= ',FRCN(I,1),FRCN(I,2),FRCN(I,3)     ! Vmol
C          ENDDO                                 ! Vmol
C        ENDIF

        DO J=1,NDIM                           ! Vmol
        DO I=1,NNUC                           ! Vmol
          FRC(I,J) = FRC(I,J) + FRCN(I,J)     ! Vmol
        ENDDO                                 ! Vmol
        ENDDO                                 ! Vmol
                                              ! Vmol
        IF(NLINK.GT.0) THEN                   ! Vmol
C         CALL EXTCORR( FRCS,FRC,SCRD )       ! Vmol
C         write(*,*) 'spltfrc:start'
C         CALL SPLTFRC( FRCS,FRC,SCRD )       ! Vmol
!         CALL SPLTFRC( FRCS,FRC,SCRD1 )      ! Vmol
C         write(*,*) 'spltfrc:end'
        ENDIF                                 ! Vmol
                                              ! Vmol
      ENDIF

      IF(MYID.EQ.0) THEN
      WRITE(*,*)                     
C     WRITE(*,995)                     
      ENDIF

197   FORMAT( 1X,A22,3F12.6 )      
199   FORMAT( 3X,A6,3F12.6 )      

994   FORMAT( '!!!!!!!!!!!!!!!!!!!!!!
     *!!!!!! Self Consistent Field
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!' )
995   FORMAT( '!!!!!!!!!!!!!!!!!!!!!!!!
     *!!!!!!!!!!!!!!!!!!!!!!!!
     *!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' )
996   FORMAT( '!!!!!!!!!!!!!!!!!!!!!!
     *!!!!!!!! Density Converged
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' )
997   FORMAT( '====================
     *======== Orbital Energy [a.u.]
     * =============================' )
998   FORMAT( '==========================
     * Residual for each Orbital
     * ==========================' )
999   FORMAT( '=====================
     *==============================
     *============================' )

      RETURN
      END

C---------------------------------------------
C   subroutine output orbital and residual
C---------------------------------------------

      SUBROUTINE OROUT1(ORBEA,NOR,MOR,KAZU,KAZU2)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NORA =  49 )
      
      DIMENSION ORBEA(  NORA )
 
      CHARACTER NO1*1,NO2*1,KAZU*2,KAZU1*6,KAZU2*8

C----- output orbital -----

      WRITE(*,*) 
      WRITE(*,*) ' ---- Orbital EigenValues ---- '

      NA = NOR - MOR

      IF(MOR.LE.5) THEN                             ! for alpha spin
        KAZU1 = KAZU//'1   '
        WRITE(*,99) KAZU2,KAZU1,( ORBEA(I),I=1,MOR )
      ELSE
        KAZU1 = KAZU//'1   '
        WRITE(*,99) KAZU2,KAZU1,( ORBEA(I),I=1,5 )
        IF(MOD(MOR,5).EQ.0) THEN
          ILOOP = MOR/5
        ELSE
          ILOOP = MOR/5 + 1
        ENDIF
        DO N=2,ILOOP
          IS = 5*N-4
          IE = 5*N
          IF(IE.GT.MOR) IE=MOR
          IF(N.LT.10) THEN
            NO1 = CHAR(N+48)
            KAZU1 = KAZU//NO1//'   '
            WRITE(*,96) KAZU1,( ORBEA(I),I=IS,IE )
          ELSEIF(N.LT.100) THEN
            IN = N/10
            NO2 = CHAR(IN+48)
            IN = MOD(N,10)
            NO1 = CHAR(IN+48)
            KAZU1 = KAZU//NO2//NO1//'  '
            WRITE(*,96) KAZU1,( ORBEA(I),I=IS,IE )
          ENDIF
        ENDDO
      ENDIF

      IF(NA.GT.0) THEN
        IF(NA.LE.5) THEN
          KAZU1 = KAZU//'1 V '
          WRITE(*,97) KAZU1,( ORBEA(I),I=MOR+1,NOR )
        ELSE
          KAZU1 = KAZU//'1 V '
          WRITE(*,97) KAZU1,( ORBEA(I),I=MOR+1,MOR+5 )
          IF(MOD(NA,5).EQ.0) THEN
            ILOOP = NA/5
          ELSE
            ILOOP = NA/5 + 1
          ENDIF
          DO N=2,ILOOP
            IS = 5*N-4+MOR
            IE = 5*N+MOR
            IF(IE.GT.NOR) IE=NOR
            IF(N.LT.10) THEN
              NO1 = CHAR(N+48)
              KAZU1 = KAZU//NO1//' V '
              WRITE(*,96) KAZU1,( ORBEA(I),I=IS,IE )
            ELSEIF(N.LT.100) THEN
              IN = N/10
              NO2 = CHAR(IN+48)
              IN = MOD(N,10)
              NO1 = CHAR(IN+48)
              KAZU1 = KAZU//NO2//NO1//' V'
              WRITE(*,96) KAZU1,( ORBEA(I),I=IS,IE )
            ENDIF
          ENDDO
        ENDIF
      ENDIF

96    FORMAT(             8X,A6,5F12.6 )      
97    FORMAT( 'virtual ',    A6,5F12.6 )      
99    FORMAT(             A8,A6,5F12.6 )      

      RETURN

      END

      SUBROUTINE OROUT2(A1,NOR,MOR,KAZU,KAZU2)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NORA =  49 )
      
      DIMENSION A1( NORA )
 
      CHARACTER NO1*1,NO2*1,KAZU*2,KAZU1*6,KAZU2*8

C----- output residual -----

      WRITE(*,*) ' ---- Residuals of Orbitals --- '

      NA = NOR - MOR

      IF(MOR.LE.5) THEN                             ! for alpha spin
        WRITE(*,89) KAZU2,( A1(I)   ,I=1,MOR )
      ELSE
        WRITE(*,89) KAZU2,( A1(I),I=1,5 )
        WRITE(*,86) ( A1(I),I=6,MOR )
      ENDIF

      IF(NA.GT.0) THEN
        IF(NA.LE.5) THEN
          WRITE(*,87) ( A1(I),I=MOR+1,NOR )
        ELSE
          WRITE(*,87) ( A1(I),I=MOR+1,MOR+5 )
          WRITE(*,86) ( A1(I),I=MOR+6,NOR )
        ENDIF
      ENDIF

86    FORMAT(               14X,5E12.4 )      
87    FORMAT( 'virtual ',    6X,5E12.4 )      
89    FORMAT(             A8,6X,5E12.4 )      

      RETURN

      END

C-----------------------------------
C  SUBROUTINE OPERATE HAMILTONIAN
C-----------------------------------

      SUBROUTINE OPERATE(RWF,RWFN,NOR,VEFF,TGNALM,DGNALM)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol
      include 'nlocd.i'                         ! Vmol   non-local d

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NORA =  49 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( RCUT =  2.5D0)
C     PARAMETER ( NSL  =  3000 )

      PARAMETER ( PI   = 3.14159265358979323D0 )

      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / DBLG  / WNLOCS( NNUC,NSL ),WNLOCPX( NNUC,NSL ),
     *                 WNLOCPY(NNUC,NSL ),WNLOCPZ( NNUC,NSL )

C     DIMENSION RWF( NMAX,NMAX,NMAX,NORA)
C     DIMENSION RWFN(NMAX,NMAX,NMAX,NORA)
      DIMENSION RWF( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NOR)
      DIMENSION RWFN(NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NOR)
      DIMENSION VEFF(NMAXX/NX,NMAXY/NY,NMAXZ/NZ)
      DIMENSION TGNALM(NNUC,2,3,NORA)            ! num of elctrns,atoms,angular momentum quantum numbers
      DIMENSION DGNALM(NNUC,  5,NORA)            ! for non-local d pseudo pot.
      DIMENSION NCT(NNUC)

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      COM  = 1.D0 / DSQRT( 4.D0*PI )

      RCUT2 = RCUT*RCUT
      DO NA=1,NNUC                               ! initialize each counter
        NCT(NA) = 0
      ENDDO

C----- kinetic and local potential

      DO L = 1, NOR
      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

        RWFN( I,J,K,L ) = RWFN( I,J,K,L ) 
     *                  + VEFF(I,J,K) * RWF( I,J,K,L )

      ENDDO
      ENDDO
      ENDDO
      ENDDO

C----- non local ps potential -----

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO 30 NA = 1, NNUC

        NJ = NPD(NA)

      SELECT CASE (NJ)

      CASE ( 0 )    ! NJ = 0: for s and p non-local

      DO 40 K = 1, JZ
        RZ = DZ*( K-NNZ ) - PNUC(NA,3)
      DO 50 J = 1, JY
        RY = DY*( J-NNY ) - PNUC(NA,2)
        R2T = RZ**2 + RY**2
C       IF( R2T >= RCUT2 ) CYCLE
      DO 60 I = 1, JX

        RX = DX*( I-NNX ) - PNUC(NA,1) 
        R2 = RX**2 + R2T
C       R = DSQRT(R2)

          IF( R2 .LT. RCUT2 ) THEN
            NCT(NA) = NCT(NA) + 1

            DO 70 L  = 1, NOR
              VPSUM = 0.D0

C             ----- s-component ----- 

              VPSUM = VPSUM  
     *          +     COM * WNLOCS( NA,NCT(NA))*TGNALM( NA,1,1,L )

C             ----- p-component ----- 

              VPSUM = VPSUM                                             ! Px component
     *          +  3.D0*COM * WNLOCPX(NA,NCT(NA))*TGNALM( NA,2,1,L )
              VPSUM = VPSUM                                             ! Py component
     *          +  3.D0*COM * WNLOCPY(NA,NCT(NA))*TGNALM( NA,2,2,L )
              VPSUM = VPSUM                                             ! Pz component
     *          +  3.D0*COM * WNLOCPZ(NA,NCT(NA))*TGNALM( NA,2,3,L )

              RWFN( I,J,K,L ) = RWFN( I,J,K,L ) 
     *                        + VPSUM

70          CONTINUE

          ENDIF

60    CONTINUE
50    CONTINUE
40    CONTINUE

      CASE ( 1 )    ! NJ = 1: for s,p and d non-local

      DO 45 K = 1, JZ
        RZ = DZ*( K-NNZ ) - PNUC(NA,3)
      DO 55 J = 1, JY
        RY = DY*( J-NNY ) - PNUC(NA,2)
        R2T = RZ**2 + RY**2
C       IF( R2T >= RCUT2 ) CYCLE
      DO 65 I = 1, JX

        RX = DX*( I-NNX ) - PNUC(NA,1) 
        R2 = RX**2 + R2T
C       R = DSQRT(R2)

          IF( R2 .LT. RCUT2 ) THEN
            NCT(NA) = NCT(NA) + 1

            DO 75 L  = 1, NOR
              VPSUM = 0.D0

C             ----- s-component ----- 

              VPSUM = VPSUM  
     *          +     COM * WNLOCS( NA,NCT(NA))*TGNALM( NA,1,1,L )

C             ----- p-component ----- 

              VPSUM = VPSUM                                             ! Px component
     *          +  3.D0*COM * WNLOCPX(NA,NCT(NA))*TGNALM( NA,2,1,L )
              VPSUM = VPSUM                                             ! Py component
     *          +  3.D0*COM * WNLOCPY(NA,NCT(NA))*TGNALM( NA,2,2,L )
              VPSUM = VPSUM                                             ! Pz component
     *          +  3.D0*COM * WNLOCPZ(NA,NCT(NA))*TGNALM( NA,2,3,L )

C             ----- d-component ----- 

              VPSUM = VPSUM                                             ! Dxy component
     *          + 15.D0*COM * WNLOCDXY(NA,NCT(NA))*DGNALM( NA,1,L )
              VPSUM = VPSUM                                             ! Dyz component
     *          + 15.D0*COM * WNLOCDYZ(NA,NCT(NA))*DGNALM( NA,2,L )
              VPSUM = VPSUM                                             ! Dzx component
     *          + 15.D0*COM * WNLOCDZX(NA,NCT(NA))*DGNALM( NA,3,L )
              VPSUM = VPSUM                                             ! Dz2 component
     *          + 5.D0/4.D0 * COM * WNLOCDZ2(NA,NCT(NA))*DGNALM( NA,4,L )
              VPSUM = VPSUM                                             ! Dx2-y2 component
     *          +15.D0/4.D0 * COM * WNLOCDX2(NA,NCT(NA))*DGNALM( NA,5,L )

              RWFN( I,J,K,L ) = RWFN( I,J,K,L ) 
     *                        + VPSUM

75          CONTINUE

          ENDIF

65    CONTINUE
55    CONTINUE
45    CONTINUE

      END SELECT

30    CONTINUE

      RETURN
      END

C----------------------------
C   SUBROUTINE RAYLEIGH-RITZ 
C----------------------------

      SUBROUTINE RITZ( RWF,RWFN,NOR,VEFF,ORBE ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX =  80 )
C     PARAMETER ( NDIM =   3 )
C     PARAMETER ( NORA =  49 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( RCUT =  2.5D0)
C     PARAMETER ( NSL  =  3000 )
C     PARAMETER ( NCYCLE  = 100 )

      PARAMETER ( PI   = 3.14159265358979323D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ
      COMMON / NCLR  / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)

C     DIMENSION RWF(   NMAX,NMAX,NMAX,NORA )            ! trial wave function
C     DIMENSION RWFN(  NMAX,NMAX,NMAX,NORA )            ! true wave function
      DIMENSION RWF(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NOR )    ! trial wave function
      DIMENSION RWFN(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NOR )    ! true wave function
      DIMENSION VEFF(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION ORBE(  NORA )
      DIMENSION ORBE1( NORA )
      DIMENSION HOR(   NOR,NOR )
      DIMENSION HORA(  NOR,NOR )
      DIMENSION FAI(   NOR,NOR )
      DIMENSION TGNALM(NNUC,2,3,NORA )         ! NUM OF ELCTRNS,ATOMS,ANGULAR MOMENTUM QUANTUM NUMBERS
      DIMENSION DGNALM(NNUC,  5,NORA )         ! NUM OF ELCTRNS,ATOMS,ANGULAR MOMENTUM QUANTUM NUMBERS
      DIMENSION NCT(   NNUC )

      DIMENSION AL(NOR)
      DIMENSION BE(NOR)
      DIMENSION CO(NOR)
      DIMENSION W( NOR)
      DIMENSION P( NOR)
      DIMENSION Q( NOR)
      DIMENSION B2(NOR)
      DIMENSION DI(NOR)
      DIMENSION BL(NOR)
      DIMENSION BU(NOR)
      DIMENSION BV(NOR)
      DIMENSION CM(NOR)
      DIMENSION LEX(NOR)

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      DV   = DX*DY*DZ

C----- ESTIMATE KINETIC ENERGY TERM -----

      CALL KNTC( RWF,NOR,RWFN )

C----- ESTIMATE EFFECTIVE POTENTIAL (NON-LOCAL POTENTIAL) -----

      CALL GNALM( RWF,TGNALM,DGNALM,NOR )

C----- OPERATE HAMILTONIAN --------

      CALL OPERATE(RWF,RWFN,NOR,VEFF,TGNALM,DGNALM)

C===== CALCULATE MATRIX ELEMENT =====

C----- INITIALIZE -----

      DO NO = 1,NOR
      DO MO = 1,NOR

      HORA(NO,MO) = 0.D0

      ENDDO
      ENDDO

      DO NO = 1,NOR
      DO MO = 1,NO
      DO L  = 1,JZ
      DO M  = 1,JY
      DO N  = 1,JX

      HORA(NO,MO) = HORA(NO,MO)
     *            + RWF(N,M,L,NO)*RWFN(N,M,L,MO)*DV

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      NNDAT = NOR*NOR
      CALL MPI_REDUCE(HORA(1,1),HOR(1,1),NNDAT,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
C     CALL MPI_BCAST(HOR(1,1),NNDAT,MPI_DOUBLE_PRECISION,
C    *                0,MPI_COMM_WORLD,IERR)

C===== GENERATE RITZ VECTOR =====

      IF(MYID.EQ.0) THEN
      CALL HOUSE(NOR,HOR,AL,BE,CO,W,P,Q)

      NLS = -1
      CALL BISEC(NOR,NLS,ORBE,AL,BE,B2)

      CALL INVITR(NOR,HOR,ORBE,FAI,AL,BE,CO,
     *            DI,BL,BU,BV,CM,LEX)
      ENDIF

      CALL MPI_BCAST(FAI(1,1),NNDAT,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)
      

C===== GENERATE NEW WAVE FUNCTIONS =====

      DO NO = 1,NOR
      DO L  = 1,JZ
      DO M  = 1,JY
      DO N  = 1,JX
        RWFN(N,M,L,NO) = 0.D0
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO MO = 1,NOR
      DO NO = 1,NOR
      DO L  = 1,JZ
      DO M  = 1,JY
      DO N  = 1,JX

      RWFN(N,M,L,NO) = RWFN(N,M,L,NO)
     *               + FAI(MO,NO)*RWF(N,M,L,MO)

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO NO = 1,NOR
      ORBE(NO) = 0.D0
      ENDDO

      DO NO = 1,NOR
      DO L  = 1,JZ
      DO M  = 1,JY
      DO N  = 1,JX

      ORBE(NO) = ORBE(NO) + RWFN(N,M,L,NO)*RWFN(N,M,L,NO)*DV

      ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE(ORBE(1),ORBE1(1),NOR,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)

      DO NO = 1,NOR
      ORBE(NO) = DSQRT(ORBE1(NO))
      ENDDO

      CALL MPI_BCAST(ORBE(1),NOR,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

      DO NO = 1,NOR
      DO L  = 1,JZ
      DO M  = 1,JY
      DO N  = 1,JX

      RWFN(N,M,L,NO) = RWFN(N,M,L,NO)/ORBE(NO)

      ENDDO
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END


C----------------------------------
C   SUBROUTINE REORDER ORBITALS
C----------------------------------

      SUBROUTINE RORD( ORBE,NRD,NOR ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NORA = 49 )

      DIMENSION ORBE( NORA )
C     DIMENSION NRD(  NOR  )
      DIMENSION NRD(  NORA )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      EMAX = ORBE( 1 )

      DO 10 I = 1, NOR
      EMAX = DMAX1( EMAX, ORBE(I) )
10    CONTINUE

      DO 20 K=1,NOR

      EMIN = EMAX

        DO 30 I=1,NOR

        EMIN = DMIN1( EMIN, ORBE(I) )

        IF( ORBE(I) .EQ. EMIN ) THEN 
        ITMP = I    
        ENDIF

30      CONTINUE

        ORBE(ITMP) = EMAX + 1.D0
        NRD(K)     = ITMP

20    CONTINUE

      RETURN
      END

C-------------------------------------------
C   SUBROUTINE COMPUTING ORBITAL ENERGIES
C-------------------------------------------

      SUBROUTINE OENGY( RWF,RWFN,ORBE,NOR ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NORA = 49 )
   
      COMMON / MMPI4  / JX,JY,JZ
      COMMON / GRID2  / DX, DY, DZ
      COMMON / ELCTRN / AX1, AX2, AY1, AY2, AZ1, AZ2, AO,
     *                  AX3, AX4, AY3, AY4, AZ3, AZ4

      DIMENSION ORBE1( NORA )
      DIMENSION ORBE(  NORA )

C     DIMENSION RWF(  NMAX,NMAX,NMAX,NORA )
C     DIMENSION RWFN( NMAX,NMAX,NMAX,NORA )
      DIMENSION RWF(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NOR )
      DIMENSION RWFN( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NOR )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      DV = DX*DY*DZ

C----- initialize -----

      DO K = 1, NOR
        ORBE1( K ) = 0.D0
      ENDDO

C----- orbital energies -----

      DO K = 1, NOR
      DO L = 1, JZ
      DO M = 1, JY
      DO N = 1, JX

        ORBE1( K ) = ORBE1( K ) + RWF( N,M,L,K )  
     *             * RWFN( N,M,L,K )*DV

      ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE(ORBE1(1),ORBE(1),NOR,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST( ORBE(1),NOR,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

C     IF (MYID .EQ. 0) THEN
C     WRITE(*,'(6E13.5)') (ORBE(N),N=1,5)
C     WRITE(*,*) ORBE(1),ORBE(NORA)
C     ENDIF

      RETURN
      END


C------------------------------------
C   SUBROUTINE KINETIC ( Simple Version ) 
C------------------------------------

C     SUBROUTINE KNTC( RWF,NOR,TWF )
C     SUBROUTINE KNTCS( RWF,NOR,TWF )

C     IMPLICIT REAL*8 ( A-H,O-Z )      
C     IMPLICIT INTEGER*4 ( I-N )      

C     include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NORA = 49 )
C     PARAMETER ( INUM =  8 )
C  
C     COMMON / ELCTRN / AX1, AX2, AY1, AY2, AZ1, AZ2, AO,
C    *                  AX3, AX4, AY3, AY4, AZ3, AZ4

C     DIMENSION RWF( NMAX,NMAX,NMAX,NORA )
C     DIMENSION TWF( NMAX,NMAX,NMAX,NORA )
C     DIMENSION RWF( NMAX,NMAX,NMAX,NOR )
C     DIMENSION TWF( NMAX,NMAX,NMAX,NOR )
C     DIMENSION TWF1( NMAX+INUM )

C     ZERO    = 0.D0
C     TWF1(1) = ZERO
C     TWF1(2) = ZERO
C     TWF1(3) = ZERO
C     TWF1(4) = ZERO
C     TWF1(NMAX+5) = ZERO
C     TWF1(NMAX+6) = ZERO
C     TWF1(NMAX+7) = ZERO
C     TWF1(NMAX+8) = ZERO

C     NS = 5
C     NE = NMAX+4

C----- Z-AXIS -----
C     
C     DO 210 LL = 1, NOR
C     DO 10 J = 1, NMAX
C       DO 20 I = 1, NMAX

C       DO K = NS, NE
C         TWF1(K) = RWF( I,J,K-4,LL )
C       ENDDO

C       DO 30 K = NS, NE

C       TWF( I,J,K-4,LL ) = AZ1 * ( TWF1( K-1 )
C    *                    +         TWF1( K+1 ))
C    *                    + AZ2 * ( TWF1( K-2 )
C    *                    +         TWF1( K+2 ))
C    *                    + AZ3 * ( TWF1( K-3 )
C    *                    +         TWF1( K+3 ))
C    *                    + AZ4 * ( TWF1( K-4 )
C    *                    +         TWF1( K+4 ))
C    *                    + AO  *   TWF1( K )

C30      CONTINUE

C20      CONTINUE
C10    CONTINUE
C210   CONTINUE

C----- Y-AXIS -----
C     
C     DO 220 LL = 1, NOR
C     DO 40 K = 1, NMAX
C       DO 50 I = 1, NMAX

C       DO J = NS, NE
C         TWF1(J) = RWF( I,J-4,K,LL )
C       ENDDO

C       DO 60 J = NS, NE

C       TWF( I,J-4,K,LL ) = AY1 * ( TWF1( J-1 )
C    *                    +         TWF1( J+1 ))
C    *                    + AY2 * ( TWF1( J-2 )
C    *                    +         TWF1( J+2 ))
C    *                    + AY3 * ( TWF1( J-3 )
C    *                    +         TWF1( J+3 ))
C    *                    + AY4 * ( TWF1( J-4 )
C    *                    +         TWF1( J+4 ))
C    *                    + TWF( I,J-4,K,LL )

C60      CONTINUE

C50      CONTINUE
C40    CONTINUE
C220   CONTINUE

C----- X-AXIS -----
C     
C     DO 230 LL = 1, NOR
C     DO 70 K = 1, NMAX
C       DO 80 J = 1, NMAX

C       DO I = NS, NE
C         TWF1(I) = RWF( I-4,J,K,LL )
C       ENDDO

C       DO 90 I = NS, NE

C       TWF( I-4,J,K,LL ) = AX1 * ( TWF1( I-1 )
C    *                    +         TWF1( I+1 ))
C    *                    + AX2 * ( TWF1( I-2 )
C    *                    +         TWF1( I+2 ))
C    *                    + AX3 * ( TWF1( I-3 )
C    *                    +         TWF1( I+3 ))
C    *                    + AX4 * ( TWF1( I-4 )
C    *                    +         TWF1( I+4 ))
C    *                    + TWF( I-4,J,K,LL )

C90      CONTINUE

C80      CONTINUE
C70    CONTINUE
C230   CONTINUE

C     RETURN
C     END

C---------------------------------------
C   SUBROUTINE COMPUTING TOTAL ENERGY  
C---------------------------------------

      SUBROUTINE TENGY0( ORBEA,HFEX,TOTE ) ! for restricted orbitals

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NORA =  49 )
C     PARAMETER ( MORA =  49 )
     
      COMMON / CLMB / PCLMB, PEXC, PEXT

      DIMENSION ORBEA( NORA )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      SUM = 0.D0

C----- for alpha spin -----

      DO 10 I = 1, MORA
        SUM = SUM + ORBEA( I )
10    CONTINUE

      TOTE = 2.D0*SUM - PCLMB + PEXC - PEXT - 2.D0*HFEX 

      IF (MYID .EQ. 0) THEN
C     WRITE(*,'(6E13.5)') (ORBEA(N),N=1,5)
C     WRITE(*,*) ORBE(1),ORBE(NORA)

C     write(*,*) 'orb_e=', SUM*2.D0
C     write(*,*) 'pclmb=', PCLMB 
C     write(*,*) 'pexc=',  PEXC 
C     write(*,*) 'pext=', -PEXT 
C     write(*,*) 'pe=',    TOTE 
      ENDIF

      RETURN
      END

      SUBROUTINE TENGY1( ORBEA,ORBEB,HFEX,TOTE ) ! for unrestricted orbitals

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NORA =  49 )
C     PARAMETER ( MORA =  49 )
C     PARAMETER ( MORB =  1 )
     
      COMMON / CLMB / PCLMB, PEXC, PEXT

      DIMENSION ORBEA( NORA )
      DIMENSION ORBEB( NORA )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      SUM = 0.D0

C----- for alpha spin -----

      DO 10 I = 1, MORA
        SUM = SUM + ORBEA( I )
10    CONTINUE

C----- for beta spin -----

      DO 15 I = 1, MORB
        SUM = SUM + ORBEB( I )
15    CONTINUE

      TOTE = SUM - PCLMB + PEXC - PEXT - HFEX

C     write(*,*) 'ek=',   EK 
C     write(*,*) 'pclmb=',PCLMB 
C     write(*,*) 'pe=',   PE 
C     write(*,*) 'pexc=', PEXC 
C     write(*,*) 'pext=', PEXT 

      RETURN
      END

C--------------------------------------------------------
C   SUBROUTINE PS_DATAA ( Read Pseudo-potential data NCPS97  
C   Kleinman and Bylander form )                       
C   Phys. Rev. Lett. 48, 1425( 1982 )
C   NCPS 97
C--------------------------------------------------------

      SUBROUTINE PSDAT

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'QMpara.i'                        ! Vmol
      include 'nlocd.i'                         ! Vmol   non-local d

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NORB =  1 )
C     PARAMETER ( NOCC =  1 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NSL  = 3000  )

      PARAMETER ( PI = 3.14159265358979323D0 )

      DIMENSION NREAD( 100 )

      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)

C--------------------------------------
C      Read number of grid points 
C      dummy                           
C      delta r in log scale from HKB.DAT          ! read from HKB.DAT
C--------------------------------------

      READ(41,*) NREAD(1),DUMMY,DRLOG(1)
      DRLOG(1) = DRLOG(1) * 3.D0                  ! true grid spacing in log-mesh

C----- Read denominators -----

      READ(41,*) SVS(1),PVP(1)                    ! for s and p pseudowavefunction
C
C----- Read grid points in real space -----
C
      K = 1
      DO I=1,INT(NREAD(1)/4)
      READ(41,*) RAD(1,K), RAD(1,K+1), RAD(1,K+2), RAD(1,K+3)
      K = K + 4
      ENDDO
      READ(41,*) RAD(1,K)
C
C----- Read |Vnls|Ys> -----                      ! caution:Ys is multiplied by r
C                                                ! Vnls = Vs - Vd
      K = 1
      DO I=1,INT(NREAD(1)/4)
      READ(41,*) VNLS(1,K), VNLS(1,K+1), VNLS(1,K+2), VNLS(1,K+3)
      K = K + 4
      ENDDO
      READ(41,*) VNLS(1,K)
C
C----- Read |Vnlp|Yp> -----                      ! caution:Yp is multiplied by r
C                                                ! Vnls = Vs - Vd
      K = 1
      DO I=1,INT(NREAD(1)/4)
      READ(41,*) VNLP(1,K), VNLP(1,K+1), VNLP(1,K+2), VNLP(1,K+3)
      K = K + 4
      ENDDO
      READ(41,*) VNLP(1,K)

      DO I=1,NREAD(1)
      VNLS(1,I)=VNLS(1,I)/RAD(1,I)
      VNLP(1,I)=VNLP(1,I)/RAD(1,I)
      ENDDO

C--------------------------------------
C      Read number of grid points                
C      dummy                          
C      delta r in log scale from HD.DAT          ! read from HD.DAT
C--------------------------------------

      READ(40,*) NREAD(1),DUMMY,DRLOG(1)
      DRLOG(1) = DRLOG(1) * 3.D0                 ! true grid spacing in log-mesh     
C
C----- Read grid points in real space -----
C
      K = 1
      DO I=1,INT(NREAD(1)/4)
      READ(40,*) RAD(1,K), RAD(1,K+1), RAD(1,K+2), RAD(1,K+3)
      K = K + 4
      ENDDO
      READ(40,*) RAD(1,K)
C
C----- Read Vld -----                            ! Vlocd
C                                                
      K = 1
      DO I=1,INT(NREAD(1)/4)
      READ(40,*) VLD(1,K), VLD(1,K+1), VLD(1,K+2), VLD(1,K+3)
      K = K + 4
      ENDDO
      READ(40,*) VLD(1,K)
C
C----- Read Vldc -----                           ! Vlocd - Vcore
C                                                
      K = 1
      DO I=1,INT(NREAD(1)/4)
      READ(40,*) VLDC(1,K), VLDC(1,K+1), VLDC(1,K+2), VLDC(1,K+3)
      K = K + 4
      ENDDO
      READ(40,*) VLDC(1,K)

C--------------------------------------
C      Read number of grid points 
C      dummy                           
C      delta r in log scale from CKB.DAT          ! read from CKB.DAT
C--------------------------------------

      READ(43,*) NREAD(6),DUMMY,DRLOG(6)
      DRLOG(6) = DRLOG(6) * 3.D0                  ! true grid spacing in log-mesh

C----- Read denominators -----

      READ(43,*) SVS(6),PVP(6)                    ! for s and p pseudowavefunction
C
C----- Read grid points in real space -----
C
      K = 1
      DO I=1,INT(NREAD(6)/4)
      READ(43,*) RAD(6,K), RAD(6,K+1), RAD(6,K+2), RAD(6,K+3)
      K = K + 4
      ENDDO
      READ(43,*) RAD(6,K)
C
C----- Read |Vnls|Ys> -----                      ! caution:Ys is multiplied by r
C                                                ! Vnls = Vs - Vd
      K = 1
      DO I=1,INT(NREAD(6)/4)
      READ(43,*) VNLS(6,K), VNLS(6,K+1), VNLS(6,K+2), VNLS(6,K+3)
      K = K + 4
      ENDDO
      READ(43,*) VNLS(6,K)
C
C----- Read |Vnlp|Yp> -----                      ! caution:Yp is multiplied by r
C                                                ! Vnls = Vs - Vd
      K = 1
      DO I=1,INT(NREAD(6)/4)
      READ(43,*) VNLP(6,K), VNLP(6,K+1), VNLP(6,K+2), VNLP(6,K+3)
      K = K + 4
      ENDDO
      READ(43,*) VNLP(6,K)

      DO I=1,NREAD(6)
      VNLS(6,I)=VNLS(6,I)/RAD(6,I)
      VNLP(6,I)=VNLP(6,I)/RAD(6,I)
      ENDDO

C--------------------------------------
C      Read number of grid points                
C      dummy                          
C      delta r in log scale from CD.DAT          ! read from CD.DAT
C--------------------------------------

      READ(42,*) NREAD(6),DUMMY,DRLOG(6)
      DRLOG(6) = DRLOG(6) * 3.D0                 ! true grid spacing in log-mesh     
C
C----- Read grid points in real space -----
C
      K = 1
      DO I=1,INT(NREAD(6)/4)
      READ(42,*) RAD(6,K), RAD(6,K+1), RAD(6,K+2), RAD(6,K+3)
      K = K + 4
      ENDDO
      READ(42,*) RAD(6,K)
C
C----- Read Vld -----                            ! Vlocd
C                                                
      K = 1
      DO I=1,INT(NREAD(6)/4)
      READ(42,*) VLD(6,K), VLD(6,K+1), VLD(6,K+2), VLD(6,K+3)
      K = K + 4
      ENDDO
      READ(42,*) VLD(6,K)
C
C----- Read Vldc -----                           ! Vlocd - Vcore
C                                                
      K = 1
      DO I=1,INT(NREAD(6)/4)
      READ(42,*) VLDC(6,K), VLDC(6,K+1), VLDC(6,K+2), VLDC(6,K+3)
      K = K + 4
      ENDDO
      READ(42,*) VLDC(6,K)
      
C--------------------------------------
C      Read number of grid points 
C      dummy                           
C      delta r in log scale from OKB.DAT          ! read from OKB.DAT
C--------------------------------------

      READ(45,*) NREAD(8),DUMMY,DRLOG(8)
      DRLOG(8) = DRLOG(8) * 3.D0                  ! true grid spacing in log-mesh

C----- Read denominators -----

      READ(45,*) SVS(8),PVP(8)                    ! for s and p pseudowavefunction
C
C----- Read grid points in real space -----
C
      K = 1
      DO I=1,INT(NREAD(8)/4)
      READ(45,*) RAD(8,K), RAD(8,K+1), RAD(8,K+2), RAD(8,K+3)
      K = K + 4
      ENDDO
      READ(45,*) RAD(8,K)
C
C----- Read |Vnls|Ys> -----                      ! caution:Ys is multiplied by r
C                                                ! Vnls = Vs - Vd
      K = 1
      DO I=1,INT(NREAD(8)/4)
      READ(45,*) VNLS(8,K), VNLS(8,K+1), VNLS(8,K+2), VNLS(8,K+3)
      K = K + 4
      ENDDO
      READ(45,*) VNLS(8,K)
C
C----- Read |Vnlp|Yp> -----                      ! caution:Yp is multiplied by r
C                                                ! Vnls = Vs - Vd
      K = 1
      DO I=1,INT(NREAD(8)/4)
      READ(45,*) VNLP(8,K), VNLP(8,K+1), VNLP(8,K+2), VNLP(8,K+3)
      K = K + 4
      ENDDO
      READ(45,*) VNLP(8,K)

      DO I=1,NREAD(8)
      VNLS(8,I)=VNLS(8,I)/RAD(8,I)
      VNLP(8,I)=VNLP(8,I)/RAD(8,I)
      ENDDO

C--------------------------------------
C      Read number of grid points                
C      dummy                          
C      delta r in log scale from OD.DAT          ! read from OD.DAT
C--------------------------------------

      READ(44,*) NREAD(8),DUMMY,DRLOG(8)
      DRLOG(8) = DRLOG(8) * 3.D0                 ! true grid spacing in log-mesh     
C
C----- Read grid points in real space -----
C
      K = 1
      DO I=1,INT(NREAD(8)/4)
      READ(44,*) RAD(8,K), RAD(8,K+1), RAD(8,K+2), RAD(8,K+3)
      K = K + 4
      ENDDO
      READ(44,*) RAD(8,K)
C
C----- Read Vld -----                            ! Vlocd
C                                                
      K = 1
      DO I=1,INT(NREAD(8)/4)
      READ(44,*) VLD(8,K), VLD(8,K+1), VLD(8,K+2), VLD(8,K+3)
      K = K + 4
      ENDDO
      READ(44,*) VLD(8,K)
C
C----- Read Vldc -----                           ! Vlocd - Vcore
C                                                
      K = 1
      DO I=1,INT(NREAD(8)/4)
      READ(44,*) VLDC(8,K), VLDC(8,K+1), VLDC(8,K+2), VLDC(8,K+3)
      K = K + 4
      ENDDO
      READ(44,*) VLDC(8,K)

      DO I=1,NREAD(8)   
      WRITE(46,*) RAD(1,I),VNLS(1,I),VNLP(1,I),VLD(1,I)
      WRITE(47,*) RAD(6,I),VNLS(6,I),VNLP(6,I),VLD(6,I)
      WRITE(48,*) RAD(8,I),VNLS(8,I),VNLP(8,I),VLD(8,I)
      ENDDO

C--------------------------------------
C      Read number of grid points
C      dummy
C      delta r in log scale from FKB.DAT         ! read from FKB.DAT
C--------------------------------------

      READ(67,*) NREAD(9),DUMMY,DRLOG(9),
     *                     SVS(9),PVP(9)         ! for s and p pseudowavefunction

      DRLOG(9) = DRLOG(9) * 3.D0                 ! true grid spacing in log-mesh

C
C----- Read grid points in real space -----
C
      READ(67,*) (RAD( 9,K),K=1,NREAD(9)),
     *           (VNLS(9,K),K=1,NREAD(9)),
     *           (VNLP(9,K),K=1,NREAD(9))

      DO I=1,NREAD(9)
      VNLS(9,I)=VNLS(9,I)/RAD(9,I)
      VNLP(9,I)=VNLP(9,I)/RAD(9,I)
      ENDDO

C--------------------------------------
C      Read number of grid points
C      dummy
C      delta r in log scale from FD.DAT          ! read from FD.DAT
C--------------------------------------

      READ(66,*) NREAD(9),DUMMY,DRLOG(9),
     *           (RAD(9,K), K=1,NREAD(9)),
     *           (VLD(9,K), K=1,NREAD(9)),
     *           (VLDC(9,K),K=1,NREAD(9))

      DRLOG(9) = DRLOG(9) * 3.D0                 ! true grid spacing in log-mesh

C-------------------------------------------
C      Read number of grid points 
C      dummy                           
C      delta r in log scale from CLKBHS.DAT       ! read from OKB.DAT
C-------------------------------------------

      READ(51,*) NREAD(17),DUMMY,DRLOG(17),
     *           SVS(17),PVP(17)                    ! for s and p pseudowavefunction
      DRLOG(17) = DRLOG(17) * 3.D0                  ! true grid spacing in log-mesh

C----- Read grid points in real space -----
C
      K = 1
      DO I=1,72
      READ(51,*) RAD(17,K), RAD(17,K+1), RAD(17,K+2)
      K = K + 3
      ENDDO
      DO I=1,51
      READ(51,*) RAD(17,K), RAD(17,K+1), RAD(17,K+2), RAD(17,K+3)
      K = K + 4
      ENDDO
      READ(51,*) RAD(17,K),VNLS(17,1), VNLS(17,2)
C
C----- Read |Vnls|Ys> -----                      ! caution:Ys is multiplied by r
C                                                ! Vnls = Vs - Vd
      K = 3
      DO I=1,56
      READ(51,*) VNLS(17,K), VNLS(17,K+1), VNLS(17,K+2) 
      K = K + 3
      ENDDO
      DO I=1,30
      READ(51,*) VNLS(17,K), VNLS(17,K+1), VNLS(17,K+2), VNLS(17,K+3)
      K = K + 4
      ENDDO
      DO I=1,43
      READ(51,*) VNLS(17,K), VNLS(17,K+1), VNLS(17,K+2) 
      K = K + 3
      ENDDO
      READ(51,*) VNLS(17,K),VNLS(17,K+1),VNLP(17,1)
C
C----- Read |Vnlp|Yp> -----                      ! caution:Yp is multiplied by r
C                                                ! Vnls = Vs - Vd
      K = 2
      DO I=1,75
      READ(51,*) VNLP(17,K), VNLP(17,K+1), VNLP(17,K+2)
      K = K + 3
      ENDDO
      DO I=1,17
      READ(51,*) VNLP(17,K), VNLP(17,K+1), VNLP(17,K+2), VNLP(17,K+3)
      K = K + 4
      ENDDO
      DO I=1,42
      READ(51,*) VNLP(17,K), VNLP(17,K+1), VNLP(17,K+2)
      K = K + 3
      ENDDO
      READ(51,*) VNLP(17,K)

      DO I=1,NREAD(17)
      VNLS(17,I)=VNLS(17,I)/RAD(17,I)
      VNLP(17,I)=VNLP(17,I)/RAD(17,I)
      ENDDO

C-------------------------------------------
C      Read number of grid points                
C      dummy                          
C      delta r in log scale from CLDBHS.DAT      ! read from OD.DAT
C-------------------------------------------

      READ(50,*) NREAD(17),DUMMY,DRLOG(17),
     *           RAD(17,1),RAD(17,2)
      DRLOG(17) = DRLOG(17) * 3.D0               ! true grid spacing in log-mesh     
C
C----- Read grid points in real space -----
C
      K = 3
      DO I=1,72
      READ(50,*) RAD(17,K), RAD(17,K+1), RAD(17,K+2)
      K = K + 3
      ENDDO
      DO I=1,50
      READ(50,*) RAD(17,K), RAD(17,K+1), RAD(17,K+2), RAD(17,K+3)
      K = K + 4
      ENDDO
      READ(50,*) RAD(17,K), RAD(17,K+1), RAD(17,K+2), VLD(17,1)
C
C----- Read Vld -----                            ! Vlocd
C                                                
      K = 2
      DO I=1,105
      READ(50,*) VLD(17,K), VLD(17,K+1), VLD(17,K+2), VLD(17,K+3)
      K = K + 4
      ENDDO
C
C----- Read Vldc -----                           ! Vlocd - Vcore
C                                                
      K = 1
      DO I=1,75
      READ(50,*) VLDC(17,K), VLDC(17,K+1), VLDC(17,K+2), VLDC(17,K+3)
      K = K + 4
      ENDDO
      DO I=1,40
      READ(50,*) VLDC(17,K), VLDC(17,K+1), VLDC(17,K+2)
      K = K + 3
      ENDDO
      READ(50,*) VLDC(17,K)

      DO I=1,NREAD(17)   
      WRITE(52,*) RAD(17,I),VNLS(17,I),VNLP(17,I),VLD(17,I)
      ENDDO

C--------------------------------------
C      Read number of grid points 
C      dummy                           
C      delta r in log scale from NKB.DAT          ! read from OKB.DAT
C--------------------------------------

      READ(53,*) NREAD(7),DUMMY,DRLOG(7)
      DRLOG(7) = DRLOG(7) * 3.D0                  ! true grid spacing in log-mesh

C----- Read denominators -----

      READ(53,*) SVS(7),PVP(7)                    ! for s and p pseudowavefunction
C
C----- Read grid points in real space -----
C
      K = 1
      DO I=1,INT(NREAD(7)/4)
      READ(53,*) RAD(7,K), RAD(7,K+1), RAD(7,K+2), RAD(7,K+3)
      K = K + 4
      ENDDO
      READ(53,*) RAD(7,K)
C
C----- Read |Vnls|Ys> -----                      ! caution:Ys is multiplied by r
C                                                ! Vnls = Vs - Vd
      K = 1
      DO I=1,INT(NREAD(7)/4)
      READ(53,*) VNLS(7,K), VNLS(7,K+1), VNLS(7,K+2), VNLS(7,K+3)
      K = K + 4
      ENDDO
      READ(53,*) VNLS(7,K)
C
C----- Read |Vnlp|Yp> -----                      ! caution:Yp is multiplied by r
C                                                ! Vnls = Vs - Vd
      K = 1
      DO I=1,INT(NREAD(7)/4)
      READ(53,*) VNLP(7,K), VNLP(7,K+1), VNLP(7,K+2), VNLP(7,K+3)
      K = K + 4
      ENDDO
      READ(53,*) VNLP(7,K)

      DO I=1,NREAD(7)
      VNLS(7,I)=VNLS(7,I)/RAD(7,I)
      VNLP(7,I)=VNLP(7,I)/RAD(7,I)
      ENDDO

C--------------------------------------
C      Read number of grid points                
C      dummy                          
C      delta r in log scale from ND.DAT          ! read from OD.DAT
C--------------------------------------

      READ(49,*) NREAD(7),DUMMY,DRLOG(7)
      DRLOG(7) = DRLOG(7) * 3.D0                 ! true grid spacing in log-mesh     
C
C----- Read grid points in real space -----
C
      K = 1
      DO I=1,INT(NREAD(7)/4)
      READ(49,*) RAD(7,K), RAD(7,K+1), RAD(7,K+2), RAD(7,K+3)
      K = K + 4
      ENDDO
      READ(49,*) RAD(7,K)
C
C----- Read Vld -----                            ! Vlocd
C                                                
      K = 1
      DO I=1,INT(NREAD(7)/4)
      READ(49,*) VLD(7,K), VLD(7,K+1), VLD(7,K+2), VLD(7,K+3)
      K = K + 4
      ENDDO
      READ(49,*) VLD(7,K)
C
C----- Read Vldc -----                           ! Vlocd - Vcore
C                                                
      K = 1
      DO I=1,INT(NREAD(7)/4)
      READ(49,*) VLDC(7,K), VLDC(7,K+1), VLDC(7,K+2), VLDC(7,K+3)
      K = K + 4
      ENDDO
      READ(49,*) VLDC(7,K)

      DO I=1,NREAD(7)   
      WRITE(54,*) RAD(7,I),VNLS(7,I),VNLP(7,I),VLD(7,I)
      ENDDO

C--------------------------------------
C      Read number of grid points
C      dummy
C      delta r in log scale from PKB.DAT         ! read from PKB.DAT
C--------------------------------------

      READ(69,*) NREAD(15),DUMMY,DRLOG(15),
     *                     SVS(15),PVP(15)       ! for s and p pseudowavefunction

      DRLOG(15) = DRLOG(15) * 3.D0               ! true grid spacing in log-mesh

      READ(69,*) (RAD( 15,K),K=1,NREAD(15)),
     *           (VNLS(15,K),K=1,NREAD(15)),
     *           (VNLP(15,K),K=1,NREAD(15))

C
C----- Read grid points in real space -----
C
C     READ(69,*) (RAD(15,K),K=1,NREAD(15))

C
C----- Read |Vnls|Ys> -----                      ! caution:Ys is multiplied by r
C                                                ! Vnls = Vs - Vd
C     READ(69,*) (VNLS(15,K),K=1,NREAD(15))

C
C----- Read |Vnlp|Yp> -----                      ! caution:Yp is multiplied by r
C                                                ! Vnls = Vs - Vd
C     READ(69,*) (VNLP(15,K),K=1,NREAD(15))

      DO I=1,NREAD(15)
      VNLS(15,I)=VNLS(15,I)/RAD(15,I)
      VNLP(15,I)=VNLP(15,I)/RAD(15,I)
      ENDDO

C--------------------------------------
C      Read number of grid points
C      dummy
C      delta r in log scale from PD.DAT          ! read from PD.DAT
C--------------------------------------

      READ(68,*) NREAD(15),DUMMY,DRLOG(15),
     *           (RAD(15,K), K=1,NREAD(15)),
     *           (VLD(15,K), K=1,NREAD(15)),
     *           (VLDC(15,K),K=1,NREAD(15))

      DRLOG(15) = DRLOG(15) * 3.D0               ! true grid spacing in log-mesh


C--------------------------------------
C      Read number of grid points                 ! for s,p,d non-local 
C      dummy
C      delta r in log scale from CAKB.DAT         ! read from CAKB.DAT
C--------------------------------------

      READ(76,*) NREAD(20),DUMMY,DRLOG(20),
     *                SVS(20),PVP(20),DVD(20)    ! for s,p, and d pseudowavefunction

      DRLOG(20) = DRLOG(20) * 3.D0               ! true grid spacing in log-mesh

      READ(76,*) (RAD( 20,K),K=1,NREAD(20)),
     *           (VNLS(20,K),K=1,NREAD(20)),
     *           (VNLP(20,K),K=1,NREAD(20)),
     *           (VNLD(20,K),K=1,NREAD(20))

C
C----- Read |Vnlp|Yp> -----                      ! caution:Yp is multiplied by r
C                                                ! Vnls = Vs - Vd

      DO I=1,NREAD(20)
      VNLS(20,I)=VNLS(20,I)/RAD(20,I)
      VNLP(20,I)=VNLP(20,I)/RAD(20,I)
      VNLD(20,I)=VNLD(20,I)/RAD(20,I)
      ENDDO

      DO I=1,NREAD(20)   
      WRITE(57,'(4F15.8)') RAD(20,I),VNLS(20,I),VNLP(20,I),VNLD(20,I)
      ENDDO

C--------------------------------------
C      Read number of grid points                 ! for s,p,d non-local 
C      dummy
C      delta r in log scale from CAKB.DAT         ! read from CAKB.DAT
C--------------------------------------

      READ(77,*) NREAD(25),DUMMY,DRLOG(25),
     *                SVS(25),PVP(25),DVD(25)    ! for s,p, and d pseudowavefunction

      DRLOG(25) = DRLOG(25) * 3.D0               ! true grid spacing in log-mesh

      READ(77,*) (RAD( 25,K),K=1,NREAD(25)),
     *           (VNLS(25,K),K=1,NREAD(25)),
     *           (VNLP(25,K),K=1,NREAD(25)),
     *           (VNLD(25,K),K=1,NREAD(25))

C
C----- Read |Vnlp|Yp> -----                      ! caution:Yp is multiplied by r
C                                                ! Vnls = Vs - Vd

      DO I=1,NREAD(20)
      VNLS(25,I)=VNLS(25,I)/RAD(25,I)
      VNLP(25,I)=VNLP(25,I)/RAD(25,I)
      VNLD(25,I)=VNLD(25,I)/RAD(25,I)
      ENDDO

      DO I=1,NREAD(25)   
      WRITE(58,'(4F15.8)') RAD(25,I),VNLS(25,I),VNLP(25,I),VNLD(25,I)
      ENDDO

      RETURN 
      END

C--------------------------------------------------------
C   SUBROUTINE Gnalm  
C   Coefficients of nonlocal Pseudo-potential for
C   Kleinman and Bylander form                         
C   Phys. Rev. Lett. 48, 1425( 1982 )
C   Phys. Rev. B.    50, 12234(1994 )
C   NCPS 97
C--------------------------------------------------------

      SUBROUTINE GNALM( RWF,TGNALM,DGNALM,NOR )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol
      include 'nlocd.i'                         ! Vmol   non-local d

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NORA = 49 )
C     PARAMETER ( NOCC =  1 )
C     PARAMETER ( NNUC = 36 )
C     PARAMETER ( RCUT =  2.5D0)
C     PARAMETER ( NSL  =  3000 )

      PARAMETER ( PI   = 3.14159265358979323D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ
      COMMON / ACELL / XL, YL, ZL
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)
      COMMON / DBLG  / WNLOCS( NNUC,NSL ),WNLOCPX( NNUC,NSL ),
     *                 WNLOCPY(NNUC,NSL ),WNLOCPZ( NNUC,NSL )

C     DIMENSION RWF( NMAX,NMAX,NMAX,NORA )
      DIMENSION RWF( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NOR )

      DIMENSION TGNALM(  NNUC,2,3,NORA )          ! num of elctrns,atoms,angular momentum quantum numbers
      DIMENSION TGNALM1( NNUC,2,3,NORA )          ! num of elctrns,atoms,angular momentum quantum numbers

      DIMENSION DGNALM(  NNUC,5,NORA )            ! num of elctrns,atoms,angular momentum quantum numbers
      DIMENSION DGNALM1( NNUC,5,NORA )            ! num of elctrns,atoms,angular momentum quantum numbers

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      COM  = 1.D0 / DSQRT( 4.D0*PI )
      DV   = DX*DY*DZ

C----- initialize -----

      DO 1 I =1,NOR
      DO 1 M =1,3
      DO 1 L =1,2
      DO 1 NA=1,NNUC
      TGNALM1( NA,L,M,I ) = 0.D0
1     CONTINUE
      DGNALM1( :,:,: ) = 0.D0

C----- calculate < Y(r)Vps | rwf > -----

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO 10 NA=1, NNUC
        NUMZ = INT( ZA( NA ) )
        NCOUNT = 0
        NJ = NPD(NA)

      SELECT CASE (NJ)

      CASE ( 0 )    ! NJ = 0: for s and p non-local

      DO 20 K=1, JZ
        DO 30 J=1, JY
          DO 40 I=1, JX

          RX = DX*( I-NNX ) - PNUC(NA,1) 
          RY = DY*( J-NNY ) - PNUC(NA,2)
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)

          R2 = RX**2+RY**2+RZ**2
          R = DSQRT(R2)

          IF( R .LT. RCUT ) THEN
          NCOUNT = NCOUNT + 1

          DO 50 N=1, NOR

C           ----- for s-component -----

            TGNALM1( NA,1,1,N ) = TGNALM1( NA,1,1,N ) +
     *      COM * WNLOCS(NA,NCOUNT)*RWF( I,J,K,N )*DV
     *      / ( SVS(NUMZ) ) 

C           ----- for p-component -----

            TGNALM1( NA,2,1,N ) = TGNALM1( NA,2,1,N )               ! Px component
     *      +  COM * WNLOCPX(NA,NCOUNT)*RWF( I,J,K,N )*DV
     *      / ( PVP(NUMZ) ) 

            TGNALM1( NA,2,2,N ) = TGNALM1( NA,2,2,N )               ! Py component
     *      +  COM * WNLOCPY(NA,NCOUNT)*RWF( I,J,K,N )*DV
     *      / ( PVP(NUMZ) ) 

            TGNALM1( NA,2,3,N ) = TGNALM1( NA,2,3,N )               ! Pz component
     *      +  COM * WNLOCPZ(NA,NCOUNT)*RWF( I,J,K,N )*DV
     *      / ( PVP(NUMZ) ) 

50        CONTINUE

          ENDIF

40        CONTINUE
30      CONTINUE
20    CONTINUE

      CASE ( 1 )    ! NJ = 1: for s,p and d non-local

      DO 25 K=1, JZ
        DO 35 J=1, JY
          DO 45 I=1, JX

          RX = DX*( I-NNX ) - PNUC(NA,1) 
          RY = DY*( J-NNY ) - PNUC(NA,2)
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)

          R2 = RX**2+RY**2+RZ**2
          R = DSQRT(R2)

          IF( R .LT. RCUT ) THEN
          NCOUNT = NCOUNT + 1

          DO 55 N=1, NOR

C           ----- for s-component -----

            TGNALM1( NA,1,1,N ) = TGNALM1( NA,1,1,N ) +
     *      COM * WNLOCS(NA,NCOUNT)*RWF( I,J,K,N )*DV
     *      / ( SVS(NUMZ) ) 

C           ----- for p-component -----

            TGNALM1( NA,2,1,N ) = TGNALM1( NA,2,1,N )               ! Px component
     *      +  COM * WNLOCPX(NA,NCOUNT)*RWF( I,J,K,N )*DV
     *      / ( PVP(NUMZ) ) 

            TGNALM1( NA,2,2,N ) = TGNALM1( NA,2,2,N )               ! Py component
     *      +  COM * WNLOCPY(NA,NCOUNT)*RWF( I,J,K,N )*DV
     *      / ( PVP(NUMZ) ) 

            TGNALM1( NA,2,3,N ) = TGNALM1( NA,2,3,N )               ! Pz component
     *      +  COM * WNLOCPZ(NA,NCOUNT)*RWF( I,J,K,N )*DV
     *      / ( PVP(NUMZ) ) 

C           ----- for d-component -----

            DGNALM1( NA,1,N ) = DGNALM1( NA,1,N )               ! dxy component
     *      +  COM * WNLOCDXY(NA,NCOUNT)*RWF( I,J,K,N )*DV
     *      / ( DVD(NUMZ) ) 

            DGNALM1( NA,2,N ) = DGNALM1( NA,2,N )               ! dyz component
     *      +  COM * WNLOCDYZ(NA,NCOUNT)*RWF( I,J,K,N )*DV
     *      / ( DVD(NUMZ) ) 

            DGNALM1( NA,3,N ) = DGNALM1( NA,3,N )               ! dzx component
     *      +  COM * WNLOCDZX(NA,NCOUNT)*RWF( I,J,K,N )*DV
     *      / ( DVD(NUMZ) ) 

            DGNALM1( NA,4,N ) = DGNALM1( NA,4,N )               ! dz2 component
     *      +  COM * WNLOCDZ2(NA,NCOUNT)*RWF( I,J,K,N )*DV
     *      / ( DVD(NUMZ) ) 

            DGNALM1( NA,5,N ) = DGNALM1( NA,5,N )               ! dx2-y2 component
     *      +  COM * WNLOCDX2(NA,NCOUNT)*RWF( I,J,K,N )*DV
     *      / ( DVD(NUMZ) ) 

55        CONTINUE

          ENDIF

45        CONTINUE
35      CONTINUE
25    CONTINUE

      END SELECT

10    CONTINUE

      NUMDAT = NNUC*6*NOR
      CALL MPI_REDUCE(TGNALM1(1,1,1,1),TGNALM(1,1,1,1),NUMDAT,
     *          MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(TGNALM(1,1,1,1),NUMDAT,
     *          MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

      NUMDAT = NNUC*5*NOR
      CALL MPI_REDUCE(DGNALM1(1,1,1),DGNALM(1,1,1),NUMDAT,
     *          MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DGNALM(1,1,1),NUMDAT,
     *          MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

C     WRITE(*,*) 'rank ',MYID  ,'TGNALM =',TGNALM(5,1,1,5)

      RETURN
      END

C--------------------------------------------------------
C   SUBROUTINE NUCLEAR POTENSHAL FOR NON-PERIODIC SYSTEM
C--------------------------------------------------------
   
      SUBROUTINE NPOT

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC =  36 )

      COMMON / MPIFRC / FRCM( NNUC,NDIM )

      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / FRCE  / FRC( NNUC,NDIM )
C     COMMON / LATTICE / POT
      COMMON / LTC     / POT                    ! Vmol

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      POT = 0.D0

      IF(MYID.EQ.0) THEN
      DO I = 1,NNUC-1
      DO J = I+1,NNUC

      RX = PNUC( I,1 ) - PNUC( J,1 ) 
      RY = PNUC( I,2 ) - PNUC( J,2 ) 
      RZ = PNUC( I,3 ) - PNUC( J,3 ) 
      R2 = RX**2 + RY**2 + RZ**2
      R  = DSQRT( R2 )

      POT = POT + ZVAL( I ) * ZVAL( J ) / R

      ENDDO
      ENDDO

      DO NA = 1,NNUC
      DO I = 1,NNUC

      IF( NA .NE. I ) THEN

      RX = PNUC( NA,1 ) - PNUC( I,1 ) 
      RY = PNUC( NA,2 ) - PNUC( I,2 ) 
      RZ = PNUC( NA,3 ) - PNUC( I,3 ) 
      R2 = RX**2 + RY**2 + RZ**2
      R  = DSQRT( R2 )

      FRCM(NA,1) = FRCM(NA,1) + ZVAL( NA ) * ZVAL( I )*RX / (R2*R)
      FRCM(NA,2) = FRCM(NA,2) + ZVAL( NA ) * ZVAL( I )*RY / (R2*R)
      FRCM(NA,3) = FRCM(NA,3) + ZVAL( NA ) * ZVAL( I )*RZ / (R2*R)

      ENDIF

      ENDDO
      ENDDO

C     IF(MYID.EQ.0) THEN
      WRITE(*,*) 'Nuclear Repulsion Energy [a.u.]'
      WRITE(*,*)  POT
      WRITE(*,*)
      ENDIF

      RETURN
      END


C ---------------------------------------------------
      SUBROUTINE PUTPS
C ---------------------------------------------------

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'QMpara.i'                       ! Vmol
      include 'nlocd.i'                        ! Vmol

C     PARAMETER ( NSL  = 3000  )

      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),          ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)

      DIMENSION A(2108)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      DO 20 K = 1,100

      NATOM = K

      DO NA = 1,NNUC
       NUMZ = INT( ZA(NA) )
       IF(NUMZ.EQ.NATOM) GOTO 10
      ENDDO

      GOTO 20

 10   CONTINUE

      NJ = NPD(NA)

      SELECT CASE (NJ)

      CASE ( 0 )    ! NJ = 0: for s and p non-local

      A(:) = 0.D0

      IF(MYID.EQ.0)THEN

      A(1)=DRLOG(NATOM)
      A(2)=SVS(NATOM)
      A(3)=PVP(NATOM)

      NCOUNT=3

      DO I=1,421
       NCOUNT=NCOUNT+1
       A(NCOUNT)=RAD(NATOM,I)
      ENDDO

      DO I=1,421
       NCOUNT=NCOUNT+1
       A(NCOUNT)=VNLS(NATOM,I)
      ENDDO

      DO I=1,421
       NCOUNT=NCOUNT+1
       A(NCOUNT)=VNLP(NATOM,I)
      ENDDO

      DO I=1,421
       NCOUNT=NCOUNT+1
       A(NCOUNT)=VLD(NATOM,I)
      ENDDO

      DO I=1,421
       NCOUNT=NCOUNT+1
       A(NCOUNT)=VLDC(NATOM,I)
      ENDDO

      ENDIF

      CALL MPI_BCAST(A(1),2108,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)

      IF(MYID.NE.0) THEN

      DRLOG(NATOM)=A(1)
      SVS(NATOM)=A(2)
      PVP(NATOM)=A(3)

      NCOUNT=3

      DO I=1,421
       NCOUNT=NCOUNT+1
       RAD(NATOM,I)=A(NCOUNT)
      ENDDO

      DO I=1,421
       NCOUNT=NCOUNT+1
       VNLS(NATOM,I)=A(NCOUNT)
      ENDDO

      DO I=1,421
       NCOUNT=NCOUNT+1
       VNLP(NATOM,I)=A(NCOUNT)
      ENDDO

      DO I=1,421
       NCOUNT=NCOUNT+1
       VLD(NATOM,I)=A(NCOUNT)
      ENDDO

      DO I=1,421
       NCOUNT=NCOUNT+1
       VLDC(NATOM,I)=A(NCOUNT)
      ENDDO

      ENDIF

      CASE ( 1 )    ! NJ = 1: for s, p and d non-local

      A(:) = 0.D0

      IF(MYID.EQ.0)THEN

      A(1)=DRLOG(NATOM)
      A(2)=SVS(  NATOM)
      A(3)=PVP(  NATOM)
      A(4)=DVD(  NATOM)

      NCOUNT=4

      DO I=1,421
       NCOUNT=NCOUNT+1
       A(NCOUNT)=RAD(NATOM,I)
      ENDDO

      DO I=1,421
       NCOUNT=NCOUNT+1
       A(NCOUNT)=VNLS(NATOM,I)
      ENDDO

      DO I=1,421
       NCOUNT=NCOUNT+1
       A(NCOUNT)=VNLP(NATOM,I)
      ENDDO

      DO I=1,421
       NCOUNT=NCOUNT+1
       A(NCOUNT)=VNLD(NATOM,I)
      ENDDO

      ENDIF

      CALL MPI_BCAST(A(1),1688,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)

      IF(MYID.NE.0) THEN

      DRLOG(NATOM)=A(1)
      SVS(  NATOM)=A(2)
      PVP(  NATOM)=A(3)
      DVD(  NATOM)=A(4)

      NCOUNT=4

      DO I=1,421
       NCOUNT=NCOUNT+1
       RAD(NATOM,I)=A(NCOUNT)
      ENDDO

      DO I=1,421
       NCOUNT=NCOUNT+1
       VNLS(NATOM,I)=A(NCOUNT)
      ENDDO

      DO I=1,421
       NCOUNT=NCOUNT+1
       VNLP(NATOM,I)=A(NCOUNT)
      ENDDO

      DO I=1,421
       NCOUNT=NCOUNT+1
       VNLD(NATOM,I)=A(NCOUNT)
      ENDDO

      ENDIF

      END SELECT 

 20   CONTINUE

      RETURN
      END


C--------------------------------------------------------
C   SUBROUTINE Hellmann Feynman Force  
C   Force from electrons and nuclei              
C   Kleinman and Bylander form                         
C   Phys. Rev. B. 50, 12234( 1994 )
C   NCPS 97
C--------------------------------------------------------

      SUBROUTINE LFRC( RHO ) ! local part for non-periodic system

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol
      include 'nlocd.i'

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NSL  = 3000 )

      PARAMETER ( PI   = 3.14159265358979323D0 )
C     PARAMETER ( RCUT = 2.5D0      )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / MPIFRC / FRCM( NNUC,NDIM )

      COMMON / FRCE  / FRC( NNUC,NDIM )
      COMMON / GRID2 / DX, DY, DZ
      COMMON / NCLR  / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PLOC / VLOC( NNUC,NSL )
      COMMON / BHS / C1(100),C2(100),AL1(100),AL2(100)

      DIMENSION DWFX( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFY( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFZ( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHO(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      DV  = DX*DY*DZ

C
C----- contribution from local part -----
C

      CALL DWF2( RHO,DWFX,DWFY,DWFZ )                   ! derivative of density 
     
      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO 1000 NA=1,NNUC

      NZA = INT(ZA(NA))

      NCOUNT = 0

      DO 20 K=1,JZ
      DO 30 J=1,JY
      DO 40 I=1,JX
        
          RX = DX*( I-NNX ) - PNUC(NA,1) 
          RY = DY*( J-NNY ) - PNUC(NA,2)
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)
          R2 = RX**2 + RY**2 + RZ**2
          R  = DSQRT(R2) 

          IF( R .LT. RCUT ) THEN
          NCOUNT = NCOUNT + 1

          FRCM(NA,1) = FRCM(NA,1) - DWFX(I,J,K)
     *                    * VLOC(NA,NCOUNT) * DV
          FRCM(NA,2) = FRCM(NA,2) - DWFY(I,J,K)
     *                    * VLOC(NA,NCOUNT) * DV
          FRCM(NA,3) = FRCM(NA,3) - DWFZ(I,J,K)
     *                    * VLOC(NA,NCOUNT) * DV

          ELSE

          VT = -ZVAL(NA)*( C1(NZA)*ERF( DSQRT(AL1(NZA)*R2) )                  ! analytical function of BHS
     *       +             C2(NZA)*ERF( DSQRT(AL2(NZA)*R2) ) ) 
     *       /  R           
            
          FRCM(NA,1) = FRCM(NA,1) - DWFX(I,J,K)
     *                    * VT * DV
          FRCM(NA,2) = FRCM(NA,2) - DWFY(I,J,K)
     *                    * VT * DV
          FRCM(NA,3) = FRCM(NA,3) - DWFZ(I,J,K)
     *                    * VT * DV

          ENDIF

40    CONTINUE
30    CONTINUE
20    CONTINUE

1000  CONTINUE

      RETURN
      END

      SUBROUTINE NLFRC( RWF,GNALM,DGNALM,MOR,BF ) ! non-local part

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                       ! Vmol
      include 'nlocd.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NORA =  49 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NSL  = 3000 )

      PARAMETER ( PI   = 3.14159265358979323D0 )
C     PARAMETER ( RCUT = 2.5D0      )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / MPIFRC / FRCM( NNUC,NDIM )

      COMMON / FRCE  / FRC( NNUC,NDIM )
      COMMON / GRID2 / DX, DY, DZ
      COMMON / NCLR  / PNUC( NNUC,NDIM )
      COMMON / DBLG  / WNLOCS( NNUC,NSL ),WNLOCPX( NNUC,NSL ),
     *                 WNLOCPY(NNUC,NSL ),WNLOCPZ( NNUC,NSL )

      DIMENSION RWF(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION DWFX( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFY( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFZ( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION TWF(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      DIMENSION GNALM(  NNUC,2,3,NORA )             ! num of elctrns,atoms,angular momentum quantum numbers
      DIMENSION DGNALM( NNUC,  5,NORA )             ! num of elctrns,atoms,angular momentum quantum numbers
      DIMENSION FRNL(  2,3,3 ) 
      DIMENSION FRNLD( 5,3 )       ! non-local d 

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      DV     = DX*DY*DZ
      COM    = 1.D0 / DSQRT( 4.D0*PI )
      COM3   =  3.D0*COM
      COM15  = 15.D0*COM
      COM54  =  5.0D0/4.D0*COM
      COM154 = 15.0D0/4.D0*COM
      CF     = 2.D0*BF

C
C----- contribution from non-local part -----
C

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO 2000 NA=1,NNUC

      NJ = NPD(NA) 

      SELECT CASE (NJ)  

      CASE ( 0 )    ! NJ = 0: for s and p non-local
      
      DO 90 N=1,MOR                                   ! loop over alpha orbitals
      
      DO 100 K=1,JZ
      DO 100 J=1,JY
      DO 100 I=1,JX
      TWF( I,J,K ) = RWF( I,J,K,N )
100   CONTINUE

      DO K=1,3
      DO J=1,3
      DO I=1,2
      FRNL( I,J,K ) = 0.D0                            ! initialize force
      ENDDO
      ENDDO
      ENDDO
      
      CALL DWF2( TWF,DWFX,DWFY,DWFZ )                 ! derivative of alpha wavefunction
            
      NCOUNT = 0                                      ! reset counter
      DO 110 K=1,JZ
      DO 120 J=1,JY
      DO 130 I=1,JX
            
        RX = DX*( I-NNX ) - PNUC(NA,1) 
        RY = DY*( J-NNY ) - PNUC(NA,2)
        RZ = DZ*( K-NNZ ) - PNUC(NA,3)

        R2 = RX**2 + RY**2 + RZ**2
        R = DSQRT(R2)

        IF( R .LT. RCUT ) THEN     
        NCOUNT = NCOUNT + 1

C----- for s-component -----

        FRNL( 1,1,1 ) = FRNL( 1,1,1 ) +
     *    COM*WNLOCS(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNL( 1,1,2 ) = FRNL( 1,1,2 ) +
     *    COM*WNLOCS(NA,NCOUNT)*DWFY( I,J,K )*DV
        FRNL( 1,1,3 ) = FRNL( 1,1,3 ) +
     *    COM*WNLOCS(NA,NCOUNT)*DWFZ( I,J,K )*DV

C----- for p-component -----

        FRNL( 2,1,1 ) = FRNL( 2,1,1 ) +                 ! dx 
     *    COM3*WNLOCPX(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNL( 2,2,1 ) = FRNL( 2,2,1 ) +                 
     *    COM3*WNLOCPY(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNL( 2,3,1 ) = FRNL( 2,3,1 ) +                 
     *    COM3*WNLOCPZ(NA,NCOUNT)*DWFX( I,J,K )*DV 

        FRNL( 2,1,2 ) = FRNL( 2,1,2 ) +                 ! dy 
     *    COM3*WNLOCPX(NA,NCOUNT)*DWFY( I,J,K )*DV 
        FRNL( 2,2,2 ) = FRNL( 2,2,2 ) +                 
     *    COM3*WNLOCPY(NA,NCOUNT)*DWFY( I,J,K )*DV 
        FRNL( 2,3,2 ) = FRNL( 2,3,2 ) +                 
     *    COM3*WNLOCPZ(NA,NCOUNT)*DWFY( I,J,K )*DV 

        FRNL( 2,1,3 ) = FRNL( 2,1,3 ) +                 ! dz 
     *    COM3*WNLOCPX(NA,NCOUNT)*DWFZ( I,J,K )*DV 
        FRNL( 2,2,3 ) = FRNL( 2,2,3 ) +                 
     *    COM3*WNLOCPY(NA,NCOUNT)*DWFZ( I,J,K )*DV 
        FRNL( 2,3,3 ) = FRNL( 2,3,3 ) +                 
     *    COM3*WNLOCPZ(NA,NCOUNT)*DWFZ( I,J,K )*DV 

        ENDIF

130   CONTINUE
120   CONTINUE
110   CONTINUE

      FRCM(NA,1) = FRCM(NA,1) - CF*( FRNL(1,1,1)*GNALM(NA,1,1,N)
     *          + FRNL(2,1,1)*GNALM(NA,2,1,N) 
     *          + FRNL(2,2,1)*GNALM(NA,2,2,N)
     *          + FRNL(2,3,1)*GNALM(NA,2,3,N) )
      FRCM(NA,2) = FRCM(NA,2) - CF*( FRNL(1,1,2)*GNALM(NA,1,1,N)
     *          + FRNL(2,1,2)*GNALM(NA,2,1,N) 
     *          + FRNL(2,2,2)*GNALM(NA,2,2,N)
     *          + FRNL(2,3,2)*GNALM(NA,2,3,N) )  
      FRCM(NA,3) = FRCM(NA,3) - CF*( FRNL(1,1,3)*GNALM(NA,1,1,N)
     *          + FRNL(2,1,3)*GNALM(NA,2,1,N) 
     *          + FRNL(2,2,3)*GNALM(NA,2,2,N)
     *          + FRNL(2,3,3)*GNALM(NA,2,3,N) )  

90    CONTINUE
  
      CASE ( 1 )    ! NJ = 1: for s, p, and d non-local

      DO 95 N=1,MOR                                   ! loop over alpha orbitals
      
      DO 105 K=1,JZ
      DO 105 J=1,JY
      DO 105 I=1,JX
      TWF( I,J,K ) = RWF( I,J,K,N )
105   CONTINUE

      DO K=1,3
      DO J=1,3
      DO I=1,2
      FRNL( I,J,K ) = 0.D0                            ! initialize force
      ENDDO
      ENDDO
      ENDDO

      FRNLD( :,: ) = 0.D0                             ! non-local d 
      
      CALL DWF2( TWF,DWFX,DWFY,DWFZ )                 ! derivative of alpha wavefunction
            
      NCOUNT = 0                                      ! reset counter
      DO 115 K=1,JZ
      DO 125 J=1,JY
      DO 135 I=1,JX
            
        RX = DX*( I-NNX ) - PNUC(NA,1) 
        RY = DY*( J-NNY ) - PNUC(NA,2)
        RZ = DZ*( K-NNZ ) - PNUC(NA,3)

        R2 = RX**2 + RY**2 + RZ**2
        R = DSQRT(R2)

        IF( R .LT. RCUT ) THEN     
        NCOUNT = NCOUNT + 1

C----- for s-component -----

        FRNL( 1,1,1 ) = FRNL( 1,1,1 ) +
     *    COM*WNLOCS(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNL( 1,1,2 ) = FRNL( 1,1,2 ) +
     *    COM*WNLOCS(NA,NCOUNT)*DWFY( I,J,K )*DV
        FRNL( 1,1,3 ) = FRNL( 1,1,3 ) +
     *    COM*WNLOCS(NA,NCOUNT)*DWFZ( I,J,K )*DV

C----- for p-component -----

        FRNL( 2,1,1 ) = FRNL( 2,1,1 ) +                 ! dx 
     *    COM3*WNLOCPX(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNL( 2,2,1 ) = FRNL( 2,2,1 ) +                 
     *    COM3*WNLOCPY(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNL( 2,3,1 ) = FRNL( 2,3,1 ) +                 
     *    COM3*WNLOCPZ(NA,NCOUNT)*DWFX( I,J,K )*DV 

        FRNL( 2,1,2 ) = FRNL( 2,1,2 ) +                 ! dy 
     *    COM3*WNLOCPX(NA,NCOUNT)*DWFY( I,J,K )*DV 
        FRNL( 2,2,2 ) = FRNL( 2,2,2 ) +                 
     *    COM3*WNLOCPY(NA,NCOUNT)*DWFY( I,J,K )*DV 
        FRNL( 2,3,2 ) = FRNL( 2,3,2 ) +                 
     *    COM3*WNLOCPZ(NA,NCOUNT)*DWFY( I,J,K )*DV 

        FRNL( 2,1,3 ) = FRNL( 2,1,3 ) +                 ! dz 
     *    COM3*WNLOCPX(NA,NCOUNT)*DWFZ( I,J,K )*DV 
        FRNL( 2,2,3 ) = FRNL( 2,2,3 ) +                 
     *    COM3*WNLOCPY(NA,NCOUNT)*DWFZ( I,J,K )*DV 
        FRNL( 2,3,3 ) = FRNL( 2,3,3 ) +                 
     *    COM3*WNLOCPZ(NA,NCOUNT)*DWFZ( I,J,K )*DV 

C----- for d-component -----

        FRNLD( 1,1 ) = FRNLD( 1,1 ) +                 ! dxy
     *    COM15*WNLOCDXY(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNLD( 2,1 ) = FRNLD( 2,1 ) +                 ! dyz
     *    COM15*WNLOCDYZ(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNLD( 3,1 ) = FRNLD( 3,1 ) +                 ! dzx
     *    COM15*WNLOCDZX(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNLD( 4,1 ) = FRNLD( 4,1 ) +                 ! dz2
     *    COM54*WNLOCDZ2(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNLD( 5,1 ) = FRNLD( 5,1 ) +                 ! dx2-y2
     *    COM154*WNLOCDX2(NA,NCOUNT)*DWFX( I,J,K )*DV 

        FRNLD( 1,2 ) = FRNLD( 1,2 ) +                 ! dxy
     *    COM15*WNLOCDXY(NA,NCOUNT)*DWFY( I,J,K )*DV 
        FRNLD( 2,2 ) = FRNLD( 2,2 ) +                 ! dyz
     *    COM15*WNLOCDYZ(NA,NCOUNT)*DWFY( I,J,K )*DV 
        FRNLD( 3,2 ) = FRNLD( 3,2 ) +                 ! dzx
     *    COM15*WNLOCDZX(NA,NCOUNT)*DWFY( I,J,K )*DV 
        FRNLD( 4,2 ) = FRNLD( 4,2 ) +                 ! dz2
     *    COM54*WNLOCDZ2(NA,NCOUNT)*DWFY( I,J,K )*DV 
        FRNLD( 5,2 ) = FRNLD( 5,2 ) +                 ! dx2-y2
     *    COM154*WNLOCDX2(NA,NCOUNT)*DWFY( I,J,K )*DV 

        FRNLD( 1,3 ) = FRNLD( 1,3 ) +                 ! dxy
     *    COM15*WNLOCDXY(NA,NCOUNT)*DWFZ( I,J,K )*DV 
        FRNLD( 2,3 ) = FRNLD( 2,3 ) +                 ! dyz
     *    COM15*WNLOCDYZ(NA,NCOUNT)*DWFZ( I,J,K )*DV 
        FRNLD( 3,3 ) = FRNLD( 3,3 ) +                 ! dzx
     *    COM15*WNLOCDZX(NA,NCOUNT)*DWFZ( I,J,K )*DV 
        FRNLD( 4,3 ) = FRNLD( 4,3 ) +                 ! dz2
     *    COM54*WNLOCDZ2(NA,NCOUNT)*DWFZ( I,J,K )*DV 
        FRNLD( 5,3 ) = FRNLD( 5,3 ) +                 ! dx2-y2
     *    COM154*WNLOCDX2(NA,NCOUNT)*DWFZ( I,J,K )*DV 

        ENDIF

135   CONTINUE
125   CONTINUE
115   CONTINUE

      FRCM(NA,1) = FRCM(NA,1) - CF*( FRNL(1,1,1)*GNALM(NA,1,1,N)
     *          + FRNL(2,1,1)*GNALM(NA,2,1,N) 
     *          + FRNL(2,2,1)*GNALM(NA,2,2,N)
     *          + FRNL(2,3,1)*GNALM(NA,2,3,N)  
     *          + FRNLD(1,1)*DGNALM(NA,1,N) 
     *          + FRNLD(2,1)*DGNALM(NA,2,N)
     *          + FRNLD(3,1)*DGNALM(NA,3,N)  
     *          + FRNLD(4,1)*DGNALM(NA,4,N)  
     *          + FRNLD(5,1)*DGNALM(NA,5,N) )

      FRCM(NA,2) = FRCM(NA,2) - CF*( FRNL(1,1,2)*GNALM(NA,1,1,N)
     *          + FRNL(2,1,2)*GNALM(NA,2,1,N) 
     *          + FRNL(2,2,2)*GNALM(NA,2,2,N)
     *          + FRNL(2,3,2)*GNALM(NA,2,3,N)    
     *          + FRNLD(1,2)*DGNALM(NA,1,N) 
     *          + FRNLD(2,2)*DGNALM(NA,2,N)
     *          + FRNLD(3,2)*DGNALM(NA,3,N)  
     *          + FRNLD(4,2)*DGNALM(NA,4,N)  
     *          + FRNLD(5,2)*DGNALM(NA,5,N) )

      FRCM(NA,3) = FRCM(NA,3) - CF*( FRNL(1,1,3)*GNALM(NA,1,1,N)
     *          + FRNL(2,1,3)*GNALM(NA,2,1,N) 
     *          + FRNL(2,2,3)*GNALM(NA,2,2,N)
     *          + FRNL(2,3,3)*GNALM(NA,2,3,N) 
     *          + FRNLD(1,3)*DGNALM(NA,1,N) 
     *          + FRNLD(2,3)*DGNALM(NA,2,N)
     *          + FRNLD(3,3)*DGNALM(NA,3,N)  
     *          + FRNLD(4,3)*DGNALM(NA,4,N)  
     *          + FRNLD(5,3)*DGNALM(NA,5,N) )

95    CONTINUE

      END SELECT

2000  CONTINUE

      NUMDAT = NNUC*NDIM
      CALL MPI_REDUCE( FRCM(1,1),FRC(1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

      RETURN
      END

      SUBROUTINE LFRC1( RHO ) ! local part for non-periodic system
                              ! For link atoms
      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
      include 'nlocd.i'

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NSL  = 3000 )

      PARAMETER ( PI   = 3.14159265358979323D0 )
C     PARAMETER ( RCUT = 2.5D0      )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / MPIFRC / FRCM( NNUC,NDIM )

      COMMON / FRCE  / FRC( NNUC,NDIM )
      COMMON / GRID2 / DX, DY, DZ
      COMMON / NCLR  / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PLOC / VLOC( NNUC,NSL )
      COMMON / BHS / C1(100),C2(100),AL1(100),AL2(100)

      DIMENSION DWFX( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFY( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFZ( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHO(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      DV = DX*DY*DZ
      ns = NNUC - NLINK + 1

C
C----- contribution from local part -----
C

      CALL DWF2( RHO,DWFX,DWFY,DWFZ )                   ! derivative of density 
     
      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO 1000 NA=ns,NNUC

      NZA = INT(ZA(NA))

      NCOUNT = 0

      DO 20 K=1,JZ
      DO 30 J=1,JY
      DO 40 I=1,JX
        
          RX = DX*( I-NNX ) - PNUC(NA,1) 
          RY = DY*( J-NNY ) - PNUC(NA,2)
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)
          R2 = RX**2 + RY**2 + RZ**2
          R  = DSQRT(R2) 

          IF( R .LT. RCUT ) THEN
          NCOUNT = NCOUNT + 1

          FRCM(NA,1) = FRCM(NA,1) - DWFX(I,J,K)
     *                    * VLOC(NA,NCOUNT) * DV
          FRCM(NA,2) = FRCM(NA,2) - DWFY(I,J,K)
     *                    * VLOC(NA,NCOUNT) * DV
          FRCM(NA,3) = FRCM(NA,3) - DWFZ(I,J,K)
     *                    * VLOC(NA,NCOUNT) * DV

          ELSE

          VT = -ZVAL(NA)*( C1(NZA)*ERF( DSQRT(AL1(NZA)*R2) )                  ! analytical function of BHS
     *       +             C2(NZA)*ERF( DSQRT(AL2(NZA)*R2) ) ) 
     *       /  R           
            
          FRCM(NA,1) = FRCM(NA,1) - DWFX(I,J,K)
     *                    * VT * DV
          FRCM(NA,2) = FRCM(NA,2) - DWFY(I,J,K)
     *                    * VT * DV
          FRCM(NA,3) = FRCM(NA,3) - DWFZ(I,J,K)
     *                    * VT * DV

          ENDIF

40    CONTINUE
30    CONTINUE
20    CONTINUE

1000  CONTINUE

      RETURN
      END

      SUBROUTINE NLFRC1( RWF,GNALM,MOR,BF ) ! non-local part
                                            ! For link atoms

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
      include 'nlocd.i'

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NORA =  49 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NSL  = 3000 )

      PARAMETER ( PI   = 3.14159265358979323D0 )
C     PARAMETER ( RCUT = 2.5D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / MPIFRC / FRCM( NNUC,NDIM )

      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / FRCE  / FRC( NNUC,NDIM )
      COMMON / GRID2 / DX, DY, DZ
      COMMON / NCLR  / PNUC( NNUC,NDIM )
      COMMON / DBLG  / WNLOCS( NNUC,NSL ),WNLOCPX( NNUC,NSL ),
     *                 WNLOCPY(NNUC,NSL ),WNLOCPZ( NNUC,NSL )

      DIMENSION RWF(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION DWFX( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFY( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFZ( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION TWF(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      DIMENSION GNALM( NNUC,2,3,NORA )             ! num of elctrns,atoms,angular momentum quantum numbers
      DIMENSION FRNL(   2,3,3 ) 

      DIMENSION AFRC( NLINK1*NDIM ) 
      DIMENSION BFRC( NLINK1*NDIM ) 

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      DV  = DX*DY*DZ
      COM = 1.D0 / DSQRT( 4.D0*PI )
      CF  = 2.D0*BF
      ns = NNUC - NLINK + 1

C
C----- contribution from non-local part -----
C

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO 2000 NA=ns,NNUC
      
      DO 90 N=1,MOR                                   ! loop over alpha orbitals
      
      DO 100 K=1,JZ
      DO 100 J=1,JY
      DO 100 I=1,JX
      TWF( I,J,K ) = RWF( I,J,K,N )
100   CONTINUE

      DO K=1,3
      DO J=1,3
      DO I=1,2
      FRNL( I,J,K ) = 0.D0                            ! initialize force
      ENDDO
      ENDDO
      ENDDO
      
      CALL DWF2( TWF,DWFX,DWFY,DWFZ )                 ! derivative of alpha wavefunction
            
      NCOUNT = 0                                      ! reset counter
      DO 110 K=1,JZ
      DO 120 J=1,JY
      DO 130 I=1,JX
            
        RX = DX*( I-NNX ) - PNUC(NA,1) 
        RY = DY*( J-NNY ) - PNUC(NA,2)
        RZ = DZ*( K-NNZ ) - PNUC(NA,3)

        R2 = RX**2 + RY**2 + RZ**2
        R = DSQRT(R2)

        IF( R .LT. RCUT ) THEN     
        NCOUNT = NCOUNT + 1

C----- for s-component -----

        FRNL( 1,1,1 ) = FRNL( 1,1,1 ) +
     *    COM*WNLOCS(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNL( 1,1,2 ) = FRNL( 1,1,2 ) +
     *    COM*WNLOCS(NA,NCOUNT)*DWFY( I,J,K )*DV
        FRNL( 1,1,3 ) = FRNL( 1,1,3 ) +
     *    COM*WNLOCS(NA,NCOUNT)*DWFZ( I,J,K )*DV

C----- for p-component -----

        FRNL( 2,1,1 ) = FRNL( 2,1,1 ) +                 ! dx 
     *    3.D0*COM*WNLOCPX(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNL( 2,2,1 ) = FRNL( 2,2,1 ) +                 
     *    3.D0*COM*WNLOCPY(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNL( 2,3,1 ) = FRNL( 2,3,1 ) +                 
     *    3.D0*COM*WNLOCPZ(NA,NCOUNT)*DWFX( I,J,K )*DV 

        FRNL( 2,1,2 ) = FRNL( 2,1,2 ) +                 ! dy 
     *    3.D0*COM*WNLOCPX(NA,NCOUNT)*DWFY( I,J,K )*DV 
        FRNL( 2,2,2 ) = FRNL( 2,2,2 ) +                 
     *    3.D0*COM*WNLOCPY(NA,NCOUNT)*DWFY( I,J,K )*DV 
        FRNL( 2,3,2 ) = FRNL( 2,3,2 ) +                 
     *    3.D0*COM*WNLOCPZ(NA,NCOUNT)*DWFY( I,J,K )*DV 

        FRNL( 2,1,3 ) = FRNL( 2,1,3 ) +                 ! dz 
     *    3.D0*COM*WNLOCPX(NA,NCOUNT)*DWFZ( I,J,K )*DV 
        FRNL( 2,2,3 ) = FRNL( 2,2,3 ) +                 
     *    3.D0*COM*WNLOCPY(NA,NCOUNT)*DWFZ( I,J,K )*DV 
        FRNL( 2,3,3 ) = FRNL( 2,3,3 ) +                 
     *    3.D0*COM*WNLOCPZ(NA,NCOUNT)*DWFZ( I,J,K )*DV 

        ENDIF

130   CONTINUE
120   CONTINUE
110   CONTINUE

      FRCM(NA,1) = FRCM(NA,1) - CF*( FRNL(1,1,1)*GNALM(NA,1,1,N)
     *          + FRNL(2,1,1)*GNALM(NA,2,1,N) 
     *          + FRNL(2,2,1)*GNALM(NA,2,2,N)
     *          + FRNL(2,3,1)*GNALM(NA,2,3,N) )
      FRCM(NA,2) = FRCM(NA,2) - CF*( FRNL(1,1,2)*GNALM(NA,1,1,N)
     *          + FRNL(2,1,2)*GNALM(NA,2,1,N) 
     *          + FRNL(2,2,2)*GNALM(NA,2,2,N)
     *          + FRNL(2,3,2)*GNALM(NA,2,3,N) )  
      FRCM(NA,3) = FRCM(NA,3) - CF*( FRNL(1,1,3)*GNALM(NA,1,1,N)
     *          + FRNL(2,1,3)*GNALM(NA,2,1,N) 
     *          + FRNL(2,2,3)*GNALM(NA,2,2,N)
     *          + FRNL(2,3,3)*GNALM(NA,2,3,N) )  

90    CONTINUE
  
2000  CONTINUE

      NCOUNT = 0
      DO J = 1, NDIM
      DO I = ns, NNUC
        NCOUNT = NCOUNT+1
        AFRC( NCOUNT ) = FRCM(I,J)
      ENDDO
      ENDDO

      NUMDAT = NLINK*NDIM
      CALL MPI_REDUCE( AFRC(1),BFRC(1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

      NCOUNT = 0
      DO J = 1, NDIM
      DO I = ns, NNUC
        NCOUNT = NCOUNT+1
        FRC(I,J) = BFRC( NCOUNT )
      ENDDO
      ENDDO

      RETURN
      END


      SUBROUTINE NLFRC2( RWF,GNALM,MOR,BF,NR ) ! non-local part
                                               ! For link atoms

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
      include 'nlocd.i'

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NORA =  49 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NSL  = 3000 )

      PARAMETER ( PI   = 3.14159265358979323D0 )
C     PARAMETER ( RCUT = 2.5D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / MPIFRC / FRCM( NNUC,NDIM )

      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / FRCE  / FRC( NNUC,NDIM )
      COMMON / GRID2 / DX, DY, DZ
      COMMON / NCLR  / PNUC( NNUC,NDIM )
      COMMON / DBLG  / WNLOCS( NNUC,NSL ),WNLOCPX( NNUC,NSL ),
     *                 WNLOCPY(NNUC,NSL ),WNLOCPZ( NNUC,NSL )

      DIMENSION RWF(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION DWFX( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFY( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFZ( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION TWF(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      DIMENSION GNALM( NNUC,2,3,NORA )             ! num of elctrns,atoms,angular momentum quantum numbers
      DIMENSION FRNL(   2,3,3 ) 

      DIMENSION AFRC( NLINK1*NDIM ) 
      DIMENSION BFRC( NLINK1*NDIM ) 

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      DV  = DX*DY*DZ
      COM = 1.D0 / DSQRT( 4.D0*PI )
      CF  = 2.D0*BF
      ns = NNUC - NLINK + 1

C
C----- contribution from non-local part -----
C

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1
      
      DO 90 N=1,MOR                                   ! loop over alpha orbitals
      
      DO 100 K=1,JZ
      DO 100 J=1,JY
      DO 100 I=1,JX
      TWF( I,J,K ) = RWF( I,J,K,N )
100   CONTINUE

      CALL DWF2( TWF,DWFX,DWFY,DWFZ )                 ! derivative of alpha wavefunction
            
      DO 2000 NA=ns,NNUC

      NCOUNT = 0                                      ! reset counter
      DO K=1,3
      DO J=1,3
      DO I=1,2
      FRNL( I,J,K ) = 0.D0                            ! initialize force
      ENDDO
      ENDDO
      ENDDO
      
      DO 110 K=1,JZ
      DO 120 J=1,JY
      DO 130 I=1,JX
            
        RX = DX*( I-NNX ) - PNUC(NA,1) 
        RY = DY*( J-NNY ) - PNUC(NA,2)
        RZ = DZ*( K-NNZ ) - PNUC(NA,3)

        R2 = RX**2 + RY**2 + RZ**2
        R = DSQRT(R2)

        IF( R .LT. RCUT ) THEN     
        NCOUNT = NCOUNT + 1

C----- for s-component -----

        FRNL( 1,1,1 ) = FRNL( 1,1,1 ) +
     *    COM*WNLOCS(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNL( 1,1,2 ) = FRNL( 1,1,2 ) +
     *    COM*WNLOCS(NA,NCOUNT)*DWFY( I,J,K )*DV
        FRNL( 1,1,3 ) = FRNL( 1,1,3 ) +
     *    COM*WNLOCS(NA,NCOUNT)*DWFZ( I,J,K )*DV

C----- for p-component -----

        FRNL( 2,1,1 ) = FRNL( 2,1,1 ) +                 ! dx 
     *    3.D0*COM*WNLOCPX(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNL( 2,2,1 ) = FRNL( 2,2,1 ) +                 
     *    3.D0*COM*WNLOCPY(NA,NCOUNT)*DWFX( I,J,K )*DV 
        FRNL( 2,3,1 ) = FRNL( 2,3,1 ) +                 
     *    3.D0*COM*WNLOCPZ(NA,NCOUNT)*DWFX( I,J,K )*DV 

        FRNL( 2,1,2 ) = FRNL( 2,1,2 ) +                 ! dy 
     *    3.D0*COM*WNLOCPX(NA,NCOUNT)*DWFY( I,J,K )*DV 
        FRNL( 2,2,2 ) = FRNL( 2,2,2 ) +                 
     *    3.D0*COM*WNLOCPY(NA,NCOUNT)*DWFY( I,J,K )*DV 
        FRNL( 2,3,2 ) = FRNL( 2,3,2 ) +                 
     *    3.D0*COM*WNLOCPZ(NA,NCOUNT)*DWFY( I,J,K )*DV 

        FRNL( 2,1,3 ) = FRNL( 2,1,3 ) +                 ! dz 
     *    3.D0*COM*WNLOCPX(NA,NCOUNT)*DWFZ( I,J,K )*DV 
        FRNL( 2,2,3 ) = FRNL( 2,2,3 ) +                 
     *    3.D0*COM*WNLOCPY(NA,NCOUNT)*DWFZ( I,J,K )*DV 
        FRNL( 2,3,3 ) = FRNL( 2,3,3 ) +                 
     *    3.D0*COM*WNLOCPZ(NA,NCOUNT)*DWFZ( I,J,K )*DV 

        ENDIF

130   CONTINUE
120   CONTINUE
110   CONTINUE

      FRCM(NA,1) = FRCM(NA,1) - CF*( FRNL(1,1,1)*GNALM(NA,1,1,N)
     *          + FRNL(2,1,1)*GNALM(NA,2,1,N) 
     *          + FRNL(2,2,1)*GNALM(NA,2,2,N)
     *          + FRNL(2,3,1)*GNALM(NA,2,3,N) )
      FRCM(NA,2) = FRCM(NA,2) - CF*( FRNL(1,1,2)*GNALM(NA,1,1,N)
     *          + FRNL(2,1,2)*GNALM(NA,2,1,N) 
     *          + FRNL(2,2,2)*GNALM(NA,2,2,N)
     *          + FRNL(2,3,2)*GNALM(NA,2,3,N) )  
      FRCM(NA,3) = FRCM(NA,3) - CF*( FRNL(1,1,3)*GNALM(NA,1,1,N)
     *          + FRNL(2,1,3)*GNALM(NA,2,1,N) 
     *          + FRNL(2,2,3)*GNALM(NA,2,2,N)
     *          + FRNL(2,3,3)*GNALM(NA,2,3,N) )  

2000  CONTINUE

90    CONTINUE

      NCOUNT = 0
      DO J = 1, NDIM
      DO I = ns, NNUC
        NCOUNT = NCOUNT+1
        AFRC( NCOUNT ) = FRCM(I,J)
      ENDDO
      ENDDO

      NUMDAT = NLINK*NDIM
      CALL MPI_REDUCE( AFRC(1),BFRC(1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

      NCOUNT = 0
      DO J = 1, NDIM
      DO I = ns, NNUC
        NCOUNT = NCOUNT+1
        IF( NR .EQ. 0 ) THEN
         FRC(I,J) = BFRC( NCOUNT )
        ELSEIF( NR .EQ. 1 ) THEN
         FRC(I,J) = FRC(I,J) + BFRC( NCOUNT )
        ENDIF 
      ENDDO
      ENDDO

      IF(MYID.EQ.0) THEN
      WRITE(*,*) ' NLFRC2: '
      DO I = ns, NNUC
        WRITE(*,'(3F12.6)') FRC(I,1:3)
      ENDDO
      ENDIF

      RETURN
      END

C---------------------------------------------------------
C    block data for BHS parameters
C---------------------------------------------------------

      BLOCK DATA BHSLOC
      IMPLICIT REAL*8 ( A-H,O-Z )      
      COMMON / BHS / C1(100),C2(100),AL1(100),AL2(100)

      DATA C1(1),C2(1),AL1(1),AL2(1)                ! H
     &    / 1.1924D0,-0.1924D0,16.22D0,5.55D0 /
      DATA C1(3),C2(3),AL1(3),AL2(3)                ! Li 2008.03.21
     &    / 2.9081D0,-1.9081D0,1.84D0,0.73D0 /
      DATA C1(6),C2(6),AL1(6),AL2(6)                ! C
     &    / 1.5222D0,-0.5222D0,9.28D0,3.69D0 /
      DATA C1(7),C2(7),AL1(7),AL2(7)                ! N  2003.01.21
     &    / 1.4504D0,-0.4504D0,12.87D0,5.12D0 /
      DATA C1(8),C2(8),AL1(8),AL2(8)                ! O
     &    / 1.4224D0,-0.4224D0,18.09D0,7.19D0 /
      DATA C1(9),C2(9),AL1(9),AL2(9)                ! F  2008.03.21
     &    / 1.3974D0,-0.3974D0,23.78D0,9.45D0 /
      DATA C1(14),C2(14),AL1(14),AL2(14)            ! Si 2008.03.21
     &    / 1.6054D0,-0.6054D0,2.16D0,0.86D0 /
      DATA C1(15),C2(15),AL1(15),AL2(15)            ! P  2008.03.21
     &    / 1.4995D0,-0.4995D0,2.59D0,1.03D0 /
      DATA C1(17),C2(17),AL1(17),AL2(17)            ! Cl
     &    / 1.3860D0,-0.3860D0,3.48D0,1.38D0 /
C     DATA C1(20),C2(20),AL1(20),AL2(20)            ! Ca  BHS local
C    &    / 4.8360D0,-3.8360D0,1.61D0,0.45D0 /
      DATA C1(20),C2(20),AL1(20),AL2(20)            ! Ca  determined by Dr. K.Kobayashi 
     &    / 0.5000D0, 0.5000D0,1.00D0,0.50D0 /
C     DATA C1(25),C2(25),AL1(25),AL2(25)            ! Mn KBA  determined by Dr. K.Kobayashi 
C    &    / 0.5000D0, 0.5000D0,1.75D0,0.75D0 /
      DATA C1(25),C2(25),AL1(25),AL2(25)            ! Mn KBB  determined by Dr. K.Kobayashi 
     &    / 0.5000D0, 0.5000D0,1.00D0,0.50D0 /


      END

CCC MPI kokomade OK !! CCC
C-------------------------------------------
C     SUBROUTINE POINT CHARGE 
C-------------------------------------------

      SUBROUTINE POCH( SCRD )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"

      include 'QMpara.i'                        ! Vmol
      include "mpi.i"
!     include 'sizes.i'                   ! Vmol

C     PARAMETER ( NMAX =  80 )
C     PARAMETER ( NDIM =   3 )
C     PARAMETER ( NNUC =  36 )
      PARAMETER ( NSP  = 255 )
      PARAMETER ( ALPHA = 1.0D0    )

      PARAMETER ( EOO   = 6.614D0  )
      PARAMETER ( ECLO  = 7.559D0  )

      PARAMETER ( PI   = 3.14159265358979323D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / QCHA / CHARGE(10)
C     COMMON / SLT / VPCE(NMAX,NMAX,NMAX),VPCZ,PLJ
      COMMON / SLT / VPCE(NMAXX/NX,NMAXY/NY,NMAXZ/NZ),VPCZ,PLJ

C     DIMENSION SCRD( NDIM,4*NSP )
      DIMENSION SCRD( NDIM,maxatm )       ! Vmol

      CHARACTER ATOM*2

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      CHARGE(2) = 0.520D0
      CHARGE(3) = 0.520D0
      CHARGE(4) =-1.040D0

      DO K=1,JZ
      DO J=1,JY
      DO I=1,JX

        VPCE( I,J,K )=0.0D0
            
      ENDDO
      ENDDO
      ENDDO

C     NNMAX = NMAX/2 + 1
      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO 70 N=1,NSP

       DO K=1,JZ
        RZ=DZ*( K-NNZ ) - SCRD( 3,4*N-2 )             !H1 density-site energy
       DO J=1,JY
        RY=DY*( J-NNY ) - SCRD( 2,4*N-2 )
       DO I=1,JX
        RX=DX*( I-NNX ) - SCRD( 1,4*N-2 )         
        R =DSQRT(RX**2+RY**2+RZ**2)
        VPCE( I,J,K )=VPCE( I,J,K ) - CHARGE( 2 )*(ERF(ALPHA*R)/R)
       ENDDO
       ENDDO
       ENDDO

       DO K=1,JZ
        RZ=DZ*( K-NNZ ) - SCRD( 3,4*N-1 )             !H2 density-site energy
       DO J=1,JY
        RY=DY*( J-NNY ) - SCRD( 2,4*N-1 )
       DO I=1,JX
        RX=DX*( I-NNX ) - SCRD( 1,4*N-1 )         
        R =DSQRT(RX**2+RY**2+RZ**2)
        VPCE( I,J,K )=VPCE( I,J,K ) - CHARGE( 3 )*(ERF(ALPHA*R)/R)
       ENDDO
       ENDDO
       ENDDO

       DO K=1,JZ
        RZ=DZ*( K-NNZ ) - SCRD( 3,4*N )               !O1(Msite) density-site energy
       DO J=1,JY
        RY=DY*( J-NNY ) - SCRD( 2,4*N )
       DO I=1,JX
        RX=DX*( I-NNX ) - SCRD( 1,4*N )         
        R =DSQRT(RX**2+RY**2+RZ**2)
        VPCE( I,J,K )=VPCE( I,J,K ) - CHARGE( 4 )*(ERF(ALPHA*R)/R)
       ENDDO
       ENDDO
       ENDDO

70    CONTINUE

      VPCZ=0.0
      PLJ=0.0

      IF(MYID.EQ.0) THEN

      DO 80 N=1,NSP
        DO 90 M=1,NNUC 

        RX=PNUC( M,1 ) - SCRD( 1,4*N-2 )      !H1 nuclear-site energy
        RY=PNUC( M,2 ) - SCRD( 2,4*N-2 )
        RZ=PNUC( M,3 ) - SCRD( 3,4*N-2 )
        R =DSQRT(RX**2+RY**2+RZ**2)
        VPCZ=VPCZ + ZVAL( M )*CHARGE( 2 )/R
         
        RX=PNUC( M,1 ) - SCRD( 1,4*N-1 )      !H2 nuclear-site energy
        RY=PNUC( M,2 ) - SCRD( 2,4*N-1 )
        RZ=PNUC( M,3 ) - SCRD( 3,4*N-1 )
        R =DSQRT(RX**2+RY**2+RZ**2)
        VPCZ=VPCZ + ZVAL( M )*CHARGE( 3 )/R
         
        RX=PNUC( M,1 ) - SCRD( 1,4*N )        !O1(Msite) nuclear-site energy
        RY=PNUC( M,2 ) - SCRD( 2,4*N )
        RZ=PNUC( M,3 ) - SCRD( 3,4*N )
        R =DSQRT(RX**2+RY**2+RZ**2)
        VPCZ=VPCZ + ZVAL( M )*CHARGE( 4 )/R

90      CONTINUE
80    CONTINUE

      DO 100 M=1,NNUC
      DO 110 N=1,NSP

        RX=PNUC( M,1 ) - SCRD( 1,4*N-3 )      !nuclear-site LJ potential
        RY=PNUC( M,2 ) - SCRD( 2,4*N-3 )
        RZ=PNUC( M,3 ) - SCRD( 3,4*N-3 )
        R =DSQRT(RX**2+RY**2+RZ**2)
        PLJ=PLJ+4.*EPSQM(M)*((SIG(M)/R)**12.-(SIG(M)/R)**6.) ! Vmol
       
        IF(R.LT.2) WRITE(*,*) 'N7=',N,' R=',R 

110   CONTINUE
100   CONTINUE

      WRITE(*,*)
      WRITE(*,*)'Nuclear - Site coulomb energy'
      WRITE(*,*)'  VPCZ = ' ,VPCZ
      WRITE(*,*)
      WRITE(*,*)'Nuclear - Site LJ potential energy'
      WRITE(*,*)'  PLJ  = ' ,PLJ

      ENDIF

      CALL MPI_BCAST( VPCZ,1,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )
      CALL MPI_BCAST( PLJ,1,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )

      RETURN
      END

C-------------------------------------------
C     SUBROUTINE QM/MM FORCE
C-------------------------------------------

      SUBROUTINE QMF( SCRD,FRCS,FRCN,RHO )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'                   ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC =  36 )
      PARAMETER ( NSP = 255 )
      PARAMETER ( DCUT  = 1.D-7 ) 
      PARAMETER ( ALPHA = 1.0D0 )

      PARAMETER ( PI    = 3.14159265358979323D0 )
      
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / QCHA  /   CHARGE(10)

C     DIMENSION SCRD( NDIM,4*NSP )
      DIMENSION SCRD( NDIM,maxatm )       ! Vmol

      DIMENSION FDSS(  4*NSP,NDIM )               ! density-site force on site
      DIMENSION FODSS( 4*NSP,NDIM )               ! density-site force on site
      DIMENSION FNSS(  4*NSP,NDIM )               ! nuclear-site force on site
      DIMENSION FQMLJS(4*NSP,NDIM )               ! LJ force on site
C     DIMENSION FRCS(  4*NSP,NDIM )               ! QM/MM force on site
      DIMENSION FRCS( maxatm,NDIM )       ! Vmol

      DIMENSION FNSN(   NNUC,NDIM )               ! nuclear-site force on nuclear
      DIMENSION FQMLJN( NNUC,NDIM )               ! LJ force on nuclear
      DIMENSION FRCN(   NNUC,NDIM )               ! QM/MM force on nuclear

      DIMENSION RHO( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
C     CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      DV=DX*DY*DZ

      DO L=1,3
      DO N=1,4*NSP

          FODSS( N,L )  = 0.D0
          FNSS(  N,L )  = 0.D0
          FQMLJS(N,L )  = 0.D0
          FRCS(  N,L )  = 0.D0
            
      ENDDO
      ENDDO
     
      DO L=1,3
      DO N=1,NNUC
     
          FNSN( N,L )   = 0.D0
          FQMLJN( N,L ) = 0.D0
          FRCN( N,L )   = 0.D0
            
      ENDDO
      ENDDO

C-----density-site force CALCULATION
     
C     NNMAX   = NMAX/2 + 1
      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      C2DV   = CHARGE(2) * DV
      C3DV   = CHARGE(3) * DV
      C4DV   = CHARGE(4) * DV
      SQPII  = 1.D0/DSQRT(PI)
      A2SQPI = 2.D0*ALPHA*SQPII

      DO 60 K=1,JZ
        ADZ   = DZ*( K-NNZ )
      DO 70 J=1,JY
        ADY   = DY*( J-NNY )
      DO 80 I=1,JX
        ADX   = DX*( I-NNX )

       IF( RHO(I,J,K).GT.DCUT ) THEN

       DO 50 N=1,NSP

       RX    = ADX - SCRD( 1,4*N-2 )    !H1 density-site force 
       RY    = ADY - SCRD( 2,4*N-2 )
       RZ    = ADZ - SCRD( 3,4*N-2 )
       R2    = RX**2+RY**2+RZ**2
       R     = DSQRT( R2 )
       R3    = R*R2
       APR   = ALPHA*R
       TFDSS = C2DV*(A2SQPI*DEXP(-APR**2)/R2-ERF(APR)/R3)*RHO(I,J,K)

       FODSS(4*N-2,1) = FODSS(4*N-2,1)-RX*TFDSS
       FODSS(4*N-2,2) = FODSS(4*N-2,2)-RY*TFDSS
       FODSS(4*N-2,3) = FODSS(4*N-2,3)-RZ*TFDSS

       RX    = ADX - SCRD( 1,4*N-1 )    !H2 density-site force 
       RY    = ADY - SCRD( 2,4*N-1 )
       RZ    = ADZ - SCRD( 3,4*N-1 )
       R2    = RX**2+RY**2+RZ**2
       R     = DSQRT( R2 )
       R3    = R*R2
       APR   = ALPHA*R
       TFDSS = C3DV*(A2SQPI*DEXP(-APR**2)/R2-ERF(APR)/R3)*RHO(I,J,K)
        
       FODSS(4*N-1,1) = FODSS(4*N-1,1)-RX*TFDSS
       FODSS(4*N-1,2) = FODSS(4*N-1,2)-RY*TFDSS
       FODSS(4*N-1,3) = FODSS(4*N-1,3)-RZ*TFDSS

       RX    = ADX - SCRD( 1,4*N )      !O1 density-site force 
       RY    = ADY - SCRD( 2,4*N )
       RZ    = ADZ - SCRD( 3,4*N )
       R2    = RX**2+RY**2+RZ**2
       R     = DSQRT( R2 )
       R3    = R*R2
       APR   = ALPHA*R
       TFDSS = C4DV*(A2SQPI*DEXP(-APR**2)/R2-ERF(APR)/R3)*RHO(I,J,K)

       FODSS(4*N,1) = FODSS(4*N,1)-RX*TFDSS
       FODSS(4*N,2) = FODSS(4*N,2)-RY*TFDSS
       FODSS(4*N,3) = FODSS(4*N,3)-RZ*TFDSS

50     CONTINUE

       ENDIF

80    CONTINUE
70    CONTINUE
60    CONTINUE

      NUMDAT = 4*NSP*3
      CALL MPI_REDUCE( FODSS(1,1),FDSS(1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

C-----nuclear-site force CALCULATION

      IF(MYID.EQ.0) THEN

      DO 90 N=1,NSP
      DO 100 M=1,NNUC 

        RX=PNUC( M,1 ) - SCRD( 1,4*N-2 )      !H1 nuclear-site force
        RY=PNUC( M,2 ) - SCRD( 2,4*N-2 )
        RZ=PNUC( M,3 ) - SCRD( 3,4*N-2 )
        R2=RX**2+RY**2+RZ**2
        R =DSQRT(R2)
        R3=R*R2
        CHRX = CHARGE(2)*RX/R3
        CHRY = CHARGE(2)*RY/R3
        CHRZ = CHARGE(2)*RZ/R3

        FNSS(4*N-2,1)=FNSS(4*N-2,1)+ZVAL( M )*CHRX
        FNSS(4*N-2,2)=FNSS(4*N-2,2)+ZVAL( M )*CHRY
        FNSS(4*N-2,3)=FNSS(4*N-2,3)+ZVAL( M )*CHRZ

        FNSN(M,1) = FNSN(M,1)+ZVAL( M )*CHRX
        FNSN(M,2) = FNSN(M,2)+ZVAL( M )*CHRY
        FNSN(M,3) = FNSN(M,3)+ZVAL( M )*CHRZ

        RX=PNUC( M,1 ) - SCRD( 1,4*N-1 )      !H2 nuclear-site force
        RY=PNUC( M,2 ) - SCRD( 2,4*N-1 )
        RZ=PNUC( M,3 ) - SCRD( 3,4*N-1 )
        R2=RX**2+RY**2+RZ**2
        R =DSQRT(R2)
        R3=R*R2
        CHRX = CHARGE(3)*RX/R3
        CHRY = CHARGE(3)*RY/R3
        CHRZ = CHARGE(3)*RZ/R3

        FNSS(4*N-1,1)=FNSS(4*N-1,1)+ZVAL( M )*CHRX
        FNSS(4*N-1,2)=FNSS(4*N-1,2)+ZVAL( M )*CHRY
        FNSS(4*N-1,3)=FNSS(4*N-1,3)+ZVAL( M )*CHRZ

        FNSN(M,1) = FNSN(M,1)+ZVAL( M )*CHRX
        FNSN(M,2) = FNSN(M,2)+ZVAL( M )*CHRY
        FNSN(M,3) = FNSN(M,3)+ZVAL( M )*CHRZ

        RX=PNUC( M,1 ) - SCRD( 1,4*N )        !O1 nuclear-site force
        RY=PNUC( M,2 ) - SCRD( 2,4*N )
        RZ=PNUC( M,3 ) - SCRD( 3,4*N )
        R2=RX**2+RY**2+RZ**2
        R =DSQRT(R2)
        R3=R*R2
        CHRX = CHARGE(4)*RX/R3
        CHRY = CHARGE(4)*RY/R3
        CHRZ = CHARGE(4)*RZ/R3

        FNSS(4*N,1)=FNSS(4*N,1)+ZVAL( M )*CHRX
        FNSS(4*N,2)=FNSS(4*N,2)+ZVAL( M )*CHRY
        FNSS(4*N,3)=FNSS(4*N,3)+ZVAL( M )*CHRZ

        FNSN(M,1) = FNSN(M,1)+ZVAL( M )*CHRX
        FNSN(M,2) = FNSN(M,2)+ZVAL( M )*CHRY
        FNSN(M,3) = FNSN(M,3)+ZVAL( M )*CHRZ
         
100   CONTINUE
90    CONTINUE
      
C---- LJ force CALCULATION
      
      TWO   = 2.D0
      SIX   = 6.D0
 
      TF    = 24.D0

      DO 160 M=1,NNUC

C         TFEPS = TF*EPS(M)
          TFEPS = TF*EPSQM(M) ! Vmol

          DO N=1,NSP

            RX = PNUC( M,1 ) - SCRD( 1,4*N-3 )      !LJ force
            RY = PNUC( M,2 ) - SCRD( 2,4*N-3 )
            RZ = PNUC( M,3 ) - SCRD( 3,4*N-3 )
            R2 = RX**2+RY**2+RZ**2
            R  = DSQRT(R2)
            ASIG   = SIG(M)/R
            ASIG6  = ASIG**6
            ASIG12 = ASIG6*ASIG6
            BSIG   = 2.D0*ASIG12 - ASIG6
            CREPS  = TFEPS/R2*BSIG
            CREPSX = CREPS*RX
            CREPSY = CREPS*RY
            CREPSZ = CREPS*RZ

            FQMLJS(4*N-3,1)=FQMLJS(4*N-3,1)+CREPSX
            FQMLJS(4*N-3,2)=FQMLJS(4*N-3,2)+CREPSY
            FQMLJS(4*N-3,3)=FQMLJS(4*N-3,3)+CREPSZ

            FQMLJN(M,1)=FQMLJN(M,1)+CREPSX
            FQMLJN(M,2)=FQMLJN(M,2)+CREPSY
            FQMLJN(M,3)=FQMLJN(M,3)+CREPSZ
        
          ENDDO

160   CONTINUE

      DO L=1,3
      DO N=1,4*NSP
        FRCS(N,L) = FDSS(N,L) - FNSS(N,L) - FQMLJS(N,L)  ! QM/MM force on site
      ENDDO
      ENDDO

      DO L=1,3
      DO N=1,NNUC
        FRCN(N,L) = FNSN(N,L) + FQMLJN(N,L)              ! QM/MM force on nuclear
      ENDDO
      ENDDO

      REWIND(61)
      DO N=1,NNUC
      DO L=1,1
       WRITE(61,*) 'FNSN  = ',FNSN(N,L)
       WRITE(61,*) 'FQMLJN= ',FQMLJN(N,L) 
      ENDDO
      ENDDO

      REWIND(60)
      DO N=1,4*NSP
      DO L=1,1
       WRITE(60,*) 'FDSS=',  FDSS(N,L)  
       WRITE(60,*) 'FNSS=',  FNSS(N,L)  
       WRITE(60,*) 'FQMLJS=',FQMLJS(N,L)  
      ENDDO
      ENDDO

      ENDIF

      RETURN
      END

C---------------------------------------------------------

      SUBROUTINE TRANS( SCRD,FRCS,NCU )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include 'QMpara.i'                  ! Vmol
!     include 'sizes.i'                   ! Vmol

      PARAMETER ( NSP = 255 )
      PARAMETER ( AUL = 0.5291771D-10 )     
      PARAMETER ( AUM = 9.109534D-31 )
      PARAMETER ( AUT = 2.418884D-17 )

      COMMON/CMFACT/FAC(4),FLJ(4),FKT(3),CLR(6)

C     DIMENSION SCRD( 3,4*NSP )
      DIMENSION SCRD( 3,maxatm )          ! Vmol
C     DIMENSION FRCS( 4*NSP,3 )
      DIMENSION FRCS( maxatm,3 )          ! Vmol

      IF ( NCU.EQ.1 ) THEN         !translation to a.u.

      DO 10 N=1,4*NSP
        DO 20 M=1,3

          SCRD( M,N )=FAC(1)*SCRD( M,N )/AUL

C         FRCS( N,M )=FAC(2)*FRCS( N,M )/FAC(1)*AUT*AUT/AUM/AUL

          FRCS( N,M )=FAC(2)*FRCS( N,M )/FAC(1)/8.23885992D-8


20      CONTINUE
10    CONTINUE

      ELSEIF( NCU.EQ.2 ) THEN      !translation from a.u.

      DO 50 N=1,4*NSP
        DO 60 M=1,3

          SCRD( M,N )=AUL*SCRD( M,N )/FAC(1)

C         FRCS( N,M )=FAC(1)*FRCS( N,M )/FAC(2)*AUM*AUL/AUT/AUT

          FRCS( N,M )=FAC(1)*FRCS( N,M )/FAC(2)*8.23885992D-8

60      CONTINUE
50    CONTINUE

      ENDIF

      RETURN
      END

C-------------------------------------------
C
C     SUBROUTINE SOLUTE-SOLVENT INTERACTION
C
C-------------------------------------------

      SUBROUTINE SSINT( RWFA,RWFB,IMD,SSE )

      IMPLICIT REAL*8 ( A-H,O-Z )
      IMPLICIT INTEGER*4 ( I-N )

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NORA =  49 )
C     PARAMETER ( NORB =  1 )
C     PARAMETER ( MORA =  49 )
C     PARAMETER ( MORB =  1 )
C     PARAMETER ( NSP = 255 )
 
      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / GRID2 / DX, DY, DZ
C     COMMON / SLT / VPCE(NMAX,NMAX,NMAX),VPCZ,PLJ
      COMMON / SLT / VPCE(NMAXX/NX,NMAXY/NY,NMAXZ/NZ),VPCZ,PLJ
      COMMON / REST / PDPM(3),TSSE

      DIMENSION RWFA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION RWFB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
C     DIMENSION RWFB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,1 )

C     DIMENSION QESPM( NSP )
C     DIMENSION SSEM(  NSP )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      SSE  = 0.D0
      QQESP = 0.D0

      DV=DX*DY*DZ

      DO L=1,MORA
      DO K=1,JZ
      DO J=1,JY
      DO I=1,JX
      
        QQESP = QQESP + RWFA( I,J,K,L )
     *       * VPCE( I,J,K )*RWFA( I,J,K,L )*DV

      ENDDO
      ENDDO
      ENDDO
      ENDDO

      IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
     *                 .OR.EXC.EQ.'UXalpha') THEN     

        DO L=1,MORB
        DO K=1,JZ
        DO J=1,JY
        DO I=1,JX
      
          QQESP = QQESP + RWFB( I,J,K,L )
     *         * VPCE( I,J,K )*RWFB( I,J,K,L )*DV

        ENDDO
        ENDDO
        ENDDO
        ENDDO

      ELSEIF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
     *                     .OR.EXC.EQ.'RXalpha') THEN     

        QQESP = 2.D0*QQESP

      ENDIF

      CALL MPI_REDUCE( QQESP,QESP,1,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

      IF(MYID.EQ.0) THEN

      SSE  = QESP + VPCZ + PLJ      

      write(*,*)
      write(*,100) 'QESP,VPCZ,PLJ,SSE = ',QESP,VPCZ,PLJ,SSE

      WRITE(63,*) SSE

      IF( IMD .EQ. 5000 ) THEN
       TSSE=0.0D0
      ENDIF

      IF( IMD .GT. 5000 ) THEN

       TSSE = TSSE + SSE
       ASSE = TSSE/(IMD-5000) 
       WRITE(*,*) '  AVERAGE OF SOLUTE-SOLVENT INT. = ',ASSE

      ENDIF

      ENDIF

100   FORMAT(A,4F15.8)

      RETURN
      END

C-------------------------------------------
C
C     SUBROUTINE SOLUTE-SOLVENT INTERACTION
C
C-------------------------------------------

      SUBROUTINE SSINT1( RHO,IMD )

      IMPLICIT REAL*8 ( A-H,O-Z )
      IMPLICIT INTEGER*4 ( I-N )

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NORA = 49 )
C     PARAMETER ( NORB =  1 )
C     PARAMETER ( MORA = 49 )
C     PARAMETER ( MORB =  1 )
C     PARAMETER ( NSP = 255 )
 
      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / GRID2 / DX, DY, DZ
C     COMMON / SLT / VPCE(NMAX,NMAX,NMAX),VPCZ,PLJ
      COMMON / SLT / VPCE(NMAXX/NX,NMAXY/NY,NMAXZ/NZ),VPCZ,PLJ
      COMMON / REST / PDPM(3),TSSE

      DIMENSION RHO( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

C     DIMENSION QESPM( NSP )
C     DIMENSION SSEM(  NSP )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      SSE   = 0.D0
      QQESP = 0.D0

      DV=DX*DY*DZ

      DO K=1,JZ
      DO J=1,JY
      DO I=1,JX
      
        QQESP = QQESP + 
     *          VPCE( I,J,K )*RHO( I,J,K )*DV

      ENDDO
      ENDDO
      ENDDO

C     IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
C    *                 .OR.EXC.EQ.'UXalpha') THEN     

C       DO L=1,MORB
C       DO K=1,JZ
C       DO J=1,JY
C       DO I=1,JX
C     
C         QQESP = QQESP + RWFB( I,J,K,L )
C    *         * VPCE( I,J,K )*RWFB( I,J,K,L )*DV

C       ENDDO
C       ENDDO
C       ENDDO
C       ENDDO

C     ELSEIF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
C    *                     .OR.EXC.EQ.'RXalpha') THEN     

C       QQESP = 2.D0*QQESP

C     ENDIF

      CALL MPI_REDUCE( QQESP,QESP,1,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

      IF(MYID.EQ.0) THEN

      SSE  = QESP + VPCZ + PLJ      

      WRITE(*,*)
      WRITE(*,100) 'QESP,VPCZ,PLJ,SSE =',QESP,VPCZ,PLJ,SSE
      WRITE(*,*)

      WRITE(63,*) SSE

      IF( IMD .EQ. 5000 ) THEN
       TSSE=0.0D0
      ENDIF

      IF( IMD .GT. 5000 ) THEN

       TSSE = TSSE + SSE
       ASSE = TSSE/(IMD-5000) 
       WRITE(*,*) '  AVERAGE OF SOLUTE-SOLVENT INT. = ',ASSE

      ENDIF

      ENDIF

100   FORMAT(X,A19,4F11.6)

      RETURN
      END

C--------------------------------------------------      
C
      SUBROUTINE CS (SCRD)
C
C     warning : available only for water solvent
C
C--------------------------------------------------      

      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)

      include 'QMpara.i'                  ! Vmol
!     include 'sizes.i'                   ! Vmol

      PARAMETER ( NSP = 255 )
C     PARAMETER ( NNUC = 36 )
C     PARAMETER ( NDIM =  3 )

      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV

C     DIMENSION SCRD(NDIM,4*NSP)
      DIMENSION SCRD(NDIM,maxatm)         ! Vmol
      DIMENSION SSCRD(NDIM,4*NSP)
      
C---- output CS.DAT ----

      REWIND(99)
      WRITE(99,*) NSP*3
   
      DO I=1,NSP
      DO J=1,NDIM
      SSCRD(J,4*I-3)=SCRD(J,4*I-3)*0.529177
      SSCRD(J,4*I-2)=SCRD(J,4*I-2)*0.529177
      SSCRD(J,4*I-1)=SCRD(J,4*I-1)*0.529177
      ENDDO
      ENDDO

      DO I=1,NSP
       WRITE(99,90) 3*I-2,SSCRD(1,4*I-3),
     *             SSCRD(2,4*I-3),SSCRD(3,4*I-3),3*I-1,3*I
       WRITE(99,91) 3*I-1,SSCRD(1,4*I-2),
     *             SSCRD(2,4*I-2),SSCRD(3,4*I-2),3*I-2
       WRITE(99,91) 3*I  ,SSCRD(1,4*I-1),
     *             SSCRD(2,4*I-1),SSCRD(3,4*I-1),3*I-2
      ENDDO


90    FORMAT('O ',1X,I3,1X,F8.4,1X,F8.4,1X,F8.4,1X,'82 ',1X,I3,1X,I3)
91    FORMAT('H ',1X,I3,1X,F8.4,1X,F8.4,1X,F8.4,1X,'11 ',1X,I3)

      RETURN
      END


C-------------------------------------------------------
C     Write coordinates,velocities,forces
C
      SUBROUTINE  WCRD(CRD,SCRD,VLC,FRC1,IMD)
C-------------------------------------------------------

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
        
      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'                   ! Vmol
!     include 'atoms.i'                   ! Vmol

      PARAMETER ( NSP = 255 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NLINK1 = 10 )           ! Vmol
     
      DIMENSION CRD(NNUC,3)
C     DIMENSION SCRD(3,4*NSP)
      DIMENSION SCRD(3,maxatm)            ! Vmol
      DIMENSION VLC( NNUC,3 )
      DIMENSION FRC1(NNUC,3)

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)

      IF(MOD(IMD,10000).EQ.1) THEN ! for test MM calculation
      WRITE(32,*) IMD
      WRITE(34,*) IMD
      DO NA=1,NNUC                                             ! write coordinates
      WRITE(32, 99) ZA(NA),ZVAL(NA),
     *              CRD(NA,1),CRD(NA,2),CRD(NA,3)
      WRITE(34,999) ZA(NA),FRC1(NA,1),FRC1(NA,2),FRC1(NA,3)    ! write forces
      ENDDO
      IF(NOPT.EQ.'MD') THEN
      WRITE(33,*) IMD
      DO NA=1,NNUC
      WRITE(33,999) ZA(NA),VLC(NA,1),VLC(NA,2),VLC(NA,3)       ! write velocities
      ENDDO
      ENDIF
      ENDIF ! for test MM calculation

      IF(NQMMM.EQ.'QMMM') THEN

      OXY1  = 8.0
      OXY2  = 6.0
      HYD   = 1.0
      SITEM = 0.0

      IF(IMD.GT.1000) THEN

      WRITE(70,*) IMD
      DO I = 1,NSP
      WRITE(70,99) OXY1,OXY2,
     *             SCRD(1,4*I-3),SCRD(2,4*I-3),SCRD(3,4*I-3)
      WRITE(70,99) HYD,HYD,
     *             SCRD(1,4*I-2),SCRD(2,4*I-2),SCRD(3,4*I-2)
      WRITE(70,99) HYD,HYD,
     *             SCRD(1,4*I-1),SCRD(2,4*I-1),SCRD(3,4*I-1)
      WRITE(70,99) SITEM,SITEM,
     *             SCRD(1,4*I),SCRD(2,4*I),SCRD(3,4*I)
      ENDDO

      ENDIF
      
      ELSEIF(NQMMM.EQ.'LINK') THEN               ! Vmol

      IF( MOD(IMD,10) .EQ. 0) THEN

      REWIND(70)
C                                                ! Vmol
C     number of loop in QM and MM subsystems     ! Vmol
C                                                ! Vmol
C     nlmm = n-NMM2                              ! Vmol
      nlmm = n-NNUC+NLINK                        ! Vmol
                                                 ! Vmol
      AZERO = 0.d0                               ! Vmol
                                                 ! Vmol
      DO I = 1, nlmm                             ! Vmol

C     WRITE(70,99) AZERO,AZERO,                  ! Vmol
C    *             SCRD(1,I),SCRD(2,I),SCRD(3,I) ! Vmol

      SX = SCRD(1,I)*0.529177D0
      SY = SCRD(2,I)*0.529177D0
      SZ = SCRD(3,I)*0.529177D0

      IF( MOD(I,3) .EQ. 1) THEN                  ! specific treatment for water solvent 
        WRITE(70,999) 'O',SX,SY,SZ
      ELSE
        WRITE(70,999) 'H',SX,SY,SZ
      ENDIF

      ENDDO                                      ! Vmol
                                                 ! Vmol
      ENDIF
      ENDIF

99    FORMAT( F6.1,F6.1,3F12.6 )
999   FORMAT( X,A1,7X,3F12.6 )

      RETURN
      END

C--------------------------------------------------------
C     SUBROUTINE FOR RESTART CALCULATION
C--------------------------------------------------------

      SUBROUTINE AGAIN(NK,KMD)

      IMPLICIT REAL*8 ( A-H,O-Z )
      IMPLICIT INTEGER*4 ( I-N )

      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NNUC =  36 )
      PARAMETER ( NSP = 255 )
      PARAMETER ( NDX = 100 )
C     PARAMETER ( NMAX = 80 )

      COMMON/CMCRDS/XYZ(3,NSP),QTN(4,NSP),ROT(9,NSP),BOXL,HBOX,WP(3,NSP)
      COMMON/CMVELO/VMOM(3,NSP),ANGM(3,NSP)
      COMMON / REST / PDPM(3),TSSE
      COMMON / NRDF1 / NMRDF,NRDF(NNUC)
      COMMON / NNRDF / NGO(NNUC,NDX),NGH(NNUC,NDX)

      IF(NK.EQ.1) THEN

C------- READ RESTART FILE -------

        OPEN(17,FILE='chk.inp',STATUS='OLD')

        WRITE(*,*)
        WRITE(*,*)'--- READING RESTART FILE ---'

        DO 10 I = 1,NSP
          READ(17,*) XYZ(1,I),XYZ(2,I),XYZ(3,I)
          READ(17,*) QTN(1,I),QTN(2,I),QTN(3,I),QTN(4,I)
          READ(17,*) VMOM(1,I),VMOM(2,I),VMOM(3,I)
          READ(17,*) ANGM(1,I),ANGM(2,I),ANGM(3,I)
10      CONTINUE  

        DO 20 I = 1,3
          READ(17,*) PDPM(I)
20      CONTINUE  

        IF(NMRDF.NE.0) THEN
        DO I = 1,NMRDF
        DO J = 1,NDX
          READ(17,*) NGO(I,J),NGH(I,J)
        ENDDO
        ENDDO
        ENDIF

        READ(17,*) TSSE,KMD

        CLOSE(17)

      ENDIF

C------- WRITE RESTART FILE -------

      IF(NK.EQ.2) THEN
       
        REWIND(18)

        DO 30 I = 1,NSP
          WRITE(18,*) XYZ(1,I),XYZ(2,I),XYZ(3,I)
          WRITE(18,*) QTN(1,I),QTN(2,I),QTN(3,I),QTN(4,I)
          WRITE(18,*) VMOM(1,I),VMOM(2,I),VMOM(3,I)
          WRITE(18,*) ANGM(1,I),ANGM(2,I),ANGM(3,I)
30     CONTINUE

       DO 40 I = 1,3
         WRITE(18,*) PDPM(I)
40     CONTINUE

        IF(NMRDF.NE.0) THEN
        DO I = 1,NMRDF
        DO J = 1,NDX
          WRITE(18,*) NGO(I,J),NGH(I,J)
        ENDDO
        ENDDO
        ENDIF

       WRITE(18,*) TSSE,KMD

      ENDIF

      RETURN
      END

C ----------------------------------------------------------------
C     CALCULATION OF DIPOLE MOMENT
C ---------------------------------------------------------------

      SUBROUTINE DPL ( RHO,IMD )
     
      IMPLICIT REAL*8 ( A-H,O-Z )
      IMPLICIT INTEGER*4 ( I-N )

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NNUC = 36 )
C     PARAMETER ( NDIM =  3 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV

      COMMON / GRID2 / DX, DY, DZ            
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / REST / PDPM(3),TSSE

      DIMENSION RHO ( NMAXX/NX,NMAXY/NY,NMAXZ/NZ)
      DIMENSION DPM ( NDIM)
      DIMENSION DPM1 (NDIM)
      DIMENSION ADPM (NDIM)

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
C     CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      DV = DX*DY*DZ

      DO I = 1,NDIM
      DPM1(I) = 0.D0
      ENDDO

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      RRHO = 0.D0     

      DO K = 1,JZ
      DO J = 1,JY
      DO I = 1,JX
      DPM1(1) = DPM1(1) - RHO(I,J,K)*DX*(I-NNX)*DV 
      DPM1(2) = DPM1(2) - RHO(I,J,K)*DY*(J-NNY)*DV 
      DPM1(3) = DPM1(3) - RHO(I,J,K)*DZ*(K-NNZ)*DV 
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE( DPM1(1),DPM(1),3,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

      IF(MYID.EQ.0) THEN
  
        DO J = 1,3    
        DO I = 1,NNUC    
        DPM(J) = DPM(J) + ZVAL(I)*PNUC(I,J)
        ENDDO
        ENDDO
  
        DB = 2.54154D0
        WRITE(*,*) 
        WRITE(*,*) ' -------- Dipole Moment [Debye] ---------'
        WRITE(*,*) 
        WRITE(*,*) ' Dipole (X)   = ',  DPM(1)*DB
        WRITE(*,*) ' Dipole (Y)   = ',  DPM(2)*DB
        WRITE(*,*) ' Dipole (Z)   = ',  DPM(3)*DB
  
        TDPM = DSQRT( DPM(1)**2 + DPM(2)**2
     *                          + DPM(3)**2 ) 
  
        WRITE(*,*) ' Total Dipole = ',  TDPM*DB
  
        IF (IMD.EQ.1000) THEN
        DO I = 1,NDIM
        PDPM (I) = 0.D0
        ADPM (I) = 0.D0
        ENDDO
        ENDIF
  
        IF (IMD.GT.1000) THEN
  
        NIM = IMD-1000
        DO I = 1,NDIM
        PDPM(I) = PDPM(I) + DPM(I)
        ADPM(I) = PDPM(I)/NIM
        ENDDO
  
        TADPM = DSQRT( ADPM(1)**2+ADPM(2)**2
     *                           +ADPM(3)**2 )
  
        WRITE(*,*)                     
        WRITE(*,197) '  QM dipole(vector) = ',ADPM(1),ADPM(2),ADPM(3)
        WRITE(*,*)   '  QM total dipole   = ',TADPM
  
        ENDIF
  
      ENDIF

197   FORMAT( 1X,A23,3F12.6 )      

      RETURN
      END

C ----------------------------------------------------------------
C     RADIAL DITRIBUTION FUNCTION
C ---------------------------------------------------------------

      SUBROUTINE RDF(CRD,SCRD,IMD)

      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'                   ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NNUC =  36 )
      PARAMETER ( NSP = 255 )
      PARAMETER ( NDX = 100 )
      PARAMETER ( AUL  = 0.5291771D-10 )     
      PARAMETER ( PI   = 3.14159265358979323D0 )

      DIMENSION CRD (NNUC,3)
C     DIMENSION SCRD( 3,4*NSP )
      DIMENSION SCRD( 3,maxatm )          ! Vmol

      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / CMCRDS / XYZ(3,NSP),QTN(4,NSP),
     *                ROT(9,NSP),BOXL,HBOX,WP(3,NSP)
      COMMON / CMFACT / FAC(4),FLJ(4),FKT(3),CLR(6)
      COMMON / NRDF1 / NMRDF,NRDF(NNUC)
      COMMON / NNRDF / NGO(NNUC,NDX),NGH(NNUC,NDX)
     
      DX = 0.188972687D0
      TBOXL = FAC(1)*BOXL/AUL
      TVOL  = TBOXL**3
    
      IF (IMD.EQ.1) THEN
      DO I = 1 ,NNUC
      DO J = 1 ,NDX 
        NGO(I,J) = 0
        NGH(I,J) = 0
      ENDDO
      ENDDO
      ENDIF

      DO J = 1 ,NMRDF
 
      K=NRDF(J)
 
      DO I = 1 ,NSP

      RO = (CRD(K,1)-SCRD(1,4*I-3))**2
     *   + (CRD(K,2)-SCRD(2,4*I-3))**2
     *   + (CRD(K,3)-SCRD(3,4*I-3))**2
     
      SRO = DSQRT(RO) / DX
      NRO = INT(SRO)

      IF(NRO.LE.NDX) THEN

      NGO(J,NRO) = NGO(J,NRO) + 1

      ENDIF

      RH = (CRD(K,1)-SCRD(1,4*I-2))**2
     *   + (CRD(K,2)-SCRD(2,4*I-2))**2
     *   + (CRD(K,3)-SCRD(3,4*I-2))**2
     
      SRH = DSQRT(RH) / DX
      NRH = INT(SRH)

      IF(NRH.LE.NDX) THEN

      NGH(J,NRH) = NGH(J,NRH) + 1

      ENDIF

      RH = (CRD(K,1)-SCRD(1,4*I-1))**2
     *   + (CRD(K,2)-SCRD(2,4*I-1))**2
     *   + (CRD(K,3)-SCRD(3,4*I-1))**2
     
      SRH = DSQRT(RH) / DX
      NRH = INT(SRH)

      IF(NRH.LE.NDX) THEN

      NGH(J,NRH) = NGH(J,NRH) + 1

      ENDIF

      ENDDO  
      
      ENDDO  

      DDX  = DX*0.529177D0
      SDDX = 0.5D0*DDX

      IF (MOD(IMD,1000).EQ.0) THEN

      DO J = 1,NMRDF
      N = J+100
      WRITE(N,*) IMD 
      DO I = 1,NDX
      R  = DX*DBLE(I)
      R2 = R**2
      SDV = 4.D0*PI*R2*DX+4.D0*PI*R*DX*DX+4.D0/3.D0*PI*DX**3
      OTMP = NGO(J,I)*TVOL/(      NSP*IMD*SDV )    
      HTMP = NGH(J,I)*TVOL/( 2.D0*NSP*IMD*SDV )
      RDX = DDX * I + SDDX 
      WRITE(N,*) RDX,OTMP,HTMP
      ENDDO
      ENDDO

      ENDIF
      
      RETURN
      END


C----------------------------------------------------------------

      SUBROUTINE WSS ( SSE )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      PARAMETER ( NSP = 255 )
      
      COMMON/MMDATA/STPPEN,TEMMMT,TEMMMR,TEMMM

      WRITE(*,*)
      WRITE(*,*) '-----------------------------------------------------
     *-----------------'
      WRITE(*,*)
      WRITE(*,*) 'Solute - Solvent interaction energy     [a.u.]'
      WRITE(*,*) '    SSE        = ',SSE
      WRITE(*,*)
      WRITE(*,*) 'Solvent - Solvent interaction energy     [kcal/mol]'
      STPPEN = STPPEN *DBLE( NSP ) / 4.184D3
      WRITE(*,*) '    STPPE      = ',STPPEN
      WRITE(*,*)      
      WRITE(*,*) 'Translational Temp.  = ',TEMMMT
      WRITE(*,*) 'Rotational Temp.     = ',TEMMMR
      WRITE(*,*) 'Total Temp.          = ',TEMMM
      WRITE(*,*)      
      WRITE(*,*) '-----------------------------------------------------
     *-----------------'

      RETURN
      END

C--------------------------------------------------
C
C   SUBROUTINE: 
C
C     Optimizing fractional charges on each atom
C
C     by least square fittings
C
C--------------------------------------------------

      SUBROUTINE OPTFC( IMD )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NNUC  = 36 )
C     PARAMETER ( NORA  = 49 )
C     PARAMETER ( NORB  =  1 )
C     PARAMETER ( MORA =  49 )
C     PARAMETER ( MORB =   1 )
C     PARAMETER ( NDIM  =  3 )
      PARAMETER ( RCUT  = 4.D0 )
C     PARAMETER ( NMAX  = 80 )
C     PARAMETER ( NTRY  = 3000  )
      PARAMETER ( NTRY  = 6400*NX*NY*NZ )
      PARAMETER ( ALPHA = 0.1D0 )
   
      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          FREEZE*5,PRINT*5,DGF*3

C     DIMENSION COORD( NTRY,3 )
      DIMENSION COORD( 3,NTRY )
      DIMENSION VALUE( NTRY )
      DIMENSION SIGMA( NTRY )
      DIMENSION TCOORD( 3,NTRY/(NX*NY*NZ) )
      DIMENSION TVALUE( NTRY/(NX*NY*NZ) )
      DIMENSION TSIGMA( NTRY/(NX*NY*NZ) )
C     DIMENSION A(     NNUC )

      DIMENSION IRCNT(  0:NX*NY*NZ-1 )
      DIMENSION IRCNT1( 0:NX*NY*NZ-1 )
      DIMENSION IDISP(  0:NX*NY*NZ-1 )
      DIMENSION IDISP1( 0:NX*NY*NZ-1 )
      
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / NCLR / PNUC( NNUC,NDIM )
C     COMMON / OFC / VCOU( NMAX,NMAX,NMAX )
      COMMON / OFC  / VCOU( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      COMMON / OFC1 / A( NNUC )                 ! Vmol
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,                 ! Vmol
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,                 ! Vmol
C    *                  WNUC(NFUZZY,NMAX,NMAX,NMAX),RSIZE(NNUC),           ! Vmol
C    *                  WNUC(NFUZZY,NMAX/NX,NMAX/NY,NMAX/NZ),
     *                  RSIZE(NNUC),
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF                       ! Vmol

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      STIME = MPI_WTIME()

      IF(MYID.EQ.0 ) THEN
        WRITE(*,*) 
        WRITE(*,*) 'Optimization of Fractional Charges'
      ENDIF

      RDX2  = (2.D0*DX)**2*3
      HUGE  = DX*NMAXX

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      NCOUNT = 0

      IF( IMD.EQ.0 ) THEN
        IF(MYID.EQ.0) write(*,*) 'Initializing start point of RANFQ'
        DO II = 1,MYID
          NX1 = INT(RANFQ(DUMMY))
          NY1 = INT(RANFQ(DUMMY))
          NZ1 = INT(RANFQ(DUMMY))
        ENDDO
      ENDIF

      NTRY1 = NTRY/NUMPROCS
C     DO 10 I = 2, NTRY
      DO 10 I = 1, NTRY1

        NX1 = INT( RANFQ(DUMMY)*JX )+1  
        NY1 = INT( RANFQ(DUMMY)*JY )+1  
        NZ1 = INT( RANFQ(DUMMY)*JZ )+1  

        XX = DX*( NX1-NNX )
        YY = DY*( NY1-NNY )
        ZZ = DZ*( NZ1-NNZ )

        ZCLM = 0.D0
        RMIN = HUGE 
        DO 11 J = 1, NNUC
          IF( SIG(J) .EQ. 0.D0 ) THEN
          RCUT2=RCUT**2
          ELSE
C         RCUT2 = SIG(J)**2
          RCUT2=RCUT**2
          ENDIF
          RX = XX - PNUC( J,1 ) 
          RY = YY - PNUC( J,2 ) 
          RZ = ZZ - PNUC( J,3 ) 
          RR = RX**2 + RY**2 + RZ**2
C         IF( RR .LT. RCUT2 ) GOTO 10
          IF( RR .LT. RCUT2 ) GOTO 14
          R  = DSQRT( RR ) 
          ZCLM = ZCLM + ZVAL(J)/R
          IF( R .LT. RMIN  ) RMIN = R
11      CONTINUE      

C       DO 12 K = 1, I-1 
        DO 12 K = 1, NCOUNT
          RX = XX - TCOORD(1,K)
          RY = YY - TCOORD(2,K)
          RZ = ZZ - TCOORD(3,K)
          RR = RX**2 + RY**2 + RZ**2
C         IF( RR .LT. RDX2 ) GOTO 10
          IF( RR .LT. RDX2 ) GOTO 14
12      CONTINUE      

      NCOUNT = NCOUNT + 1  
      TCOORD( 1,NCOUNT ) = XX 
      TCOORD( 2,NCOUNT ) = YY
      TCOORD( 3,NCOUNT ) = ZZ

      TVALUE( NCOUNT )   = ZCLM - VCOU( NX1,NY1,NZ1 ) 
      TSIGMA( NCOUNT )   = 1/(1-DEXP( -ALPHA*RMIN ))

14    CONTINUE      

      DO II = 1,NUMPROCS-1
        NX1 = INT(RANFQ(DUMMY))
        NY1 = INT(RANFQ(DUMMY))
        NZ1 = INT(RANFQ(DUMMY))
      ENDDO

10    CONTINUE      

C     NDATA = NCOUNT

      DO ID = 0, NUMPROCS-1
        IRCNT1(ID) = 0
      ENDDO
      CALL MPI_GATHER( NCOUNT,1,MPI_INTEGER,
     *                 IRCNT1(0),1,MPI_INTEGER,
     *                 0,MPI_COMM_WORLD,IERR )

      IF(MYID.EQ.0) THEN
        NDATA = IRCNT1(0)
        IRCNT(0) = 3*IRCNT1(0)
        IDISP(0) = 0 
        IDISP1(0) = 0 
        DO ID = 1, NUMPROCS-1
          IRCNT1(ID) = IRCNT1(ID)
          IRCNT(ID)  = 3*IRCNT1(ID)
          IDISP(ID)  = IDISP(ID-1)+IRCNT(ID-1)
          IDISP1(ID) = IDISP1(ID-1)+IRCNT1(ID-1)
          NDATA      = NDATA + IRCNT1(ID)
        ENDDO
      ENDIF

      NUMDAT = NCOUNT*3
      CALL MPI_GATHERV( TCOORD,NUMDAT,MPI_DOUBLE_PRECISION,
     *                  COORD,IRCNT,IDISP,MPI_DOUBLE_PRECISION,
     *                  0,MPI_COMM_WORLD,IERR )
      CALL MPI_GATHERV( TVALUE,NCOUNT,MPI_DOUBLE_PRECISION,
     *                  VALUE,IRCNT1,IDISP1,MPI_DOUBLE_PRECISION,
     *                  0,MPI_COMM_WORLD,IERR )
      CALL MPI_GATHERV( TSIGMA,NCOUNT,MPI_DOUBLE_PRECISION,
     *                  SIGMA,IRCNT1,IDISP1,MPI_DOUBLE_PRECISION,
     *                  0,MPI_COMM_WORLD,IERR )

      IF(MYID.EQ.0 ) THEN
 
      WRITE(*,90) NDATA,NTRY
      
90    FORMAT(1X,I4,' points have been accepted'/
     *       1X,'out of ',I4,' trials.')

      CALL LSFIT( COORD,VALUE,SIGMA,NDATA,A,CHISQ,NTRY )

      SUM = 0.D0
      DO 20 I=1,NNUC
        SUM = SUM + A(I)
20    CONTINUE
      WRITE(*,*) 

      IF( EXC .EQ. 'RPZ' .OR. EXC .EQ. 'RBLYP' .OR. EXC.EQ.'RHF' 
     *                   .OR. EXC .EQ. 'RXalpha' ) THEN
C       QSUM = DBLE( 2.D0 * NORA )
        QSUM = DBLE( 2.D0 * MORA )
      ELSE
        QSUM = DBLE( MORA + MORB )
      ENDIF
      ZSUM = 0.D0
      DO 30 I=1,NNUC
        ZSUM = ZSUM + ZVAL(I)
30    CONTINUE

      QSUM = ZSUM - QSUM 
      RES = QSUM - SUM
      DO 40 I=1,NNUC
        QSHIFT = RES * ZVAL(I)/ZSUM
        A(I) = A(I) + QSHIFT 
40    CONTINUE

      WRITE(*,*) '----- Fractional Charges -----'

      WRITE(*,*)
      WRITE(*,93)   ! Vmol
      DO 50 I=1,NNUC
C       WRITE(*,91) I,A(I)        
        WRITE(*,91) I,A(I),ZVAL(I)-ZPOP(I)        
C91    FORMAT(3X,'ATOM',I4,' : ',F14.10)
91    FORMAT(3X,'ATOM',I4,' : ',F14.10,2x,F14.10)   ! Vmol
93    FORMAT(18x,'OPTFC',9x,'FUZZYCELL')   ! Vmol
50    CONTINUE
      WRITE(*,*) 
      WRITE(*,92) CHISQ
92    FORMAT(3X,'CHISQ    = ',F14.10)
      WRITE(*,*) 
      WRITE(*,*) '------------------------------'

      ENDIF

      ETIME = MPI_WTIME()
      IF(MYID.EQ.0) THEN
        write(*,*) 'elapsed time:OPTFC =',ETIME-STIME
      ENDIF

C     WRITE(35,*) IMD
C     DO NA = 1,NNUC
C       WRITE(35,95) ZA(NA),ZVAL(NA),A(NA)
C     ENDDO

C95    FORMAT(F6.1,F6.1,F14.10)

      RETURN
      END

C***************************************************

      SUBROUTINE OPTFC1(IMD)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol

      PARAMETER ( RCUT  = 4.D0  )
      PARAMETER ( NTRY  = 6400*NX*NY*NZ )
      PARAMETER ( ALPHA = 0.1D0 )
   
      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          FREEZE*5,PRINT*5,DGF*3

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / OFC  / VCOU( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      COMMON / OFC1 / A( NNUC )                 ! Vmol
      COMMON / OFC2 / NX1(NTRY/(NX*NY*NZ)),NY1(NTRY/(NX*NY*NZ)),
     *                NZ1(NTRY/(NX*NY*NZ)),NDATA,NDATA1,
     *                TCOORD( 3,NTRY/(NX*NY*NZ) ),
     *                COORD( 3,NTRY ),CHISQ
      COMMON / OFC3 / IRCNT( 0:NX*NY*NZ-1 ),IRCNT1( 0:NX*NY*NZ-1 ),        ! Vmol
     *                IDISP( 0:NX*NY*NZ-1 ),IDISP1( 0:NX*NY*NZ-1 )
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,                 ! Vmol
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,                 ! Vmol
C    *                  WNUC(NFUZZY,NMAX,NMAX,NMAX),RSIZE(NNUC),           ! Vmol
     *                  RSIZE(NNUC),                                       ! Vmol
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF                       ! Vmol

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      RDX2  = (2.D0*DX)**2*3
      HUGE  = DX*NMAXX

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      IF( IMD.EQ.0 ) THEN
        IF(MYID.EQ.0) write(*,*) 
        IF(MYID.EQ.0) write(*,*) 'Initializing start point of RANFQ'
        DO II = 1,MYID
          NX2 = INT(RANFQ(DUMMY))
          NY2 = INT(RANFQ(DUMMY))
          NZ2 = INT(RANFQ(DUMMY))
        ENDDO
      ENDIF

      NTRY1 = NTRY/NUMPROCS    ! 20260903
!     NTRY1 = NTRY             ! 20260902
      NCOUNT = 0

      DO 10 I = 1, NTRY1

        NX2 = INT( RANFQ(DUMMY)*JX )+1  
        NY2 = INT( RANFQ(DUMMY)*JY )+1  
        NZ2 = INT( RANFQ(DUMMY)*JZ )+1  

        XX = DX*( NX2-NNX )
        YY = DY*( NY2-NNY )
        ZZ = DZ*( NZ2-NNZ )

        ZCLM = 0.D0
        RMIN = HUGE 
        DO 11 J = 1, NNUC
          IF( SIG(J) .EQ. 0.D0 ) THEN
            RCUT2=RCUT**2
          ELSE
C           RCUT2 = SIG(J)**2
            RCUT2=RCUT**2
          ENDIF
          RX = XX - PNUC( J,1 ) 
          RY = YY - PNUC( J,2 ) 
          RZ = ZZ - PNUC( J,3 ) 
          RR = RX**2 + RY**2 + RZ**2
          IF( RR .LT. RCUT2 ) GOTO 10
11      CONTINUE      

        DO 12 K = 1, NCOUNT
          RX = XX - TCOORD(1,K)
          RY = YY - TCOORD(2,K)
          RZ = ZZ - TCOORD(3,K)
          RR = RX**2 + RY**2 + RZ**2
          IF( RR .LT. RDX2 ) GOTO 10
12      CONTINUE      

      NCOUNT = NCOUNT + 1  
      TCOORD( 1,NCOUNT ) = XX 
      TCOORD( 2,NCOUNT ) = YY
      TCOORD( 3,NCOUNT ) = ZZ
      NX1(NCOUNT)       = NX2
      NY1(NCOUNT)       = NY2
      NZ1(NCOUNT)       = NZ2

10    CONTINUE      

      NDATA1 = NCOUNT

      CALL MPI_GATHER( NCOUNT,1,MPI_INTEGER,
     *                 IRCNT1(0),1,MPI_INTEGER,
     *                 0,MPI_COMM_WORLD,IERR )

      IF(MYID.EQ.0) THEN
        NDATA = IRCNT1(0)
        IRCNT(0) = 3*IRCNT1(0)
        IDISP(0) = 0 
        IDISP1(0) = 0 
        DO ID = 1, NUMPROCS-1
          IRCNT(ID)  = 3*IRCNT1(ID)
          IDISP(ID)  = IDISP(ID-1)+IRCNT(ID-1)
          IDISP1(ID) = IDISP1(ID-1)+IRCNT1(ID-1)
          NDATA      = NDATA + IRCNT1(ID)
        ENDDO
        WRITE(*,*)
        WRITE(*,90) NDATA,NTRY
      ENDIF

      NUMDAT = NCOUNT*3
      CALL MPI_GATHERV( TCOORD,NUMDAT,MPI_DOUBLE_PRECISION,
     *                  COORD,IRCNT,IDISP,MPI_DOUBLE_PRECISION,
     *                  0,MPI_COMM_WORLD,IERR )

90    FORMAT(1X,I6,' points have been accepted'/
     *       1X,'out of ',I6,' trials.')

      RETURN
      END

      SUBROUTINE OPTFC2                         ! For Poisson Equation

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol

      PARAMETER ( RCUT  = 4.D0 )
      PARAMETER ( NTRY  = 6400*NX*NY*NZ )
      PARAMETER ( ALPHA = 0.1D0 )
   
      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          FREEZE*5,PRINT*5,DGF*3

      DIMENSION VALUE(  NTRY )
      DIMENSION SIGMA(  NTRY )
      DIMENSION TVALUE( NTRY/(NX*NY*NZ) )
      DIMENSION TSIGMA( NTRY/(NX*NY*NZ) )
      
      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / OFC  / VCOU( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      COMMON / OFC1 / A( NNUC )                 ! Vmol
      COMMON / OFC2 / NX1(NTRY/(NX*NY*NZ)),NY1(NTRY/(NX*NY*NZ)),
     *                NZ1(NTRY/(NX*NY*NZ)),NDATA,NDATA1,
     *                TCOORD( 3,NTRY/(NX*NY*NZ) ),
     *                COORD( 3,NTRY ),CHISQ
      COMMON / OFC3 / IRCNT( 0:NX*NY*NZ-1 ),IRCNT1( 0:NX*NY*NZ-1 ),        ! Vmol
     *                IDISP( 0:NX*NY*NZ-1 ),IDISP1( 0:NX*NY*NZ-1 )
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,                 ! Vmol
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,                 ! Vmol
C    *                  WNUC(NFUZZY,NMAX,NMAX,NMAX),RSIZE(NNUC),           ! Vmol
     *                  RSIZE(NNUC),                                       ! Vmol
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF                       ! Vmol

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      CHISQ = 0.D0

      RDX2  = (2.D0*DX)**2*3
      HUGE  = DX*NMAXX

      DO 10 I = 1, NDATA1

        XX = TCOORD(1,I)
        YY = TCOORD(2,I)
        ZZ = TCOORD(3,I)

        ZCLM = 0.D0
        RMIN = HUGE 
        DO 11 J = 1, NNUC
          RX = XX - PNUC( J,1 ) 
          RY = YY - PNUC( J,2 ) 
          RZ = ZZ - PNUC( J,3 ) 
          RR = RX**2 + RY**2 + RZ**2
          R  = DSQRT( RR ) 
          ZCLM = ZCLM + ZVAL(J)/R
          IF( R .LT. RMIN  ) RMIN = R
11      CONTINUE      

        TVALUE( I )   = ZCLM - VCOU(NX1(I),NY1(I),NZ1(I)) 
        TSIGMA( I )   = 1/(1-DEXP(-ALPHA*RMIN))

10    CONTINUE      

      CALL MPI_GATHERV( TVALUE,NDATA1,MPI_DOUBLE_PRECISION,
     *                  VALUE,IRCNT1,IDISP1,MPI_DOUBLE_PRECISION,
     *                  0,MPI_COMM_WORLD,IERR )
      CALL MPI_GATHERV( TSIGMA,NDATA1,MPI_DOUBLE_PRECISION,
     *                  SIGMA,IRCNT1,IDISP1,MPI_DOUBLE_PRECISION,
     *                  0,MPI_COMM_WORLD,IERR )

      IF(MYID.EQ.0 ) THEN

      CALL LSFIT( COORD,VALUE,SIGMA,NDATA,A,CHISQ,NTRY )

      SUM = 0.D0
      DO 20 I=1,NNUC
        SUM = SUM + A(I)
20    CONTINUE
      WRITE(*,*) 

      IF( EXC .EQ. 'RPZ' .OR. EXC .EQ. 'RBLYP' .OR. EXC.EQ.'RHF' 
     *                   .OR. EXC .EQ. 'RXalpha' ) THEN
        QSUM = DBLE( 2.D0 * MORA )
      ELSE
        QSUM = DBLE( MORA + MORB )
      ENDIF
      ZSUM = 0.D0
      DO 30 I=1,NNUC
        ZSUM = ZSUM + ZVAL(I)
30    CONTINUE

      QSUM = ZSUM - QSUM 
      RES = QSUM - SUM
      DO 40 I=1,NNUC
        QSHIFT = RES * ZVAL(I)/ZSUM
        A(I) = A(I) + QSHIFT 
        ZPOP(I) = ZVAL(I) - A(I)
40    CONTINUE

C       WRITE(*,91) (ZPOP(I),I=1,NNUC)   ! comment 2005.11.22 takahashi
91      FORMAT(2X,'ZPOP',' : ',5F10.5)   ! Vmol
        WRITE(*,92) CHISQ
92      FORMAT(2X,'CHISQ    = ',F14.10)
      ENDIF

      CALL MPI_BCAST( ZPOP(1),NNUC,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )

      RETURN
      END

      SUBROUTINE OPTFC3                         ! Print Fractional(ESP) Charges

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol

      PARAMETER ( RCUT  = 4.D0 )
      PARAMETER ( NTRY  = 6400*NX*NY*NZ )
      PARAMETER ( ALPHA = 0.1D0 )
   
      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          FREEZE*5,PRINT*5,DGF*3

      DIMENSION VALUE( NTRY )
      DIMENSION SIGMA( NTRY )
      DIMENSION TVALUE( NTRY/(NX*NY*NZ) )
      DIMENSION TSIGMA( NTRY/(NX*NY*NZ) )
      
      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / OFC  / VCOU( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      COMMON / OFC1 / A( NNUC )                 ! Vmol
      COMMON / OFC2 / NX1(NTRY/(NX*NY*NZ)),NY1(NTRY/(NX*NY*NZ)),
     *                NZ1(NTRY/(NX*NY*NZ)),NDATA,NDATA1,
     *                TCOORD( 3,NTRY/(NX*NY*NZ) ),
     *                COORD( 3,NTRY ),CHISQ
      COMMON / OFC3 / IRCNT( 0:NX*NY*NZ-1 ),IRCNT1( 0:NX*NY*NZ-1 ),        ! Vmol
     *                IDISP( 0:NX*NY*NZ-1 ),IDISP1( 0:NX*NY*NZ-1 )
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,                 ! Vmol
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,                 ! Vmol
C    *                  WNUC(NFUZZY,NMAX,NMAX,NMAX),RSIZE(NNUC),           ! Vmol
     *                  RSIZE(NNUC),                                       ! Vmol
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF                       ! Vmol

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      CHISQ = 0.D0

      RDX2  = (2.D0*DX)**2*3
      HUGE  = DX*NMAXX

      DO 10 I = 1, NDATA1

        XX = TCOORD(1,I)
        YY = TCOORD(2,I)
        ZZ = TCOORD(3,I)

        ZCLM = 0.D0
        RMIN = HUGE 
        DO 11 J = 1, NNUC
          RX = XX - PNUC( J,1 ) 
          RY = YY - PNUC( J,2 ) 
          RZ = ZZ - PNUC( J,3 ) 
          RR = RX**2 + RY**2 + RZ**2
          R  = DSQRT( RR ) 
          ZCLM = ZCLM + ZVAL(J)/R
          IF( R .LT. RMIN  ) RMIN = R
11      CONTINUE      

        TVALUE( I )   = ZCLM - VCOU(NX1(I),NY1(I),NZ1(I)) 
        TSIGMA( I )   = 1/(1-DEXP(-ALPHA*RMIN))

10    CONTINUE      

      CALL MPI_GATHERV( TVALUE,NDATA1,MPI_DOUBLE_PRECISION,
     *                  VALUE,IRCNT1,IDISP1,MPI_DOUBLE_PRECISION,
     *                  0,MPI_COMM_WORLD,IERR )
      CALL MPI_GATHERV( TSIGMA,NDATA1,MPI_DOUBLE_PRECISION,
     *                  SIGMA,IRCNT1,IDISP1,MPI_DOUBLE_PRECISION,
     *                  0,MPI_COMM_WORLD,IERR )

      IF(MYID.EQ.0 ) THEN

      CALL LSFIT( COORD,VALUE,SIGMA,NDATA,A,CHISQ,NTRY )

      SUM = 0.D0
      DO 20 I=1,NNUC
        SUM = SUM + A(I)
20    CONTINUE
      WRITE(*,*) 

      IF( EXC .EQ. 'RPZ' .OR. EXC .EQ. 'RBLYP' .OR. EXC.EQ.'RHF' 
     *                   .OR. EXC .EQ. 'RXalpha' ) THEN
        QSUM = DBLE( 2.D0 * MORA )
      ELSE
        QSUM = DBLE( MORA + MORB )
      ENDIF
      ZSUM = 0.D0
      DO 30 I=1,NNUC
        ZSUM = ZSUM + ZVAL(I)
30    CONTINUE

      QSUM = ZSUM - QSUM 
      RES = QSUM - SUM
      DO 40 I=1,NNUC
        QSHIFT = RES * ZVAL(I)/ZSUM
        A(I) = A(I) + QSHIFT 
40    CONTINUE

        WRITE(*,*) 
        WRITE(*,*) 'Optimization of Fractional Charges'

        WRITE(*,*) '----- Fractional Charges -----'

        WRITE(*,*)
        DO 50 I=1,NNUC
          WRITE(*,91) I,A(I)
91        FORMAT(3X,'ATOM',I4,' : ',F14.10)
50      CONTINUE
        WRITE(*,*) 
        WRITE(*,92) CHISQ
92      FORMAT(3X,'CHISQ    = ',F14.10)
        WRITE(*,*) 
        WRITE(*,*) '------------------------------'
      ENDIF

      RETURN
      END


C---------------------------------------- 
C
C  SUBROUTINE Basis Functions for 
C  Least Square Fittings
C
C  Basis: Coulombic potentials( 1/r )
C 
C---------------------------------------- 

      SUBROUTINE FUNCS( NTMP,COORD,AFUNC,NTRY )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include 'QMpara.i'                        ! Vmol

C     PARAMETER (NNUC  = 36 )
C     PARAMETER (NDIM  = 3 )
C     DIMENSION COORD(NTRY,3),AFUNC(NNUC)
      DIMENSION COORD(3,NTRY),AFUNC(NNUC)

      COMMON / NCLR / PNUC( NNUC,NDIM )

      DO 10 I = 1,NNUC

C       RX = COORD( NTMP,1 ) - PNUC( I,1 )
C       RY = COORD( NTMP,2 ) - PNUC( I,2 )
C       RZ = COORD( NTMP,3 ) - PNUC( I,3 )
        RX = COORD( 1,NTMP ) - PNUC( I,1 )
        RY = COORD( 2,NTMP ) - PNUC( I,2 )
        RZ = COORD( 3,NTMP ) - PNUC( I,3 )
        R2 = RX**2 + RY**2 + RZ**2
        R  = DSQRT(R2)

        AFUNC(I) = 1.D0/R

10    CONTINUE
     
      RETURN
      END

 
C     ######################################################
C     #   qm/mm boundary soubroutines for tinker package   #
C     ######################################################

C---------------------------------------------------
C     subroutine sorting QM / MM parameters
C---------------------------------------------------

      subroutine spqmmm(NMID1)

      implicit real*8 ( a-h,o-z )
      implicit integer*4 ( i-n )

      include "mpif.h"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
!     include 'atoms.i'

C     PARAMETER ( NNUC  = 36 )
C     PARAMETER ( NLINK1 = 10 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          FREEZE*5,PRINT*5,DGF*3

      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)

      dimension NOQM(NLINK1)
      dimension NOMM(NLINK1)
      dimension NMID(  NNUC)
      dimension NMID1( NNUC)
      dimension ZA1(   NNUC)
      dimension ZVAL1( NNUC)
      dimension TSIG(  NNUC)
      dimension TEPSQM(NNUC)

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      IF(MYID.EQ.0) THEN

      write(*,*)
      write(*,*) 'sorting parameters and coordinates'
C
C     sort QM/MM parameters
C
C
C     sort MMID
C
      IF(NQMMM.EQ.'LINK') THEN
        nlqm = NNUC-NLINK
      ELSEIF(NQMMM.EQ.'MM') THEN
C       nlqm = NNUC
        nlqm = NNUC-NLINK
      ELSE
        write(*,*) 'error: subroutine spqmmm'
        STOP
      ENDIF

      do i = 1, nlqm
        NMID(i)   = MMID(i)
        NMID1(i)  = MMID(i)
        ZA1(i)    = ZA(i)
        ZVAL1(i)  = ZVAL(i)
        TSIG(i)   = SIG(i)
        TEPSQM(i) = EPSQM(i)
      enddo

C     idmax = nlqm + 1

      do j = 1, nlqm
C       nid1 = n
        nid1 = n + 1                                        ! 2004.09.24
        ii   = 0                                            ! 2004.07.01
        do i = 1, nlqm

          if ( NMID(i) .le. nid1 .and. NMID(i) .lt. n+1 ) then ! 2004.09.24
C         if ( NMID(i) .le. nid1 ) then

            MMID(j)  = NMID(i)
            ZA(j)    = ZA1(i)
            ZVAL(j)  = ZVAL1(i)
            SIG(j)   = TSIG(i)
            EPSQM(j) = TEPSQM(i)

            if( MMID(j) .eq. nid1 .and. ii.gt.0 ) then      ! 2004.07.01
C           if( MMID(j) .eq. nid1 .and. ii.lt.0 ) then
C           if( MMID(j) .eq. nid1 ) then
              write(*,*)
     *          'ERROR: Two atoms(QM-A) have identical ID number.'
              write(*,98) 'Atom(QM) = ',MMID1,i
              write(*,99) 'MMID(QM) = ',MMID(j)
              STOP
            endif

            nid1    = NMID(i)
            MMID1   = i
            ii      = ii + 1                                ! 2004.07.01

          endif

        enddo

        NMID(MMID1) = n + 1

      enddo

      ENDIF

CC
CC     sort NLAQM and NLAMM
CC
C      do i = 1, NLINK
C        NOQM(i) = NLAQM(i)
C        NOMM(i) = NLAMM(i)
C      enddo
C
C      nqmmax = NLINK + 1
C      nmmmax = n + 1
C
C      do j = 1, NLINK
C        nqm1  = nqmmax
C        nmm1  = nmmmax
C        do i = 1, NLINK
CC
CC        -- for QM part --
CC
C          if ( NOQM(i) .le. nqm1 .and. NOQM(i) .ne. nqmmax ) then
C            NLAQM(j) = NOQM(i)
C            if( NLAQM(j) .eq. nqm1 ) then
C              write(*,98)
C     *          'ERROR: Two atoms(QM-B) have identical ID number.',
C     *           NOQM1,i
C              STOP
C            endif
C            nqm1     = NOQM(i)
C            NOQM1    = i
C          endif
CC
CC        -- for MM part --
CC
C          if ( NOMM(i) .le. nmm1 .and. NOMM(i) .ne. nmmmax ) then
C            NLAMM(j) = NOMM(i)
C            if( NLAMM(j) .eq. nmm1 ) then
C              write(*,98)
C     *          'ERROR: Two atoms(MM) have identical ID number.',
C     *           NOMM1,i
C              STOP
C            endif
C            nmm1     = NOMM(i)
C            NOMM1    = i
C          endif
C
C        enddo
C
C        NOQM(NOQM1) = nqmmax
C        NOMM(NOMM1) = nmmmax
C
C      enddo

98    format(a,2i6)
99    format(a,i6)

      return
      end

C---------------------------------------------------
C     subroutine translate from angstrom to bohr
C---------------------------------------------------

!     subroutine trntin( scrd,crd )

!     implicit real*8 ( a-h,o-z )
!     implicit integer*4 ( i-n )

!     include "mpif.h"

!     include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
!     include 'atoms.i'
!     include 'units.i'

!     CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
!    *          FREEZE*5,PRINT*5,DGF*3

!     COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
!    *                 FREEZE,PRINT,DGF
!     COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
!    *                 NMMSW(maxatm),MMID(NNUC)

!     dimension scrd( ndim,maxatm )
!     DIMENSION CRD( NNUC,NDIM )

!     CALL MPI_BCAST( x(1),n,MPI_DOUBLE_PRECISION,0,
!    *                MPI_COMM_WORLD,IERR )
!     CALL MPI_BCAST( y(1),n,MPI_DOUBLE_PRECISION,0,
!    *                MPI_COMM_WORLD,IERR )
!     CALL MPI_BCAST( z(1),n,MPI_DOUBLE_PRECISION,0,
!    *                MPI_COMM_WORLD,IERR )

!     bohri = 1.d0/bohr

!     ne1 = NNUC - NLINK

!     ns2 = 1
!     ncount = 0
!     do j = 1, ne1
!       MMID1 = MMID(j)
!       ne2   = MMID1 - 1
!       do i = ns2, ne2
!         ncount = ncount + 1
!         scrd(1,ncount) = x(i)*bohri
!         scrd(2,ncount) = y(i)*bohri
!         scrd(3,ncount) = z(i)*bohri
!       enddo
!       crd(j,1) = x(MMID1)*bohri
!       crd(j,2) = y(MMID1)*bohri
!       crd(j,3) = z(MMID1)*bohri
!       ns2 = MMID1 + 1
!     enddo
!     do i = ns2, n
!       ncount = ncount + 1
!       scrd(1,ncount) = x(i)*bohri
!       scrd(2,ncount) = y(i)*bohri
!       scrd(3,ncount) = z(i)*bohri
!     enddo

C     ENDIF

!     return
!     end

C----------------------------------------------
C     translation of force on each MM site
C     from hartree/bohr to kcal/mol/Ang
C----------------------------------------------

      subroutine trntin1( FRCS,FRC )

      implicit real*8 ( a-h,o-z )
      implicit integer*4 ( i-n )

      include "mpif.h"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
!     include 'atoms.i'
!     include 'units.i'

C     PARAMETER ( NNUC  = 36 )
C     PARAMETER ( NDIM  = 3 )
C     PARAMETER ( NLINK1 = 10 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          FREEZE*5,PRINT*5,DGF*3

      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)

      DIMENSION FRC( NNUC,NDIM )
      DIMENSION FRCS( maxatm,NDIM )               !QM/MM force on site ( for tinker )
      DIMENSION FRCS1( maxatm,NDIM )              !QM/MM force on site ( for Vmol )

      coeff = hartree / bohr

      ne3 = n - NNUC + NLINK

      do i = 1, ne3 
        FRCS1(i,1) = FRCS(i,1) * coeff
        FRCS1(i,2) = FRCS(i,2) * coeff
        FRCS1(i,3) = FRCS(i,3) * coeff
      enddo

      ne1 = NNUC - NLINK

      ns2 = 1
      ncount = 0
      do j = 1, ne1
        MMID1 = MMID(j)
        ne2   = MMID1 - 1
        do i = ns2, ne2
          ncount = ncount + 1
          FRCS(i,1) = FRCS1(ncount,1)
          FRCS(i,2) = FRCS1(ncount,2)
          FRCS(i,3) = FRCS1(ncount,3)
C       write(*,11) 'trntin1',i,FRCS(i,1),FRCS(i,2),FRCS(i,3)   ! test
        enddo
C       write(*,10) 'trntin1',MMID1,j,FRC(j,1),FRC(j,2),FRC(j,3)   ! test
        FRCS(MMID1,1) = FRC(j,1) * coeff
        FRCS(MMID1,2) = FRC(j,2) * coeff
        FRCS(MMID1,3) = FRC(j,3) * coeff
        ns2 = MMID1 + 1
      enddo
      do i = ns2, n
        ncount = ncount + 1
        FRCS(i,1) = FRCS1(ncount,1)
        FRCS(i,2) = FRCS1(ncount,2)
        FRCS(i,3) = FRCS1(ncount,3)
C       write(*,11) 'trntin1',i,FRCS(i,1),FRCS(i,2),FRCS(i,3)   ! test
      enddo

10    format(a,2i6,3f12.6)
11    format(a,i6,3f12.6)

      return
      end

C------------------------------------------------------
C     determination of link atom position by SPLAM
C------------------------------------------------------

      subroutine LINKA_FIX( CRD )

      implicit real*8 ( a-h,o-z )
      implicit integer*4 ( i-n )

      include "mpif.h"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'

C     PARAMETER ( NNUC  = 36 )
C     PARAMETER ( NDIM  = 3 )
C     PARAMETER ( NLINK1 = 10 )

      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PRMT3 / blch(NLINK1),blcc(NLINK1),
     *                 bkch(NLINK1),bkcc(NLINK1)
      COMMON / MMDAT / SCRD( NDIM,maxatm ),SCRD1( NDIM,maxatm ),   ! Vmol
     *                 FRCS( maxatm,NDIM )                         ! Vmol

      DIMENSION CRD( NNUC,NDIM )
      DIMENSION CRDL( NLINK1*NDIM )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
C     CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      ns = NNUC - NLINK + 1

      IF(MYID.EQ.0) THEN

      write(*,*)
      write(*,*) 'Read link.dat file...'

      rewind(85)
      do i = ns, NNUC          ! Link atoms will be fixed at their average positions
        read(85,*)  id,inum,(crd(i,j),j=1,3)   
      enddo

      write(*,*) 'Link Atom Position'

      nc = 0

      NCOUNT = 0

      do i = ns, NNUC

C       nc = nc + 1

C       xcc = scrd(1,NLAMM(nc)) - crd(NLAQM(nc),1)
C       ycc = scrd(2,NLAMM(nc)) - crd(NLAQM(nc),2)
C       zcc = scrd(3,NLAMM(nc)) - crd(NLAQM(nc),3)
C       rcc = xcc*xcc + ycc*ycc + zcc*zcc
C       rcc = dsqrt(rcc)
C       xcc = xcc/rcc
C       ycc = ycc/rcc
C       zcc = zcc/rcc

C       rch = blch(nc) + (bkcc(nc)/bkch(nc))*(rcc-blcc(nc))

C       crd(i,1) = crd(NLAQM(nc),1) + rch*xcc
C       crd(i,2) = crd(NLAQM(nc),2) + rch*ycc
C       crd(i,3) = crd(NLAQM(nc),3) + rch*zcc

        do j = 1, NDIM
          NCOUNT = NCOUNT + 1
          CRDL(NCOUNT) = crd(i,j)
        enddo

        write(*,99) ZA(i),ZVAL(i),CRD(i,1),CRD(i,2),CRD(i,3)

      enddo

      CALL MPI_BCAST( CRDL(1),NCOUNT,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )

      ELSE

      NDAT = NLINK*NDIM  ! NDAT must be same as NCOUNT of rank0
      CALL MPI_BCAST( CRDL(1),NDAT,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )

      NCOUNT = 0
      do i = ns, NNUC
      do j = 1,  NDIM
        NCOUNT = NCOUNT + 1
        crd(i,j) = CRDL(NCOUNT)
      enddo
      enddo

      ENDIF

99    format( 2f6.1,3f12.6 )

      return
      end



C------------------------------------------------------
C     determination of link atom position by SPLAM
C------------------------------------------------------

      subroutine LINKA( CRD )

      implicit real*8 ( a-h,o-z )
      implicit integer*4 ( i-n )

      include "mpif.h"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'

C     PARAMETER ( NNUC  = 36 )
C     PARAMETER ( NDIM  = 3 )
C     PARAMETER ( NLINK1 = 10 )

      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PRMT3 / blch(NLINK1),blcc(NLINK1),
     *                 bkch(NLINK1),bkcc(NLINK1)
      COMMON / MMDAT / SCRD( NDIM,maxatm ),SCRD1( NDIM,maxatm ),   ! Vmol
     *                 FRCS( maxatm,NDIM )                         ! Vmol

      DIMENSION CRD( NNUC,NDIM )
      DIMENSION CRDL( NLINK1*NDIM )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
C     CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      ns = NNUC - NLINK + 1

      IF(MYID.EQ.0) THEN

      write(*,*)
      write(*,*) 'Link Atom Position'

C     do i = ns, NNUC          ! Link atoms will be fixed at their average positions
C       read(85,*)  id,inum,(crd(i,j),j=1,3)   
C     enddo

      nc = 0

      NCOUNT = 0

      do i = ns, NNUC

        nc = nc + 1

        xcc = scrd(1,NLAMM(nc)) - crd(NLAQM(nc),1)
        ycc = scrd(2,NLAMM(nc)) - crd(NLAQM(nc),2)
        zcc = scrd(3,NLAMM(nc)) - crd(NLAQM(nc),3)
        rcc = xcc*xcc + ycc*ycc + zcc*zcc
        rcc = dsqrt(rcc)
        xcc = xcc/rcc
        ycc = ycc/rcc
        zcc = zcc/rcc

        rch = blch(nc) + (bkcc(nc)/bkch(nc))*(rcc-blcc(nc))

        crd(i,1) = crd(NLAQM(nc),1) + rch*xcc
        crd(i,2) = crd(NLAQM(nc),2) + rch*ycc
        crd(i,3) = crd(NLAQM(nc),3) + rch*zcc

        do j = 1, NDIM
          NCOUNT = NCOUNT + 1
          CRDL(NCOUNT) = crd(i,j)
        enddo

        write(*,99) ZA(i),ZVAL(i),CRD(i,1),CRD(i,2),CRD(i,3)

      enddo

      CALL MPI_BCAST( CRDL(1),NCOUNT,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )

      ELSE

      NDAT = NLINK*NDIM  ! NDAT must be same as NCOUNT of rank0
      CALL MPI_BCAST( CRDL(1),NDAT,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )

      NCOUNT = 0
      do i = ns, NNUC
      do j = 1,  NDIM
        NCOUNT = NCOUNT + 1
        crd(i,j) = CRDL(NCOUNT)
      enddo
      enddo

      ENDIF

99    format( 2f6.1,3f12.6 )

      return
      end


C-------------------------------------------------
C
C     SUBROUTINE EDF(RHO,SCRD,PENG,IMD)
C
C     calculation of energy distribution function
C     in the solution system  
C     for the use with Tinker package 
C
C-------------------------------------------------

C     IMPLICIT REAL*8( A-H,O-Z ) 
C     IMPLICIT INTEGER*4( I-N )

C     include 'mpif.h'
C     include 'mpi.i'

C     INCLUDE 'QMpara.i'                        ! Vmol  
C     INCLUDE 'sizes.i'
C     INCLUDE 'atoms.i'
C     include 'atmtyp.i'

C     PARAMETER ( NCUT  = 0    )             ! cut step
C     PARAMETER ( ROMG  =   35.D0 )             ! solute-solvent interaction of finite range (Angstrom)

C-----/home3/takahasi/1b4v/FAD/mpi/fad/Vmol-tinker-----
C     PARAMETER ( EISO  = -131.676068149D0 )    ! energy of FAD with N electrons at isolation

C     PARAMETER ( NDXL  =  500    )             ! number of grids for the energy coordinate
C     PARAMETER ( EMIN  =  -60.D0 )             ! unit of kcal/mol 
C     PARAMETER ( EMAX  =   20.D0 )             ! unit of kcal/mol 
C     
C     PARAMETER ( ALPHA = 1.0D0  )      
C     PARAMETER ( AUE   = 627.0898D0 )      
C     PARAMETER ( AUL   = 0.529177D0 )      

C     PARAMETER ( maxsol = 3000 )
C     PARAMETER ( maxsit =   20 )

C     COMMON / MMPI1  / MX,MY,MZ
C     COMMON / MMPI4  / JX,JY,JZ

C     COMMON / GRID2  / DX,DY,DZ
C     COMMON / PRMT   / COE,TEMP,DTMD,
C    *                  NDEN,MDMAX,
C    *                  ZA(NNUC),ZVAL(NNUC),
C    *                  SIG(NNUC),EPSQM(NNUC),
C    *                  NRVLC,NCHK,CONV
C     COMMON / PRMT2  / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
C    *                  NMMSW(maxatm),MMID(NNUC)
C     COMMON / PRMT4  / nion1,iiont(maxatm)
C     COMMON / NCLR   / PNUC( NNUC,NDIM )
C     COMMON / SLT    / VPCE(NMAXX/NX,NMAXY/NY,NMAXZ/NZ),VPCZ,PLJ
C     COMMON / LJQMMM / SQMMM(NNUC,maxatm),EPQMMM(NNUC,maxatm),
C    *                  chgmm(maxatm),chgqm(NNUC)

C     COMMON / PRMTEDF / zmm(maxatm),nsol,nsite(maxsol),iiont1          ! 2005.08.31 takahashi
C     COMMON / ZMASS / ZCM(100)

C     DIMENSION iiont1(maxsol,maxsit)

C     DIMENSION RHO  ( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
C     DIMENSION SCRD ( NDIM,maxatm )
C     DIMENSION CRD  ( NNUC,3 )
C     DIMENSION NECHI( NDXL )
C     DIMENSION ETMP ( NDXL )
C     DIMENSION BLCHI( NDXL,NDXL )
C     DIMENSION NEDFQM ( NDXL )
C     DIMENSION NEDFMM ( NDXL )
C     DIMENSION NERDF  ( NDXL )

C     CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
C     CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

C     DE    = (EMAX-EMIN)/DBLE(NDXL)   ! kcal/mol
C     DE    =  DE/AUE                  ! atomic unit
C     EMINA =  EMIN/AUE                ! atomic unit
C     EMAXA =  EMAX/AUE                ! atomic unit
C     ROMGA =  ROMG/AUL                ! atomic unit

C     DV    =  DX*DY*DZ

C     BB = (DLOG10(EMAX/ECORE))/NM
C     DV = 8.D0*DX*DY*DZ

C     IF(MYID.EQ.0) THEN
C       WRITE(*,*) 'iiont1(30,1)=',iiont1(30,1)
C     ENDIF

C     NMAXHX = NMAXX/2
C     NMAXHY = NMAXY/2
C     NMAXHZ = NMAXZ/2

C     NNMAXX = (1-MX)*NMAXHX + 1
C     NNMAXY = (1-MY)*NMAXHY + 1
C     NNMAXZ = (1-MZ)*NMAXHZ + 1

C-----compute center of mass of QM atoms------

C     CX   = 0.D0
C     CY   = 0.D0
C     CZ   = 0.D0
C     ZSUM = 0.D0

C     DO NN = 1,NNUC
C       NZA  = INT(ZA(NN))
C       ZSUM = ZSUM + ZCM(NZA)
C       CX   = CX + ZCM(NZA)*PNUC(NN,1)
C       CY   = CY + ZCM(NZA)*PNUC(NN,2)
C       CZ   = CZ + ZCM(NZA)*PNUC(NN,3)
C     ENDDO
C     CMASX = CX/ZSUM
C     CMASY = CY/ZSUM
C     CMASZ = CZ/ZSUM

C-----initialize---------------

C     IF( IMD .EQ. NCUT ) THEN
C     
C     DO I = 1,NDXL
C      NEDFQM(I) = 0          ! energy distribution for QM energy
C      NEDFMM(I) = 0          ! energy distribution for QMMM interaction
C      NERDF( I) = 0          ! total energy distribution (QM + QM/MM)
C     ENDDO

C     ENDIF

C------------------------------

C     IF( IMD .GT. NCUT ) THEN

C     EPC  = 0.D0             ! for QM-energy coordinate
C     DO K = 1,JZ
C     DO J = 1,JY
C     DO I = 1,JX
C       EPC = EPC + VPCE(I,J,K)*RHO(I,J,K)*DV
C     ENDDO
C     ENDDO
C     ENDDO

C     CALL MPI_REDUCE(EPC,SEPC,1,MPI_DOUBLE_PRECISION,
C    *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
C     EPC = SEPC
C   
C     CALL MPI_BCAST(EPC,1,MPI_DOUBLE_PRECISION,0,
C    *               MPI_COMM_WORLD,IERR)


C-----solute-solvent interaction energy calculation----------------

C     SSE = 0.D0
C     DO 10 NSL  = 1,nsol

C     SUMM = 0.D0 
C     CX   = 0.D0
C     CY   = 0.D0
C     CZ   = 0.D0

C     DO NS = 1,nsite(NSL)

C       in   = iiont1(NSL,NS)   
C       IF(NMMSW(in).EQ.1) GOTO 10 

C       mn = in + 22
C       SUMM = SUMM + zmm(mn)
C       CX   = CX   + zmm(mn)*SCRD(1,in)
C       CY   = CY   + zmm(mn)*SCRD(2,in)
C       CZ   = CZ   + zmm(mn)*SCRD(3,in)

C     ENDDO

C     CMMX = CX/SUMM
C     CMMY = CY/SUMM
C     CMMZ = CZ/SUMM

C     RX = CMASX - CMMX
C     RY = CMASY - CMMY
C     RZ = CMASZ - CMMZ
C     R  = DSQRT(RX**2+RY**2+RZ**2)

C-----solute-solvent interaction of finite range-------------------

C     IF(R .LT. ROMGA) THEN

C     CNS   = 0.D0
C     POTLJ = 0.D0
C     CELS  = 0.D0

C     DO 15 NS = 1,nsite(NSL)

C     in = iiont1(NSL,NS)    

C-----nuclear-site interaction-------------------------------------

C     IF(MYID.EQ.0) THEN

C     DO 20 M = 1,NNUC

C      RX = PNUC(M,1) - SCRD(1,in)
C      RY = PNUC(M,2) - SCRD(2,in)
C      RZ = PNUC(M,3) - SCRD(3,in)
C      R  = DSQRT(RX**2+RY**2+RZ**2)
C      RI = 1.D0/R
C      CNS = CNS + ZVAL( M )*chgmm( in )*RI   ! Coulomb (nuclear-site)

C20    CONTINUE

C-----LJ potential--------------------------------------------------

C     DO 30 M = 1,NNUC

C      RX = PNUC(M,1) - SCRD(1,in)
C      RY = PNUC(M,2) - SCRD(2,in)
C      RZ = PNUC(M,3) - SCRD(3,in)
C      R  = DSQRT(RX**2+RY**2+RZ**2)
C      RI = 1.D0/R

C      sigr6  = (SQMMM(M,in)*RI*RI)**3
C      sigr12 =  sigr6*sigr6 

C      POTLJ  = POTLJ + EPQMMM(M,in) * (sigr12 - sigr6)

C30    CONTINUE

C     ENDIF

C-----density-site interaction-------------------------------------

C      DO 40 K=1,JZ
C        ADZ  = DZ*(K-NNMAXZ) 
C        ARZ  = ADZ - SCRD(3,in)
C        ARZ2 = ARZ*ARZ
C      DO 50 J=1,JY
C        ADY  = DY*(J-NNMAXY)
C        ARY  = ADY - SCRD(2,in)
C        ARY2 = ARY*ARY
C      DO 60 I=1,JX
C        ADX  = DX*(I-NNMAXX)
C        ARX  = ADX - SCRD(1,in)
C        ARX2 = ARX*ARX

C        AR   = DSQRT(ARX2+ARY2+ARZ2)         
C        AEF  = DERF(ALPHA*AR)/AR
C        VSSE = chgmm( in )*AEF    

C        CELS = CELS - RHO(I,J,K)*VSSE*DV

C60     CONTINUE
C50     CONTINUE
C40     CONTINUE

C15     CONTINUE

C      CALL MPI_REDUCE(CELS,SCELS,1,MPI_DOUBLE_PRECISION,
C    *                 MPI_SUM,0,MPI_COMM_WORLD,IERR)
C      CELS = SCELS
C   
C      CALL MPI_BCAST(CELS,1,MPI_DOUBLE_PRECISION,0,
C    *                MPI_COMM_WORLD,IERR)

C      ENDIF

C      IF(MYID.EQ.0) THEN

C      EQMMM = POTLJ + CNS + CELS
C      SSE   = SSE + EQMMM

C----- energy distribution function calculation------------------

C      IF(EQMMM.GT.EMINA) THEN
C      IF(EQMMM.LT.EMAXA) THEN

C        A=EQMMM-EMINA
C        NESS=1+JIDINT(A/DE)

C        NEDFMM(NESS) = NEDFMM(NESS) + 1
C        NERDF( NESS) = NERDF( NESS) + 1

C      ENDIF
C      ENDIF

C      ENDIF      ! close 'if(myid.eq.0) then'

C10    CONTINUE

C     IF(MYID.EQ.0) THEN

C     WRITE(64,*) SSE
C     WRITE(*,*) POTLJ
C     WRITE(*,*) CNS
C     WRITE(*,*) CELS

C-----energy distribution for the QM-electron interaction-----
C   
C     VQM = PENG-EPC-VPCZ-PLJ-EISO
C     IF(VQM.GT.EMINA) THEN
C     IF(VQM.LT.EMAXA) THEN
C       A = VQM-EMINA                       ! 2005.12.01 takahashi
C       NVQM = 1+JIDINT(A/DE)
C       NEDFQM(NVQM) = NEDFQM(NVQM) + 1
C       NERDF( NVQM) = NERDF( NVQM) + 1
C     ENDIF
C     ENDIF

C     WRITE(*,*) 'VQM =',VQM

C-----the average of distribution function--------------------

C     IF(MOD(IMD,100).EQ.0) THEN
C       WRITE(65,*) IMD
C       WRITE(66,*) IMD
C       WRITE(67,*) IMD
C       DENOM = 1.D0/DBLE(IMD-NCUT)
C       DO I = 1,NDXL
C         ETMP = NERDF( I)*DENOM
C         ETMP1= NEDFQM(I)*DENOM
C         ETMP2= NEDFMM(I)*DENOM
C         WRITE(65,*) ETMP
C         WRITE(66,*) ETMP1
C         WRITE(67,*) ETMP2
C       ENDDO
C     ENDIF

C-----CHI(I,J) calculation------------------

C      IF(MOD(IMD,10000).EQ.NCUT) THEN
C       WRITE(66,*) IMD 
C      ENDIF

C      REWIND(115)
C      REWIND(120)
C      REWIND(130)
C      DO I = 1,NDXL
C      DO J = 1,NDXL

C        CHI(I,J) = CHI(I,J) + NECHI(I)*NECHI(J)
C        BLCHI(I,J) = BLCHI(I,J) + NECHI(I)*NECHI(J)

C       IF(MOD(IMD,10000).EQ.NCUT) THEN
C        AVCHI = CHI(I,J)/DBLE(IMD-NCUT)
C        ACHI(I,J) = ACHI(I,J) + BLCHI(I,J)/1000 
C        AVBL = ACHI(I,J)/(NNST/1000)
C        BLCHI(I,J) = 0.D0
C	 WCHI = AVCHI-ETMP(I)*ETMP(J)
C        WRITE(66,*) AVCHI 
C        WRITE(120,*) AVBL
C        WRITE(130,*) WCHI
C       ENDIF

C      ENDDO
C      ENDDO

C     ENDIF

C     IF( NSTP .EQ. MMST ) THEN
C      NSTP=0
C     ENDIF

C     ENDIF       ! close 'if(myid.eq.0) then'

C     ENDIF

C     RETURN
C     END

C-------------------------------------------------
C
      SUBROUTINE EDF(RHO,SCRD,PENG,IMD)
C
C     calculation of energy distribution function
C     in the solution system  
C     for the use with Tinker package 
C
C-------------------------------------------------

      IMPLICIT REAL*8( A-H,O-Z ) 
      IMPLICIT INTEGER*4( I-N )

      include 'mpif.h'
      include 'mpi.i'

      INCLUDE 'QMpara.i'                        ! Vmol  
!     INCLUDE 'sizes.i'
!     INCLUDE 'atoms.i'
!     include 'atmtyp.i'

      PARAMETER ( NCUT  = 0  )               ! cut step
      PARAMETER ( ROMG  = 15.D0 )               ! solute-solvent interaction of finite range (Angstrom)

C-----/home3/takahasi/1b4v/FAD/mpi/fad/Vmol-tinker-----
      PARAMETER ( EISO  = -131.676068149D0 )    ! energy of FAD with N electrons at isolation

      PARAMETER ( NDXL  =  500    )             ! number of grids for the energy coordinate
      PARAMETER ( EMIN  =  -60.D0 )             ! unit of kcal/mol 
      PARAMETER ( EMAX  =   20.D0 )             ! unit of kcal/mol 
      
      PARAMETER ( ALPHA =   1.0D0  )      
      PARAMETER ( AUE   = 627.0898D0 )      
      PARAMETER ( AUL   = 0.529177D0 )      

      PARAMETER ( maxsol = 3000 )
      PARAMETER ( maxsit =   20 )

      COMMON / MMPI1  / MX,MY,MZ
      COMMON / MMPI4  / JX,JY,JZ

      COMMON / GRID2  / DX,DY,DZ
      COMMON / PRMT   / COE,TEMP,DTMD,
     *                  NDEN,MDMAX,
     *                  ZA(NNUC),ZVAL(NNUC),
     *                  SIG(NNUC),EPSQM(NNUC),
     *                  NRVLC,NCHK,CONV
      COMMON / PRMT2  / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                  NMMSW(maxatm),MMID(NNUC)
      COMMON / PRMT4  / nion1,iiont(maxatm)
      COMMON / NCLR   / PNUC( NNUC,NDIM )
      COMMON / SLT    / VPCE(NMAXX/NX,NMAXY/NY,NMAXZ/NZ),VPCZ,PLJ
      COMMON / CEDF   / DHOMO(NMAXX/NX,NMAXY/NY,NMAXZ/NZ)
      COMMON / LJQMMM / SQMMM(NNUC,maxatm),EPQMMM(NNUC,maxatm),
     *                  chgmm(maxatm),chgqm(NNUC)
      COMMON / AVRHO  / RHOAV(NMAXX/NX,NMAXY/NY,NMAXZ/NZ)

      COMMON / PRMTEDF / zmm(maxatm),nsol,nsite(maxsol),iiont1          ! 2005.08.31 takahashi
      COMMON / ZMASS   / ZCM(100)

      DIMENSION iiont1(maxsol,maxsit)

      DIMENSION RHO  ( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION SCRD ( NDIM,maxatm )
      DIMENSION CRD  ( NNUC,3 )
      DIMENSION NECHI( NDXL )
C     DIMENSION ETMP ( NDXL )                                           ! 2005.12.03 takahashi
      DIMENSION BLCHI( NDXL,NDXL )
      DIMENSION NEDFQM ( NDXL )
      DIMENSION NEDFMM ( NDXL )
      DIMENSION NERDF  ( NDXL )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      DE    = (EMAX-EMIN)/DBLE(NDXL)   ! kcal/mol
      DE    =  DE/AUE                  ! atomic unit
      EMINA =  EMIN/AUE                ! atomic unit
      EMAXA =  EMAX/AUE                ! atomic unit
      ROMGA =  ROMG/AUL                ! atomic unit
      ROMGA2 = ROMGA**2

      DV    =  DX*DY*DZ

C     IF(MYID.EQ.0) THEN
C       WRITE(*,*) 'iiont1(30,1)=',iiont1(30,1)
C     ENDIF

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

C-----compute center of mass of QM atoms------

      CX   = 0.D0
      CY   = 0.D0
      CZ   = 0.D0
      ZSUM = 0.D0

      DO NN = 1,NNUC
        NZA  = INT(ZA(NN))
        ZSUM = ZSUM + ZCM(NZA)
        CX   = CX + ZCM(NZA)*PNUC(NN,1)
        CY   = CY + ZCM(NZA)*PNUC(NN,2)
        CZ   = CZ + ZCM(NZA)*PNUC(NN,3)
      ENDDO
      CMASX = CX/ZSUM
      CMASY = CY/ZSUM
      CMASZ = CZ/ZSUM

      IF( IMD.EQ.0) THEN
      IF(MYID.EQ.0) THEN
       WRITE(*,*) 'C.M. = ',CMASX,CMASY,CMASZ
      ENDIF
      ENDIF

C-----initialize---------------

      IF( IMD .EQ. NCUT ) THEN
      
      DO I = 1,NDXL
       NEDFQM(I) = 0          ! energy distribution for QM energy
       NEDFMM(I) = 0          ! energy distribution for QMMM interaction
       NERDF( I) = 0          ! total energy distribution (QM + QM/MM)
      ENDDO

      ENDIF

C------------------------------

      IF( IMD .GT. NCUT ) THEN

      EPC  = 0.D0             ! for QM-energy coordinate
      DO K = 1,JZ
      DO J = 1,JY
      DO I = 1,JX
        EPC = EPC + VPCE(I,J,K)*RHO(I,J,K)*DV
        DHOMO(I,J,K) = RHO(I,J,K) - RHOAV(I,J,K)
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE(EPC,SEPC,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      EPC = SEPC
    
C     CALL MPI_BCAST(EPC,1,MPI_DOUBLE_PRECISION,0,
C    *               MPI_COMM_WORLD,IERR)


C-----solute-solvent interaction energy calculation----------------

      SSE = 0.D0
      DO 10 NSL  = 1,nsol

      SUMM = 0.D0 
      CX   = 0.D0
      CY   = 0.D0
      CZ   = 0.D0

      DO NS = 1,nsite(NSL)

        in   = iiont1(NSL,NS)   
        IF(NMMSW(in).EQ.1) GOTO 10 

        mn = in + 22                     ! 2005.11.26 takahashi
        SUMM = SUMM + zmm(mn)
        CX   = CX   + zmm(mn)*SCRD(1,in)
        CY   = CY   + zmm(mn)*SCRD(2,in)
        CZ   = CZ   + zmm(mn)*SCRD(3,in)

C     IF(MYID.EQ.3) THEN
C     write(*,*) zmm(mn),SCRD(1,in),SCRD(2,in),SCRD(3,in),chgmm(in)
C     ENDIF

      ENDDO

      CMMX = CX/SUMM
      CMMY = CY/SUMM
      CMMZ = CZ/SUMM

      RX = CMASX - CMMX
      RY = CMASY - CMMY
      RZ = CMASZ - CMMZ
      R2 = RX**2+RY**2+RZ**2

C-----solute-solvent interaction of finite range-------------------

      CELS  = 0.D0

      IF(R2 .LT. ROMGA2) THEN

      DO 15 NS = 1,nsite(NSL)

      in = iiont1(NSL,NS)    

C----- excess charge - site interaction ---------------------------

       DO 40 K=1,JZ
         ADZ  = DZ*(K-NNZ) 
         ARZ  = ADZ - SCRD(3,in)
         ARZ2 = ARZ*ARZ
       DO 50 J=1,JY
         ADY  = DY*(J-NNY)
         ARY  = ADY - SCRD(2,in)
         ARY2 = ARY*ARY
       DO 60 I=1,JX
         ADX  = DX*(I-NNX)
         ARX  = ADX - SCRD(1,in)
         ARX2 = ARX*ARX

         AR   = DSQRT(ARX2+ARY2+ARZ2)         
         AEF  = DERF(ALPHA*AR)/AR
         VSSE = chgmm( in )*AEF    

         CELS = CELS - DHOMO(I,J,K)*VSSE*DV

60     CONTINUE
50     CONTINUE
40     CONTINUE

15     CONTINUE

       CALL MPI_REDUCE(CELS,SCELS,1,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR)
       CELS = SCELS

       IF(MYID.EQ.0) THEN

       EQMMM = CELS
       SSE   = SSE + EQMMM

C----- energy distribution function calculation------------------

       IF(EQMMM.GT.EMINA) THEN
       IF(EQMMM.LT.EMAXA) THEN

         A=EQMMM-EMINA
         NESS=1+ IDINT(A/DE)

         NEDFMM(NESS) = NEDFMM(NESS) + 1
         NERDF( NESS) = NERDF( NESS) + 1

       ENDIF
       ENDIF

       ENDIF      ! close 'if(myid.eq.0) then'

       ENDIF

10    CONTINUE

      IF(MYID.EQ.0) THEN

      WRITE(64,*) SSE
C     WRITE(*,*) POTLJ
C     WRITE(*,*) CNS
C     WRITE(*,*) CELS

C-----energy distribution for the QM-electron interaction-----
    
      VQM = PENG-EPC-VPCZ-PLJ-EISO
      IF(VQM.GT.EMINA) THEN
      IF(VQM.LT.EMAXA) THEN
        A = VQM-EMINA                       ! 2005.12.01 takahashi
        NVQM = 1+ IDINT(A/DE)
        NEDFQM(NVQM) = NEDFQM(NVQM) + 1
        NERDF( NVQM) = NERDF( NVQM) + 1
      ENDIF
      ENDIF

      WRITE(*,*) 'VQM =',VQM

C-----the average of distribution function--------------------

      IF(MOD(IMD,100).EQ.0) THEN
        WRITE(65,*) IMD
        WRITE(66,*) IMD
        WRITE(67,*) IMD
        DENOM = 1.D0/DBLE(IMD-NCUT)
        DO I = 1,NDXL
          ETMP = NERDF( I)*DENOM
          ETMP1= NEDFQM(I)*DENOM
          ETMP2= NEDFMM(I)*DENOM
          WRITE(65,*) ETMP
          WRITE(66,*) ETMP1
          WRITE(67,*) ETMP2
        ENDDO
      ENDIF

C-----CHI(I,J) calculation------------------

C      IF(MOD(IMD,10000).EQ.NCUT) THEN
C       WRITE(66,*) IMD 
C      ENDIF

C      REWIND(115)
C      REWIND(120)
C      REWIND(130)
C      DO I = 1,NDXL
C      DO J = 1,NDXL

C        CHI(I,J) = CHI(I,J) + NECHI(I)*NECHI(J)
C        BLCHI(I,J) = BLCHI(I,J) + NECHI(I)*NECHI(J)

C       IF(MOD(IMD,10000).EQ.NCUT) THEN
C        AVCHI = CHI(I,J)/DBLE(IMD-NCUT)
C        ACHI(I,J) = ACHI(I,J) + BLCHI(I,J)/1000 
C        AVBL = ACHI(I,J)/(NNST/1000)
C        BLCHI(I,J) = 0.D0
C	 WCHI = AVCHI-ETMP(I)*ETMP(J)
C        WRITE(66,*) AVCHI 
C        WRITE(120,*) AVBL
C        WRITE(130,*) WCHI
C       ENDIF

C      ENDDO
C      ENDDO

C     ENDIF

C     IF( NSTP .EQ. MMST ) THEN
C      NSTP=0
C     ENDIF

      ENDIF       ! close 'if(myid.eq.0) then'

      ENDIF

      RETURN
      END

C----------------------------------
C   SUBROUTINE COMPUTING DENSITY 
C----------------------------------

      SUBROUTINE DNSHOMO( RWFAB ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ
      COMMON / NCLR  / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / CEDF  / DHOMO(NMAXX/NX,NMAXY/NY,NMAXZ/NZ)

      DIMENSION RWFAB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )


      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      DV = DX*DY*DZ


C----- density of excess electron -----

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

C        DHOMO( I,J,K ) = RWFAB(I,J,K,MORA)**2

      ENDDO
      ENDDO
      ENDDO

C----- normalize -----

      ANORM = 0.D0

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

C         ANORM = ANORM + DHOMO( I,J,K )*DV

      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE(ANORM,BNORM,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST( BNORM,1,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

C     IF( myid .eq. 0 ) then
C      write(*,*) BNORM
C     endif

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

C         DHOMO( I,J,K ) = DHOMO( I,J,K ) / BNORM  

      ENDDO
      ENDDO
      ENDDO

      RETURN 
      END

C-------------------------------------------
C     SUBROUTINE POINT CHARGE 
C-------------------------------------------

      SUBROUTINE POCH1( SCRD )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
!     include 'atoms.i'
C     include 'charge.i'

C     PARAMETER ( NMAX =  80 )
C     PARAMETER ( NDIM =   3 )
C     PARAMETER ( NNUC =  36 )
      PARAMETER ( ALPHA = 1.0D0   )
C     PARAMETER ( NLINK1 = 10 )
      PARAMETER ( RRCUT = 4.724D0 )
      PARAMETER ( RRCUT1= 5.669D0 )
      PARAMETER ( DRCUT = 0.945D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PRMT4 / nion1,iiont(maxatm)
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / SLT / VPCE(NMAXX/NX,NMAXY/NY,NMAXZ/NZ),VPCZ,PLJ
      COMMON / LJQMMM / SQMMM(NNUC,maxatm),EPQMMM(NNUC,maxatm),
     *                  chgmm(maxatm),chgqm(NNUC)

      DIMENSION SCRD( NDIM,maxatm )
      DIMENSION ADX(  NMAXX/NX )
      DIMENSION ADY(  NMAXY/NY )
      DIMENSION ADZ(  NMAXZ/NZ )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

C
C     number of loop in QM and MM subsystems
C

      IF(MYID.EQ.0) THEN
C     OPEN(111,FILE='CHG.DAT',STATUS='UNKNOWN')            ! open output files 
        WRITE(*,*) 'Start poch1'
      qsum = 0.d0
      DO iin = 1, nion1
        in   = iiont(iin)
C       write(111,*) chgmm(in)
        qsum = qsum + chgmm(in)
      ENDDO
      write(*,*) 'qsum = ', qsum
      ENDIF

      TIME1 = MPI_WTIME()

      nlmm = n-NNUC+NLINK
      nlqm = NNUC-NMM2-NLINK

      ONE    = 1.D0
      RRCUT2 = ONE/(RRCUT*RRCUT)
C     NNMAX  = NMAX/2 + 1

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO K=1,JZ
        ADZ(K) = DZ*( K-NNZ )
      ENDDO
      DO J=1,JY
        ADY(J) = DY*( J-NNY )
      ENDDO
      DO I=1,JX
        ADX(I) = DX*( I-NNX )
      ENDDO

      DO K=1,JZ
      DO J=1,JY
      DO I=1,JX
        VPCE( I,J,K )=0.0D0            ! initialize
      ENDDO
      ENDDO
      ENDDO

      DO 70 iin=1,nion1
        in = iiont(iin)

       IF(NMMSW(in).EQ.1) GOTO 70

       DO    K=1,JZ
        RZ    = ADZ(K) - SCRD( 3,in )
       DO    J=1,JY
        RY    = ADY(J) - SCRD( 2,in )
       DO 75 I=1,JX
        RX    = ADX(I) - SCRD( 1,in )

        R2    = RX**2+RY**2+RZ**2
        R     = DSQRT( R2 )

        IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
          ESW  = ONE
        ELSEIF(R.GT.DRCUT) THEN
          ARD  = R - DRCUT
          ARD2 = ARD*ARD*RRCUT2
          ESW1 = ONE - ARD2
          ESW2 = ESW1*ESW1
          ESW  = ONE - ESW2
        ELSE
          GOTO 75
        ENDIF

        VPCE( I,J,K )= VPCE( I,J,K )
     *               - ESW*chgmm( in )*(ERF(ALPHA*R)/R)

75     CONTINUE
       ENDDO
       ENDDO

70    CONTINUE

      TIME2 = MPI_WTIME()

      IF(MYID.EQ.0) THEN
        WRITE(*,99) 'Elapsed Time (vpce) = ',TIME2-TIME1,' scnds'
      ENDIF

      VPCZ=0.D0
      PLJ =0.D0

      IF(MYID.EQ.0) THEN

      DO 80 in=1,nlmm
      DO 90 M=1,NNUC 

        RX=PNUC( M,1 ) - SCRD( 1,in )
        RY=PNUC( M,2 ) - SCRD( 2,in )
        RZ=PNUC( M,3 ) - SCRD( 3,in )
        R =DSQRT(RX**2+RY**2+RZ**2)
        RI = 1.D0/R

        sigr6  = (SQMMM(M,in)*RI*RI)**3
        sigr12 = sigr6*sigr6

        ETMP = EPQMMM(M,in) * ( sigr12 - sigr6 )
        IF( ETMP .GT. 1.0D-2 ) THEN
          WRITE(*,*) 'PLJ in poch1.f: ', ETMP, in 
        ENDIF
        PLJ  = PLJ + EPQMMM(M,in) * ( sigr12 - sigr6 )    ! LJ

90    CONTINUE
80    CONTINUE
      
      TIME3 = MPI_WTIME()
      WRITE(*,99) 'Elapsed Time (plj)  = ',TIME3-TIME2,' scnds'

      DO 85 iin=1,nion1
        in = iiont(iin)
      DO 95 M=1,NNUC 

        RX=PNUC( M,1 ) - SCRD( 1,in )
        RY=PNUC( M,2 ) - SCRD( 2,in )
        RZ=PNUC( M,3 ) - SCRD( 3,in )
        R =DSQRT(RX**2+RY**2+RZ**2)
        RI = 1.D0/R

        IF(NMMSW(in).EQ.1) GOTO 95

        IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
          ESW = 1.D0
        ELSEIF(R.GT.DRCUT) THEN
          ARD  = R - DRCUT
          ARD2 = ARD*ARD*RRCUT2
          ESW1 = 1.D0 - ARD2
          ESW2 = ESW1*ESW1
          ESW  = 1.D0 - ESW2
        ELSE
          GOTO 95
        ENDIF

        VPCZ = VPCZ + ZVAL( M )*ESW*chgmm( in )*RI   ! Coulomb (nuclear-site)

95    CONTINUE
85    CONTINUE

      TIME4 = MPI_WTIME()
      WRITE(*,99) 'Elapsed Time (vpcz) = ',TIME4-TIME3,' scnds'
      
      WRITE(*,*)
      WRITE(*,*)'Nuclear - Site coulomb energy'
      WRITE(*,*)'  VPCZ = ' ,VPCZ
      WRITE(*,*)
      WRITE(*,*)'Nuclear - Site LJ potential energy'
      WRITE(*,*)'  PLJ  = ' ,PLJ

      ENDIF

C     WRITE(*,*) 'rank',myid,'is alive.'

      CALL MPI_BCAST( VPCZ,1,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )
      CALL MPI_BCAST( PLJ,1,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )

99    FORMAT( X,A22,F10.6,3X,A6 )

      RETURN
      END

C-------------------------------------------
C     SUBROUTINE QM/MM FORCE
C-------------------------------------------

      SUBROUTINE QMF1( SCRD,FRCS,FRCN,RHO )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
!     include 'atoms.i'

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC =  36 )
      PARAMETER ( DCUT  = 1.D-9 ) 
      PARAMETER ( ALPHA = 1.0D0 )
C     PARAMETER ( NLINK1 = 10 )
      PARAMETER ( RRCUT = 4.724D0 )
      PARAMETER ( RRCUT1= 5.669D0 )
      PARAMETER ( DRCUT = 0.945D0 )

      PARAMETER ( PI    = 3.14159265358979323D0 )
      
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PRMT4 / nion1,iiont(maxatm)
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / LJQMMM / SQMMM(NNUC,maxatm),EPQMMM(NNUC,maxatm),
     *                  chgmm(maxatm),chgqm(NNUC)

      DIMENSION SCRD( NDIM,maxatm )

      DIMENSION FDSS(  maxatm,NDIM )              ! density-site force on site
      DIMENSION FODSS( maxatm,NDIM )              ! density-site force on site
      DIMENSION FNSS(  maxatm,NDIM )              ! nuclear-site and LJ force on site
      DIMENSION FRCS(  maxatm,NDIM )              ! QM/MM force on site

      DIMENSION FRCN( NNUC,NDIM )                 ! QM/MM force on nuclear

      DIMENSION RHO( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      DIMENSION ADX( NMAXX/NX )
      DIMENSION ADY( NMAXY/NY )
      DIMENSION ADZ( NMAXZ/NZ )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      IF(MYID.EQ.0) THEN
        WRITE(*,*) 'Start qmf1'
      ENDIF

      DV=DX*DY*DZ

C
C     number of loop in QM and MM subsystems
C
      nlmm = n-NNUC+NLINK

      DO L=1,3
      DO in=1,n

        FDSS(  in,L ) = 0.D0
        FODSS( in,L ) = 0.D0
        FNSS(  in,L ) = 0.D0
        FRCS(  in,L ) = 0.D0
            
      ENDDO
      ENDDO
     
      DO L=1,3
      DO in=1,NNUC
     
        FRCN( in,L )   = 0.D0
            
      ENDDO
      ENDDO

C-----density-site force CALCULATION
     
      TIME1 = MPI_WTIME()

      ONE    = 1.D0
      TWO    = 2.D0
      FOUR   = 4.D0
      RRCUT2 = ONE/(RRCUT*RRCUT)
C     NNMAX  = NMAX/2 + 1
      SQPII  = ONE/DSQRT(PI)
      A2SQPI = TWO*ALPHA*SQPII

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO K=1,JZ
        ADZ(K) = DZ*( K-NNZ )
      ENDDO
      DO J=1,JY
        ADY(J) = DY*( J-NNY )
      ENDDO
      DO I=1,JX
        ADX(I) = DX*( I-NNX )
      ENDDO

      DO 50 iin=1,nion1
       in = iiont(iin)

      IF(NMMSW(in).EQ.1) THEN
        IF(MYID.EQ.0) WRITE(*,*) 'enter 1',in
        GOTO 50
      ENDIF

      DO 60 K=1,JZ
       RZ = ADZ(K) - SCRD( 3,in )
      DO 70 J=1,JY
       RY = ADY(J) - SCRD( 2,in )
      DO 80 I=1,JX

       ARHO = RHO(I,J,K)
       IF( ARHO.GT.DCUT ) THEN

       RX = ADX(I) - SCRD( 1,in )

       R2 = RX**2+RY**2+RZ**2
       R  = DSQRT( R2 )

       IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
         RI    = ONE/R
         R2I   = RI*RI
         R3I   = RI*R2I
         APR   = ALPHA*R
         C2RDV = ARHO * chgmm( in ) * DV
         TFDSS = C2RDV*(A2SQPI*DEXP(-APR**2)*R2I-ERF(APR)*R3I)
       ELSEIF(R.GT.DRCUT) THEN
         ARD   = R - DRCUT
         ARD2  = ARD*ARD*RRCUT2
         ESW1  = ONE - ARD2
         ESW2  = ESW1*ESW1
         ESW   = ONE - ESW2
         ESWA  = FOUR*ESW1*ARD*RRCUT2
         RI    = ONE/R
         R2I   = RI*RI
         APR   = ALPHA*R
         C2RDV = ARHO * chgmm( in ) * DV
         DER2I = ERF(APR)*R2I
         TFDSS = ESW*C2RDV*(A2SQPI*DEXP(-APR**2)*R2I-DER2I*RI)
     *         + ESWA*C2RDV*DER2I
       ELSE
         IF(MYID.EQ.0) WRITE(*,*) 'enter 2'
         GOTO 50
       ENDIF

       FODSS(in,1) = FODSS(in,1)-RX*TFDSS
       FODSS(in,2) = FODSS(in,2)-RY*TFDSS
       FODSS(in,3) = FODSS(in,3)-RZ*TFDSS

       ENDIF

80    CONTINUE
70    CONTINUE
60    CONTINUE

50    CONTINUE

C     NUMDAT = n*NDIM
      NUMDAT = maxatm*NDIM
      CALL MPI_REDUCE( FODSS(1,1),FDSS(1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )
C     WRITE(*,*) 'MYID:qmf1=',MYID,RHO(5,5,5)

      TIME2 = MPI_WTIME()

      IF(MYID.EQ.0) THEN
        WRITE(*,99) 'Elapsed Time (qmf1) = ',TIME2-TIME1,' scnds'
      ENDIF

C-----nuclear-site force CALCULATION
C-----LJ force CALCULATION

      IF(MYID.EQ.0) THEN

      SIX   = 6.D0

      DO 90 in=1,nlmm
      DO 100 M=1,NNUC 

        RX=PNUC( M,1 ) - SCRD( 1,in )
        RY=PNUC( M,2 ) - SCRD( 2,in )
        RZ=PNUC( M,3 ) - SCRD( 3,in )
        R2=RX**2+RY**2+RZ**2
        R   = DSQRT(R2)
        RI  = ONE/R
        R2I = RI*RI

        ASIG2  = SQMMM(M,in)*R2I
        ASIG6  = ASIG2*ASIG2*ASIG2
        ASIG12 = ASIG6*ASIG6
        BSIG   = TWO*ASIG12 - ASIG6
        TFEPS  = SIX*EPQMMM(M,in)
        CREPS  = TFEPS*R2I*BSIG
        CREPSX = CREPS*RX
        CREPSY = CREPS*RY
        CREPSZ = CREPS*RZ

        FNSS(in,1) = FNSS(in,1) + CREPSX   ! LJ force on site
        FNSS(in,2) = FNSS(in,2) + CREPSY
        FNSS(in,3) = FNSS(in,3) + CREPSZ

        FRCN(M,1)  = FRCN(M,1) + CREPSX    ! LJ force on nuclear
        FRCN(M,2)  = FRCN(M,2) + CREPSY
        FRCN(M,3)  = FRCN(M,3) + CREPSZ

100   CONTINUE
90    CONTINUE
      
      TIME3 = MPI_WTIME()
        WRITE(*,99) 'Elapsed Time (qmf2) = ',TIME3-TIME2,' scnds'

      DO 95 iin=1,nion1
        in = iiont(iin)
      DO 105 M=1,NNUC 

        RX=PNUC( M,1 ) - SCRD( 1,in )
        RY=PNUC( M,2 ) - SCRD( 2,in )
        RZ=PNUC( M,3 ) - SCRD( 3,in )
        R2=RX**2+RY**2+RZ**2
        R   = DSQRT(R2)
        RI  = ONE/R
        R2I = RI*RI

       IF(NMMSW(in).EQ.1) GOTO 105

       IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
         R3I   = RI*R2I
         ACHR  = ZVAL( M )*chgmm( in )*R3I
       ELSEIF(R.GT.DRCUT) THEN
         ARD   = R - DRCUT
         ARD2  = ARD*ARD*RRCUT2
         ESW1  = ONE - ARD2
         ESW2  = ESW1*ESW1
         ESW   = ONE - ESW2
         ESWA  = FOUR*ESW1*ARD*RRCUT2
         ZVPCH = ZVAL( M )*chgmm( in )*R2I
         ACHR  = ESW*ZVPCH*RI - ESWA*ZVPCH
       ELSE
         GOTO 105
       ENDIF

        CHRX = ACHR*RX
        CHRY = ACHR*RY
        CHRZ = ACHR*RZ

        FNSS(in,1) = FNSS(in,1) + CHRX      ! Coulomb force on site
        FNSS(in,2) = FNSS(in,2) + CHRY
        FNSS(in,3) = FNSS(in,3) + CHRZ

        FRCN(M,1)  = FRCN(M,1) + CHRX       ! Coulomb force on nuclear
        FRCN(M,2)  = FRCN(M,2) + CHRY
        FRCN(M,3)  = FRCN(M,3) + CHRZ

105   CONTINUE
95    CONTINUE

      TIME4 = MPI_WTIME()
        WRITE(*,99) 'Elapsed Time (qmf3) = ',TIME4-TIME3,' scnds'

      DO L=1,3
      DO in=1,nlmm
        FRCS(in,L) = FDSS(in,L) - FNSS(in,L)                 !QM/MM force on site
      ENDDO
      ENDDO

      ENDIF

99    FORMAT( X,A22,F10.6,3X,A6 )

      RETURN
      END

C-----------------------------------------------
C     SUBROUTINE remove artificial electrostatic
C     energy and force between artificial bond
C     including link atom and MM subsystem
C-----------------------------------------------

      SUBROUTINE EXTCORR( FRCS,FRC,SCRD )

      implicit real*8 ( a-h,o-z )
      implicit integer*4 ( i-n )

      include "mpi.i"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
!     include 'atoms.i'
!     include 'charge.i'

C     PARAMETER ( NMAX  = 80 )
C     PARAMETER ( NNUC  = 36 )
C     PARAMETER ( NDIM  = 3 )
C     PARAMETER ( NLINK1 = 10 )
      PARAMETER ( RRCUT = 4.724D0 )
      PARAMETER ( RRCUT1= 5.669D0 )
      PARAMETER ( DRCUT = 0.945D0 )

      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / OFC1 / A( NNUC )
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
C     COMMON / SLT / VPCE(NMAX,NMAX,NMAX),VPCZ,PLJ
      COMMON / SLT / VPCE(NMAXX/NX,NMAXY/NY,NMAXZ/NZ),VPCZ,PLJ
      COMMON / LJQMMM / SQMMM(NNUC,maxatm),EPQMMM(NNUC,maxatm),
     *                  chgmm(maxatm),chgqm(NNUC)

      DIMENSION FRCS( maxatm,NDIM )
      DIMENSION FRC( NNUC,NDIM )
      DIMENSION SCRD( NDIM,maxatm )

C
C     number of loop in QM and MM subsystems
C
      nlmm   = n-NMM2
      ncount = NNUC - NLINK

      VPCZC  = 0.D0
      ONE    = 1.D0
      FOUR   = 4.D0
      RRCUT2 = ONE/(RRCUT*RRCUT)

      DO 80 in=1,nlmm
        IF(NMMSW(in).EQ.1) GOTO 80
      DO 90 M=1,NLINK 

        ncount = ncount + 1

        RX=PNUC( NLAQM(M),1 ) - SCRD( 1,in )
        RY=PNUC( NLAQM(M),2 ) - SCRD( 2,in )
        RZ=PNUC( NLAQM(M),3 ) - SCRD( 3,in )
        R =DSQRT(RX**2+RY**2+RZ**2)
        RI = ONE/R
        R2I = RI*RI

        IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
          ESW   = ONE
          R3I   = RI*R2I
          ACHR  = A( ncount )*chgmm( in )*R3I
        ELSEIF(R.GT.DRCUT) THEN
          ARD   = R - DRCUT
          ARD2  = ARD*ARD*RRCUT2
          ESW1  = ONE - ARD2
          ESW2  = ESW1*ESW1
          ESW   = ONE - ESW2
          ESWA  = FOUR*ESW1*ARD*RRCUT2
          ZVPCH = A( ncount )*chgmm( in )*R2I
          ACHR  = ESW*ZVPCH*RI - ESWA*ZVPCH
        ELSE
          GOTO 95
        ENDIF

        VPCZC = VPCZC + A( ncount )*ESW*chgmm( in )*RI   ! Coulomb (nuclear-site)

        CHRX = ACHR*RX
        CHRY = ACHR*RY
        CHRZ = ACHR*RZ

        FRCS(in,1) = FRCS(in,1) - CHRX      ! Coulomb force on site
        FRCS(in,2) = FRCS(in,2) - CHRY
        FRCS(in,3) = FRCS(in,3) - CHRZ

        FRC(NLAQM(M),1)  = FRC(NLAQM(M),1) + CHRX       ! Coulomb force on nuclear
        FRC(NLAQM(M),2)  = FRC(NLAQM(M),2) + CHRY
        FRC(NLAQM(M),3)  = FRC(NLAQM(M),3) + CHRZ

95      CONTINUE

        RX=PNUC( ncount,1 ) - SCRD( 1,in )
        RY=PNUC( ncount,2 ) - SCRD( 2,in )
        RZ=PNUC( ncount,3 ) - SCRD( 3,in )
        R =DSQRT(RX**2+RY**2+RZ**2)
        RI = ONE/R
        R2I = RI*RI

        IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
          ESW   = ONE
          R3I   = RI*R2I
          ACHR  = A( ncount )*chgmm( in )*R3I
        ELSEIF(R.GT.DRCUT) THEN
          ARD   = R - DRCUT
          ARD2  = ARD*ARD*RRCUT2
          ESW1  = ONE - ARD2
          ESW2  = ESW1*ESW1
          ESW   = ONE - ESW2
          ESWA  = FOUR*ESW1*ARD*RRCUT2
          ZVPCH = A( ncount )*chgmm( in )*R2I
          ACHR  = ESW*ZVPCH*RI - ESWA*ZVPCH
        ELSE
          GOTO 90
        ENDIF

        VPCZC = VPCZC - A( ncount )*ESW*chgmm( in )*RI   ! Coulomb (nuclear-site)

        CHRX = ACHR*RX
        CHRY = ACHR*RY
        CHRZ = ACHR*RZ

        FRCS(in,1) = FRCS(in,1) + CHRX              ! Coulomb force on site
        FRCS(in,2) = FRCS(in,2) + CHRY
        FRCS(in,3) = FRCS(in,3) + CHRZ

        FRC(ncount,1)  = FRC(ncount,1) - CHRX       ! Coulomb force on nuclear
        FRC(ncount,2)  = FRC(ncount,2) - CHRY
        FRC(ncount,3)  = FRC(ncount,3) - CHRZ

90    CONTINUE

80    CONTINUE
      
      WRITE(*,*)
      WRITE(*,*)'Artificial QM - MM coulomb energy'
      WRITE(*,*)'  VPCZC = ' ,VPCZC

      VPCZ = VPCZ + VPCZC

      RETURN
      END

C-----------------------------------------------
C     SUBROUTINE split QM/MM boundary forse
C     for covalent bonds
C-----------------------------------------------

      SUBROUTINE SPLTFRC( FRCS,FRC,SCRD )

      implicit real*8 ( a-h,o-z )
      implicit integer*4 ( i-n )

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
!     include 'atoms.i'

C     PARAMETER ( NNUC  = 36 )
C     PARAMETER ( NDIM  = 3 )
C     PARAMETER ( NLINK1 = 10 )

      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)

      DIMENSION FRCS( maxatm,NDIM )
      DIMENSION FRC( NNUC,NDIM )
      DIMENSION SCRD( NDIM,maxatm )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      IF( MYID .EQ. 0 ) THEN

      ns = NNUC - NLINK + 1
      nc = 0

      do i = ns, NNUC

        nc = nc + 1

        xcc = scrd(1,NLAMM(nc)) - PNUC(NLAQM(nc),1)
        ycc = scrd(2,NLAMM(nc)) - PNUC(NLAQM(nc),2)
        zcc = scrd(3,NLAMM(nc)) - PNUC(NLAQM(nc),3)
        rcc = xcc*xcc + ycc*ycc + zcc*zcc
        rcc = dsqrt(rcc)
        xcc = xcc/rcc
        ycc = ycc/rcc
        zcc = zcc/rcc

        fba = FRC(i,1)*xcc + FRC(i,2)*ycc + FRC(i,3)*zcc
        fbx = fba*xcc
        fby = fba*ycc
        fbz = fba*zcc
        fpx = FRC(i,1) - fbx
        fpy = FRC(i,2) - fby
        fpz = FRC(i,3) - fbz

        xch = PNUC(NLAQM(nc),1) - PNUC(i,1)
        ych = PNUC(NLAQM(nc),2) - PNUC(i,2)
        zch = PNUC(NLAQM(nc),3) - PNUC(i,3)
        rch = xch*xch + ych*ych + zch*zch
        rch = dsqrt(rch)

        coeff = rch/rcc
        fpmx  = coeff * fpx
        fpmy  = coeff * fpy
        fpmz  = coeff * fpz

        FRCS( NLAMM(nc),1 ) = FRCS( NLAMM(nc),1 ) + fpmx + fbx
        FRCS( NLAMM(nc),2 ) = FRCS( NLAMM(nc),2 ) + fpmy + fby
        FRCS( NLAMM(nc),3 ) = FRCS( NLAMM(nc),3 ) + fpmz + fbz

        fpqx = fpx - fpmx
        fpqy = fpy - fpmy
        fpqz = fpz - fpmz

        FRC( NLAQM(nc),1 ) = FRC( NLAQM(nc),1 ) + fpqx
        FRC( NLAQM(nc),2 ) = FRC( NLAQM(nc),2 ) + fpqy
        FRC( NLAQM(nc),3 ) = FRC( NLAQM(nc),3 ) + fpqz

      enddo

      ENDIF

C     nlmm = n-NMM2
C     nlqm = NNUC-NMM2-NLINK

C     DO L = 1, 3
C       NC = nlqm
C     DO in = nlmm+1, n
C       NC = NC + 1
C       FRCS(in,L) = FRC(NC,L)
C     ENDDO
C     ENDDO

      RETURN
      END

C ----------------------------------------------------------------
C     RADIAL DITRIBUTION FUNCTION 
C ---------------------------------------------------------------

      SUBROUTINE RDFL(CRD,SCRD,IMD)

      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT INTEGER*4 (I-N) 

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
!     include 'atoms.i'
!     include 'boxes.i'

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NNUC =  36 )
      PARAMETER ( NDX = 100 )
      PARAMETER ( AUL  = 0.5291771D0 )     
      PARAMETER ( PI   = 3.14159265358979323D0 )

      DIMENSION CRD (NNUC,3)
      DIMENSION SCRD( 3,maxatm )

      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / NRDF1 / NMRDF,NRDF(NNUC)
      COMMON / NNRDF / NGO(NNUC,NDX),NGH(NNUC,NDX)
     
      DX = 0.188972687D0
      TBOXL = xbox/AUL
      TVOL  = TBOXL**3
    
      IF (IMD.EQ.1.OR.IMD.EQ.1000) THEN
      NMD = 0
      DO I = 1 ,NNUC
      DO J = 1 ,NDX 
        NGO(I,J) = 0
        NGH(I,J) = 0
      ENDDO
      ENDDO
      ENDIF

      N3I = n/3

      DO J = 1 ,NMRDF
 
      K=NRDF(J)
 
      DO I = 1 ,N3I

      RO = (CRD(K,1)-SCRD(1,3*I-2))**2
     *   + (CRD(K,2)-SCRD(2,3*I-2))**2
     *   + (CRD(K,3)-SCRD(3,3*I-2))**2
     
      SRO = DSQRT(RO) / DX
      NRO = INT(SRO)

      IF(NRO.LE.NDX) THEN

      NGO(J,NRO) = NGO(J,NRO) + 1

      ENDIF

      RH = (CRD(K,1)-SCRD(1,3*I-1))**2
     *   + (CRD(K,2)-SCRD(2,3*I-1))**2
     *   + (CRD(K,3)-SCRD(3,3*I-1))**2
     
      SRH = DSQRT(RH) / DX
      NRH = INT(SRH)

      IF(NRH.LE.NDX) THEN

      NGH(J,NRH) = NGH(J,NRH) + 1

      ENDIF

      RH = (CRD(K,1)-SCRD(1,3*I))**2
     *   + (CRD(K,2)-SCRD(2,3*I))**2
     *   + (CRD(K,3)-SCRD(3,3*I))**2
     
      SRH = DSQRT(RH) / DX
      NRH = INT(SRH)

      IF(NRH.LE.NDX) THEN

      NGH(J,NRH) = NGH(J,NRH) + 1

      ENDIF

      ENDDO  
      
      ENDDO  

      DDX  = DX*0.529177D0
      SDDX = 0.5D0*DDX

      NMD = NMD + 1

      IF (MOD(IMD,1000).EQ.0.AND.IMD.GT.1000) THEN

      DO J = 1,NMRDF
      IN = J+100
      WRITE(IN,*) IMD 
      DO I = 1,NDX
      R  = DX*DBLE(I)
      R2 = R**2
      SDV = 4.D0*PI*R2*DX+4.D0*PI*R*DX*DX+4.D0/3.D0*PI*DX**3
      OTMP = NGO(J,I)*TVOL/(      N3I*NMD*SDV )    
      HTMP = NGH(J,I)*TVOL/( 2.D0*N3I*NMD*SDV )
      RDX = DDX * I + SDDX 
      WRITE(IN,*) RDX,OTMP,HTMP
      ENDDO
      ENDDO

      ENDIF
      
      RETURN
      END

C--------------------------------------------------------
C     SUBROUTINE POINT CHARGE for Molecular Mechanics
C--------------------------------------------------------

      SUBROUTINE POCH2( SCRD )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpi.i"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
!     include 'atoms.i'
C     include 'charge.i'

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NLINK1 = 10 )
      PARAMETER ( RRCUT = 4.724D0 )
      PARAMETER ( RRCUT1= 5.669D0 )
      PARAMETER ( DRCUT = 0.945D0 )

      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PRMT4 / nion1,iiont(maxatm)
      COMMON / NCLR / PNUC( NNUC,NDIM )
C     COMMON / SLT / VPCE(NMAX,NMAX,NMAX),VPCZ,PLJ
      COMMON / SLT / VPCE(NMAXX/NX,NMAXY/NY,NMAXZ/NZ),VPCZ,PLJ
      COMMON / LJQMMM / SQMMM(NNUC,maxatm),EPQMMM(NNUC,maxatm),
     *                  chgmm(maxatm),chgqm(NNUC)

      DIMENSION SCRD( NDIM,maxatm )

C
C     number of loop in QM and MM subsystems
C
      nlqm = NNUC-NLINK
      nlqm1 = nlqm+1
      nlmm = n-NNUC+NLINK

      ONE    = 1.D0
      RRCUT2 = ONE/(RRCUT*RRCUT)

      VPCZ=0.D0
      TVPCZ=0.D0
      PLJ=0.D0

      DO 80 in=1,nlmm
      DO 90 M=1,nlqm 

        RX=PNUC( M,1 ) - SCRD( 1,in )
        RY=PNUC( M,2 ) - SCRD( 2,in )
        RZ=PNUC( M,3 ) - SCRD( 3,in )
        R =DSQRT(RX**2+RY**2+RZ**2)
        RI = 1.D0/R

        sigr6  = (SQMMM(M,in)*RI*RI)**3
        sigr12 = sigr6*sigr6

        PLJ  = PLJ + EPQMMM(M,in) * ( sigr12 - sigr6 )    ! LJ

        IF(NMMSW(in).EQ.1) GOTO 90

        IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
          ESW = ONE
        ELSEIF(R.GT.DRCUT) THEN
          ARD  = R - DRCUT
          ARD2 = ARD*ARD*RRCUT2
          ESW1 = ONE - ARD2
          ESW2 = ESW1*ESW1
          ESW  = ONE - ESW2
        ELSE
          GOTO 90
        ENDIF
        VPCZ = VPCZ + chgqm( M )*ESW*chgmm( in )*RI   ! Coulomb (nuclear-site)
        TVPCZ = TVPCZ + ZVAL( M )*ESW*chgmm( in )*RI   ! Coulomb (nuclear-site)

90    CONTINUE
      DO 91 M=nlqm1,NNUC

        RX=PNUC( M,1 ) - SCRD( 1,in )
        RY=PNUC( M,2 ) - SCRD( 2,in )
        RZ=PNUC( M,3 ) - SCRD( 3,in )
        R =DSQRT(RX**2+RY**2+RZ**2)
        RI = 1.D0/R

        IF(NMMSW(in).EQ.1) GOTO 91

        IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
          ESW = ONE
        ELSEIF(R.GT.DRCUT) THEN
          ARD  = R - DRCUT
          ARD2 = ARD*ARD*RRCUT2
          ESW1 = ONE - ARD2
          ESW2 = ESW1*ESW1
          ESW  = ONE - ESW2
        ELSE
          GOTO 91
        ENDIF
        TVPCZ = TVPCZ + ZVAL( M )*ESW*chgmm( in )*RI   ! Coulomb (nuclear-site)

91    CONTINUE
80    CONTINUE

      ETOT = VPCZ + PLJ
      
      WRITE(*,*)
      WRITE(*,*)'  Eelec = ',VPCZ
      WRITE(*,*)'  Elj   = ',PLJ
      WRITE(*,*)'  Etot  = ',ETOT
      WRITE(*,*)'  VPCZ  = ',TVPCZ

      RETURN
      END

C--------------------------------------------------------
C     SUBROUTINE QM/MM FORCE for Molecular Mechanics
C--------------------------------------------------------

      SUBROUTINE QMF2( SCRD,FRCS,FRCN )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
!     include 'atoms.i'

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC =  36 )
      PARAMETER ( DCUT  = 1.D-7 ) 
      PARAMETER ( ALPHA = 1.0D0    )
C     PARAMETER ( NLINK1 = 10 )
      PARAMETER ( RRCUT = 4.724D0 )
      PARAMETER ( RRCUT1= 5.669D0 )
      PARAMETER ( DRCUT = 0.945D0 )

      PARAMETER ( PI    = 3.14159265358979323D0 )
      
      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PRMT4 / nion1,iiont(maxatm)
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / LJQMMM / SQMMM(NNUC,maxatm),EPQMMM(NNUC,maxatm),
     *                  chgmm(maxatm),chgqm(NNUC)

      DIMENSION SCRD( NDIM,maxatm )

C     DIMENSION FNSS( maxatm,NDIM )              ! nuclear-site and LJ force on site
      DIMENSION FRCS( maxatm,NDIM )              ! QM/MM force on site

      DIMENSION FRCN( NNUC,NDIM )                ! QM/MM force on nuclear

      DV=DX*DY*DZ

C
C     number of loop in QM and MM subsystems
C
      nlqm = NNUC-NLINK
      nlqm1 = nlqm+1
C     nlmm = n-NNUC
      nlmm = n-NNUC+NLINK

      DO L=1,3
      DO in=1,n

          FRCS( in,L )   = 0.D0
            
      ENDDO
      ENDDO
     
      DO L=1,3
      DO in=1,NNUC
     
          FRCN( in,L )   = 0.D0
            
      ENDDO
      ENDDO

C-----nuclear-site force CALCULATION
C-----LJ force CALCULATION

      ONE    = 1.D0
      TWO    = 2.D0
      FOUR   = 4.D0
      SIX    = 6.D0
      RRCUT2 = ONE/(RRCUT*RRCUT)

      DO 90 in=1,nlmm
      DO 100 M=1,nlqm 

               RX=PNUC( M,1 ) - SCRD( 1,in )
        RY=PNUC( M,2 ) - SCRD( 2,in )
        RZ=PNUC( M,3 ) - SCRD( 3,in )
        R2=RX**2+RY**2+RZ**2
        R   = DSQRT(R2)
        RI  = ONE/R
        R2I = RI*RI

        ASIG2  = SQMMM(M,in)*R2I
        ASIG6  = ASIG2*ASIG2*ASIG2
        ASIG12 = ASIG6*ASIG6
        BSIG   = TWO*ASIG12 - ASIG6
        TFEPS  = SIX*EPQMMM(M,in)
        CREPS  = TFEPS*R2I*BSIG
        IF(CREPS.GT.0.1D0) THEN
          WRITE(*,*) 'QM,MM,CREPS(OLD) = ',M,in,CREPS
C         CREPS = 0.1D0
          CREPS = 0.01D0
C         WRITE(*,*) 'QM,MM,CREPS(NEW) = ',M,in,CREPS
        ENDIF
        CREPSX = CREPS*RX
        CREPSY = CREPS*RY
        CREPSZ = CREPS*RZ

        FRCS(in,1) = FRCS(in,1) - CREPSX   ! LJ force on site
        FRCS(in,2) = FRCS(in,2) - CREPSY
        FRCS(in,3) = FRCS(in,3) - CREPSZ

        FRCN(M,1)  = FRCN(M,1) + CREPSX    ! LJ force on nuclear
        FRCN(M,2)  = FRCN(M,2) + CREPSY
        FRCN(M,3)  = FRCN(M,3) + CREPSZ

       IF(NMMSW(in).EQ.1) GOTO 100

       IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
         R3I   = RI*R2I
         ACHR  = chgqm( M )*chgmm( in )*R3I
       ELSEIF(R.GT.DRCUT) THEN
         ARD   = R - DRCUT
         ARD2  = ARD*ARD*RRCUT2
         ESW1  = ONE - ARD2
         ESW2  = ESW1*ESW1
         ESW   = ONE - ESW2
         ESWA  = FOUR*ESW1*ARD*RRCUT2
         ZVPCH = chgqm( M )*chgmm( in )*R2I
         ACHR  = ESW*ZVPCH*RI - ESWA*ZVPCH
       ELSE
         GOTO 100
       ENDIF

       ABSCHR = ABS(ACHR)
       IF(ABSCHR.GT.1.0D0) THEN
         WRITE(*,*) 'QM,MM,ACHR(OLD) = ',M,in,ACHR
         ACHR = 0.01D0*ACHR/ABSCHR
C        WRITE(*,*) 'QM,MM,ACHR(NEW) = ',M,in,ACHR
       ENDIF

        CHRX = ACHR*RX
        CHRY = ACHR*RY
        CHRZ = ACHR*RZ

        FRCS(in,1) = FRCS(in,1) - CHRX      ! Coulomb force on site
        FRCS(in,2) = FRCS(in,2) - CHRY
        FRCS(in,3) = FRCS(in,3) - CHRZ

        FRCN(M,1)  = FRCN(M,1) + CHRX       ! Coulomb force on nuclear
        FRCN(M,2)  = FRCN(M,2) + CHRY
        FRCN(M,3)  = FRCN(M,3) + CHRZ

100   CONTINUE
90    CONTINUE
      
C     DO L=1,3
C     DO in=1,nlmm
C       FRCS(in,L) = FDSS(in,L) - FNSS(in,L)                 !QM/MM force on site
C     ENDDO
C     ENDDO

      RETURN
      END


C--------------------------------------------------------
C     SUBROUTINE SPLITTING FUZZY CELLS
C--------------------------------------------------------
 
      SUBROUTINE FUZZYC
 
      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
 
      include "mpif.h"
      include "mpi.i"
 
      include 'QMpara.i'                        ! Vmol
 
C     PARAMETER ( NMAX = 80 )
      PARAMETER ( PI = 3.14159265358979323D0 )
 
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
 
      COMMON / GRID2 / DX, DY, DZ
      COMMON / NCLR  / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / FCELL / WNUC(NFUZZY,NMAXX/NX,NMAXY/NY,NMAXZ/NZ)
C     COMMON / ELCTRN / AX1, AX2, AY1, AY2, AZ1, AZ2, AO,
C    *                  AX3, AX4, AY3, AY4, AZ3, AZ4
      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,
C    *                  WNUC(NFUZZY,NMAX,NMAX,NMAX),RSIZE(NNUC),
C    *                  WNUC(NFUZZY,NMAX/NX,NMAX/NY,NMAX/NZ),
     *                  RSIZE(NNUC),
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF
C
CC     DIMENSION PGRID( NNUC )
CC     DIMENSION ASIZE( 1:NNUC-1,2:NNUC )
       DIMENSION PGRID( NFUZZY )
       DIMENSION ASIZE( NFUZZY*NFUZZY )
 
       CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
       CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )
C
CC------------------
C      TWO  = 2.d0
C      TAO  = TWO*AO
C      TAX1 = TWO*AX1
C      TAX2 = TWO*AX2
C      TAX3 = TWO*AX3
C      TAX4 = TWO*AX4
C      TAY1 = TWO*AY1
C      TAY2 = TWO*AY2
C      TAY3 = TWO*AY3
C      TAY4 = TWO*AY4
C      TAZ1 = TWO*AZ1
C      TAZ2 = TWO*AZ2
C      TAZ3 = TWO*AZ3
C      TAZ4 = TWO*AZ4
C      PI4  = 4.d0*PI
CC------------------
 
      ONE    = 1.D0
C     NNMAXX = NMAXX/2+1
C     NNMAXY = NMAXY/2+1
C     NNMAXZ = NMAXZ/2+1
 
C      NNMAXP = NMAX/2
C      NNX    = (1-MX)*NNMAXP + 1
C      NNY    = (1-MY)*NNMAXP + 1
C      NNZ    = (1-MZ)*NNMAXP + 1
 
      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

       DO K  = 1, JZ
       DO J  = 1, JY
       DO I  = 1, JX
       DO NA = 1, NFUZZY
         WNUC(NA,I,J,K)=ONE
       ENDDO
       ENDDO
       ENDDO
       ENDDO
 
       NF = 0
       DO NA = 1, NNUC
C        IF( ZA(NA) .GT. 1.D0 ) THEN
           NF = NF + 1
           NPOP(NF) = NA
C        ENDIF
       ENDDO
 
       NF1 = 0
       DO NA1 = 1, NF-1
       DO NA2 = NA1+1, NF
         RU = (RSIZE(NPOP(NA1))-RSIZE(NPOP(NA2)))
     *     / (RSIZE(NPOP(NA1))+RSIZE(NPOP(NA2)))
         NF1 = NF1 + 1
         ASIZE(NF1) = RU/(RU*RU-1.D0)
         AS    = ASIZE(NF1)
         ABSAS = ABS(AS)
         IF(ABSAS.GT.0.5D0) ASIZE(NF1) = 0.5D0*AS/ABSAS
C        write(*,*) 'ASIZE=',ASIZE(NA1,NA2)
       ENDDO
       ENDDO
 
       DO 10 K = 1, JZ
 
         RZ1 = DZ*( K-NNZ )
 
         DO 10 J = 1, JY
 
           RY1 = DY*( J-NNY )
 
           DO 10 I = 1, JX
 
             RX1 = DX*( I-NNX )
 
             NF1 = 0
             DO 10 NA1 = 1, NF-1
 
               NNA1 = NPOP(NA1)
               RX = RX1 - PNUC(NNA1,1) 
               RY = RY1 - PNUC(NNA1,2)
               RZ = RZ1 - PNUC(NNA1,3)
               R2 = RX**2 + RY**2 + RZ**2
               RNA1 = DSQRT(R2)
 
               DO 10 NA2 = NA1+1, NF
 
                 NNA2 = NPOP(NA2)
                 RX = RX1 - PNUC(NNA2,1) 
                 RY = RY1 - PNUC(NNA2,2)
                 RZ = RZ1 - PNUC(NNA2,3)
                 R2 = RX**2 + RY**2 + RZ**2
                 RNA2 = DSQRT(R2)
 
                 RX = PNUC(NNA1,1) - PNUC(NNA2,1) 
                 RY = PNUC(NNA1,2) - PNUC(NNA2,2)
                 RZ = PNUC(NNA1,3) - PNUC(NNA2,3)
                 R2 = RX**2 + RY**2 + RZ**2
                 RZA = DSQRT(R2)
 
                 RU = (RNA1-RNA2)/RZA
                 NF1 = NF1 + 1
                 RV = RU + ASIZE(NF1)*(1.D0-RU*RU)
 
                 FV = 1.5D0*RV-0.5D0*RV*RV*RV
                 FV = 1.5D0*FV-0.5D0*FV*FV*FV
                 FV = 1.5D0*FV-0.5D0*FV*FV*FV
 
                 HFV = 0.5D0*FV
                 SRU12 = 0.5D0-HFV
                 SRU21 = 0.5D0+HFV
 
                 WNUC(NA1,I,J,K)=WNUC(NA1,I,J,K)*SRU12
                 WNUC(NA2,I,J,K)=WNUC(NA2,I,J,K)*SRU21
                 IF(SRU12.LT.0.d0) THEN
                   write(*,*) 'SRU12=',SRU12
                 ENDIF
                 IF(SRU21.LT.0.d0) THEN
                   write(*,*) 'SRU21=',SRU21
                 ENDIF
 
10    CONTINUE
 
       DO K  = 1, JZ
       DO J  = 1, JY
       DO I  = 1, JX
         SPGRID = 0.D0
         DO NA = 1, NF
           PGRID(NA) = WNUC(NA,I,J,K)
           SPGRID = SPGRID + PGRID(NA)
         ENDDO
         SPGI = 1.D0/SPGRID
         DO NA = 1, NF
           WNUC(NA,I,J,K) = PGRID(NA)*SPGI
         ENDDO
       ENDDO
       ENDDO
       ENDDO
 
       RETURN
       END
 
C--------------------------------------------------------
C     SUBROUTINE Effective Charges by FUZZY CELL method 
C--------------------------------------------------------
 
      SUBROUTINE EFTC(RHOA,RHOB)
 
      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
 
      include "mpif.h"
      include "mpi.i"
 
      include 'QMpara.i'                        ! Vmol
 
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
 
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / GRID2 / DX, DY, DZ
      COMMON / FCELL / WNUC(NFUZZY,NMAXX/NX,NMAXY/NY,NMAXZ/NZ)
C     COMMON / ELCTRN / AX1, AX2, AY1, AY2, AZ1, AZ2, AO,
      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,
C    *                  WNUC(NFUZZY,NMAX,NMAX,NMAX),RSIZE(NNUC),
C    *                  WNUC(NFUZZY,NMAX/NX,NMAX/NY,NMAX/NZ),
     *                  RSIZE(NNUC),
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF
 
       DIMENSION TZPOP( NFUZZY )
 
       DIMENSION RHO(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
       DIMENSION RHOA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
       DIMENSION RHOB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
 
       CHARACTER ATOM*2

       CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
       CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )
 
       RHO(:,:,:) = RHOA(:,:,:) - RHOB(:,:,:)    ! spin density

       DO NA = 1, NF
         TZPOP(NA) = 0.D0
       ENDDO
 
       DV=DX*DY*DZ
 
       DO K  = 1, JZ
       DO J  = 1, JY
       DO I  = 1, JX
         TMPRHO = RHO(I,J,K)
         DO NA = 1, NF
           TZPOP(NA) =  TZPOP(NA) + WNUC(NA,I,J,K)*TMPRHO  
         ENDDO
       ENDDO
       ENDDO
       ENDDO
 
       CALL MPI_ALLREDUCE( TZPOP(1),ZPOP(1),NF,MPI_DOUBLE_PRECISION,
     *                    MPI_SUM,MPI_COMM_WORLD,IERR )
 
       DO NA = 1, NF
         ZPOP(NA) = ZPOP(NA)*DV
       ENDDO
 
      IF( MYID .EQ. 0 ) THEN

C     WRITE(*,*)
C     WRITE(*,*) '----- Spin density on site (Fuzzy cell) -----'
      DO K = 1, NF
        NNA = NPOP(K)
      IF(ZA(K).EQ.1.D0) THEN
        ATOM = 'H '
      ELSEIF(ZA(K).EQ.6.D0) THEN
        ATOM = 'C '
      ELSEIF(ZA(K).EQ.7.D0) THEN
        ATOM = 'N '
      ELSEIF(ZA(K).EQ.8.D0) THEN
        ATOM = 'O '
      ELSEIF(ZA(K).EQ.17.D0) THEN
        ATOM = 'Cl'
      ELSEIF(ZA(K).EQ.20.D0) THEN
        ATOM = 'Ca'
      ELSEIF(ZA(K).EQ.25.D0) THEN
        ATOM = 'Mn'
      ENDIF

C     WRITE(*,91) K,ATOM,ZPOP(K),RSIZE(K)

      ENDDO

C     WRITE(*,*) '---------------------------------------------'

      ENDIF
 
91    FORMAT(2X,I4,2X,A2,' : ',F10.6,5X,F6.4)

      RETURN
      END
 
C----------------------------------------------
C   SUBROUTINE COULOMB POTENTIAL IN REAL SPACE
C----------------------------------------------

C----- RVCLMB -----

      SUBROUTINE RVCLMB( RHO )

      IMPLICIT REAL*8( A-H,O-Z )
      IMPLICIT INTEGER*4( I-N )

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'

      PARAMETER( INUM = 8 )
      PARAMETER( INUM2 = INUM/2 )
      PARAMETER( PI = 3.14159265358979323D0 )
      PARAMETER( NITMAX = 100  )
      PARAMETER( NCG    = 10  )
      PARAMETER( EPS1 = 1.0D0 )
C     PARAMETER( EPS2 = 1.0D-05 )
      PARAMETER( EPS2 = 1.0D-04 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          FREEZE*5,PRINT*5,DGF*3

      DIMENSION TA(  3 )
      DIMENSION TA1( 3 )

      DIMENSION RHO(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DMMY(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION BRHO(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION PVEC(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RVEC(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION APVEC( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
C     DIMENSION VCOUX( INUM,NMAX/NY,NMAX/NZ ) 
C     DIMENSION VCOUY( NMAX/NX,INUM,NMAX/NZ ) 

      DIMENSION VCOUX( NMAXY/NY,NMAXZ/NZ,INUM ) ! note: Dimension(Y,Z,X)
      DIMENSION VCOUY( NMAXX/NX,NMAXZ/NZ,INUM ) ! note: Dimension(X,Z,Y)
      DIMENSION VCOUZ( NMAXX/NX,NMAXY/NY,INUM )
      DIMENSION VBC  ( 1-INUM2:NMAXX/NX+INUM2,1-INUM2:NMAXY/NY+INUM2,
     *                 1-INUM2:NMAXZ/NZ+INUM2)

      COMMON / MMPI4 / JX,JY,JZ

      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / CLMB  / PCLMB, PEXC, PEXT
      COMMON / GRID2 / DX,DY,DZ
C     COMMON / OFC / VCOU( NMAX,NMAX,NMAX )
      COMMON / OFC / VCOU( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,
C    *                  WNUC(NFUZZY,NMAX,NMAX,NMAX),RSIZE(NNUC),
C    *                  WNUC(NFUZZY,NMAX/NX,NMAX/NY,NMAX/NZ),
     *                  RSIZE(NNUC),
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      DV = DX*DY*DZ

C----- BOUNDARY CONDITION -----

C     DMMY(:,:,:) = 0.D0
C     CALL EFTC(RHO,DMMY)
      CALL OPTFC2

      CALL RCLMB2( VBC )
      PVEC(:,:,:) = 0.D0
      CALL RCLMB5( PVEC,VBC,BRHO )
      BRHO( :,:,: ) = PI4 * RHO( :,:,: ) - BRHO( :,:,: )

C------- takahashi MPI check 05.08.01-----------------
      ZSUM = 0.D0
      DO I = 1, JX
      DO J = 1, JY 
      DO K = 1, JZ
        ZSUM = ZSUM + BRHO( I,J,K )**2*DV
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE( ZSUM,SZSUM,1,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

C     IF(MYID.EQ.0) THEN
C     write(6,*) 'zsum = ',SZSUM
C     ENDIF
C-----------------------------------------------------

C----- STEEPEST DESCENT -----

10    CONTINUE

C     write(*,*) 'steepest descent'

      VBC(:,:,:) = 0.D0
      CALL RCLMB5( VCOU,VBC,APVEC )

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        RVEC( I,J,K ) = BRHO( I,J,K ) - APVEC( I,J,K )
        PVEC( I,J,K ) = RVEC( I,J,K ) 
      ENDDO
      ENDDO
      ENDDO

C----- CONJUGATE GRADIENT -----
C----- ITERATION -----

C     write(*,*) 'conjugate gradient'

      DO 20 L = 1, NITMAX
  
        CALL RCLMB5( PVEC,VBC,APVEC )
  
        TA1(1) = 0.D0
        TA1(2) = 0.D0
        DO K = 1, JZ
        DO J = 1, JY
        DO I = 1, JX
          TA1(1) = TA1(1) + RVEC( I,J,K )*RVEC( I,J,K )
          TA1(2) = TA1(2) + PVEC( I,J,K )*APVEC( I,J,K )
        ENDDO
        ENDDO
        ENDDO
        
        CALL MPI_ALLREDUCE( TA1(1),TA(1),2,MPI_DOUBLE_PRECISION,
     *                      MPI_SUM,MPI_COMM_WORLD,IERR )

        ALPHA = TA(1) / TA(2) 
  
        TA1(3) = 0.D0
        DO K = 1, JZ
        DO J = 1, JY
        DO I = 1, JX
          VCOU( I,J,K ) = VCOU( I,J,K ) + ALPHA*PVEC(  I,J,K )
          RVEC( I,J,K ) = RVEC( I,J,K ) - ALPHA*APVEC( I,J,K )
          TA1(3) = TA1(3) + RVEC( I,J,K )**2
        ENDDO
        ENDDO
        ENDDO
  
        CALL MPI_ALLREDUCE( TA1(3),TA(3),1,MPI_DOUBLE_PRECISION,
     *                      MPI_SUM,MPI_COMM_WORLD,IERR )

        DIST = DSQRT(TA(3)*DV)
C       IF(MOD(L,10).EQ.1) WRITE(*,*) 'DIST = ',L,DIST
C       IF(MOD(L,20).EQ.1) write(*,*) 'alpha,DIST=',ALPHA,DIST
C       IF( DIST .LT. EPS2 .OR. DIST .GE. EPS1 ) THEN
        IF( DIST .LT. EPS2 ) THEN
          GOTO 60
        ENDIF
       
        BETA = TA(3)/TA(1)
  
        DO K = 1, JZ
        DO J = 1, JY
        DO I = 1, JX
          PVEC( I,J,K ) = RVEC( I,J,K ) + BETA*PVEC( I,J,K )
        ENDDO
        ENDDO
        ENDDO
  
20    CONTINUE
C     WRITE(*,*) 'DIST = ',L,DIST
      L = L - 1

60    CONTINUE
C     IF( DIST .GE. EPS1 ) GOTO 10
      IF(MYID.EQ.0) THEN
      IF(PRINT.EQ.'LARGE') THEN
        WRITE(*,*)                     
        WRITE(*,*) ' Poisson It,Dist = ',L,DIST
      ENDIF
      ENDIF

C----- classical coulomb potential -----   

      TPCLMB = 0.D0

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        TPCLMB = TPCLMB + VCOU( I,J,K )*RHO( I,J,K )
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE( TPCLMB,PCLMB,1,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

      PCLMB = 0.5D0*PCLMB*DV

C     IF(MYID.EQ.0) THEN
C     WRITE(*,*) 'PCLMB=',PCLMB
C     ENDIF

C     REWIND( 77 )

C     DO I = 1, NMAX
C     DO J = 1, NMAX 
C     WRITE(77,*) VCOU( I,J,NMAX/2 )
C     ENDDO
C     ENDDO

      RETURN
      END

 
C----------------------------------------------
C   SUBROUTINE COULOMB POTENTIAL IN REAL SPACE
C   TO BE USED SPECIFICALLY IN SUB. OEP_YW.F 
C----------------------------------------------

C----- RVCLMB -----

      SUBROUTINE RVCLMB_OEP( RHO,VCOU )

      IMPLICIT REAL*8( A-H,O-Z )
      IMPLICIT INTEGER*4( I-N )

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'

      PARAMETER( INUM = 8 )
      PARAMETER( INUM2 = INUM/2 )
      PARAMETER( PI = 3.14159265358979323D0 )
      PARAMETER( NITMAX = 100  )
      PARAMETER( NCG    = 10  )
      PARAMETER( EPS1 = 1.0D0 )
C     PARAMETER( EPS2 = 1.0D-05 )
      PARAMETER( EPS2 = 1.0D-04 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          FREEZE*5,PRINT*5,DGF*3

      DIMENSION TA(  3 )
      DIMENSION TA1( 3 )

      DIMENSION RHO(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DMMY( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION BRHO( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION PVEC( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RVEC( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION APVEC(NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VCOU( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

C     DIMENSION VCOUX( INUM,NMAX/NY,NMAX/NZ ) 
C     DIMENSION VCOUY( NMAX/NX,INUM,NMAX/NZ ) 

      DIMENSION VCOUX( NMAXY/NY,NMAXZ/NZ,INUM ) ! note: Dimension(Y,Z,X)
      DIMENSION VCOUY( NMAXX/NX,NMAXZ/NZ,INUM ) ! note: Dimension(X,Z,Y)
      DIMENSION VCOUZ( NMAXX/NX,NMAXY/NY,INUM )
      DIMENSION VBC  ( 1-INUM2:NMAXX/NX+INUM2,1-INUM2:NMAXY/NY+INUM2,
     *                 1-INUM2:NMAXZ/NZ+INUM2)

      COMMON / MMPI4 / JX,JY,JZ

      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / CLMB  / PCLMB, PEXC, PEXT
      COMMON / GRID2 / DX,DY,DZ
C     COMMON / OFC / VCOU( NMAX,NMAX,NMAX )
C     COMMON / OFC / VCOU( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )    ! should not be in the common block 
      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,
C    *                  WNUC(NFUZZY,NMAX,NMAX,NMAX),RSIZE(NNUC),
C    *                  WNUC(NFUZZY,NMAX/NX,NMAX/NY,NMAX/NZ),
     *                  RSIZE(NNUC),
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      DV = DX*DY*DZ

C----- BOUNDARY CONDITION -----

C     DMMY(:,:,:) = 0.D0
C     CALL EFTC(RHO,DMMY)
      CALL OPTFC2

      CALL RCLMB2( VBC )
      PVEC(:,:,:) = 0.D0
      CALL RCLMB5( PVEC,VBC,BRHO )
      BRHO( :,:,: ) = PI4 * RHO( :,:,: ) - BRHO( :,:,: )

C------- takahashi MPI check 05.08.01-----------------
      ZSUM = 0.D0
      DO I = 1, JX
      DO J = 1, JY 
      DO K = 1, JZ
        ZSUM = ZSUM + BRHO( I,J,K ) **2*DV
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE( ZSUM,SZSUM,1,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

C     IF(MYID.EQ.0) THEN
C     write(6,*) 'zsum = ',SZSUM
C     ENDIF
C-----------------------------------------------------

C----- STEEPEST DESCENT -----

10    CONTINUE

C     write(*,*) 'steepest descent'

      VBC(:,:,:) = 0.D0
      CALL RCLMB5( VCOU,VBC,APVEC )

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        RVEC( I,J,K ) = BRHO( I,J,K ) - APVEC( I,J,K )
        PVEC( I,J,K ) = RVEC( I,J,K ) 
      ENDDO
      ENDDO
      ENDDO

C----- CONJUGATE GRADIENT -----
C----- ITERATION -----

C     write(*,*) 'conjugate gradient'

      DO 20 L = 1, NITMAX
  
        CALL RCLMB5( PVEC,VBC,APVEC )
  
        TA1(1) = 0.D0
        TA1(2) = 0.D0
        DO K = 1, JZ
        DO J = 1, JY
        DO I = 1, JX
          TA1(1) = TA1(1) + RVEC( I,J,K )*RVEC( I,J,K )
          TA1(2) = TA1(2) + PVEC( I,J,K )*APVEC( I,J,K )
        ENDDO
        ENDDO
        ENDDO
        
        CALL MPI_ALLREDUCE( TA1(1),TA(1),2,MPI_DOUBLE_PRECISION,
     *                      MPI_SUM,MPI_COMM_WORLD,IERR )

        ALPHA = TA(1) / TA(2) 
  
        TA1(3) = 0.D0
        DO K = 1, JZ
        DO J = 1, JY
        DO I = 1, JX
          VCOU( I,J,K ) = VCOU( I,J,K ) + ALPHA*PVEC( I,J,K )
          RVEC( I,J,K ) = RVEC( I,J,K ) - ALPHA*APVEC( I,J,K )
          TA1(3) = TA1(3) + RVEC( I,J,K )**2
        ENDDO
        ENDDO
        ENDDO
  
        CALL MPI_ALLREDUCE( TA1(3),TA(3),1,MPI_DOUBLE_PRECISION,
     *                      MPI_SUM,MPI_COMM_WORLD,IERR )

        DIST = DSQRT(TA(3)*DV)
C       IF(MOD(L,10).EQ.1) WRITE(*,*) 'DIST = ',L,DIST
C       IF(MOD(L,20).EQ.1) write(*,*) 'alpha,DIST=',ALPHA,DIST
C       IF( DIST .LT. EPS2 .OR. DIST .GE. EPS1 ) THEN
        IF( DIST .LT. EPS2 ) THEN
          GOTO 60
        ENDIF
       
        BETA = TA(3)/TA(1)
  
        DO K = 1, JZ
        DO J = 1, JY
        DO I = 1, JX
          PVEC( I,J,K ) = RVEC( I,J,K ) + BETA*PVEC( I,J,K )
        ENDDO
        ENDDO
        ENDDO
  
20    CONTINUE
C     WRITE(*,*) 'DIST = ',L,DIST
      L = L - 1

60    CONTINUE
C     IF( DIST .GE. EPS1 ) GOTO 10
      IF(MYID.EQ.0) THEN
      IF(PRINT.EQ.'LARGE') THEN
        WRITE(*,*)                     
        WRITE(*,*) ' Poisson It,Dist = ',L,DIST
      ENDIF
      ENDIF

C----- classical coulomb potential -----   

      TPCLMB = 0.D0

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        TPCLMB = TPCLMB + VCOU( I,J,K )*RHO( I,J,K )
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE( TPCLMB,PCLMB,1,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

      PCLMB = 0.5D0*PCLMB*DV

C     IF(MYID.EQ.0) THEN
C     WRITE(*,*) 'PCLMB=',PCLMB
C     ENDIF

C     REWIND( 77 )

C     DO I = 1, NMAX
C     DO J = 1, NMAX 
C     WRITE(77,*) VCOU( I,J,NMAX/2 )
C     ENDDO
C     ENDDO

      RETURN
      END

C----- RCLMB2 -----
C Make Boundary Condition of VCOU
C     SUBROUTINE RCLMB2( VBC )

C     IMPLICIT REAL*8( A-H,O-Z )
C     IMPLICIT INTEGER*4( I-N )

C     include "mpif.h"
C     include "mpi.i"
C     include 'QMpara.i'

C     PARAMETER( INUM  = 8 )
C     PARAMETER( INUM2 = INUM/2 )

C     COMMON / MMPI1 / MX,MY,MZ
C     COMMON / MMPI4 / JX,JY,JZ

C     COMMON / GRID2 / DX,DY,DZ
C     COMMON / ACELL  / XL,YL,ZL
C     COMMON / NCLR / PNUC( NNUC,NDIM )
C     COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,
C    *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,
C    *                  RSIZE(NNUC),
C    *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF

C     DIMENSION VBC( 1-INUM2:NMAXX/NX+INUM2,1-INUM2:NMAXY/NY+INUM2,
C    *               1-INUM2:NMAXZ/NZ+INUM2 )

C     CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
C     CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C----- INITIALIZE -----

C     VBC(:,:,:) = 0.D0

C     NNMAXPX = NMAXX/2
C     NNMAXPY = NMAXY/2
C     NNMAXPZ = NMAXZ/2

C     NNX     = -JX*MX + NNMAXPX + 1
C     NNY     = -JY*MY + NNMAXPY + 1
C     NNZ     = -JZ*MZ + NNMAXPZ + 1

C     DO K = 1-INUM2, JZ+INUM2
C     DO J = 1-INUM2, JY+INUM2
C     DO I = 1-INUM2, JX+INUM2
C     DO L = 1, NF

C       RX1 = DX*( I-NNX ) 
C       RY1 = DY*( J-NNY )
C       RZ1 = DZ*( K-NNZ )
C       
C       NNA = NPOP(L)
C       RX = RX1 - PNUC( NNA,1 )
C       RY = RY1 - PNUC( NNA,2 )
C       RZ = RZ1 - PNUC( NNA,3 )
C 
C       R = RX**2 + RY**2 + RZ**2
C       R = DSQRT( R )
C       VBC( I,J,K ) = VBC( I,J,K ) + ZPOP( NNA ) / R

C     ENDDO
C     ENDDO
C     ENDDO
C     ENDDO

C     VBC(1:JX,1:JY,1:JZ) = 0.D0
C     IF(MX.NE.0   ) VBC(1-INUM2:0    ,:,:) = 0.D0
C     IF(MX.NE.NX-1) VBC(JX+1:JX+INUM2,:,:) = 0.D0
C     IF(MY.NE.0   ) VBC(:,1-INUM2:0    ,:) = 0.D0
C     IF(MY.NE.NY-1) VBC(:,JY+1:JY+INUM2,:) = 0.D0
C     IF(MZ.NE.0   ) VBC(:,:,1-INUM2:0    ) = 0.D0
C     IF(MZ.NE.NZ-1) VBC(:,:,JZ+1:JZ+INUM2) = 0.D0

C     RETURN
C     END

C----- RCLMB2 -----
C Make Boundary Condition of VCOU
      SUBROUTINE RCLMB2( VBC )

      IMPLICIT REAL*8( A-H,O-Z )
      IMPLICIT INTEGER*4( I-N )

      include "mpif.h"
      include "mpi.i"
      include 'QMpara.i'

      PARAMETER( INUM  = 8 )
      PARAMETER( INUM2 = INUM/2 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX,DY,DZ
      COMMON / ACELL  / XL,YL,ZL
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,
     *                  RSIZE(NNUC),
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF

      DIMENSION VBC( 1-INUM2:NMAXX/NX+INUM2,1-INUM2:NMAXY/NY+INUM2,
     *               1-INUM2:NMAXZ/NZ+INUM2 )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C----- INITIALIZE -----

      VBC(:,:,:) = 0.D0

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2

      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO IBOUND = 1, 6
         SELECT CASE (IBOUND)
         CASE(1); IF(MX.NE.0   ) CYCLE
         CASE(2); IF(MY.NE.0   ) CYCLE
         CASE(3); IF(MZ.NE.0   ) CYCLE
         CASE(4); IF(MX.NE.NX-1) CYCLE
         CASE(5); IF(MY.NE.NY-1) CYCLE
         CASE(6); IF(MZ.NE.NZ-1) CYCLE
         END SELECT

         IBGN = 1; IEND = JX
         JBGN = 1; JEND = JY
         KBGN = 1; KEND = JZ
         SELECT CASE (IBOUND)
         CASE(1); IBGN = 1-INUM2; IEND = 0
         CASE(2); JBGN = 1-INUM2; JEND = 0
         CASE(3); KBGN = 1-INUM2; KEND = 0
         CASE(4); IBGN = JX+1;    IEND = JX+INUM2
         CASE(5); JBGN = JY+1;    JEND = JY+INUM2
         CASE(6); KBGN = JZ+1;    KEND = JZ+INUM2
         END SELECT

         DO K = KBGN, KEND
         DO J = JBGN, JEND
         DO I = IBGN, IEND
         DO L = 1, NF

            RX1 = DX*( I-NNX ) 
            RY1 = DY*( J-NNY )
            RZ1 = DZ*( K-NNZ )
            
            NNA = NPOP(L)
            RX = RX1 - PNUC( NNA,1 )
            RY = RY1 - PNUC( NNA,2 )
            RZ = RZ1 - PNUC( NNA,3 )
            
            R = RX**2 + RY**2 + RZ**2
            R = DSQRT( R )
            VBC( I,J,K ) = VBC( I,J,K ) + ZPOP( NNA ) / R

         ENDDO
         ENDDO
         ENDDO
         ENDDO

      ENDDO

      RETURN
      END

C----- RCLMB4 -----

C      SUBROUTINE RCLMB4( PVEC,APVEC )
C
C      IMPLICIT REAL*8( A-H,O-Z )
C      IMPLICIT INTEGER*4( I-N )
C
C      include 'QMpara.i'
C
C      PARAMETER( INUM = 8 )
C
C      PARAMETER( PI = 3.14159265358979323D0 )
C
C      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,
C     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,
CC     *                  WNUC(NFUZZY,NMAX,NMAX,NMAX),RSIZE(NNUC),
CC     *                  WNUC(NFUZZY,NMAX/NX,NMAX/NY,NMAX/NZ),
C     *                  RSIZE(NNUC),
C     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF
C
C      DIMENSION PVEC( NMAX,NMAX,NMAX )
C      DIMENSION PVEC1( NMAX+INUM )
C      DIMENSION APVEC( NMAX,NMAX,NMAX )
C     
C      ZERO     = 0.D0
C      PVEC1(1) = ZERO
C      PVEC1(2) = ZERO
C      PVEC1(3) = ZERO
C      PVEC1(4) = ZERO
C      PVEC1(NMAX+5) = ZERO
C      PVEC1(NMAX+6) = ZERO
C      PVEC1(NMAX+7) = ZERO
C      PVEC1(NMAX+8) = ZERO
C
C      NS = 5
C      NE = NMAX+4
C
CC----- Z-AXIS -----
C      
C      DO 10 J = 1, NMAX
C        DO 20 I = 1, NMAX
C
C        DO K = NS, NE
C          PVEC1(K) = PVEC( I,J,K-4 )
C        ENDDO
C
C        DO 30 K = NS, NE
C
C        APVEC( I,J,K-4 ) = TAZ1 * ( PVEC1( K-1 )
C     *                   +         PVEC1( K+1 ))
C     *                   + TAZ2 * ( PVEC1( K-2 )
C     *                   +         PVEC1( K+2 ))
C     *                   + TAZ3 * ( PVEC1( K-3 )
C     *                   +         PVEC1( K+3 ))
C     *                   + TAZ4 * ( PVEC1( K-4 )
C     *                   +         PVEC1( K+4 ))
C     *                   + TAO  *   PVEC1( K )
C
C30      CONTINUE
C
C20      CONTINUE
C10    CONTINUE
C
CC----- Y-AXIS -----
C      
C      DO 40 K = 1, NMAX
C        DO 50 I = 1, NMAX
C
C        DO J = NS, NE
C          PVEC1(J) = PVEC( I,J-4,K )
C        ENDDO
C
C        DO 60 J = NS, NE
C
C        APVEC( I,J-4,K ) = TAY1 * ( PVEC1( J-1 )
C     *                   +         PVEC1( J+1 ))
C     *                   + TAY2 * ( PVEC1( J-2 )
C     *                   +         PVEC1( J+2 ))
C     *                   + TAY3 * ( PVEC1( J-3 )
C     *                   +         PVEC1( J+3 ))
C     *                   + TAY4 * ( PVEC1( J-4 )
C     *                   +         PVEC1( J+4 ))
C     *                   + APVEC( I,J-4,K )
C
C60      CONTINUE
C
C50      CONTINUE
C40    CONTINUE
C
CC----- X-AXIS -----
C      
C      DO 70 K = 1, NMAX
C        DO 80 J = 1, NMAX
C
C        DO I = NS, NE
C          PVEC1(I) = PVEC( I-4,J,K )
C        ENDDO
C
C        DO 90 I = NS, NE
C
C        APVEC( I-4,J,K ) = TAX1 * ( PVEC1( I-1 )
C     *                   +         PVEC1( I+1 ))
C     *                   + TAX2 * ( PVEC1( I-2 )
C     *                   +         PVEC1( I+2 ))
C     *                   + TAX3 * ( PVEC1( I-3 )
C     *                   +         PVEC1( I+3 ))
C     *                   + TAX4 * ( PVEC1( I-4 )
C     *                   +         PVEC1( I+4 ))
C     *                   + APVEC( I-4,J,K )
C
C90      CONTINUE
C
C80      CONTINUE
C70    CONTINUE
C
C      RETURN
C      END

C----- RCLMB4 -----

      SUBROUTINE RCLMB4( PVEC,APVEC )

      IMPLICIT REAL*8 ( A-H,O-Z )
      IMPLICIT INTEGER*4 ( I-N )

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'
 
      PARAMETER( PI = 3.14159265358979323D0 )

      COMMON / MMPI4 / JX,JY,JZ

      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,
C    *                  WNUC(NFUZZY,NMAX,NMAX,NMAX),RSIZE(NNUC),
C    *                  WNUC(NFUZZY,NMAX/NX,NMAX/NY,NMAX/NZ),
     *                  RSIZE(NNUC),
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF
 
      DIMENSION PVEC(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION APVEC( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
     
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C----- Z-AXIS -----
      
      DO K = 5, JZ-4
      DO J = 1, JY
      DO I = 1, JX

        APVEC( I,J,K ) = TAZ1 * ( PVEC( I,J,K-1 )
     *                 +         PVEC( I,J,K+1 ))
     *                 + TAZ2 * ( PVEC( I,J,K-2 )
     *                 +         PVEC( I,J,K+2 ))
     *                 + TAZ3 * ( PVEC( I,J,K-3 )
     *                 +         PVEC( I,J,K+3 ))
     *                 + TAZ4 * ( PVEC( I,J,K-4 )
     *                 +         PVEC( I,J,K+4 ))
     *                 + TAO  *   PVEC( I,J,K )

       ENDDO
       ENDDO
       ENDDO

      DO 10 J = 1, JY
        DO 20 I = 1, JX

        APVEC( I,J,1 ) = TAZ1 * PVEC( I,J,2 )
     *                 + TAZ2 * PVEC( I,J,3 )
     *                 + TAZ3 * PVEC( I,J,4 )
     *                 + TAZ4 * PVEC( I,J,5 )
     *                 + TAO  * PVEC( I,J,1 )

        APVEC( I,J,2 ) = TAZ1 * ( PVEC( I,J,1 )
     *                 +         PVEC( I,J,3 ))
     *                 + TAZ2 * PVEC( I,J,4 )
     *                 + TAZ3 * PVEC( I,J,5 )
     *                 + TAZ4 * PVEC( I,J,6 )
     *                 + TAO  * PVEC( I,J,2 )

        APVEC( I,J,3 ) = TAZ1 * ( PVEC( I,J,2 )
     *                 +         PVEC( I,J,4 ))
     *                 + TAZ2 * ( PVEC( I,J,1 )
     *                 +         PVEC( I,J,5 ))
     *                 + TAZ3 * PVEC( I,J,6 )
     *                 + TAZ4 * PVEC( I,J,7 )
     *                 + TAO  * PVEC( I,J,3 )

        APVEC( I,J,4 ) = TAZ1 * ( PVEC( I,J,3 )
     *                 +         PVEC( I,J,5 ))
     *                 + TAZ2 * ( PVEC( I,J,2 )
     *                 +         PVEC( I,J,6 ))
     *                 + TAZ3 * ( PVEC( I,J,1 )
     *                 +         PVEC( I,J,7 ))
     *                 + TAZ4 * PVEC( I,J,8 )
     *                 + TAO  * PVEC( I,J,4 )

        APVEC( I,J,JZ-3 ) = TAZ1 * ( PVEC( I,J,JZ-4 )
     *                    +         PVEC( I,J,JZ-2 ))
     *                    + TAZ2 * ( PVEC( I,J,JZ-5 )
     *                    +         PVEC( I,J,JZ-1 ))
     *                    + TAZ3 * ( PVEC( I,J,JZ-6 )
     *                    +         PVEC( I,J,JZ ))
     *                    + TAZ4 * PVEC( I,J,JZ-7 )
     *                    + TAO  * PVEC( I,J,JZ-3 )
        
        APVEC( I,J,JZ-2 ) = TAZ1 * ( PVEC( I,J,JZ-3 )
     *                    +         PVEC( I,J,JZ-1 ))
     *                    + TAZ2 * ( PVEC( I,J,JZ-4 )
     *                    +         PVEC( I,J,JZ ))
     *                    + TAZ3 * PVEC( I,J,JZ-5 )
     *                    + TAZ4 * PVEC( I,J,JZ-6 )
     *                    + TAO  * PVEC( I,J,JZ-2 )

        APVEC( I,J,JZ-1 ) = TAZ1 * ( PVEC( I,J,JZ-2 )
     *                    +         PVEC( I,J,JZ ))
     *                    + TAZ2 * PVEC( I,J,JZ-3 )
     *                    + TAZ3 * PVEC( I,J,JZ-4 )
     *                    + TAZ4 * PVEC( I,J,JZ-5 )
     *                    + TAO  * PVEC( I,J,JZ-1 ) 

        APVEC( I,J,JZ )   = TAZ1 * PVEC( I,J,JZ-1 )
     *                    + TAZ2 * PVEC( I,J,JZ-2 )
     *                    + TAZ3 * PVEC( I,J,JZ-3 )
     *                    + TAZ4 * PVEC( I,J,JZ-4 )
     *                    + TAO  * PVEC( I,J,JZ ) 

20      CONTINUE
10    CONTINUE

C----- Y-AXIS -----
      
      DO K = 1, JZ
      DO J = 5, JY-4
      DO I = 1, JX

        APVEC( I,J,K ) = TAY1 * ( PVEC( I,J-1,K )
     *                 +         PVEC( I,J+1,K ))
     *                 + TAY2 * ( PVEC( I,J-2,K )
     *                 +         PVEC( I,J+2,K ))
     *                 + TAY3 * ( PVEC( I,J-3,K )
     *                 +         PVEC( I,J+3,K ))
     *                 + TAY4 * ( PVEC( I,J-4,K )
     *                 +         PVEC( I,J+4,K ))
     *                 +         APVEC( I,J,K )

      ENDDO
      ENDDO
      ENDDO

      DO 40 K = 1, JZ
        DO 50 I = 1, JX

        APVEC( I,1,K ) = TAY1 * PVEC( I,2,K )
     *                 + TAY2 * PVEC( I,3,K )
     *                 + TAY3 * PVEC( I,4,K )
     *                 + TAY4 * PVEC( I,5,K )
     *                 +        APVEC( I,1,K )

        APVEC( I,2,K ) = TAY1 * ( PVEC( I,1,K )
     *                 +         PVEC( I,3,K ))
     *                 + TAY2 * PVEC( I,4,K )
     *                 + TAY3 * PVEC( I,5,K )
     *                 + TAY4 * PVEC( I,6,K )
     *                 +        APVEC( I,2,K )

        APVEC( I,3,K ) = TAY1 * ( PVEC( I,2,K )
     *                 +         PVEC( I,4,K ))
     *                 + TAY2 * ( PVEC( I,1,K )
     *                 +         PVEC( I,5,K ))
     *                 + TAY3 * PVEC( I,6,K )
     *                 + TAY4 * PVEC( I,7,K )
     *                 +        APVEC( I,3,K )

        APVEC( I,4,K ) = TAY1 * ( PVEC( I,3,K )
     *                 +         PVEC( I,5,K ))
     *                 + TAY2 * ( PVEC( I,2,K )
     *                 +         PVEC( I,6,K ))
     *                 + TAY3 * ( PVEC( I,1,K )
     *                 +         PVEC( I,7,K ))
     *                 + TAY4 * PVEC( I,8,K )
     *                 +        APVEC( I,4,K )

        APVEC( I,JY-3,K ) = TAY1 * ( PVEC( I,JY-4,K )
     *                    +         PVEC( I,JY-2,K ))
     *                    + TAY2 * ( PVEC( I,JY-5,K )
     *                    +         PVEC( I,JY-1,K ))
     *                    + TAY3 * ( PVEC( I,JY-6,K )
     *                    +         PVEC( I,JY,K ))
     *                    + TAY4 * PVEC( I,JY-7,K )
     *                    +        APVEC( I,JY-3,K )
        
        APVEC( I,JY-2,K ) = TAY1 * ( PVEC( I,JY-3,K )
     *                    +         PVEC( I,JY-1,K ))
     *                    + TAY2 * ( PVEC( I,JY-4,K )
     *                    +         PVEC( I,JY,K ))
     *                    + TAY3 * PVEC( I,JY-5,K )
     *                    + TAY4 * PVEC( I,JY-6,K )
     *                    +        APVEC( I,JY-2,K )

        APVEC( I,JY-1,K ) = TAY1 * ( PVEC( I,JY-2,K )
     *                    +         PVEC( I,JY,K ))
     *                    + TAY2 * PVEC( I,JY-3,K )
     *                    + TAY3 * PVEC( I,JY-4,K )
     *                    + TAY4 * PVEC( I,JY-5,K )
     *                    +        APVEC( I,JY-1,K ) 

        APVEC( I,JY,K )   = TAY1 * PVEC( I,JY-1,K )
     *                    + TAY2 * PVEC( I,JY-2,K )
     *                    + TAY3 * PVEC( I,JY-3,K )
     *                    + TAY4 * PVEC( I,JY-4,K )
     *                    +        APVEC( I,JY,K ) 

50      CONTINUE
40    CONTINUE

C----- X-AXIS -----
      
      DO K = 1, JZ
      DO J = 1, JY
      DO I = 5, JX-4

        APVEC( I,J,K ) = TAX1 * ( PVEC( I-1,J,K )
     *                 +         PVEC( I+1,J,K ))
     *                 + TAX2 * ( PVEC( I-2,J,K )
     *                 +         PVEC( I+2,J,K ))
     *                 + TAX3 * ( PVEC( I-3,J,K )
     *                 +         PVEC( I+3,J,K ))
     *                 + TAX4 * ( PVEC( I-4,J,K )
     *                 +         PVEC( I+4,J,K ))
     *                 +         APVEC( I,J,K )

      ENDDO
      ENDDO
      ENDDO

      DO 70 K = 1, JZ
        DO 80 J = 1, JY

        APVEC( 1,J,K ) = TAX1 * PVEC( 2,J,K )
     *                 + TAX2 * PVEC( 3,J,K )
     *                 + TAX3 * PVEC( 4,J,K )
     *                 + TAX4 * PVEC( 5,J,K )
     *                 +        APVEC( 1,J,K )

        APVEC( 2,J,K ) = TAX1 * ( PVEC( 1,J,K )
     *                 +         PVEC( 3,J,K ))
     *                 + TAX2 * PVEC( 4,J,K )
     *                 + TAX3 * PVEC( 5,J,K )
     *                 + TAX4 * PVEC( 6,J,K )
     *                 +        APVEC( 2,J,K )

        APVEC( 3,J,K ) = TAX1 * ( PVEC( 2,J,K )
     *                 +         PVEC( 4,J,K ))
     *                 + TAX2 * ( PVEC( 1,J,K )
     *                 +         PVEC( 5,J,K ))
     *                 + TAX3 * PVEC( 6,J,K )
     *                 + TAX4 * PVEC( 7,J,K )
     *                 +        APVEC( 3,J,K )

        APVEC( 4,J,K ) = TAX1 * ( PVEC( 3,J,K )
     *                 +         PVEC( 5,J,K ))
     *                 + TAX2 * ( PVEC( 2,J,K )
     *                 +         PVEC( 6,J,K ))
     *                 + TAX3 * ( PVEC( 1,J,K )
     *                 +         PVEC( 7,J,K ))
     *                 + TAX4 * PVEC( 8,J,K )
     *                 +        APVEC( 4,J,K )

        APVEC( JX-3,J,K ) = TAX1 * ( PVEC( JX-4,J,K )
     *                    +         PVEC( JX-2,J,K ))
     *                    + TAX2 * ( PVEC( JX-5,J,K )
     *                    +         PVEC( JX-1,J,K ))
     *                    + TAX3 * ( PVEC( JX-6,J,K )
     *                    +         PVEC( JX,J,K ))
     *                    + TAX4 * PVEC( JX-7,J,K )
     *                    +        APVEC( JX-3,J,K )
        
        APVEC( JX-2,J,K ) = TAX1 * ( PVEC( JX-3,J,K )
     *                    +         PVEC( JX-1,J,K ))
     *                    + TAX2 * ( PVEC( JX-4,J,K )
     *                    +         PVEC( JX,J,K ))
     *                    + TAX3 * PVEC( JX-5,J,K )
     *                    + TAX4 * PVEC( JX-6,J,K )
     *                    +        APVEC( JX-2,J,K )

        APVEC( JX-1,J,K ) = TAX1 * ( PVEC( JX-2,J,K )
     *                    +         PVEC( JX,J,K ))
     *                    + TAX2 * PVEC( JX-3,J,K )
     *                    + TAX3 * PVEC( JX-4,J,K )
     *                    + TAX4 * PVEC( JX-5,J,K )
     *                    +        APVEC( JX-1,J,K ) 

        APVEC( JX,J,K )   = TAX1 * PVEC( JX-1,J,K )
     *                    + TAX2 * PVEC( JX-2,J,K )
     *                    + TAX3 * PVEC( JX-3,J,K )
     *                    + TAX4 * PVEC( JX-4,J,K )
     *                    +        APVEC( JX,J,K ) 

80      CONTINUE
70    CONTINUE

      RETURN
      END


C-------------------------------------------------
C
C   SUBROUTINE BLYP 
C   electronic exchange and correlation  
C
C   exchange energy proposed by Becke
C   Phys. Rev. A, 38, 3098(1988)
C
C   correlation energy proposed by Lee-Yang-Parr
C   Phys. Rev. B, 37, 785(1988)
C
C-------------------------------------------------

      SUBROUTINE BLYP( RRHO,RRHOA,RRHOB,VEXA,VEXB,PEXC,PEXT )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol
      include 'nlocd.i'                         ! Vmol   non-local d
      include 'pc_crr.i'                        ! Vmol   pcc
      include 'BHHpara.i'                       ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NOCC =  1 )
      PARAMETER ( PI   = 3.14159265358979323D0 )
      PARAMETER ( BETA = 0.0042D0 )
C     PARAMETER ( EPSX = 1.D-4 )
C     PARAMETER ( EPSC = 1.D-6 )
      PARAMETER ( EPSX = 1.D-6 )
      PARAMETER ( EPSC = 1.D-8 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI3 / NIDX(1:2),NIDY(1:2),NIDZ(1:2)
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ            
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF

      DIMENSION RHO(  -3:NMAXX/NX+4,-3:NMAXY/NY+4,-3:NMAXZ/NZ+4 )
      DIMENSION RHOA( -3:NMAXX/NX+4,-3:NMAXY/NY+4,-3:NMAXZ/NZ+4 )
      DIMENSION RHOB( -3:NMAXX/NX+4,-3:NMAXY/NY+4,-3:NMAXZ/NZ+4 ) 
      DIMENSION RRHO(    NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RRHOA(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RRHOB(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DRHOAX(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DRHOAY(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DRHOAZ(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DRHOBX(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DRHOBY(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DRHOBZ(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEXA(    NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEXB(    NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION F1A(  -3:NMAXX/NX+4,-3:NMAXY/NY+4,-3:NMAXZ/NZ+4 )
      DIMENSION F1B(  -3:NMAXX/NX+4,-3:NMAXY/NY+4,-3:NMAXZ/NZ+4 )
      DIMENSION G2(   -3:NMAXX/NX+4,-3:NMAXY/NY+4,-3:NMAXZ/NZ+4 )
C     DIMENSION AZ( NMAX*NMAX*4*3)
C     DIMENSION BZ( NMAX*NMAX*4*3)
      DIMENSION AZ( NMAXX*NMAXX*4*3)
      DIMENSION BZ( NMAXX*NMAXX*4*3)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      RHEXR = 1.D0 - RHEX       ! RHEX = mixing ratio of HF exchange

      DV = DX*DY*DZ

C----- COEFFICIENTS -----

      A  = -0.4582D0            ! EXCHANGE
      A1 =  1.0D0               ! RS >= 1.
      A2 =  1.0529D0
      A3 =  0.3334D0
      A4 = -0.1423D0

      AA = 0.04918D0
      BB = 0.132D0
      CC = 0.2533D0
      DD = 0.349D0
      CF = (3.D0/10.D0)*((3.D0*PI*PI)**(2.D0/3.D0)) 
      
C----- initialize boxes for Send and Receive -----
    
      DO K = -3,JZ+4
      DO J = -3,JY+4
      DO I = -3,JX+4
      RHO(I,J,K) =0.D0
      RHOA(I,J,K)=0.D0
      RHOB(I,J,K)=0.D0
      G2(I,J,K)  =0.D0
      F1A(I,J,K) =0.D0
      F1B(I,J,K) =0.D0
      ENDDO
      ENDDO
      ENDDO
      
      IF(MYID .EQ. 41) THEN
        WRITE(*,*) PCDNSB(4,4,4)
      ENDIF

      DO K = 1,JZ
      DO J = 1,JY
      DO I = 1,JX

C----- add partial charge density -------

      IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN
        RHO( I,J,K)=RRHO( I,J,K) + PCDNSA(I,J,K) + PCDNSB(I,J,K) 
      ELSE
        RHO( I,J,K)=RRHO( I,J,K) + 2.D0*PCDNSA(I,J,K)
      ENDIF

      RHOA(I,J,K)=RRHOA(I,J,K) + PCDNSA(I,J,K)
      RHOB(I,J,K)=RRHOB(I,J,K) + PCDNSB(I,J,K)

C     RHO( I,J,K)=RRHO( I,J,K) 
C     RHOA(I,J,K)=RRHOA(I,J,K) 
C     RHOB(I,J,K)=RRHOB(I,J,K) 

      ENDDO
      ENDDO
      ENDDO

      CALL SENDRECV_BUF(RHO ,JX,JY,JZ,4)
      CALL SENDRECV_BUF(RHOA,JX,JY,JZ,4)
      CALL SENDRECV_BUF(RHOB,JX,JY,JZ,4)

C----- CALCULATE EXCHANGE AND CORRELATION -----

      PEXC = 0.D0
      PEXT = 0.D0

      CALL DWF( RHOA,DRHOAX,DRHOAY,DRHOAZ ) 
      
      IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN
      
      CALL DWF( RHOB,DRHOBX,DRHOBY,DRHOBZ ) 
      
      ENDIF

      DO 1 L = 1, JZ
      DO 1 M = 1, JY
      DO 1 N = 1, JX
          
C====== CALCULATION FOR EXCHANGE ======      

C       FOR ALPHA SPIN

        IF(RHOA(N,M,L) .LT. 1.D-10) THEN

        F1A(N,M,L) = 0.D0

        ELSE

        RHOA43 = RHOA(N,M,L)**(4.D0/3.D0)

        DRHOA2 = DRHOAX(N,M,L)*DRHOAX(N,M,L)
     *         + DRHOAY(N,M,L)*DRHOAY(N,M,L)
     *         + DRHOAZ(N,M,L)*DRHOAZ(N,M,L)
        
        DRHOA = DSQRT(DRHOA2)

        XA = DRHOA/RHOA43
        XA2 = XA*XA 
        XA215 = DSQRT(1.D0+XA2)
        ASINHA = DLOG(XA+XA215)
        DENOMA = 1.D0+6.D0*BETA*XA*ASINHA

        F1A(N,M,L) = (BETA/RHOA43)*((3.D0*BETA*XA
     *             * ((XA/XA215)-ASINHA)-1.D0)
     *             / (DENOMA**2.D0))

        ENDIF

C       FOR BETA SPIN

        IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN

         IF(RHOB(N,M,L) .LT. 1.D-10) THEN

        F1B(N,M,L) = 0.D0

         ELSE

        RHOB43 = RHOB(N,M,L)**(4.D0/3.D0)

        DRHOB2 = DRHOBX(N,M,L)*DRHOBX(N,M,L)
     *         + DRHOBY(N,M,L)*DRHOBY(N,M,L)
     *         + DRHOBZ(N,M,L)*DRHOBZ(N,M,L)
        
         DRHOB = DSQRT(DRHOB2)

        XB = DRHOB/RHOB43
        XB2 = XB*XB 
        XB215 = DSQRT(1.D0+XB2)
        ASINHB = DLOG(XB+XB215)
        DENOMB = 1.D0+6.D0*BETA*XB*ASINHB

        F1B(N,M,L) = (BETA/RHOB43)*((3.D0*BETA*XB
     *             * ((XB/XB215)-ASINHB)-1.D0)
     *             / (DENOMB**2.D0))

        ENDIF

        ENDIF

C       CALCULATION FOR CORRELATION       

        RHOI13 = RHO(N,M,L)**(-1.D0/3.D0) 
        RHOI53 = RHO(N,M,L)**(-5.D0/3.D0) 
        RHO2 = RHO(N,M,L)*RHO(N,M,L)
        RHOA2 = RHOA(N,M,L)*RHOA(N,M,L)
        P1 = 1.D0 + DD*RHOI13
        P2 = -CC*RHOI13

        IF(RHO2 .LT. 1.D-20) THEN

        GAMMA = 0.D0
        
        ELSEIF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' ) THEN
       
        GAMMA = 2.D0*(1.D0-((2.D0*RHOA2)/RHO2))
        
        ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN
         
        RHOB2 = RHOB(N,M,L)*RHOB(N,M,L)
        GAMMA = 2.D0*(1.D0-((RHOA2+RHOB2)/RHO2))
        
        ENDIF

        F2 = GAMMA/P1

        G2(N,M,L) = F2*RHOI53*DEXP(P2)

1     CONTINUE

      CALL SENDRECV_BUF(G2 ,JX,JY,JZ,4)
      CALL SENDRECV_BUF(F1A,JX,JY,JZ,4)
      CALL SENDRECV_BUF(F1B,JX,JY,JZ,4)

C--------------------------------------

      DO 10 L = 1,JZ
      DO 10 M = 1,JY
      DO 10 N = 1,JX
          
        RSA= (3.D0/(4.D0*PI))**(1.D0/3.D0)
     *     * (1.D0/(2.D0*RHOA( N,M,L )))**(1.D0/3.D0)

       IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN
         
        RSB= (3.D0/(4.D0*PI))**(1.D0/3.D0)
     *     * (1.D0/(2.D0*RHOB( N,M,L )))**(1.D0/3.D0)

       ENDIF

        IF ( ( RHO(N,M,L) .LT. EPSC )
     *  .OR. ((N.LE.4   ).AND.(MX.EQ.0   ))
     *  .OR. ((N.GE.JX-3).AND.(MX.EQ.NX-1))
     *  .OR. ((M.LE.4   ).AND.(MY.EQ.0   ))
     *  .OR. ((M.GE.JY-3).AND.(MY.EQ.NY-1))
     *  .OR. ((L.LE.4   ).AND.(MZ.EQ.0   ))
     *  .OR. ((L.GE.JZ-3).AND.(MZ.EQ.NZ-1)) ) THEN

            DENA = ( A1+A2*RSA**(1.D0/2.D0)+A3*RSA ) 

            VEXA( N,M,L ) = A/RSA + A4/DENA
     *                    - (RSA/3.D0)*( -A/RSA**2.D0 
     *                    - A4*( (A2/2.D0)*(RSA)**(-1.D0/2.D0)+A3 )
     *                    / DENA**2.D0)
          
            PVEXA  = A/RSA + A4/DENA

         IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN
         
            DENB = ( A1+A2*RSB**(1.D0/2.D0)+A3*RSB ) 

            VEXB( N,M,L ) = A/RSB + A4/DENB
     *                    - (RSB/3.D0)*( -A/RSB**2.D0 
     *                    - A4*( (A2/2.D0)*(RSB)**(-1.D0/2.D0)+A3 )
     *                    / DENB**2.D0)

            PVEXB  = A/RSB + A4/DENB

          
            PEXC = PEXC + ( PVEXA*RHOA( N,M,L )                        ! total exchange and correlation
     *                  +   PVEXB*RHOB( N,M,L ) )*DV             

         ELSEIF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' ) THEN

            PEXC = PEXC + 2.D0*( PVEXA*RHOA( N,M,L ))*DV               ! total exchange and correlation

         ENDIF

        ELSE

C------- FIRST DRIVATIVE -------

        DRHOX = DFX(RHO(N+1,M,L),RHO(N-1,M,L),
     *              RHO(N+2,M,L),RHO(N-2,M,L),
     *              RHO(N+3,M,L),RHO(N-3,M,L),
     *              RHO(N+4,M,L),RHO(N-4,M,L))

        DRHOY = DFY(RHO(N,M+1,L),RHO(N,M-1,L),
     *              RHO(N,M+2,L),RHO(N,M-2,L),
     *              RHO(N,M+3,L),RHO(N,M-3,L),
     *              RHO(N,M+4,L),RHO(N,M-4,L))

        DRHOZ = DFZ(RHO(N,M,L+1),RHO(N,M,L-1),
     *              RHO(N,M,L+2),RHO(N,M,L-2),
     *              RHO(N,M,L+3),RHO(N,M,L-3),
     *              RHO(N,M,L+4),RHO(N,M,L-4))

        DRHO2 = DRHOX*DRHOX
     *        + DRHOY*DRHOY
     *        + DRHOZ*DRHOZ
        
C       DRHO = DSQRT(DRHO2)

        DRHOA2 = DRHOAX(N,M,L)*DRHOAX(N,M,L)
     *         + DRHOAY(N,M,L)*DRHOAY(N,M,L)
     *         + DRHOAZ(N,M,L)*DRHOAZ(N,M,L)
        
        DRHOA = DSQRT(DRHOA2)

        DG2X = DFX(G2(N+1,M,L),G2(N-1,M,L),
     *             G2(N+2,M,L),G2(N-2,M,L),
     *             G2(N+3,M,L),G2(N-3,M,L),
     *             G2(N+4,M,L),G2(N-4,M,L))

        DG2Y = DFY(G2(N,M+1,L),G2(N,M-1,L),
     *             G2(N,M+2,L),G2(N,M-2,L),
     *             G2(N,M+3,L),G2(N,M-3,L),
     *             G2(N,M+4,L),G2(N,M-4,L))

        DG2Z = DFZ(G2(N,M,L+1),G2(N,M,L-1),
     *             G2(N,M,L+2),G2(N,M,L-2),
     *             G2(N,M,L+3),G2(N,M,L-3),
     *             G2(N,M,L+4),G2(N,M,L-4))

C------- SECOND DERIVATIVE --------

         SDRHO = SDF(RHO(N,M,L),
     *               RHO(N+1,M,L),RHO(N-1,M,L),
     *               RHO(N+2,M,L),RHO(N-2,M,L),
     *               RHO(N+3,M,L),RHO(N-3,M,L),
     *               RHO(N+4,M,L),RHO(N-4,M,L),
     *               RHO(N,M+1,L),RHO(N,M-1,L),
     *               RHO(N,M+2,L),RHO(N,M-2,L),
     *               RHO(N,M+3,L),RHO(N,M-3,L),
     *               RHO(N,M+4,L),RHO(N,M-4,L),
     *               RHO(N,M,L+1),RHO(N,M,L-1),
     *               RHO(N,M,L+2),RHO(N,M,L-2),
     *               RHO(N,M,L+3),RHO(N,M,L-3),
     *               RHO(N,M,L+4),RHO(N,M,L-4))

         SDRHOA = SDF(RHOA(N,M,L),
     *               RHOA(N+1,M,L),RHOA(N-1,M,L),
     *               RHOA(N+2,M,L),RHOA(N-2,M,L),
     *               RHOA(N+3,M,L),RHOA(N-3,M,L),
     *               RHOA(N+4,M,L),RHOA(N-4,M,L),
     *               RHOA(N,M+1,L),RHOA(N,M-1,L),
     *               RHOA(N,M+2,L),RHOA(N,M-2,L),
     *               RHOA(N,M+3,L),RHOA(N,M-3,L),
     *               RHOA(N,M+4,L),RHOA(N,M-4,L),
     *               RHOA(N,M,L+1),RHOA(N,M,L-1),
     *               RHOA(N,M,L+2),RHOA(N,M,L-2),
     *               RHOA(N,M,L+3),RHOA(N,M,L-3),
     *               RHOA(N,M,L+4),RHOA(N,M,L-4))

         SDG2 = SDF(G2(N,M,L),
     *              G2(N+1,M,L),G2(N-1,M,L),
     *              G2(N+2,M,L),G2(N-2,M,L),
     *              G2(N+3,M,L),G2(N-3,M,L),
     *              G2(N+4,M,L),G2(N-4,M,L),
     *              G2(N,M+1,L),G2(N,M-1,L),
     *              G2(N,M+2,L),G2(N,M-2,L),
     *              G2(N,M+3,L),G2(N,M-3,L),
     *              G2(N,M+4,L),G2(N,M-4,L),
     *              G2(N,M,L+1),G2(N,M,L-1),
     *              G2(N,M,L+2),G2(N,M,L-2),
     *              G2(N,M,L+3),G2(N,M,L-3),
     *              G2(N,M,L+4),G2(N,M,L-4))

       IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN

        DRHOB2 = DRHOBX(N,M,L)*DRHOBX(N,M,L)
     *         + DRHOBY(N,M,L)*DRHOBY(N,M,L)
     *         + DRHOBZ(N,M,L)*DRHOBZ(N,M,L)
        
         DRHOB = DSQRT(DRHOB2)

         SDRHOB = SDF(RHOB(N,M,L),
     *                RHOB(N+1,M,L),RHOB(N-1,M,L),
     *                RHOB(N+2,M,L),RHOB(N-2,M,L),
     *                RHOB(N+3,M,L),RHOB(N-3,M,L),
     *                RHOB(N+4,M,L),RHOB(N-4,M,L),
     *                RHOB(N,M+1,L),RHOB(N,M-1,L),
     *                RHOB(N,M+2,L),RHOB(N,M-2,L),
     *                RHOB(N,M+3,L),RHOB(N,M-3,L),
     *                RHOB(N,M+4,L),RHOB(N,M-4,L),
     *                RHOB(N,M,L+1),RHOB(N,M,L-1),
     *                RHOB(N,M,L+2),RHOB(N,M,L-2),
     *                RHOB(N,M,L+3),RHOB(N,M,L-3),
     *                RHOB(N,M,L+4),RHOB(N,M,L-4))

       ENDIF

C------- exchange for alpha spin -------

        IF( RHOA(N,M,L) .LT. EPSX ) THEN                          ! AND (RSA .GT. 1)

          VEXA( N,M,L ) = (4.D0/3.D0)*(A/RSA)                     ! EXCHANGE (LDA)
          VEXA( N,M,L ) = RHEXR*VEXA( N,M,L )                     ! treatment for BHandHLYP

C------- total exchange -------
      
         IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' ) THEN

C           PEXC = PEXC + 2.D0*(A/RSA)*RHOA( N,M,L )*DV
            PEXC = PEXC + RHEXR*2.D0*(A/RSA)*RHOA( N,M,L )*DV      ! treatment for BHandHLYP

         ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN
          
C           PEXC = PEXC + (A/RSA)*RHOA( N,M,L )*DV
            PEXC = PEXC + RHEXR*(A/RSA)*RHOA( N,M,L )*DV           ! treatment for BHandHLYP
       
         ENDIF

        ELSE

C------- FIRST DERIVATIVE -------

        DF1AX = DFX(F1A(N+1,M,L),F1A(N-1,M,L),
     *              F1A(N+2,M,L),F1A(N-2,M,L),
     *              F1A(N+3,M,L),F1A(N-3,M,L),
     *              F1A(N+4,M,L),F1A(N-4,M,L))

        DF1AY = DFY(F1A(N,M+1,L),F1A(N,M-1,L),
     *              F1A(N,M+2,L),F1A(N,M-2,L),
     *              F1A(N,M+3,L),F1A(N,M-3,L),
     *              F1A(N,M+4,L),F1A(N,M-4,L))

        DF1AZ = DFZ(F1A(N,M,L+1),F1A(N,M,L-1),
     *              F1A(N,M,L+2),F1A(N,M,L-2),
     *              F1A(N,M,L+3),F1A(N,M,L-3),
     *              F1A(N,M,L+4),F1A(N,M,L-4))

        RHOA13 = RHOA(N,M,L)**(1.D0/3.D0)
        RHOA43 = RHOA13*RHOA(N,M,L)

        XA = DRHOA/RHOA43
        XA2 = XA*XA 
        XA215 = DSQRT(1.D0+XA2)
        ASINHA = DLOG(XA+XA215)
        DENOMA = 1.D0+6.D0*BETA*XA*ASINHA

          VEXA( N,M,L ) = (4.D0/3.D0)*(A/RSA)                     ! EXCHANGE (LDA)
     *                  - BETA*(4.D0/3.D0)*RHOA13                 ! GRADIENT CORRECTION FOR EXCHANGE (BECKE)
     *                  * ((XA2/DENOMA)
     *                  + 2.D0*XA2*(RHOA43/BETA)*F1A(N,M,L))
     *                  - 2.D0*(F1A(N,M,L)*SDRHOA
     *                  + DRHOAX(N,M,L)*DF1AX
     *                  + DRHOAY(N,M,L)*DF1AY
     *                  + DRHOAZ(N,M,L)*DF1AZ)
          VEXA( N,M,L ) = RHEXR*VEXA( N,M,L )                     ! treatment for BHandHLYP

C------- total exchange -------
      
           IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' ) THEN

           PEXC = PEXC +RHEXR*2.D0*((A/RSA)*RHOA( N,M,L )        ! treatment for BHandHLYP
     *         - BETA*RHOA43*(XA2/DENOMA))*DV                    ! gradient correction for exchange

           ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN

           PEXC = PEXC +RHEXR*((A/RSA)*RHOA( N,M,L )             ! treatment for BHandHLYP
     *         - BETA*RHOA43*(XA2/DENOMA))*DV                    ! gradient correction for exchange

           ENDIF

        ENDIF

       IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN

C------- exchange for beta spin -------

        IF( RHOB(N,M,L) .LT. EPSX ) THEN                          ! AND (RSB .GT. 1)

          VEXB( N,M,L ) = (4.D0/3.D0)*(A/RSB)                     ! EXCHANGE (LDA)
          VEXB( N,M,L ) = RHEXR*VEXB( N,M,L )                     ! treatment for BHandHLYP

C------- total exchange -------
      
C         PEXC = PEXC + (A/RSB)*RHOB( N,M,L )*DV
          PEXC = PEXC + RHEXR*(A/RSB)*RHOB( N,M,L )*DV            ! treatment for BHandHLYP

        ELSE

C------- FIRST DERIVATIVE -------

        DF1BX = DFX(F1B(N+1,M,L),F1B(N-1,M,L),
     *              F1B(N+2,M,L),F1B(N-2,M,L),
     *              F1B(N+3,M,L),F1B(N-3,M,L),
     *              F1B(N+4,M,L),F1B(N-4,M,L))

        DF1BY = DFY(F1B(N,M+1,L),F1B(N,M-1,L),
     *              F1B(N,M+2,L),F1B(N,M-2,L),
     *              F1B(N,M+3,L),F1B(N,M-3,L),
     *              F1B(N,M+4,L),F1B(N,M-4,L))

        DF1BZ = DFZ(F1B(N,M,L+1),F1B(N,M,L-1),
     *              F1B(N,M,L+2),F1B(N,M,L-2),
     *              F1B(N,M,L+3),F1B(N,M,L-3),
     *              F1B(N,M,L+4),F1B(N,M,L-4))

        RHOB13 = RHOB(N,M,L)**(1.D0/3.D0)
        RHOB43 = RHOB13*RHOB(N,M,L)

        XB  = DRHOB/RHOB43
        XB2 = XB*XB 
        XB215  = DSQRT(1.D0+XB2)
        ASINHB = DLOG(XB+XB215)
        DENOMB = 1.D0+6.D0*BETA*XB*ASINHB

          VEXB( N,M,L ) = (4.D0/3.D0)*(A/RSB)                     ! EXCHANGE (LDA)
     *                  - BETA*(4.D0/3.D0)*RHOB13                 ! GRADIENT CORRECTION FOR EXCHANGE (BECKE)
     *                  * ((XB2/DENOMB)
     *                  + 2.D0*XB2*(RHOB43/BETA)*F1B(N,M,L))
     *                  - 2.D0*(F1B(N,M,L)*SDRHOB
     *                  + DRHOBX(N,M,L)*DF1BX
     *                  + DRHOBY(N,M,L)*DF1BY
     *                  + DRHOBZ(N,M,L)*DF1BZ)
          VEXB( N,M,L ) = RHEXR*VEXB( N,M,L )                     ! treatment for BHandHLYP

C------- total exchange -------
      
          PEXC = PEXC +RHEXR*((A/RSB)*RHOB( N,M,L )               ! treatment for BHandHLYP
     *         - BETA*RHOB43*(XB2/DENOMB))*DV                     ! gradient correction for exchange

        ENDIF

       ENDIF

C------ calculation for correlation ------

         RHOA53 = RHOA(N,M,L)**(5.D0/3.D0)
         RHOA83 = RHOA(N,M,L)**(8.D0/3.D0)
         RHOI = 1.D0/RHO(N,M,L) 
         RHOI13 = RHO(N,M,L)**(-1.D0/3.D0) 
         RHOI43 = RHO(N,M,L)**(-4.D0/3.D0) 
         RHOI53 = RHO(N,M,L)**(-5.D0/3.D0) 
         RHO2 = RHO(N,M,L)*RHO(N,M,L)
         RHO3 = RHO(N,M,L)*RHO2
         RHOA2 = RHOA(N,M,L)*RHOA(N,M,L)

         P1 = 1.D0 + DD*RHOI13
         P2 = -CC*RHOI13
         P3 = RHO(N,M,L)*SDG2 
         P4 = 4.D0*(DG2X*DRHOX+DG2Y*DRHOY+DG2Z*DRHOZ)
         P5 = 4.D0*G2(N,M,L)*SDRHO 
         P6 = RHO(N,M,L)*SDRHO-DRHO2

         PA1 = 3.D0*RHOA(N,M,L)*SDG2
         PA2 = 4.D0*(DRHOAX(N,M,L)*DG2X
     *       + DRHOAY(N,M,L)*DG2Y+DRHOAZ(N,M,L)*DG2Z)
         PA3 = 4.D0*G2(N,M,L)*SDRHOA

       IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN

         RHOB53 = RHOB(N,M,L)**(5.D0/3.D0)
         RHOB83 = RHOB(N,M,L)**(8.D0/3.D0)
           RHOB2 = RHOB(N,M,L)*RHOB(N,M,L)
          P7 = 3.D0*(RHOA(N,M,L)*SDRHOA+RHOB(N,M,L)*SDRHOB)
     *      + (DRHOA2+DRHOB2)
          PB1 = 3.D0*RHOB(N,M,L)*SDG2
          PB2 = 4.D0*(DRHOBX(N,M,L)*DG2X
     *       + DRHOBY(N,M,L)*DG2Y+DRHOBZ(N,M,L)*DG2Z)
          PB3 = 4.D0*G2(N,M,L)*SDRHOB
         GAMMA = 2.D0*(1.D0-((RHOA2+RHOB2)/RHO2))
         GAMMAA = 4.D0*(((RHOA2+RHOB2)-RHOA(N,M,L)*RHO(N,M,L))/RHO3)
         GAMMAB = 4.D0*(((RHOA2+RHOB2)-RHOB(N,M,L)*RHO(N,M,L))/RHO3)
         F2 = GAMMA/P1
         F2A = (GAMMAA*P1+(1.D0/3.D0)*DD*GAMMA*RHOI43)/(P1*P1)
         F2B = (GAMMAB*P1+(1.D0/3.D0)*DD*GAMMA*RHOI43)/(P1*P1)
         G2A = RHOI53*DEXP(P2)*(F2A-(5.D0/3.D0)*RHOI*F2
     *       + (1.D0/3.D0)*CC*RHOI43*F2)
         G2B = RHOI53*DEXP(P2)*(F2B-(5.D0/3.D0)*RHOI*F2
     *       + (1.D0/3.D0)*CC*RHOI43*F2)

       ELSEIF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' ) THEN       

          P7 = 6.D0*(RHOA(N,M,L)*SDRHOA)
     *      + (2.D0*DRHOA2)
         GAMMA = 2.D0*(1.D0-((2.D0*RHOA2)/RHO2))
         GAMMAA = 4.D0*((2.D0*(RHOA2)-RHOA(N,M,L)*RHO(N,M,L))/RHO3)
          F2 = GAMMA/P1
         F2A = (GAMMAA*P1+(1.D0/3.D0)*DD*GAMMA*RHOI43)/(P1*P1)
         G2A = RHOI53*DEXP(P2)*(F2A-(5.D0/3.D0)*RHOI*F2
     *       + (1.D0/3.D0)*CC*RHOI43*F2)

       ENDIF

C------ calculate exchange-correlation potential ------

       IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN

          VEXA( N,M,L ) = VEXA( N,M,L )                           ! EXCHANGE
     *                  - AA*(F2A*RHO(N,M,L)+F2)                  ! CORRELATION ( LEE-YANG-PARR )
     *                  - (2.D0**(5.D0/3.D0))*AA*BB*CF 
C    *                  * (G2A*(RHOA83*RHOB83)                    ! Error
     *                  * (G2A*(RHOA83+RHOB83)
     *                  + (8.D0/3.D0)*G2(N,M,L)*RHOA53)
     *                  - ((AA*BB)/4.D0)*(P3+P4+P5+G2A*P6)
     *                  - ((AA*BB)/36.D0)*(PA1+PA2+PA3+G2A*P7)

          VEXB( N,M,L ) = VEXB( N,M,L )                           ! EXCHANGE
     *                  - AA*(F2B*RHO(N,M,L)+F2)                  ! CORRELATION ( LEE-YANG-PARR )
     *                  - (2.D0**(5.D0/3.D0))*AA*BB*CF 
     *                  * (G2B*(RHOA83+RHOB83)
     *                  + (8.D0/3.D0)*G2(N,M,L)*RHOB53)
     *                  - ((AA*BB)/4.D0)*(P3+P4+P5+G2B*P6)
     *                  - ((AA*BB)/36.D0)*(PB1+PB2+PB3+G2B*P7)

          TW = (1.D0/8.D0)*((DRHO2/RHO(N,M,L))-SDRHO)
          TWA = (1.D0/8.D0)*((DRHOA2/RHOA(N,M,L))-SDRHOA)
          TWB = (1.D0/8.D0)*((DRHOB2/RHOB(N,M,L))-SDRHOB)
 
          Q1 = 2.D0**(2.D0/3.D0)*CF*(RHOA83+RHOB83)
          Q2 = -RHO(N,M,L)*TW
     *       + (1.D0/9.D0)*(RHOA(N,M,L)*TWA+RHOB(N,M,L)*TWB)
          Q3 = (1.D0/18.D0)*(RHOA(N,M,L)*SDRHOA+RHOB(N,M,L)*SDRHOB)

          PEXC = PEXC - (AA*(F2*RHO(N,M,L)                        ! CORRELATION ( LEE-YANG-PARR )
     *         + 2.D0*BB*G2(N,M,L)*(Q1+Q2+Q3)))*DV

       ELSEIF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' ) THEN

          VEXA( N,M,L ) = VEXA( N,M,L )                           ! EXCHANGE
     *                  - AA*(F2A*RHO(N,M,L)+F2)                  ! CORRELATION ( LEE-YANG-PARR )
     *                  - (2.D0**(5.D0/3.D0))*AA*BB*CF 
     *                  * (G2A*2.D0*(RHOA83)
     *                  + (8.D0/3.D0)*G2(N,M,L)*RHOA53)
     *                  - ((AA*BB)/4.D0)*(P3+P4+P5+G2A*P6)
     *                  - ((AA*BB)/36.D0)*(PA1+PA2+PA3+G2A*P7)

           TW = (1.D0/8.D0)*((DRHO2/RHO(N,M,L))-SDRHO)
           TWA = (1.D0/8.D0)*((DRHOA2/RHOA(N,M,L))-SDRHOA)

           Q1 = 2.D0**(2.D0/3.D0)*CF*(2.D0*RHOA83)
           Q2 = -RHO(N,M,L)*TW
     *       + (1.D0/9.D0)*(2.D0*RHOA(N,M,L)*TWA)
           Q3 = (1.D0/18.D0)*(2.D0*RHOA(N,M,L)*SDRHOA)

          PEXC = PEXC - (AA*(F2*RHO(N,M,L)                        ! CORRELATION ( LEE-YANG-PARR )
     *         + 2.D0*BB*G2(N,M,L)*(Q1+Q2+Q3)))*DV

       ENDIF

      ENDIF

C----- external exchange and correlation -----
    
      IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN
C      PEXT = PEXT + ( RHOA(N,M,L)*VEXA(N,M,L)
C    *      +          RHOB(N,M,L)*VEXB(N,M,L) )*DV
       PEXT = PEXT + ( (RHOA(N,M,L)-PCDNSA(N,M,L))*VEXA(N,M,L)
     *      +          (RHOB(N,M,L)-PCDNSB(N,M,L))*VEXB(N,M,L) )*DV
      ELSEIF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' ) THEN
       PEXT = PEXT + (2.D0*(RHOA(N,M,L)-PCDNSA(N,M,L))*VEXA(N,M,L))*DV
      ENDIF

C---------------------------------------------

10    CONTINUE

      CALL MPI_REDUCE(PEXC,PPEXC,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_REDUCE(PEXT,PPEXT,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)

      PEXC = PPEXC
      PEXT = PPEXT
    
      CALL MPI_BCAST(PEXC,1,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(PEXT,1,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)

      IF(MYID.EQ.0) THEN
      write(*,*) 'pexc =',pexc
      write(*,*) 'pext =',pext
      ENDIF

      RETURN
      END

C---------------------------------------------
C   SUBROUTINE DERIVATIVES OF WAVEFUNCTION 
C--------------------------------------------

      SUBROUTINE DWF( RWFN,DWFX,DWFY,DWFZ ) 
      
      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

      DIMENSION RWFN( -3:NMAXX/NX+4,-3:NMAXY/NY+4,-3:NMAXZ/NZ+4 )
      DIMENSION DWFX( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFY( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFZ( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
 
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ

      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)

      A1  = 449.D0 / 360.D0
      A2  = -11.D0 /  72.D0
      A3  =  11.D0 / 504.D0
      A4  =  -1.D0 / 560.D0

      AX1 = A1 / ( 2.D0*DX )
      AX2 = A2 / ( 2.D0*DX )
      AX3 = A3 / ( 2.D0*DX )
      AX4 = A4 / ( 2.D0*DX )

      AY1 = A1 / ( 2.D0*DY )
      AY2 = A2 / ( 2.D0*DY )
      AY3 = A3 / ( 2.D0*DY )
      AY4 = A4 / ( 2.D0*DY )

      AZ1 = A1 / ( 2.D0*DZ )
      AZ2 = A2 / ( 2.D0*DZ )
      AZ3 = A3 / ( 2.D0*DZ )
      AZ4 = A4 / ( 2.D0*DZ )

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX 
        DWFZ( I,J,K ) = AZ1 * ( RWFN( I,J,K+1 ) - RWFN( I,J,K-1 ))
     *                + AZ2 * ( RWFN( I,J,K+2 ) - RWFN( I,J,K-2 ))
     *                + AZ3 * ( RWFN( I,J,K+3 ) - RWFN( I,J,K-3 ))
     *                + AZ4 * ( RWFN( I,J,K+4 ) - RWFN( I,J,K-4 ))

        DWFY( I,J,K ) = AY1 * ( RWFN( I,J+1,K ) - RWFN( I,J-1,K ))
     *                + AY2 * ( RWFN( I,J+2,K ) - RWFN( I,J-2,K ))
     *                + AY3 * ( RWFN( I,J+3,K ) - RWFN( I,J-3,K ))
     *                + AY4 * ( RWFN( I,J+4,K ) - RWFN( I,J-4,K ))

        DWFX( I,J,K ) = AX1 * ( RWFN( I+1,J,K ) - RWFN( I-1,J,K ))
     *                + AX2 * ( RWFN( I+2,J,K ) - RWFN( I-2,J,K ))
     *                + AX3 * ( RWFN( I+3,J,K ) - RWFN( I-3,J,K ))
     *                + AX4 * ( RWFN( I+4,J,K ) - RWFN( I-4,J,K ))

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END

C---------------------------------------------
C   SUBROUTINE DERIVATIVES OF WAVEFUNCTION 
C--------------------------------------------

      SUBROUTINE DWF2( RWFN,DWFX,DWFY,DWFZ ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

      PARAMETER ( NPTS  =  7 )    ! number of sample points for interpolation ( It must be an odd number! ) 
      PARAMETER ( NPTS2 =  3 )    ! NPTS/2(INTEGER)

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI3 / NIDX(1:2),NIDY(1:2),NIDZ(1:2)
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ

      DIMENSION RWFN1( 1-NPTS2:NMAXX/NX+NPTS2,1-NPTS2:NMAXY/NY+NPTS2,
     *                 1-NPTS2:NMAXZ/NZ+NPTS2 )
      DIMENSION RWFN( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFX( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFY( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DWFZ( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      DIMENSION Y(  NPTS )
      DIMENSION RX( NPTS )
      DIMENSION RY( NPTS )
      DIMENSION RZ( NPTS )
      DIMENSION COE( NPTS )
      DIMENSION PD( 2 )
      
      DIMENSION A( NPTS2*2*NMAXX*NMAXX/(NX*NY*NZ) )
      DIMENSION B( NPTS2*2*NMAXX*NMAXX/(NX*NY*NZ) )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      DO INC=1,NPTS
        RZ(INC) = DBLE(-NPTS2+INC-1)*DZ
        RY(INC) = DBLE(-NPTS2+INC-1)*DY
        RX(INC) = DBLE(-NPTS2+INC-1)*DX
      ENDDO

      RWFN1(:,:,:) = 0.D0
      RWFN1(1:JX,1:JY,1:JZ) = RWFN(1:JX,1:JY,1:JZ)
      CALL SENDRECV_BUF(RWFN1,JX,JY,JZ,NPTS2)

C----- z-axis -----

      DO I = 1, JX 
      DO J = 1, JY
      DO K = 1, JZ

C----- z-axis -----

        DO INC=1,NPTS
          Y(INC) = RWFN1(I,J,K-NPTS2+INC-1)
        ENDDO

        CALL POLCOE( RZ,Y,NPTS,COE )
        CALL DDPOLY( COE,NPTS,0.D0,PD,2 )              

        DWFZ( I,J,K ) = PD(2)

C----- y-axis -----

        DO INC=1,NPTS
          Y(INC) = RWFN1(I,J-NPTS2+INC-1,K)
        ENDDO

        CALL POLCOE( RY,Y,NPTS,COE )
        CALL DDPOLY( COE,NPTS,0.D0,PD,2 )              

        DWFY( I,J,K ) = PD(2)

C----- x-axis -----

        DO INC=1,NPTS
         Y(INC) = RWFN1(I-NPTS2+INC-1,J,K)
        ENDDO

        CALL POLCOE( RX,Y,NPTS,COE )
        CALL DDPOLY( COE,NPTS,0.D0,PD,2 )              

        DWFX( I,J,K ) = PD(2)

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END

CCC MPI kokokara CCC
C--------------------------------------------------------
C     SUBROUTINE Averaging Electronic Field formed by MM molecules
C--------------------------------------------------------

      SUBROUTINE ELECF( IMD )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
      INTEGER*4 STAT

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
!     include 'atoms.i'
C     include 'charge.i'

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC = 36 )
      PARAMETER ( ALPHA = 1.0D0    )
C     PARAMETER ( NLINK1 = 10 )
      PARAMETER ( RRCUT = 4.724D0 )
      PARAMETER ( RRCUT1= 5.669D0 )
      PARAMETER ( DRCUT = 0.945D0 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          FREEZE*5,PRINT*5,DGF*3

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI2 / LX(0:NX*NY*NZ-1),LY(0:NX*NY*NZ-1),
     *                 LZ(0:NX*NY*NZ-1)
      COMMON / MPIBUF / TBUF(  NMAXX,NMAXY,NMAXZ ),
     *                  TBUF1( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PRMT4 / nion1,iiont(maxatm)
C     COMMON / SLT / VPCE(NMAX,NMAX,NMAX),VPCZ,PLJ
      COMMON / SLT / VPCE(NMAXX/NX,NMAXY/NY,NMAXZ/NZ),VPCZ,PLJ
      COMMON / LJQMMM / SQMMM(NNUC,maxatm),EPQMMM(NNUC,maxatm),
     *                  chgmm(maxatm),chgqm(NNUC)
      COMMON / MMDAT / SCRD( NDIM,maxatm ),SCRD1( NDIM,maxatm ),   ! Vmol
     *                 FRCS( maxatm,NDIM )                         ! Vmol

      DIMENSION STAT(MPI_STATUS_SIZE)

      DIMENSION ADX( NMAXX/NX )
      DIMENSION ADY( NMAXY/NY )
      DIMENSION ADZ( NMAXZ/NZ )
      DIMENSION HADX( NMAXX/2 )
      DIMENSION HADY( NMAXY/2 )
      DIMENSION HADZ( NMAXZ/2 )
      DIMENSION HVPCE( NMAXX/2,NMAXY/2,NMAXZ/2 )
      DIMENSION AVPCE( NMAXX/2,NMAXY/2,NMAXZ/2 )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      ONE    = 1.D0
      RRCUT2 = ONE/(RRCUT*RRCUT)

      NMAXXH  = NMAXX/2
      NMAXYH  = NMAXY/2
      NMAXZH  = NMAXZ/2

      MNMAXX  = NMAXX/2 + 2
      MNMAXY  = NMAXY/2 + 2
      MNMAXZ  = NMAXZ/2 + 2

      NNMAXX  = NMAXX/2 + 1
      NNMAXY  = NMAXY/2 + 1
      NNMAXZ  = NMAXZ/2 + 1

      NNX     = (1-MX)*NMAXXH + 1
      NNY     = (1-MY)*NMAXYH + 1
      NNZ     = (1-MZ)*NMAXZH + 1

C     IF(IMD.NE.0) go to 300
      IF(IMD.NE.1) go to 300
C
C     write electronic field ( for dense grid )
C
      OPEN(unit=200,file='VPCE-1.dat',status='unknown')
      OPEN(unit=210,file='VPCE-2.dat',status='unknown')

      IF(NQMMM .EQ. 'MM') THEN      ! Vmol

      DO K=1,JZ
        ADZ(K) = DZ*( K-NNMAXZ )
      ENDDO
      DO J=1,JY
        ADY(J) = DY*( J-NNMAXY )
      ENDDO
      DO I=1,JX
        ADX(I) = DX*( I-NNMAXX )
      ENDDO

      DO K=1,JZ
      DO J=1,JY
      DO I=1,JX
        VPCE( I,J,K )=0.0D0            ! initialize
      ENDDO
      ENDDO
      ENDDO

      DO 70 iin=1,nion1
        in = iiont(iin)

       IF(NMMSW(in).EQ.1) GOTO 70

       DO    K=1,JZ
       DO    J=1,JY
       DO 75 I=1,JX

        RX    = ADX(I) - SCRD1( 1,in )
        RY    = ADY(J) - SCRD1( 2,in )
        RZ    = ADZ(K) - SCRD1( 3,in )
        R2    = RX**2+RY**2+RZ**2
        R     = DSQRT( R2 )

        IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
          ESW  = ONE
        ELSEIF(R.GT.DRCUT) THEN
          ARD  = R - DRCUT
          ARD2 = ARD*ARD*RRCUT2
          ESW1 = ONE - ARD2
          ESW2 = ESW1*ESW1
          ESW  = ONE - ESW2
        ELSE
          GOTO 75
        ENDIF

        VPCE( I,J,K )= VPCE( I,J,K )
     *               - ESW*chgmm( in )*(ERF(ALPHA*R)/R)

75     CONTINUE
       ENDDO
       ENDDO

70    CONTINUE

      ENDIF

      IF( MYID .EQ. 0 ) THEN

C     NNMAXP = NMAX/2
      NUMDAT = NMAXX*NMAXY*NMAXZ/(NX*NY*NZ)

      I1 = MX*NMAXXH
      J1 = MY*NMAXYH
      K1 = MZ*NMAXZH

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        TBUF(I+I1,J+J1,K+K1) = VPCE(I,J,K)
      ENDDO
      ENDDO
      ENDDO

      DO ID = 1, NUMPROCS - 1

        CALL MPI_RECV(TBUF1(1,1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                ID,MPI_ANY_TAG,MPI_COMM_WORLD,STAT,IERR)

        I1 = LX(ID)*NMAXXH
        J1 = LY(ID)*NMAXYH
        K1 = LZ(ID)*NMAXZH

        DO K = 1, JZ
        DO J = 1, JY
        DO I = 1, JX
          TBUF(I+I1,J+J1,K+K1) = TBUF1(I,J,K)
        ENDDO
        ENDDO
        ENDDO

      ENDDO

      WRITE(200,*) IMD
      DO K=1,NMAXZ
      DO J=1,NMAXY
      DO I=1,NMAXX
        WRITE(200,*) TBUF( I,J,K )
      ENDDO
      ENDDO
      ENDDO

C     CLOSE(unit=200)

      ELSE                                                 ! For slave CPU

        NUMDAT = NMAXX*NMAXY*NMAXZ/(NX*NY*NZ)
        CALL MPI_SEND(VPCE(1,1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                0,1,MPI_COMM_WORLD,IERR)

      ENDIF

      NVPCE = 0
      DO K=1,NMAXZH
      DO J=1,NMAXYH
      DO I=1,NMAXXH
        HVPCE( I,J,K ) = 0.d0
      ENDDO
      ENDDO
      ENDDO

300   CONTINUE
C
C     averaging electronic field ( for coarse grid )
C
      IF( MYID .EQ. 0 ) THEN  ! For master CPU

      NVPCE = NVPCE + 1
      HDX = 0.5d0*DX
      HDY = 0.5d0*DY
      HDZ = 0.5d0*DZ
      DO K=1,NMAXZH
        HADZ(K) = DZ*( 2*K-MNMAXZ ) + HDZ
      ENDDO
      DO J=1,NMAXYH
        HADY(J) = DY*( 2*J-MNMAXY ) + HDY
      ENDDO
      DO I=1,NMAXXH
        HADX(I) = DX*( 2*I-MNMAXX ) + HDX
      ENDDO

      DO 80 iin=1,nion1
        in = iiont(iin)

       IF(NMMSW(in).EQ.1) GOTO 80

       DO    K=1,NMAXZH
       DO    J=1,NMAXYH
       DO 85 I=1,NMAXXH

        RX    = HADX(I) - SCRD1( 1,in )
        RY    = HADY(J) - SCRD1( 2,in )
        RZ    = HADZ(K) - SCRD1( 3,in )
        R2    = RX**2+RY**2+RZ**2
        R     = DSQRT( R2 )

        IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
          ESW  = ONE
        ELSEIF(R.GT.DRCUT) THEN
          ARD  = R - DRCUT
          ARD2 = ARD*ARD*RRCUT2
          ESW1 = ONE - ARD2
          ESW2 = ESW1*ESW1
          ESW  = ONE - ESW2
        ELSE
          GOTO 85
        ENDIF

        HVPCE( I,J,K )= HVPCE( I,J,K )
     *                - ESW*chgmm( in )*(ERF(ALPHA*R)/R)

85     CONTINUE
       ENDDO
       ENDDO

80    CONTINUE
C
C     write electronic field ( for coarse grid )
C
      IF(MOD(IMD,5000).EQ.0) THEN
        ANVPCE = 1.d0/DBLE(NVPCE)
        WRITE(210,*) IMD
        DO K=1,NMAXZH
        DO J=1,NMAXYH
        DO I=1,NMAXXH
          AVPCE( I,J,K ) = ANVPCE * HVPCE( I,J,K )
          WRITE(210,*) AVPCE( I,J,K )
        ENDDO
        ENDDO
        ENDDO
      ENDIF

      ENDIF  ! For master CPU
C
C     close file
C
C     IF(IMD.EQ.MDMAX) THEN
C       CLOSE(unit=210)
C     ENDIF

      RETURN
      END

C--------------------------------------------------------
C     SUBROUTINE INITIAL HARTREE POTENTIAL
C--------------------------------------------------------

      SUBROUTINE IVCOU( CRD )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
      INTEGER*4 STAT

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol

      PARAMETER( EPS = 1.D-1 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI2 / LX(0:NX*NY*NZ-1),LY(0:NX*NY*NZ-1),
     *                 LZ(0:NX*NY*NZ-1)
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / MPIBUF / TBUF(  NMAXX,NMAXY,NMAXZ ),
     *                  TBUF1( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      COMMON / GRID2 / DX, DY, DZ
C     COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
C     COMMON / OFC / VCOU( NMAX,NMAX,NMAX )
      COMMON / OFC / VCOU( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      DIMENSION STAT(MPI_STATUS_SIZE)

      DIMENSION CRD(   NNUC,NDIM )
      DIMENSION TVCOU( NMAXX,NMAXY,NMAXZ )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

C----- read -----

      IF(MYID.EQ.0 ) THEN

C     READ(23,*)
C     READ(23,*)
C     READ(23,*) N0,X0,Y0,Z0
C     READ(23,*) N1,X1,Y1,Z1
C     READ(23,*) N2,X2,Y2,Z2
C     READ(23,*) N3,X3,Y3,Z3

      DO 10  L=1,NNUC
C     READ(23,*)
   10 CONTINUE 

      DO 30 I=1,N1
      DO 30 J=1,N2
C       READ(23,'(6E13.5)') (TVCOU(I,J,K),K=1,N3)
   30 CONTINUE

C-----------------------------------------------------------

      NNMAXX = NMAXX/2+1
      NNMAXY = NMAXY/2+1
      NNMAXZ = NMAXZ/2+1

      DO 40 K = 1, NMAXZ

        RZ1 = DZ*( K-NNMAXZ )

        DO 40 J = 1, NMAXY

          RY1 = DY*( J-NNMAXY )

          DO 40 I = 1, NMAXX

            RX1 = DX*( I-NNMAXX )

            DO 40 NA = 1, NNUC

              RX = RX1 - CRD(NA,1) 
              RY = RY1 - CRD(NA,2)
              RZ = RZ1 - CRD(NA,3)
              R2 = RX**2 + RY**2 + RZ**2
              R  = DSQRT(R2)

              IF( R .GT. EPS ) THEN
                TVCOU( I,J,K ) = TVCOU(I,J,K)+ZVAL(NA)*DERF(R)/R
              ENDIF

40    CONTINUE

      NCOUNT = 0

      DO K = 1, NMAXZ
        RZ1 = DZ*( K-NNMAXZ )
      DO J = 1, NMAXY
        RY1 = DY*( J-NNMAXY )
      DO I = 1, NMAXX
        RX1 = DX*( I-NNMAXX )
        TBUF( I,J,K ) = TVCOU( I,J,K )
      DO NA = 1, NNUC
        RX = RX1 - CRD(NA,1) 
        RY = RY1 - CRD(NA,2)
        RZ = RZ1 - CRD(NA,3)
        R2 = RX**2 + RY**2 + RZ**2
        R  = DSQRT(R2)
        IF( R .LE. EPS ) THEN
          TBUF( I,J,K ) = (1.d0/6.d0) 
     *                  * ( TVCOU( I-1,J,K ) + TVCOU( I+1,J,K )
     *                  +   TVCOU( I,J-1,K ) + TVCOU( I,J+1,K )
     *                  +   TVCOU( I,J,K-1 ) + TVCOU( I,J,K+1 ) )
          write(*,*) 'VCOU:',TVCOU(I,J,K),' => ',TBUF(I,J,K)
          write(*,*) '     ',I,J,K
          NCOUNT = NCOUNT+1
          GOTO 45
        ENDIF
      ENDDO
45      CONTINUE
        IF(TBUF(I,J,K).LT.0.d0) write(*,*) 'Warning: IVCOU'
      ENDDO
      ENDDO
      ENDDO

      IF(NCOUNT.EQ.1) THEN
        write(*,100) 'Warning(IVCOU)!: There is ',NCOUNT,
     *             ' grid point which is too close to a nucleus site.  '
        write(*,*)
      ELSEIF(NCOUNT.GT.1) THEN
        write(*,100) 'Warning(IVCOU)!: There are',NCOUNT,
     *             ' grid points which are too close to a nucleus site.'
        write(*,*)
      ENDIF

100   FORMAT(1x,a26,i4,a51)

      NUMDAT = JX*JY*JZ

      I1 = MX*JX
      J1 = MY*JY 
      K1 = MZ*JZ 

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        VCOU(I,J,K) = TBUF(I+I1,J+J1,K+K1)
      ENDDO
      ENDDO
      ENDDO

      DO ID = 1, NUMPROCS - 1

        I1 = LX(ID)*JX
        J1 = LY(ID)*JY 
        K1 = LZ(ID)*JZ 

        DO K = 1, JZ
        DO J = 1, JY
        DO I = 1, JX
          TBUF1(I,J,K) = TBUF(I+I1,J+J1,K+K1)
        ENDDO
        ENDDO
        ENDDO

        CALL MPI_SEND(TBUF1(1,1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                ID,1,MPI_COMM_WORLD,IERR)

      ENDDO

      ELSE                           ! For slave CPU

        NUMDAT = JX*JY*JZ
        CALL MPI_RECV(VCOU(1,1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                 0,MPI_ANY_TAG,MPI_COMM_WORLD,STAT,IERR)

      ENDIF

      RETURN
      END

C------------------------------------
C   SUBROUTINE KINETIC  
C------------------------------------

      SUBROUTINE KNTC( RWF,NOR,TWF )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'
      include 'QMpara.i'                        ! Vmol

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI3 / NIDX(1:2),NIDY(1:2),NIDZ(1:2)
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / ELCTRN / AX1, AX2, AY1, AY2, AZ1, AZ2, AO,
     *                  AX3, AX4, AY3, AY4, AZ3, AZ4

      DIMENSION RWF( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NOR )
      DIMENSION TWF( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NOR )

      DIMENSION RWF1( -3:NMAXX/NX+4,-3:NMAXY/NY+4,-3:NMAXZ/NZ+4 )
      DIMENSION A( 8*NMAXX*NMAXX/(NX*NY*NZ) )
      DIMENSION B( 8*NMAXX*NMAXX/(NX*NY*NZ) )

      DO LL = 1, NOR

      RWF1(:,:,:) = 0.D0
      RWF1(1:JX,1:JY,1:JZ) = RWF(1:JX,1:JY,1:JZ,LL)
      CALL SENDRECV_BUF(RWF1,JX,JY,JZ,4)

      DO  K = 1, JZ
      DO  J = 1, JY
      DO  I = 1, JX

        TWF( I,J,K,LL ) = AZ1 * ( RWF1(I,J,K-1) + RWF1(I,J,K+1) )
     *                  + AZ2 * ( RWF1(I,J,K-2) + RWF1(I,J,K+2) )
     *                  + AZ3 * ( RWF1(I,J,K-3) + RWF1(I,J,K+3) )
     *                  + AZ4 * ( RWF1(I,J,K-4) + RWF1(I,J,K+4) )
     *                  + AY1 * ( RWF1(I,J-1,K) + RWF1(I,J+1,K) )
     *                  + AY2 * ( RWF1(I,J-2,K) + RWF1(I,J+2,K) )
     *                  + AY3 * ( RWF1(I,J-3,K) + RWF1(I,J+3,K) )
     *                  + AY4 * ( RWF1(I,J-4,K) + RWF1(I,J+4,K) )
     *                  + AX1 * ( RWF1(I-1,J,K) + RWF1(I+1,J,K) )
     *                  + AX2 * ( RWF1(I-2,J,K) + RWF1(I+2,J,K) )
     *                  + AX3 * ( RWF1(I-3,J,K) + RWF1(I+3,J,K) )
     *                  + AX4 * ( RWF1(I-4,J,K) + RWF1(I+4,J,K) )
     *                  + AO  *   RWF1(I,J,K)

      ENDDO
      ENDDO
      ENDDO

      ENDDO

      RETURN
      END


C---------------------------------------------
C   SUBROUTINE PRECONDITIONING 
C
C   REFERENCE : Phys. Rev. B 52 5459( 1995 ) 
C--------------------------------------------

      SUBROUTINE PC( RVEC,TRVEC,NOR ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX =  80 )
C     PARAMETER ( NORA =  49 )
      PARAMETER ( C1 = 0.5D0 )

      DIMENSION RVEC1(0:NMAXX/NX+1,0:NMAXY/NY+1,0:NMAXZ/NZ+1)
      DIMENSION RVEC( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION TRVEC(NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
C     DIMENSION A( 2*NMAX*NMAX/(NX*NY*NZ) )
C     DIMENSION B( 2*NMAX*NMAX/(NX*NY*NZ) )
      PARAMETER( JMAX = MAX(NMAXX/NX,NMAXY/NY,NMAXZ/NZ) )
      DIMENSION A( JMAX*JMAX )
      DIMENSION B( JMAX*JMAX )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI3 / NIDX(1:2),NIDY(1:2),NIDZ(1:2)
      COMMON / MMPI4 / JX,JY,JZ

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      AZ = 1.D0 / 6.D0
      AY = 1.D0 / 6.D0
      AX = 1.D0 / 6.D0

      DO L = 1, NOR

      RVEC1(:,:,:) = 0D0
      RVEC1(1:JX,1:JY,1:JZ) = RVEC(1:JX,1:JY,1:JZ,L)
      CALL SENDRECV_BUF(RVEC1,JX,JY,JZ,1)

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        TRVEC(I,J,K,L) = AZ*( RVEC1(I,J,K-1) + RVEC1(I,J,K+1) )
     *                 + AY*( RVEC1(I,J-1,K) + RVEC1(I,J+1,K) )
     *                 + AX*( RVEC1(I-1,J,K) + RVEC1(I+1,J,K) )
        RVEC (I,J,K,L) = C1*( RVEC (I,J,K,L) + TRVEC(I,J,K,L) )
      ENDDO
      ENDDO
      ENDDO

      ENDDO

      RETURN
      END

      SUBROUTINE RCLMB5( PVEC,VBC,APVEC )

      IMPLICIT REAL*8 ( A-H,O-Z )
      IMPLICIT INTEGER*4 ( I-N )

      include "mpif.h"
      include "mpi.i"
      include 'QMpara.i'
 
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,
     *                  RSIZE(NNUC),
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF
 
      DIMENSION PVEC ( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION APVEC( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VBC  ( -3:NMAXX/NX+4,-3:NMAXY/NY+4,-3:NMAXZ/NZ+4 )
      DIMENSION PVEC1( -3:NMAXX/NX+4,-3:NMAXY/NY+4,-3:NMAXZ/NZ+4 )
     
      PVEC1(:,:,:) = 0.D0
      PVEC1(1:JX,1:JY,1:JZ) = PVEC(1:JX,1:JY,1:JZ)
      CALL SENDRECV_BUF(PVEC1,JX,JY,JZ,4)
      PVEC1(:,:,:) = PVEC1(:,:,:) + VBC(:,:,:)

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX

        APVEC( I,J,K ) = TAZ1 * ( PVEC1( I,J,K-1 ) + PVEC1( I,J,K+1 ) )
     *                 + TAZ2 * ( PVEC1( I,J,K-2 ) + PVEC1( I,J,K+2 ) )
     *                 + TAZ3 * ( PVEC1( I,J,K-3 ) + PVEC1( I,J,K+3 ) )
     *                 + TAZ4 * ( PVEC1( I,J,K-4 ) + PVEC1( I,J,K+4 ) )
     *                 + TAY1 * ( PVEC1( I,J-1,K ) + PVEC1( I,J+1,K ) )
     *                 + TAY2 * ( PVEC1( I,J-2,K ) + PVEC1( I,J+2,K ) )
     *                 + TAY3 * ( PVEC1( I,J-3,K ) + PVEC1( I,J+3,K ) )
     *                 + TAY4 * ( PVEC1( I,J-4,K ) + PVEC1( I,J+4,K ) )
     *                 + TAX1 * ( PVEC1( I-1,J,K ) + PVEC1( I+1,J,K ) )
     *                 + TAX2 * ( PVEC1( I-2,J,K ) + PVEC1( I+2,J,K ) )
     *                 + TAX3 * ( PVEC1( I-3,J,K ) + PVEC1( I+3,J,K ) )
     *                 + TAX4 * ( PVEC1( I-4,J,K ) + PVEC1( I+4,J,K ) )
     *                 + TAO  *   PVEC1( I,J,K )

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END
C--------------------------------------------------------
C     SUBROUTINE Read Electronic Field formed by MM molecules
C--------------------------------------------------------

      SUBROUTINE RVPCE( IMD )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include 'QMpara.i'                        ! Vmol
      include "mpi.i"

C     PARAMETER ( NMAX = 80 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          FREEZE*5,PRINT*5,DGF*3

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI2 / LX(0:NX*NY*NZ-1),LY(0:NX*NY*NZ-1),
     *                 LZ(0:NX*NY*NZ-1)
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / MPIBUF / TBUF(  NMAXX,NMAXY,NMAXZ ),
     *                  TBUF1( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
C     COMMON / SLT / VPCE(NMAX,NMAX,NMAX),VPCZ,PLJ
      COMMON / SLT / VPCE(NMAXX/NX,NMAXY/NY,NMAXZ/NZ),VPCZ,PLJ

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      IF( IMD .NE. 0 ) RETURN

      IF( MYID .NE. 0 ) THEN

C     For dense-grid

      IF( NQMMM.EQ.'QMVD' ) THEN

        WRITE(*,*)
        WRITE(*,*) 'Read VPCE(dense) data'

        READ(220,*) NO
        DO K = 1, NMAXZ
        DO J = 1, NMAXY
        DO I = 1, NMAXX
C         READ(220,*) VPCE(I,J,K)
          READ(220,*) TBUF(I,J,K)
        ENDDO
        ENDDO
        ENDDO

C     For coarse-grid

      ELSEIF( NQMMM.EQ.'QMVC' ) THEN

        WRITE(*,*)
        WRITE(*,*) 'Read VPCE(coarse) data'
  
        NMAXHX  = NMAXX/2
        NMAXHY  = NMAXY/2
        NMAXHZ  = NMAXZ/2
  
        READ(230,*) NO
        DO K = 1, NMAXHZ
        DO J = 1, NMAXHY
        DO I = 1, NMAXHX
  
          READ(230,*) BVPCE
  
          DO KK = 2*K-1, 2*K
          DO JJ = 2*J-1, 2*J
          DO II = 2*I-1, 2*I
C           VPCE(II,JJ,KK) = BVPCE
            TBUF(II,JJ,KK) = BVPCE
          ENDDO
          ENDDO
          ENDDO
  
        ENDDO
        ENDDO
        ENDDO

      ENDIF

      NMAXHX  = NMAXX/2
      NMAXHY  = NMAXY/2
      NMAXHZ  = NMAXZ/2
      NUMDAT  = NMAXX*NMAXY*NMAXZ/(NX*NY*NZ)

      I1 = MX*NMAXHX  
      J1 = MY*NMAXHY  
      K1 = MZ*NMAXHZ  

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        VPCE(I,J,K) = TBUF(I+I1,J+J1,K+K1)
      ENDDO
      ENDDO
      ENDDO

      DO ID = 1, NUMPROCS - 1

        I1 = LX(ID)*NMAXHX  
        J1 = LY(ID)*NMAXHY  
        K1 = LZ(ID)*NMAXHZ  

        DO K = 1, JZ
        DO J = 1, JY
        DO I = 1, JX
          TBUF1(I,J,K) = TBUF(I+I1,J+J1,K+K1)
        ENDDO
        ENDDO
        ENDDO

        CALL MPI_SEND(TBUF1(1,1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                ID,1,MPI_COMM_WORLD,IERR)

      ENDDO

      ELSE                                                 ! For slave CPU

        NUMDAT = NMAXX*NMAXY*NMAXZ/(NX*NY*NZ)
        CALL MPI_RECV(VPCE(1,1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                0,MPI_ANY_TAG,MPI_COMM_WORLD,STAT,IERR)

      ENDIF

      RETURN
      END

C-------------------------------------
C     SUBROUTINE  SEND,RECV
C-------------------------------------

      SUBROUTINE SEND(A,NUMDAT,NDEST)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"

      DIMENSION A(NUMDAT)

      CALL MPI_SEND(A,NUMDAT,MPI_DOUBLE_PRECISION,
     *              NDEST,NDEST,MPI_COMM_WORLD,IERR)

      RETURN
      END

      SUBROUTINE RECV(B,NUMDAT,NROOT)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
      INTEGER*4 STAT

      include "mpif.h"

      DIMENSION B(NUMDAT)
      DIMENSION STAT(MPI_STATUS_SIZE)

      CALL MPI_RECV(B,NUMDAT,MPI_DOUBLE_PRECISION,
     *              NROOT,MPI_ANY_TAG,
     *              MPI_COMM_WORLD,STAT,IERR)

      RETURN
      END

      SUBROUTINE SENDRECV(A,B,NUMDAT,NDEST,NROOT)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
      INTEGER*4 STAT

      include "mpif.h"

      DIMENSION A(NUMDAT), B(NUMDAT)
      DIMENSION STAT(MPI_STATUS_SIZE)

      IF(NDEST.NE.-1 .AND. NROOT.NE.-1) THEN
         CALL MPI_SENDRECV(
     *        A,NUMDAT,MPI_DOUBLE_PRECISION,NDEST,NDEST,
     *        B,NUMDAT,MPI_DOUBLE_PRECISION,NROOT,MPI_ANY_TAG,
     *        MPI_COMM_WORLD,STAT,IERR)
      ELSE
         IF(NDEST.NE.-1) CALL SEND(A,NUMDAT,NDEST)
         IF(NROOT.NE.-1) CALL RECV(B,NUMDAT,NROOT)
      ENDIF

      RETURN
      END

      SUBROUTINE SENDRECV_BUF(ARRAY,JX,JY,JZ,NB)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

      PARAMETER( JMAX = MAX(NMAXX/NX,NMAXY/NY,NMAXZ/NZ) )
      PARAMETER( NBMAX = 4 )

      COMMON / MMPI3 / NIDX(1:2),NIDY(1:2),NIDZ(1:2)

      DIMENSION ARRAY(1-NB:JX+NB,1-NB:JY+NB,1-NB:JZ+NB)
      DIMENSION A( NBMAX*JMAX*JMAX )
      DIMENSION B( NBMAX*JMAX*JMAX )

      IF( NB > NBMAX ) THEN
         WRITE(*,*) "Error: NB > NBMAX in SENDRECV_BUF"
         STOP
      ENDIF

C----- z-axis -----

      NCOUNT = NB*JX*JY
      A(1:NCOUNT) = RESHAPE( ARRAY(1:JX,1:JY,JZ-NB+1:JZ), (/NCOUNT/) )
      B(1:NCOUNT) = 0D0
      CALL SENDRECV(A,B,NCOUNT,NIDZ(2),NIDZ(1))
      ARRAY(1:JX,1:JY,1-NB:0) = RESHAPE( B(1:NCOUNT), (/JX,JY,NB/) )

      A(1:NCOUNT) = RESHAPE( ARRAY(1:JX,1:JY,1:NB), (/NCOUNT/) )
      B(1:NCOUNT) = 0D0
      CALL SENDRECV(A,B,NCOUNT,NIDZ(1),NIDZ(2))
      ARRAY(1:JX,1:JY,JZ+1:JZ+NB) = RESHAPE( B(1:NCOUNT), (/JX,JY,NB/) )

C----- y-axis -----

      NCOUNT = NB*JX*JZ
      A(1:NCOUNT) = RESHAPE( ARRAY(1:JX,JY-NB+1:JY,1:JZ), (/NCOUNT/) )
      B(1:NCOUNT) = 0D0
      CALL SENDRECV(A,B,NCOUNT,NIDY(2),NIDY(1))
      ARRAY(1:JX,1-NB:0,1:JZ) = RESHAPE( B(1:NCOUNT), (/JX,NB,JZ/) )

      A(1:NCOUNT) = RESHAPE( ARRAY(1:JX,1:NB,1:JZ), (/NCOUNT/) )
      B(1:NCOUNT) = 0D0
      CALL SENDRECV(A,B,NCOUNT,NIDY(1),NIDY(2))
      ARRAY(1:JX,JY+1:JY+NB,1:JZ) = RESHAPE( B(1:NCOUNT), (/JX,NB,JZ/) )

C----- x-axis -----

      NCOUNT = NB*JY*JZ
      A(1:NCOUNT) = RESHAPE( ARRAY(JX-NB+1:JX,1:JY,1:JZ), (/NCOUNT/) )
      B(1:NCOUNT) = 0D0
      CALL SENDRECV(A,B,NCOUNT,NIDX(2),NIDX(1))
      ARRAY(1-NB:0,1:JY,1:JZ) = RESHAPE( B(1:NCOUNT), (/NB,JY,JZ/) )

      A(1:NCOUNT) = RESHAPE( ARRAY(1:NB,1:JY,1:JZ), (/NCOUNT/) )
      B(1:NCOUNT) = 0D0
      CALL SENDRECV(A,B,NCOUNT,NIDX(1),NIDX(2))
      ARRAY(JX+1:JX+NB,1:JY,1:JZ) = RESHAPE( B(1:NCOUNT), (/NB,JY,JZ/) )

      RETURN
      END

      FUNCTION ASUM(D)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"

      REAL*8 DSUM

      CALL MPI_ALLREDUCE(D,ASUM,1,MPI_DOUBLE_PRECISION,
     *     MPI_SUM,MPI_COMM_WORLD,IERR)
      
      RETURN
      END

      SUBROUTINE PRINTALL(A)

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

      PARAMETER( JMAX = MAX(NMAXX/NX,NMAXY/NY,NMAXZ/NZ) )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      DIMENSION A(JX,JY,JZ)
      DIMENSION B(NMAXX,NMAXY,NMAXZ)
      DIMENSION C(NMAXX,NMAXY,NMAXZ)

      CHARACTER FILENAME*9

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      NNX = MX*JX
      NNY = MY*JY
      NNZ = MZ*JZ

      B = 0d0
      B(NNX+1:NNX+JX,NNY+1:NNY+JY,NNZ+1:NNZ+JZ) = A(1:JX,1:JY,1:JZ)

      NUMDAT = NMAXX*NMAXY*NMAXZ

      C = 0d0
      CALL MPI_ALLREDUCE(B,C,NUMDAT,MPI_DOUBLE_PRECISION,
     *     MPI_SUM,MPI_COMM_WORLD,IERR)

      IF(MYID.eq.0) THEN
         WRITE(FILENAME,'("np",3i1,".txt")') NX,NY,NZ
         OPEN(999,FILE=FILENAME)
         DO K=1,NMAXZ
         DO J=1,NMAXY
         DO I=1,NMAXX
            WRITE(999,*)I,J,K,C(I,J,K)
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      RETURN
      END
C     Subroutine send/receive data for Poisson Equation
      SUBROUTINE SRPSN( PVEC,VCOUX,VCOUY,VCOUZ )

      IMPLICIT REAL*8( A-H,O-Z )
      IMPLICIT INTEGER*4( I-N )
      INTEGER*4 STAT

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'

      PARAMETER( INUM = 8 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI3 / NIDX(1:2),NIDY(1:2),NIDZ(1:2)
      COMMON / MMPI4 / JX,JY,JZ

      DIMENSION STAT(MPI_STATUS_SIZE)

      DIMENSION PVEC(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VCOUX( NMAXY/NY,NMAXZ/NZ,INUM ) ! note: Dimension(Y,Z,X)
      DIMENSION VCOUY( NMAXX/NX,NMAXZ/NZ,INUM ) ! note: Dimension(X,Z,Y)
      DIMENSION VCOUZ( NMAXX/NX,NMAXY/NY,INUM )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      INUM2  = INUM/2
      INUM3  = INUM2 + 1

C----- X-AXIS -----

      IF( NX .EQ. 2 ) THEN
      IF( MX .EQ. 0 ) THEN

        DO K = 1, JZ
        DO J = 1, JY
        DO I = JX-3, JX 
          I1 = JX + 1 - I
          VCOUX( J,K,I1 ) = PVEC( I,J,K )
        ENDDO
        ENDDO
        ENDDO
  
        NUMDAT = JY*JZ*INUM2
        CALL MPI_SEND(VCOUX(1,1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                NIDX(2),1,MPI_COMM_WORLD,IERR)
        CALL MPI_RECV(VCOUX(1,1,INUM3),NUMDAT,MPI_DOUBLE_PRECISION,
     *                NIDX(2),MPI_ANY_TAG,MPI_COMM_WORLD,STAT,IERR)

      ELSEIF( MX .EQ. 1 ) THEN

        DO K = 1, JZ
        DO J = 1, JY
        DO I = 1, INUM2
          I1 = I + INUM2
          VCOUX( J,K,I1 ) = PVEC( I,J,K )
        ENDDO
        ENDDO
        ENDDO

        NUMDAT = JY*JZ*INUM2
        CALL MPI_RECV(VCOUX(1,1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                NIDX(1),MPI_ANY_TAG,MPI_COMM_WORLD,STAT,IERR)
        CALL MPI_SEND(VCOUX(1,1,INUM3),NUMDAT,MPI_DOUBLE_PRECISION,
     *                NIDX(1),1,MPI_COMM_WORLD,IERR)

      ENDIF
      ENDIF

C----- Y-AXIS -----

      IF( NY .EQ. 2 ) THEN
      IF( MY .EQ. 0 ) THEN

        DO K = 1, JZ
        DO J = JY-3, JY
          J1 = JY + 1 - J
        DO I = 1, JX 
          VCOUY( I,K,J1 ) = PVEC( I,J,K )
        ENDDO
        ENDDO
        ENDDO
  
        NUMDAT = JX*JZ*INUM2
        CALL MPI_SEND(VCOUY(1,1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                NIDY(2),1,MPI_COMM_WORLD,IERR)
        CALL MPI_RECV(VCOUY(1,1,INUM3),NUMDAT,MPI_DOUBLE_PRECISION,
     *                NIDY(2),MPI_ANY_TAG,MPI_COMM_WORLD,STAT,IERR)

      ELSEIF( MY .EQ. 1 ) THEN

        DO K = 1, JZ
        DO J = 1, INUM2
          J1 = J + INUM2
        DO I = 1, JX
          VCOUY( I,K,J1 ) = PVEC( I,J,K )
        ENDDO
        ENDDO
        ENDDO

        NUMDAT = JX*JZ*INUM2
        CALL MPI_RECV(VCOUY(1,1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                NIDY(1),MPI_ANY_TAG,MPI_COMM_WORLD,STAT,IERR)
        CALL MPI_SEND(VCOUY(1,1,INUM3),NUMDAT,MPI_DOUBLE_PRECISION,
     *                NIDY(1),1,MPI_COMM_WORLD,IERR)

      ENDIF
      ENDIF

C----- Z-AXIS -----

      IF( NZ .EQ. 2 ) THEN
      IF( MZ .EQ. 0 ) THEN

        DO K = JZ-3, JZ
          K1 = JZ + 1 - K
        DO J = 1, JY
        DO I = 1, JX 
          VCOUZ( I,J,K1 ) = PVEC( I,J,K )
        ENDDO
        ENDDO
        ENDDO
  
        NUMDAT = JX*JY*INUM2
        CALL MPI_SEND(VCOUZ(1,1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                NIDZ(2),1,MPI_COMM_WORLD,IERR)
        CALL MPI_RECV(VCOUZ(1,1,INUM3),NUMDAT,MPI_DOUBLE_PRECISION,
     *                NIDZ(2),MPI_ANY_TAG,MPI_COMM_WORLD,STAT,IERR)

      ELSEIF( MZ .EQ. 1 ) THEN

        DO K = 1, INUM2
          K1 = K + INUM2
        DO J = 1, JY
        DO I = 1, JX
          VCOUZ( I,J,K1 ) = PVEC( I,J,K )
        ENDDO
        ENDDO
        ENDDO

        NUMDAT = JX*JY*INUM2
        CALL MPI_RECV(VCOUZ(1,1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                NIDZ(1),MPI_ANY_TAG,MPI_COMM_WORLD,STAT,IERR)
        CALL MPI_SEND(VCOUZ(1,1,INUM3),NUMDAT,MPI_DOUBLE_PRECISION,
     *                NIDZ(1),1,MPI_COMM_WORLD,IERR)

      ENDIF
      ENDIF

      RETURN
      END

C--------------------------------------------------------
C     SUBROUTINE Write Electronic Density of the QM Subsytem
C--------------------------------------------------------

      SUBROUTINE WDNS( IMD,RHO )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
      INTEGER*4 STAT

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol

      COMMON / MMPI4 / JX,JY,JZ
      COMMON / MMPI2 / LX(0:NX*NY*NZ-1),LY(0:NX*NY*NZ-1),
     *                 LZ(0:NX*NY*NZ-1)
      COMMON / MPIBUF / TBUF( NMAXX,NMAXY,NMAXZ ),
     *                  TBUF1( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

C     PARAMETER ( NMAX = 80 )

      DIMENSION STAT(MPI_STATUS_SIZE)

      DIMENSION RHO( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      IF(MYID.EQ.0 ) THEN

      NUMDAT = JX*JY*JZ

      DO ID = 1, NUMPROCS - 1

        CALL MPI_RECV(TBUF1(1,1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                ID,MPI_ANY_TAG,MPI_COMM_WORLD,STAT,IERR)

        I1 = LX(ID)*JX
        J1 = LY(ID)*JY 
        K1 = LZ(ID)*JZ 

        DO K = 1, JZ
        DO J = 1, JY
        DO I = 1, JX
          TBUF(I+I1,J+J1,K+K1) = TBUF1(I,J,K)
        ENDDO
        ENDDO
        ENDDO

      ENDDO

      REWIND(31)
      WRITE(31,*) IMD
      DO K = 1, NMAXZ
      DO J = 1, NMAXY
      DO I = 1, NMAXX
C       WRITE(31,*) RHO(I,J,K)
        WRITE(31,*) TBUF(I,J,K)
      ENDDO
      ENDDO
      ENDDO

      ELSE                                                 ! For slave CPU

        NUMDAT = JX*JY*JZ
        CALL MPI_SEND(RHO(1,1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                0,1,MPI_COMM_WORLD,IERR)

      ENDIF

      RETURN
      END


C--------------------------------------------------------
C     SUBROUTINE INITIAL HARTREE POTENTIAL
C     specific for use of the system which consists
C     of N water molecules placed on uniform
C     grids.
C--------------------------------------------------------

      SUBROUTINE IVCOU_NW( CRD )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
      INTEGER*4 STAT

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol

      PARAMETER( EPS   = 1.D-1 )
      PARAMETER( RCUT  = 3.5D0 )
      PARAMETER( ALPHA = 1.0D0 )
      PARAMETER( DELTA = 0.8D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI2 / LX(0:NX*NY*NZ-1),LY(0:NX*NY*NZ-1),
     *                 LZ(0:NX*NY*NZ-1)
      COMMON / MMPI4 / JX,JY,JZ
C     COMMON / MPIBUF / TBUF(  NMAXX,NMAXY,NMAXZ ),
C    *                  TBUF1( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      COMMON / GRID2 / DX, DY, DZ
C     COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
C     COMMON / OFC / VCOU( NMAX,NMAX,NMAX )
      COMMON / OFC / VCOU( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      DIMENSION STAT(MPI_STATUS_SIZE)

      DIMENSION CRD(   NNUC,NDIM )
C     DIMENSION TVCOU( NMAXX,NMAXY,NMAXZ )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

C     IF(MYID.EQ.11 ) THEN
C     OPEN(77,FILE='VCOU.DAT',STATUS='UNKNOWN')
C     ENDIF

C----- read -----

C     IF(MYID.EQ.0 ) THEN

C     READ(23,*)
C     READ(23,*)
C     READ(23,*) N0,X0,Y0,Z0
C     READ(23,*) N1,X1,Y1,Z1
C     READ(23,*) N2,X2,Y2,Z2
C     READ(23,*) N3,X3,Y3,Z3

C     DO 10  L=1,NNUC
C     READ(23,*)
C  10 CONTINUE 

C     DO 30 I=1,N1
C     DO 30 J=1,N2
C       READ(23,'(6E13.5)') (TVCOU(I,J,K),K=1,N3)
C  30 CONTINUE

C-----------------------------------------------------------

      VCOU(:,:,:) = 0.D0

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1
      RCUT2   = RCUT*RCUT
      SQAL    = DSQRT(ALPHA)

      DO 40 K = 1, JZ
        RZ1 = DZ*( K-NNZ ) 
      DO 40 J = 1, JY
        RY1 = DY*( J-NNY ) 
      DO 40 I = 1, JX
        RX1 = DX*( I-NNX ) 
      DO 40 NA = 1, NNUC

        IF( ZVAL(NA) .EQ. 6.D0) THEN
          DEV = DELTA
        ELSE
          DEV = -0.5*DELTA
        ENDIF

        RX = RX1 - CRD(NA,1) 
        RY = RY1 - CRD(NA,2)
        RZ = RZ1 - CRD(NA,3)
        R2 = RX**2 + RY**2 + RZ**2
        R  = DSQRT(R2)
        RR = 1/R
C       VH = (ZVAL(NA)+DEV)*DERF(SQAL*R)*RR         ! error function for all atoms
        VH = ZVAL(NA)*DERF(SQAL*R)*RR               ! error function for all atoms
        VCOU(I,J,K) = VCOU(I,J,K) + VH 

40    CONTINUE

C     IF(MYID .EQ. 11) THEN
C       WRITE(*,*) MX,MY,MZ,JX,JY,JZ
C       DO 30 I=1,JX
C       DO 30 J=1,JZ
C         WRITE(77,*) VCOU(I,16,J)
C         WRITE(77,'(6E13.5)') ((VCOU(I,J,K)),K=1,JZ)
C  30   CONTINUE
C     ENDIF

      RETURN
      END



C--------------------------------------------------------
C
C   SUBROUTINE Modified Double Grid 
C   employing fourth order Lagrange interpolation
C
C   Phys. Rev. Lett. 82, 5016( 1999 )
C   
C--------------------------------------------------------

      SUBROUTINE DG4_MOD(VNUC)                  ! for non-periodic system

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
      real*4 tim,ta(2)

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NNUC = 36 )
C     PARAMETER ( NDIM =  3 )

      PARAMETER ( NP   = 4     )
      PARAMETER ( RCUT = 3.5D0 )
      PARAMETER ( NSL  = 9000  )

      PARAMETER ( PI    = 3.14159265358979323D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ
      COMMON / ACELL / XL, YL, ZL
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)
      COMMON / DBLG  / WNLOCS( NNUC,NSL ),WNLOCPX( NNUC,NSL ),
     *                 WNLOCPY(NNUC,NSL ),WNLOCPZ( NNUC,NSL )
      COMMON / PLOC / VLOC( NNUC,NSL )
      COMMON / BHS / C1(100),C2(100),AL1(100),AL2(100)

      DIMENSION VNUC( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      DIMENSION XA(  NP )
      DIMENSION YS(  NP )
      DIMENSION YP(  NP )
      DIMENSION YLD( NP )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) 'Double Grid(4th-order Lagrange Interpolation): Start'
      STIME = MPI_WTIME()
      ENDIF

      RT3  = 1.D0/3.D0
      
C----- initialize -----

      DO J=1, NSL 
      DO I=1, NNUC  

        WNLOCS( I,J) = 0.D0
        WNLOCPX(I,J) = 0.D0
        WNLOCPY(I,J) = 0.D0
        WNLOCPZ(I,J) = 0.D0

        VLOC(I,J) = 0.D0

      ENDDO 
      ENDDO 

      DO K=1,JZ  
      DO J=1,JY 
      DO I=1,JX 

        VNUC( I,J,K ) = 0.D0                                       

      ENDDO 
      ENDDO 
      ENDDO 

C----- compute weight factor for each coarse grid -----  

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO 100 NA = 1, NNUC                  ! loop over atoms
        NUMZ = INT(  ZA( NA ) )
        RLN  = DLOG( RAD( NUMZ,1 ) )
        NCOUNT = 0
      DO 10 K=1,JZ
        DO 20 J=1,JY
          DO 30 I=1,JX

          RX = DX*( I-NNX ) - PNUC(NA,1) 
          RY = DY*( J-NNY ) - PNUC(NA,2)
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)
          R2 = RX**2 + RY**2 + RZ**2
          R  = DSQRT(R2) 

          IF( R .LT. RCUT ) THEN
          NCOUNT = NCOUNT + 1
           IF( NCOUNT .GT. NSL ) THEN
           WRITE(*,*) ' dg4.f : size of dimension too small '
           STOP
           ENDIF

          DO 40 KD = -2*NDEN,2*NDEN
          DO 50 JD = -2*NDEN,2*NDEN
          DO 60 ID = -2*NDEN,2*NDEN

          DID = DABS(DBLE(ID)/DBLE(NDEN))         
          DJD = DABS(DBLE(JD)/DBLE(NDEN))         
          DKD = DABS(DBLE(KD)/DBLE(NDEN))         

          RXX = ( RX + DX*DBLE(ID)/DBLE(NDEN)) 
          RYY = ( RY + DY*DBLE(JD)/DBLE(NDEN)) 
          RZZ = ( RZ + DZ*DBLE(KD)/DBLE(NDEN)) 
          RD2 = RXX**2 + RYY**2 + RZZ**2
          RD  = DSQRT( RD2 )

          IF( RD .LT. RAD( NUMZ,3 ) ) THEN
            
            DO IPOL = 1,NP
              XA( IPOL) = RAD( NUMZ,IPOL)
              YS( IPOL) = VNLS(NUMZ,IPOL)
              YP( IPOL) = VNLP(NUMZ,IPOL)
              YLD(IPOL) = VLDC(NUMZ,IPOL)
            ENDDO

          ELSE
           
            NRAD = INT( ( 0.5D0*DLOG( RD2 ) - RLN )  
     *           / DRLOG(NUMZ) ) + 1 

            IT = NRAD - 2
            DO IPOL = 1,NP
              XA( IPOL) = RAD( NUMZ,IT+IPOL)
              YS( IPOL) = VNLS(NUMZ,IT+IPOL)
              YP( IPOL) = VNLP(NUMZ,IT+IPOL)
              YLD(IPOL) = VLDC(NUMZ,IT+IPOL)
            ENDDO

          ENDIF

          CALL POLINT(XA,YS, NP,RD,Y1,DELY)
          CALL POLINT(XA,YP, NP,RD,Y2,DELY)
          CALL POLINT(XA,YLD,NP,RD,Y3,DELY)

          VNLOCS = Y1
          VNLOCP = Y2
          VLOCD2 = Y3

          IF( (IABS(ID) .LE. NDEN)   .AND.                        ! 1
     *        (IABS(JD) .LE. NDEN)   .AND.
     *        (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 2
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 3
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 4
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 5
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 6
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 7
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 8
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ENDIF

C--------------------------
C      non-local part      
C--------------------------

C----- for s-component -----

          WNLOCS( NA,NCOUNT )  = WNLOCS( NA,NCOUNT ) 
     *                         + WIJ*VNLOCS 

C----- for p-component -----

          IF( RD .LT. 1.0D-06 ) THEN

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 

C         NW = NW + 1

          ELSE

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RXX/RD 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RYY/RD
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RZZ/RD

          ENDIF

C----------------------------------------------
C      local part ( BHS + local-d )
C----------------------------------------------

          IF( RD .LT. 1.0D-06 ) THEN

            VT = -ZVAL(NA)*( C1(NUMZ)*2.D0*DSQRT(AL1(NUMZ))                        ! analytical function ( limit zero ) of BHS
     *         +             C2(NUMZ)*2.D0*DSQRT(AL2(NUMZ)) )  
     *         /  DSQRT(PI)           
            
          ELSE

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*RD2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*RD2) ) ) 
     *         /  RD           
            
          ENDIF

            VLOC( NA,NCOUNT ) = VLOC( NA,NCOUNT )
     *                        + WIJ*(VLOCD2+VT)    

60        CONTINUE
50        CONTINUE
40        CONTINUE

            VNUC( I,J,K ) = VNUC( I,J,K ) 
     *                          + VLOC( NA,NCOUNT )

          ELSE

C----------------------------------------------
C      local part ( BHS )
C----------------------------------------------

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*R2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*R2) ) ) 
     *         /  R           
            
            VNUC( I,J,K ) = VNUC( I,J,K ) + VT                       

          ENDIF

30        CONTINUE
20      CONTINUE
10    CONTINUE

100   CONTINUE

      write(*,*) 'done! rank =',MYID

      CALL MPI_REDUCE(NCOUNT,NCOUNT1,1,MPI_INTEGER,MPI_SUM,0,
     *                MPI_COMM_WORLD,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) '  ncount = ', NCOUNT1
      ETIME = MPI_WTIME()
      write(*,*) 'Elapsed Time (dg4) = ',ETIME-STIME,' scnds'
      ENDIF

C     write(*,*) 'NR = ',NR 
C     write(*,*) 'NW = ',NW 

      RETURN
      END


C--------------------------------------------------------
C
C   SUBROUTINE Modified Double Grid 
C   employing bilinear interpolation
C
C   Phys. Rev. Lett. 82, 5016( 1999 )
C   
C--------------------------------------------------------

      SUBROUTINE DG4_MOD2(VNUC)                  ! for non-periodic system

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
      real*4 tim,ta(2)

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NNUC = 36 )
C     PARAMETER ( NDIM =  3 )

      PARAMETER ( NP   = 4     )
      PARAMETER ( RCUT = 3.5D0 )
      PARAMETER ( NSL  = 9000  )

      PARAMETER ( PI    = 3.14159265358979323D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ
      COMMON / ACELL / XL, YL, ZL
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)
      COMMON / DBLG  / WNLOCS( NNUC,NSL ),WNLOCPX( NNUC,NSL ),
     *                 WNLOCPY(NNUC,NSL ),WNLOCPZ( NNUC,NSL )
      COMMON / PLOC / VLOC( NNUC,NSL )
      COMMON / BHS / C1(100),C2(100),AL1(100),AL2(100)

      DIMENSION VNUC( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      DIMENSION XA(  NP )
      DIMENSION YS(  NP )
      DIMENSION YP(  NP )
      DIMENSION YLD( NP )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) 'Double Grid(4th-order Lagrange Interpolation): Start'
      STIME = MPI_WTIME()
      ENDIF

      RT3  = 1.D0/3.D0
      
C----- initialize -----

      DO J=1, NSL 
      DO I=1, NNUC  

        WNLOCS( I,J) = 0.D0
        WNLOCPX(I,J) = 0.D0
        WNLOCPY(I,J) = 0.D0
        WNLOCPZ(I,J) = 0.D0

        VLOC(I,J) = 0.D0

      ENDDO 
      ENDDO 

      DO K=1,JZ  
      DO J=1,JY 
      DO I=1,JX 

        VNUC( I,J,K ) = 0.D0                                       

      ENDDO 
      ENDDO 
      ENDDO 

C----- compute weight factor for each coarse grid -----  

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO 100 NA = 1, NNUC                  ! loop over atoms
        NUMZ = INT(  ZA( NA ) )
        RLN  = DLOG( RAD( NUMZ,1 ) )
        NCOUNT = 0
      DO 10 K=1,JZ
        DO 20 J=1,JY
          DO 30 I=1,JX

          RX = DX*( I-NNX ) - PNUC(NA,1) 
          RY = DY*( J-NNY ) - PNUC(NA,2)
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)
          R2 = RX**2 + RY**2 + RZ**2
          R  = DSQRT(R2) 

          IF( R .LT. RCUT ) THEN
          NCOUNT = NCOUNT + 1
           IF( NCOUNT .GT. NSL ) THEN
           WRITE(*,*) ' dg4.f : size of dimension too small '
           STOP
           ENDIF

          DO 40 KD = -NDEN,NDEN             ! double grid with bilinear interpolations
          DO 50 JD = -NDEN,NDEN             ! 2013.04.27 takahashi
          DO 60 ID = -NDEN,NDEN

          DI = DBLE(ID)/DBLE(NDEN)   ! parameter s in the note      
          DJ = DBLE(JD)/DBLE(NDEN)   ! parameter t in the note     
          DK = DBLE(KD)/DBLE(NDEN)   ! parameter r in the note     

          RXX = ( RX + DX*DI ) 
          RYY = ( RY + DY*DJ ) 
          RZZ = ( RZ + DZ*DK ) 
          RD2 = RXX**2 + RYY**2 + RZZ**2
          RD  = DSQRT( RD2 )

          IF( RD .LT. RAD( NUMZ,3 ) ) THEN       ! '3' is the initial number.
            
            DO IPOL = 1,NP
              XA( IPOL) = RAD( NUMZ,IPOL)
              YS( IPOL) = VNLS(NUMZ,IPOL)
              YP( IPOL) = VNLP(NUMZ,IPOL)
              YLD(IPOL) = VLDC(NUMZ,IPOL)
            ENDDO

          ELSE
           
            NRAD = INT( ( 0.5D0*DLOG( RD2 ) - RLN )  
     *           / DRLOG(NUMZ) ) + 1 

            IT = NRAD - 2
            DO IPOL = 1,NP
              XA( IPOL) = RAD( NUMZ,IT+IPOL)
              YS( IPOL) = VNLS(NUMZ,IT+IPOL)
              YP( IPOL) = VNLP(NUMZ,IT+IPOL)
              YLD(IPOL) = VLDC(NUMZ,IT+IPOL)
            ENDDO

          ENDIF

          CALL POLINT(XA,YS, NP,RD,Y1,DELY)
          CALL POLINT(XA,YP, NP,RD,Y2,DELY)
          CALL POLINT(XA,YLD,NP,RD,Y3,DELY)

          VNLOCS = Y1
          VNLOCP = Y2
          VLOCD2 = Y3

          DID = DABS(DI)  
          DJD = DABS(DJ) 
          DKD = DABS(DK) 

          WIJ = (1.D0-DID)*(1.D0-DJD)*(1.D0-DKD)         ! This expression is common to eight cases.
     *        / DBLE(NDEN**3)

C         IF( (IABS(ID) .LE. NDEN)   .AND.                        ! 1
C    *        (IABS(JD) .LE. NDEN)   .AND.
C    *        (IABS(KD) .LE. NDEN) ) THEN 
C         WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
C    *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
C    *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
C    *           / DBLE(NDEN**3)

C         ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 2
C    *            (IABS(JD) .LE. NDEN)   .AND.
C    *            (IABS(KD) .LE. NDEN) ) THEN 
C         WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
C    *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
C    *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
C    *           / DBLE(NDEN**3)

C         ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 3
C    *            (IABS(JD) .GT. NDEN)   .AND.
C    *            (IABS(KD) .LE. NDEN) ) THEN 
C         WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
C    *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
C    *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
C    *           / DBLE(NDEN**3)

C         ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 4
C    *            (IABS(JD) .LE. NDEN)   .AND.
C    *            (IABS(KD) .GT. NDEN) ) THEN 
C         WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
C    *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
C    *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
C    *           / DBLE(NDEN**3)

C         ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 5
C    *            (IABS(JD) .GT. NDEN)   .AND.
C    *            (IABS(KD) .GT. NDEN) ) THEN 
C         WIJ  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
C    *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
C    *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
C    *           / DBLE(NDEN**3)

C         ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 6
C    *            (IABS(JD) .LE. NDEN)   .AND.
C    *            (IABS(KD) .GT. NDEN) ) THEN 
C         WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
C    *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
C    *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
C    *           / DBLE(NDEN**3)

C         ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 7
C    *            (IABS(JD) .GT. NDEN)   .AND.
C    *            (IABS(KD) .LE. NDEN) ) THEN 
C         WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
C    *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
C    *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
C    *           / DBLE(NDEN**3)

C         ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 8
C    *            (IABS(JD) .GT. NDEN)   .AND.
C    *            (IABS(KD) .GT. NDEN) ) THEN 
C         WIJ  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
C    *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
C    *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
C    *           / DBLE(NDEN**3)

C         ENDIF

C--------------------------
C      non-local part      
C--------------------------

C----- for s-component -----

          WNLOCS( NA,NCOUNT )  = WNLOCS( NA,NCOUNT ) 
     *                         + WIJ*VNLOCS 

C----- for p-component -----

          IF( RD .LT. 1.0D-06 ) THEN

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 

C         NW = NW + 1

          ELSE

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RXX/RD 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RYY/RD
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RZZ/RD

          ENDIF

C----------------------------------------------
C      local part ( BHS + local-d )
C----------------------------------------------

          IF( RD .LT. 1.0D-06 ) THEN

            VT = -ZVAL(NA)*( C1(NUMZ)*2.D0*DSQRT(AL1(NUMZ))                        ! analytical function ( limit zero ) of BHS
     *         +             C2(NUMZ)*2.D0*DSQRT(AL2(NUMZ)) )  
     *         /  DSQRT(PI)           
            
          ELSE

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*RD2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*RD2) ) ) 
     *         /  RD           
            
          ENDIF

            VLOC( NA,NCOUNT ) = VLOC( NA,NCOUNT )
     *                        + WIJ*(VLOCD2+VT)    

60        CONTINUE
50        CONTINUE
40        CONTINUE

            VNUC( I,J,K ) = VNUC( I,J,K ) 
     *                          + VLOC( NA,NCOUNT )

          ELSE

C----------------------------------------------
C      local part ( BHS )
C----------------------------------------------

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*R2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*R2) ) ) 
     *         /  R           
            
            VNUC( I,J,K ) = VNUC( I,J,K ) + VT                       

          ENDIF

30        CONTINUE
20      CONTINUE
10    CONTINUE

100   CONTINUE

      write(*,*) 'done! rank =',MYID

      CALL MPI_REDUCE(NCOUNT,NCOUNT1,1,MPI_INTEGER,MPI_SUM,0,
     *                MPI_COMM_WORLD,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) '  ncount = ', NCOUNT1
      ETIME = MPI_WTIME()
      write(*,*) 'Elapsed Time (dg4) = ',ETIME-STIME,' scnds'
      ENDIF

C     write(*,*) 'NR = ',NR 
C     write(*,*) 'NW = ',NW 

      RETURN
      END

C--------------------------------------------------------
C
C   SUBROUTINE Modified Double Grid 
C   employing cubic interpolation
C
C   Phys. Rev. Lett. 82, 5016( 1999 )
C   
C--------------------------------------------------------

      SUBROUTINE DG4_CUB(VNUC)                  ! for non-periodic system

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
      real*4 tim,ta(2)

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol
      include 'nlocd.i'                         ! Vmol   non-local d
!     include 'sizes.i'                         ! Vmol   non-local d

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NNUC = 36 )
C     PARAMETER ( NDIM =  3 )

      PARAMETER ( NP   = 4     )
C     PARAMETER ( RCUT = 2.5D0 )
C     PARAMETER ( NSL  = 3000  )

      PARAMETER ( PI    = 3.14159265358979323D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ
      COMMON / ACELL / XL, YL, ZL
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)
      COMMON / DBLG  / WNLOCS( NNUC,NSL ),WNLOCPX( NNUC,NSL ),
     *                 WNLOCPY(NNUC,NSL ),WNLOCPZ( NNUC,NSL )
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PLOC / VLOC( NNUC,NSL )
      COMMON / BHS / C1(100),C2(100),AL1(100),AL2(100)

      DIMENSION VNUC( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      DIMENSION XA(  NP )
      DIMENSION YS(  NP )
      DIMENSION YP(  NP )
      DIMENSION YD(  NP )       ! non-local d
      DIMENSION YLD( NP )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) 'Double Grid (Cubic Interpolation): Start'
      STIME = MPI_WTIME()

      WRITE(*,*)                             
      WRITE(*,*) '    QM Coordinates'
      DO J=1,NNUC                                ! write coordinates
      WRITE(*,99) ZA(J),ZVAL(J),
     *            PNUC(J,1),PNUC(J,2),PNUC(J,3)
      ENDDO

      ENDIF

99    FORMAT( F6.1,F6.1,3F12.6 )

      RT3   = 1.D0/3.D0
      RCUT2 = RCUT*RCUT
      ne = NNUC - NLINK
      
C----- initialize -----

      DO J=1, NSL 
      DO I=1, ne  
C     DO I=1, NNUC  

        WNLOCS( I,J) = 0.D0
        WNLOCPX(I,J) = 0.D0
        WNLOCPY(I,J) = 0.D0
        WNLOCPZ(I,J) = 0.D0

        VLOC(I,J) = 0.D0

      ENDDO 
      ENDDO 

      WNLOCDXY(:,:) = 0.D0
      WNLOCDYZ(:,:) = 0.D0
      WNLOCDZX(:,:) = 0.D0
      WNLOCDZ2(:,:) = 0.D0
      WNLOCDX2(:,:) = 0.D0

      DO K=1,JZ  
      DO J=1,JY 
      DO I=1,JX 

        VNUC( I,J,K ) = 0.D0                                       

      ENDDO 
      ENDDO 
      ENDDO 

C----- compute weight factor for each coarse grid -----  

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

C     DO 100 NA = 1, NNUC                  ! loop over atoms
      DO 100 NA = 1, ne                    ! loop over atoms
        NUMZ = INT(  ZA( NA ) )
        RLN  = DLOG( RAD( NUMZ,1 ) )
        NCOUNT = 0
        NJ = NPD(NA)

      SELECT CASE (NJ)

      CASE ( 0 )    ! NJ = 0: for s and p non-local

      DO 10 K=1,JZ
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)
        DO 20 J=1,JY
          RY = DY*( J-NNY ) - PNUC(NA,2)
C         RT2 = RZ**2 + RY**2
C         IF( RT2 >= RCUT2 ) CYCLE
          DO 30 I=1,JX

          RX = DX*( I-NNX ) - PNUC(NA,1) 
C         R2 = RX**2 + RT2
          R2 = RX**2 + RY**2 + RZ**2
          R  = DSQRT(R2) 

          IF( R < RCUT ) THEN
          NCOUNT = NCOUNT + 1
           IF( NCOUNT .GT. NSL ) THEN
           WRITE(*,*) ' dg4.f : size of dimension too small '
           STOP
           ENDIF

          DO 40 KD = -2*NDEN,2*NDEN
          RZZ = ( RZ + DZ*DBLE(KD)/DBLE(NDEN)) 
          DKD = DABS(DBLE(KD)/DBLE(NDEN))         
          DO 50 JD = -2*NDEN,2*NDEN
          RYY = ( RY + DY*DBLE(JD)/DBLE(NDEN)) 
          DJD = DABS(DBLE(JD)/DBLE(NDEN))         
          RDT2 = RYY**2 + RZZ**2
          DO 60 ID = -2*NDEN,2*NDEN

          DID = DABS(DBLE(ID)/DBLE(NDEN))         

          RXX = ( RX + DX*DBLE(ID)/DBLE(NDEN)) 
          RD2 = RXX**2 + RYY**2 + RZZ**2
C         RD2 = RXX**2 + RDT2 
          RD  = DSQRT( RD2 )

          IF( RD .LT. RAD( NUMZ,3 ) ) THEN
            
            DO IPOL = 1,NP
              XA( IPOL) = RAD( NUMZ,IPOL)
              YS( IPOL) = VNLS(NUMZ,IPOL)
              YP( IPOL) = VNLP(NUMZ,IPOL)
              YLD(IPOL) = VLDC(NUMZ,IPOL)
            ENDDO

          ELSE
           
            NRAD = INT( ( 0.5D0*DLOG( RD2 ) - RLN )  
     *           / DRLOG(NUMZ) ) + 1 

            IT = NRAD - 2
            DO IPOL = 1,NP
              XA( IPOL) = RAD( NUMZ,IT+IPOL)
              YS( IPOL) = VNLS(NUMZ,IT+IPOL)
              YP( IPOL) = VNLP(NUMZ,IT+IPOL)
              YLD(IPOL) = VLDC(NUMZ,IT+IPOL)
            ENDDO

          ENDIF

          CALL POLINT(XA,YS, NP,RD,Y1,DELY)
          CALL POLINT(XA,YP, NP,RD,Y2,DELY)
          CALL POLINT(XA,YLD,NP,RD,Y3,DELY)

          VNLOCS = Y1
          VNLOCP = Y2
          VLOCD2 = Y3

          IF( (IABS(ID) .LE. NDEN)   .AND.                        ! 1
     *        (IABS(JD) .LE. NDEN)   .AND.
     *        (IABS(KD) .LE. NDEN) ) THEN 

          AI  = DID - 1.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 1.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 1.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID 
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD 
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD 
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = ( 2.D0*AI3 - 0.5D0*BI3 + 3.D0*AI2 + 0.5D0*BI2 )
     *        * ( 2.D0*AJ3 - 0.5D0*BJ3 + 3.D0*AJ2 + 0.5D0*BJ2 )
     *        * ( 2.D0*AK3 - 0.5D0*BK3 + 3.D0*AK2 + 0.5D0*BK2 )
     *        /   DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 2
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 

          AI  = DID - 2.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 1.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 1.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID - 1.D0 
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD 
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD 
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = ( -0.5D0*AI3 - 0.5D0*AI2 )
     *        * (   2.D0*AJ3 - 0.5D0*BJ3 + 3.D0*AJ2 + 0.5D0*BJ2 )
     *        * (   2.D0*AK3 - 0.5D0*BK3 + 3.D0*AK2 + 0.5D0*BK2 )
     *        /     DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 3
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 

          AI  = DID - 1.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 2.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 1.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID 
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD - 1.D0
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD 
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = (   2.D0*AI3 - 0.5D0*BI3 + 3.D0*AI2 + 0.5D0*BI2 )
     *        * ( -0.5D0*AJ3 - 0.5D0*AJ2 )
     *        * (   2.D0*AK3 - 0.5D0*BK3 + 3.D0*AK2 + 0.5D0*BK2 )
     *        /     DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 4
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 

          AI  = DID - 1.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 1.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 2.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID 
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD 
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD - 1.D0
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = (   2.D0*AI3 - 0.5D0*BI3 + 3.D0*AI2 + 0.5D0*BI2 )
     *        * (   2.D0*AJ3 - 0.5D0*BJ3 + 3.D0*AJ2 + 0.5D0*BJ2 )
     *        * ( -0.5D0*AK3 - 0.5D0*AK2 )
     *        /     DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 5
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 

          AI  = DID - 1.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 2.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 2.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID 
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD - 1.D0
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD - 1.D0
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = (   2.D0*AI3 - 0.5D0*BI3 + 3.D0*AI2 + 0.5D0*BI2 )
     *        * ( -0.5D0*AJ3 - 0.5D0*AJ2 )
     *        * ( -0.5D0*AK3 - 0.5D0*AK2 )
     *        /     DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 6
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 

          AI  = DID - 2.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 1.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 2.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID - 1.D0
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD 
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD - 1.D0
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = ( -0.5D0*AI3 - 0.5D0*AI2 )
     *        * (   2.D0*AJ3 - 0.5D0*BJ3 + 3.D0*AJ2 + 0.5D0*BJ2 )
     *        * ( -0.5D0*AK3 - 0.5D0*AK2 )
     *        /     DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 7
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 

          AI  = DID - 2.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 2.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 1.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID - 1.D0
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD - 1.D0
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD 
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = ( -0.5D0*AI3 - 0.5D0*AI2 )
     *        * ( -0.5D0*AJ3 - 0.5D0*AJ2 )
     *        * (   2.D0*AK3 - 0.5D0*BK3 + 3.D0*AK2 + 0.5D0*BK2 )
     *        /     DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 8
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 

          AI  = DID - 2.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 2.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 2.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID - 1.D0
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD - 1.D0
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD - 1.D0
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = ( -0.5D0*AI3 - 0.5D0*AI2 )
     *        * ( -0.5D0*AJ3 - 0.5D0*AJ2 )
     *        * ( -0.5D0*AK3 - 0.5D0*AK2 )
     *        /     DBLE(NDEN**3)

          ENDIF

C--------------------------
C      non-local part      
C--------------------------

C----- for s-component -----

          WNLOCS( NA,NCOUNT )  = WNLOCS( NA,NCOUNT ) 
     *                         + WIJ*VNLOCS 

C----- for p-component -----

          IF( RD .LT. 1.0D-06 ) THEN

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 

C         NW = NW + 1

          ELSE

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RXX/RD 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RYY/RD
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 
     *                         + WIJ*VNLOCP*RZZ/RD

          ENDIF

C----------------------------------------------
C      local part ( BHS + local-d )
C----------------------------------------------

          IF( RD .LT. 1.0D-06 ) THEN

            VT = -ZVAL(NA)*( C1(NUMZ)*2.D0*DSQRT(AL1(NUMZ))                        ! analytical function ( limit zero ) of BHS
     *         +             C2(NUMZ)*2.D0*DSQRT(AL2(NUMZ)) )  
     *         /  DSQRT(PI)           
            
          ELSE

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*RD2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*RD2) ) ) 
     *         /  RD           
            
          ENDIF

            VLOC( NA,NCOUNT ) = VLOC( NA,NCOUNT )
     *                        + WIJ*(VLOCD2+VT)    

60        CONTINUE
50        CONTINUE
40        CONTINUE

            VNUC( I,J,K ) = VNUC( I,J,K ) 
     *                          + VLOC( NA,NCOUNT )

          ELSE

C----------------------------------------------
C      local part ( BHS )
C----------------------------------------------

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*R2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*R2) ) ) 
     *         /  R           
            
            VNUC( I,J,K ) = VNUC( I,J,K ) + VT                       

          ENDIF

30        CONTINUE
20      CONTINUE
10    CONTINUE

C     IF( MYID .EQ. 10) WRITE(*,*) NA,AL1(NUMZ),AL2(NUMZ),R,VT,VNUC(1,1,1) 

      CASE ( 1 )    ! NJ = 1: for s, p, and d non-local

      DO 15 K=1,JZ
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)
        DO 25 J=1,JY
          RY = DY*( J-NNY ) - PNUC(NA,2)
          RT2 = RZ**2 + RY**2
          DO 35 I=1,JX

          RX = DX*( I-NNX ) - PNUC(NA,1) 
C         R2 = RX**2 + RT2
          R2 = RX**2 + RY**2 + RZ**2
          R  = DSQRT(R2) 

          IF( R < RCUT ) THEN
          NCOUNT = NCOUNT + 1
           IF( NCOUNT .GT. NSL ) THEN
           WRITE(*,*) ' dg4.f : size of dimension too small '
           STOP
           ENDIF

          DO 45 KD = -2*NDEN,2*NDEN
          RZZ = ( RZ + DZ*DBLE(KD)/DBLE(NDEN)) 
          DKD = DABS(DBLE(KD)/DBLE(NDEN))         
          DO 55 JD = -2*NDEN,2*NDEN
          RYY = ( RY + DY*DBLE(JD)/DBLE(NDEN)) 
          DJD = DABS(DBLE(JD)/DBLE(NDEN))         
C         RDT2 = RYY**2 + RZZ**2
          DO 65 ID = -2*NDEN,2*NDEN

          DID = DABS(DBLE(ID)/DBLE(NDEN))         

          RXX = ( RX + DX*DBLE(ID)/DBLE(NDEN)) 
          RD2 = RXX**2 + RYY**2 + RZZ**2
C         RD2 = RXX**2 + RDT2 
          RD  = DSQRT( RD2 )

          IF( RD .LT. RAD( NUMZ,3 ) ) THEN
            
            DO IPOL = 1,NP
              XA( IPOL) = RAD( NUMZ,IPOL)
              YS( IPOL) = VNLS(NUMZ,IPOL)
              YP( IPOL) = VNLP(NUMZ,IPOL)
              YD( IPOL) = VNLD(NUMZ,IPOL)
            ENDDO

          ELSE
           
            NRAD = INT( ( 0.5D0*DLOG( RD2 ) - RLN )  
     *           / DRLOG(NUMZ) ) + 1 

            IT = NRAD - 2
            DO IPOL = 1,NP
              XA( IPOL) = RAD( NUMZ,IT+IPOL)
              YS( IPOL) = VNLS(NUMZ,IT+IPOL)
              YP( IPOL) = VNLP(NUMZ,IT+IPOL)
              YD( IPOL) = VNLD(NUMZ,IT+IPOL)
            ENDDO

          ENDIF

          CALL POLINT(XA,YS, NP,RD,Y1,DELY)
          CALL POLINT(XA,YP, NP,RD,Y2,DELY)
          CALL POLINT(XA,YD, NP,RD,Y3,DELY)

          VNLOCS = Y1
          VNLOCP = Y2
          VNLOCD = Y3

          IF( (IABS(ID) .LE. NDEN)   .AND.                        ! 1
     *        (IABS(JD) .LE. NDEN)   .AND.
     *        (IABS(KD) .LE. NDEN) ) THEN 

          AI  = DID - 1.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 1.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 1.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID 
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD 
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD 
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = ( 2.D0*AI3 - 0.5D0*BI3 + 3.D0*AI2 + 0.5D0*BI2 )
     *        * ( 2.D0*AJ3 - 0.5D0*BJ3 + 3.D0*AJ2 + 0.5D0*BJ2 )
     *        * ( 2.D0*AK3 - 0.5D0*BK3 + 3.D0*AK2 + 0.5D0*BK2 )
     *        /   DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 2
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 

          AI  = DID - 2.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 1.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 1.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID - 1.D0 
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD 
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD 
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = ( -0.5D0*AI3 - 0.5D0*AI2 )
     *        * (   2.D0*AJ3 - 0.5D0*BJ3 + 3.D0*AJ2 + 0.5D0*BJ2 )
     *        * (   2.D0*AK3 - 0.5D0*BK3 + 3.D0*AK2 + 0.5D0*BK2 )
     *        /     DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 3
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 

          AI  = DID - 1.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 2.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 1.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID 
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD - 1.D0
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD 
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = (   2.D0*AI3 - 0.5D0*BI3 + 3.D0*AI2 + 0.5D0*BI2 )
     *        * ( -0.5D0*AJ3 - 0.5D0*AJ2 )
     *        * (   2.D0*AK3 - 0.5D0*BK3 + 3.D0*AK2 + 0.5D0*BK2 )
     *        /     DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 4
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 

          AI  = DID - 1.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 1.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 2.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID 
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD 
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD - 1.D0
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = (   2.D0*AI3 - 0.5D0*BI3 + 3.D0*AI2 + 0.5D0*BI2 )
     *        * (   2.D0*AJ3 - 0.5D0*BJ3 + 3.D0*AJ2 + 0.5D0*BJ2 )
     *        * ( -0.5D0*AK3 - 0.5D0*AK2 )
     *        /     DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 5
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 

          AI  = DID - 1.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 2.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 2.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID 
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD - 1.D0
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD - 1.D0
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = (   2.D0*AI3 - 0.5D0*BI3 + 3.D0*AI2 + 0.5D0*BI2 )
     *        * ( -0.5D0*AJ3 - 0.5D0*AJ2 )
     *        * ( -0.5D0*AK3 - 0.5D0*AK2 )
     *        /     DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 6
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 

          AI  = DID - 2.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 1.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 2.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID - 1.D0
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD 
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD - 1.D0
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = ( -0.5D0*AI3 - 0.5D0*AI2 )
     *        * (   2.D0*AJ3 - 0.5D0*BJ3 + 3.D0*AJ2 + 0.5D0*BJ2 )
     *        * ( -0.5D0*AK3 - 0.5D0*AK2 )
     *        /     DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 7
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 

          AI  = DID - 2.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 2.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 1.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID - 1.D0
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD - 1.D0
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD 
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = ( -0.5D0*AI3 - 0.5D0*AI2 )
     *        * ( -0.5D0*AJ3 - 0.5D0*AJ2 )
     *        * (   2.D0*AK3 - 0.5D0*BK3 + 3.D0*AK2 + 0.5D0*BK2 )
     *        /     DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 8
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 

          AI  = DID - 2.D0
          AI2 = AI *AI
          AI3 = AI2*AI
          AJ  = DJD - 2.D0
          AJ2 = AJ *AJ
          AJ3 = AJ2*AJ
          AK  = DKD - 2.D0
          AK2 = AK *AK
          AK3 = AK2*AK

          BI  = DID - 1.D0
          BI2 = BI *BI
          BI3 = BI2*BI
          BJ  = DJD - 1.D0
          BJ2 = BJ *BJ
          BJ3 = BJ2*BJ
          BK  = DKD - 1.D0
          BK2 = BK *BK
          BK3 = BK2*BK

          WIJ = ( -0.5D0*AI3 - 0.5D0*AI2 )
     *        * ( -0.5D0*AJ3 - 0.5D0*AJ2 )
     *        * ( -0.5D0*AK3 - 0.5D0*AK2 )
     *        /     DBLE(NDEN**3)

          ENDIF

C--------------------------
C      non-local part ( s, p, and d )     
C--------------------------

          RDINV  = 1.D0/RD
          RDINV2 = RDINV*RDINV

C----- for s-component -----

          WNLOCS( NA,NCOUNT )  = WNLOCS( NA,NCOUNT ) 
     *                         + WIJ*VNLOCS 

C----- for p-component -----

          IF( RD .LT. 1.0D-06 ) THEN

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 

C         NW = NW + 1

          ELSE

          FACP = WIJ*RDINV
          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
     *                         + VNLOCP*RXX*FACP
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
     *                         + VNLOCP*RYY*FACP
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 
     *                         + VNLOCP*RZZ*FACP

          ENDIF

C----- for d-component -----

          IF( RD .LT. 1.0D-05 ) THEN

          WNLOCDXY( NA,NCOUNT ) = WNLOCDXY( NA,NCOUNT ) 
          WNLOCDYZ( NA,NCOUNT ) = WNLOCDYZ( NA,NCOUNT ) 
          WNLOCDZX( NA,NCOUNT ) = WNLOCDZX( NA,NCOUNT ) 
          WNLOCDZ2( NA,NCOUNT ) = WNLOCDZ2( NA,NCOUNT ) 
          WNLOCDX2( NA,NCOUNT ) = WNLOCDX2( NA,NCOUNT ) 

          ELSE

          FACD = WIJ*RDINV2
          WNLOCDXY( NA,NCOUNT ) = WNLOCDXY( NA,NCOUNT ) 
     *                          + VNLOCD*RXX*RYY*FACD
          WNLOCDYZ( NA,NCOUNT ) = WNLOCDYZ( NA,NCOUNT ) 
     *                          + VNLOCD*RYY*RZZ*FACD
          WNLOCDZX( NA,NCOUNT ) = WNLOCDZX( NA,NCOUNT ) 
     *                          + VNLOCD*RZZ*RXX*FACD
          WNLOCDZ2( NA,NCOUNT ) = WNLOCDZ2( NA,NCOUNT ) 
     *                          + VNLOCD*(3.D0*RZZ**2-RD2)*FACD
          WNLOCDX2( NA,NCOUNT ) = WNLOCDX2( NA,NCOUNT ) 
     *                          + VNLOCD*(RXX**2-RYY**2)*FACD

          ENDIF

C----------------------------------------------
C      local part ( BHS + local-d )
C----------------------------------------------

          IF( RD .LT. 1.0D-06 ) THEN

            VT = -ZVAL(NA)*( C1(NUMZ)*2.D0*DSQRT(AL1(NUMZ))                        ! analytical function ( limit zero ) of BHS
     *         +             C2(NUMZ)*2.D0*DSQRT(AL2(NUMZ)) )  
     *         /  DSQRT(PI)           
            
          ELSE

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*RD2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*RD2) ) ) 
     *         *  RDINV           

          ENDIF

            VLOC( NA,NCOUNT ) = VLOC( NA,NCOUNT )
     *                        + WIJ*VT    

65        CONTINUE
55        CONTINUE
45        CONTINUE

            VNUC( I,J,K ) = VNUC( I,J,K ) 
     *                    + VLOC( NA,NCOUNT )

          ELSE

C----------------------------------------------
C      local part ( BHS )
C----------------------------------------------

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*R2) )                  ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*R2) ) ) 
     *         /  R           
            
            VNUC( I,J,K ) = VNUC( I,J,K ) + VT                       

          ENDIF

35        CONTINUE
25      CONTINUE
15    CONTINUE

C     IF( MYID .EQ. 20) WRITE(*,*) NA,ZVAL(NA),C1(NUMZ),C2(NUMZ),VNUC(1,1,1) 

      END SELECT

100   CONTINUE

      write(*,*) 'done! rank =',MYID

      CALL MPI_REDUCE(NCOUNT,NCOUNT1,1,MPI_INTEGER,MPI_SUM,0,
     *                MPI_COMM_WORLD,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) '  ncount = ', NCOUNT1
      ETIME = MPI_WTIME()
      write(*,*) 'Elapsed Time (dg4) = ',ETIME-STIME,' scnds'
      ENDIF

C     write(*,*) 'NR = ',NR 
C     write(*,*) 'NW = ',NW 

C     REWIND(24)
      DO K = 1,JZ  
      DO J = 1,JY
      DO I = 1,JX
C       WRITE(24) VNUC(I,J,K)
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END

C------------------------------------------
C   SUBROUTINE DIAGONALIZE WAVEFUNCTIONS
C   BY modified GRAM-SCHMIDT METHOD
C   with OpenMP
C------------------------------------------

      SUBROUTINE MDIAG( NRD,NOR,TWF,RWF ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                            ! mpi
      include 'QMpara.i'                         ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NORA = 49 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NOCC =  1 )

      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ

      DIMENSION NRD( NOR )

      DIMENSION TWF( NMAXX/NX*NMAXY/NY*NMAXZ/NZ,NOR )
      DIMENSION RWF( NMAXX/NX*NMAXY/NY*NMAXZ/NZ,NOR )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C!$ call omp_set_num_threads(8)

      IF(MYID.EQ.0) THEN
      STIME = MPI_WTIME()
      ENDIF
C     call second(dtim0)

      DV = DX*DY*DZ
      nmaxxyz=NMAXX*NMAXY*NMAXZ/(NX*NY*NZ)

C----- modified gram schmidt -----

      Do NTMP = 1, NOR        

         NRDK=NRD( NTMP )

         If( NTMP.gt.1 ) then
            Do K = 1,NTMP-1        

               sumdot = 0.d0
!$omp parallel default(shared)
!$omp& private(N)
!$omp do reduction(+:sumdot)
               Do N = 1, nmaxxyz
                  sumdot = sumdot + RWF(N,K)*TWF(N,NRDK)
               End do
!$omp end do
!$omp single
               effk1 = sumdot*DV
!$omp end single
!$omp end parallel

      CALL MPI_REDUCE(effk1,effk,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(effk,1,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

!$omp parallel default(shared)
!$omp& private(N)
!$omp do
CDEC$ IVDEP
               Do N = 1, nmaxxyz
                  TWF( N,NRDK ) = TWF( N,NRDK ) - effk*RWF( N,K )
               End do
!$omp end do nowait
!$omp end parallel

           End do
         End if

C----- normalize -----

      ANORM = 0.D0
!$omp parallel default(shared)
!$omp& private(N)
!$omp do reduction(+:ANORM)
         Do N = 1, nmaxxyz
            ANORM = ANORM + TWF( N,NRDK )**2
         Enddo
!$omp end do
!$omp end parallel
      ANORM = ANORM*DV

      CALL MPI_REDUCE(ANORM,BNORM,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)

      IF(MYID.EQ.0) THEN
      ANORM = DSQRT( BNORM )
      ENDIF

      CALL MPI_BCAST(ANORM,1,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

      ANORMiv = 1.d0/DSQRT( ANORM )

!$omp parallel default(shared)
!$omp& private(N)
!$omp do
CDEC$ IVDEP
         Do N = 1, nmaxxyz
            TWF( N,NRDK )  = TWF( N,NRDK )*ANORMiv
            RWF( N,NTMP )  = TWF( N,NRDK )
         Enddo
!$omp end do nowait
!$omp end parallel
      Enddo

      IF(MYID.EQ.0) THEN
      ETIME = MPI_WTIME()
C     write(*,'("   mDIAG  cpu=",f15.4)') ETIME-STIME
      ENDIF
C     call second(dtim1)

      RETURN
      END

C------------------------------------------
C   SUBROUTINE DIAGONALIZE WAVEFUNCTIONS
C   BY modified GRAM-SCHMIDT METHOD
C   with OpenMP
C------------------------------------------

      SUBROUTINE MDIAG_BLK( NRD,NOR,TWF,RWF ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                            ! mpi
      include 'QMpara.i'                         ! Vmol

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NORA = 49 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NOCC =  1 )

      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ

      DIMENSION NRD( NOR )

      DIMENSION TWF( NMAXX/NX*NMAXY/NY*NMAXZ/NZ,NOR )
      DIMENSION RWF( NMAXX/NX*NMAXY/NY*NMAXZ/NZ,NOR )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C!$ call omp_set_num_threads(4)

      IF(MYID.EQ.0) THEN
      STIME = MPI_WTIME()
      ENDIF
C     call second(dtim0)

      DV = DX*DY*DZ
      nmaxxyz=NMAXX*NMAXY*NMAXZ/(NX*NY*NZ)
      lblk = 40960

C----- modified gram schmidt -----

      Do NTMP = 1, NOR        

         NRDK=NRD( NTMP )

         If( NTMP.gt.1 ) then
            Do K = 1,NTMP-1        

            sumdot = 0.d0
            Do ib = 1,nmaxxyz, lblk
               ibend = min(ib+lblk-1,nmaxxyz)
!$omp parallel default(shared)
!$omp& private(N)
!$omp do reduction(+:sumdot)
               Do N = ib, ibend
                  sumdot = sumdot + RWF(N,K)*TWF(N,NRDK)
               End do
!$omp end do
!$omp end parallel
            End do
            effk1 = sumdot*DV

      CALL MPI_REDUCE(effk1,effk,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(effk,1,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

!$omp parallel default(shared)
!$omp& private(N)
!$omp do
CDEC$ IVDEP
               Do N = 1, nmaxxyz
                  TWF( N,NRDK ) = TWF( N,NRDK ) - effk*RWF( N,K )
               End do
!$omp end do nowait
!$omp end parallel

           End do
         End if

C----- normalize -----

      ANORM = 0.D0
!$omp parallel default(shared)
!$omp& private(N)
!$omp do reduction(+:ANORM)
         Do N = 1, nmaxxyz
            ANORM = ANORM + TWF( N,NRDK )**2
         Enddo
!$omp end do
!$omp end parallel
      ANORM = ANORM*DV

      CALL MPI_REDUCE(ANORM,BNORM,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)

      IF(MYID.EQ.0) THEN
      ANORM = DSQRT( BNORM )
      ENDIF

      CALL MPI_BCAST(ANORM,1,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

      ANORMiv = 1.d0/DSQRT( ANORM )

!$omp parallel default(shared)
!$omp& private(N)
!$omp do
CDEC$ IVDEP
         Do N = 1, nmaxxyz
            TWF( N,NRDK )  = TWF( N,NRDK )*ANORMiv
            RWF( N,NTMP )  = TWF( N,NRDK )
         Enddo
!$omp end do nowait
!$omp end parallel
      Enddo

      IF(MYID.EQ.0) THEN
      ETIME = MPI_WTIME()
      write(*,'(" mDIAG  cpu=",f15.4)') ETIME-STIME
      ENDIF
C     call second(dtim1)

      RETURN
      END


C-----------------------------------------------
C
C     COMPUTE AVERAGED ELECTRON DENSITY
C
C-----------------------------------------------

      SUBROUTINE AVDNST ( RHO,IMD )      

      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4( I-N ) 

      include 'mpif.h'
      include 'mpi.i'
      include 'QMpara.i'                        ! Vmol

      PARAMETER ( NCUT = 5000 ) 

      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX,DY,DZ
C     COMMON / RST / PDPM(3),SUMRHO(NMAX,NMAX,NMAX),
C    *               TSSE,NERDF(NDXL),CHI(NDXL,NDXL)      

      DIMENSION RWF(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA)
      DIMENSION RHO(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ)
      DIMENSION AVRHO( NMAXX/NX,NMAXY/NY,NMAXZ/NZ)
      DIMENSION SUMRHO(NMAXX/NX,NMAXY/NY,NMAXZ/NZ)

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)


      DV = DX*DY*DZ
      REDV = 1.0 / DV

      IF(IMD.EQ.1) THEN

       DO 10 I = 1,JX
       DO 20 J = 1,JY
       DO 30 K = 1,JZ
         SUMRHO(I,J,K) = 0.0 
30     CONTINUE
20     CONTINUE
10     CONTINUE

      ICOUNT = 0

      ENDIF 

      IF(IMD.GT.NCUT) THEN

        ICOUNT = ICOUNT + 1

        DO 60 K = 1,JZ  
        DO 50 J = 1,JY
        DO 40 I = 1,JX

          SUMRHO(I,J,K) = SUMRHO(I,J,K) 
     *                  + RHO(I,J,K)*DV
          AVRHO(I,J,K)  = SUMRHO(I,J,K)*REDV/DBLE(ICOUNT)

40      CONTINUE
50      CONTINUE
60      CONTINUE

      IF( MOD(IMD,100) .EQ. 0) THEN
C       REWIND(24)
        DO 65 K = 1,JZ  
        DO 55 J = 1,JY
        DO 45 I = 1,JX
C         WRITE(24) RWF(I,J,K,NORA-1)
C         WRITE(24) AVRHO(I,J,K)
45      CONTINUE
55      CONTINUE
65      CONTINUE
      ENDIF

      ENDIF 

      RETURN
      END   


C--------------------------------------------------------
C    Reduce electron densities defined over 64^3 grids
C     to 32^3 ones
C--------------------------------------------------------

      SUBROUTINE REDDNS(RHO,IMD)

      IMPLICIT INTEGER*4 ( I-N )
      IMPLICIT REAL*8 ( A-H,O-Z )      
      
      include 'mpif.h'
      include 'mpi.i'
      include 'QMpara.i'                        ! Vmol

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5
      
      DIMENSION RHO( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI2 / LX(0:NX*NY*NZ-1),LY(0:NX*NY*NZ-1),
     *                 LZ(0:NX*NY*NZ-1)
      COMMON / MMPI3 / NIDX(1:2),NIDY(1:2),NIDZ(1:2)
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX, DY, DZ
      COMMON / RDDNST /RRHO(NMAXX/(2*NX),NMAXY/(2*NY),NMAXZ/(2*NZ))

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      IF ( MYID .EQ. 0 ) THEN
      NERR = 0
      IF( MOD(JX,2) .NE. 0 ) THEN
       write(*,*) 'Error in sub. reddns: jx cannot be divided by 2.' 
       NERR = NERR + 1
      ENDIF
      IF( MOD(JY,2) .NE. 0 ) THEN
       write(*,*) 'Error in sub. reddns: jy cannot be divided by 2.' 
       NERR = NERR + 1
      ENDIF
      IF( MOD(JZ,2) .NE. 0 ) THEN
       write(*,*) 'Error in sub. reddns: jz cannot be divided by 2.' 
       NERR = NERR + 1
      ENDIF
      IF( NERR .NE. 0 ) STOP
      ENDIF

      DV  = DX*DY*DZ

      JX2 = JX/2
      JY2 = JY/2
      JZ2 = JZ/2

C----- initialize -----

      SUM    = 0.D0
      PSUM   = 0.D0
      SUM1   = 0.D0
      SUMRHO = 0.D0

      RRHO(:,:,:) = 0.D0

C----- Normalize ----------------

C     IF(IMD .EQ. 0) THEN

      DO K=1,JZ
      DO J=1,JY
      DO I=1,JX
        PSUM = PSUM + RHO(I,J,K)*DV
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE(PSUM,SUM,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST( SUM,1,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

      IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
     *                 .OR.EXC.EQ.'RXalpha') THEN         ! total charge of electrons
        ZSUM = DBLE(2*MORA)
      ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
     *                     .OR.EXC.EQ.'UXalpha') THEN   
        ZSUM = DBLE(MORA+MORB)
      ENDIF

      PSUM = 0.D0
      DO K = 1,JZ
      DO J = 1,JY
      DO I = 1,JX
        RHO(I,J,K) = RHO(I,J,K)*ZSUM/SUM
        PSUM = PSUM + RHO(I,J,K)*DV
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE(PSUM,SUM,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST( SUM,1,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

      IF( MYID .EQ. 0 ) THEN
        WRITE(*,*) 'REDDNS: SUM = ',SUM
      ENDIF

C     ENDIF

C--------------------------------

       DO K = 1,JZ2
       DO J = 1,JY2
       DO I = 1,JX2
         
        RRHO(I,J,K)=RHO(2*I-1,2*J-1,2*K-1)
     *             +RHO(2*I-1,2*J-1,2*K  )
     *             +RHO(2*I-1,2*J,  2*K-1)
     *             +RHO(2*I-1,2*J,  2*K  )
     *             +RHO(2*I,  2*J-1,2*K-1)
     *             +RHO(2*I,  2*J-1,2*K  )
     *             +RHO(2*I,  2*J,  2*K-1)
     *             +RHO(2*I,  2*J,  2*K  )

        RRHO(I,J,K)=RRHO(I,J,K)/8.D0

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END


C-------------------------------------------
C     SUBROUTINE QM/MM FORCE
C-------------------------------------------

      SUBROUTINE RQMF1( SCRD,FRCS,FRCN )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
!     include 'atoms.i'

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC = 36 )
      PARAMETER ( DCUT  = 1.D-9 ) 
      PARAMETER ( ALPHA = 1.0D0 )
C     PARAMETER ( NLINK1 = 10 )
      PARAMETER ( RRCUT = 4.724D0 )
      PARAMETER ( RRCUT1= 5.669D0 )
      PARAMETER ( DRCUT = 0.945D0 )

      PARAMETER ( PI    = 3.14159265358979323D0 )
      
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PRMT4 / nion1,iiont(maxatm)
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / LJQMMM / SQMMM(NNUC,maxatm),EPQMMM(NNUC,maxatm),
     *                  chgmm(maxatm),chgqm(NNUC)

      COMMON / RDDNST /RRHO(NMAXX/(2*NX),NMAXY/(2*NY),NMAXZ/(2*NZ))

      DIMENSION SCRD( NDIM,maxatm )

      DIMENSION FDSS(  maxatm,NDIM )              ! density-site force on site
      DIMENSION FODSS( maxatm,NDIM )              ! density-site force on site
      DIMENSION FNSS(  maxatm,NDIM )              ! nuclear-site and LJ force on site
      DIMENSION FRCS(  maxatm,NDIM )              ! QM/MM force on site

      DIMENSION FRCN( NNUC,NDIM )                 ! QM/MM force on nuclear

      DIMENSION ADX( NMAXX/(NX*2) )
      DIMENSION ADY( NMAXY/(NY*2) )
      DIMENSION ADZ( NMAXZ/(NZ*2) )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      IF(MYID.EQ.0) THEN
        WRITE(*,*) 'Start qmf1 using reduced density'
      ENDIF

      DV = 8.D0*DX*DY*DZ

      JX2 = JX/2
      JY2 = JY/2
      JZ2 = JZ/2

C
C     number of loop in QM and MM subsystems
C
      nlmm = n-NNUC+NLINK

      DO L  = 1,3
      DO in = 1,n

        FDSS(  in,L ) = 0.D0
        FODSS( in,L ) = 0.D0
        FNSS(  in,L ) = 0.D0
        FRCS(  in,L ) = 0.D0
            
      ENDDO
      ENDDO
     
      DO L=1,3
      DO in=1,NNUC
     
        FRCN( in,L )   = 0.D0
            
      ENDDO
      ENDDO

C-----density-site force CALCULATION
     
      TIME1 = MPI_WTIME()

      ONE    = 1.D0
      TWO    = 2.D0
      FOUR   = 4.D0
      RRCUT2 = ONE/(RRCUT*RRCUT)
      SQPII  = ONE/DSQRT(PI)
      A2SQPI = TWO*ALPHA*SQPII

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO K=1,JZ2
        ADZ(K) = DZ*( 2*K-NNZ ) - DZ/2.D0
      ENDDO
      DO J=1,JY2
        ADY(J) = DY*( 2*J-NNY ) - DY/2.D0
      ENDDO
      DO I=1,JX2
        ADX(I) = DX*( 2*I-NNX ) - DX/2.D0
      ENDDO

      DO 50 iin=1,nion1
        in = iiont(iin)

      IF(NMMSW(in).EQ.1) THEN
C       WRITE(*,*) 'enter 1'
        IF(MYID.EQ.0) WRITE(*,*) 'enter 1',in
        GOTO 50
      ENDIF

      DO 60 K=1,JZ2
       RZ = ADZ(K) - SCRD( 3,in )
      DO 70 J=1,JY2
       RY = ADY(J) - SCRD( 2,in )
      DO 80 I=1,JX2

       ARHO = RRHO(I,J,K)

       IF( ARHO.GT.DCUT ) THEN

       RX = ADX(I) - SCRD( 1,in )

       R2 = RX**2+RY**2+RZ**2
       R  = DSQRT( R2 )

       IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
         RI    = ONE/R
         R2I   = RI*RI
         R3I   = RI*R2I
         APR   = ALPHA*R
         C2RDV = ARHO * chgmm( in ) * DV
         TFDSS = C2RDV*(A2SQPI*DEXP(-APR**2)*R2I-ERF(APR)*R3I)
       ELSEIF(R.GT.DRCUT) THEN
         ARD   = R - DRCUT
         ARD2  = ARD*ARD*RRCUT2
         ESW1  = ONE - ARD2
         ESW2  = ESW1*ESW1
         ESW   = ONE - ESW2
         ESWA  = FOUR*ESW1*ARD*RRCUT2
         RI    = ONE/R
         R2I   = RI*RI
         APR   = ALPHA*R
         C2RDV = ARHO * chgmm( in ) * DV
         DER2I = ERF(APR)*R2I
         TFDSS = ESW*C2RDV*(A2SQPI*DEXP(-APR**2)*R2I-DER2I*RI)
     *         + ESWA*C2RDV*DER2I
       ELSE
C        WRITE(*,*) 'enter 2', MYID
         GOTO 50
       ENDIF

       FODSS(in,1) = FODSS(in,1)-RX*TFDSS
       FODSS(in,2) = FODSS(in,2)-RY*TFDSS
       FODSS(in,3) = FODSS(in,3)-RZ*TFDSS

       ENDIF

80    CONTINUE
70    CONTINUE
60    CONTINUE

50    CONTINUE

C     NUMDAT = n*NDIM
      NUMDAT = maxatm*NDIM
      CALL MPI_REDUCE( FODSS(1,1),FDSS(1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )
C     WRITE(*,*) 'MYID:qmf1=',MYID,RHO(5,5,5)

      TIME2 = MPI_WTIME()

      IF(MYID.EQ.0) THEN
        WRITE(*,99) 'Elapsed Time (qmf1) = ',TIME2-TIME1,' scnds'
      ENDIF

C-----nuclear-site force CALCULATION
C-----LJ force CALCULATION

      IF(MYID.EQ.0) THEN

      SIX   = 6.D0

      DO 90 in=1,nlmm
      DO 100 M=1,NNUC 

        RX=PNUC( M,1 ) - SCRD( 1,in )
        RY=PNUC( M,2 ) - SCRD( 2,in )
        RZ=PNUC( M,3 ) - SCRD( 3,in )
        R2=RX**2+RY**2+RZ**2
        R   = DSQRT(R2)
        RI  = ONE/R
        R2I = RI*RI

        ASIG2  = SQMMM(M,in)*R2I
        ASIG6  = ASIG2*ASIG2*ASIG2
        ASIG12 = ASIG6*ASIG6
        BSIG   = TWO*ASIG12 - ASIG6
        TFEPS  = SIX*EPQMMM(M,in)
        CREPS  = TFEPS*R2I*BSIG
        CREPSX = CREPS*RX
        CREPSY = CREPS*RY
        CREPSZ = CREPS*RZ

        FNSS(in,1) = FNSS(in,1) + CREPSX   ! LJ force on site
        FNSS(in,2) = FNSS(in,2) + CREPSY
        FNSS(in,3) = FNSS(in,3) + CREPSZ

        FRCN(M,1)  = FRCN(M,1) + CREPSX    ! LJ force on nuclear
        FRCN(M,2)  = FRCN(M,2) + CREPSY
        FRCN(M,3)  = FRCN(M,3) + CREPSZ

100   CONTINUE
90    CONTINUE
      
      TIME3 = MPI_WTIME()
        WRITE(*,99) 'Elapsed Time (qmf2) = ',TIME3-TIME2,' scnds'

      DO 95 iin=1,nion1
        in = iiont(iin)
      DO 105 M=1,NNUC 

        RX=PNUC( M,1 ) - SCRD( 1,in )
        RY=PNUC( M,2 ) - SCRD( 2,in )
        RZ=PNUC( M,3 ) - SCRD( 3,in )
        R2=RX**2+RY**2+RZ**2
        R   = DSQRT(R2)
        RI  = ONE/R
        R2I = RI*RI

       IF(NMMSW(in).EQ.1) GOTO 105

       IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
         R3I   = RI*R2I
         ACHR  = ZVAL( M )*chgmm( in )*R3I
       ELSEIF(R.GT.DRCUT) THEN
         ARD   = R - DRCUT
         ARD2  = ARD*ARD*RRCUT2
         ESW1  = ONE - ARD2
         ESW2  = ESW1*ESW1
         ESW   = ONE - ESW2
         ESWA  = FOUR*ESW1*ARD*RRCUT2
         ZVPCH = ZVAL( M )*chgmm( in )*R2I
         ACHR  = ESW*ZVPCH*RI - ESWA*ZVPCH
       ELSE
         GOTO 105
       ENDIF

        CHRX = ACHR*RX
        CHRY = ACHR*RY
        CHRZ = ACHR*RZ

        FNSS(in,1) = FNSS(in,1) + CHRX      ! Coulomb force on site
        FNSS(in,2) = FNSS(in,2) + CHRY
        FNSS(in,3) = FNSS(in,3) + CHRZ

        FRCN(M,1)  = FRCN(M,1) + CHRX       ! Coulomb force on nuclear
        FRCN(M,2)  = FRCN(M,2) + CHRY
        FRCN(M,3)  = FRCN(M,3) + CHRZ

105   CONTINUE
95    CONTINUE

      TIME4 = MPI_WTIME()
        WRITE(*,99) 'Elapsed Time (qmf3) = ',TIME4-TIME3,' scnds'

      DO L=1,3
      DO in=1,nlmm
        FRCS(in,L) = FDSS(in,L) - FNSS(in,L)                 !QM/MM force on site
      ENDDO
      ENDDO

      ENDIF

99    FORMAT( X,A22,F10.6,3X,A6 )

      RETURN
      END

C----------------------------------------------------------
C     Partial Charge Correction for pseudopotentials 
C     S.G. Louie, S.Froyen, and M.L. Cohen,
C     Phys. Rev. B 26, 1738, 1982

      SUBROUTINE PCC

C     for the use with RBLYP or UBLYP
C----------------------------------------------------------

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      
      real*4 tim,ta(2)

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol
      include 'nlocd.i'                         ! Vmol  non-local d
      include 'pc_crr.i'                        ! Vmol  pcc

      PARAMETER ( NP = 4 )
      PARAMETER ( PI = 3.14159265358979323D0 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          FREEZE*5,PRINT*5,DGF*3

      DIMENSION DPC(100,421)
      DIMENSION XA( NP )
      DIMENSION YS( NP )
      DIMENSION A( 421 )

      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ 
      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)


      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      DPC(:,:) = 0.D0

C---- Partial Charge Density for Core Electrons of Mn -----

      IF(MYID.EQ.0) THEN

      WRITE(*,*) ' MYID = 0 is now in pcc.f '

      OPEN(78,FILE=PS_DIR//'/MN/MN_PCC/MNPC.DAT'    ! 2016.01.28
     *              ,STATUS='OLD')                  ! 2016.01.28

      READ(78,*) (INDX,R_DUMMY  ,K=1,421),
     *             DUMMY,
     *           (INDX,DPC(25,K),K=1,421)

      ENDIF

C----- Put PCC data on other cpus -------------------

      DO 20 K = 1,100     ! loop over elements

      NATOM = K

      DO NA = 1,NNUC
        NUMZ = INT( ZA(NA) )
        IF(NUMZ.EQ.NATOM) GOTO 10
      ENDDO 

      GOTO 20

10    CONTINUE

      IF( NPC(NA) == 0 ) CYCLE

      IF(MYID.EQ.0) THEN

      DO I=1,421
        A(I)=DPC(NATOM,I)
      ENDDO 

      ENDIF

      CALL MPI_BCAST(A(1),421,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)

      IF(MYID.NE.0) THEN

      DO I=1,421
        DPC(NATOM,I)=A(I)
      ENDDO

      ENDIF

20    CONTINUE

C----------------------------------------------------

      FAC = 1.D0/(4.D0*PI)

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      PCDNSA(:,:,:) = 0.D0
      PCDNSB(:,:,:) = 0.D0

      DO NA = 1,NNUC

      IF( NPC(NA) == 0 ) CYCLE

      NUMZ = INT( ZA( NA ) )
      RLN  = DLOG(RAD(NUMZ,1))

      DO KK = 1,JZ
      DO JJ = 1,JY
      DO II = 1,JX

        RX = DX*( II-NNX ) - PNUC(NA,1)
        RY = DY*( JJ-NNY ) - PNUC(NA,2)
        RZ = DZ*( KK-NNZ ) - PNUC(NA,3)
        R2 = RX**2 + RY**2 + RZ**2
        R  = DSQRT(R2) 

        IF( R .LT. RAD( NUMZ,3 ) ) THEN
            
          DO IPOL = 1,NP
            XA(IPOL) = RAD(NUMZ,IPOL)
            YS(IPOL) = DPC(NUMZ,IPOL)
          ENDDO

          VAL = 0.D0

        ELSE
           
          NRAD = INT( ( 0.5D0*DLOG(R2) - RLN )  
     *         / DRLOG(NUMZ) ) + 1 

          IT = NRAD - 2
          DO IPOL = 1,NP
            XA(IPOL) = RAD(NUMZ,IT+IPOL)
            YS(IPOL) = DPC(NUMZ,IT+IPOL)
          ENDDO

          CALL POLINT(XA,YS,NP,R,Y1,DELY)
          VAL = FAC*Y1/(R**2)       

        ENDIF

        PCDNSA(  II,JJ,KK) = PCDNSA(II,JJ,KK) + 0.5D0*VAL

        IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN
          PCDNSB(II,JJ,KK) = PCDNSB(II,JJ,KK) + 0.5D0*VAL
        ENDIF

      ENDDO
      ENDDO
      ENDDO
        
C     IF(MYID .EQ. 17) THEN
C       WRITE(*,*) NA,DPC(25,30)
C     ENDIF

      ENDDO

C     REWIND(24)
      DO 65 K = 1,JZ  
      DO 55 J = 1,JY
      DO 45 I = 1,JX
C       WRITE(24) PCDNSA(I,J,K)
45    CONTINUE
55    CONTINUE
65    CONTINUE

      RETURN
      END


C---------------------------------------------------

      SUBROUTINE TRWF2( CRD,WFA,WFB,WFNA,WFNB )
       
      IMPLICIT REAL*8 ( A-H,O-Z )
      IMPLICIT INTEGER*4 ( I-N )
      
      include "mpif.h"                          ! mpi
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX =  80 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NORA =  49 )
C     PARAMETER ( NORB =   1 )
C     PARAMETER ( MORA =  49 )
C     PARAMETER ( MORB =   1 )
      PARAMETER ( PI   =   3.14159265358979323D0 )
C     PARAMETER ( ALPHA = 3.D-1 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5

      CHARACTER ATMORB*5
      CHARACTER ATMINX*2
      
      CHARACTER NAME1*11,NAME2*13

      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / AVRHO / RHOAV(NMAXX/NX,NMAXY/NY,NMAXZ/NZ)

      DIMENSION WFA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION WFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
C     DIMENSION WFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,1 )
      DIMENSION WFNA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION WFNB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION NRDA( NORA )
      DIMENSION NRDB( NORB )
      
      DIMENSION CRD(    NNUC,NDIM )
      DIMENSION COEMOA( NBASIS,NBASIS )
      DIMENSION COEMOB( NBASIS,NBASIS )
      DIMENSION BASIS(  NBASIS )
      DIMENSION BASINF( 5,20,3 )
      DIMENSION NPRIM( 20)
      DIMENSION NTYP(  20)

      NAMELIST/INIDAT/COE,NDEN,DTMD,MDMAX,TEMP,
     *                NRST,NOPT,EXC,NQMMM,NRVLC,
     *                NCHK,CONV,FREEZE,PRINT,
C    *                DGF,NMRDF                  
     *                DGF,NMRDF,NMM2,NLINK          ! Vmol

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      NDAT = NNUC*NDIM

      CALL MPI_BCAST( CRD(1,1), NDAT,MPI_DOUBLE_PRECISION,0,   
     *                               MPI_COMM_WORLD,IERR )

C     if(myid.eq.5) then
C       write(*,*) crd(10,1) 
C     endif

C----- open avrho.dat file -----

      write (NAME2, '("avrho.dat", i4.4)') MYID                 ! takahasi 2013.07.30
      OPEN(24,FILE=NAME2,STATUS='UNKNOWN',FORM='UNFORMATTED')   ! takahasi 2013.07.30

C---- COMPUTE PARAMETERS ----

      DV = DX*DY*DZ

      IF(MYID .EQ. 0) THEN

      OPEN(75,FILE ='basis.dat',STATUS='OLD')       ! information for the basis set given by gfinput option
C     OPEN(86,FILE ='fort.7',   STATUS='OLD')       ! lcao coefficients generated by punch=mo option
      OPEN(86,FILE ='valence.dat',STATUS='OLD')     ! lcao coefficients generated by punch=mo option

      IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                 .OR. EXC.EQ.'RXalpha') THEN

         READ(86,*)

         DO J=1,NORA
            READ(86,'()') 
            READ(86,'(5D15.8)') (COEMOA(I,J),I=1,NBASIS)
         ENDDO

      ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                     .OR. EXC.EQ.'Uxalpha') THEN

         READ(86,*)

         DO J=1,NORA
            READ(86,'()') 
            READ(86,'(5D15.8)') (COEMOA(I,J),I=1,NBASIS)
          ENDDO

         READ(86,*)

         DO J=1,NORB
            READ(86,'()') 
            READ(86,'(5D15.8)') (COEMOB(I,J),I=1,NBASIS)
          ENDDO

      ENDIF

      ENDIF


      IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                 .OR. EXC.EQ.'RXalpha') THEN
        NUM = NORA*NBASIS
        CALL MPI_BCAST(COEMOA(1,1),NUM,MPI_DOUBLE_PRECISION,0,
     *                 MPI_COMM_WORLD,IERR)
      ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                      .OR. EXC.EQ.'Uxalpha') THEN
        NUM = NORA*NBASIS
        CALL MPI_BCAST(COEMOA(1,1),NUM,MPI_DOUBLE_PRECISION,0,
     *                 MPI_COMM_WORLD,IERR)
        NUM = NORB*NBASIS
        CALL MPI_BCAST(COEMOB(1,1),NUM,MPI_DOUBLE_PRECISION,0,
     *                 MPI_COMM_WORLD,IERR)
      ENDIF

C     if(myid .eq. 3) then
C     write(*,*) COEMOA(1,1)
C     endif

      WFA(:,:,:,:) = 0.0D0
      WFB(:,:,:,:) = 0.0D0

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      INDX = 0
      DO N=1,NNUC

       BASINF(:,:,:) = 0.D0
       NSHELL = 0

      IF(MYID .EQ. 0) THEN

       READ(75,*) NATOM,IDUM
       IF( NATOM .NE. N ) THEN
         WRITE(*,*) 'Warining: natom is not coincident with n.'
         STOP
       ENDIF

       DO     ! loop over primitive basis in each atom

           READ(75,*) ATMORB
           IF(ATMORB .EQ. '****') THEN
             EXIT  
           ELSEIF(ATMORB .EQ. 'G') THEN
             STOP
           ENDIF

           BACKSPACE(75)
           READ(75,*) ATMORB,NGAUSS
           NSHELL = NSHELL + 1 
           NPRIM(NSHELL) = NGAUSS

           IF(ATMORB .EQ. 'S') THEN
             NTYP(NSHELL) = 1
             DO N1=1,NGAUSS
               READ(75,*) FACTOR,COEF
               BASINF(1,N1,1) = FACTOR
               BASINF(1,N1,2) = COEF
             ENDDO
           ELSEIF(ATMORB .EQ. 'SP') THEN
             NTYP(NSHELL) = 2
             DO N1=1,NGAUSS
               READ(75,*) FACTOR,COEFS,COEFP
               BASINF(2,N1,1) = FACTOR
               BASINF(2,N1,2) = COEFS
               BASINF(2,N1,3) = COEFP
             ENDDO
           ELSEIF(ATMORB .EQ. 'P') THEN
             NTYP(NSHELL) = 3
             DO N1=1,NGAUSS
               READ(75,*) FACTOR,COEF
               BASINF(3,N1,1) = FACTOR
               BASINF(3,N1,2) = COEF
             ENDDO
           ELSEIF(ATMORB .EQ. 'D') THEN
             NTYP(NSHELL) = 4
             DO N1=1,NGAUSS
               READ(75,*) FACTOR,COEF
               BASINF(4,N1,1) = FACTOR
               BASINF(4,N1,2) = COEF
             ENDDO
           ELSEIF(ATMORB .EQ. 'F') THEN
             NTYP(NSHELL) = 5
             DO N1=1,NGAUSS
               READ(75,*) FACTOR,COEF
               BASINF(5,N1,1) = FACTOR
               BASINF(5,N1,2) = COEF
             ENDDO
           ENDIF

       ENDDO

      ENDIF

      CALL MPI_BCAST(BASINF(1,1,1),300,MPI_DOUBLE_PRECISION,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NTYP(1),20,MPI_INTEGER,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NPRIM(1),20,MPI_INTEGER,0,
     *               MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NSHELL,1,MPI_INTEGER,0,
     *               MPI_COMM_WORLD,IERR)

      if(myid .eq. 3) then
      write(*,*) 'trwf2'
      write(*,*) BASINF(1,1,1)
      write(*,*) NTYP(1),NPRIM(1),NSHELL
      endif

      DO N2 = 1,NSHELL

      INDX  = INDX+1
      NNTYP = NTYP(N2)

C     write(*,*) NTYP(1),NPRIM(1),NSHELL

      DO K=1,JZ
         RZA = DZ*( K-NNZ ) - CRD(N,3)
      DO J=1,JY
         RYA = DY*( J-NNY ) - CRD(N,2)
      DO I=1,JX
         RXA = DX*( I-NNX ) - CRD(N,1)

         SQR = RXA**2 + RYA**2 + RZA**2

!------- for s-type atomic orb. -----------

          IF(NNTYP .EQ. 1) THEN

            PNORM=(8.0D0/PI**3)**0.25D0
              
            BASIS(:) = 0
            DO N1 = 1,NPRIM(N2)
            
              FACTOR =  BASINF(1,N1,1) 
              COEF   =  BASINF(1,N1,2) 

              PRES=FACTOR**0.75D0*DEXP(-FACTOR*SQR)
              PRES=COEF*PNORM*PRES
              BASIS(INDX)=BASIS(INDX)+PRES

            ENDDO

            IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                  .OR. EXC.EQ.'RXalpha') THEN
             DO NO = 1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX)
     *                                    *COEMOA(INDX,NVAL) 
             ENDDO
            ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                      .OR. EXC.EQ.'Uxalpha') THEN
             DO NO=1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX)*COEMOA(INDX,NVAL) 
             ENDDO
             DO NO=1,NORB 
               NVAL =         NO
               WFB(I,J,K,NO)=WFB(I,J,K,NO)+BASIS(INDX)*COEMOB(INDX,NVAL) 
             ENDDO
            ENDIF
     
          ELSEIF(NNTYP .EQ. 2) THEN

            PNORMS=(8.0D0/PI**3)**0.25D0
            PNORMP=(128.0D0/PI**3)**0.25D0

            BASIS(:) = 0
            DO N1 = 1,NPRIM(N2)

                FACTOR =  BASINF(2,N1,1) 
                COEFS  =  BASINF(2,N1,2) 
                COEFP  =  BASINF(2,N1,3) 

                PRESP=FACTOR**0.75D0*DEXP(-FACTOR*SQR)
                PRESP=COEFS*PNORMS*PRESP
                BASIS(INDX)=BASIS(INDX)+PRESP

                PRESPX=FACTOR**1.25D0*(RXA)*DEXP(-FACTOR*SQR)
                PRESPX=COEFP*PNORMP*PRESPX
                BASIS(INDX+1)=BASIS(INDX+1)+PRESPX

                PRESPY=FACTOR**1.25D0*(RYA)*DEXP(-FACTOR*SQR)
                PRESPY=COEFP*PNORMP*PRESPY
                BASIS(INDX+2)=BASIS(INDX+2)+PRESPY

                PRESPZ=FACTOR**1.25D0*(RZA)*DEXP(-FACTOR*SQR)
                PRESPZ=COEFP*PNORMP*PRESPZ
                BASIS(INDX+3)=BASIS(INDX+3)+PRESPZ

            ENDDO
         
            IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                  .OR. EXC.EQ.'RXalpha') THEN
             DO NO = 1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,3
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                      .OR. EXC.EQ.'Uxalpha') THEN
             DO NO=1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,3
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
             DO NO=1,NORB 
               NVAL =         NO
               DO INC = 0,3
               WFB(I,J,K,NO)=WFB(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOB(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ENDIF

          ELSEIF(NNTYP .EQ. 3) THEN

            PNORMP=(128.0D0/PI**3)**0.25D0

            BASIS(:) = 0
            DO N1 = 1,NPRIM(N2)

              FACTOR =  BASINF(3,N1,1) 
              COEFP  =  BASINF(3,N1,2) 

                PREPX=FACTOR**1.25D0*(RXA)*DEXP(-FACTOR*SQR)
                PREPX=COEFP*PNORMP*PREPX
                BASIS(INDX)=BASIS(INDX)+PREPX

                PREPY=FACTOR**1.25D0*(RYA)*DEXP(-FACTOR*SQR)
                PREPY=COEFP*PNORMP*PREPY
                BASIS(INDX+1)=BASIS(INDX+1)+PREPY

                PREPZ=FACTOR**1.25D0*(RZA)*DEXP(-FACTOR*SQR)
                PREPZ=COEFP*PNORMP*PREPZ
                BASIS(INDX+2)=BASIS(INDX+2)+PREPZ

            ENDDO
         
            IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                  .OR. EXC.EQ.'RXalpha') THEN
             DO NO = 1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,2
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                      .OR. EXC.EQ.'Uxalpha') THEN
             DO NO=1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,2
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
             DO NO=1,NORB 
               NVAL =         NO
               DO INC = 0,2
               WFB(I,J,K,NO)=WFB(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOB(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ENDIF

          ELSEIF(NNTYP .EQ. 4) THEN

              PNORMD1=(2048.0D0/(9.0D0*PI**3))**0.25D0
              PNORMD2=(2048.0D0/(PI**3))**0.25D0
              PNORMD3=(2048.0D0/(16.0D0*PI**3))**0.25D0
              PNORMD4=(2048.0D0/(144.0D0*PI**3))**0.25D0

            BASIS(:) = 0
            DO N1 = 1,NPRIM(N2)

              FACTOR =  BASINF(4,N1,1) 
              COEFD  =  BASINF(4,N1,2) 

                PRED0=(2.0D0*RZA**2-RXA**2-RYA**2)*DEXP(-FACTOR*SQR)
                PRED0=FACTOR**1.75D0*COEFD*PNORMD4*PRED0
                BASIS(INDX)=BASIS(INDX)+PRED0

                PRED1=FACTOR**1.75D0*(RZA)*(RXA)*DEXP(-FACTOR*SQR)
                PRED1=COEFD*PNORMD2*PRED1
                BASIS(INDX+1)=BASIS(INDX+1)+PRED1

                PRED2=FACTOR**1.75D0*(RYA)*(RZA)*DEXP(-FACTOR*SQR)
                PRED2=COEFD*PNORMD2*PRED2
                BASIS(INDX+2)=BASIS(INDX+2)+PRED2

                PRED3=FACTOR**1.75D0*(RXA**2-RYA**2)*DEXP(-FACTOR*SQR)
                PRED3=COEFD*PNORMD3*PRED3
                BASIS(INDX+3)=BASIS(INDX+3)+PRED3

                PRED4=FACTOR**1.75D0*(RXA)*(RYA)*DEXP(-FACTOR*SQR)
                PRED4=COEFD*PNORMD2*PRED4
                BASIS(INDX+4)=BASIS(INDX+4)+PRED4

            ENDDO

            IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                  .OR. EXC.EQ.'RXalpha') THEN
             DO NO = 1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,4
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                      .OR. EXC.EQ.'Uxalpha') THEN
             DO NO=1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,4
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
             DO NO=1,NORB 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,4
               WFB(I,J,K,NO)=WFB(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOB(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ENDIF

          ELSEIF(NNTYP .EQ. 5) THEN
     
               PNORMF1=(32768.0D0/(3600.0D0*PI**3))**0.25D0
               PNORMF2=(32768.0D0/(1600.0D0*PI**3))**0.25D0
               PNORMF3=(32768.0D0/(16.0D0*PI**3))**0.25D0
               PNORMF4=(32768.0D0/PI**3)**0.25D0
               PNORMF5=(32768.0D0/(576.0D0*PI**3))**0.25D0

            BASIS(:) = 0
            DO N1 = 1,NPRIM(N2)

              FACTOR =  BASINF(5,N1,1) 
              COEFF  =  BASINF(5,N1,2) 

                 PREF0=RZA*(5.0D0*RZA**2-3.0D0*SQR)*DEXP(-FACTOR*SQR)
                 PREF0=FACTOR**2.25D0*COEFF*PNORMF1*PREF0
                 BASIS(INDX)=BASIS(INDX)+PREF0

                 PREF1=RXA*(5.0D0*RZA**2-SQR)*DEXP(-FACTOR*SQR)
                 PREF1=FACTOR**2.25D0*COEFF*PNORMF2*PREF1
                 BASIS(INDX+1)=BASIS(INDX+1)+PREF1

                 PREF2=RYA*(5.0D0*RZA**2-SQR)*DEXP(-FACTOR*SQR)
                 PREF2=FACTOR**2.25D0*COEFF*PNORMF2*PREF2
                 BASIS(INDX+2)=BASIS(INDX+2)+PREF2

                 PREF3=RZA*(RXA**2-RYA**2)*DEXP(-FACTOR*SQR)
                 PREF3=FACTOR**2.25D0*COEFF*PNORMF3*PREF3
                 BASIS(INDX+3)=BASIS(INDX+3)+PREF3

                 PREF4=RXA*RYA*RZA*DEXP(-FACTOR*SQR)
                 PREF4=FACTOR**2.25D0*COEFF*PNORMF4*PREF4
                 BASIS(INDX+4)=BASIS(INDX+4)+PREF4

                 PREF5=RXA*(RXA**2-3.0D0*RYA**2)*DEXP(-FACTOR*SQR)
                 PREF5=FACTOR**2.25D0*COEFF*PNORMF5*PREF5
                 BASIS(INDX+5)=BASIS(INDX+5)+PREF5

                 PREF6=RYA*(3.0D0*RXA**2-RYA**2)*DEXP(-FACTOR*SQR)
                 PREF6=FACTOR**2.25D0*COEFF*PNORMF5*PREF6
                 BASIS(INDX+6)=BASIS(INDX+6)+PREF6

            ENDDO

            IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                  .OR. EXC.EQ.'RXalpha') THEN
             DO NO = 1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,6
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                      .OR. EXC.EQ.'Uxalpha') THEN
             DO NO=1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,6
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
             DO NO=1,NORB 
               NVAL =         NO
               DO INC = 0,6
               WFB(I,J,K,NO)=WFB(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOB(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ENDIF

          ENDIF

       ENDDO
       ENDDO
       ENDDO

          IF    (NNTYP .EQ. 1) THEN
            NADD = 0
          ELSEIF(NNTYP .EQ. 2) THEN
            NADD = 3
          ELSEIF(NNTYP .EQ. 3) THEN
            NADD = 2
          ELSEIF(NNTYP .EQ. 4) THEN
            NADD = 4
          ELSEIF(NNTYP .EQ. 5) THEN
            NADD = 6
          ENDIF

          INDX = INDX + NADD

       ENDDO

       ENDDO

C     if(myid.eq.0) then
C       write(*,*) wfa(10,10,10,10)
C       write(*,*) basis(10)
C     endif

C----- Normalize for alpha spin -----

      DO 40 K = 1, NORA      

        ANORM = 0.D0

        DO 50 L = 1,JZ      
        DO 50 M = 1,JY      
        DO 50 N = 1,JX      
          ANORM = ANORM + WFA( N,M,L,K )**2*DV
50      CONTINUE

        CALL MPI_REDUCE(ANORM,BNORM,1,MPI_DOUBLE_PRECISION,
     *                  MPI_SUM,0,MPI_COMM_WORLD,IERR)

      if(myid.eq.0) then
      write(*,*) 'anorm=',BNORM
      endif
        ANORM = DSQRT(BNORM)

        CALL MPI_BCAST(ANORM,1,MPI_DOUBLE_PRECISION,
     *                 0,MPI_COMM_WORLD,IERR)

        DO 60 L = 1,JZ      
        DO 60 M = 1,JY      
        DO 60 N = 1,JX      
          WFA(  N,M,L,K ) = WFA( N,M,L,K ) / ANORM 
          WFNA( N,M,L,K ) = WFA( N,M,L,K ) 
60      CONTINUE

40    CONTINUE

      DO K = 1, NORA
        NRDA(K) = K
      ENDDO

C     CALL  DIAG(      NRDA,NORA,WFNA,WFA )
      CALL MDIAG (     NRDA,NORA,WFNA,WFA )
C     CALL MDIAG_BLK ( NRDA,NORA,WFNA,WFA )

C----- normalize for beta spin -----

      IF(EXC.EQ.'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                .OR. EXC.EQ.'UXalpha') THEN

        DO 70 K = 1, NORB      

          ANORM = 0.D0

          DO 80 L = 1,JZ      
          DO 80 M = 1,JY      
          DO 80 N = 1,JX      
            ANORM = ANORM + WFB( N,M,L,K )**2*DV
80        CONTINUE

        CALL MPI_REDUCE(ANORM,BNORM,1,MPI_DOUBLE_PRECISION,
     *                  MPI_SUM,0,MPI_COMM_WORLD,IERR)

        ANORM = DSQRT(BNORM)

C       IF(MYID.EQ.0) THEN
C         write(*,*) 'anorm =',ANORM
C       ENDIF

        CALL MPI_BCAST(ANORM,1,MPI_DOUBLE_PRECISION,
     *                 0,MPI_COMM_WORLD,IERR)

          DO 90 L = 1,JZ           ! 2005.07.25 takahashi
          DO 90 M = 1,JY      
          DO 90 N = 1,JX      
            WFB(  N,M,L,K ) = WFB( N,M,L,K ) / ANORM
            WFNB( N,M,L,K ) = WFB( N,M,L,K ) 
90        CONTINUE

70      CONTINUE

        DO K = 1, NORB
          NRDB(K) = K
        ENDDO

C       CALL  DIAG(      NRDB,NORB,WFNB,WFB )
        CALL MDIAG (     NRDB,NORB,WFNB,WFB )
C       CALL MDIAG_BLK ( NRDB,NORB,WFNB,WFB )

      ENDIF

C      IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
C    *                  .OR. EXC.EQ.'RXalpha') THEN

C        DO I=1,NMAX
C        DO J=1,NMAX
C            WRITE(3,'(6E13.5)') ((WFA(I,J,K,NO),NO=2,5),K=1,NMAX)
C        ENDDO
C        ENDDO

C       ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
C    *                       .OR. EXC.EQ.'Uxalpha') THEN

C        DO I=1,NMAX
C        DO J=1,NMAX
C            WRITE(3,'(6E13.5)') ((WFA(I,J,K,NO),NO=2,5),
C    *                           (WFB(I,J,K,NO),NO=2,5),K=1,NMAX)
C        ENDDO
C        ENDDO
C      
C       ENDIF

       RETURN
       END  


      SUBROUTINE DG4L_Y(VNUCN,VNUC)          ! for non-periodic system
                                             ! for link atoms
      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'

      include 'QMpara.i'                     ! Vmol
!     include 'sizes.i'

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NNUC = 36 )
C     PARAMETER ( NDIM =  3 )
      PARAMETER ( RCUT =  3.5D0 )
      PARAMETER ( NSL  = 9000 )

      PARAMETER ( PI    = 3.14159265358979323D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX, DY, DZ
      COMMON / ACELL / XL, YL, ZL
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)
      COMMON / DBLG  / WNLOCS( NNUC,NSL ),WNLOCPX( NNUC,NSL ),
     *                 WNLOCPY(NNUC,NSL ),WNLOCPZ( NNUC,NSL )
C     COMMON / DBLG1 / VNUCQM( NMAX,NMAX,NMAX )
      COMMON / PLOC / VLOC( NNUC,NSL )
      COMMON / BHS / C1(100),C2(100),AL1(100),AL2(100)

      DIMENSION VNUC(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VNUCL( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VNUCN( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      ! yamacta
      dimension WIJ_tmp(-2*NDEN:2*NDEN,-2*NDEN:2*NDEN,-2*NDEN:2*NDEN)
      dimension VLOCD2_tmp(-2*NDEN:2*NDEN,-2*NDEN:2*NDEN,-2*NDEN:2*NDEN)
      dimension VNLOCP_tmp(-2*NDEN:2*NDEN,-2*NDEN:2*NDEN,-2*NDEN:2*NDEN)
      dimension VNLOCS_tmp(-2*NDEN:2*NDEN,-2*NDEN:2*NDEN,-2*NDEN:2*NDEN)
      ! yamacta

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) 'Double Grid(4th-order Lagrange Interpolation): Start'
      WRITE(*,*) '--- for Link Atoms ---'
      STIME = MPI_WTIME()
      ENDIF

      RT3  = 1.D0/3.D0
      ns = NNUC - NLINK + 1
      
C----- initialize -----

      DO J=1, NSL 
      DO I=ns, NNUC  

        WNLOCS( I,J) = 0.D0
        WNLOCPX(I,J) = 0.D0
        WNLOCPY(I,J) = 0.D0
        WNLOCPZ(I,J) = 0.D0

        VLOC(I,J) = 0.D0

      ENDDO 
      ENDDO 

      DO K=1, JZ  
      DO J=1, JY 
      DO I=1, JX 

        VNUCL( I,J,K ) = 0.0D0     ! VNUC for Link atoms 2016.03.18 takahashi
        VNUC( I,J,K )  = 0.0D0

      ENDDO 
      ENDDO 
      ENDDO 

      NCOUNT = 0
      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

C----- compute weight factor for each coarse grid -----  

      DO 100 NA = ns, NNUC                  ! loop over atoms
        NUMZ = INT(  ZA( NA ) )
        RLN  = DLOG( RAD( NUMZ,1 ) )
        NCOUNT = 0
      DO 10 K=1, JZ
        DO 20 J=1, JY
          DO 30 I=1, JX

          RX = DX*( I-NNX ) - PNUC(NA,1) 
          RY = DY*( J-NNY ) - PNUC(NA,2)
          RZ = DZ*( K-NNZ ) - PNUC(NA,3)
          R2 = RX**2 + RY**2 + RZ**2
          R  = DSQRT(R2) 

          IF( R .LT. RCUT ) THEN
          NCOUNT = NCOUNT + 1
           IF( NCOUNT .GT. NSL ) THEN
           WRITE(*,*) ' dg4.f : size of dimension too small '
           STOP
          ENDIF

          DO 40 KD = -2*NDEN,2*NDEN
          DO 50 JD = -2*NDEN,2*NDEN
          DO 60 ID = -2*NDEN,2*NDEN

          DID = DABS(DBLE(ID)/DBLE(NDEN))         
          DJD = DABS(DBLE(JD)/DBLE(NDEN))         
          DKD = DABS(DBLE(KD)/DBLE(NDEN))         

          RXX = ( RX + DX*DBLE(ID)/DBLE(NDEN)) 
          RYY = ( RY + DY*DBLE(JD)/DBLE(NDEN)) 
          RZZ = ( RZ + DZ*DBLE(KD)/DBLE(NDEN)) 
          RD2 = RXX**2 + RYY**2 + RZZ**2
          RD  = DSQRT( RD2 )

          IF( RD .LT. RAD( NUMZ,1 ) ) THEN
          WGT1   =  RAD( NUMZ,1 ) -  RD 
          WGT2   =  RD
          VNLOCS_tmp(ID,JD,KD) = VNLS( NUMZ,1 )  
          VNLOCP_tmp(ID,JD,KD) = WGT2*VNLP( NUMZ,1 )  / ( WGT1 + WGT2 )     
          VLOCD2_tmp(ID,JD,KD) = VLDC( NUMZ,1 )  
C         NR = NR + 1
          ELSE
          NRAD   =  INT( ( 0.5D0*DLOG( RD2 ) - RLN )  
     *           / DRLOG(NUMZ) ) + 1 
          WGT1   =  RAD( NUMZ,NRAD+1 ) - RD
          WGT2   = -RAD( NUMZ,NRAD )   + RD
          VNLOCS_tmp(ID,JD,KD) = ( WGT1*VNLS( NUMZ,NRAD ) 
     *                 +   WGT2*VNLS( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )     
          VNLOCP_tmp(ID,JD,KD) = ( WGT1*VNLP( NUMZ,NRAD ) 
     *                 +   WGT2*VNLP( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )     
          VLOCD2_tmp(ID,JD,KD) = ( WGT1*VLDC( NUMZ,NRAD ) 
     *           +   WGT2*VLDC( NUMZ,NRAD+1 ) ) 
     *           / ( WGT1 + WGT2 )                
          ENDIF

          IF( (IABS(ID) .LE. NDEN)   .AND.                        ! 1
     *        (IABS(JD) .LE. NDEN)   .AND.
     *        (IABS(KD) .LE. NDEN) ) THEN 
          WIJ_tmp(ID,JD,KD)  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 2
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ_tmp(ID,JD,KD)  = ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 3
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ_tmp(ID,JD,KD)  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 4
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ_tmp(ID,JD,KD)  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .LE. NDEN)   .AND.                    ! 5
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ_tmp(ID,JD,KD)  =  ( 1.0-DID )*( 1.0+DID )*( 1.0-0.5*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 6
     *            (IABS(JD) .LE. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ_tmp(ID,JD,KD)  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0+DJD )*( 1.0-0.5*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 7
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .LE. NDEN) ) THEN 
          WIJ_tmp(ID,JD,KD)  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0+DKD )*( 1.0-0.5*DkD )    
     *           / DBLE(NDEN**3)

          ELSEIF( (IABS(ID) .GT. NDEN)   .AND.                    ! 8
     *            (IABS(JD) .GT. NDEN)   .AND.
     *            (IABS(KD) .GT. NDEN) ) THEN 
          WIJ_tmp(ID,JD,KD)  =  ( 1.0-DID )*( 1.0-0.5*DID )*( 1.0-RT3*DID )    
     *         *  ( 1.0-DJD )*( 1.0-0.5*DJD )*( 1.0-RT3*DJD )    
     *         *  ( 1.0-DKD )*( 1.0-0.5*DKD )*( 1.0-RT3*DKD )    
     *           / DBLE(NDEN**3)

          ENDIF

C--------------------------
C      non-local part      
C--------------------------

C----- for s-component -----

          WNLOCS( NA,NCOUNT )  = WNLOCS( NA,NCOUNT ) 
     *                         + WIJ_tmp(ID,JD,KD)*VNLOCS_tmp(ID,JD,KD)

C----- for p-component -----

          IF( RD .LT. 1.0D-06 ) THEN

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 

          NW = NW + 1

          ELSE

          WNLOCPX( NA,NCOUNT ) = WNLOCPX( NA,NCOUNT ) 
     *                         + WIJ_tmp(ID,JD,KD)*VNLOCP_tmp(ID,JD,KD)*RXX/RD 
          WNLOCPY( NA,NCOUNT ) = WNLOCPY( NA,NCOUNT ) 
     *                         + WIJ_tmp(ID,JD,KD)*VNLOCP_tmp(ID,JD,KD)*RYY/RD
          WNLOCPZ( NA,NCOUNT ) = WNLOCPZ( NA,NCOUNT ) 
     *                         + WIJ_tmp(ID,JD,KD)*VNLOCP_tmp(ID,JD,KD)*RZZ/RD

          ENDIF

C----------------------------------------------
C      local part ( BHS + local-d )
C----------------------------------------------

          IF( RD .LT. 1.0D-06 ) THEN

            VT = -ZVAL(NA)*( C1(NUMZ)*2.D0*DSQRT(AL1(NUMZ))            ! analytical function ( limit zero ) of BHS
     *         +             C2(NUMZ)*2.D0*DSQRT(AL2(NUMZ)) )  
     *         /  DSQRT(PI)           
            
          ELSE

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*RD2) )        ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*RD2) ) ) 
     *         /  RD           
            
          ENDIF

            VLOC( NA,NCOUNT ) = VLOC( NA,NCOUNT )
     *                        + WIJ_tmp(ID,JD,KD)*(VLOCD2_tmp(ID,JD,KD)+VT)    

60        CONTINUE
50        CONTINUE
40        CONTINUE

            VNUCL( I,J,K ) = VNUCL( I,J,K ) 
     *                     + VLOC( NA,NCOUNT )

          ELSE

C----------------------------------------------
C      local part ( BHS )
C----------------------------------------------

            VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*R2) )         ! analytical function of BHS
     *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*R2) ) ) 
     *         /  R           
            
            VNUCL( I,J,K ) = VNUCL( I,J,K ) + VT                       

          ENDIF

30        CONTINUE
20      CONTINUE
10    CONTINUE

100   CONTINUE

      CALL MPI_REDUCE(NCOUNT,NCOUNT1,1,MPI_INTEGER,MPI_SUM,0,
     *                MPI_COMM_WORLD,IERR)

      IF(MYID.EQ.0) THEN
      WRITE(*,*) '  ncount = ', NCOUNT1
      ETIME = MPI_WTIME()
      write(*,*) 'Elapsed Time (dg4l) = ',ETIME-STIME,' scnds'
C     write(*,*) 'VNUCN(5,5,5) =', VNUCN(5,5,5)
      ENDIF

      VNUC(:,:,:) = VNUCN(:,:,:) + VNUCL(:,:,:)      ! add vnucl due to Link atoms 2016.03.18 takahashi

C     write(*,*) 'NR = ',NR 
C     write(*,*) 'NW = ',NW 

      ne = NNUC - NLINK
      DO 105 NA = 1, ne                  ! loop over atoms
        NUMZ = INT( ZA( NA ) )
        NCOUNT = 0
      DO 15 K=1, JZ
        DO 25 J=1, JY
          DO 35 I=1, JX

C         RX = DX*( I-NNX ) - PNUC(NA,1) 
C         RY = DY*( J-NNY ) - PNUC(NA,2)
C         RZ = DZ*( K-NNZ ) - PNUC(NA,3)
C         R2 = RX**2 + RY**2 + RZ**2
C         R  = DSQRT(R2) 

          IF( R .LT. RCUT ) THEN
C         NCOUNT = NCOUNT + 1
C         VNUC( I,J,K ) = VNUC( I,J,K ) 
C    *                  + VLOC( NA,NCOUNT )

          ELSE

C----------------------------------------------
C      local part ( BHS )
C----------------------------------------------

C           VT = -ZVAL(NA)*( C1(NUMZ)*ERF( DSQRT(AL1(NUMZ)*R2) )                  ! analytical function of BHS
C    *         +             C2(NUMZ)*ERF( DSQRT(AL2(NUMZ)*R2) ) ) 
C    *         /  R           
            
C           VNUC( I,J,K ) = VNUC( I,J,K ) + VT                       

          ENDIF

35        CONTINUE
25      CONTINUE
15    CONTINUE

105   CONTINUE

      RETURN
      END


C-------------------------------------------
C     SUBROUTINE POINT CHARGE 
C-------------------------------------------

      SUBROUTINE POCH1_S( SCRD )      ! tuned by S.Sakuraba Mar.2015

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
!     include 'atoms.i'
C     include 'charge.i'

C     PARAMETER ( NMAX =  80 )
C     PARAMETER ( NDIM =   3 )
C     PARAMETER ( NNUC =  36 )
      PARAMETER ( ALPHA = 1.0D0    )
C     PARAMETER ( NLINK1 = 10 )
      PARAMETER ( RRCUT = 4.724D0 )
      PARAMETER ( RRCUT1= 5.669D0 )
      PARAMETER ( DRCUT = 0.945D0 )


      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PRMT4 / nion1,iiont(maxatm)
      COMMON / NCLR / PNUC( NNUC,NDIM )
C     COMMON / SLT / VPCE(NMAX,NMAX,NMAX),VPCZ,PLJ
      COMMON / SLT / VPCE(NMAXX/NX,NMAXY/NY,NMAXZ/NZ),VPCZ,PLJ
      COMMON / LJQMMM / SQMMM(NNUC,maxatm),EPQMMM(NNUC,maxatm),
     *                  chgmm(maxatm),chgqm(NNUC)

      REAL*8 CHARGE0, CHARGE2, CCRD0, CCRD2, RINV, PVPCE
      REAL*8 DRCUT_SQ, RRCUT1_SQ, R2T
      INTEGER NATM0, NATM2

      DIMENSION SCRD( NDIM,maxatm )
      DIMENSION ADX(  NMAXX/NX )
      DIMENSION ADY(  NMAXY/NY )
      DIMENSION ADZ(  NMAXZ/NZ )

      DIMENSION CHARGE0(maxatm), CHARGE2(maxatm)
      DIMENSION CCRD0( maxatm, NDIM ), CCRD2( maxatm, NDIM )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

C
C     number of loop in QM and MM subsystems
C

      IF(MYID.EQ.0) THEN
      OPEN(111,FILE='CHG.DAT',STATUS='UNKNOWN')            ! open output files 
        WRITE(*,*) 'Start poch1'
      qsum = 0.d0
      DO iin = 1, nion1
        in   = iiont(iin)
C       write(111,*) chgmm(in)
        qsum = qsum + chgmm(in)
      ENDDO
      write(*,*) 'qsum = ', qsum
      write(*,*) 'nion1 = ', nion1
      ENDIF

      TIME1 = MPI_WTIME()

      nlmm = n-NNUC+NLINK
      nlqm = NNUC-NMM2-NLINK

      ONE    = 1.D0
      RRCUT2 = ONE/(RRCUT*RRCUT)
      RRCUT1_SQ = RRCUT1 * RRCUT1
      DRCUT_SQ = DRCUT * DRCUT
C     NNMAX  = NMAX/2 + 1

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO K=1,JZ
        ADZ(K) = DZ*( K-NNZ )
      ENDDO
      DO J=1,JY
        ADY(J) = DY*( J-NNY )
      ENDDO
      DO I=1,JX
        ADX(I) = DX*( I-NNX )
      ENDDO

C     count atoms up for each type
      NATM0 = 0
      NATM2 = 0
      DO iin = 1, nion1
        in = iiont(iin)
        IF(NMMSW(in).EQ.0) THEN
          NATM0 = NATM0 + 1
          CCRD0(NATM0, 1) = SCRD(1, in)
          CCRD0(NATM0, 2) = SCRD(2, in)
          CCRD0(NATM0, 3) = SCRD(3, in)
          CHARGE0(NATM0) = chgmm(in)
        ELSEIF(NMMSW(in).EQ.2) THEN
          NATM2 = NATM2 + 1
          CCRD2(NATM2, 1) = SCRD(1, in)
          CCRD2(NATM2, 2) = SCRD(2, in)
          CCRD2(NATM2, 3) = SCRD(3, in)
          CHARGE2(NATM2) = chgmm(in)
        ENDIF
      ENDDO

!$omp parallel do
!$omp&default(none) 
!$omp&private(I, J, K, PVPCE, iin)
!$omp&private(RX, RY, RZ, R2, RINV)
!$omp&shared(JX, JY, JZ, NATM0, NATM2, CCRD0, CCRD2)
!$omp&shared(ADX, ADY, ADZ)
!$omp&shared(CHARGE0, CHARGE2, VPCE, RRCUT1_SQ)
      DO K=1,JZ
      DO J=1,JY
      DO I=1,JX
       PVPCE = 0.0D0
!OCL SIMD
       DO iin = 1, NATM0
        RZ = ADZ(K) - CCRD0(iin, 3)
        RY = ADY(J) - CCRD0(iin, 2)
        RX = ADX(I) - CCRD0(iin, 1)

        R2    = RX**2+RY**2+RZ**2
        RINV  = 1.0D0 / DSQRT( R2 )

        IF(R2.LT.RRCUT1_SQ) THEN
          RINV = 0.0D0
        ENDIF
        PVPCE = PVPCE - CHARGE0(iin) * RINV
       ENDDO
!OCL SIMD
       DO iin = 1, NATM2
        RZ = ADZ(K) - CCRD2(iin, 3)
        RY = ADY(J) - CCRD2(iin, 2)
        RX = ADX(I) - CCRD2(iin, 1)

        R2    = RX**2+RY**2+RZ**2
        RINV  = 1.0D0 / DSQRT( R2 )

        IF(R2.LT.RRCUT1_SQ) THEN
          RINV = 0.0D0
        ENDIF
        PVPCE = PVPCE - CHARGE2(iin) * RINV
       ENDDO
       VPCE(I, J, K) = PVPCE
      ENDDO
      ENDDO
      ENDDO

      DO iin = 1, NATM0
       DO K = 1,JZ
       DO J = 1,JY
        RZ = ADZ(K) - CCRD0(iin, 3)
        RY = ADY(J) - CCRD0(iin, 2)
        R2T = RY**2+RZ**2
        IF(R2T >= RRCUT1_SQ) THEN
          CYCLE
        ENDIF
        DO I = 1, JX
         RX = ADX(I) - CCRD0(iin, 1)
         R2 = R2T + RX**2
         RINV = 1.0D0 / DSQRT(R2)
         R = RINV * R2
         IF(R2 >= RRCUT1_SQ) THEN
           CYCLE
         ENDIF
         VPCE( I,J,K ) = VPCE( I,J,K )
     *                   - CHARGE0(iin)*(ERF(ALPHA*R)*RINV)
        ENDDO
       ENDDO
       ENDDO
      ENDDO

      DO iin = 1, NATM2
       DO K = 1,JZ
       DO J = 1,JY
        RZ = ADZ(K) - CCRD2(iin, 3)
        RY = ADY(J) - CCRD2(iin, 2)
        R2T = RY**2+RZ**2
        IF(R2T >= RRCUT1_SQ) THEN
          CYCLE
        ENDIF
        DO I = 1, JX
         RX = ADX(I) - CCRD2(iin, 1)
         R2 = R2T + RX**2
         RINV = 1.0D0 / DSQRT(R2)
         R = RINV * R2
         IF(R2 >= RRCUT1_SQ .OR. R2 < DRCUT_SQ) THEN
           CYCLE
         ENDIF
         ARD  = R - DRCUT
         ARD2 = ARD*ARD*RRCUT2
         ESW1 = ONE - ARD2
         ESW2 = ESW1*ESW1
         ESW  = ONE - ESW2
         VPCE( I,J,K ) = VPCE( I,J,K )
     *                   - ESW * CHARGE2(iin)*(ERF(ALPHA*R)*RINV)
        ENDDO
       ENDDO
       ENDDO
      ENDDO

      TIME2 = MPI_WTIME()

      IF(MYID.EQ.0) THEN
        WRITE(*,99) 'Elapsed Time (vpce) = ',TIME2-TIME1,' scnds'
      ENDIF

      VPCZ=0.D0
      PLJ =0.D0

      IF(MYID.EQ.0) THEN

      DO 80 in=1,nlmm
      DO 90 M=1,NNUC 

        RX=PNUC( M,1 ) - SCRD( 1,in )
        RY=PNUC( M,2 ) - SCRD( 2,in )
        RZ=PNUC( M,3 ) - SCRD( 3,in )
        R =DSQRT(RX**2+RY**2+RZ**2)
        RI = 1.D0/R

        sigr6  = (SQMMM(M,in)*RI*RI)**3
        sigr12 = sigr6*sigr6

        PLJ  = PLJ + EPQMMM(M,in) * ( sigr12 - sigr6 )    ! LJ

90    CONTINUE
80    CONTINUE
      
      TIME3 = MPI_WTIME()
      WRITE(*,99) 'Elapsed Time (plj)  = ',TIME3-TIME2,' scnds'

      DO 85 iin=1,nion1
        in = iiont(iin)
      DO 95 M=1,NNUC 

        RX=PNUC( M,1 ) - SCRD( 1,in )
        RY=PNUC( M,2 ) - SCRD( 2,in )
        RZ=PNUC( M,3 ) - SCRD( 3,in )
        R =DSQRT(RX**2+RY**2+RZ**2)
        RI = 1.D0/R

        IF(NMMSW(in).EQ.1) GOTO 95

        IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
          ESW = 1.D0
        ELSEIF(R.GT.DRCUT) THEN
          ARD  = R - DRCUT
          ARD2 = ARD*ARD*RRCUT2
          ESW1 = 1.D0 - ARD2
          ESW2 = ESW1*ESW1
          ESW  = 1.D0 - ESW2
        ELSE
          GOTO 95
        ENDIF

        VPCZ = VPCZ + ZVAL( M )*ESW*chgmm( in )*RI   ! Coulomb (nuclear-site)

95    CONTINUE
85    CONTINUE

      TIME4 = MPI_WTIME()
      WRITE(*,99) 'Elapsed Time (vpcz) = ',TIME4-TIME3,' scnds'
      
      WRITE(*,*)
      WRITE(*,*)'Nuclear - Site coulomb energy'
      WRITE(*,*)'  VPCZ = ' ,VPCZ
      WRITE(*,*)
      WRITE(*,*)'Nuclear - Site LJ potential energy'
      WRITE(*,*)'  PLJ  = ' ,PLJ

      ENDIF

C     WRITE(*,*) 'rank',myid,'is alive.'

      CALL MPI_BCAST( VPCZ,1,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )
      CALL MPI_BCAST( PLJ,1,MPI_DOUBLE_PRECISION,0,
     *                MPI_COMM_WORLD,IERR )

99    FORMAT( X,A22,F10.6,3X,A6 )

      RETURN
      END



C----------------------------
C   SUBROUTINE RAYLEIGH-RITZ 
C----------------------------

      SUBROUTINE RITZ_S( RWF,RWFN,NOR,VEFF,ORBE ) 

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX =  80 )
C     PARAMETER ( NDIM =   3 )
C     PARAMETER ( NORA =  49 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( RCUT =  3.5D0)
C     PARAMETER ( NSL  =  9000 )
C     PARAMETER ( NCYCLE  = 100 )

      PARAMETER ( PI   = 3.14159265358979323D0 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX, DY, DZ
      COMMON / NCLR  / PNUC( NNUC,NDIM )
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PSPOT / DRLOG(100),SVS(100),PVP(100),
     *                 RAD( 100,421),
     *                 VNLS(100,421),VNLP(100,421),
     *                 VLD( 100,421),VLDC(100,421)

C     DIMENSION RWF(   NMAX,NMAX,NMAX,NORA )            ! trial wave function
C     DIMENSION RWFN(  NMAX,NMAX,NMAX,NORA )            ! true wave function
      DIMENSION RWF(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NOR )    ! trial wave function
      DIMENSION RWFN(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NOR )    ! true wave function
      DIMENSION VEFF(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION ORBE(  NORA )
      DIMENSION ORBE1( NORA )
      DIMENSION HOR(   NOR,NOR )
      DIMENSION HORA(  NOR,NOR )
      DIMENSION FAI(   NOR,NOR )
      DIMENSION TGNALM(NNUC,2,3,NORA )         ! NUM OF ELCTRNS,ATOMS,ANGULAR MOMENTUM QUANTUM NUMBERS
      DIMENSION DGNALM(NNUC,  5,NORA )         ! NUM OF ELCTRNS,ATOMS,ANGULAR MOMENTUM QUANTUM NUMBERS
      DIMENSION NCT(   NNUC )

      DIMENSION AL(NOR)
      DIMENSION BE(NOR)
      DIMENSION CO(NOR)
      DIMENSION W( NOR)
      DIMENSION P( NOR)
      DIMENSION Q( NOR)
      DIMENSION B2(NOR)
      DIMENSION DI(NOR)
      DIMENSION BL(NOR)
      DIMENSION BU(NOR)
      DIMENSION BV(NOR)
      DIMENSION CM(NOR)
      DIMENSION LEX(NOR)

      DIMENSION WORK(NOR * 26)                  ! LAPACK work space, 26n
      DIMENSION IWORK(NOR * 10)                 ! ibid (10n)
      DIMENSION ISUPPZ(NOR * 2)

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      DV   = DX*DY*DZ

C----- ESTIMATE KINETIC ENERGY TERM -----

      CALL KNTC( RWF,NOR,RWFN )

C----- ESTIMATE EFFECTIVE POTENTIAL (NON-LOCAL POTENTIAL) -----

C     CALL GNALM( RWF,TGNALM,NOR )
      CALL GNALM( RWF,TGNALM,DGNALM,NOR )

C----- OPERATE HAMILTONIAN --------

      CALL OPERATE(RWF,RWFN,NOR,VEFF,TGNALM,DGNALM)

C===== CALCULATE MATRIX ELEMENT =====

C----- INITIALIZE -----

      DO NO = 1,NOR
      DO MO = 1,NOR

      HORA(NO,MO) = 0.D0

      ENDDO
      ENDDO

      CALL DGEMM('T', 'N', NOR, NOR, JX * JY * JZ,
     *      DV, RWF(1, 1, 1, 1), JX * JY * JZ,
     *          RWFN(1, 1, 1, 1), JX * JY * JZ,
     *   0.0D0, HORA(1, 1), NOR)

      NNDAT = NOR*NOR
      CALL MPI_REDUCE(HORA(1,1),HOR(1,1),NNDAT,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
C     CALL MPI_BCAST(HOR(1,1),NNDAT,MPI_DOUBLE_PRECISION,
C    *                0,MPI_COMM_WORLD,IERR)

C===== GENERATE RITZ VECTOR =====

      ! DSYEVR may return different values if thread parallel version is
      ! used. Thus the diagonalization is performed on rank0 only.
      IF(MYID.EQ.0) THEN
         CALL DSYEVR('V', 'A', 'U', NOR, HOR, NOR,
     *      0.0D0, 0.0D0, 0, 0, 0.0D0, 
     *      MFOUND, ORBE, FAI, NOR, ISUPPZ, 
     *      WORK, SIZE(WORK, 1), IWORK, SIZE(IWORK, 1), INFO)
         IF(INFO /= 0) THEN
            print *, "info = ", info
            STOP "ERROR(RITZ): DSYEVR failed"
         ENDIF
         IF(MFOUND /= NOR) STOP "ERROR: DSYEVR"
      ENDIF

      CALL MPI_BCAST(FAI(1,1),NNDAT,MPI_DOUBLE_PRECISION,
     *                0,MPI_COMM_WORLD,IERR)

C===== GENERATE NEW WAVE FUNCTIONS =====

      CALL DGEMM('N', 'N', JX * JY * JZ, NOR, NOR,
     *   1.0D0, RWF(1, 1, 1, 1), JX * JY * JZ,
     *          FAI(1, 1), NOR,
     *   0.0D0, RWFN(1, 1, 1, 1), JX * JY * JZ)

C     RWFN is already normalized (because FAI is unitary and RWF is
C     normalized). No need to renormalize.

      RETURN
      END


C-------------------------------------------
C     SUBROUTINE QM/MM FORCE
C-------------------------------------------

      SUBROUTINE RQMF1_S( SCRD,FRCS,FRCN )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include "mpi.i"

      include 'QMpara.i'                        ! Vmol
!     include 'sizes.i'
!     include 'atoms.i'

C     PARAMETER ( NMAX = 80 )
C     PARAMETER ( NDIM =  3 )
C     PARAMETER ( NNUC = 36 )
      PARAMETER ( DCUT  = 1.D-9 ) 
      PARAMETER ( ALPHA = 1.0D0 )
C     PARAMETER ( NLINK1 = 10 )
      PARAMETER ( RRCUT = 4.724D0 )
      PARAMETER ( RRCUT1= 5.669D0 )
      PARAMETER ( DRCUT = 0.945D0 )

      PARAMETER ( PI    = 3.14159265358979323D0 )
      
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ

      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT2 / NMM2,NLINK,NLAQM(NLINK1),NLAMM(NLINK1),
     *                 NMMSW(maxatm),MMID(NNUC)
      COMMON / PRMT4 / nion1,iiont(maxatm)
      COMMON / NCLR / PNUC( NNUC,NDIM )
      COMMON / LJQMMM / SQMMM(NNUC,maxatm),EPQMMM(NNUC,maxatm),
     *                  chgmm(maxatm),chgqm(NNUC)

      COMMON / RDDNST /RRHO(NMAXX/(2*NX),NMAXY/(2*NY),NMAXZ/(2*NZ))

      REAL*8 CHARGE0, CHARGE2, CCRD0, CCRD2, RINV, PVPCE
      REAL*8 DRCUT_SQ, RRCUT1_SQ, R2T
      INTEGER NATM0, NATM2

      DIMENSION SCRD( NDIM,maxatm )

      DIMENSION FDSS(  maxatm,NDIM )              ! density-site force on site
      DIMENSION FODSS( maxatm,NDIM )              ! density-site force on site
      DIMENSION FNSS(  maxatm,NDIM )              ! nuclear-site and LJ force on site
      DIMENSION FRCS(  maxatm,NDIM )              ! QM/MM force on site

      DIMENSION FRCN( NNUC,NDIM )                 ! QM/MM force on nuclear

      DIMENSION ADX( NMAXX/(NX*2) )
      DIMENSION ADY( NMAXY/(NY*2) )
      DIMENSION ADZ( NMAXZ/(NZ*2) )

C     DIMENSION CHARGE0(maxatm), CHARGE2(maxatm)
C     DIMENSION CCRD0( maxatm, NDIM ), CCRD2( maxatm, NDIM )

      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MYID,IERR )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NUMPROCS,IERR )

      IF(MYID.EQ.0) THEN
        WRITE(*,*) 'Start qmf1 using reduced density'
      ENDIF

      DV = 8.D0*DX*DY*DZ

      RRCUT1_SQ = RRCUT1*RRCUT1
      DRCUT_SQ  = DRCUT *DRCUT

      JX2 = JX/2
      JY2 = JY/2
      JZ2 = JZ/2

C
C     number of loop in QM and MM subsystems
C
      nlmm = n-NNUC+NLINK

      DO L  = 1,3
      DO in = 1,n

        FDSS(  in,L ) = 0.D0
        FODSS( in,L ) = 0.D0
        FNSS(  in,L ) = 0.D0
        FRCS(  in,L ) = 0.D0
            
      ENDDO
      ENDDO
     
      DO L=1,3
      DO in=1,NNUC
     
        FRCN( in,L )   = 0.D0
            
      ENDDO
      ENDDO

C-----density-site force CALCULATION
     
      TIME1 = MPI_WTIME()

      ONE    = 1.D0
      TWO    = 2.D0
      FOUR   = 4.D0
      RRCUT2 = ONE/(RRCUT*RRCUT)
      SQPII  = ONE/DSQRT(PI)
      A2SQPI = TWO*ALPHA*SQPII

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO K=1,JZ2
        ADZ(K) = DZ*( 2*K-NNZ ) - DZ/2.D0
      ENDDO
      DO J=1,JY2
        ADY(J) = DY*( 2*J-NNY ) - DY/2.D0
      ENDDO
      DO I=1,JX2
        ADX(I) = DX*( 2*I-NNX ) - DX/2.D0
      ENDDO

C---- count atoms up for each type -----

C     NATM0 = 0
C     NATM2 = 0
C     DO iin = 1, nion1
C       in = iiont(iin)
C       IF(NMMSW(in).EQ.0) THEN
C         NATM0 = NATM0 + 1
C         CCRD0(NATM0, 1) = SCRD(1, in)
C         CCRD0(NATM0, 2) = SCRD(2, in)
C         CCRD0(NATM0, 3) = SCRD(3, in)
C         CHARGE0(NATM0) = chgmm(in)
C       ELSEIF(NMMSW(in).EQ.2) THEN
C         NATM2 = NATM2 + 1
C         CCRD2(NATM2, 1) = SCRD(1, in)
C         CCRD2(NATM2, 2) = SCRD(2, in)
C         CCRD2(NATM2, 3) = SCRD(3, in)
C         CHARGE2(NATM2) = chgmm(in)
C       ENDIF
C     ENDDO

C---------------------------------------

      DO 50 iin=1,nion1

      in = iiont(iin)
      IF(NMMSW(in).EQ.1) CYCLE

      DO 60 K=1,JZ2
       RZ = ADZ(K) - SCRD( 3,in )
      DO 70 J=1,JY2
       RY = ADY(J) - SCRD( 2,in )
      DO 80 I=1,JX2

       ARHO = RRHO(I,J,K)

C     IF( ARHO.GT.DCUT ) THEN

       RX = ADX(I) - SCRD( 1,in )
       R2 = RX**2+RY**2+RZ**2

C      IF( NMMSW(in).EQ.0 .OR. R2.GE.RRCUT1_SQ ) THEN
       IF( NMMSW(in).EQ.0 ) THEN

       IF( R2 >= RRCUT1_SQ ) THEN

         R  = DSQRT( R2 )

         RI  = 1.D0/R
         R2I = RI*RI
         R3I = RI*R2I

         TFDSS = -ARHO*chgmm(in)*R3I

       ELSEIF( R2 >= DRCUT_SQ .AND. R2 <= RRCUT1_SQ ) THEN

         R     = DSQRT( R2 )
         RI    = ONE/R
         R2I   = RI*RI
         R3I   = RI*R2I
         APR   = ALPHA*R
         C2RDV = ARHO * chgmm( in ) 

         TFDSS = C2RDV*(A2SQPI*DEXP(-APR**2)*R2I-ERF(APR)*R3I)

       ENDIF

       ELSEIF( NMMSW(in).EQ.2 ) THEN

       IF( R2 >= RRCUT1_SQ ) THEN

         R  = DSQRT( R2 )

         RI  = 1.D0/R
         R2I = RI*RI
         R3I = RI*R2I

         TFDSS = -ARHO*chgmm(in)*R3I

       ELSEIF( R2 >= DRCUT_SQ .AND. R2 <= RRCUT1_SQ ) THEN
         R     = DSQRT( R2 )
         ARD   = R - DRCUT
         ARD2  = ARD*ARD*RRCUT2
         ESW1  = ONE - ARD2
         ESW2  = ESW1* ESW1
         ESW   = ONE - ESW2
         ESWA  = FOUR*ESW1*ARD*RRCUT2
         RI    = ONE/R
         R2I   = RI*RI
         APR   = ALPHA*R
         C2RDV = ARHO * chgmm(in) 
         DER2I = ERF(APR)*R2I
         TFDSS = ESW* C2RDV*(A2SQPI*DEXP(-APR**2)*R2I-DER2I*RI)
     *         + ESWA*C2RDV*DER2I
       ENDIF

       ELSE
        CYCLE
       ENDIF

       FODSS(in,1) = FODSS(in,1) - RX*TFDSS*DV
       FODSS(in,2) = FODSS(in,2) - RY*TFDSS*DV
       FODSS(in,3) = FODSS(in,3) - RZ*TFDSS*DV

C     ENDIF

80    CONTINUE
70    CONTINUE
60    CONTINUE

50    CONTINUE
 
C     NUMDAT = n*NDIM
      NUMDAT = maxatm*NDIM
      CALL MPI_REDUCE( FODSS(1,1),FDSS(1,1),NUMDAT,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )
C     WRITE(*,*) 'MYID:qmf1=',MYID,RHO(5,5,5)

      TIME2 = MPI_WTIME()

      IF(MYID.EQ.0) THEN
        WRITE(*,99) 'Elapsed Time (qmf1) = ',TIME2-TIME1,' scnds'
      ENDIF

C-----nuclear-site force CALCULATION
C-----LJ force CALCULATION

      IF(MYID.EQ.0) THEN

      SIX   = 6.D0

      DO 90 in=1,nlmm
      DO 100 M=1,NNUC 

        RX=PNUC( M,1 ) - SCRD( 1,in )
        RY=PNUC( M,2 ) - SCRD( 2,in )
        RZ=PNUC( M,3 ) - SCRD( 3,in )
        R2=RX**2+RY**2+RZ**2
        R   = DSQRT(R2)
        RI  = ONE/R
        R2I = RI*RI

        ASIG2  = SQMMM(M,in)*R2I
        ASIG6  = ASIG2*ASIG2*ASIG2
        ASIG12 = ASIG6*ASIG6
        BSIG   = TWO*ASIG12 - ASIG6
        TFEPS  = SIX*EPQMMM(M,in)
        CREPS  = TFEPS*R2I*BSIG
        CREPSX = CREPS*RX
        CREPSY = CREPS*RY
        CREPSZ = CREPS*RZ

        FNSS(in,1) = FNSS(in,1) + CREPSX   ! LJ force on site
        FNSS(in,2) = FNSS(in,2) + CREPSY
        FNSS(in,3) = FNSS(in,3) + CREPSZ

        FRCN(M,1)  = FRCN(M,1) + CREPSX    ! LJ force on nuclear
        FRCN(M,2)  = FRCN(M,2) + CREPSY
        FRCN(M,3)  = FRCN(M,3) + CREPSZ

100   CONTINUE
90    CONTINUE
      
      TIME3 = MPI_WTIME()
        WRITE(*,99) 'Elapsed Time (qmf2) = ',TIME3-TIME2,' scnds'

      DO 95 iin=1,nion1
        in = iiont(iin)
      DO 105 M=1,NNUC 

        RX=PNUC( M,1 ) - SCRD( 1,in )
        RY=PNUC( M,2 ) - SCRD( 2,in )
        RZ=PNUC( M,3 ) - SCRD( 3,in )
        R2=RX**2+RY**2+RZ**2
        R   = DSQRT(R2)
        RI  = ONE/R
        R2I = RI*RI

       IF(NMMSW(in).EQ.1) GOTO 105

       IF(NMMSW(in).EQ.0 .OR. R.GE.RRCUT1) THEN
         R3I   = RI*R2I
         ACHR  = ZVAL( M )*chgmm( in )*R3I
       ELSEIF(R.GT.DRCUT) THEN
         ARD   = R - DRCUT
         ARD2  = ARD*ARD*RRCUT2
         ESW1  = ONE - ARD2
         ESW2  = ESW1*ESW1
         ESW   = ONE - ESW2
         ESWA  = FOUR*ESW1*ARD*RRCUT2
         ZVPCH = ZVAL( M )*chgmm( in )*R2I
         ACHR  = ESW*ZVPCH*RI - ESWA*ZVPCH
       ELSE
         GOTO 105
       ENDIF

        CHRX = ACHR*RX
        CHRY = ACHR*RY
        CHRZ = ACHR*RZ

        FNSS(in,1) = FNSS(in,1) + CHRX      ! Coulomb force on site
        FNSS(in,2) = FNSS(in,2) + CHRY
        FNSS(in,3) = FNSS(in,3) + CHRZ

        FRCN(M,1)  = FRCN(M,1) + CHRX       ! Coulomb force on nuclear
        FRCN(M,2)  = FRCN(M,2) + CHRY
        FRCN(M,3)  = FRCN(M,3) + CHRZ

105   CONTINUE
95    CONTINUE

      TIME4 = MPI_WTIME()
        WRITE(*,99) 'Elapsed Time (qmf3) = ',TIME4-TIME3,' scnds'

      DO L=1,3
      DO in=1,nlmm
        FRCS(in,L) = FDSS(in,L) - FNSS(in,L)                 !QM/MM force on site
      ENDDO
      ENDDO

      ENDIF

99    FORMAT( X,A22,F10.6,3X,A6 )

      RETURN
      END



C----------------------------------------------
C
C   SUBROUTINE HARTREE-FOCK EXCHANGE IN REAL SPACE
C   USING POISSON EQUATION FOR PARALLEL IMPLEMENTATION
C
C   2018.03.09 
C
C----------------------------------------------

C----- RVCLMB is modified -----

      SUBROUTINE HF_EXCH_PSN( RWF,RWFT,NOR,MOR,VHIJX,HFEX )
      USE POISSON_SOLVER,
     *     PSOLVER_RUN => RUN

      IMPLICIT REAL*8( A-H,O-Z )
      IMPLICIT INTEGER*4( I-N )

      include 'mpif.h'
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'

      PARAMETER( INUM  = 8      )
      PARAMETER( INUM2 = INUM/2 )
      PARAMETER( PI    = 3.14159265358979323D0 )
      PARAMETER( NITMAX = 100 )
      PARAMETER( NCG    =  10 )
      PARAMETER( NLR    =   5 )
      PARAMETER( EPS1 = 1.0D0 )
C     PARAMETER( EPS2 = 1.0D-05 )
      PARAMETER( EPS2 = 1.0D-04 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          FREEZE*5,PRINT*5,DGF*3

      DIMENSION TA(  3 )
      DIMENSION TA1( 3 )

      DIMENSION CTR( 2,3,NORA*NORA ) 
      DIMENSION SW(  2,NORA*NORA ) 

      DIMENSION RWF(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION RWFT(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
!     REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: RWF
!     REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: RWFT
      DIMENSION RHOIJ( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DMMY(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VHIJ(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VHIJX( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
!     REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: VHIJX

      DIMENSION RHO(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION BRHO(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VCOUX( NMAXY/NY,NMAXZ/NZ,INUM ) ! note: Dimension(Y,Z,X)
      DIMENSION VCOUY( NMAXX/NX,NMAXZ/NZ,INUM ) ! note: Dimension(X,Z,Y)
      DIMENSION VCOUZ( NMAXX/NX,NMAXY/NY,INUM )
      DIMENSION VBC  ( 1-INUM2:NMAXX/NX+INUM2,1-INUM2:NMAXY/NY+INUM2,
     *                 1-INUM2:NMAXZ/NZ+INUM2)
      DIMENSION PVEC(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RVEC(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION APVEC( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHOL(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION PHIL(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / CLMB  / PCLMB, PEXC, PEXT
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX,DY,DZ
      COMMON / OFC / VCOU( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

C     OPEN(75,FILE='TMP.DAT',STATUS='UNKNOWN')

      DV = DX*DY*DZ

C     write(*,*)  MYID, ' enter: hf_exch_psn '

C     CALL CNTROID( RWF,CTR,NOR,SW )

      RWFT(:,:,:,:) = 0.D0
      DMMY( :,:,:)  = 0.D0

C     NCNT = 0
C     ID   = 0

C     CALL LRCLMB0

      DO 50 IOC = 1,  NOR     ! loop over orbital pairs
C     DO 70 JOC = IOC,NOR
      DO 70 JOC = 1,  MOR     ! Do loop is changed so that it can be used for virtual orbitals.

C     NCNT = NCNT + 1
C     IF(MOD(NCNT,NUMPROCS) .EQ. MYID) THEN
C     ID = ID + 1

C----- BOUNDARY CONDITION -----

C     IF( IOC .EQ. JOC ) THEN
C       CALL BNDRY( CTR,SW,VBC,NCNT )
C     ELSE
C       CALL BNDRY( CTR,NOR,IOC,VBC,0 )
C     ENDIF

C     SUM = 0.D0
      DO K = 1,JZ
      DO J = 1,JY
      DO I = 1,JX
        RHOIJ(I,J,K) = RWF(I,J,K,IOC)*RWF(I,J,K,JOC)
C       SUM = SUM + RHOIJ(I,J,K)*DV
      ENDDO
      ENDDO
      ENDDO

      DMMY(:,:,:) = 0.D0
      CALL EFTC(RHOIJ,DMMY) 

C     IF(IOC /= JOC) THEN
      IF( .TRUE. ) THEN      ! for inc_i=j

C      IF( MOD(NJ,NLR) == 1 ) THEN
         CALL LRCLMB(RHOL,PHIL)    ! for long-range part of density and potential
C      ENDIF
C        RHO(:,:,:) = RHO(:,:,:) - 4.D0*PI*RHOL(:,:,:)
         PHI_OR_RHO(1:JX, 1:JY, 1:JZ) = -4.D0*PI*(RHOIJ(:, :, :) - RHOL(:,:,:))
         CALL PSOLVER_RUN()
         VHIJ(:,:,:) = PHI_OR_RHO(1:JX, 1:JY, 1:JZ)
         VHIJ(:,:,:) = VHIJ(:,:,:) + PHIL(:,:,:)  ! sum the long-range part
         GOTO 60
      ENDIF

C     CALL BNDRY(  VBC )
      CALL RCLMB2( VBC )
      PVEC(:,:,:) = 0.D0
      CALL RCLMB5( PVEC,VBC,BRHO )
      BRHO(:,:,:) = 4.D0*PI*RHOIJ(:,:,:) - BRHO(:,:,:)
   
C     IF( IOC .NE. JOC ) THEN
C       VBC(:,:,:) = 0.D0
C     ENDIF

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        TMP = VHIJX(I,J,K,IOC) 
C       TMP = VHIJX(I,J,K,ID) 
        VHIJ(I,J,K) = TMP                 ! initial guess 
      ENDDO
      ENDDO
      ENDDO

C----- STEEPEST DESCENT -----

10    CONTINUE

C     write(*,*) 'steepest descent'

      VBC(:,:,:) = 0.D0
      CALL RCLMB5( VHIJ,VBC,APVEC )

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        RVEC(I,J,K) = BRHO(I,J,K) - APVEC(I,J,K)
        PVEC(I,J,K) = RVEC( I,J,K ) 
      ENDDO
      ENDDO
      ENDDO
  
C----- CONJUGATE GRADIENT -----
C----- ITERATION -----

C     write(*,*) 'conjugate gradient'

C     IF( IOC .NE. JOC ) THEN
C       VBC(:,:,:) = 0.D0
C     ENDIF

      DO 20 L = 1, NITMAX
  
C       CALL RCLMB4( PVEC,APVEC )
        CALL RCLMB5( PVEC,VBC,APVEC )
  
        TA1(1) = 0.D0
        TA1(2) = 0.D0
        DO K = 1, JZ
        DO J = 1, JY
        DO I = 1, JX
          TA1(1) = TA1(1) + RVEC( I,J,K )*RVEC(  I,J,K )
          TA1(2) = TA1(2) + PVEC( I,J,K )*APVEC( I,J,K )
C         TA2 = A2 + RVEC( I,J,K )*APVEC( I,J,K )       ! 2015.07.16 takahasi
        ENDDO
        ENDDO
        ENDDO
        
        CALL MPI_ALLREDUCE( TA1(1),TA(1),2,MPI_DOUBLE_PRECISION,
     *                      MPI_SUM,MPI_COMM_WORLD,IERR )

        ALPHA = TA(1) / TA(2) 
  
        TA1(3) = 0.D0
        DO K = 1, JZ
        DO J = 1, JY
        DO I = 1, JX
          VHIJ( I,J,K ) = VHIJ( I,J,K ) + ALPHA*PVEC(  I,J,K )
          RVEC( I,J,K ) = RVEC( I,J,K ) - ALPHA*APVEC( I,J,K )
          TA1(3) = TA1(3) + RVEC( I,J,K )**2
        ENDDO
        ENDDO
        ENDDO
  
        CALL MPI_ALLREDUCE( TA1(3),TA(3),1,MPI_DOUBLE_PRECISION,
     *                      MPI_SUM,MPI_COMM_WORLD,IERR )

        DIST = DSQRT(TA(3)*DV)
C       IF(MOD(L,10).EQ.1) WRITE(*,*) 'DIST = ',L,DIST
C       IF( DIST .LT. EPS2 .OR. DIST .GE. EPS1 ) THEN
        IF( DIST .LT. EPS2 ) THEN
          GOTO 60
        ENDIF
       
        BETA = TA(3)/TA(1)
  
        DO K = 1, JZ
        DO J = 1, JY
        DO I = 1, JX
          PVEC( I,J,K ) = RVEC( I,J,K ) + BETA*PVEC( I,J,K )
        ENDDO
        ENDDO
        ENDDO
  
20    CONTINUE
C     WRITE(*,*) 'DIST = ',L,DIST

60    CONTINUE
C     IF( DIST .GE. EPS1 ) GOTO 10
C     IF(PRINT.EQ.'LARGE') THEN
C       WRITE(*,*)                     
C       WRITE(*,*) ' Poisson It,Dist = ',L-1,DIST
C     ENDIF

      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        RWFT(I,J,K,IOC) = RWFT(I,J,K,IOC)
     *                  - VHIJ(I,J,K)*RWF(I,J,K,JOC)
      ENDDO
      ENDDO
      ENDDO
  
      IF( JOC .GT. IOC ) THEN
      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
C       RWFT(I,J,K,JOC) = RWFT(I,J,K,JOC)
C    *                  - VHIJ(I,J,K)*RWF(I,J,K,IOC)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
     
      IF( IOC == JOC ) THEN
      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        TMP = VHIJ(I,J,K) 
        VHIJX(I,J,K,IOC) = TMP        ! replace the guess for next scf step
      ENDDO
      ENDDO
      ENDDO
      ENDIF

70    CONTINUE
50    CONTINUE

C     DO J = 1, NMAXY
C     DO I = 1, NMAXX
C       WRITE(75,*) VHIJ(I,J,NMAXZ/2)
C     ENDDO
C     ENDDO

C     NDAT = NMAXX*NMAXY*NMAXZ*NORA
C     CALL MPI_REDUCE(RWFT1(1,1,1,1),RWFT(1,1,1,1),NDAT,
C    *         MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)

C----- TOTAL HF EXCHANGE ENERGY -----   

      THFEX = 0.D0
C     DO IOC = 1,NOR      ! Do loop should run only over the occ. orbs. 
C     DO IOC = 1,MORA     ! 20230620
      DO IOC = 1,MOR      ! 20230804
      DO K   = 1,JZ
      DO J   = 1,JY
      DO I   = 1,JX
        THFEX = THFEX + RWFT(I,J,K,IOC)*RWF(I,J,K,IOC)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_ALLREDUCE( THFEX,HFEX,1,MPI_DOUBLE_PRECISION,
     *                      MPI_SUM,MPI_COMM_WORLD,IERR )

      HFEX = 0.5D0*HFEX*DV

      IF(MYID .EQ. 0) THEN
        WRITE(*,*) ' HFEX  = ', HFEX
      ENDIF

C     REWIND( 77 )

C     DO I = 1, NMAX
C     DO J = 1, NMAX 
C     WRITE(77,*) VCOU( I,J,NMAX/2 )
C     ENDDO
C     ENDDO

      RETURN
      END


C---------------------------------------------------
C
C  Construct initial guess w.f. from the lcao coefficents
C  (fort.7) yielded by Gaussian package.
C  input: basis.dat    <-- #gfinput option
C         ( extracted from the output file )
C         valence.dat = fort.7 <-- #Punch(MO) option
C         ( core w.f. eliminated )
C
C  orginal:   H. Kambe
C  bug fixed and parallelized: H. Takahashi 
C   
C---------------------------------------------------

      SUBROUTINE TRWF3( CRD,WFA,WFB,WFNA,WFNB )
       
      IMPLICIT REAL*8 ( A-H,O-Z )
      IMPLICIT INTEGER*4 ( I-N )
      
      include "mpif.h"                          ! mpi
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'                        ! Vmol

C     PARAMETER ( NMAX =  80 )
C     PARAMETER ( NNUC =  36 )
C     PARAMETER ( NORA =  49 )
C     PARAMETER ( NORB =   1 )
C     PARAMETER ( MORA =  49 )
C     PARAMETER ( MORB =   1 )
      PARAMETER ( PI   =   3.14159265358979323D0 )
C     PARAMETER ( ALPHA = 3.D-1 )

      CHARACTER NRST*7,NOPT*3,EXC*7,NQMMM*4,
     *          PRINT*5,DGF*3,FREEZE*5

      CHARACTER ATMORB*5
      CHARACTER ATMINX*2
      
      CHARACTER NAME1*11,NAME2*13

      COMMON / GRID2 / DX, DY, DZ
      COMMON / PRMT  / COE,TEMP,DTMD,
     *                 NDEN,MDMAX,
     *                 ZA(NNUC),ZVAL(NNUC),
     *                 SIG(NNUC),EPSQM(NNUC),   ! Vmol
     *                 NRVLC,NCHK,CONV
      COMMON / PRMT1 / NRST,NOPT,EXC,NQMMM,
     *                 FREEZE,PRINT,DGF
      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / AVRHO / RHOAV(NMAXX/NX,NMAXY/NY,NMAXZ/NZ)

      DIMENSION WFA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION WFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
C     DIMENSION WFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,1 )
      DIMENSION WFNA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION WFNB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION NRDA( NORA )
      DIMENSION NRDB( NORB )
      
      DIMENSION CRD(    NNUC,NDIM )
      DIMENSION COEMOA( NBASIS,NBASIS )
      DIMENSION COEMOB( NBASIS,NBASIS )
      DIMENSION BASIS(  NBASIS )
      DIMENSION BASINF( 5,20,3 )
      DIMENSION NPRIM( 20)
      DIMENSION NTYP(  20)

      NAMELIST/INIDAT/COE,NDEN,DTMD,MDMAX,TEMP,
     *                NRST,NOPT,EXC,NQMMM,NRVLC,
     *                NCHK,CONV,FREEZE,PRINT,
C    *                DGF,NMRDF                  
     *                DGF,NMRDF,NMM2,NLINK          ! Vmol

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      NDAT = NNUC*NDIM

      CALL MPI_BCAST( CRD(1,1), NDAT,MPI_DOUBLE_PRECISION,0,   
     *                               MPI_COMM_WORLD,IERR )

C     if(myid.eq.5) then
C       write(*,*) crd(3,1) 
C     endif

C----- open avrho.dat file -----

      write (NAME2, '("avrho.dat", i4.4)') MYID                 ! takahasi 2013.07.30
      OPEN(24,FILE=NAME2,STATUS='UNKNOWN',FORM='UNFORMATTED')   ! takahasi 2013.07.30

C---- COMPUTE PARAMETERS ----

      DV = DX*DY*DZ

      OPEN(75,FILE ='basis.dat',STATUS='OLD')       ! information for the basis set given by gfinput option
C     OPEN(86,FILE ='fort.7',   STATUS='OLD')       ! lcao coefficients generated by punch=mo option
      OPEN(86,FILE ='valence.dat',STATUS='OLD')     ! lcao coefficients generated by punch=mo option

      IF(MYID .EQ. 0) THEN

      IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                 .OR. EXC.EQ.'RXalpha') THEN

         READ(86,*)

         DO J=1,NORA
            READ(86,'()') 
            READ(86,'(5D15.8)') (COEMOA(I,J),I=1,NBASIS)
         ENDDO

      ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                     .OR. EXC.EQ.'Uxalpha') THEN

         READ(86,*)

         DO J=1,NORA
            READ(86,'()') 
            READ(86,'(5D15.8)') (COEMOA(I,J),I=1,NBASIS)
         ENDDO

         DO J=1,NBASIS-NORA-NCORE     ! 20230802 takahasi
            READ(86,'()')    
            READ(86,'(5D15.8)') (DUMMY,I=1,NBASIS)     ! The orbitals \phi_n with n > nora will be discarded.
         ENDDO

         READ(86,*)

         DO J=1,NORB
            READ(86,'()') 
            READ(86,'(5D15.8)') (COEMOB(I,J),I=1,NBASIS)
         ENDDO

      ENDIF

      ENDIF

      IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                 .OR. EXC.EQ.'RXalpha') THEN
        NUM = NORA*NBASIS
        CALL MPI_BCAST(COEMOA(1,1),NUM,MPI_DOUBLE_PRECISION,0,
     *                 MPI_COMM_WORLD,IERR)
      ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                      .OR. EXC.EQ.'Uxalpha') THEN
        NUM = NORA*NBASIS
        CALL MPI_BCAST(COEMOA(1,1),NUM,MPI_DOUBLE_PRECISION,0,
     *                 MPI_COMM_WORLD,IERR)
        NUM = NORB*NBASIS
        CALL MPI_BCAST(COEMOB(1,1),NUM,MPI_DOUBLE_PRECISION,0,
     *                 MPI_COMM_WORLD,IERR)
      ENDIF

C     if(myid .eq. 3) then
C       write(*,*) COEMOA(10,4)
C     endif

      WFA(:,:,:,:) = 0.0D0
      WFB(:,:,:,:) = 0.0D0

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      INDX = 0
      DO N = 1,NNUC

       BASINF(:,:,:) = 0.D0
       NSHELL = 0

C     IF(MYID .EQ. 0) THEN

       READ(75,*) NATOM,IDUM
       IF( NATOM .NE. N ) THEN
         WRITE(*,*) 'Warining: natom is not coincident with n.'
         STOP
       ENDIF

C     ENDIF

       DO     ! loop over primitive basis in a shell

C      IF(MYID .EQ. 0) THEN

           BASINF(:,:,:) = 0.D0
           NTYP( :)      = 0
           NPRIM(:)      = 0

           READ(75,*) ATMORB
           IF(ATMORB .EQ. '****') THEN
             EXIT  
           ELSEIF(ATMORB .EQ. 'G') THEN
             STOP
           ENDIF

           BACKSPACE(75)
           READ(75,*) ATMORB,NGAUSS
           NSHELL = NSHELL + 1 
           NPRIM(NSHELL) = NGAUSS

           IF(ATMORB .EQ. 'S') THEN
             NTYP(NSHELL) = 1
             DO N1=1,NGAUSS
               READ(75,*) FACTOR,COEF
               BASINF(1,N1,1) = FACTOR
               BASINF(1,N1,2) = COEF
             ENDDO
           ELSEIF(ATMORB .EQ. 'SP') THEN
             NTYP(NSHELL) = 2
             DO N1=1,NGAUSS
               READ(75,*) FACTOR,COEFS,COEFP
               BASINF(2,N1,1) = FACTOR
               BASINF(2,N1,2) = COEFS
               BASINF(2,N1,3) = COEFP
             ENDDO
           ELSEIF(ATMORB .EQ. 'P') THEN
             NTYP(NSHELL) = 3
             DO N1=1,NGAUSS
               READ(75,*) FACTOR,COEF
               BASINF(3,N1,1) = FACTOR
               BASINF(3,N1,2) = COEF
             ENDDO
           ELSEIF(ATMORB .EQ. 'D') THEN
             NTYP(NSHELL) = 4
             DO N1=1,NGAUSS
               READ(75,*) FACTOR,COEF
               BASINF(4,N1,1) = FACTOR
               BASINF(4,N1,2) = COEF
             ENDDO
           ELSEIF(ATMORB .EQ. 'F') THEN
             NTYP(NSHELL) = 5
             DO N1=1,NGAUSS
               READ(75,*) FACTOR,COEF
               BASINF(5,N1,1) = FACTOR
               BASINF(5,N1,2) = COEF
             ENDDO
           ENDIF

C     ENDIF
C      ENDDO

C     CALL MPI_BCAST(BASINF(1,1,1),300,MPI_DOUBLE_PRECISION,0,
C    *               MPI_COMM_WORLD,IERR)
C     CALL MPI_BCAST(NTYP(1),20,MPI_INTEGER,0,
C    *               MPI_COMM_WORLD,IERR)
C     CALL MPI_BCAST(NPRIM(1),20,MPI_INTEGER,0,
C    *               MPI_COMM_WORLD,IERR)
C     CALL MPI_BCAST(NSHELL,1,MPI_INTEGER,0,
C    *               MPI_COMM_WORLD,IERR)

      if(myid .eq. 3) then
      write(*,*) BASINF(3,1,1), BASINF(3,2,1)
      write(*,*) N,NSHELL,NTYP(NSHELL),NPRIM(NSHELL)
      endif

C     DO N2 = 1,NSHELL
      N2 = NSHELL

      INDX  = INDX+1
      NNTYP = NTYP(N2)

C     write(*,*) NTYP(1),NPRIM(1),NSHELL

      DO K=1,JZ
         RZA = DZ*( K-NNZ ) - CRD(N,3) - 0.5D0*DZ      ! 2018.03.16 Takahashi
C        RZA = DZ*( K-NNZ ) - CRD(N,3) 
      DO J=1,JY                                        ! shifted by half of the grid size to fit raw orb.dat
         RYA = DY*( J-NNY ) - CRD(N,2) - 0.5D0*DY      ! This lowers the energy of the first step of the SCF cycle.
C        RYA = DY*( J-NNY ) - CRD(N,2) 
      DO I=1,JX
         RXA = DX*( I-NNX ) - CRD(N,1) - 0.5D0*DX
C        RXA = DX*( I-NNX ) - CRD(N,1) 

         SQR = RXA**2 + RYA**2 + RZA**2

!------- for s-type atomic orb. -----------

          IF(NNTYP .EQ. 1) THEN

            PNORM=(8.0D0/PI**3)**0.25D0
              
            BASIS(:) = 0
            DO N1 = 1,NPRIM(N2)
            
              FACTOR =  BASINF(1,N1,1) 
              COEF   =  BASINF(1,N1,2) 

              PRES=FACTOR**0.75D0*DEXP(-FACTOR*SQR)
              PRES=COEF*PNORM*PRES
              BASIS(INDX)=BASIS(INDX)+PRES

            ENDDO

            IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                  .OR. EXC.EQ.'RXalpha') THEN
             DO NO = 1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX)
     *                                    *COEMOA(INDX,NVAL) 
             ENDDO
            ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                      .OR. EXC.EQ.'Uxalpha') THEN
             DO NO=1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX)*COEMOA(INDX,NVAL) 
             ENDDO
             DO NO=1,NORB 
               NVAL =         NO
               WFB(I,J,K,NO)=WFB(I,J,K,NO)+BASIS(INDX)*COEMOB(INDX,NVAL) 
             ENDDO
            ENDIF
     
          ELSEIF(NNTYP .EQ. 2) THEN

            PNORMS=(8.0D0/PI**3)**0.25D0
            PNORMP=(128.0D0/PI**3)**0.25D0

            BASIS(:) = 0
            DO N1 = 1,NPRIM(N2)

                FACTOR =  BASINF(2,N1,1) 
                COEFS  =  BASINF(2,N1,2) 
                COEFP  =  BASINF(2,N1,3) 

                PRESP=FACTOR**0.75D0*DEXP(-FACTOR*SQR)
                PRESP=COEFS*PNORMS*PRESP
                BASIS(INDX)=BASIS(INDX)+PRESP

                PRESPX=FACTOR**1.25D0*(RXA)*DEXP(-FACTOR*SQR)
                PRESPX=COEFP*PNORMP*PRESPX
                BASIS(INDX+1)=BASIS(INDX+1)+PRESPX

                PRESPY=FACTOR**1.25D0*(RYA)*DEXP(-FACTOR*SQR)
                PRESPY=COEFP*PNORMP*PRESPY
                BASIS(INDX+2)=BASIS(INDX+2)+PRESPY

                PRESPZ=FACTOR**1.25D0*(RZA)*DEXP(-FACTOR*SQR)
                PRESPZ=COEFP*PNORMP*PRESPZ
                BASIS(INDX+3)=BASIS(INDX+3)+PRESPZ

            ENDDO
         
            IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                  .OR. EXC.EQ.'RXalpha') THEN
             DO NO = 1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,3
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                      .OR. EXC.EQ.'Uxalpha') THEN
             DO NO=1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,3
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
             DO NO=1,NORB 
               NVAL =         NO
               DO INC = 0,3
               WFB(I,J,K,NO)=WFB(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOB(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ENDIF

          ELSEIF(NNTYP .EQ. 3) THEN

            PNORMP=(128.0D0/PI**3)**0.25D0

            BASIS(:) = 0
            DO N1 = 1,NPRIM(N2)

              FACTOR =  BASINF(3,N1,1) 
              COEFP  =  BASINF(3,N1,2) 

                PREPX=FACTOR**1.25D0*(RXA)*DEXP(-FACTOR*SQR)
                PREPX=COEFP*PNORMP*PREPX
                BASIS(INDX)=BASIS(INDX)+PREPX

                PREPY=FACTOR**1.25D0*(RYA)*DEXP(-FACTOR*SQR)
                PREPY=COEFP*PNORMP*PREPY
                BASIS(INDX+1)=BASIS(INDX+1)+PREPY

                PREPZ=FACTOR**1.25D0*(RZA)*DEXP(-FACTOR*SQR)
                PREPZ=COEFP*PNORMP*PREPZ
                BASIS(INDX+2)=BASIS(INDX+2)+PREPZ

            ENDDO
         
            IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                  .OR. EXC.EQ.'RXalpha') THEN
             DO NO = 1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,2
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                      .OR. EXC.EQ.'Uxalpha') THEN
             DO NO=1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,2
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
             DO NO=1,NORB 
               NVAL =         NO
               DO INC = 0,2
               WFB(I,J,K,NO)=WFB(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOB(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ENDIF

          ELSEIF(NNTYP .EQ. 4) THEN

              PNORMD1=(2048.0D0/(9.0D0*PI**3))**0.25D0
              PNORMD2=(2048.0D0/(PI**3))**0.25D0
              PNORMD3=(2048.0D0/(16.0D0*PI**3))**0.25D0
              PNORMD4=(2048.0D0/(144.0D0*PI**3))**0.25D0

            BASIS(:) = 0
            DO N1 = 1,NPRIM(N2)

              FACTOR =  BASINF(4,N1,1) 
              COEFD  =  BASINF(4,N1,2) 

                PRED0=(2.0D0*RZA**2-RXA**2-RYA**2)*DEXP(-FACTOR*SQR)
                PRED0=FACTOR**1.75D0*COEFD*PNORMD4*PRED0
                BASIS(INDX)=BASIS(INDX)+PRED0

                PRED1=FACTOR**1.75D0*(RZA)*(RXA)*DEXP(-FACTOR*SQR)
                PRED1=COEFD*PNORMD2*PRED1
                BASIS(INDX+1)=BASIS(INDX+1)+PRED1

                PRED2=FACTOR**1.75D0*(RYA)*(RZA)*DEXP(-FACTOR*SQR)
                PRED2=COEFD*PNORMD2*PRED2
                BASIS(INDX+2)=BASIS(INDX+2)+PRED2

                PRED3=FACTOR**1.75D0*(RXA**2-RYA**2)*DEXP(-FACTOR*SQR)
                PRED3=COEFD*PNORMD3*PRED3
                BASIS(INDX+3)=BASIS(INDX+3)+PRED3

                PRED4=FACTOR**1.75D0*(RXA)*(RYA)*DEXP(-FACTOR*SQR)
                PRED4=COEFD*PNORMD2*PRED4
                BASIS(INDX+4)=BASIS(INDX+4)+PRED4

            ENDDO

            IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                  .OR. EXC.EQ.'RXalpha') THEN
             DO NO = 1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,4
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                      .OR. EXC.EQ.'Uxalpha') THEN
             DO NO=1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,4
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
             DO NO=1,NORB 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,4
               WFB(I,J,K,NO)=WFB(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOB(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ENDIF

          ELSEIF(NNTYP .EQ. 5) THEN
     
               PNORMF1=(32768.0D0/(3600.0D0*PI**3))**0.25D0
               PNORMF2=(32768.0D0/(1600.0D0*PI**3))**0.25D0
               PNORMF3=(32768.0D0/(16.0D0*PI**3))**0.25D0
               PNORMF4=(32768.0D0/PI**3)**0.25D0
               PNORMF5=(32768.0D0/(576.0D0*PI**3))**0.25D0

            BASIS(:) = 0
            DO N1 = 1,NPRIM(N2)

              FACTOR =  BASINF(5,N1,1) 
              COEFF  =  BASINF(5,N1,2) 

                 PREF0=RZA*(5.0D0*RZA**2-3.0D0*SQR)*DEXP(-FACTOR*SQR)
                 PREF0=FACTOR**2.25D0*COEFF*PNORMF1*PREF0
                 BASIS(INDX)=BASIS(INDX)+PREF0

                 PREF1=RXA*(5.0D0*RZA**2-SQR)*DEXP(-FACTOR*SQR)
                 PREF1=FACTOR**2.25D0*COEFF*PNORMF2*PREF1
                 BASIS(INDX+1)=BASIS(INDX+1)+PREF1

                 PREF2=RYA*(5.0D0*RZA**2-SQR)*DEXP(-FACTOR*SQR)
                 PREF2=FACTOR**2.25D0*COEFF*PNORMF2*PREF2
                 BASIS(INDX+2)=BASIS(INDX+2)+PREF2

                 PREF3=RZA*(RXA**2-RYA**2)*DEXP(-FACTOR*SQR)
                 PREF3=FACTOR**2.25D0*COEFF*PNORMF3*PREF3
                 BASIS(INDX+3)=BASIS(INDX+3)+PREF3

                 PREF4=RXA*RYA*RZA*DEXP(-FACTOR*SQR)
                 PREF4=FACTOR**2.25D0*COEFF*PNORMF4*PREF4
                 BASIS(INDX+4)=BASIS(INDX+4)+PREF4

                 PREF5=RXA*(RXA**2-3.0D0*RYA**2)*DEXP(-FACTOR*SQR)
                 PREF5=FACTOR**2.25D0*COEFF*PNORMF5*PREF5
                 BASIS(INDX+5)=BASIS(INDX+5)+PREF5

                 PREF6=RYA*(3.0D0*RXA**2-RYA**2)*DEXP(-FACTOR*SQR)
                 PREF6=FACTOR**2.25D0*COEFF*PNORMF5*PREF6
                 BASIS(INDX+6)=BASIS(INDX+6)+PREF6

            ENDDO

            IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
     *                  .OR. EXC.EQ.'RXalpha') THEN
             DO NO = 1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,6
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                      .OR. EXC.EQ.'Uxalpha') THEN
             DO NO=1,NORA 
C              NVAL = NCORE + NO
               NVAL =         NO
               DO INC = 0,6
               WFA(I,J,K,NO)=WFA(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOA(INDX+INC,NVAL) 
               ENDDO
             ENDDO
             DO NO=1,NORB 
               NVAL =         NO
               DO INC = 0,6
               WFB(I,J,K,NO)=WFB(I,J,K,NO)+BASIS(INDX+INC)
     *                                    *COEMOB(INDX+INC,NVAL) 
               ENDDO
             ENDDO
            ENDIF

          ENDIF

       ENDDO
       ENDDO
       ENDDO

          IF    (NNTYP .EQ. 1) THEN
            NADD = 0
          ELSEIF(NNTYP .EQ. 2) THEN
            NADD = 3
          ELSEIF(NNTYP .EQ. 3) THEN
            NADD = 2
          ELSEIF(NNTYP .EQ. 4) THEN
            NADD = 4
          ELSEIF(NNTYP .EQ. 5) THEN
            NADD = 6
          ENDIF

          INDX = INDX + NADD

      if(myid .eq. 3) then
      write(*,*) ' INDX = ', INDX
      endif

       ENDDO

       ENDDO

C     if(myid.eq.0) then
C       write(*,*) wfa(10,10,10,10)
C       write(*,*) basis(10)
C     endif

C     REWIND(24)
      DO K = 1,JZ  
      DO J = 1,JY
      DO I = 1,JX
C       WRITE(24) WFA(I,J,K,1)
      ENDDO
      ENDDO
      ENDDO

C----- Normalize for alpha spin -----

      DO 40 K = 1, NORA      

        ANORM = 0.D0

        DO 50 L = 1,JZ      
        DO 50 M = 1,JY      
        DO 50 N = 1,JX      
          ANORM = ANORM + WFA( N,M,L,K )**2*DV
50      CONTINUE

        CALL MPI_REDUCE(ANORM,BNORM,1,MPI_DOUBLE_PRECISION,
     *                  MPI_SUM,0,MPI_COMM_WORLD,IERR)

      if(myid.eq.0) then
      write(*,*) 'anorm=',BNORM
      endif
        ANORM = DSQRT(BNORM)

        CALL MPI_BCAST(ANORM,1,MPI_DOUBLE_PRECISION,
     *                 0,MPI_COMM_WORLD,IERR)

        DO 60 L = 1,JZ      
        DO 60 M = 1,JY      
        DO 60 N = 1,JX      
          WFA(  N,M,L,K ) = WFA( N,M,L,K ) / ANORM 
          WFNA( N,M,L,K ) = WFA( N,M,L,K ) 
60      CONTINUE

40    CONTINUE

      DO K = 1, NORA
        NRDA(K) = K
      ENDDO

C     CALL  DIAG(      NRDA,NORA,WFNA,WFA )
      CALL MDIAG (     NRDA,NORA,WFNA,WFA )
C     CALL MDIAG_BLK ( NRDA,NORA,WFNA,WFA )

C----- normalize for beta spin -----

      IF(EXC.EQ.'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
     *                .OR. EXC.EQ.'UXalpha') THEN

        DO 70 K = 1, NORB      

          ANORM = 0.D0

          DO 80 L = 1,JZ      
          DO 80 M = 1,JY      
          DO 80 N = 1,JX      
            ANORM = ANORM + WFB( N,M,L,K )**2*DV
80        CONTINUE

        CALL MPI_REDUCE(ANORM,BNORM,1,MPI_DOUBLE_PRECISION,
     *                  MPI_SUM,0,MPI_COMM_WORLD,IERR)

        ANORM = DSQRT(BNORM)

C       IF(MYID.EQ.0) THEN
C         write(*,*) 'anorm =',ANORM
C       ENDIF

        CALL MPI_BCAST(ANORM,1,MPI_DOUBLE_PRECISION,
     *                 0,MPI_COMM_WORLD,IERR)

          DO 90 L = 1,JZ           ! 2005.07.25 takahashi
          DO 90 M = 1,JY      
          DO 90 N = 1,JX      
            WFB(  N,M,L,K ) = WFB( N,M,L,K ) / ANORM
            WFNB( N,M,L,K ) = WFB( N,M,L,K ) 
90        CONTINUE

70      CONTINUE

        DO K = 1, NORB
          NRDB(K) = K
        ENDDO

C       CALL  DIAG(      NRDB,NORB,WFNB,WFB )
        CALL MDIAG (     NRDB,NORB,WFNB,WFB )
C       CALL MDIAG_BLK ( NRDB,NORB,WFNB,WFB )

      ENDIF

C      IF(EXC.EQ. 'RPZ' .OR. EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' 
C    *                  .OR. EXC.EQ.'RXalpha') THEN

C        DO I=1,NMAX
C        DO J=1,NMAX
C            WRITE(3,'(6E13.5)') ((WFA(I,J,K,NO),NO=2,5),K=1,NMAX)
C        ENDDO
C        ENDDO

C       ELSEIF(EXC.EQ. 'UPZ' .OR. EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' 
C    *                       .OR. EXC.EQ.'Uxalpha') THEN

C        DO I=1,NMAX
C        DO J=1,NMAX
C            WRITE(3,'(6E13.5)') ((WFA(I,J,K,NO),NO=2,5),
C    *                           (WFB(I,J,K,NO),NO=2,5),K=1,NMAX)
C        ENDDO
C        ENDDO
C      
C       ENDIF

       RETURN
       END  


C----- LONG-RANGE PART OF THE RHOIJ and VHIJ -----

      SUBROUTINE LRCLMB0

      IMPLICIT REAL*8( A-H,O-Z )
      IMPLICIT INTEGER*4( I-N )

      include 'mpif.h'
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'

      PARAMETER( PI = 3.14159265358979323D0 )
C     PARAMETER( ALPHA = 1.0D0 )
C     PARAMETER( ALPHA = 0.5D0 )
C     PARAMETER( ALPHA = 0.3D0 )
C     PARAMETER( ALPHA = 0.1D0 )
      PARAMETER( ALPHA = 0.5D-1 )
C     PARAMETER( ALPHA = 0.1D-1 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX,DY,DZ
      COMMON / ACELL / XL,YL,ZL
      COMMON / NCLR  / PNUC( NNUC,NDIM )
      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,
     *                  RSIZE(NNUC),
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF

      COMMON / PFFT/ RHOLL( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NNUC ),
     *               PHILL( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NNUC )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      RHOLL(:,:,:,:) = 0.D0
      PHILL(:,:,:,:) = 0.D0
      
      SQAL   = DSQRT(ALPHA)
      SQPI   = DSQRT(PI)
      ALPI32 = (SQAL/SQPI)**3.D0

      NNMAXPX = NMAXX/2
      NNMAXPY = NMAXY/2
      NNMAXPZ = NMAXZ/2
      NNX     = -JX*MX + NNMAXPX + 1
      NNY     = -JY*MY + NNMAXPY + 1
      NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO 10 L = 1, NF
        NNA = NPOP(L)
      DO 10 K = 1, JZ
        RZ = DZ*( K-NNZ ) - PNUC( NNA,3 )
      DO 10 J = 1, JY
        RY = DY*( J-NNY ) - PNUC( NNA,2 )
      DO 10 I = 1, JX
        RX = DX*( I-NNX ) - PNUC( NNA,1 )

C         RX = RX2 - PNUC( NNA,1 )
C         RY = RY2 - PNUC( NNA,2 )
C         RZ = RZ2 - PNUC( NNA,3 )
          R2 = RX**2 + RY**2 + RZ**2
          R  = DSQRT( R2 ) 

	  IF( R .LT. 1.0D-06 ) THEN
            VT = 2.D0*DSQRT(ALPHA)/SQPI         ! analytical function ( limit zero ) 
            PHILL(I,J,K,L) = VT
          ELSE
            PHILL(I,J,K,L) = DERF(SQAL*R)/R
          ENDIF

          RHOLL(I,J,K,L) = ALPI32*DEXP(-ALPHA*R2) 

10    CONTINUE

      RETURN
      END


C----- LONG-RANGE PART OF THE RHOIJ and VHIJ -----

      SUBROUTINE LRCLMB(RHOL,PHIL)

      IMPLICIT REAL*8( A-H,O-Z )
      IMPLICIT INTEGER*4( I-N )

      include 'mpif.h'
      include 'mpi.i'                           ! mpi
      include 'QMpara.i'

      PARAMETER( PI = 3.14159265358979323D0 )
C     PARAMETER( ALPHA = 1.0D0 )
C     PARAMETER( ALPHA = 0.5D0 )
C     PARAMETER( ALPHA = 0.3D0 )
C     PARAMETER( ALPHA = 0.1D0 )
      PARAMETER( ALPHA = 0.5D-1 )
C     PARAMETER( ALPHA = 0.1D-1 )

      COMMON / MMPI1 / MX,MY,MZ
      COMMON / MMPI4 / JX,JY,JZ
      COMMON / GRID2 / DX,DY,DZ
      COMMON / ACELL / XL,YL,ZL
      COMMON / NCLR  / PNUC( NNUC,NDIM )
      COMMON / PSN    / TAX1,TAX2,TAY1,TAY2,TAZ1,TAZ2,TAO,
     *                  TAX3,TAX4,TAY3,TAY4,TAZ3,TAZ4,PI4,
     *                  RSIZE(NNUC),
     *                  ZPOP(NFUZZY),NPOP(NFUZZY),NF

      COMMON / PFFT/ RHOLL( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NNUC ),
     *               PHILL( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NNUC )

      DIMENSION RHOL( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION PHIL( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      RHOL(:,:,:) = 0.D0
      PHIL(:,:,:) = 0.D0
      
C     SQAL   = DSQRT(ALPHA)
C     SQPI   = DSQRT(PI)
C     ALPI32 = (SQAL/SQPI)**3.D0

C     NNMAXPX = NMAXX/2
C     NNMAXPY = NMAXY/2
C     NNMAXPZ = NMAXZ/2
C     NNX     = -JX*MX + NNMAXPX + 1
C     NNY     = -JY*MY + NNMAXPY + 1
C     NNZ     = -JZ*MZ + NNMAXPZ + 1

      DO 10 L = 1, NF
        NNA = NPOP(L)
      DO 10 K = 1, JZ
C       RZ = DZ*( K-NNZ ) - PNUC( NNA,3 )
      DO 10 J = 1, JY
C       RY = DY*( J-NNY ) - PNUC( NNA,2 )
      DO 10 I = 1, JX
C       RX = DX*( I-NNX ) - PNUC( NNA,1 )

        PHIL(I,J,K) = PHIL(I,J,K) + ZPOP(L)*PHILL(I,J,K,L)    
        RHOL(I,J,K) = RHOL(I,J,K) + ZPOP(L)*RHOLL(I,J,K,L)    

C         RX = RX2 - PNUC( NNA,1 )
C         RY = RY2 - PNUC( NNA,2 )
C         RZ = RZ2 - PNUC( NNA,3 )
C         R2 = RX**2 + RY**2 + RZ**2
C         R  = DSQRT( R2 ) 

C        IF( R .LT. 1.0D-06 ) THEN
C           VT = ZPOP(L)*2.D0*DSQRT(ALPHA)          ! analytical function ( limit zero ) 
C    *         / SQPI           
C           PHIL(I,J,K) = PHIL(I,J,K) + VT
C         ELSE
C           PHIL(I,J,K) = PHIL(I,J,K) + ZPOP(L)*DERF(SQAL*R)/R
C         ENDIF

C         RHOL(I,J,K) = RHOL(I,J,K) + ZPOP(L)*ALPI32*DEXP(-ALPHA*R2) 

10    CONTINUE

      RETURN
      END

!     include 'qmmmlj.f'
!     include 'mm_obsolete.f'
      include 'ext_routines.f'


