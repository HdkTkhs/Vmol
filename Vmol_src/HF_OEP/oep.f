
! Copyright 2026 Hideaki Takahashi

! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C---------------------------------------------------------------
C
C     SUBROUTINE LOCAL EXCHANGE POTENTIAL PROPOSED BY J.C.Slater
C     This potential will be used as an initial guess of the OEP
C     SCF.
C
C     July,2023, H. Takahashi
C
C---------------------------------------------------------------

      SUBROUTINE OEP_LOCX( RWFA,RWFX,MOR,VXLOC )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'
      include 'QMpara.i'                     ! Vmol
!     include 'sizes.i'                      ! Vmol
      include 'BHHpara.i'                    ! Vmol

      PARAMETER ( EPS = 1.00D0 )
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
      COMMON / GRID2 / DX, DY, DZ
      COMMON / LTC   / POT                                         ! Vmol

      DIMENSION VXLOC( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )    ! Local exchange potential by J.C.Slater
!----- VXLOC ===> VEXA in sub. SCF.f -----
      DIMENSION RHO(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHOA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      DIMENSION WGT(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION VEXAL( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION RWFA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )    ! This is both for alpha and beta spins
      DIMENSION RWFX(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
 
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      CALL DNST( RWFA,RHOA,MOR )          ! spin  density for orbitals

      DO IORB = 1,MOR
        DO K = 1,JZ
        DO J = 1,JY
        DO I = 1,JX
          WGT(I,J,K,IORB) = RWFA(I,J,K,IORB)**2/RHOA(I,J,K)       
        ENDDO
        ENDDO
        ENDDO
      ENDDO

      DO IORB = 1,MOR
        DO K = 1,JZ
        DO J = 1,JY
        DO I = 1,JX
          VEXAL(I,J,K,IORB) = RWFX(I,J,K,IORB)/RWFA(I,J,K,IORB)
        ENDDO
        ENDDO
        ENDDO
      ENDDO

      DO K = 1,JZ
      DO J = 1,JY
      DO I = 1,JX
        SUM = 0.D0
        DO IORB = 1,MOR
          SUM = SUM + WGT(I,J,K,IORB)*VEXAL(I,J,K,IORB)
        ENDDO
        VXLOC(I,J,K) = SUM
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END


C------------------------------------
C
C   SUBROUTINE to optimize the wavefunctions 
C   on the local oep exchange potential
C
C   June 2023, H.Takahashi
C
C------------------------------------

      SUBROUTINE OEP_SCF( RWFA,RWFB,RWFNA,RWFNB,ORBEA,ORBEB,VNUC,VCOU_FA,VEXA,VEXB,RHO ) 
!     SUBROUTINE OEP_SCF( RWFA,RWFB,RWFNA,RWFNB,ORBEA,      VNUC,VCOU_FA,VEXA,VEXB,RHO ) 
!     SUBROUTINE OEP_SCF( RWFA,RWFB,RWFNA,RWFNB,ORBEA,VNUC,VCOU_FA,VEXA,VEXB     ) 

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
C     COMMON / OFC   / VCOU( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )     

      DIMENSION ORBEA( NORA )
      DIMENSION ORBEB( NORA )
      DIMENSION NRDA(  NORA )
      DIMENSION NRDB(  NORA )
      DIMENSION VNUC(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VCOU(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEXA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEXB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEFFA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEFFB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHO(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION TRHO(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHOA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHOB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
C     DIMENSION VHIJA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA,NORA)
C     DIMENSION VHIJB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA,NORA)
      DIMENSION VHIJA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA)
      DIMENSION VHIJB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA)
      DIMENSION VCOU_TMP( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VCOU_FA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )        ! Fermi-Amaldi potential as reference potential 20230623

      DIMENSION GNALMA(  NNUC,2,3,NORA )
      DIMENSION GNALMB(  NNUC,2,3,NORA )
      DIMENSION DGNALMA( NNUC,5,NORA )
      DIMENSION DGNALMB( NNUC,5,NORA )

      DIMENSION RWFA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ, NORA )
      DIMENSION RWFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ, NORA )
      DIMENSION RWFX(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ, NORA )
C     DIMENSION RWFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ, 1    )
      DIMENSION RWFNA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ, NORA )
      DIMENSION RWFNB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ, NORA )

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
        CALL POCH ( SCRD )
        APOT = APOT + VPCZ + PLJ                          ! add nuclear - site and LJ potential
      ELSEIF(NQMMM.EQ.'LINK') THEN                        ! Vmol
C       CALL POCH1( SCRD  )                               ! Vmol
        CALL POCH1_S( SCRD1 )                             ! Vmol
        APOT = APOT + VPCZ + PLJ                          ! Vmol
        IF( MOD(IMD,NQM) .NE. 0 ) GOTO 777
C       WRITE(*,*) 'rank',myid,'APOT=',APOT
      ELSEIF(NQMMM.EQ.'MM') THEN                          ! Vmol
C       CALL POCH2( SCRD )                                ! Vmol
        CALL POCH2( SCRD1 )                               ! Vmol
C       CALL QMF2( SCRD, FRCS,FRCN )                      ! Vmol
        CALL QMF2( SCRD1,FRCS,FRCN )                      ! Vmol
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
C     WRITE(*,*)                     
C     WRITE(*,994)                     
      ENDIF

      TOTE = 0.D0
      VHIJA( :,:,:,: ) = 0.D0       ! initialize the exchange pot. 20180309
      VHIJB( :,:,:,: ) = 0.D0       ! initialize the exchange pot. 20180309

C     IOEP = 0    ! 20230623 initialize step number of OEP SCF

10    CONTINUE

      STIME = MPI_WTIME()

      NJ   = NJ + 1
      PREE = TOTE

C----- For Restricted orbitals -----

      IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
     *                 .OR.EXC.EQ.'RXalpha') THEN     

C      compute density       

C       CALL DNST(  RWFA,RHOA,MORA )          ! spin density for orbitals
C       WRITE(*,*) 'rank ',MYID,'RHOA =',RHOA(1,1,1)
C       CALL RDNST(  RHO,RHOA,ZSUM )          ! total density for restricted orbitals
C       WRITE(*,*) 'rank ',MYID,'RHO =',RHO(1,1,1)

C      coulomb potential

C       CALL VCLMB(  RHO )
C       TIME1 = MPI_WTIME()
C       CALL RVCLMB( RHO )
C       TIME2 = MPI_WTIME()
C       IF(MYID.EQ.0) THEN
C         write(*,*) 'elapsed time:RVCLMB =',TIME2-TIME1
C       ENDIF
C   
C      exchange and correlation      

C       IF(EXC.EQ.'RXalpha') THEN
C       
C         CALL VXALP(RHO,RHOA,RHOB,VEXA,VEXB,PEXC,PEXT)
C       
C       ELSEIF(EXC.EQ.'RPZ') THEN
C       
C         PEXC = 0.D0
C         CALL PZ( RHOA,VEXA,PEXC )           ! Perdew & Zunger for spin density
C         PEXC = 2.D0*PEXC
C     
C         PEXT = 0.D0                         ! external exchange and correlation
C         DO L = 1, JZ                
C         DO M = 1, JY
C         DO N = 1, JX
C           PEXT = PEXT + RHOA(N,M,L)*VEXA(N,M,L)*DV
C         ENDDO
C         ENDDO
C         ENDDO

C         CALL MPI_REDUCE(PEXT,PEXT1,1,MPI_DOUBLE_PRECISION,
C    *                  MPI_SUM,0,MPI_COMM_WORLD,IERR)
C         PEXT = 2.D0*PEXT1

C         CALL MPI_BCAST(PEXT,1,MPI_DOUBLE_PRECISION,
C    *                 0,MPI_COMM_WORLD,IERR)

C       ELSEIF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' ) THEN
C     
C         TIME1 = MPI_WTIME()
!         CALL BLYP(RHO,RHOA,RHOB,VEXA,VEXB,PEXC,PEXT)    ! comment for pure HF calculation 20230619
C         TIME2 = MPI_WTIME()
C         IF(MYID.EQ.0) THEN
C           write(*,*) 'elapsed time:BLYP =',TIME2-TIME1
C         ENDIF

C       ENDIF

C15    CONTINUE     ! loop for OEP-SCF   20230623

C      IF ( IOEP > 0 ) THEN                  ! increment step number for OEP 20230623
       CALL DNST(  RWFA,RHOA,MORA )          ! spin density for orbitals
       CALL RDNST(  RHO,RHOA,ZSUM )          ! total density for restricted orbitals
C      ENDIF

C      estimate effective potential (local potential)

        DO L = 1, JZ
        DO M = 1, JY
        DO N = 1, JX

          VEFFA(N,M,L) = VNUC(    N,M,L )
     *                 + VCOU_FA( N,M,L )
     *                 + VEXA(    N,M,L ) 

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
!       CALL RVCLMB( RHO )

C      exchange and correlation      
    
!       IF(EXC.EQ.'UXalpha') THEN
!       
!         CALL VXALP(RHO,RHOA,RHOB,VEXA,VEXB,PEXC,PEXT)
!       
!       ELSEIF(EXC.EQ.'UPZ') THEN
!       
!       PEXC = 0.D0
!       CALL PZ( RHOA,VEXA,PEXC ) ! Perdew & Zunger for spin density
!       CALL PZ( RHOB,VEXB,PEXC ) ! Perdew & Zunger for spin density
!     
!       PEXT = 0.D0                                       ! external exchange and correlation
!       DO L = 1, JZ                
!       DO M = 1, JY
!       DO N = 1, JX
!         PEXT = PEXT + ( RHOA(N,M,L)*VEXA(N,M,L)
!    *         +          RHOB(N,M,L)*VEXB(N,M,L) )*DV
!       ENDDO
!       ENDDO
!       ENDDO
!       CALL MPI_REDUCE(PEXT,PEXT1,1,MPI_DOUBLE_PRECISION,
!    *                  MPI_SUM,0,MPI_COMM_WORLD,IERR)
!       PEXT = PEXT1

!       CALL MPI_BCAST(PEXT,1,MPI_DOUBLE_PRECISION,
!    *                 0,MPI_COMM_WORLD,IERR)

!       ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN
!     
!         CALL BLYP(RHO,RHOA,RHOB,VEXA,VEXB,PEXC,PEXT)

!       ENDIF

C      estimate effective potential (local potential)

        DO L = 1, JZ
        DO M = 1, JY
        DO N = 1, JX

          VEFFA(N,M,L) = VNUC(N,M,L)
     *                 + VCOU_FA( N,M,L )
     *                 + VEXA(N,M,L) 
          VEFFB(N,M,L) = VNUC(N,M,L)
     *                 + VCOU_FA( N,M,L )
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

!     CALL HF_EXCH_PSN( RWFA,RWFX,NORA,VHIJA,HFEXA,NJ )     ! MPI parallel
!     HFEXA = RHEX*HFEXA
!     RWFNA(:,:,:,:) = RWFNA(:,:,:,:) + RHEX*RWFX(:,:,:,:)     ! kinetic + hf-exchange

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
C----- estimate effective potential (non-local potential) -----
C----- operate hamiltonian --------
C----- orbital energy -----
C----- total energy -----

!     HFEXB = 0.D0
!     IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' ) THEN
!       CALL HF_EXCH_PSN( RWFB,RWFX,NORB,VHIJB,HFEXB )          ! MPI parallel
!     ENDIF
!     HFEXB = RHEX*HFEXB

      CALL KNTC(   RWFB,NORB,RWFNB )
!     RWFNB(:,:,:,:) = RWFNB(:,:,:,:) + RHEX*RWFX(:,:,:,:)     ! kinetic + hf-exchange
      CALL GNALM(  RWFB,GNALMB,DGNALMB,NORB )
      CALL OPERATE(RWFB,RWFNB,NORB,VEFFB,GNALMB,DGNALMB)
      CALL OENGY(  RWFB,RWFNB,ORBEB,NORB ) 
!     HFEX = HFEXA + HFEXB
!     CALL TENGY1( ORBEA,ORBEB,HFEX,TOTE )            ! for unrestricted orbitals

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
C       WRITE(*,*) ' SCF   Iteration = ', NJ   
C       WRITE(*,*) ' Total Energy    = ', PENG
        IF(NJ.NE.1) THEN
          DE = TOTE - PREE
C         WRITE(*,*) ' Delta E         = ', DE
C         WRITE(*,*) ' Residual(rho)   = ', DERHO
        ENDIF
      ENDIF

      ENDIF

      ETIME = MPI_WTIME()
      IF(MYID.EQ.0) THEN
C       write(*,*) ' elapsed time    =', ETIME-STIME
      ENDIF

      GOTO 10    ! SCF for the reference density

100   CONTINUE

C     CALL DNSHOMO( RWFA )           ! 2005.12.02 takahashi

      CALL KNTC(   RWFA,NORA,RWFNA )

C----- Hartree-Fock exchange potential -----

!     CALL HF_EXCH_PSN( RWFA,RWFX,NORA,VHIJA,HFEXA )     ! MPI parallel
!     HFEXA = RHEX*HFEXA
!     RWFNA(:,:,:,:) = RWFNA(:,:,:,:) + RHEX*RWFX(:,:,:,:)     ! kinetic + hf-exchange

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

      IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
     *                 .OR.EXC.EQ.'RXalpha') THEN     

C       CALL TENGY0( ORBEA,TOTE )              ! for restricted orbitals
        CALL TENGY0( ORBEA,HFEXA,TOTE )        ! for restricted orbitals

      ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
     *                     .OR.EXC.EQ.'UXalpha') THEN     

        CALL KNTC( RWFB,NORB,RWFNB )

C----- Hartree-Fock exchange potential -----

!     CALL HF_EXCH_PSN( RWFB,RWFX,NORB,VHIJB,HFEXB )     ! MPI parallel
!     HFEXB = RHEX*HFEXB
!     RWFNB(:,:,:,:) = RWFNB(:,:,:,:) + RHEX*RWFX(:,:,:,:)     ! kinetic + hf-exchange

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

      ENDIF

      PENG = TOTE + APOT
      DE = TOTE - PREE

C     CALL RVCLMB_OEP( RHO,VCOU_TMP )        

C     NEL = 2.D0*MORA
C     REWIND(24)
C     DO K = 1,JZ  
C     DO J = 1,JY
C     DO I = 1,JX
C       WRITE(24) VEXA(I,J,K) - VCOU_TMP(I,J,K) + DBLE(NEL)/DBLE(NEL-1)*VCOU(I,J,K)
C     ENDDO
C     ENDDO
C     ENDDO

      IF(MYID.EQ.0) THEN

C     WRITE(*,*) ' SCF   Iteration = ', NJ   
C     WRITE(*,*) ' Total Energy    = ', PENG
C     WRITE(*,*) ' Delta E         = ', DE
C     WRITE(*,*) ' Residual(rho)   = ', DERHO
C     WRITE(*,*)                     
C     WRITE(*,996)

      KAZU = ' A'
      KAZU2= ' Alpha  ' 
      
C     WRITE(*,*)                     
C     WRITE(*,997)

      CALL OROUT1(ORBEA,NORA,MORA,KAZU,KAZU2)

      IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
     *                 .OR.EXC.EQ.'UXalpha') THEN     

C     WRITE(*,999)                     

      KAZU = ' B'
      KAZU2= ' Beta   ' 
      
      CALL OROUT1(ORBEB,NORB,MORB,KAZU,KAZU2)

C     WRITE(*,999)                     
C     WRITE(*,*)                     
      
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
      WRITE(*,*) '  SCF   Iteration = ', NJ   
      WRITE(*,*) '  Final Energy    = ', PENG
C     WRITE(30,*)  PENG                    

      ENDIF                                         ! IF(MYID.EQ.0) THEN

C----- output dipole moment -----

C     CALL DPL( RHO,IMD )                           ! for QM subsystem

      CALL FLUSH(06)                                ! for SX-ACE
      CALL FLUSH(30)                                ! for SX-ACE

C----- compute energy distribution function -----

C     CALL EDF( RHO,SCRD1,PENG,IMD )                ! takahashi 2005.09.13

C----- output electronic density -----

C     IF(MYID.EQ.0) THEN
C     CALL WDNS( IMD,RHO )
C     ENDIF

C----------------------------------

C     IF(FREEZE.EQ.'FALSE') THEN

C     CALL LFRC( RHO ) ! Force (local part) for non-periodic system

C     IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
C    *                 .OR.EXC.EQ.'RXalpha') THEN     
C       BF = 2.D0
C       CALL NLFRC( RWFA,GNALMA,DGNALMA,MORA,BF ) ! Force (non-local part)
C     ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
C    *                 .OR.EXC.EQ.'UXalpha') THEN     
C       BF = 1.D0
C       CALL NLFRC( RWFA,GNALMA,DGNALMA,MORA,BF ) ! Force (non-local part)
C       CALL NLFRC( RWFB,GNALMB,DGNALMB,MORB,BF ) ! Force (non-local part)
C     ENDIF

C     ELSEIF(NQMMM.EQ.'LINK') THEN

C     CALL LFRC1( RHO ) ! Force (local part) for non-periodic system

C     IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
C    *                 .OR.EXC.EQ.'RXalpha') THEN     
C       BF = 2.D0
C       CALL NLFRC1( RWFA,GNALMA,MORA,BF ) ! Force (non-local part)
C       CALL NLFRC2( RWFA,GNALMA,MORA,BF,0 ) ! Force (non-local part)
C     ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
C    *                 .OR.EXC.EQ.'UXalpha') THEN     
C       BF = 1.D0
C       CALL NLFRC1( RWFA,GNALMA,MORA,BF ) ! Force (non-local part)
C       CALL NLFRC1( RWFB,GNALMB,MORB,BF ) ! Force (non-local part)
C       CALL NLFRC2( RWFA,GNALMA,MORA,BF,0 ) ! Force (non-local part)
C       CALL NLFRC2( RWFB,GNALMB,MORB,BF,1 ) ! Force (non-local part)
C     ENDIF

C     ENDIF

C     CALL OPTFC
C     CALL OPTFC( IMD )
C     CALL OPTFC3

C     CALL AVDNST(RHO, IMD)
C     CALL AVDNST(RWFA,IMD)

C     CALL FLUSH(24)                                ! for SX-ACE

C     CALL REDDNS(RHO,IMD)

777   CONTINUE

C     CALL FUZZYC                             ! Vmol
C     CALL EFTC(RHOA,RHOB)                          ! Vmol
C     WRITE(*,*) 'PCLMB =',PCLMB              ! Vmol
C----------------------------------------------------------------------
C     IF(NQMMM.EQ.'QMMM') THEN

C       CALL QMF( SCRD,FRCS,FRCN,RHO )

C       DO J=1,NDIM
C       DO I=1,NNUC
C         FRC(I,J) = FRC(I,J) + FRCN(I,J)
C       ENDDO
C       ENDDO

C     ELSEIF(NQMMM.EQ.'LINK') THEN            ! Vmol
C                                             ! Vmol
C       CALL SSINT1( RHO,IMD )                ! Vmol  2013.08.05 takahashi
C       CALL QMF1(   SCRD, FRCS,FRCN,RHO )    ! Vmol
C       CALL QMF1(   SCRD1,FRCS,FRCN,RHO )    ! Vmol
C       CALL RQMF1_S(SCRD1,FRCS,FRCN     )    ! Vmol
C                                             ! Vmol
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

C       DO J=1,NDIM                           ! Vmol
C       DO I=1,NNUC                           ! Vmol
C         FRC(I,J) = FRC(I,J) + FRCN(I,J)     ! Vmol
C       ENDDO                                 ! Vmol
C       ENDDO                                 ! Vmol
C                                             ! Vmol
C       IF(NLINK.GT.0) THEN                   ! Vmol
C         CALL EXTCORR( FRCS,FRC,SCRD )       ! Vmol
C         write(*,*) 'spltfrc:start'
C         CALL SPLTFRC( FRCS,FRC,SCRD )       ! Vmol
C         CALL SPLTFRC( FRCS,FRC,SCRD1 )      ! Vmol
C         write(*,*) 'spltfrc:end'
C       ENDIF                                 ! Vmol
C                                             ! Vmol
C     ENDIF

      IF(MYID.EQ.0) THEN
      WRITE(*,*)                     
C     WRITE(*,995)                     
      ENDIF

197   FORMAT( 1X,A22,3F12.6 )      
199   FORMAT( 3X,A6,3F12.6 )      

994   FORMAT( '!!!!!!!!!!!!!!!!!!!!!!
     *!!!!!! Self Consistent Field for OEP
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



C---------------------------------------------------------------
C
C     SUBROUTINE OPTIMIZED EFFECTIVE POTENTIAL BASED ON
C     THE Yang-Wu APPROACH (PRL,DOI: 10.1103/PhysRevLett.89.143002).
C
C     Note: This routine can be used both in rhf and uhf calculations. 
C
C     June,2023, H. Takahashi
C
C---------------------------------------------------------------

      SUBROUTINE OEP_YW_UHF( RWFA,RWFB,ORBELCA,ORBELCB,VNUC,VHIJA,VHIJB,VEXA,VEXB,DELV,IOEP,RHO )

      IMPLICIT REAL*8 ( A-H,O-Z )      
      IMPLICIT INTEGER*4 ( I-N )      

      include "mpif.h"
      include 'mpi.i'
      include 'QMpara.i'                     ! Vmol
!     include 'sizes.i'                      ! Vmol
      include 'BHHpara.i'                    ! Vmol

      PARAMETER ( EPS = 5.0D0 )

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
      COMMON / LTC   / POT                                         ! Vmol
!     COMMON / OFC   / VCOU( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )          ! should not be in the common block
      COMMON / OFC   / VCOU( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )          ! 20230702

      DIMENSION ORBELCA( NORA )     ! This is for the local potential Veff
      DIMENSION ORBELCB( NORA )     ! This is for the local potential Veff
      DIMENSION ORBEA(   NORA )     ! This is for the HF hamiltonian
      DIMENSION ORBEB(   NORA )
      DIMENSION NRDA(    NORA )
      DIMENSION NRDB(    NORA )

!     DIMENSION HIA(   1:MORA,MORA+1:NORA )       ! matrix: <\phi_i|H_eff|\phi_a>
!     DIMENSION HIA1(  1:MORA,MORA+1:NORA )       ! for MPI reduction

      DIMENSION HIA_A( 1:NORA,1:NORA )            ! matrix: <\phi_i|H_eff|\phi_a>
      DIMENSION HIA1_A(1:NORA,1:NORA )            ! for MPI reduction

      DIMENSION HIA_B( 1:NORA,1:NORA )            ! matrix: <\phi_i|H_eff|\phi_a>
      DIMENSION HIA1_B(1:NORA,1:NORA )            ! for MPI reduction

      DIMENSION VNUC(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEXA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEXB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION VEFF(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHO(   NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION TRHO(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHOA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION RHOB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ )
      DIMENSION DVEFF( NMAXX/NX,NMAXY/NY,NMAXZ/NZ )

      DIMENSION GNALMA(  NNUC,2,3,NORA )
      DIMENSION GNALMB(  NNUC,2,3,NORB )
      DIMENSION DGNALMA( NNUC,5,NORA )
      DIMENSION DGNALMB( NNUC,5,NORB )

      DIMENSION RWFA(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION RWFB(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION RWFX(  NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION RWFNA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION RWFNB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION VHIJA( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
      DIMENSION VHIJB( NMAXX/NX,NMAXY/NY,NMAXZ/NZ,NORA )
 
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,IERR)

      PEXC = 0.D0
      PEXT = 0.D0

      IOEP = IOEP + 1
      DV   = DX*DY*DZ

C---- Density calculation should be performed in the oep-scf.f routine. -----

!     ZSUM = DBLE(2*MORA)
!     CALL DNST(  RWFA,RHOA,MORA )          ! spin  density for orbitals
!     CALL RDNST( RHO, RHOA,ZSUM )          ! total density for restricted orbitals

!     CALL DPL( RHO,IMD )                   ! This is for check.            

C     CALL RVCLMB_OEP( RHO,VCOU )        
      CALL RVCLMB( RHO )        

      CALL KNTC(  RWFA,NORA,RWFNA )
      CALL GNALM( RWFA,GNALMA,DGNALMA,NORA )

      CALL HF_EXCH_PSN( RWFA,RWFX,NORA,MORA,VHIJA,HFEXA )   ! MPI parallel
      RWFNA(:,:,:,:) = RWFNA(:,:,:,:) + RWFX(:,:,:,:)       ! kinetic + hf-exchange

      VEFF(:,:,:) = VNUC(:,:,:) + VCOU(:,:,:)

      CALL OPERATE(RWFA,RWFNA,NORA,VEFF,GNALMA,DGNALMA)
      CALL OENGY(  RWFA,RWFNA,ORBEA,NORA ) 

      IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'
     *                 .OR.EXC.EQ.'UXalpha') THEN     
        CALL KNTC(  RWFB,NORB,RWFNB )
        CALL GNALM( RWFB,GNALMB,DGNALMB,NORB )
        CALL HF_EXCH_PSN( RWFB,RWFX,NORB,MORB,VHIJB,HFEXB )   ! MPI parallel
        RWFNB(:,:,:,:) = RWFNB(:,:,:,:) + RWFX(:,:,:,:)       ! kinetic + hf-exchange
        CALL OPERATE(RWFB,RWFNB,NORB,VEFF,GNALMB,DGNALMB)
        CALL OENGY(  RWFB,RWFNB,ORBEB,NORB ) 
      ENDIF

      IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
     *                 .OR.EXC.EQ.'RXalpha') THEN

        CALL TENGY0( ORBEA,HFEXA,TOTE )        ! for restricted orbitals

      ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'    ! UBLYP ==> UHF
     *                     .OR.EXC.EQ.'UXalpha') THEN

        HFEX = HFEXA + HFEXB
        CALL TENGY1( ORBEA,ORBEB,HFEX,TOTE ) ! for unrestricted orbitals

      ENDIF

      OLD_ENE = OEP_ENE
      OEP_ENE = TOTE + POT
      DEL_OEP = OEP_ENE - OLD_ENE

      IF ( MYID == 0 ) THEN 
        WRITE(*,*) 
        WRITE(*,80) ' Number of OEP Cycles    = ', IOEP    
        WRITE(*,90) ' Total OEP Energy (a.u.) = ', OEP_ENE
        WRITE(*,90) ' Delta E          (a.u.) = ', DEL_OEP
      ENDIF

C---- build HIA matrix for alpha spins -----

      HIA1_A(:,:) = 0.D0
      DO IVIR = MORA+1,NORA
      DO IOCC = 1,MORA

        VAL = 0.D0
        DO K = 1,JZ
        DO J = 1,JY
        DO I = 1,JX
          VAL = VAL + RWFA(I,J,K,IOCC)*RWFNA(I,J,K,IVIR)
        ENDDO
        ENDDO
        ENDDO
        HIA1_A(IOCC,IVIR) = VAL*DV

      ENDDO
      ENDDO

      NUM_DAT = NORA*NORA

      CALL MPI_REDUCE( HIA1_A(1,1),HIA_A(1,1),NUM_DAT,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

      CALL MPI_BCAST( HIA_A(1,1),NUM_DAT,MPI_DOUBLE_PRECISION,0,
     *                               MPI_COMM_WORLD,IERR )


C---- derivative of EHF with respect to rsg for alpha spins -----

      DO K = 1,JZ
      DO J = 1,JY
      DO I = 1,JX
 
        SUM = 0.D0
        DO IVIR = MORA+1,NORA
        DO IOCC = 1,MORA
          EDIFF  = ORBELCA(IOCC) - ORBELCA(IVIR)    ! Caution! ORBELCA is not equal to ORBEA!
          REVEDF = 1.D0/EDIFF 
          WF1 = RWFA(I,J,K,IVIR)
          WF2 = RWFA(I,J,K,IOCC)
          SUM = SUM + HIA_A(IOCC,IVIR)*WF1*WF2*REVEDF
        ENDDO
        ENDDO

        DVEFF(I,J,K) = SUM*DV

      ENDDO
      ENDDO
      ENDDO
  
      VEXA(:,:,:) = VEXA(:,:,:) - EPS*DVEFF(:,:,:)    ! initial guess of vexa is given in sub. oep_locx.f

      DELV1 = 0.D0      
      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        DELV1 = DELV1 + DVEFF(I,J,K)**2*DV
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE(DELV1,DELV,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      
      CALL MPI_BCAST(DELV,1,MPI_DOUBLE_PRECISION,
     *               0,MPI_COMM_WORLD,IERR)

      DELV = DSQRT(DELV)
      DELVA = DELV

C------------------------------------------------------------------------------

      IF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'             ! for beta spin
     *                 .OR.EXC.EQ.'UXalpha') THEN     

C---- build HIA matrix for beta spins -----

      HIA1_B(:,:) = 0.D0
      DO IVIR = MORB+1,NORB
      DO IOCC = 1,MORB

        VAL = 0.D0
        DO K = 1,JZ
        DO J = 1,JY
        DO I = 1,JX
          VAL = VAL + RWFB(I,J,K,IOCC)*RWFNB(I,J,K,IVIR)
        ENDDO
        ENDDO
        ENDDO
        HIA1_B(IOCC,IVIR) = VAL*DV

      ENDDO
      ENDDO

      NUM_DAT = NORA*NORA

      CALL MPI_REDUCE( HIA1_B(1,1),HIA_B(1,1),NUM_DAT,MPI_DOUBLE_PRECISION,
     *                 MPI_SUM,0,MPI_COMM_WORLD,IERR )

      CALL MPI_BCAST( HIA_B(1,1),NUM_DAT,MPI_DOUBLE_PRECISION,0,
     *                               MPI_COMM_WORLD,IERR )

C---- derivative of EHF with respect to rsg for beta spins -----

      DO K = 1,JZ
      DO J = 1,JY
      DO I = 1,JX
 
        SUM = 0.D0
        DO IVIR = MORB+1,NORB
        DO IOCC = 1,MORB
          EDIFF  = ORBELCB(IOCC) - ORBELCB(IVIR)    ! Caution! ORBELCA is not equal to ORBEA!
          REVEDF = 1.D0/EDIFF 
          WF1 = RWFB(I,J,K,IVIR)
          WF2 = RWFB(I,J,K,IOCC)
          SUM = SUM + HIA_B(IOCC,IVIR)*WF1*WF2*REVEDF
        ENDDO
        ENDDO

        DVEFF(I,J,K) = SUM*DV

      ENDDO
      ENDDO
      ENDDO
  
      VEXB(:,:,:) = VEXB(:,:,:) - EPS*DVEFF(:,:,:)    ! initial guess of vexb is given in sub. oep_locx.f

      DELV1 = 0.D0      
      DO K = 1, JZ
      DO J = 1, JY
      DO I = 1, JX
        DELV1 = DELV1 + DVEFF(I,J,K)**2*DV
      ENDDO
      ENDDO
      ENDDO

      CALL MPI_REDUCE(DELV1,DELV,1,MPI_DOUBLE_PRECISION,
     *                MPI_SUM,0,MPI_COMM_WORLD,IERR)
      
      CALL MPI_BCAST(DELV,1,MPI_DOUBLE_PRECISION,
     *               0,MPI_COMM_WORLD,IERR)

      DELV  = DSQRT(DELV)
      DELVB = DELV

      ENDIF     ! for beta spins

C------------------------------------------------------------------------------

!     REWIND(24)
      DO K = 1,JZ  
      DO J = 1,JY
      DO I = 1,JX
!       WRITE(24) VEXA(I,J,K)
      ENDDO
      ENDDO
      ENDDO

      IF ( MYID == 0 ) THEN 
        IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
     *                   .OR.EXC.EQ.'RXalpha') THEN
          WRITE(*,*) ' Norm of delta Veff = ', DELVA
        ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'             ! for beta spin
     *                 .OR.EXC.EQ.'UXalpha') THEN     
          WRITE(*,*) ' Norm of delta Veff for alpha spin = ', DELVA
          WRITE(*,*) ' Norm of delta Veff for beta  spin = ', DELVB
        ENDIF
      ENDIF

      IF(EXC.EQ.'RBLYP' .OR. EXC.EQ.'RHF' .OR.EXC.EQ.'RPZ'
     *                 .OR.EXC.EQ.'RXalpha') THEN
        DELV = DELVA
      ELSEIF(EXC.EQ.'UBLYP' .OR. EXC.EQ.'UHF' .OR.EXC.EQ.'UPZ'        ! UBLYP ==> UHF
     *                     .OR.EXC.EQ.'UXalpha') THEN
        DELV = DELVA*DBLE(MORA) + DELVB*DBLE(MORB)
        DELV = DELV/DBLE(MORA*MORB)
      ENDIF

80    FORMAT(X,A27,I4)
90    FORMAT(X,A27,F14.7)

      RETURN
      END


