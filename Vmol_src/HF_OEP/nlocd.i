C     parameters for non-local d pseudopotentials

C     NNUC:    number of atoms in QM region
C     NPD:     NPD = 0: with no non-local d
C              NPD = 1: with    non-local d

C     IMPLICIT REAL*8 ( A-H,O-Z )      
C     IMPLICIT INTEGER*4 ( I-N )      

C     integer*4 NPD

C     real*8 DVD,VNLD
C     real*8 WNLOCDXY,WNLOCDYZ,WNLOCDZX
C     real*8 WNLOCDZ2,WNLOCDX2

C     parameter ( NNUC =    10 )
      parameter ( NSL  =  9000 )
      PARAMETER ( RCUT = 3.5D0 )

      COMMON / NLOCD / NPD(NNUC),DVD(100),VNLD(100,421),
     *                 WNLOCDXY(NNUC,NSL),WNLOCDYZ(NNUC,NSL),
     *                 WNLOCDZX(NNUC,NSL),WNLOCDZ2(NNUC,NSL),
     *                 WNLOCDX2(NNUC,NSL)



