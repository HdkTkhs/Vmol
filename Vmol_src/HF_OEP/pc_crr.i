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

      COMMON / PS_PCC / NPC(NNUC),
     *                  PCDNSA(NMAXX/NX,NMAXY/NY,NMAXZ/NZ),
     *                  PCDNSB(NMAXX/NX,NMAXY/NY,NMAXZ/NZ)



