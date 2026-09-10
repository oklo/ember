      PROGRAM PROBE
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      INCLUDE 'ARRAY1.FOR'
      COMMON/RYBMTX/RA(MDEPTH),RB(MDEPTH),RC(MDEPTH),VR(MDEPTH),
     *              UA(MDEPTH),UB(MDEPTH),UC(MDEPTH),
     *              VA(MDEPTH),VB(MDEPTH),VC(MDEPTH),WR(MDEPTH),
     *              WM(MDEPTH,MDEPTH)
      COMMON/DERIDT/DERT
      REAL*8 BASE(3),JAC(3),FD(3)
      ND=3
      GRAV=1.D5
      TEFF=2800.D0
      HMIX0=1.9D0
      ICONV=1
      IDCONZ=2
      ICBEGP=2
      ICENTR=0
      IOPTAB=-1
      ACONML=0.125D0
      BCONML=0.5D0
      CCONML=16.D0
      DERT=1.D-5
      BASE(1)=2800.D0
      BASE(2)=3000.D0
      BASE(3)=3350.D0
      DO IC=1,3
         SCALE=10.D0**(IC+3)
         PTOTAL(1)=SCALE
         PTOTAL(2)=1.25D0*SCALE
         PTOTAL(3)=1.75D0*SCALE
         DO ID=1,ND
            TEMP(ID)=BASE(ID)
            DENS(ID)=PTOTAL(ID)*2.3D0*1.66053906660D-24/
     *               (1.38054D-16*TEMP(ID))
            DM(ID)=PTOTAL(ID)/GRAV
            ABROSD(ID)=GRAV/HMIX0/PTOTAL(ID)
            TAURS(ID)=1.D0
            REINT(ID)=0.D0
            REDIF(ID)=0.D0
         END DO
         REINT(2)=1.D0
         CALL EVALUE(V0)
         DO J=1,3
            JAC(J)=WM(2,J)
         END DO
         DO J=1,3
            STEP=1.D-6*BASE(J)
            TEMP(J)=BASE(J)+STEP
            CALL EVALUE(VP)
            TEMP(J)=BASE(J)-STEP
            CALL EVALUE(VM)
            TEMP(J)=BASE(J)
            FD(J)=-(VP-VM)/(2.D0*STEP)
            WRITE(*,100) IC,J,JAC(J),FD(J),JAC(J)/FD(J)-1.D0
         END DO
      END DO
  100 FORMAT(2I3,3ES24.15)
      END

      SUBROUTINE EVALUE(V)
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      COMMON/RYBMTX/RA(MDEPTH),RB(MDEPTH),RC(MDEPTH),VR(MDEPTH),
     *              UA(MDEPTH),UB(MDEPTH),UC(MDEPTH),
     *              VA(MDEPTH),VB(MDEPTH),VC(MDEPTH),WR(MDEPTH),
     *              WM(MDEPTH,MDEPTH)
      DO I=1,ND
         WR(I)=0.D0
         DO J=1,ND
            WM(I,J)=0.D0
         END DO
      END DO
      CALL RYBENE
      V=WR(2)
      END

C     Analytic ideal molecular gas isolates convection matrix assembly.
C     Pressure, opacity and the divergence density prefactor are held
C     fixed; CONVEC still differentiates its thermodynamic face state.
      SUBROUTINE TRMDER(ID,T,PG,PRAD,TAU,CP,DLR,AD,RHO)
      IMPLICIT REAL*8(A-H,O-Z)
      CP=3.5D0*1.38054D-16/(2.3D0*1.66053906660D-24)
      AD=2.D0/7.D0
      DLR=-1.D0
      RHO=PG*2.3D0*1.66053906660D-24/(1.38054D-16*T)
      END

      SUBROUTINE TRMDRT(ID,T,P,CP,DLR,AD,RHO)
      IMPLICIT REAL*8(A-H,O-Z)
      STOP 'UNUSED TABLE-EOS BRANCH'
      END
