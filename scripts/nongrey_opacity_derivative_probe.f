C     Exercise native OPACTR with an analytic, temperature-dependent
C     molecular weight and opacity. Check the prescribed 1% secant
C     derivative at fixed gas pressure and restoration of the state.
      PROGRAM PROBE
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      COMMON/GRDPRA/GRD(MDEPTH),PRA(MDEPTH),PGS0(MDEPTH),
     *              ANTP(MDEPTH)
      REAL*8 TBASE(3),RBASE(3),MBASE(3),PBASE(3)
      ND=3
      NFREQ=2
      LTE=.TRUE.
      IOPTAB=-1
      IDISK=0
      IFPRAD=0
      IFMOL=1
      FREQ(1)=1.D14
      FREQ(2)=2.D14
      BNUE(1)=1.D-5
      BNUE(2)=2.D-5
      DO ID=1,ND
         TEMP(ID)=1000.D0+1000.D0*ID
         PGS(ID)=1.D5*ID
         AN=PGS(ID)/BOLK/TEMP(ID)
         CALL ELDENS(ID,TEMP(ID),AN,ANE,ENRG,ENTT,WM,1)
         TBASE(ID)=TEMP(ID)
         RBASE(ID)=DENS(ID)
         MBASE(ID)=WMM(ID)
         PBASE(ID)=PGS(ID)
      END DO
      DO IJ=1,NFREQ
         CALL OPACTR(IJ)
         DO ID=1,ND
            AK=IJ*(TBASE(ID)/3000.D0)**2*
     *            SQRT(RBASE(ID)/1.D-6)
C           At fixed P, rho scales as T**(-1.2), so kappa scales
C           as T**1.4 for the analytic EOS/opacity pair below.
            EXPECT=AK*(1.01D0**1.4D0-1.D0)/(.01D0*TBASE(ID))
            WRITE(*,100) IJ,ID,DABT1(ID)/EXPECT-1.D0,
     *       TEMP(ID)/TBASE(ID)-1.D0,DENS(ID)/RBASE(ID)-1.D0,
     *       WMM(ID)/MBASE(ID)-1.D0,PGS(ID)/PBASE(ID)-1.D0
         END DO
      END DO
  100 FORMAT(2I3,5ES24.15)
      END

      SUBROUTINE ELDENS(ID,T,AN,ANE,ENRG,ENTT,WM,IPRI)
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      ANE=1.D-7*AN
      WMM(ID)=2.3D0*1.66053906660D-24*(3000.D0/T)**.2D0
      DENS(ID)=WMM(ID)*(AN-ANE)
      ELEC(ID)=ANE
      ENRG=0.D0
      ENTT=0.D0
      WM=WMM(ID)
      END

      SUBROUTINE OPACF1(IJ)
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      DO ID=1,ND
         AK=IJ*(TEMP(ID)/3000.D0)**2*SQRT(DENS(ID)/1.D-6)
         ABSO1(ID)=DENS(ID)*AK
         SCAT1(ID)=0.D0
         EMIS1(ID)=ABSO1(ID)
      END DO
      END

      SUBROUTINE TDPINI
      INCLUDE 'IMPLIC.FOR'
      INCLUDE 'BASICS.FOR'
      INCLUDE 'MODELQ.FOR'
      INCLUDE 'ALIPAR.FOR'
      DO ID=1,ND
         HKT1(ID)=HK/TEMP(ID)
      END DO
      END

      SUBROUTINE WNSTOR(ID)
      END
      SUBROUTINE STEQEQ(ID,POP,MODE)
      END
      SUBROUTINE OPAINI(MODE)
      END
      SUBROUTINE SABOLF(ID)
      STOP 'UNUSED NON-LTE BRANCH'
      END
      SUBROUTINE RATMAL(ID,AES,BES)
      STOP 'UNUSED NON-LTE BRANCH'
      END
      SUBROUTINE LEVSOL(AES,BES,POPLTE,IIFOR,NLEVEL,MODE)
      STOP 'UNUSED NON-LTE BRANCH'
      END
      SUBROUTINE PGSET(MODE)
      STOP 'UNUSED DISK BRANCH'
      END
