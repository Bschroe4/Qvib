      SUBROUTINE MATREL(NKOORD,KDIM,KPOT,M1,H1,QP1,GEW1,M2,H2,QP2,GEW2,
     &                  M3,H3,QP3,GEW3,M4,H4,QP4,GEW4,M5,H5,QP5,GEW5,
     &                  JROT)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL ROTAT
      INTEGER COUNTER
      PARAMETER (NNMAX=3000)
      DIMENSION H1(50,50,3),QP1(50),GEW1(50)
      DIMENSION H2(50,50,3),QP2(50),GEW2(50)
      DIMENSION H3(50,50,3),QP3(50),GEW3(50)
      DIMENSION H4(50,50,3),QP4(50),GEW4(50)
      DIMENSION H5(50,50,3),QP5(50),GEW5(50)
      DIMENSION Q(10),PROD(1001)
      DIMENSION ABC(1001),ABD(1001),PABC(1001)
      DIMENSION PLAM(5)
      DIMENSION NK(5)
      DIMENSION RINT(10)
      COMMON/TRANS/TRA(10,10),YLAM(10)
      COMMON/SKA/RE(10),NZ,NPROP,IP(NNMAX,5),KPRO,MAX
      COMMON/POTE/NU(300,5),CV(300),NPOT,NKOART
      COMMON/MATR/XK(NNMAX,NNMAX),EVAL(NNMAX)
      COMMON/DRUCK/NDRU
      DATA EFAKT/219474.631/,Z05/0.5D0/
      GOTO (10,20,30,40,50,999,999,999),NKOORD
   10 WRITE(6,1000) M1,M1
 1000 FORMAT(/1X,'INTEGRATION POINTS PER NORMAL COORDINATE',
     11X,I3,1X,'TOTAL',1X,I6)
      if(NDRU.NE.0) WRITE(6,1003)
      DO 11 I=1,NKOORD
   11 PLAM(I)=1.D0/YLAM(I)**2
      DO 12 I1=1,M1
      Q(1)=QP1(I1)*YLAM(1)
      DO 13 I=1,NKOORD
      SU=0.
      DO 14 K=1,NKOORD
   14 SU=SU+TRA(I,K)*Q(K)
      IF(NKOART.NE.0) SU=QP1(I1)
   13 RINT(I)=SU
      AX=GEW1(I1)
      DO 15 I=1,KDIM
      N1=IP(I,1)
      ABC(I)=H1(N1,I1,1)
   15 ABD(I)=ABC(I)*AX
      VPOT=0.D0
      DO 16 I=1,KPOT
      V=CV(I)
      DO 17 J=1,NKOORD
   17 V=V*RINT(J)**NU(I,J)
   16 VPOT=VPOT+V
      IF(NDRU.NE.0) THEN
        WRITE(6,1001) (Q(I),I=1,NKOORD),(RINT(I),I=1,NKOORD),VPOT*EFAKT
      ENDIF
 1001 FORMAT(1x,2F10.5,D20.10)
 1002 FORMAT(1x,F10.5,D20.10)
 1003 FORMAT(/1x,8x,"Q1",8x,"q1",10x,"POT [CM-1]")
      DO 18 IU=1,KDIM
      Y=ABC(IU)
      N1=IP(IU,1)
      Z2=-0.5D0*PLAM(1)*H1(N1,I1,3)+Y*VPOT
      DO 18 IV=1,IU
      XK(IU,IV)=XK(IU,IV)+Z2*ABD(IV)
   18 XK(IV,IU)=XK(IU,IV)
   12 CONTINUE
      GOTO 999
   20 WRITE(6,2000) M1,M2,M1*M2
 2000 FORMAT(/1X,'INTEGRATION POINTS PER NORMAL COORDINATE',
     11X,2I3,1X,'TOTAL',1X,I6)
      DO 21 I=1,NKOORD
   21 PLAM(I)=1.D0/YLAM(I)**2
      DO 22 I1=1,M1
      Q(1)=QP1(I1)*YLAM(1)
      DO 22 I2=1,M2
      Q(2)=QP2(I2)*YLAM(2)
      DO 23 I=1,NKOORD
      SU=0.
      DO 24 K=1,NKOORD
   24 SU=SU+TRA(I,K)*Q(K)
      IF(NKOART.NE.0) SU=Q(I)/YLAM(I)
   23 RINT(I)=SU
      AX=GEW1(I1)*GEW2(I2)
      DO 25 I=1,KDIM
      N1=IP(I,1)
      N2=IP(I,2)
      ABC(I)=H1(N1,I1,1)*H2(N2,I2,1)
   25 ABD(I)=ABC(I)*AX
      VPOT=0.D0
      DO 26 I=1,KPOT
      V=CV(I)
      DO 27 J=1,NKOORD
   27 V=V*RINT(J)**NU(I,J)
   26 VPOT=VPOT+V
      IF(NDRU.NE.0) THEN
        WRITE(6,2001) (Q(I),I=1,NKOORD),(RINT(J),J=1,NKOORD),VPOT*EFAKT
      ENDIF
 2001 FORMAT(1x,4F10.5,F10.2)
      DO 28 IU=1,KDIM
      Y=ABC(IU)
      N1=IP(IU,1)
      N2=IP(IU,2)
      Z2=-Z05*(H1(N1,I1,3)*H2(N2,I2,1)*PLAM(1)
     &        +H1(N1,I1,1)*H2(N2,I2,3)*PLAM(2))
      Z2=Z2+Y*VPOT
      DO 28 IV=1,IU
      XK(IU,IV)=XK(IU,IV)+Z2*ABD(IV)
   28 XK(IV,IU)=XK(IU,IV)
   22 CONTINUE
      GOTO 999
   30 WRITE(6,3000) M1,M2,M3,M1*M2*M3
 3000 FORMAT(/1X,'INTEGRATION POINTS PER NORMAL COORDINATE',
     11X,3I3,1X,'TOTAL',1X,I6)
      DO 31 I=1,NKOORD
   31 PLAM(I)=1.D0/YLAM(I)**2
      DO 32 I1=1,M1
      Q(1)=QP1(I1)*YLAM(1)
      DO 32 I2=1,M2
      Q(2)=QP2(I2)*YLAM(2)
      DO 32 I3=1,M3
      Q(3)=QP3(I3)*YLAM(3)
      DO 33 I=1,NKOORD
      SU=0.
      DO 34 K=1,NKOORD
   34 SU=SU+TRA(I,K)*Q(K)
      IF(NKOART.NE.0) SU=Q(I)/YLAM(I)
   33 RINT(I)=SU
      AX=GEW1(I1)*GEW2(I2)*GEW3(I3)
      DO 35 I=1,KDIM
      N1=IP(I,1)
      N2=IP(I,2)
      N3=IP(I,3)
      ABC(I)=H1(N1,I1,1)*H2(N2,I2,1)*H3(N3,I3,1)
   35 ABD(I)=ABC(I)*AX
      VPOT=0.D0
      DO 36 I=1,KPOT
      V=CV(I)
      DO 37 J=1,NKOORD
   37 V=V*RINT(J)**NU(I,J)
   36 VPOT=VPOT+V
      IF(NDRU.NE.0) THEN
        WRITE(6,3001) (Q(I),I=1,NKOORD),(RINT(J),J=1,NKOORD),VPOT*EFAKT
        IF(VPOT.LT.0) WRITE(6,'(1x,"!!!WARNING!!! HOLE IN THE PEF!")')
      ENDIF
 3001 FORMAT(1x,6F10.5,F10.2)
      DO 38 IU=1,KDIM
      Y=ABC(IU)
      N1=IP(IU,1)
      N2=IP(IU,2)
      N3=IP(IU,3)
      Z2=-0.5D0*(H1(N1,I1,3)*H2(N2,I2,1)*H3(N3,I3,1)*PLAM(1)
     &          +H1(N1,I1,1)*H2(N2,I2,3)*H3(N3,I3,1)*PLAM(2)
     &          +H1(N1,I1,1)*H2(N2,I2,1)*H3(N3,I3,3)*PLAM(3))
      Z2=Z2+Y*VPOT
      DO 38 IV=1,IU
      XK(IU,IV)=XK(IU,IV)+Z2*ABD(IV)
   38 XK(IV,IU)=XK(IU,IV)
   32 CONTINUE
      GOTO 999
   40 WRITE(6,4000) M1,M2,M3,M4,M1*M2*M3*M4
 4000 FORMAT(/1X,'INTEGRATION POINTS PER NORMAL COORDINATE',
     11X,4I3,1X,'TOTAL',1X,I6)
      DO 41 I=1,NKOORD
   41 PLAM(I)=1.D0/YLAM(I)**2
      DO 42 I1=1,M1
      Q(1)=QP1(I1)*YLAM(1)
      DO 42 I2=1,M2
      Q(2)=QP2(I2)*YLAM(2)
      DO 42 I3=1,M3
      Q(3)=QP3(I3)*YLAM(3)
      DO 42 I4=1,M4
      Q(4)=QP4(I4)*YLAM(4)
      DO 43 I=1,NKOORD
      SU=0.
      DO 44 K=1,NKOORD
   44 SU=SU+TRA(I,K)*Q(K)
      IF(NKOART.NE.0) SU=Q(I)/YLAM(I)
   43 RINT(I)=SU
      AX=GEW1(I1)*GEW2(I2)*GEW3(I3)
      DO 45 I=1,KDIM
      N1=IP(I,1)
      N2=IP(I,2)
      N3=IP(I,3)
      N4=IP(I,4)
      ABC(I)=H1(N1,I1,1)*H2(N2,I2,1)*H3(N3,I3,1)*H4(N4,I4,1)
      ABD(I)=ABC(I)*AX
      DO 45 J=1,NKOORD
      DO 46 K=1,NKOORD
      NK(J)=1
      IF(K.EQ.J) NK(J)=NK(J)+2
   46 CONTINUE
      PABC(I)=PABC(I)+H1(N1,I1,NK(1))
     &               *H2(N2,I2,NK(2))
     &               *H3(N3,I3,NK(3))
     &               *H4(N4,I4,NK(4))*PLAM(J)
   45 CONTINUE
      VPOT=0.D0
      DO 47 I=1,KPOT
      V=CV(I)
      DO 48 J=1,NKOORD
   48 V=V*RINT(J)**NU(I,J)
   47 VPOT=VPOT+V
      IF(NDRU.NE.0) THEN
        WRITE(6,4001) (Q(I),I=1,NKOORD),(RINT(J),J=1,NKOORD),VPOT*EFAKT
      ENDIF
 4001 FORMAT(1x,8F10.5,F10.2)
      DO 49 IU=1,KDIM
      Y=ABC(IU)
      Z2=-Z05*PABC(IU)
      Z2=Z2+Y*VPOT
      DO 49 IV=1,IU
   49 XK(I,J)=XK(I,J)+Z2*ABD(IV)
   42 CONTINUE
      GOTO 999
   50 CONTINUE
      WRITE(6,5000) M1,M2,M3,M4,M5,M1*M2*M3*M4*M5
 5000 FORMAT(/1X,'INTEGRATION POINTS PER NORMAL COORDINATE',
     11X,4I3,1X,'TOTAL',1X,I6)
      DO 51 I=1,NKOORD
   51 PLAM(I)=1.D0/YLAM(I)**2
      DO 52 I1=1,M1
      Q(1)=QP1(I1)*YLAM(1)
      DO 52 I2=1,M2
      Q(2)=QP2(I2)*YLAM(2)
      DO 52 I3=1,M3
      Q(3)=QP3(I3)*YLAM(3)
      DO 52 I4=1,M4
      Q(4)=QP4(I4)*YLAM(4)
      DO 52 I5=1,M5
      Q(5)=QP5(I5)*YLAM(5)
      DO 53 I=1,NKOORD
      SU=0.
      DO 54 K=1,NKOORD
   54 SU=SU+TRA(I,K)*Q(K)
      IF(NKOART.NE.0) SU=Q(I)/YLAM(I)
   53 RINT(I)=SU
      AX=GEW1(I1)*GEW2(I2)*GEW3(I3)*GEW4(I4)*GEW5(I5)
      DO 55 I=1,KDIM
      N1=IP(I,1)
      N2=IP(I,2)
      N3=IP(I,3)
      N4=IP(I,4)
      N5=IP(I,5)
      ABC(I)=H1(N1,I1,1)*H2(N2,I2,1)*H3(N3,I3,1)*H4(N4,I4,1)*H5(N5,I5,1)
      ABD(I)=ABC(I)*AX
      DO 55 J=1,NKOORD
      DO 56 K=1,NKOORD
      NK(J)=1
      IF(K.EQ.J) NK(J)=NK(J)+2
   56 CONTINUE
      PABC(I)=PABC(I)+H1(N1,I1,NK(1))
     &               *H2(N2,I2,NK(2))
     &               *H3(N3,I3,NK(3))
     &               *H4(N4,I4,NK(4))
     &               *H5(N5,I5,NK(5))*PLAM(J)
   55 CONTINUE
      VPOT=0.D0
      DO 57 I=1,KPOT
      V=CV(I)
      DO 58 J=1,NKOORD
   58 V=V*RINT(J)**NU(I,J)
   57 VPOT=VPOT+V
      IF(NDRU.NE.0) THEN
        WRITE(6,5001) (Q(I),I=1,NKOORD),(RINT(J),J=1,NKOORD),VPOT*EFAKT
      ENDIF
 5001 FORMAT(1x,8F10.5,F10.2)
      DO 59 IU=1,KDIM
      Y=ABC(IU)
      Z2=-Z05*PABC(IU)
      Z2=Z2+Y*VPOT
      DO 59 IV=1,IU
   59 XK(I,J)=XK(I,J)+Z2*ABD(IV)
   52 CONTINUE
  999 CONTINUE
      RETURN
      END
      SUBROUTINE MOMENTUM(NKOORD,NN,YLAM,IP,A)
      IMPLICIT NONE
      INTEGER,PARAMETER ::NNMAX=3000
      INTEGER,PARAMETER ::I4B =SELECTED_INT_KIND(9)
      INTEGER,PARAMETER ::DP =KIND(1.0D0)
      INTEGER(I4B), INTENT(IN) :: NKOORD,NN
      INTEGER(I4B), INTENT(IN) :: IP(NNMAX,5)
      REAL(DP), INTENT(IN)     :: YLAM(NKOORD)
      REAL(DP), INTENT(INOUT)  :: A(NNMAX,NNMAX)
      INTEGER(I4B) :: I, J, K, IU, IV
      INTEGER(I4B) :: VK, VI(NKOORD,2), DELV(NKOORD)
      INTEGER(I4B) :: BRAKET(2)
      INTEGER(I4B) :: OVL,POVL
      REAL(DP) :: P, P2, ZZ
      REAL(DP) :: PLAM(NKOORD)
      REAL(DP), SAVE :: Z05=0.5D0
      DO I=1,NKOORD
        PLAM(I)=(1.D0/YLAM(I))**2
      ENDDO
      DO IU=1,NN
        DO IV=1,IU
          BRAKET(1)=IU
          BRAKET(2)=IV
          DO I=1,2
            K=BRAKET(I)
            DO J=1,NKOORD
              VI(J,I)=IP(K,J)-1
            ENDDO
          ENDDO
          DO I=1,NKOORD
            DELV(I)=VI(I,2)-VI(I,1)
          ENDDO
          P2=0.D0
          DO I=1,NKOORD
            POVL=1.D0
            IF(ABS(DELV(I)).NE.2.AND.DELV(I).NE.0) CYCLE
            DO J=1,NKOORD
              IF(J.NE.I) THEN
                OVL=0
                IF(DELV(J).EQ.0) OVL=1
                POVL=POVL*OVL
              ENDIF
            ENDDO
            IF(POVL.NE.1) CYCLE
            VK=VI(I,1)
            ZZ=PLAM(I)
            SELECT CASE (DELV(I))
              CASE (2)
                P=(VK+1)*(VK+2)
                ZZ=-ZZ*Z05*SQRT(P)
              CASE (0)
                P=(VK+Z05)
                ZZ=ZZ*P
              CASE (-2)
                P=VK*(VK-1)
                ZZ=-ZZ*Z05*SQRT(P)
            END SELECT
            P2=P2+ZZ
          ENDDO
          A(IU,IV)=A(IU,IV)+Z05*P2
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE
