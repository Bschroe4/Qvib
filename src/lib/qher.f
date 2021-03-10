      SUBROUTINE HERPOL(MM1,NDIM,H,XQ,XJ,NQ)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NNMAX=3000)
      DIMENSION H(50,50,3),XQ(50),XJ(50),X(3000),C(3000,3000)
      COMMON/DRUCK/NDRU
      COMMON/TRANS/TRA(10,10),YLAM(10)
      COMMON/SKA/RE(10),NZ,NPROP,IP(NNMAX,5),KPRO,MAX
      COMMON/FAKUL/AZ(3000)
      DATA EFAKT/219474.631/
      PIQ=1.33133536380038D0                                            0000127
      Z05=0.5D0
      AZ(1)=1.0D0                                                       0000128
      DO 13 I=1,2999                                                    0000128
      AZ(I+1)=AZ(I)*I
   13 CONTINUE                                                          0000132
      CALL HERMIT(MM1,XQ,XJ)
      DO 4 I=1,NDIM
      DO 4 J=1,NDIM
   4  C(I,J)=0
      C(1,1)=1.D0
      C(2,1)=0
      C(2,2)=2.D0
      DO 3 I=3,NDIM
      IM1=I-1
      IM2=I-2
      DO 3 J=1,NDIM
      HMAT=0
      JM1=J-1
      IF(JM1.LT.1)GOTO 3
      HMAT=2.D0*C(IM1,JM1)
    3 C(I,J)=HMAT-2*IM2*C(IM2,J)
      X(1)=1.D0
      DO 10 I=1,MM1                                                     0000151
      DO 1 L=2,NDIM
      LL=L-1                                                            0000153
    1 X(L)=XQ(I)*X(LL)                                                  0000154
      DO 10 J=1,NDIM
      HMAT=0
      DO 5 K=1,NDIM
   5  HMAT=HMAT+C(J,K)*X(K)
   10 H(J,I,1)=HMAT
      DO 15 J=1,MM1
      H(1,J,2)=0
      H(2,J,2)=2
      H(1,J,3)=0
      H(2,J,3)=0
   15 CONTINUE
      DO 14 I=3,NDIM
      I1=I-1
      I2=I-2
      DO 14 J=1,MM1
      H(I,J,2)=2*I1*H(I1,J,1)
      H(I,J,3)=4*I1*I2*H(I2,J,1)
   14 CONTINUE
      DO 11 I=1,MM1
      DO 11 J=1,NDIM
      H(J,I,3)=H(J,I,3)-2*XQ(I)*H(J,I,2)+(XQ(I)*XQ(I)-1)
     +         *H(J,I,1)
   11 H(J,I,2)=H(J,I,2)-XQ(I)*H(J,I,1)
      DO 12 I=1,NDIM
      I1=I-1                                                            0000206
      A1=SQRT(2.D0**I1)
      A1=PIQ*SQRT(AZ(I))*A1
      DO 12 J=1,MM1                                                     0000208
      DO 12 K=1,3
   12 H(I,J,K)=H(I,J,K)/A1
      IF(NDRU.EQ.0) RETURN
      WRITE(6,200)                                      
  200 FORMAT(//,1X,'INTEGRATION POINTS AND WEIGHTS',/)                  0000149
      WRITE(6,101)(XQ(J),XJ(J),J=1,MM1)
  101 FORMAT (F20.15,D25.12)                                            0000148
      IF(NDRU.LT.2) RETURN
      WRITE(6,123) NDIM, MM1
  123 FORMAT(/1X,I3,' x',I3,' H(1)-MATRIX'/)
      DO 21 I=1,NDIM
      SU=0.D0
      DO J=1,MM1
        SU=SU+H(I,J,1)*H(I,J,1)*XJ(J)
      ENDDO
      WRITE(6,102) (H(I,J,1),J=1,MM1)   
   21 WRITE(6,103) SU      
      DO I=1,NDIM
      ENDDO
      WRITE(6,124) NDIM, MM1
  124 FORMAT(/1X,I3,' x',I3,' H(2)-MATRIX'/)
      DO 22 I=1,NDIM
      SU=0.D0
      DO J=1,MM1
        SU=SU+H(I,J,2)*H(I,J,2)*XJ(J)
      ENDDO
      WRITE(6,102) (H(I,J,2),J=1,MM1)   
   22 WRITE(6,103) SU      
      WRITE(6,125) NDIM, MM1
  125 FORMAT(/1X,I3,' x',I3,' H(3)-MATRIX'/)
      PLAM=(1.D0/YLAM(NQ))**2
      WRITE(6,'(1x,F10.3)') PLAM*EFAKT
      DO 23 I=1,NDIM
      SU=0.D0
      DO J=1,MM1
        SU=SU-Z05*H(I,J,1)*PLAM*H(I,J,3)*XJ(J)
      ENDDO
      WRITE(6,102) (H(I,J,3),J=1,MM1)   
   23 WRITE(6,103) SU*EFAKT,0.5*((I-1)+0.5)*PLAM*EFAKT      
  102 FORMAT(1X,10D14.6)
  103 FORMAT(1X,2F10.2)
   99 RETURN                                                            0000211
      END                                                               0000212
      SUBROUTINE EBASIS(M1,M2,M3,NMAX,NNE)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NNMAX=3000)
      COMMON/SKA/RE(10),NZ,NPROP,IP(NNMAX,5),KPRO,MAX
      DIMENSION OM(3)
      READ(9,1000) (OM(I),I=1,3),EMAX
 1000 FORMAT(4F10.5)
      MM1=M1-1
      MM2=M2-1
      MM3=M3-1
      IZ=0
      DO 1 I=1,MM1
      DO 1 J=1,MM2
      DO 1 K=1,MM3
      IJK=I+J+K-3
      IM1=I-1
      JM1=J-1
      KM1=K-1
      TH=0.5D0*(OM(1)*(2*I-1)+OM(2)*(2*J-1)+OM(3)*(2*K-1))
      IF(TH.GT.EMAX)GOTO 1
      IZ=IZ+1
      IP(IZ,1)=IM1
      IP(IZ,2)=JM1
      IP(IZ,3)=KM1
    1 CONTINUE
      NNE=IZ
      IF(IZ.LE.NMAX)RETURN
      WRITE(6,2000) 
 2000 FORMAT(/1X,'ZUVIELE BASISFUNKTIONEN UEBER ENERGIESELEKTION'/)
      STOP
      END
      SUBROUTINE BASIS(NKOORD,MAXIM,MAXDIM,NNA,NSYM)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NNMAX=3000)
      COMMON/SKA/RE(10),NZ,NPROP,IP(NNMAX,5),KPRO,MAX
      DIMENSION MAXIM(10)
      MAXP1=MAXIM(1)+1
      MAXP2=MAXIM(2)+1
      MAXP3=MAXIM(3)+1
      MAXP4=MAXIM(4)+1
      MAXP5=MAXIM(5)+1
      IF(NSYM.NE.0) WRITE(6,1000) NSYM
      IF(NNA.NE.9999) WRITE(6,1001) NNA
 1000 FORMAT(1X,'SYMMETRY FACTOR:  ',I10)
 1001 FORMAT(1X,'SUM RESTRICTION:  ',I10)
      GOTO (10,20,30,40,50,999,999,999,999,999),NKOORD
   10 DO 11 I=1,MAXP1
      IP(I,1)=I-1
   11 CONTINUE
      IZ=MAXIM(1)
      GOTO 999
   20 IZ=0
      DO 21 I=1,MAXP1
      DO 21 J=1,MAXP2
      IJ=I+J-NKOORD
      IJP=0
      IF(I-1.EQ.0) IJP=IJP+1
      IF(J-1.EQ.0) IJP=IJP+1
      IF(IJ.GT.NNA.AND.IJP.LT.(NKOORD-1)) GOTO 21
      IZ=IZ+1
      IF(IZ.GT.MAXDIM) CALL ERROR(1,IZ,MAXDIM)
      IP(IZ,1)=I-1
      IP(IZ,2)=J-1
   21 CONTINUE
      GOTO 999
   30 IZ=0
      DO 31 I=1,MAXP1
      DO 31 J=1,MAXP2
      DO 31 K=1,MAXP3
      IJK=I+J+K-NKOORD
      IJKP=0
      IF(I-1.EQ.0) IJKP=IJKP+1
      IF(J-1.EQ.0) IJKP=IJKP+1
      IF(K-1.EQ.0) IJKP=IJKP+1
      IF(IJK.GT.NNA.AND.IJKP.LT.(NKOORD-1)) GOTO 31
      JH2=J/2
      IF((NSYM.EQ.-1).AND.(JH2*2.EQ.J)) GOTO 31
      IF((NSYM.EQ.1).AND.(JH2*2.NE.J)) GOTO 31
      IZ=IZ+1
      IF(IZ.GT.MAXDIM) CALL ERROR(1,IZ,MAXDIM)
      IP(IZ,1)=I-1
      IP(IZ,2)=J-1
      IP(IZ,3)=K-1
   31 CONTINUE
      GOTO 999
   40 IZ=0
      DO 41 I=1,MAXP1
      DO 41 J=1,MAXP2
      DO 41 K=1,MAXP3
      DO 41 L=1,MAXP4
      IJKL=I+J+K+L-NKOORD
      IJKLP=0
      IF(I-1.EQ.0) IJKLP=IJKLP+1
      IF(J-1.EQ.0) IJKLP=IJKLP+1
      IF(K-1.EQ.0) IJKLP=IJKLP+1
      IF(L-1.EQ.0) IJKLP=IJKLP+1
      IF(IJKL.GT.NNA.AND.IJKLP.LT.(NKOORD-1)) GOTO 41
      IZ=IZ+1
      IP(IZ,1)=I-1
      IP(IZ,2)=J-1
      IP(IZ,3)=K-1
      IP(IZ,4)=L-1
   41 CONTINUE
      IF(IZ.GT.MAXDIM) CALL ERROR(1,IZ,MAXDIM)
      GOTO 999
   50 IZ=0
      DO 51 I=1,MAXP1
      DO 51 J=1,MAXP3
      DO 51 K=1,MAXP4
      DO 51 L=1,MAXP5
      IL=I+L
      IL2=IL/2
      IF((NSYM.EQ.-1).AND.(IL2*2.NE.IL)) GOTO 51
      IF((NSYM.EQ. 1).AND.(IL2*2.EQ.IL)) GOTO 51
      DO 52 M=1,MAXP1
      IJKLM=I+J+K+L+M-NKOORD
      IJKLMP=0
      IF(I-1.EQ.0) IJKLMP=IJKLMP+1
      IF(J-1.EQ.0) IJKLMP=IJKLMP+1
      IF(K-1.EQ.0) IJKLMP=IJKLMP+1
      IF(L-1.EQ.0) IJKLMP=IJKLMP+1
      IF(M-1.EQ.0) IJKLMP=IJKLMP+1
      IF(IJKLM.GT.NNA.AND.IJKLMP.LT.(NKOORD-1)) GOTO 51
      IZ=IZ+1
      IF(IZ.GT.MAXDIM) CALL ERROR(1,IZ,MAXDIM)
      IP(IZ,1)=I-1
      IP(IZ,2)=J-1
      IP(IZ,3)=K-1
      IP(IZ,4)=L-1
      IP(IZ,5)=M-1
   52 CONTINUE
   51 CONTINUE
  999 NNA=IZ
      RETURN
      END
      SUBROUTINE ERROR(N,NVAL,MAX)
      IMPLICIT REAL*8(A-H,O-Z)
      GOTO (10,20),N
   10 WRITE(6,11), NVAL,MAX
   11 FORMAT(//1X,'!!!!!ZUVIELE OSZILLATORPRODUKTFUNKTIONEN:',I4,
     &            ' MAX:',I4/)
      STOP
   20 WRITE(6,21), VAL,MAX
   21 FORMAT(//1X,'!!!!!ZAHL DER INTEGRATIONSPUNKTE ZU GROSS:',I4,
     &            'MAX:',I4/)
      END
      SUBROUTINE HERMIT(NN,X,A)                                          100
      IMPLICIT REAL*8(A-H,O-Z)                                               11
C     NN: GRAD DES HERMITE-POLYNOMS                                          12
C     X:  NULLSTELLEN DES HERMITE-POLYNOMS                                   13
C     A:  GEWICHTSFAKTOREN                                                   14
      DIMENSION X(NN),A(NN),FAK(9)
      COMMON/FAKUL/GAMMA(3000)
      DATA FAK/1.7724538509D0,.16667D0,1.85575D0,1.14D0,0.426D0,             16
     11.86D0,0.86D0,1.91D0,0.91D0/                                           17
      EPS=1.D-20
      Z2=2.D0                                                                18
      FN=NN                                                                  19
      N1=NN-1                                                                20
      N2=(NN+1)/2                                                            21
      CC=FAK(1)*GAMMA(NN)/(Z2**N1)
c      CC=CC*GAMMA(NN)
      S=(2*NN+1)**FAK(2)                                                     24
      DO 10 I=1,N2                                                           25
      IF(I-1) 10,1,2                                                         26
C     GROESSTE NULLSTELLE                                                    27
    1 XT=S**3-FAK(3)/S                                                       28
      GOTO 9                                                                 29
    2 IF(I-2) 10,3,4                                                         30
C     ZWEITE NULLSTELLE                                                      31
    3 XT=XT-FAK(4)*FN**FAK(5)/XT                                             32
      GOTO 9                                                                 33
    4 IF(I-3) 10,5,6                                                         34
C     DRITTE NULLSTELLE                                                      35
    5 XT=FAK(6)*XT-FAK(7)*X(1)                                               36
      GOTO 9                                                                 37
    6 IF(I-4) 10,7,8                                                         38
C     VIERTE NULLSTELLE                                                      39
    7 XT=FAK(8)*XT-FAK(9)*X(2)                                               40
      GOTO 9                                                                 41
C     ALLE UEBRIGEN NULLSTELLEN                                              42
    8 XT=Z2*XT-X(I-2)                                                        43
    9 CALL HROOT(XT,NN,DPN,PN1,EPS)                                          44
      X(I)=XT                                                                45
      A(I)=CC/DPN/PN1                                                        46
      NI=NN-I+1                                                              47
      X(NI)=-XT                                                              47
   10 A(NI)=A(I)                                                             47
      RETURN                                                                 49
      END                                                                    50
      SUBROUTINE HROOT(X,NN,DPN,PN1,EPS)                                     51
      IMPLICIT REAL*8(A-H,O-Z)                                               52
C     X: NAEHERUNGSWERT FUER NULLSTELLE                                      53
C     DPN: 1. ABLEITUNG VON H(N) BEI X                                       54
C     PN1: WERT VON H(N-1) BEI X                                             55
      ITER=0                                                                 56
    1 ITER=ITER+1                                                            57
      CALL HRECUR(P,DP,PN1,X,NN)                                             58
      D=P/DP                                                                 59
      X=X-D
      IF(DABS(D)-EPS) 3,3,2
    2 IF(ITER-10) 1,3,3                                                      61
    3 DPN=DP                                                                 62
      RETURN                                                                 63
      END                                                                    64
      SUBROUTINE HRECUR(PN,DPN,PN1,X,NN)                                     65
      IMPLICIT REAL*8(A-H,O-Z)                                               66
      DATA Z1,Z2/1.D0,2.D0/                                                  67
      P1=Z1                                                                  68
      P=X                                                                    69
      DP1=0.D0                                                               70
      DP=Z1                                                                  71
      DO 1 J=2,NN                                                            72
      FJ=J                                                                   73
      FJ2=(FJ-Z1)/Z2                                                         74
      Q=X*P-FJ2*P1                                                           75
      DQ=X*DP+P-FJ2*DP1                                                      76
      P1=P                                                                   77
      P=Q                                                                    78
      DP1=DP                                                                 79
    1 DP=DQ                                                                  80
      PN=P                                                                   81
      DPN=DP                                                                 82
      PN1=P1                                                                 83
      RETURN                                                                 84
      END                                                                    85
