      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NPU=5000,NC=500,NKO=10)
      DIMENSION S(NPU,NKO),E(NPU),ECEPA(NPU),W(NPU),C(NC),CN(NC)
      DIMENSION F(NPU,NC),TEXT(10),DAT(8),SCALE(NKO),NDEZ(NPU)
      DIMENSION COORD(NKO),EF(NPU),CREF(NKO),BET(3)
      DIMENSION CSTO2(NC),IPOT2(NC,NKO),IPOT(NC,NKO)
      DIMENSION SCFAK(NKO),REFKOAU(NKO),REFKO(NKO)
      COMMON/COEF/CSTO1(NC),IPTOT(NC,NKO)
      CHARACTER(LEN=90) IOString, IOString2
      CHARACTER(LEN=30) CART
      W12=1/DSQRT(2.D0)
      PI=DACOS(-1.D0)
      NQ=NPU
      READ(5,1000) TEXT
 1000 FORMAT(10A8)
      WRITE(6,1001) TEXT
 1001 FORMAT(5X,5H*****,10A8,5H*****//)
      READ(5,1002) DAT
 1002 FORMAT(8F10.8)
 1022 FORMAT(6F10.5,f20.12,2f10.5)
      NKOORD=DAT(1)
      IF(NKOORD.GT.NKO) STOP'NKOORD TOO LARGE!'
C     NUMBER OF COORDINATES
      NFU=(DAT(1)-NKOORD)*10+0.1
      WRITE(IOString,'(A1,I1,A19,A1)') '(',NKOORD,
     &'F10.5,F18.12,f10.2',')'
      SELECT CASE (NFU)
        CASE (3)
          WRITE(CART,'("NORMAL",22x)')
        CASE (7)
          WRITE(CART,'("INTERNAL WITH TORSION",8x)')
        CASE DEFAULT
          WRITE(CART,'(I3,27x)') NFU
      END SELECT
C     NUMBER OF POINTS
      NP=DAT(2)
      NTYP=(DAT(2)-NP)*10+0.1
C     NUMBER OF PARAMETERS
      NF=DAT(3)
      NFCEPA=(DAT(3)-NF)*1000
C     SHIFT IN WEIGHTING FACTORS
      NW=DAT(4)
C     SCALING FACTORS FOR COORDINATE VALUES
      NSCALE=DAT(5)
      NADD=(DAT(5)-NSCALE)*10+0.1
C     SCF REFERENCE ENERGY
!     EREF=DAT(6)
C     CEPA REFERENCE ENERGY
!     EREFCEP=DAT(7)
C     SCALING OF ENERGY VALUES
      ESCALE=DAT(8)
C     WRITE INPUT DATA...
      WRITE(6,'(1x,"INPUT")')
      WRITE(6,'(1x,"COORDINATES:  ",I3)') NKOORD
      WRITE(6,'(1x,"TYPE:           ",A30)') CART
      WRITE(6,'(1x,"POINTS:       ",I3)') NP
      WRITE(6,'(1x,"COEFFICIENTS: ",I3)') NF
C     SHIFT FOR REFERENZKOORDINATES
!     READ(5,1002) (CREF(I),I=1,NKOORD)
!     WRITE(6,'(/1x,"SHIFT FOR REFERENZKOORDINATES")')
!     WRITE(6,'(1x,8F9.5)') (CREF(I),I=1,NKOORD)
C     REFERENCE COORDINATES           
!     READ(5,1002) (REFKO(I),I=1,NKOORD)
!     WRITE(6,'(/1x,"REFERENCE COORDINATES")')
!     WRITE(6,'(1x,8F9.5)') (REFKO(I),I=1,NKOORD)
C     IF(MOD(NFU,2).EQ.0) READ(5,1002) BET
C     IF(MOD(NFU,2).EQ.0) WRITE(6,1010) BET
c1010 FORMAT(/'FIT MIT MORSE-POLYNOM(EN), BETAS=',6F10.3/)
      IF(NSCALE.NE.0) THEN
        READ(5,1002) (SCALE(I),I=1,NKOORD)
        WRITE(6,'(/1x,"SCALINGFACTORS FOR COORDINATE VALUES")')
        WRITE(6,'(1x,8F12.8)') (SCALE(I),I=1,NKOORD)
      ENDIF
      DO 89 I=1,NKOORD
      REFKOAU(I)=REFKO(I)*SCALE(I)
   89 SCFAK(I)=SCALE(I)
      if(nfu.ne.3.and.nfu.ne.9.and.NFU.ne.7) then                   
        read(5,1002) re1,re2
        WRITE(6,1002) re1,re2
      endif
      if(nfu.eq.4.or.nfu.eq.5) then                  
        read(5,1002) beta1,beta2
        WRITE(6,1002) beta1,beta2
      endif
      DO 11 I=1,NF
   11   READ(5,1005) (IPOT(I,J),J=1,NKOORD)
 1005 FORMAT(10I3)
      DO 1 I=1,NP
      READ(5,IOString) (S(I,J),J=1,NKOORD),E(I),W(I)
!     WRITE(6,IOString) (S(I,J),J=1,NKOORD),E(I),W(I)
      IF(NFU.EQ.1) THEN               
      S(I,1)=1-RE1/(S(I,1)+CREF(1))
      S(I,2)=1-RE2/(S(I,2)+CREF(2))
      S(I,3)=1-DCOS(S(I,3)*PI/180.D0)
      ENDIF
      IF(NFU.EQ.2) THEN               
      S(I,1)=S(I,1)+CREF(1)-re1
      S(I,2)=S(I,2)+CREF(2)-re2    
      S(I,3)=1.0d0*(1-DCOS(S(I,3)*PI/180.D0))
      ENDIF
C     READ(5,1002) (S(I,J),J=1,NKOORD),E(I),ECEPA(I),W(I)
C     IF(NFU.EQ.2) THEN
C     IBET=NKOORD
C     IF(BET(3).EQ.0) IBET=2
C     DO 802 J=1,MIN(NKOORD,IBET)
C 802 S(I,J)=(1-DEXP(-(S(I,J)-CREF(J))*BET(J)/CREF(J)))/BET(J)
C     ENDIF
      IF(NFU.EQ.3) THEN
      DO 33 J=1,NKOORD
   33 S(I,J)=S(I,J)-CREF(J)
      ENDIF
      IF(NFU.EQ.4) THEN
      r1=s(i,1)+cref(1)        
      r2=s(i,2)+cref(2)
      s(i,1)=(1-dexp(-beta1*(r1/re1-1)))/beta1
      s(i,2)=(1-dexp(-beta2*(r2/re2-1)))/beta2
      S(I,3)=1.0d0*(1-DCOS(S(I,3)*PI/180.D0))
      ENDIF
      IF(NFU.EQ.5) THEN
      r1=s(i,1)+cref(1)        
      r2=s(i,2)+cref(2)
      s(i,1)=(1-dexp(-beta1*(r1/re1-1)))/beta1
      s(i,2)=(1-dexp(-beta2*(r2/re2-1)))/beta2
      S(I,3)=s(i,3)*PI/180.d0                    
      ENDIF
      IF(NFU.EQ.6) THEN
      DO 803 J=1,2
  803 S(I,J)=(1-DEXP(-(S(I,J)-CREF(J))*BET(J)/CREF(J)))/BET(J)
      S(I,1)=W12*(S(I,1)+S(I,2))
      S(I,2)=S(I,1)-2*W12*S(I,2)
      S(I,3)=S(I,3)-CREF(3)
      ENDIF
      IF(NSCALE.NE.0) THEN
      DO 2 J=1,NKOORD
    2 S(I,J)=S(I,J)*SCALE(J)
      ENDIF
      IF(NFU .EQ. 8) THEN
      S(I,1)=1-DEXP(-(S(I,1)-CREF(1))*BET(1))
      S(I,2)=1-DEXP(-(S(I,2)-CREF(2))*BET(2))
      ENDIF
      IF(NFU.EQ.9) THEN
      DO 34 J=1,3
   34 S(I,J)=S(I,J)-CREF(J)*scale(J)
      S(I,4)=SIN(S(I,4))
      S(I,5)=SIN(S(I,5))
      ENDIF
      DO 3 K=1,NKOORD
    3 COORD(K)=S(I,K)
C     IF(NTYP.EQ.1) E(I)=ECEPA(I)
      E(I)=E(I)-EREF
      IF(NW.EQ.0) W(I)=1
C      IF(SHIFT.NE.0) W(I)=1/(E(I)+SHIFT)**2
      IF(NFU.EQ.7.OR.NFU.EQ.9) THEN
        CALL FITFU6T(IPOT,I,NKOORD,COORD,F,NF,NQ)
      ELSE
        CALL FITFU(IPOT,I,NKOORD,COORD,F,NF,NQ)
      ENDIF
    1 CONTINUE
C.....Counting first NDEZ zero decimal places
      NCOUNT=0
      DO 122 J=1,16
      IF(J.eq.1) then
      NDEZ(J)=E(1)*10**(J-1)
      ELSE
      NPREV=0
      DO K=1,J 
      NPREV=NPREV+NDEZ(K)*10**(J-K)
      ENDDO
      NDEZ(J)=E(1)*10**(J-1)-NPREV
      ENDIF
      IF(NDEZ(J).eq.0) THEN
      IF(NCOUNT.eq.0) NZERO=J-1
      NCOUNT=NCOUNT+1
      IF(NCOUNT.ge.3.and.NVAL.eq.1) GOTO 123
      ELSE
      NVAL=1
      NCOUNT=0
      NZERO=0
      ENDIF
      IF(J.eq.16.and.NZERO.eq.0) NZERO=10
  122 continue
C............................................
  123 NZ=0
      CALL SFIT(E,F,W,C,EF,NP,NF,NQ,SIG)
    8 IF(ESCALE.NE.0) THEN
      SIGN=SIG*ESCALE
      DO 6 I=1,NF
      IF(NZ.EQ.0) CSTO1(I)=C(I)
    6 CN(I)=C(I)*ESCALE
      ENDIF
      IF(NZ.EQ.0) WRITE(6,1003) SIG,SIGN
 1003 FORMAT(/1X,'STANDARD DEVIATION [A.U.] AND [CM-1]:',5X,2(D14.7)/)
 1008 FORMAT(1X,I4,5X,2(D13.7,5X),10I4)
      DO 21 I=1,NF
      IF(NZ.EQ.0)WRITE(6,1008) I,C(I),CN(I),(IPOT(I,J),J=1,NKOORD)
   21 IF(NZ.EQ.1)WRITE(6,1008) I,C(I),CN(I),(IPOT2(I,J),J=1,NKOORD)
      IF(NZ.EQ.0) WRITE(6,1020)
 1020 FORMAT(/1X,'COMPARISON'/)
      IF(NZ.EQ.0) THEN
      IF(NZERO.lt.10.and.NZERO+4.lt.10) THEN
      WRITE(IOString2,'(A1,A6,I1,A8,I1,A1,I1,A5,A1)') '(','1X,I4,',
     & NKOORD,'F10.6,3F',NZERO+4,'.',NZERO,',F6.3',')'
      ENDIF
      IF(NZERO.lt.10.and.NZERO+4.gt.10) THEN
      WRITE(IOString2,'(A1,A6,I1,A8,I2,A1,I1,A5,A1)') '(','1X,I4,',
     & NKOORD,'F10.6,3F',NZERO+4,'.',NZERO,',F6.3',')'
      ENDIF
      IF(NZERO.ge.10.and.NZERO+4.ge.10) THEN
      WRITE(IOString2,'(A1,A6,I1,A8,I2,A1,I2,A5,A1)') '(','1X,I4,',
     & NKOORD,'F10.6,3F',NZERO+4,'.',NZERO,',F6.3',')'
      ENDIF
      DO 4 I=1,NP
      DE=EF(I)-E(I)
    4 WRITE(6,IOString2) I,(S(I,J)/SCALE(J),J=1,NKOORD),E(I),EF(I),DE,
     &W(I)
      ENDIF
      IF(NZ.GT.0) GOTO 99
      IZ=0
      SIG=0
      NFSCF=NF
      NF=NFCEPA
      DO 12 I=1,NF
   12 READ(5,1005) (IPOT2(I,J),J=1,NKOORD)
      DO 7 I=1,NP
      DO 9 K=1,NKOORD
    9 COORD(K)=S(I,K)
      IF(ECEPA(I).EQ.0) THEN
      W(I)=0
      ELSE
      IZ=IZ+1
      ECEPA(I)=ECEPA(I)-EREFCEP
      IF(NTYP.EQ.0) E(I)=ECEPA(I)
      IF(NTYP.EQ.1) E(I)=E(I)+ECEPA(I)
      ENDIF
    7 CALL FITFU(IPOT2,I,NKOORD,COORD,F,NF,NQ)
      CALL SFIT(E,F,W,C,EF,NP,NF,NQ,SIG)
      DO 16 I=1,NF
   16 CSTO2(I)=C(I)
      SIG=SIG*DSQRT(DFLOAT(NP-1)/DFLOAT(IZ-1))
      NZ=NZ+1
      GOTO 8
 1004 FORMAT(1X,I3,6F14.10,4(F14.10,3X),F8.4)
   99 CONTINUE
      DO 22 I=1,NFSCF
      DO 22 J=1,NKOORD
   22 IPTOT(I,J)=IPOT(I,J)
      IF(NADD.EQ.0) GOTO 999
      WRITE(6,1007)
 1007 FORMAT(/1X,'KOEFFIZIENTENADDITION:'/)
      DO 17 I=1,NFSCF
      DO 18 J=1,NFCEPA
      DO 19 K=1,NKOORD
      IDIFF=IPOT2(J,K)-IPOT(I,K)
      IF(IDIFF.NE.0) GOTO 18
   19 CONTINUE
      CSTO1(I)=CSTO1(I)+CSTO2(J)
   18 CONTINUE
   17 CONTINUE
      DO 20 I=1,NFSCF
   20 WRITE(6,1006) I,CSTO1(I),(IPOT(I,J),J=1,NKOORD)
 1006 FORMAT(1X,I5,5X,D13.7,3X,10I4)
  999 CALL SHIFTE(NKOORD,NFSCF,CREF,REFKO,REFKOAU,SCFAK)
      STOP
      END
      SUBROUTINE SFIT(Y,G,W,C,YF,NP,NF,NQ,SIG)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(NP),YF(NP),G(NQ,NF),U(500,500),C(NF),W(NP),F(5000,500)
      DO 1 L=1,NP
      YF(L)=Y(L)
      DO 1 J=1,NF
    1 F(L,J)=G(L,J)
      DO 10 J=1,NF
      DO 10 K=1,J
      SO=0
      DO 11 L=1,NP
   11 SO=SO+F(L,J)**2*W(L)
      U(K,J)=0.
      S=0.
      DO 12 L=1,NP
   12 S=S+F(L,J)*F(L,K)*W(L)
      IF(K.EQ.J)GOTO 15
      DO 13 L=1,NP
   13 F(L,J)=F(L,J)-S*F(L,K)
      DO 14 L=1,K
   14 U(L,J)=U(L,J)-S*U(L,K)
      GOTO 10
   15 IF(S.GT.SO*1.D-6)S=1.D0/DSQRT(S)
      U(J,J)=1.D0
      DO 16 L=1,NP
   16 F(L,J)=F(L,J)*S
      DO 17 L=1,J
   17 U(L,J)=U(L,J)*S
   10 CONTINUE
      DO 20 J=1,NF
      UC=0.
      C(J)=0.
      DO 21 L=1,NP
   21 UC=UC+YF(L)*F(L,J)*W(L)
      DO 22 K=1,J
   22 C(K)=C(K)+UC*U(K,J)
      DO 23 L=1,NP
   23 YF(L)=YF(L)-F(L,J)*UC
   20 CONTINUE
      DO 30 L=1,NP
      SIG=SIG+YF(L)**2*W(L)
   30 YF(L)=Y(L)-YF(L)
      SIG=DSQRT(SIG/(NP-1))
      RETURN
      END
      SUBROUTINE FITFU6T(IPOT,I,NKOORD,COORD,F,NF,NQ)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NPU=5000,NC=500,NKO=10)
      DIMENSION COORD(1),F(NPU,NC),IPOT(NC,NKO),HILF(10)
      DO 1 J=1,NF
      DO 2 K=1,NKOORD
      IF((COORD(K).EQ.0).AND.(IPOT(J,K).EQ.0))THEN
        HILF(K)=1
      ELSE IF(K.EQ.6) THEN
          HILF(K)=DCOS(IPOT(J,K)*COORD(K))
C       IF(MOD(IPOT(J,4)*IPOT(J,5),2).NE.0) THEN
C         HILF(K)=DCOS(IPOT(J,K)*COORD(K))
C       ELSE
C         HILF(K)=DCOS(IPOT(J,K)*COORD(K))**2
C       ENDIF
      ELSE    
        HILF(K)=COORD(K)**IPOT(J,K)
      ENDIF
    2 CONTINUE
      PROD=HILF(1)
      DO 3 K=2,NKOORD
    3 PROD=PROD*HILF(K)
    1 F(I,J)=PROD
      RETURN
      END
      SUBROUTINE FITFU(IPOT,I,NKOORD,COORD,F,NF,NQ)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NPU=5000,NC=500,NKO=10)
      DIMENSION COORD(NKO),F(NPU,NC),IPOT(NC,NKO),HILF(NKO)
      DO 1 J=1,NF
      DO 2 K=1,NKOORD
      IF((COORD(K).EQ.0).AND.(IPOT(J,K).EQ.0))THEN
      HILF(K)=1
      ELSE
      HILF(K)=COORD(K)**IPOT(J,K)
      ENDIF
    2 CONTINUE
      PROD=HILF(1)
      DO 3 K=2,NKOORD
    3 PROD=PROD*HILF(K)
    1 F(I,J)=PROD
      RETURN
      END
      SUBROUTINE SHIFTE(NKOOR,NKOEF,SHIFT,REF,REFAU,SC)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*8 ZC,BEFEHL1,BEFEHL2
      PARAMETER (NPU=5000,NC=500,NKO=10)
      DIMENSION REF(NKO),REFAU(NKO),SC(NKO),DEL(NKO)
      DIMENSION AUMIN(NKO), SIMIN(NKO)
      CHARACTER(LEN=80) AUString
      CHARACTER(LEN=80) SIString
      CHARACTER(LEN=80) SHIFTString
      CHARACTER(LEN=80) DELString
      CHARACTER(LEN=80) IOString,IOString2
      COMMON/FAK/FA(30)
      COMMON/COEF1/DRU
      COMMON/COEF/ CO(NC),IP(NC,NKO)
      DIMENSION IND(500,2),IZ(2),IHI(10),SHIFT(10)
      DIMENSION COH(500),NN(10)
      DATA BEFEHL1/'MIN     '/
      DATA BEFEHL2/'SHIFT   '/
      DATA NWORK,NKOMAX/500,10/
      DATA FACTE/2.1947463d5/
      write(AUString,'(A15)') 'MINIMUM IN A.U.'
      write(SIString,'(A14)') 'MINIMUM IN ANG'
      write(SHIFTString,'(A13)') 'SHIFT IN A.U.'
      write(DELString,'(A12)') 'SHIFT IN ANG'
      CALL MAFAK
  601 READ(5,1000,ERR=600,END=600) DRU,EPS,H,ZIT
!     WRITE(6,'(1x,A,4F10.5)') "DRU,EPS,H,ZIT",DRU,EPS,H,ZIT
 1000 FORMAT(8F10.5)
      ITMAX=ZIT+0.1
      IF(ITMAX.EQ.0) ITMAX=80
      IF(EPS.EQ.0.D0) EPS=10
      EPS=10**(-EPS)
      IF(H.EQ.0.D0) H=5
      H=10**(-H)
!     WRITE(6,'(1x,A,F10.1,2D10.3,F10.1)') "DRU,EPS,H,ZIT",DRU,EPS,H,ZIT
      IF(NKOEF.GT.NWORK) CALL ERROR(1)
      IF(NKOOR.GT.NKOMAX) CALL ERROR(2)
      READ(5,1005,END=600) ZC,SHIFT
 1005 FORMAT(A8,2X,7F10.5)
      IF(ZC.EQ.BEFEHL1) THEN
      CALL MINI(NKOEF,SHIFT,NKOOR,EPS,H,ITMAX)
      DO 87 I=1,NKOOR
      DEL(I)=SHIFT(I)/SC(I)
      AUMIN(I)=SHIFT(I)+REFAU(I)
   87 SIMIN(I)=SHIFT(I)/SC(I)+REF(I)
      WRITE(6,2005)
      WRITE(6,2100) SHIFTString,(SHIFT(I),I=1,NKOOR)
      WRITE(6,2100) DELString,(DEL(I),I=1,NKOOR)
      write(6,2099) AUString,(AUMIN(I),I=1,NKOOR)
      write(6,2100) SIString,(SIMIN(I),I=1,NKOOR)
      write(11,2101) (SIMIN(I),I=1,NKOOR)
      ELSE IF(ZC.EQ.BEFEHL2) THEN
      DO 88 I=1,NKOOR
      DEL(I)=SHIFT(I)
      SHIFT(I)=SHIFT(I)*SC(I)
      AUMIN(I)=SHIFT(I)+REFAU(I)
   88 SIMIN(I)=SHIFT(I)/SC(I)+REF(I)
      WRITE(6,2005)
      WRITE(6,2100) SHIFTString,(SHIFT(I),I=1,NKOOR)
      WRITE(6,2100) DELString,(DEL(I),I=1,NKOOR)
      write(6,2099) AUString,(AUMIN(I),I=1,NKOOR)
      write(6,2100) SIString,(SIMIN(I),I=1,NKOOR)
      write(11,2101) (SIMIN(I),I=1,NKOOR)
      ELSE
        DO I=1,NKOEF
          COH(I)=0.D0
        ENDDO
      GOTO 50
      ENDIF
      DO 5 I=1,NKOEF
    5 COH(I)=CO(I)
      CALL IVAL(IND,IP,NWORK,NKOOR,NKOEF)
      NKOR=NKOEF
      DO 60 I=1,NKOEF
      DO 11 J=1,NKOOR
   11 NN(J)=IP(I,J)
      L1=-1
   71 L1=L1+1
      IF(L1.GT.NN(1)) GOTO 60
      IF(NKOOR.EQ.1) GOTO 61
      L2=-1
   72 L2=L2+1
      IF(L2.GT.NN(2)) GOTO 71
      IF(NKOOR.EQ.2) GOTO 62
      L3=-1
   73 L3=L3+1
      IF(L3.GT.NN(3)) GOTO 72
      IF(NKOOR.EQ.3) GOTO 63
      L4=-1
   74 L4=L4+1
      IF(L4.GT.NN(4)) GOTO 73
      IF(NKOOR.EQ.4) GOTO 64
      L5=-1
   75 L5=L5+1
      IF(L5.GT.NN(5)) GOTO 74
      IF(NKOOR.EQ.5) GOTO 65
      L6=-1
   76 L6=L6+1
      IF(L6.GT.NN(6)) GOTO 75
      IF(NKOOR.EQ.6) GOTO 66
      L7=-1
   77 L7=L7+1
      IF(L7.GT.NN(7)) GOTO 76
      IF(NKOOR.EQ.7) GOTO 67
      L8=-1
   78 L8=L8+1
      IF(L8.GT.NN(8)) GOTO 77
      IF(NKOOR.EQ.8) GOTO 68
      L9=-1
   79 L9=L9+1
      IF(L9.GT.NN(9)) GOTO 78
      IF(NKOOR.EQ.9) GOTO 69
      L10=-1
   80 L10=L10+1
      IF(L10.GT.NN(10)) GOTO 79
      IHI(10)=L10
   69 IHI(9)=L9
   68 IHI(8)=L8
   67 IHI(7)=L7
   66 IHI(6)=L6
   65 IHI(5)=L5
   64 IHI(4)=L4
   63 IHI(3)=L3
   62 IHI(2)=L2
   61 IHI(1)=L1
      PR=1.D0
      DO 17 KK=1,NKOOR
      IF(IHI(KK).EQ.0) GOTO 17
      PR=PR*BINOM(NN(KK),IHI(KK))*SHIFT(KK)**IHI(KK)
   17 IHI(KK)=IP(I,KK)-IHI(KK)
      CALL JVAL(IZ,IHI)
      IF(IZ(1).EQ.IND(I,1).AND.IZ(2).EQ.IND(I,2)) GOTO 14
      Z=COH(I)*PR
      DO 12 KK=1,NKOR
   12 IF(IZ(1).EQ.IND(KK,1).AND.IZ(2).EQ.IND(KK,2)) GOTO 13
      NKOR=NKOR+1
      IF(NKOR.GT.NWORK) CALL ERROR(3)
      KK=NKOR
      IND(KK,1)=IZ(1)
      IND(KK,2)=IZ(2)
      CO(KK)=0
   13 CO(KK)=CO(KK)+Z
   14 GOTO (71,72,73,74,75,76,77,78,79,80) NKOOR
   60 CONTINUE
      IF(NKOR.EQ.NKOEF) GOTO 50
      KK=NKOEF+1
      DO 51 I=KK,NKOR
      CALL HILF(IND(I,1),IND(I,2),IHI,NKOOR)
      DO 51 J=1,NKOOR
   51 IP(I,J)=IHI(J)
      NKOEF=NKOR
   50 DO 20 I=1,NKOEF
      SU=0
      DO J=1,NKOOR
        SU=SU+IP(I,J)
      ENDDO
      IF(ZC.NE.BEFEHL1) THEN
        WRITE(11,2006) CO(I),(IP(I,J),J=1,NKOOR)
      ELSE
        IF(SU.GE.2.d0) THEN
          WRITE(11,2000) CO(I),(IP(I,J),J=1,NKOOR)
        ENDIF
      ENDIF
   20 CONTINUE
 2006 FORMAT(F13.10,10I2)
 2000 FORMAT(F10.8,10I2)
      IF(DRU.GT.0.D0) THEN
      WRITE(6,2001)
      WRITE(IOString,'(A1,A10,I1,A12,A1)') '(','1x,F14.10,',NKOOR,
     &'I2,F14.10',')'
      WRITE(IOString2,'(A1,A11,I1,A37,A1)') '(',"/1X,'FOR',",NKOOR,
     &"(F10.6,5X),'CALCULATED VALUE:',F20.10",')'
      DO 21 I=1,NKOEF
   21 WRITE(6,IOString) CO(I),(IP(I,J),J=1,NKOOR),COH(I)
      ENDIF
      WRITE(IOString,'(A1,A10,I2,A44,A1)') "(","1X,'FOR ',",NKOOR,
     &     "(F10.6,5X),'CALCULATED VALUE:',2F20.10,F10.2",")"
      WRITE(IOString2,'(A1,I2,A12,A1)') "(",NKOOR,"F10.6,F20.10",")"
      IX=0
  500 READ(5,IOString2,END=600,ERR=600) (SHIFT(I),I=1,NKOOR),EREF
      IF(IX.EQ.0) THEN
        IX=1
        WRITE(6,'(/1x,A)') "CHECK"
      ENDIF
      Z=0.d0
      DO I=1,NKOEF
        X=CO(I)
        DO J=1,NKOOR
          X=X*SHIFT(J)**IP(I,J)
        ENDDO
        Z=Z+X
      ENDDO
      WRITE(6,IOString) (SHIFT(I),I=1,NKOOR),Z,Z-EREF,(Z-EREF)*FACTE
      GOTO 500
  600 STOP
 2005 FORMAT(/1X,'CALCULATED MINIMUM:'/)
 2099 format(/1x,A15,1x,8F14.8)
 2100 format(1x,A15,1x,8F14.8)
 2101 format(8F10.7)
 2003 FORMAT(1X,6F20.10)
 2004 FORMAT(1X,'FOR ',6(F10.6,5X),'CALCULATED VALUE:',F20.10)
 2001 FORMAT(/1X,'NEW COEFFICIENTS / OLD COEFFICIENTS'/)
 2002 FORMAT(1X,F10.7,6I3,6X,F10.7)
      END
      FUNCTION BINOM(N,K)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FAK/FA(30)
      BINOM=FA(N+1)/(FA(K+1)*FA(N-K+1))
      RETURN
      END
      SUBROUTINE MAFAK
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FAK/FA(30)
      FA(1)=1.D0
      DO 1 I=2,30
      N=I-1
    1 FA(I)=FA(N)*N
      RETURN
      END
      SUBROUTINE ERROR(I)
      GOTO (10,20,30,40) I
   10 WRITE(6,11)
   11 FORMAT(1X,'INPUT: ZU VIELE KOEFFIZIENTEN!')
      STOP 'FEHLERHALT'
   20 WRITE(6,21)
   21 FORMAT(1X,'INPUT: ZU VIELE KOORDINATEN!')
      STOP 'FEHLERHALT'
   30 WRITE(6,31)
   31 FORMAT(1X,'WORK ZU KLEIN!')
      STOP 'FEHLERHALT'
   40 WRITE(6,41)
   41 FORMAT(1X,'KEINE KONVERGENZ?')
      STOP 'FEHLERHALT'
      END
      SUBROUTINE MINI(NKOEF,START,NKOOR,EPS,H,ITMAX)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/COEF1/DRU
      DIMENSION START(NKOOR),COOR(6),F(3)
      IF(DRU.GT.1.D0) WRITE(6,2222)
 2222 FORMAT(/1X,'ITERATIONS MINIMUM SEARCH:'/)
      FST=FUWE(START,NKOOR,NKOEF)
      IZ=0
   10 IF(DRU.GT.1.D0) WRITE(6,1010) IZ,(START(I),I=1,NKOOR),FST
      DFG=0
      DO 1 I=1,NKOOR
      DO 2 J=1,NKOOR
    2 COOR(J)=START(J)
      COOR(I)=COOR(I)-H
      DO 3 J=1,3
      F(J)=FUWE(COOR,NKOOR,NKOEF)
    3 COOR(I)=COOR(I)+H
      DF=(F(3)-F(1))*H*0.5D0/(F(1)+F(3)-2*F(2))
      DFG=DFG+DABS(DF)
    1 START(I)=START(I)-DF
      FN=FUWE(START,NKOOR,NKOEF)
      IF(DFG.LT.EPS) RETURN
      IZ=IZ+1
      IF(IZ.GT.ITMAX) CALL ERROR(4)
      FST=FN
      GOTO 10
 1010 FORMAT(1X,I4,3X,7(F15.8,2X))
      END
      FUNCTION FUWE(COORD,NKOOR,NCO)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/COEF/CO(500),IP(500,10)
      COMMON/COEF1/DRU
      DIMENSION COORD(NKOOR)
      FUWE=0.D0
      DO 1 I=1,NCO
      PR=1.D0
      DO 2 J=1,NKOOR
      IF((COORD(J).EQ.0).AND.(IP(I,J).EQ.0)) THEN
      PR1=1
      ELSE
      PR1=COORD(J)**IP(I,J)
      ENDIF
    2 PR=PR*PR1
    1 FUWE=FUWE+PR*CO(I)
      RETURN
      END
      SUBROUTINE IVAL(IND,IP,NWORK,NKOOR,NCO)
      DIMENSION IND(NWORK,1),IP(NWORK,1),IHI(1),IJ(2),IK(2),INDEX(2)
      SAVE IZ,IZMAX,INDEX
      IZ=(NKOOR-1)/5+1
      INDEX(1)=5
      INDEX(IZ)=MOD(NKOOR+IZ-1,6)
      DO 1 I=1,NCO
      IND(I,1)=0
      IND(I,2)=0
      DO 1 J=1,IZ
      DO 1 K=1,INDEX(J)
    1 IND(I,J)=IND(I,J)+IP(I,(J-1)*5+K)*10**((5-K)*2)
      RETURN
      ENTRY JVAL(IJ,IHI)
      IJ(1)=0
      IJ(2)=0
      DO 2 J=1,IZ
      DO 2 K=1,INDEX(J)
    2 IJ(J)=IJ(J)+IHI((J-1)*5+K)*10**((5-K)*2)
      RETURN
      ENTRY HILF(I1,I2,IHI,NKOOR)
      IK(1)=I1
      IK(2)=I2
      DO 3 J=1,IZ
      DO 3 I=1,INDEX(J)
      IFAK=10**((5-I)*2)
      IH=IK(J)/IFAK
      IHI((J-1)*5+I)=IH
    3 IK(J)=IK(J)-IH*IFAK
      RETURN
      END
