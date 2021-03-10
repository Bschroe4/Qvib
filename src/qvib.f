      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INERT
      CHARACTER(len=20) ENDE,BEFEHL
      LOGICAL EBAS,ROTAT
      CHARACTER(len=50) IOStr1
      CHARACTER(len=50) IOStr2
      CHARACTER(len=50) IOStr3
      CHARACTER(len=80) cout
      PARAMETER (NNMAX=3000)
      DIMENSION XQ(3),AK0(3)
      DIMENSION H1(50,50,3),XQ1(50),XJ1(50)
      DIMENSION HX1(50,50,3),XQX1(50),XJX1(50)
      DIMENSION H2(50,50,3),XQ2(50),XJ2(50)
      DIMENSION HX2(50,50,3),XQX2(50),XJX2(50)
      DIMENSION H3(50,50,3),XQ3(50),XJ3(50)
      DIMENSION HX3(50,50,3),XQX3(50),XJX3(50)
      DIMENSION H4(50,50,3),XQ4(50),XJ4(50)
      DIMENSION H5(50,50,3),XQ5(50),XJ5(50)
      DIMENSION BEFEHL(20),GI(10,10),F(10,10)
      DIMENSION KFAK(10),DAT(9),FRM(10)
      DIMENSION G(10,10),EG(10,10),FR(10),XM(11),ZKOORD(11),ZZ(11)
      DIMENSION HELP(50,50),HELP2(50,50),HELP3(50,50)
      DIMENSION FCMAT(99,99)
      DIMENSION ANH(10)
      DIMENSION TRAINV(10,10)
      DIMENSION PLAM(10)
      DIMENSION NUMAX(10)
      DIMENSION NFU(10)
      DIMENSION WFU(99),MASSI(99,5)
      DIMENSION VK(5),XIJ(10,10),DSB(NNMAX)
      DIMENSION ANKOEFF(1001),IKOEFF(1001)
      DIMENSION TRACART(15,5),B(50,50)
      COMMON/TRANS/TRA(10,10),YLAM(10)
      COMMON/SKA/RE(10),NZ,NPROP,IP(NNMAX,5),KPRO,MAX
      COMMON/POTE/NU(300,5),CV(300),NPOT,NKOART
      COMMON/MATR/XK(NNMAX,NNMAX),EVAL(NNMAX)
      COMMON/WAVE/EIVEC(1000,99)
      COMMON/DRUCK/NDRU
      COMMON/PROPER/CPROP(60),IPROP(60,5),NNE
      DATA ENDE/'---'/
      DATA FA1,FA2,EFAKT,FAC/4.3597482,0.529177249,219474.631,
     *57.29577951/
      DATA FCM,FMHZ,FMAS,FKK/16.8580,505379.05,1822.88851,15.5689412/
      DATA ZM1,CVEL/0.1D0,2.99792458D4/
      ROTAT=.FALSE.
      JA=1
      JE=1
      NMAX=50
      NDMAX1=50
      FA3=FA1/FA2**2
      FA4=FA1/FA2
  100 FORMAT(10I4)
 5555 FORMAT(F10.8,10I2)
  102 FORMAT(F10.5,10I2)
  107 FORMAT(/)                                                         0000017
  111  FORMAT(1X,10E13.6)
  112  FORMAT(1X,20F10.2)
  213 FORMAT(/,1X,'POTENTIAL TERMS IN ATOMIC UNITS'/)
  214 FORMAT(//,1X,'PROPERTY-FUNCTION IN ATOMIC UNITS'/)
      WRITE(6,3002)
 3002 FORMAT(/3x,19('-'),/3x,' QVib V1.0 (2019) ',/3x,19('-'),
     1 ///3X,'******* INPUT *******'///)
      IZ=0
  199 READ(5,1000,END=198) BEFEHL
      IF(BEFEHL(1).EQ.ENDE) GOTO 198
      IZ=IZ+1
      WRITE(9,1000) BEFEHL
      WRITE(6,3001) IZ,BEFEHL
 3001 FORMAT(1X,I5,5X,20A4)
      GOTO 199
  198 REWIND 9
    1 READ(9,1000,END=99) BEFEHL
 1000 FORMAT(20A4)
      IF(BEFEHL(1).EQ.ENDE) GOTO 99
      WRITE(6,2002) BEFEHL
 2002 FORMAT(//1X,20A4//)
      READ(9,1110) DAT
      EBAS=.FALSE.
      NKOORD=DAT(1)+ZM1
      NKP1=NKOORD+1
      NV=DAT(2)+ZM1
      NPOT=(DAT(2)-NV)*10+ZM1
      IF(NPOT.EQ.0)NPOT=7
      NEIG=DAT(3)+ZM1
      IF(NEIG.EQ.0) NEIG=20
      NKOART=DAT(4)+ZM1
      NXIJ=DAT(5)+ZM1
      NDRU=DAT(6)+ZM1
      EPS=DAT(7)
      IF(EPS.EQ.0) EPS=0.05
      READ(9,1110) DAT
      NFU1=DAT(1)+ZM1
      NFU2=DAT(2)+ZM1
      NFU3=DAT(3)+ZM1
      NFU4=DAT(4)+ZM1
      NFU5=DAT(5)+ZM1
      NNA=DAT(6)+ZM1
      IF(NNA.EQ.0) NNA=9999
      READ(9,1110) DAT
      M1=DAT(1)+ZM1
      M2=DAT(2)+ZM1
      M3=DAT(3)+ZM1
      M4=DAT(4)+ZM1
      M5=DAT(5)+ZM1
      NZ=1
      NSYM=0
 1110 FORMAT(7F10.5,2F5.2)
  110 FORMAT(10F10.8)
 1001 FORMAT(10A8)
      IF(NKOART.NE.0) GOTO 1900
      READ(9,110) (RE(I),I=1,NKOORD)
 3004 FORMAT(1X,'EQUILIBRIUM BOND LENGTHS:',5F10.5)
      WRITE(6,3004) (RE(I),I=1,NKOORD)
  777 READ(9,110) (XM(I),I=1,NKP1)
 3005 FORMAT(1X,'MASSES:',10F15.8)
      WRITE(6,3005) (XM(I),I=1,NKP1)
 1900 CONTINUE
 2902 FORMAT(/1X,'STRETCH-BEND ANHARMONICITY [CM-1]',/)
      IF(NXIJ.NE.0) THEN ! ADD S-B CORRECTION TERMS TO DIAGONAL
        WRITE(6,2902)
        ! FOR EVERY BEND (NKOORD-1) READ A NEW LINE
        IE=NKOORD-1
        IF(IE.EQ.0) IE=1 ! SPECIAL TREATMENT FOR 1D CASE
        DO I=1,IE
          !READ THE NKOORD X_sb IN [CM-1]
          READ(9,1110) DAT
          ! CONVERT TO A.U. AND SAVE IN XIJ
          DO J=1,NKOORD
            XIJ(J,I)=DAT(J)/EFAKT
          ENDDO
          ! PRINT TO STDOUT
          NB=NKOORD+I
          IF(NXIJ.LT.0) NB=I
          WRITE(6,'(1x,I3,10F10.2)') NB,(XIJ(J,I)*EFAKT,J=1,NKOORD)
        ENDDO
        IF(NXIJ.EQ.2) THEN
          READ(9,'(3D25.15)') AK0
          WRITE(6,'(/1X,"A-COEFFICIENTS:",5x,3D14.6)') (AK0(I),I=1,3)
        ENDIF
      ENDIF
      IF(NKOART.NE.0) GOTO 1901
      ZKOORD(1)=0.0
      TOTMAS=XM(1)*FMAS
      DO 42 I=1,NKOORD
      ZKOORD(I+1)=ZKOORD(I)+RE(I)/FA2
   42 TOTMAS=TOTMAS+XM(I+1)*FMAS
      ZS=0.0
      DO 66 I=1,NKP1
   66 ZS=ZS+XM(I)*ZKOORD(I)*FMAS
      ZS=ZS/TOTMAS
      DO 97 I=1,NKP1
   97 ZKOORD(I)=ZKOORD(I)-ZS
      INERT=0.0
      DO 47 I=1,NKP1
   47 INERT=INERT+XM(I)*ZKOORD(I)*ZKOORD(I)*FMAS
      BECM=1.d0/(2.d0*INERT)*EFAKT
      BEMHZ=1.d0/(2.d0*INERT)*EFAKT*CVEL
      WRITE(6,3003) BECM,BEMHZ
 3003 FORMAT(/1X,'EQUILIBRIUM ROTATIONAL',
     & ' CONSTANT [CM-1] AND [MHz]',5X,F15.7,F15.2)
 1901 CONTINUE
      READ(9,'(I5)',ERR=1903) NDEZ
 1905 CONTINUE
      IF(NDEZ.GE.10) THEN
 3900 FORMAT(A1,A1,I2,A1,I2,A1,I1,A2,A1)
        WRITE(IOStr3,3900) '(','F',NDEZ+2,'.',NDEZ
     &                        ,',',NKOORD,'I2',')'
        GOTO 1904
      ELSE
 3901 FORMAT(2A1,I2,A1,I1,A1,I1,A2,A1)
        WRITE(IOStr3,3901) '(','F',NDEZ+2,'.',NDEZ
     &                        ,',',NKOORD,'I2',')'
        GOTO 1904 
      ENDIF
 1903 BACKSPACE 9
      NDEZ=8
      GOTO 1905
 1904 CONTINUE
      DO 51 I=1,NV
   51 READ(9,IOStr3) CV(I),(NU(I,K),K=1,NKOORD) 
!     FOR POTENTIAL IN NORMAL-KOORDS EXTRACT HARMONIC
!     WAVENUMBER/DIAGONAL QUADRATIC COEFFICIENTS
      IF(NKOART.NE.0) THEN
        IZ=0
        DO I=1,NV
          NSU=0 
          DO J=1,NKOORD
            NSU=NSU+NU(I,J)
          ENDDO
          IF(NSU.EQ.2) THEN
            IZ=IZ+1
            DO J=1,NKOORD
              IF(NU(I,J).EQ.2) THEN
                FR(J)=2*CV(I)*EFAKT
                YLAM(J)=(FR(J)/EFAKT)**2
              ENDIF
            ENDDO
          ENDIF
          IF(IZ.GE.NKOORD) EXIT
        ENDDO
      ENDIF
      WRITE(6,213)
      DO 61 I=1,NV
   61 WRITE(6,104) I,CV(I),(NU(I,J),J=1,NKOORD)
  104 FORMAT(1X,I5,1x,F12.10,10I3)
 2003 FORMAT(/1X,'HARMONIC VIBRATIONAL FREQUENCIES [CM-1]'/)
 2900 FORMAT(1x,I2,F10.2)
      IF(NKOART.NE.0) THEN
        WRITE(6,2003)
        DO I=1,NKOORD
          WRITE(6,2900) I,FR(I)
        ENDDO
        DO I=1,NKOORD
          YLAM(I)=1.0D0/(YLAM(I)**0.25D0)                                   000006
        ENDDO
        ZPE=0.D0
        DO I=1,NKOORD
          ZPE=ZPE+0.5D0*FR(I)
        ENDDO
 2901 FORMAT(/1X,'HARMONIC ZERO-POINT ENERGY [CM-1]'//,1x,F10.2/)
        WRITE(6,2901) ZPE
        GOTO 1902
      ENDIF
      cout='HARMONIC PART FINISHED'
      CALL KRAFTK(NKOORD,NV,F)
      WRITE(6,1123)
 1123 FORMAT(//1X,'F-MATRIX'/)
      DO 332 I=1,NKOORD
  332 WRITE(6,2001) (F(I,J),J=1,I)
      DO 43 I=1,NKP1
   43 XM(I)=1.d0/(XM(I)*FMAS)
      WRITE(6,'(/1X,"B-MATRIX"/)')
      DO 402 I=1,NKOORD
      DO 402 J=1,3*(NKOORD+1)
  402 B(I,J)=0
      IL=3
      DO 432 I=1,NKOORD
      IR=IL+3
      B(I,IL)=-1
      B(I,IR)=1
  432 IL=IL+3
      DO I=1,NKOORD
        WRITE(6,2222) (B(I,J),J=1,3*(NKOORD+1))
      ENDDO
      DO 12 I=1,NKOORD
      DO 12 J=1,NKOORD
   12 G(I,J)=0.
      DO 13 I=1,NKOORD
      DO 13 J=1,I
      IF(IABS(I-J).GT.1) GOTO 13
      IF(I.NE.J) GOTO 14
      G(I,I)=XM(I)+XM(I+1)
      GOTO 13
   14 G(I,J)=-XM(I)
      G(J,I)=G(I,J)
   13 CONTINUE
      WRITE(6,2000)
 2000 FORMAT(/1X,'G-MATRIX'/)
      DO 331 I=1,NKOORD
  331 WRITE(6,2222) (G(I,J),J=1,I)
 2001 FORMAT(1X,10D15.7)
      CALL DIAG(10,NKOORD,G,FR,EG,1)
      DO 17 I=1,NKOORD
      DO 17 J=1,I
      SU=0.D0
      SU1=0.
      DO 18 K=1,NKOORD
      SPEI=DSQRT(FR(K))
      SU1=SU1+EG(I,K)*EG(J,K)/FR(K)
   18 SU=SU+EG(I,K)*EG(J,K)*SPEI
      HELP(I,J)=SU1
      HELP(J,I)=SU1
      G(I,J)=SU
   17 G(J,I)=SU
      WRITE(6,2005)
 2005 FORMAT(/1X,'INVERSE G-MATRIX'/)
      DO 2006 I=1,NKOORD
 2006 WRITE(6,2222) (HELP(I,J),J=1,I)
      DO 62 I=1,NKP1
      DO 62 J=1,NKOORD
      HELP2(I,J)=0.D0
      IF(I.EQ.J) HELP2(I,I)=-1.D0
      IF((I-J).EQ.1) HELP2(I,J)=1.D0
   62 CONTINUE
      DO 63 I=1,NKP1
      DO 63 J=1,NKOORD
      SU=0.D0
      DO 64 K=1,NKOORD
   64 SU=SU+HELP2(I,K)*HELP(K,J)
   63 HELP3(I,J)=SU
      DO 65 I=1,NKOORD
      SU=0.D0
      DO 67 K=1,NKP1
   67 SU=SU+ZKOORD(K)*HELP3(K,I)
   65 FRM(I)=SU*2.D0/FA2
      DO 19 I=1,NKOORD
      DO 19 J=1,I
      SU=0.D0
      DO 31 K=1,NKOORD
      DO 31 L=1,NKOORD
   31 SU=SU+G(I,K)*F(K,L)*G(L,J)
   19 EG(I,J)=SU
      CALL DIAG(10,NKOORD,EG,FR,EG,1)
      DO 37 I=1,NKOORD
      DO 37 J=1,NKOORD
      TRA(I,J)=0.
      DO 38 K=1,NKOORD
   38 TRA(I,J)=TRA(I,J)+G(I,K)*EG(K,J)
   37 CONTINUE
      WRITE(6,2004)
 2004 FORMAT(/1X,'L-MATRIX'/)
      DO 333 I=1,NKOORD
      DO 334 J=1,NKOORD
  334 TRAINV(I,J)=TRA(I,J)
  333 WRITE(6,2222) (TRA(I,J),J=1,NKOORD)
      CALL MATINV(TRAINV,NKOORD,DET)
 2222 FORMAT(1X,12(F10.5,1X))
      WRITE(6,2115)
 2115 FORMAT(/1X,'INVERSE L-MATRIX'/)
      DO 2119 I=1,NKOORD
 2119 WRITE(6,2222) (TRAINV(I,J),J=1,NKOORD)
      DO I=1,NKOORD
         DO J=1,NKOORD
           SU=0.D0
           DO K=1,NKOORD
             SU=SU+TRA(I,K)*TRAINV(K,J)
           ENDDO
           IF(SU.LT.1.D-12) SU=0.D0
           HELP(I,J)=SU
         ENDDO     
      ENDDO
 2117 FORMAT(/1X,'ORTHOGONALITY CHECK'/)
      WRITE(6,2117)
      DO 2118 I=1,NKOORD
 2118 WRITE(6,2222) (HELP(I,J),J=1,NKOORD)
      WRITE(6,2003)
      DO 32 I=1,NKOORD
      YLAM(I)=FR(I)
   32 FR(I)=EFAKT*DSQRT(FR(I))
      DO I=1,NKOORD
        YLAM(I)=1.0D0/(YLAM(I)**0.25D0)
      ENDDO
      DO I=1,NKOORD
        WRITE(6,2900) I,FR(I)
      ENDDO
      ZPE=0.D0
      DO I=1,NKOORD
        ZPE=ZPE+0.5D0*FR(I)
      ENDDO
      WRITE(6,2901) ZPE
      DO I=1,NKOORD
      DO J=1,3*(NKOORD+1)
        JZ=(J+2)/3
        B(I,J)=B(I,J)*DSQRT(XM(JZ))
      ENDDO
      ENDDO
      HELP=0.d0
      TRACART=0.d0
      DO 500 I=1,3*(NKOORD+1) 
      DO 500 J=1,NKOORD
      SU=0 
      DO 505 K=1,NKOORD
  505 SU=SU+B(K,I)*TRAINV(J,K)                           
  500 HELP(I,J)=SU
      DO 510 I=1,3*NKP1
      IZ=(I+2)/3
      DO 510 J=1,NKOORD
  510 TRACART(I,J)=HELP(I,J)*DSQRT(XM(IZ))
  520 FORMAT(/1X,'TRANSFORMATION NORMAL- AND CARTESIAN KOORD.:'/)
      WRITE(6,520)
      DO I=1,3*(NKOORD+1)
        WRITE(6,'(1x,10D13.6)') (TRACART(I,J),J=1,NKOORD)
      ENDDO
      write(6,*)
!     CALL MATOUT(TRA,KARMAX,NKART,NKOORD,1)
 1902 CONTINUE
      WRITE(6,'(1x,"HARMONIC PART FINISHED")') 
      CALL TIMER(cout)
      NFU(1)=NFU1
      NFU1=NFU1+1
      NFU(2)=NFU2
      NFU2=NFU2+1
      NFU(3)=NFU3
      NFU3=NFU3+1
      NFU(4)=NFU4
      NFU4=NFU4+1
      NFU(5)=NFU5
      NFU5=NFU5+1
      WRITE(6,'(/1x,"CONSTRUCTING BASIS SET")') 
      CALL BASIS(NKOORD,NFU,NNMAX,NNA,NSYM)
      NN=NNA
      WRITE(6,211) NN
  211 FORMAT(1X,'SIZE OF BASIS SET ',I10)
      IF(NDRU.NE.0) THEN
  201 FORMAT(/1X,'BASIS SET'/)
      WRITE(6,201)
  888 DO 34 I=1,NN 
   34 WRITE(6,35) I,(IP(I,J),J=1,NKOORD)
      ENDIF
   35 FORMAT(1x,I5,10X,10I4)
      DO J=1,NKOORD
        NUMAX(J)=0
      ENDDO
      DO I=1,NV
        DO J=1,NKOORD
          IF(NU(I,J).GT.NUMAX(J)) NUMAX(J)=NU(I,J)
        ENDDO
      ENDDO
      write(cout,'("INTEGRATION FINISHED")')
      IF(M1.EQ.0) M1=0.5*(2*(NFU1-1)+NUMAX(1))+1
      IF(M2.EQ.0) M2=0.5*(2*(NFU2-1)+NUMAX(2))+1
      IF(M3.EQ.0) M3=0.5*(2*(NFU3-1)+NUMAX(3))+1
      IF(M4.EQ.0) M4=0.5*(2*(NFU4-1)+NUMAX(4))+1
      IF(M5.EQ.0) M5=0.5*(2*(NFU5-1)+NUMAX(5))+1
      DO 6 I=1,NN                                                       0000075
      DO 6 J=1,I
    6 XK(I,J)=0.0D0                                                     0000077
      CALL HERPOL(M1,NFU1,H1,XQ1,XJ1,1)
      IF(NKOORD.GE.2) CALL HERPOL(M2,NFU2,H2,XQ2,XJ2,2)
      IF(NKOORD.GE.3) CALL HERPOL(M3,NFU3,H3,XQ3,XJ3,3)
      IF(NKOORD.GE.4) CALL HERPOL(M4,NFU4,H4,XQ4,XJ4,4)
      IF(NKOORD.GE.5) CALL HERPOL(M5,NFU5,H5,XQ5,XJ5,5)
      WRITE(6,'(/1x,"INTEGRATING HAMILTONIAN")') 
      DO 101 I=1,NN
      DO 101 K=1,NKOORD
  101 IP(I,K)=IP(I,K)+1
      CALL MATREL(NKOORD,NN,NV,M1,H1,XQ1,XJ1,M2,H2,XQ2,XJ2,
     &            M3,H3,XQ3,XJ3,M4,H4,XQ4,XJ4,M5,H5,XQ5,XJ5,JROT)
      write(6,'(/1x,"INTEGRATION FINISHED")')
      CALL TIMER(cout)
  311 FORMAT(/1X,'ADDING DIAGONAL STR-BEND CORRECTION')                               0000018
      IF(NXIJ.NE.0) THEN ! ADD CORRECTION TERM TO DIAGONAL ELEMENTS
        WRITE(6,311)
        DO I=1,NN !LOOP OVER ALL DIAGONAL ELEMENTS
          ! GET QUANTUM NUMBERS FOR DIAGONAL ELEMENT
          DO J=1,NKOORD
            VK(J)=IP(I,J)-1
          ENDDO
          SU=0.D0
          JE=NKOORD-1
          IF(JE.EQ.0) JE=1 !SPECIAL TREATMENT OF 1D CASE
          DO J=1,JE !LOOP OVER BENDS
            DO K=1,NKOORD   !LOOP OVER STRETCHES
              SU=SU+XIJ(K,J)*(VK(K)+0.5D0)
            ENDDO
            DSB(I)=SU
            XK(I,I)=XK(I,I)+DSB(I)
          ENDDO
        ENDDO
      ENDIF
      IF(NDRU.NE.0.AND.(NN.LE.15.OR.NDRU.GT.1)) THEN
  204 FORMAT(/,1X,'HAMILTONIAN'/)                                       0000018
      WRITE(6,204)                                                      0000100
      DO 3 I=1,NN                                                       0000101
    3 WRITE(6,112) (XK(I,J)*EFAKT,J=1,NN)                                0000102
      ENDIF
C 412 FORMAT(/1X,'DIAGONAL ELEMENTS WITHOUT S-B-CORRECTION [CM-1]'/)
C     IF(NDRU.NE.0.OR.NXIJ.NE.0) THEN
C     WRITE(6,412)
C     WRITE(6,410) ((XK(I,I)-DSB(I))*EFAKT,I=1,NN)
C     ENDIF
C 413 FORMAT(/1X,'DIAGONAL ELEMENTS WITH S-B-CORRECTION [CM-1]'/)
C     IF(NXIJ.NE.0) THEN
C     WRITE(6,413)
C     WRITE(6,410) (XK(I,I)*EFAKT,I=1,NN)
C     ENDIF
  312 FORMAT(/,1X,'DIAGONALIZING HAMILTONIAN')                                   0000018
      WRITE(6,312)
      write(cout,'("DIAGONALIZATION FINISHED")')
      CALL DIAG(NNMAX,NN,XK,EVAL,XK,0)
      write(6,'(/1x,"DIAGONALIZATION FINISHED")')
      CALL TIMER(cout)
  403 FORMAT(//1X,'BANDENZENTREN IN REZIPROKEN ZENTIMETERN'/)
  404 FORMAT(/1X,'EIGENVALUES [CM-1]'/)
  405 FORMAT(/1X,'EIGENVECTORS OF THE HAMILTONIAN')
  410 FORMAT(1X,10F13.3)
      DO 401 I=1,NN
  401 EVAL(I)=EVAL(I)*EFAKT
      IF(NDRU.NE.0) THEN
      WRITE(6,404)
      WRITE(6,410) (EVAL(I),I=1,NN)
      ENDIF
      E0=EVAL(1)
      DO 409 I=1,NN
  409 EVAL(I)=EVAL(I)-E0
      WRITE(6,405)
      IF(NEIG.GT.NN) NEIG=NN
      DO 408 I=1,NEIG
      CMAX=0.D0
      XQ=0.d0
      DO J=1,NN
         IF(ABS(XK(J,I)).GT.CMAX) THEN
            CMAX=ABS(XK(J,I))
            IZ=J
         ENDIF
         CC=XK(J,I)
         DO K=1,NN
           CCC=CC*XK(K,I)
           DO II=1,NKOORD
             ! CHECK V[II]
             NVI=IP(J,II)-1
             NVJ=IP(K,II)-1
             IF(ABS(NVI-NVJ).NE.1) CYCLE
             ICHK=0
             DO JJ=1,NKOORD
               IF(JJ.EQ.II) CYCLE
               IF(ABS(IP(J,JJ)-IP(K,JJ)).NE.0) ICHK=ICHK+1 
             ENDDO
             IF(ICHK.NE.0) CYCLE
             CALL Q1D(NVI-NVJ,NVI,CCC)
             XQ(II)=XQ(II)+CCC
           ENDDO
         ENDDO
      ENDDO
  406 FORMAT(A1,A24,I1,A22,A1)
  415 FORMAT(A1,A24,A13,A1)
      IF(NKOORD-1.NE.0) THEN
        WRITE(IOStr1,406) '(','/1X,I3,2X,2F11.2,2x,"(",',
     &                         NKOORD-1,
     &                        '(I3,"/"),I3,")",F11.6/',')'
      ELSE
        WRITE(IOStr1,415) '(','/1X,I3,2X,2F11.2,2x,"(",',
     &                        'I3,")",F11.6/',')'
      ENDIF
  414 FORMAT(6X,5(10F11.6))
  411 FORMAT(/1X,I3,2X,2F11.2/)
      IF(CMAX.GT.0.8D0) THEN
        WRITE(6,IOStr1) I,EVAL(I)+E0,EVAL(I),(IP(IZ,J)-1,J=1,NKOORD),
     &                  CMAX
      ELSE
        WRITE(6,411) I,EVAL(I)+E0,EVAL(I)
      ENDIF
      IZ=0
      DO J=1,NN
        IF(XK(J,I)**2.GT.EPS) THEN
          IZ=IZ+1
          WFU(IZ)=XK(J,I)
          DO K=1,NKOORD
            MASSI(IZ,K)=IP(J,K)-1
          ENDDO
        ENDIF
      ENDDO
!      WRITE(6,414) (WFU(J),J=1,IZ)
      NL=IZ/10
      IF(MOD(IZ,10).NE.0) NL=NL+1
!      WRITE(6,'(1x,A3,I3,A5,I3)') 'NL=',NL,'REST=',MOD(IZ,10)
  407 FORMAT(A1,A10,I1,A17,A1)
      WRITE(IOStr2,407) '(','6x,10("|",',NKOORD,'I2,":",F10.6),"|"',')'
      DO J=1,NL
!        WRITE(6,'(1x,A2,I3)') 'J=',J
        KA=10*(J-1)+1
        KE=KA+9
        IF(KE.GT.IZ) KE=KA+MOD(IZ,10)-1
!        WRITE(6,'(1x,A6,2I3)') 'KA,KE=',KA,KE
        WRITE(6,IOStr2) ((MASSI(K,L),L=1,NKOORD),WFU(K),K=KA,KE)
      ENDDO
      IF(NKOART.EQ.0) THEN
        WRITE(6,'(/,3(3x,"Q(",I1,"):",F10.5))')
     &    (II,XQ(II),II=1,NKOORD)
        ZZ=ZKOORD
        DO II=1,3
          XQQ=XQ(II)*YLAM(II)
          DO IJ=1,NKP1
            ZZ(IJ)=ZZ(IJ)+XQQ*TRACART(3*IJ,II)
          ENDDO
        ENDDO
        WRITE(6,'(4(3x,"Z(",I1,"):",F10.5))')
     &    (II,ZZ(II)*FA2,II=1,NKP1)
        AMUE=0.d0
        DO II=1,3
          AMUE=AMUE+AK0(II)*XQ(II)*YLAM(II)
        ENDDO
        AMUE=INERT+0.5d0*AMUE
        AMUE=0.5d0*INERT/(AMUE*AMUE)
        WRITE(6,'(3x," MUE:",F15.7,F15.2)') AMUE*EFAKT,AMUE*EFAKT*CVEL
        if(I.EQ.1) BROT=AMUE
        IF(I.NE.1) THEN
        WRITE(6,'(3x,"DMUE:",F15.7,F15.2)') (AMUE-BROT)*EFAKT,
     &  (AMUE-BROT)*EFAKT*CVEL
        ENDIF
      ENDIF
  408 CONTINUE
      DO I=1,NN
         DO J=1,NN
           XK(I,J)=0.D0
         ENDDO
      ENDDO
!     CALL MOMENTUM(NKOORD,NN,YLAM,IP,XK)
!     DO I=1,NN
!        WRITE(6,111) (XK(I,J),J=1,I)
!     ENDDO
!     WRITE(6,'()')
      write(6,'(/1x,"LINVIB FINISHED")')
      CALL TIMER(cout)
   99 STOP
      END
      SUBROUTINE KRAFTK(NAT,NV,F)
!     GENERATE THE FORCE CONSTANT MATRIX F VIA EXTRAKTION OF
!     QUADRATIC PART OF THE POTENTIAL
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/POTE/NU(300,5),CV(300),NPOT,NKOART
      DIMENSION F(10,10)
      DIMENSION IK(2)
!     INITIALIZE F
      DO I=1,NAT
        DO J=1,NAT
          F(I,J)=0.D0
        END DO
      END DO
!     LOOP OVER ALL TERMS IN THE POTENTIAL
      DO I=1,NV
        SU=0
!     LOOP OVER ALL INTERNAL COORDINATES
        DO J=1,NAT
          SU=SU+NU(I,J)
        END DO
!     DID WE FIND A QUADRATIC TERM?
        IF(SU.EQ.2) THEN
          IZ=1
          FAK=1.
!     CHECK IF WE HAVE A DIAGONAL OR COUPLING TERM
          DO J=1,NAT
            IF(NU(I,J).EQ.2) THEN
!     FOR DIAGONAL
              IK(1)=J
              IK(2)=J
              FAK=2
            ENDIF
            IF(NU(I,J).EQ.1) THEN
!     FOR COUPLING
              IK(IZ)=J
              IZ=IZ+1
            ENDIF
          END DO
!     WRITE ELEMENT TO CORRECT POSITION IN F-MATRIX
          F(IK(1),IK(2))=CV(I)*FAK
          F(IK(2),IK(1))=F(IK(1),IK(2))
        ENDIF
      END DO
      RETURN
      END
      SUBROUTINE TIMER(text)
      IMPLICIT NONE
      INTEGER,PARAMETER ::I4B =SELECTED_INT_KIND(9)
      INTEGER,PARAMETER ::DP =KIND(1.0D0)
      character(len=80),INTENT(IN) :: text
      REAL(DP), SAVE :: T=0.D0,TR=0.D0
      REAL(DP), SAVE :: TSU=0.D0,TRSU=0.D0
      REAL(DP) :: TCPU,TREAL
      INTEGER(I4B) :: CLOCK,CLOCKR,CLOCKM
      IF(TR.eq.0.D0) THEN
        call system_clock(CLOCK, CLOCKR, CLOCKM)
        TR=CLOCK/CLOCKR
      ENDIF
      call CPU_TIME(TCPU)
      call system_clock(CLOCK, CLOCKR, CLOCKM)
      TREAL=CLOCK/CLOCKR
      TSU=TSU+(TCPU-T)
      TRSU=TRSU+(TREAL-TR)
!     WRITE(6,1000) text
      WRITE(6,1001) TCPU-T,(TREAL-TR)
      WRITE(6,1002) TSU,TRSU
      T=TCPU
      TR=TREAL
      RETURN
 1000 FORMAT(1x,A80)
 1001 FORMAT(1x,"STEP  CPU:",F10.3," s  WALL:",F10.2," s")
 1002 FORMAT(1x,"TOTAL CPU:",F10.3," s  WALL:",F10.2," s")
      END SUBROUTINE TIMER
      SUBROUTINE Q1D(DELV,VS,Y)
      ! Integral over 1D Harmonic Oscillator Functions
      ! Type: <v|Q|v'>
      ! DELV = v-v'
      ! VS = v
    
      IMPLICIT NONE
      INTEGER,PARAMETER ::I4B =SELECTED_INT_KIND(9)
      INTEGER,PARAMETER ::DP =KIND(1.0D0)
      
      INTEGER(I4B), INTENT(IN)    :: DELV,VS
      REAL(DP),     INTENT(INOUT) :: Y
     
      REAL(DP) :: P
    
      SELECT CASE (DELV)
        CASE (-1)
          P=0.5d0*(VS+1)
          Y=Y*SQRT(P)
        CASE (1)
          P=0.5d0*VS
          Y=Y*SQRT(P)
      END SELECT            
      RETURN  
      END SUBROUTINE
    
