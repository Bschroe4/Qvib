      SUBROUTINE DIAG (M,N,A,D,X,NORD)                                          
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
C                                                                               
C      COMPUTATION OF ALL EIGENVALUES AND EIGENVECTORS OF A REAL                
C      SYMMETRIC MATRIX BY THE METHOD OF QR TRANSFORMATIONS.                    
C      IF THE EUCLIDEAN NORM OF THE ROWS VARIES   S T R O N G L Y               
C      MOST ACCURATE RESULTS MAY BE OBTAINED BY PERMUTING ROWS AND              
C      COLUMNS TO GIVE AN ARRANGEMENT WITH INCREASING NORMS OF ROWS.            
C                                                                               
C      TWO MACHINE CONSTANTS MUST BE ADJUSTED APPROPRIATELY,                    
C      EPS = MINIMUM OF ALL X SUCH THAT 1+X IS GREATER THAN 1 ON THE            
C            COMPUTER,                                                          
C      TOL = INF / EPS  WITH INF = MINIMUM OF ALL POSITIVE X REPRESEN-          
C            TABLE WITHIN THE COMPUTER.                                         
C      A DIMENSION STATEMENT E(160) MAY ALSO BE CHANGED APPROPRIATELY.          
C                                                                               
C      INPUT                                                                    
C                                                                               
C      (M)   NOT LARGER THAN 160,  CORRESPONDING VALUE OF THE ACTUAL            
C            DIMENSION STATEMENT A(M,M), D(M), X(M,M),                          
C      (N)   NOT LARGER THAN (M), ORDER OF THE MATRIX,                          
C      (A)   THE MATRIX TO BE DIAGONALIZED, ITS LOWER TRIANGLE HAS TO           
C            BE GIVEN AS  ((A(I,J), J=1,I), I=1,N),                             
C                                                                               
C      OUTPUT                                                                   
C                                                                               
C      (D)   COMPONENTS D(1), ..., D(N) HOLD THE COMPUTED EIGENVALUES           
C            IN ASCENDING SEQUENCE. THE REMAINING COMPONENTS OF (D) ARE         
C            UNCHANGED,                                                         
C      (X)   THE COMPUTED EIGENVECTOR CORRESPONDING TO THE J-TH EIGEN-          
C            VALUE IS STORED AS COLUMN (X(I,J), I=1,N). THE EIGENVECTORS        
C            ARE NORMALIZED AND ORTHOGONAL TO WORKING ACCURACY. THE             
C            REMAINING ENTRIES OF (X) ARE UNCHANGED.                            
C                                                                               
C      ARRAY (A) IS UNALTERED. HOWEVER, THE ACTUAL PARAMETERS                   
C      THE EIGENVECTORS ON (A).                                                 
C      LEIBNIZ-RECHENZENTRUM, MUNICH 1965                                       
C                                                                               
      DIMENSION   A(M,M), D(M), X(M,M)                                          
      DIMENSION   E(3000)                                                      
C                                                                               
C     CORRECT ADJUSTMENT FOR IBM 360/91 DOUBLE PRECISION                        
C                                                                               
      EPS=2.5D-16                                                               
      TOL=2.5D-63                                                               
C                                                                               
      IF(N.EQ.1) GO TO 400                                                      
      DO 10 I=1,N                                                               
      DO 10 J=1,I                                                               
   10 X(I,J)=A(I,J)                                                             
C                                                                               
C     SIMULATION OF LOOP DO 150 I=N,2,(-1)                                      
C                                                                               
      DO 150 NI=2,N                                                             
      I=N+2-NI                                                                  
      L=I-2                                                                     
      H=0.0                                                                     
      G=X(I,I-1)                                                                
      IF(L) 140,140,20                                                          
   20 DO 30 K=1,L                                                               
   30 H=H+X(I,K)**2                                                             
      S=H+G*G                                                                   
      IF(S.GE.TOL) GO TO 50                                                     
   40 H=0.0                                                                     
      GO TO 140                                                                 
   50 IF(H) 140,140,60                                                          
   60 L=L+1                                                                     
      F=G                                                                       
      G=DSQRT(S)                                                                
      IF(F) 75,75,70                                                            
   70 G=-G                                                                      
   75 H=S-F*G                                                                   
      X(I,I-1)=F-G                                                              
      F=0.0                                                                     
C                                                                               
      DO 110 J=1,L                                                              
      X(J,I)=X(I,J)/H                                                           
      S=0.0                                                                     
      DO 80 K=1,J                                                               
   80 S=S+X(J,K)*X(I,K)                                                         
      J1=J+1                                                                    
      IF(J1.GT.L) GO TO 100                                                     
      DO 90 K=J1,L                                                              
   90 S=S+X(K,J)*X(I,K)                                                         
  100 E(J)=S/H                                                                  
  110 F=F+S*X(J,I)                                                              
C                                                                               
      F=F/(H+H)                                                                 
C                                                                               
      DO 120 J=1,L                                                              
  120 E(J)=E(J)-F*X(I,J)                                                        
C                                                                               
      DO 130 J=1,L                                                              
      F=X(I,J)                                                                  
      S=E(J)                                                                    
      DO 130 K=1,J                                                              
  130 X(J,K)=X(J,K)-F*E(K)-X(I,K)*S                                             
C                                                                               
  140 D(I)=H                                                                    
  150 E(I-1)=G                                                                  
C                                                                               
C     ACCUMULATION OF TRANSFORMATION MATRICES                                   
C                                                                               
  160 D(1)=X(1,1)                                                               
      X(1,1)=1.0                                                                
      DO 220 I=2,N                                                              
      L=I-1                                                                     
      IF(D(I)) 200,200,170                                                      
  170 DO 190 J=1,L                                                              
      S=0.0                                                                     
      DO 180 K=1,L                                                              
  180 S=S+X(I,K)*X(K,J)                                                         
      DO 190 K=1,L                                                              
  190 X(K,J)=X(K,J)-S*X(K,I)                                                    
  200 D(I)=X(I,I)                                                               
      X(I,I)=1.0                                                                
  210 DO 220 J=1,L                                                              
      X(I,J)=0.0                                                                
  220 X(J,I)=0.0                                                                
C                                                                               
C     DIAGONALIZATION OF THE TRIDIAGONAL MATRIX                                 
C                                                                               
      B=0.0                                                                     
      F=0.0                                                                     
      E(N)=0.0                                                                  
C                                                                               
      DO 340 L=1,N                                                              
      H=EPS*(DABS(D(L))+DABS(E(L)))                                             
      IF (H.GT.B) B=H                                                           
C                                                                               
C     TEST FOR SPLITTING                                                        
C                                                                               
      DO 240 J=L,N                                                              
      IF (DABS(E(J)).LE.B) GOTO 250                                             
  240 CONTINUE                                                                  
C                                                                               
C     TEST FOR CONVERGENCE                                                      
C                                                                               
  250 IF(J.EQ.L) GO TO 340                                                      
C                                                                               
C     SHIFT FROM UPPER 2*2 MINOR                                                
C                                                                               
  260 P=(D(L+1)-D(L))*0.5/E(L)                                                  
      R=DSQRT(P*P+1.0)                                                          
      IF(P) 270,280,280                                                         
  270 P=P-R                                                                     
      GO TO 290                                                                 
  280 P=P+R                                                                     
  290 H=D(L)-E(L)/P                                                             
      DO 300 I=L,N                                                              
  300 D(I)=D(I)-H                                                               
      F=F+H                                                                     
C                                                                               
C     QR TRANSFORMATION                                                         
C                                                                               
      P=D(J)                                                                    
      C=1.0                                                                     
      S=0.0                                                                     
C                                                                               
C     SIMULATION OF LOOP DO 330 I=J-1,L,(-1)                                    
C                                                                               
      J1=J-1                                                                    
      DO 330 NI=L,J1                                                            
      I=L+J1-NI                                                                 
      G=C*E(I)                                                                  
      H=C*P                                                                     
C                                                                               
C     PROTECTION AGAINST UNDERFLOW OF EXPONENTS                                 
C                                                                               
      IF (DABS(P).LT.DABS(E(I))) GOTO 310                                       
      C=E(I)/P                                                                  
      R=DSQRT(C*C+1.0)                                                          
      E(I+1)=S*P*R                                                              
      S=C/R                                                                     
      C=1.0/R                                                                   
      GO TO 320                                                                 
  310 C=P/E(I)                                                                  
      R=DSQRT(C*C+1.0)                                                          
      E(I+1)=S*E(I)*R                                                           
      S=1.0/R                                                                   
      C=C/R                                                                     
  320 P=C*D(I)-S*G                                                              
      D(I+1)=H+S*(C*G+S*D(I))                                                   
      DO 330 K=1,N                                                              
      H=X(K,I+1)                                                                
      X(K,I+1)=X(K,I)*S+H*C                                                     
  330 X(K,I)=X(K,I)*C-H*S                                                       
C                                                                               
      E(L)=S*P                                                                  
      D(L)=C*P                                                                  
      IF (DABS(E(L)).GT.B) GO TO 260                                            
C                                                                               
C     CONVERGENCE                                                               
C                                                                               
  340 D(L)=D(L)+F                                                               
      IF(NORD.NE.0) GOTO 9988                                                   
C                                                                               
C     ORDERING OF EIGENVALUES                                                   
C                                                                               
      NI=N-1                                                                    
  350 DO 380I=1,NI                                                              
      K=I                                                                       
      P=D(I)                                                                    
      J1=I+1                                                                    
      DO 360J=J1,N                                                              
      IF(D(J).GE.P) GOTO 360                                                    
      K=J                                                                       
      P=D(J)                                                                    
  360 CONTINUE                                                                  
      IF (K.EQ.I) GOTO 380                                                      
      D(K) =D(I)                                                                
      D(I)=P                                                                    
      DO 370 J=1,N                                                              
      P=X(J,I)                                                                  
      X(J,I)=X(J,K)                                                             
  370 X(J,K)=P                                                                  
  380 CONTINUE                                                                  
C                                                                               
C     FIXING OF SIGN                                                            
C                                                                               
      DO 385 I=1,N                                                              
      PM=0                                                                      
      DO 386 J=1,N                                                              
      IF(PM.GT.DABS(X(J,I))) GOTO 386                                           
      PM =DABS(X(J,I))                                                          
      K=J                                                                       
  386 CONTINUE                                                                  
      IF(X(K,I).GE.0) GOTO 385                                                  
      DO 387 J=1,N                                                              
  387 X(J,I)=-X(J,I)                                                            
  385 CONTINUE                                                                  
 9988 CONTINUE                                                                  
  390 GO TO 410                                                                 
C                                                                               
C     SPECIAL TREATMENT OF CASE N = 1                                           
C                                                                               
  400 D(1)=A(1,1)                                                               
      X(1,1)=1.0                                                                
  410 RETURN                                                                    
      END                                                                       
      SUBROUTINE MATINV(A,N,DETERM)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION INDEX(10,2), IPIVOT(10), A(10,10), B(10,1), PIVOT(10)
      EQUIVALENCE (IROW,JROW), (ICOLUM,JCOLUM), (AMAX, T, SWAP)
C     INITIALIZATION
      IROW=0
      ICOLUM=0
      JROW=0
      JCOLUM=0
      M=1
      DETERM=1.0
      DO 20 J=1,N
   20 IPIVOT(J)=0
      DO 550 I=1,N
C     SEARCH FOR PIVOT ELEMENT
      AMAX=0.0
      DO 105 J=1,N
      IF (IPIVOT(J)-1) 60, 105, 60
   60 DO 100 K=1,N
      IF (IPIVOT(K)-1) 80, 100, 740
   80 IF(DABS(AMAX)-DABS(A(J,K)))85,100,100
   85 IROW=J
      ICOLUM=K
      AMAX=A(J,K)
  100 CONTINUE
  105 CONTINUE
      IPIVOT(ICOLUM)=IPIVOT(ICOLUM)+1
C     INTERCHANGE ROWS TO PUT PIVOT ELEMENT ON DIAGONAL
      IF (IROW-ICOLUM) 140, 260, 140
  140 DETERM=-DETERM
      DO 200 L=1,N
      SWAP=A(IROW,L)
      A(IROW,L)=A(ICOLUM,L)
  200 A(ICOLUM,L)=SWAP
      IF(M) 260, 260, 210
  210 DO 250 L=1, M
      SWAP=B(IROW,L)
      B(IROW,L)=B(ICOLUM,L)
  250 B(ICOLUM,L)=SWAP
  260 INDEX(I,1)=IROW
      INDEX(I,2)=ICOLUM
      PIVOT(I)=A(ICOLUM,ICOLUM)
      DETERM=DETERM*PIVOT(I)
C     DIVIDE PIVOT ROW BY PIVOT ELEMENT
      A(ICOLUM,ICOLUM)=1.0
      DO 350 L=1,N
  350 A(ICOLUM,L)=A(ICOLUM,L)/PIVOT(I)
      IF(M) 380, 380, 360
  360 DO 370 L=1,M
  370 B(ICOLUM,L)=B(ICOLUM,L)/PIVOT(I)
C     REDUCE NON-PIVOT ROWS
  380 DO 550 L1=1,N
      IF(L1-ICOLUM) 400, 550, 400
  400 T=A(L1,ICOLUM)
      A(L1,ICOLUM)=0.0
      DO 450 L=1,N
  450 A(L1,L)=A(L1,L)-A(ICOLUM,L)*T
      IF(M) 550, 550, 460
  460 DO 500 L=1,M
  500 B(L1,L)=B(L1,L)-B(ICOLUM,L)*T
  550 CONTINUE
C     INTERCHANGE COLUMNS
C     LSTSQ10
      DO 710 I=1,N
      L=N+1-I
      IF (INDEX(L,1)-INDEX(L,2)) 630, 710, 630
  630 JROW=INDEX(L,1)
      JCOLUM=INDEX(L,2)
      DO 705 K=1,N
      SWAP=A(K,JROW)
      A(K,JROW)=A(K,JCOLUM)
      A(K,JCOLUM)=SWAP
  705 CONTINUE
  710 CONTINUE
  740 RETURN
      END
