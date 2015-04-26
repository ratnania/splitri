MODULE bnet
CONTAINS      

        SUBROUTINE BS3DM(L,M,N,A_SIZE,A,B,BSIZE1,BSIZE2,D)                    
!*THIS SUBROUTINE PRODUCES THE B-NETS OF ANY BOX SPLINES ON A THREE      
! DIRECTIONAL MESH. THE INPUT IS THE REPETITIONS (L,M,N) OF DIRECTIONS   
! (1,0), (0,1), AND (1,1). HERE L, M, AND N MUST BE GREATER THAN OR EQUAL
! TO 1. THE OUTPUT IS AN ARRAY OF THE COLLECTION OF B-NETS OF ALL        
! POLYNOMIAL PIECES OF THE BOX SPLINE INSIDE ITS SUPPORT.                
!      INPUT PARAMETERS:                                                 
!      L, M, N ARE INTEGERS. L IS THE NUMBER OF THE REPETITION OF (1,0)  
! M IS THE NUMBER OF THE REPETITION OF (0,1). N IS THE NUMBER OF THE     
! REPETITION OF (1,1).                                                   
!      SIZE IS A PARAMETER FOR STORAGE REQUIREMENT WHICH MUST BE GREATER 
! THAN THE LARGER ONE OF (L+M+N-2)*(L+N)+1 AND (L+M+N-2)*(M+N)+1.        
!      OUTPUT PARAMETERS:                                                
!      B IS AN OUTPUT ARRAY OF SIZE BSIZE1*BSIZE2. B CONTAINS THE B-NETS 
! OF THE BOX SPLINE (L+M+N-2)!*B(L,M,N).                                 
        IMPLICIT NONE
        INTEGER, INTENT(OUT) :: BSIZE1, BSIZE2, D
        INTEGER, INTENT(IN) :: A_SIZE
        INTEGER, INTENT(IN) :: L,M,N
        INTEGER, DIMENSION(0:A_SIZE,0:A_SIZE), INTENT(OUT) :: A
        INTEGER, DIMENSION(0:A_SIZE,0:A_SIZE), INTENT(OUT) :: B
! LOCAL PARAMETERS ARE L1, M1, N1, D, NX, NY, A. NOTE THAT D IS THE      
! DEGREE OF THE BOX SPLINE.                                              
        INTEGER :: I,J
        INTEGER :: L1,M1,N1,NX,NY                                         

        B(1,1)=1                                                         
        D=1                                                              
! WE START WITH THE B-NET OF BOX SPLINE B(1,1,1) AND FIRST FIND THE      
! B-NET OF BOX SPLINE B(L,1,1) IF L IS GREATER THAN 1.                   
        IF(L.GT.1)THEN                                                   
         L1=1                                                            
50       IF(L1.GE.L)GOTO 100                                             
! FIND THE DIFFERENCE OF THE B-NETS OF B(L1,1,1) AT (0,0) AND B(L1,1,1)  
! AT (L1,0) AND WRITE IT INTO THE ARRAY A.                               
          DO 60 I=0, D*(L1+2)                                            
           DO 55 J=0, 2*D                                                
            IF(I.LT.D)THEN                                               
             A(I,J)=B(I,J)                                               
            ELSE                                                         
             A(I,J)=B(I,J)-B(I-D,J)                                      
            END IF                                                       
55         CONTINUE                                                      
60        CONTINUE                                                       
         D=D+1                                                           
! SET THE INITIAL B-NET OF B(L1+1,1,1) TO BE ZERO.                       
         DO 65 I=0, D*(L1+2)                                             
          B(I,0)=0                                                       
65       CONTINUE                                                        
         DO 70 J=0,D*2                                                   
          B(0,I)=0                                                       
70       CONTINUE                                                        
! THEN SOLVE THE DIFFERENCE EQUATIONS TO FIND THE B-NET OF B(L1+1,1,1)   
         DO 95 I=0, L1+1                                                 
          DO 90 J=0,1                                                    
           DO 85 NX=1, D                                                 
            DO 80 NY=1, D                                                
             IF(NY.GE.NX)THEN                                            
        B(NX+I*D,NY+D*J)=B(NX-1+D*I,NY+D*J)+A(NX-1+(D-1)*I,NY-1+(D-1)*J) 
             ELSE                                                        
        B(NX+I*D,NY+D*J)=B(NX-1+D*I,NY+D*J)+A(NX-1+(D-1)*I,NY+(D-1)*J)   
             END IF                                                      
80          CONTINUE                                                     
85         CONTINUE                                                      
90        CONTINUE                                                       
95       CONTINUE                                                        
         L1=L1+1                                                         
         GOTO 50                                                         
        ENDIF                                                            
! THE ABOVE LOOP HAS FOUND THE B-NET OF B(L,1,1).                        
! NOW START TO FIND THE B-NET OF BOX SPLINE B(L,M,1).                    
100     IF(M.GT.1)THEN                                                   
        M1=1                                                             
150     IF(M1.GE.M)GOTO 200                                              
! FIND THE DIFFERENCE OF THE B-NETS OF B(L,M1,1) AT (0,0) AND B(L,M1,1)  
! AT (0,M1) AND WRITE IT INTO THE ARRAY A.                               
        DO 160 J=0, D*(M1+2)                                             
          DO 155 I=0,D*(L+1)                                             
           IF(J.LT.D)THEN                                                
            A(I,J)=B(I,J)                                                
           ELSE                                                          
            A(I,J)=B(I,J)-B(I,J-D)                                       
           END IF                                                        
155       CONTINUE                                                       
160      CONTINUE                                                        
        D=D+1                                                            
! SET THE INITIAL B-NET OF B(L,M1+1,1) TO BE ZERO.                       
        DO 165 I=0, D*(L+1)                                              
         B(I,0)=0                                                        
165     CONTINUE                                                         
        DO 170 J=0, D*(M1+2)                                             
         B(0,J)=0                                                        
170     CONTINUE                                                         
! THEN SOLVE THE DIFFERENCE EQUATIONS TO FIND THE B-NET OF B(L,M1+1,1)   
        DO 195 I=0, L                                                    
         DO 190 J=0, M1+1                                                
          DO 185 NY=1, D                                                 
           DO 180 NX=1,D                                                 
            IF(NX.GE.NY)THEN                                             
        B(NX+I*D,NY+J*D)=B(NX+I*D,NY-1+J*D)+A(NX+I*(D-1)-1,NY-1+J*(D-1)) 
             ELSE                                                        
        B(NX+I*D,NY+J*D)=B(NX+I*D,NY-1+J*D)+A(NX+I*(D-1),NY-1+J*(D-1))   
            ENDIF                                                        
180        CONTINUE                                                      
185       CONTINUE                                                       
190      CONTINUE                                                        
195     CONTINUE                                                         
        M1=M1+1                                                          
        GOTO 150                                                         
       ENDIF                                                             
! FINALLY  FIND THE B-NET OF BOX SPLINE B(L,M,N)                         
200     IF(N.GT.1)THEN                                                   
        N1=1                                                             
250     IF(N1.GE.N)GOTO 300                                              
! FIND THE DIFFERENCE OF THE B-NETS OF B(L,M,N1) AT (0,0) AND B(L,M,N1)  
! AT (N1,N1) AND WRITE IT INTO THE ARRAY A.                              
        DO 260 I=0,D*(L+N1+1)                                            
         DO 255 J=0,D*(M+N1+1)                                           
          IF((I.LT.D).OR.(J.LT.D))THEN                                   
           A(I,J)=B(I,J)                                                 
          ELSE                                                           
           A(I,J)=B(I,J)-B(I-D,J-D)                                      
          ENDIF                                                          
255      CONTINUE                                                        
260     CONTINUE                                                         
        D=D+1                                                            
! SET THE INITIAL B-NET OF B(L,M,N1+1) TO BE ZERO.                       
        DO 265 I=0, D*(L+N1+1)                                           
         B(I,0)=0                                                        
265     CONTINUE                                                         
        DO 270 J=0, D*(M+N1+1)                                           
         B(0,J)=0                                                        
270     CONTINUE                                                         
! THEN SOLVE THE DIFFERENCE EQUATIONS TO FIND THE B-NET OF B(L,M,N1+1)   
        DO 295 I=0, L+N1                                                 
         DO 290 J=0, M+N1                                                
          DO 285 NX=1, D                                                 
           DO 280 NY=1, D                                                
           B(NX+I*D,NY+J*D)= B(NX+I*D-1,NY+J*D-1)    &
        &                 +  A(NX+I*(D-1)-1,NY+J*(D-1)-1)              
280        CONTINUE                                                      
285       CONTINUE                                                       
290      CONTINUE                                                        
295     CONTINUE                                                         
        N1=N1+1                                                          
        GOTO 250                                                         
        ENDIF                                                            
300     PRINT*,'WE HAVE GENERATED THE B-NETS OF BOX SPLINE '             
        PRINT*, D,'!*B(',L,M,N,')'                                       
        BSIZE1= D*(L+N)+1                                                
        BSIZE2= D*(M+N)+1                                                

        RETURN                                                           
        END SUBROUTINE BS3DM                                                             






        SUBROUTINE BS4DM(K,L,M,N,A_SIZE,A,B,BSIZE1,BSIZE2,D1,D2)                  
! THIS PROGRAM FINDS THE B-NETS OF ANY BOX SPLINE ON A FOUR DIRECTIONAL   
! MESH. THE INPUT IS THE REPETITIONS (K,L,M,N) OF DIRECTIONS (1,0),(0,1), 
! (1,1), AND (1,-1). ALL K, L, M, AND N MUST BE GREATER OR EQUAL TO 1.    
! THE OUTPUT IS AN ARRAY OF COLLECTIONS OF B-NETS OF THE BOX SPLINE       
! RESTRICTED TO EVERY TRIANGLES INSIDE ITS SUPPORT.                       
!       INPUT PARAMETERS:                                                 
! K IS THE NUMBER OF THE REPETITION OF (1,0); L IS THE NUMBER OF THE      
! REPETITION OF (0,1); M IS THE NUMBER OF REPETITION OF (1,1); AND N IS   
! THE NUMBER OF THE REPETITION OF (1,-1).                                 
! SIZE IS THE PARAMETER TO DECLARE THE STORAGE REQUIREMENT OF THE OUTPUT  
! MATRIX B. SIZE MUST BE GREATER THAN THE INTERGERS 2*(K+L+M+N-2)*(K+M+N) 
! AND 2*(K+L+M+N-2)*(L+M+N).                                              
!       OUTPUT PARAMETERS:                                                
! B IS AN OUTPUT ARRAY OF SIZE BSIZE1*BSIZE2. B CONTAINS THE B-NETS OF    
! BOX SPLINE (K+L+M+N-2)!*2^{M+N+1}*B(K,L,M,N). ACTUALLY, THE B-NETS CAN  
! BE PRINTED OUT AS FOLLOWS.                                              
!       DO J=0,2*(K+L+M+N-2)*(L+M+N),2                                    
!       WRITE(*,1) (B(I,J),I=0,2*(K+L+M+N-2)*(K+M+N),2)                   
!       WRITE(*,2) (B(I+1,J+1),I=0,2*(K+L+M+N-2)*(K+M+N)-1,2)             
! 1     FORMAT(2X,20I4)                                                   
! 2     FORMAT(4X,20I4)                                                   
! BSIZE1 AND BSIZE2 ARE DESCRIBED AS ABOVE.                               
        IMPLICIT NONE
        INTEGER, INTENT(OUT) :: BSIZE1, BSIZE2
        INTEGER, INTENT(OUT) :: D1, D2
        INTEGER, INTENT(IN) :: A_SIZE
        INTEGER, INTENT(IN) :: K,L,M,N
        INTEGER, DIMENSION(0:A_SIZE,0:A_SIZE), INTENT(OUT) :: A
        INTEGER, DIMENSION(0:A_SIZE,0:A_SIZE), INTENT(OUT) :: B
! LOCAL PARAMETERS A,D,K1,L1,M1,N1,NX,NY.                                 
        INTEGER :: I,J,D,K1,L1,M1,N1,NX,NY                                       
! SET THE INITIAL B-NETS OF 16*B(1,1,1,1).                                
        B(:,:) = 0
        B(4,4)=4                                                          
        B(4,6)=8                                                          
        B(4,8)=4                                                          
        B(3,5)=4                                                          
        B(3,7)=4                                                          
        B(2,6)=2                                                          
        B(5,3)=4                                                          
        B(5,5)=8                                                          
        B(5,7)=8                                                          
        B(5,9)=4                                                          
        B(6,2)=2                                                          
        B(6,4)=8                                                          
        B(6,6)=8                                                          
        B(6,8)=8                                                          
        B(6,10)=2                                                         
        B(7,3)=4                                                          
        B(7,5)=8                                                          
        B(7,7)=8                                                          
        B(7,9)=4                                                          
        B(8,4)=4                                                          
        B(8,6)=8                                                          
        B(8,8)=4                                                          
        B(9,5)=4                                                          
        B(9,7)=4                                                          
        B(10,6)=2                                                         
        D=2                                                               
! SET DEGREE D=2 AND START TO FIND THE B-NETS OF BOX SPLINE B(K,1,1,1).   
        IF(K.GT.1)THEN                                                    
         K1=1                                                             
50       IF(K1.GE.K)GOTO 100                                              
! FIND THE DIFFERENCE OF THE B-NETS OF B(K1,1,1,1) AT (0,0) AND           
! B(K1,1,1,1)  AT (K1,0) AND WRITE IT INTO THE ARRAY A.                   
        DO 60 I=0, 2*D*(K1+3)                                             
          DO 55 J=0, 6*D                                                  
           IF(I.LT.2*D)THEN                                               
            A(I,J)=B(I,J)                                                 
           ELSE                                                           
            A(I,J)=B(I,J)-B(I-2*D,J)                                      
           END IF                                                         
55        CONTINUE                                                        
60      CONTINUE                                                          
        D=D+1                                                             
        DO 65 I=0, 2*D*(K1+3),2                                           
         B(I,0)=0                                                         
65      CONTINUE                                                          
        DO 70 J=0, 6*D,2                                                  
         B(0,J)=0                                                         
70      CONTINUE                                                          
! THEN SOLVE THE DIFFERENCE EQUATIONS TO FIND THE B-NET OF B(K1+1,1,1,1). 
        DO 95 I=0,K1+2                                                    
         DO 94 J=0,2                                                      
          DO 75 NX=1, D                                                   
           DO 74 NY=NX, 2*D-NX,2                                          
            B(I*2*D+NX,J*2*D+NY)=B(2*I*D+NX-1,J*2*D+NY-1) &                
     &      +B(2*I*D+NX-1,2*J*D+NY+1)+A(I*2*(D-1)+NX-1,J*2*(D-1)+NY-1)    
            B(I*2*D+NX,J*2*D+NY)=B(I*2*D+NX,J*2*D+NY)/2                   
74         CONTINUE                                                       
75        CONTINUE                                                        
          DO 80 NY=1,D-1                                                  
           DO 79 NX=NY+2, 2*D-NY, 2                                       
            B(I*2*D+NX,J*2*D+NY)=B(I*2*D+NX-2,J*2*D+NY)  &                 
     &                           +A(I*2*(D-1)+NX-2,J*2*(D-1)+NY)          
79         CONTINUE                                                       
80        CONTINUE                                                        
          DO 85 NY=D+1,2*D                                                
           DO 84 NX=2*D-NY+2,NY,2                                         
            B(I*2*D+NX,J*2*D+NY)=B(I*2*D+NX-2,J*2*D+NY)   &                
     &                           +A(I*2*(D-1)+NX-2,J*2*(D-1)+NY-2)        
84         CONTINUE                                                       
85        CONTINUE                                                        
          DO 90 NX=D+2, 2*D                                               
           DO 89 NY=2*D-NX+2,NX-2, 2                                      
           B(I*2*D+NX,J*2*D+NY)=B(I*2*D+NX-1,J*2*D+NY-1)*2 &               
     &         -B(I*2*D+NX,J*2*D+NY-2)+A(I*2*(D-1)+NX-2,J*2*(D-1)+NY-2)   
89         CONTINUE                                                       
90        CONTINUE                                                        
94       CONTINUE                                                         
95      CONTINUE                                                          
        K1=K1+1                                                           
        GOTO 50                                                           
       ENDIF                                                              
! THE ABOVE LOOP HAS FOUND THE B-NET OF B(K,1,1,1).                       
! START TO FIND THE B-NETS OF BOX SPLINE B(K,L,1,1).                      
100     IF(L.GT.1)THEN                                                    
        L1=1                                                              
150     IF(L1.GE.L)GOTO 200                                               
! FIND THE DIFFERENCE OF THE B-NETS OF B(K,L1,1,1) AT (0,0) AND           
! B(K,L1,1,1)  AT (0,L1) AND WRITE IT INTO THE ARRAY A.                   
        DO 160 I=0, 2*D*(K+2)                                             
          DO 155 J=0, 2*D*(L1+3)                                          
           IF(J.LT.2*D)THEN                                               
            A(I,J)=B(I,J)                                                 
           ELSE                                                           
            A(I,J)=B(I,J)-B(I,J-2*D)                                      
           ENDIF                                                          
155       CONTINUE                                                        
160     CONTINUE                                                          
        D=D+1                                                             
! SET THE INITIAL B-NET OF B(K,L1+1,1,1) TO BE ZERO.                      
        DO 165 I=0,2*D*(K+2),2                                            
         B(I,0)=0                                                         
165     CONTINUE                                                          
        DO 170 J=0,2*D*(L1+3),2                                           
         B(0,J)=0                                                         
170     CONTINUE                                                          
! THEN SOLVE THE DIFFERENCE EQUATIONS TO FIND THE B-NET OF B(K,L1+1,1,1)  
        DO 195 I=0, K+1                                                   
          DO 194 J=0, L1+2                                                
           DO 175 NY=1,D                                                  
            DO 174 NX=NY,2*D-NY,2                                         
             B(I*2*D+NX,J*2*D+NY)=(B(I*2*D+NX-1,J*2*D+NY-1)   &            
     &     +B(I*2*D+NX+1,J*2*D+NY-1)+A(I*2*(D-1)+NX-1,J*2*(D-1)+NY-1))/2  
174         CONTINUE                                                      
175        CONTINUE                                                       
           DO 180 NX=1,D-1                                                
            DO 179 NY=NX+2,2*D-NX,2                                       
             B(I*2*D+NX,J*2*D+NY)=B(I*2*D+NX,J*2*D+NY-2)      &            
     &                            +A(I*2*(D-1)+NX,J*2*(D-1)+NY-2)         
179         CONTINUE                                                      
180        CONTINUE                                                       
           DO 185 NX=D+1,2*D                                              
            DO 184 NY=2*D-NX+2,NX,2                                       
             B(I*2*D+NX,J*2*D+NY)=B(I*D*2+NX,J*2*D+NY-2)      &            
     &                            +A(I*2*(D-1)+NX-2,J*2*(D-1)+NY-2)       
184         CONTINUE                                                      
185        CONTINUE                                                       
           DO 190 NY=D+2,2*D                                              
            DO 189 NX=2*D-NY+2, NY-2, 2                                   
             B(I*2*D+NX,J*2*D+NY)=B(I*2*D+NX-1,J*2*D+NY-1)*2  &            
     &      -B(I*2*D+NX-2,J*2*D+NY)+A(I*2*(D-1)+NX-2,J*2*(D-1)+NY-2)      
189         CONTINUE                                                      
190        CONTINUE                                                       
194       CONTINUE                                                        
195      CONTINUE                                                         
        L1=L1+1                                                           
        GOTO 150                                                          
       ENDIF                                                              
! NOW COMPUTE THE B-NETS OF BOX SPLINE B(K,L,M,1).                        
200     IF(M.GT.1)THEN                                                    
         M1=1                                                             
250      IF(M1.GE.M)GOTO 300                                              
! FIND THE DIFFERENCE OF THE B-NETS OF B(K,L,M1,1) AT (0,0) AND           
! B(K,L,M1,1) AT (M1,M1) AND WRITE IT INTO THE ARRAY A.                   
         DO 260 I=0, 2*D*(K+M1+2)                                         
          DO 255 J=0, 2*D*(L+M1+2)                                        
           IF((I.LT.2*D).OR.(J.LT.2*D))THEN                               
            A(I,J)=B(I,J)                                                 
           ELSE                                                           
            A(I,J)=B(I,J)-B(I-2*D,J-2*D)                                  
           ENDIF                                                          
255       CONTINUE                                                        
260      CONTINUE                                                         
         D=D+1                                                            
! SET THE INITIAL B-NET OF B(K,L,M1+1,1) TO BE ZERO.                      
         DO 265 I=0,2*D*(K+M1+2),2                                        
          B(I,0)=0                                                        
265      CONTINUE                                                         
         DO 270 J=0, 2*D*(L+M1+2),2                                       
          B(0,J)=0                                                        
270      CONTINUE                                                         
         DO 295 I=0, K+M1+1                                               
          DO 294 J=0, L+M1+1                                              
           DO 275 NX=1, D                                                 
            DO 274 NY=NX, 2*D-NX,2                                        
             B(I*2*D+NX,J*2*D+NY)=B(I*2*D+NX-1,J*2*D+NY-1)        &        
     &                            +A(I*2*(D-1)+NX-1,J*2*(D-1)+NY-1)       
274         CONTINUE                                                      
275        CONTINUE                                                       
          DO 280 NY=1, D-1                                                
           DO 279 NX=NY+2,2*D-NY,2                                        
            B(I*2*D+NX,J*2*D+NY)=B(I*2*D+NX-1,2*J*D+NY-1)         &        
     &                           +A(I*2*(D-1)+NX-1,J*2*(D-1)+NY-1)        
279        CONTINUE                                                       
280       CONTINUE                                                        
          DO 285 NY=D+1,2*D                                               
           DO 284 NX=2*D-NY+2,NY,2                                        
            B(I*2*D+NX,J*2*D+NY)=B(I*2*D+NX-1,2*J*D+NY-1)         &        
     &                           +A(I*2*(D-1)+NX-2,J*2*(D-1)+NY-2)        
284        CONTINUE                                                       
285       CONTINUE                                                        
          DO 290  NX=D+2,2*D                                              
           DO 289 NY=2*D-NX+2,NX-2,2                                      
            B(I*2*D+NX,J*2*D+NY)=B(I*2*D+NX-1,2*J*D+NY-1)         &        
     &                           +A(I*2*(D-1)+NX-2,J*2*(D-1)+NY-2)        
289        CONTINUE                                                       
290       CONTINUE                                                        
294      CONTINUE                                                         
295     CONTINUE                                                          
        M1=M1+1                                                           
        GOTO 250                                                          
       ENDIF                                                              
! FINALLY, COMPUTE THE B-NETS OF BOX SPLINE B(K,L,M,N).                   
300     IF(N.GT.1)THEN                                                    
        N1=1                                                              
350     IF(N1.GE.N)GOTO 400                                               
! FIND THE DIFFERENCE OF THE B-NETS OF B(K,L,M,N1) AT (0,0) AND           
! B(K,L,M,N1) AT (N1,-N1) AND WRITE IT INTO THE ARRAY A.                  
        DO 360 I=0, 2*D*(K+M+N1+1)                                        
           DO 355 J=0, 2*D*(L+M+N1+1)                                     
            IF((I.LT.2*D).AND.(J.LT.2*D))THEN                             
              A(I,J)=0                                                    
            ENDIF                                                         
            IF((I.LT.2*D).AND.(J.GE.2*D))THEN                             
              A(I,J)=-B(I,J-2*D)                                          
            ENDIF                                                         
            IF((I.GE.2*D).AND.(J.LT.2*D))THEN                             
              A(I,J)=B(I-2*D,J)                                           
            ENDIF                                                         
            IF((I.GE.2*D).AND.(J.GE.2*D))THEN                             
              A(I,J)=B(I-2*D,J)-B(I,J-2*D)                                
            ENDIF                                                         
355        CONTINUE                                                       
360      CONTINUE                                                         
         D=D+1                                                            
! SET THE INITIAL B-NET OF B(K,L,M, N1+1) TO BE ZERO.                     
         DO 365 I=0, 2*D*(K+M+N1+1),2                                     
            B(I,0)=0                                                      
365      CONTINUE                                                         
         DO 370 J=0, 2*D*(L+M+N1+1),2                                     
          B(2*D*(K+M+N1+1),J)=0                                           
370      CONTINUE                                                         
         DO 395 I=K+M+N1,0,-1                                             
          DO 394 J=0,L+M+N1                                               
           DO 375 NY=1,D                                                  
            DO 374 NX=NY,2*D-NY,2                                         
             B(I*2*D+NX,J*2*D+NY)=B(I*2*D+NX+1,J*2*D+NY-1)         &       
     &                            +A(I*2*(D-1)+NX-1,J*2*(D-1)+NY-1)       
374         CONTINUE                                                      
375        CONTINUE                                                       
           DO 380 NX=2*D-1,D+1,-1                                         
            DO 379 NY=2*D-NX+2,NX,2                                       
             B(I*2*D+NX,J*2*D+NY)=B(I*2*D+NX+1,J*2*D+NY-1)         &       
     &                              +A(I*2*(D-1)+NX-1,J*2*(D-1)+NY-1)     
379         CONTINUE                                                      
380        CONTINUE                                                       
           DO 385 NX=D-1,0,-1                                             
            DO 384 NY=NX+2,2*D-NX,2                                       
             B(I*2*D+NX,J*2*D+NY)=B(I*2*D+NX+1,J*2*D+NY-1)         &       
     &                            +A(I*2*(D-1)+NX,J*2*(D-1)+NY-2)         
384         CONTINUE                                                      
385        CONTINUE                                                       
           DO 390 NY=D+2,2*D                                              
            DO 389 NX=2*D-NY+2,NY-2,2                                     
             B(I*2*D+NX,J*2*D+NY)=B(I*2*D+NX+1,J*2*D+NY-1)         &       
     &                            +A(I*2*(D-1)+NX,J*2*(D-1)+NY-2)         
389         CONTINUE                                                      
390        CONTINUE                                                       
394       CONTINUE                                                        
395      CONTINUE                                                         
         N1=N1+1                                                          
         GOTO 350                                                         
         ENDIF                                                            
400      PRINT*, 'WE HAVE GENERATED THE B-NETS OF BOX SPLINE '            
         PRINT*, D,'!*','2^(',M+N+1,')B(',K,L,M,N,').'                    
         BSIZE1=2*D*(K+M+N)+1                                             
         BSIZE2=2*D*(L+M+N)+1                                             
         D1 = D ; D2 = 2**(M+N+1)
         RETURN                                                           
         END SUBROUTINE BS4DM                                                             

END MODULE BNET 
