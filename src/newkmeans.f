C++*********************************************************************
C
C SUBKMEANS.F
C                  LUNF, * BUG                    AUG 2009 ARDEAN LEITH
C
C **********************************************************************
C
C  NEWKMEANS.F
C
C  THIS ROUTINE WILL PARTITION A NUMBER OF OBJECTS NOBJ, INTO A NUMBER OF
C  CLASSES NCLASS. AN OBJECT "I" BELONGS TO A CLASS "J" IF THE "EUCLIDIAN
C  DISTANCE" OF THIS OBJECT "I" TO THE CENTROID OF THE CLASS "J" IS THE 
C  SMALLEST OF DISTANCES TO ALL CLASS CENTROIDS.
C  A CLASS WILL HAVE MORE OBJECTS IN IT AS ALL OBJECT ARE PROCESSED.
C  THUS A CENTROID COORDINATES ARE CONTINUALLY COMPUTED AND CHANGING.
C  THIS CAUSES AN OBJECT "I" IN ONE CLASS TO FIND ITSELF CLOSER TO 
C  A CENTROID OF ANOTHER CLASS.
C  THAT OBJECT "I" IS REMOVED FROM ONE CLASS AND PUT INTO ANOTHER. 
C
C  REFERENCE: CLUSTER ANALYSIS ALGORITHMS
C             FOR DATA REDUCTION AND CLASSIFICATION OF OBJECTS.
C  BY: HELMUTH SPATH.
C
C  PUBLISHER: JOHN WILEY & SONS
C             ELLIS HORWOOD LIMITED,  1980
C
C  PAGE 72:

C	THE INITIAL ASSIGNMENT OF THE VECTORS X(NOBJ,NFAC) TO NCLASS CLUSTERS
C	IS GIVEN  BY THE ARRAY P(), WHERE P(I) IS THE CLUSTER NUMBER
C  	OF THE I-th VECTOR; THUS EACH P(I) MUST BE SUCH THAT 
C	1 <= P(I) =< NCLASS AND FOR EACH J=1,...,N AT LEAST ONE I WITH
C	P(I) = J MUST EXIST.
C
C	LET D DENOTE THE SUM OF THE E(J), WHERE E(J) IS THE SUM
C	OF THE SQUARES OF THE DISTANCES BETWEEN THE MEMBERS OF 
C	THE J-th CLUSTER AND THEIR CENTROIDS.
C
C	THE SUBROUTINE MINIMISES D AS FAR AS POSSIBLE BY
C	REPEATED EXCHANGES OF CLUSTER MEMBERS, THE P(I)
C	ARE CORRESPONDINGLY MODIFIED WITHOUT, HOWEVER,
C	CHANGING THE NUMBER OF CLUSTERS.
C
C	THE SUBROUTINE RETURNS THE VALUES AT THE FINAL
C	CONFIGURATION FOR THE CENTROIDS S(NCLASS,NFAC), THE 
C	SUMS E(NCLASS), AND D.
C
C	IF IDR = 1, THE CURRENT VALUES OF D AND OF THE VECTOR P
C	ARE PRINTED AT EACH ITERATION.
C
C	NOBJ	       : NUMBER OF OBJECTS
C       NCLASS         : NUMBER OF CLUSTERS
C       NFAC	       : NUMBER OF DIMENSION IN THE COORD. SYSTEM
C	X(NFAC,NOBJ)   : COORDINATES OF EACH OBJECT
C	P(NOBJ)	       : CLUSTER TO WHICH OBJECT I BELONGS.
C	S(NFAC,NCLASS) : COORDINATES OF THE CENTROIDS OF EACH CLUSTER.
C	E(NCLASS)      : SUM OF THE SQUARED DISTANCES OF ALL OBJECTS 
C		         (IN A CLUSTER) TO THEIR CENTROID.
C	D	       : SUM OF E(NCLASS)
C                        (THIS HAS TO BE REDUCED AFTER EACH NEW CLUSTERING)
C	Q(NCLASS)      : NUMBER OF OBJECTS IN EACH CLUSTER.
C       INUM(NFAC)     :                                          (SENT)
C       TYPE           : TYPE OF INPUT FILE e.g. _PIX, _IMC, _SEQ
C
C	NOTE: P() HAS BEEN INITIALIZED BEFORE THE CALL TO THIS ROUTINE.
C             IF (NOBJ/NCLASS < 5) SOME GOOD OPTIMAL PARTITIONS ARE
C             TO ADD THE POSSIBILITY OF WORKING IN MEMORY OR ON DISK,
C	      I ADDED THE READ OF X() DATA HERE.
C--***************************************************************************

        SUBROUTINE NEWKMEANS(NOBJ, NFAC, X, P, NCLASS, S, E, Q, D, IDR,
     &            INUM,NIMI,WT, COO, NFTOT, X_DIM, USE_DISK,LUNF,ITYPE)

        INCLUDE 'CMBLOCK.INC'

        INTEGER :: P, Q, R, U, V, W, IDR, NFTOT, PT, X_DIM

C	NOTE: IF (USE_DISK) THEN X_DIM=1 ELSE X_DIM=NOBJ

        DIMENSION X(NFAC,X_DIM), P(NOBJ), S(NFAC,NCLASS),INUM(NFAC),
     &      E(NCLASS), Q(NCLASS), COO(NFTOT),NIMI(NOBJ),WT(NFAC)

        LOGICAL :: DONE,USE_DISK

        IF (NOBJ/NCLASS .LT. 5) THEN
          WRITE(NOUT,*)' '
          WRITE(NOUT,*)' WARNING -------'
          WRITE(NOUT,*)' NUMBER OF OBJECTS / NUMBER OF CLASSES < 5'
          WRITE(NOUT,*)' SOME GOOD OPTIMAL PARTITIONS WILL NOT BE FOUND'
          WRITE(NOUT,*)'  '
        ENDIF

C	INITIALIZE TO ZERO, THE DISTANCES TO THE CENTROIDS AND
C	THE COORDINATES OF THE CENTROIDS.

c$omp   parallel do if(nfac.ge.64), private(j,k)
        DO J = 1, NCLASS
           Q(J) = 0
           E(J) = 0.
           DO K = 1, NFAC
              S(K,J) = 0.0
           ENDDO
        ENDDO

C	FOR EACH OBJECT, FIND ITS CLUSTER 'R', THE NUMBER OF OBJECT IN THE
C       CLUSTER 'Q(R)' AND ADD ITS COORDINATES 'X(I,K)' TO THE COORDINATES
C	OF THE CENTROID OF THE CLUSTER 'S(R,K)'.

        DO I = 1, NOBJ
	  IF (USE_DISK) THEN
             PT = 1
          ELSE
             PT = I
          ENDIF

C	  IN EITHER CASE, WE NEED TO READ X() AT LEAST ONCE.
          IF (ITYPE .EQ. 1) THEN
C            UNFORMATTED FILE FOR _SEQ IMAGE DATA
             READ(LUNF)   (COO(J), J = 1, NFTOT),FIM
          ELSE
C            FORMATTED FILE FOR IMAGE OR PIXEL COOR.
             READ(LUNF,*) (COO(J), J = 1, NFTOT),D1,D2,FIM,FDUM
          ENDIF
          NIMI(I) = FIM

c$omp     parallel do if(nfac.ge.64), private(j)
          DO J = 1, NFAC
              X(J,PT) = WT(J) * COO(INUM(J))
          ENDDO
          R = P(I)

          IF (R .LT. 1 .OR. R .GT. NCLASS) THEN
            WRITE(NOUT,*)' CLUSTER = ',R,' DOES NOT BELONG TO GROUP'
            WRITE(NOUT,*)' OBJECT I=',I,
     &                   ' HAS NOT BEEN PRESET TO ANY CLUSTER'
            RETURN
          ENDIF

          Q(R) = Q(R) + 1

c$omp     parallel do if(nfac.ge.64), private(k)
          DO K = 1, NFAC
             S(K,R) = S(K,R) + X(K,PT)
          ENDDO
        ENDDO 

C	REWIND THE UNIT TO READ THE DATA AGAIN       
        IF (USE_DISK) THEN
          IF (ITYPE .EQ. 1) THEN
C            UNFORMATTED FILE FOR _SEQ IMAGE DATA
             CALL REW(LUNF,1)
          ELSE
C            FORMATTED FILE FOR IMAGE OR PIXEL COOR.
             CALL REWF(LUNF,1)
          ENDIF
        ENDIF

C	COMPUTE THE CORRECT COORDINATES OF EACH CENTROID BY
C	DIVIDING THE SUM OF COORDS OF ALL ITS OBJECTS BY THE
C	NUMBER OF OBJECTS IN IT.

        DO J = 1, NCLASS
           IR = Q(J)

C	   EMPTY CLUSTER ? RETURN
           IF (IR .EQ. 0) RETURN

           F = 1. / FLOAT(IR)
c$omp      parallel do if(nfac.ge.64), private(k)
           DO K = 1, NFAC
              S(K,J) = S(K,J) * F
           ENDDO
        ENDDO

C	FIND E(R), THE SUM OF DISTANCES OF EACH OBJECT TO THE 
C	CENTROID OF THE CLUSTER 'R' IT BELONGS TO. 

        DO I = 1, NOBJ
          R = P(I)
          F = 0.

C         COMPUTE THE DISTANCE OF OBJECT 'I' TO THE CENTROID OF 
C         ITS CLUSTER 'R' AND ADD IT TO THE TOTAL SUM OF DISTANCES
C	  'E(R)' IN THAT CLUSTER.

          IF (USE_DISK) THEN
             IF (ITYPE .EQ. 1) THEN
C              UNFORMATTED FILE FOR _SEQ IMAGE DATA
               READ(LUNF) (COO(J), J = 1, NFTOT),FDUM
             ELSE
C              FORMATTED FILE FOR IMAGE OR PIXEL COOR.
               READ(LUNF,*) (COO(J), J = 1, NFTOT),FDUM,FDUM,FDUM,FDUM
             ENDIF
             DO J = 1, NFAC
                X(J,1) = WT(J) * COO(INUM(J))
             ENDDO
             PT = 1
          ELSE
             PT = I
          ENDIF

c$omp     parallel do if(nfac.ge.64), private(k),reduction(+:f)
          DO K = 1, NFAC
             F = F + (S(K,R) - X(K,PT))**2
          ENDDO
          E(R) = E(R) + F
        ENDDO

C	REWIND THE UNIT TO READ THE DATA AGAIN       
        IF (USE_DISK) THEN
           IF (ITYPE .EQ. 1) THEN
C             UNFORMATTED FILE FOR _SEQ IMAGE DATA
              CALL REW(LUNF,1)
           ELSE
C             FORMATTED FILE FOR IMAGE OR PIXEL COOR.
              CALL REWF(LUNF,1)
           ENDIF
        ENDIF

C	COMPUTE THE SUM OF ALL THESE DISTANCES. THE NEW CLUSTERING
C	HAS TO REDUCE THIS D VARIABLE OR STOP AND RETURN.

        D = 0.
        DO J = 1, NCLASS
           D = D + E(J)
        ENDDO

C	FOR EACH OBJECT 'I', CHECK ITS DISTANCE TO THE CENTROID OF
C	THE OTHER CLUSTERS. IF THIS DISTANCE IS SMALLER THAN THE DISTANCE
C	OF THIS OBJECT 'I' TO ITS OWN CLUSTER'S CENTROID, THEN MOVE
C	THE OBJECT 'I' TO THE NEW CLUSTER. STOP WHEN NO MORE OBJECTS 
C	ARE MOVED.

        I  = 0
        IT = 0

C	LOOP UNTIL DONE AND ALL OBJECTS HAVE BEEN CHECKED (I.GT.NOBJ) 
	DONE = .TRUE.

	DO
          I = I + 1
          IF (I .GT. NOBJ) THEN
             IF (DONE) RETURN
             DONE = .TRUE.

C	     REWIND THE UNIT TO READ THE DATA AGAIN       
             IF (USE_DISK) THEN
                IF (ITYPE .EQ. 1) THEN
C                  UNFORMATTED FILE FOR _SEQ IMAGE DATA
                   CALL REW(LUNF,1)
                ELSE
C                  FORMATTED FILE FOR IMAGE OR PIXEL COOR.
                   CALL REWF(LUNF,1)
                ENDIF
              ENDIF

              I = 1
           ENDIF

           IF (USE_DISK) THEN
              IF (ITYPE .EQ. 1) THEN
C                UNFORMATTED FILE FOR _SEQ IMAGE DATA
                 READ(LUNF)   (COO(J), J = 1, NFTOT),FDUM
              ELSE
C                FORMATTED FILE FOR IMAGE OR PIXEL COOR.
                 READ(LUNF,*) (COO(J), J = 1, NFTOT),FDUM,FDUM,FDUM,FDUM
              ENDIF

c$omp         parallel do if(nfac.ge.64), private(j)
              DO J = 1, NFAC
                 X(J,1) = WT(J) * COO(INUM(J))
              ENDDO
              PT = 1
           ELSE
              PT = I
           ENDIF

C	   FIND THE CLUSTER, R, AND THE NUMBER, U, OF OBJECTS IN IT.
C          NOTE THAT IF THERE IS ONLY ONE OBJECT, IT IS CLOSEST TO ITSELF.

           R = P(I)
           U = Q(R)
           IF (U .GT. 1) THEN

C	      COMPUTE DISTANCE 'A' OF OBJECT 'I' TO ITS OWN CLUSTER 'R'

              H = FLOAT(U)
              H = H / (H - 1.)
              F = 0.
c$omp         parallel do if(nfac.ge.64), private(k),reduction(+:f)
              DO K = 1, NFAC
                  F = F + (S(K,R) - X(K,PT))**2
              ENDDO
              A = H * F
              B = 1. E30

C	      FIND DISTANCE 'F' OF THE SAME OBJECT 'I' TO THE OTHER CLUSTERS.
C	      DISTANCE 'B' IS THE SMALLEST OF THESE DISTANCES 'F' CORRESPONDING
C	      TO A CLUSTER 'V' WHICH HAS 'W' NUMBER OF OBJECTS IN IT.

              DO J = 1, NCLASS
                 IF (J .NE. R) THEN
                   U = Q(J)
                   H = FLOAT(U)
                   H = H / (H + 1.)
                   F = 0.
c$omp              parallel do if(nfac.ge.64), private(k),reduction(+:f)
                   DO K = 1, NFAC
                       F = F + (S(K,J) - X(K,PT))**2
                   ENDDO
                   F = H * F
                   IF (F .LE. B) THEN
C                      SMALLER DISTANCE FOUND, SAVE THAT CLUSTER 'V', THE NUMBER 
C	               OF OBJECTS IN IT AND ITS DISTANCE TO 'I'.

                       B = F
                       V = J
                       W = U
                   ENDIF
                 ENDIF          
              ENDDO
 
              IF (B .LT. A) THEN

C	          OBJECT 'I' IS CLOSEST TO CLUSTER 'V' INSTEAD. SO TAKE OBJECT
C	         'I' FROM CLUSTER 'R' AND PUT IT IN CLUSTER 'V'
C
	         DONE = .FALSE.

C                RECOMPUTE THEIR CENTROID AND SUM OF DISTANCES.
C                these are equations (3.2.16, 3.2.17, 3.2.20 and 3.2.21) on page
C	         71-72 of the above book reference.

                 IT   = 0
                 E(R) = E(R) - A
                 E(V) = E(V) + B
                 D    = D - A + B
                 H    = FLOAT(Q(R))
                 G    = FLOAT(W)
                 A    = 1. / (H - 1.)
                 B    = 1. / (G + 1.)

C                RECOMPUTE COORDS. OF CENTROID OF 'R' AND 'V' CLUSTERS

c$omp            parallel do if(nfac.ge.64), private(k,f)
                 DO K = 1, NFAC
                    F = X(K,PT)

C                   OBJECT 'I' REMOVED FROM CLUSTER 'R'.
                    S(K,R) = ( H * S(K,R) - F) * A

C	            OBJECT 'I' ADDED TO CLUSTER 'V'
                    S(K,V) = ( G * S(K,V) + F) * B
                 ENDDO

C	         RESET THE CLUSTER NUMBER FOR OBJECT 'I', RECOMPUTE THE
C	         NUMBER OF ELEMENTS IN CLUSTER 'R' AND 'V'.

                 P(I) = V
                 Q(R) = Q(R) - 1
                 Q(V) = Q(V) + 1

C                PRINT AN OUTPUT
                 IF (IDR .EQ. 1) WRITE(NOUT,16) I, D
16               FORMAT(2X,I4,' MIMINIZE DISTANCE D= ',1X,F12.2,/,
     &           'CLUSTER MEMBERSHIP',/)

                 IF (IDR .EQ. 1) WRITE(NOUT,17) (U,P(U), U = 1, NOBJ)
17               FORMAT(10X,20I3,/)

	      ELSE
C	         OBJECT 'I' IS CLOSEST TO CLUSTER 'R', TAKE NEXT OBJECT.
                 IT = IT + I
              ENDIF
	   ENDIF
        ENDDO
   
        END
