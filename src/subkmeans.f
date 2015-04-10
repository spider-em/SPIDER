C++*********************************************************************
C
C SUBKMEANS.F
C                  OUTPUT COMMENTS                JUN 2009 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@wadsworth.org                                        *
C=*                                                                    *
C=* SPIDER is free software; you can redistribute it and/or            *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* SPIDER is distributed in the hope that it will be useful,          *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* merchantability or fitness for a particular purpose.  See the GNU  *
C=* General Public License for more details.                           *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program. If not, see <http://www.gnu.org/licenses> *
C=*                                                                    *
C **********************************************************************
C
C SUBKMEANS
C
C	W():       WEIGHTS TO BE USED PER OBJECT.
C	INUM():	   FACTORS TO BE USED TO BE USED.
C	NFAC:      NUMBER OF FACTORS TO BE USED
C	NFTOT:     TOTAL NUMBER OF FACTORS USED FOR MULTIVARIATE 
C                  STATISTICAL ANALYSIS (MSA) I.E; CORAN OR PCA.
C	NOBJ:      NUMBER OF OBJECTS.
C       KCLASS:    NUMBER OF CLUSTERS
C       X(M,L):    COORDINATES OF EACH OBJECT
C       IC(M):     CLUSTER TO WHICH OBJECT M BELONGS.
C       CM(L,N):   COORDINATES OF THE CENTROIDS OF EACH CLUSTER.
C       WSS(N):    SUM OF THE SQUARED DISTANCES OF ALL OBJECTS
C                  (IN A CLUSTER) TO THEIR CENTROID.
C       NC(N):     NUMBER OF OBJECTS IN EACH CLUSTER.
C       D:         SUM OF E(N)
C                  (THIS HAS TO BE REDUCED AFTER EACH NEW CLUSTERING)
C       TYPE:      TYPE OF INPUT FILE e.g. _PIX, _IMC, _SEQ
C
C	CALLS:
C	  NEWKMEANS, FILSEQP, FILGET, SAVDN1, SAVD, SAVDC
C
C	NOTE: FILE ON UNIT=LUNF IS ALREADY OPEN.
C             CLOSE THE FILE ON RETURN FROM THIS ROUTINE.
C
C--*********************************************************************

        SUBROUTINE SUBKMEANS(X, CM, NIMI,IC, NC, WSS, AT, NOBJ, NFAC,
     &                       KCLASS, W, NFTOT, COO, INUM, BB, WW, DB,
     &                       IX_DIM,USE_DISK,ISEED,LUNF,LUND,ITYPE)

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 

        REAL                  :: X(NFAC,IX_DIM), CM(NFAC,KCLASS)
        INTEGER               :: IC(NOBJ), NIMI(NOBJ), NC(KCLASS)
        REAL                  :: WSS(KCLASS)
        REAL                  :: AT(NFAC), W(NFAC), COO(NFTOT)
        INTEGER               :: INUM(NFAC)
        REAL                  :: DLIST(3)

        LOGICAL               :: USE_DISK,NEWFILE
        CHARACTER(LEN=MAXNAM) :: FINPAT,FINP,DOCNAM
        INTEGER               :: IDUM(1)

C	SQUARE ROOTS OF THE WEIGHTS.
        DO J = 1, NFAC
           W(J) = SQRT(W(J))
        ENDDO

	IF (ISEED .EQ. 0)  THEN
C	   IC(NOBJ): ASSIGN A CLUSTER TO EACH OBJECT. 
           I = 1
           DO WHILE (I .LE. NOBJ)
              K = 1
              DO WHILE (K .LE. KCLASS)
                 IC(I) = K
                 K     = K + 1
                 I     = I + 1
              ENDDO
           ENDDO
	ELSE
C          RANDOM ASSIGNMENT
	   T = 1.0/FLOAT(KCLASS)
	   I = 1
	   DO WHILE (I .LE. NOBJ)
	      Q  = RAND_P(ISEED)
	      QT = T
	      K  = 1
	      DO WHILE (QT .LE. Q)
	         QT = QT +T 
	         K  = K + 1
	      ENDDO
	      IC(I) = K
	      I     = I + 1
	   ENDDO
	ENDIF

        IDR = 0
        CALL NEWKMEANS(NOBJ,NFAC,X,IC,KCLASS,CM,WSS,NC,D,IDR,
     &               INUM,NIMI,W,COO,NFTOT,IX_DIM,USE_DISK,LUNF,ITYPE)

C       CALCULATE CITERIA ...

c$omp   parallel do if(nfac.ge.64), private(i)
	DO I = 1,NFAC
           AT(I) = 0.0
        ENDDO

        DO L = 1,NOBJ
c$omp     parallel do if(nfac.ge.64), private(i)
          DO I = 1,NFAC
            AT(I) = AT(I) + X(I,L)
          ENDDO
        ENDDO

c$omp   parallel do if(nfac.ge.64), private(i)
        DO I = 1,NFAC
          AT(I) = AT(I)/NOBJ
        ENDDO

        BB = 0.0
c$omp   parallel do if(nfac.ge.64), private(i,j),reduction(+:bb)
        DO I = 1,NFAC
           DO J = 1,KCLASS
              BB = BB + (CM(I,J) - AT(I))**2
           ENDDO
        ENDDO

C	REWIND THE UNIT=LUNF

        IF (USE_DISK) THEN
           IF (ITYPE .EQ. 1) THEN
C             UNFORMATTED FILE FOR _SEQ IMAGE DATA
              CALL REW(LUNF,1)
           ELSE
C             FORMATTED FILE FOR IMAGE OR PIXEL COOR.
              CALL REWF(LUNF,1)
           ENDIF
        ENDIF

	DO I=1,KCLASS
 	   WSS(I) = 0.0
	ENDDO

        DO L = 1, NOBJ
          IF (USE_DISK) THEN
             IF (ITYPE .EQ. 1) THEN
C               UNFORMATTED FILE FOR _SEQ IMAGE DATA
                READ(LUNF)   (COO(J), J = 1, NFTOT),FDUM
             ELSE
C               FORMATTED FILE FOR IMAGE OR PIXEL COOR.
                READ(LUNF,*) (COO(J), J = 1, NFTOT),FDUM,FDUM,FDUM,FDUM
             ENDIF

             DO J = 1, NFAC
                X(J,IX_DIM) = W(J) * COO(INUM(J))
             ENDDO
             K = IX_DIM
          ELSE
             K = L
          ENDIF

          WW = 0.0
c$omp     parallel do if(nfac.ge.64), private(i),reduction(+:ww)
          DO I = 1, NFAC
              WW = WW + (X(I,K) - CM(I,IC(L)))**2
          ENDDO
	  WSS(IC(L)) = WSS(IC(L))+WW
        ENDDO
        WW = 0.0
	DO I=1,KCLASS
 	   WW = WW + WSS(I)
	ENDDO

C       Davis-Bouldin criterium
	DB = 0.0
	DO J=1,KCLASS
	   RK = -1.0
	   DO I=1,KCLASS
	      IF (I.NE.J)  THEN
	         SMK = 0.0
c$omp            parallel do if(nfac.ge.64), private(l),reduction(+:smk)
	         DO L=1,NFAC
 	            SMK = SMK + (CM(L,J) - CM(L,I))**2
	         ENDDO
c$omp            end parallel do 
	         RJK = (WSS(J) + WSS(I) ) / SMK
	         RK  = AMAX1(RK,RJK)
	      ENDIF
	   ENDDO
 	   DB = DB + RK
	ENDDO

	DB = DB / K

C	WRITE THE OBJECTS IN EACH CLASS TO A FILE.
        NMAX = 0
        CALL FILSEQP(FINPAT,NLET,IDUM,NMAX,NIMA,
     &        'SELECTION DOC FILE TEMPLATE (e.g.: SEL***)',IRTFLG)

C	LOOK THROUGH IC(), PUT ALL OBJECTS BELONGING TO CLASS IN FILE. 

        DO ICLASS = 1, KCLASS

C	  CREATE A FILENAME FOR CLASS ICLASS
          CALL FILGET(FINPAT,FINP,NLET,ICLASS,IRTFLG)

C	  NUMBER OF OBJECT IN THIS GROUP
          IGRP = 0
          NLS  = 2
	  DO I = 1, NOBJ
            IF (ICLASS .EQ. IC(I)) THEN

C             THE OBJECT BELONG TO CLASS ICLASS, STORE IT.
              IGRP     = IGRP + 1
              DLIST(1) = IGRP
              DLIST(2) = NIMI(I)
              IAP      = 0
              CALL SAVDN1(LUND,FINP,DLIST,NLS,IGRP-1,IAP)
            ENDIF
	  ENDDO

C         CLOSE EACH CLASS FILE OPENED BY SAVDN1 ROUTINE.
          CLOSE(LUND)
	ENDDO

        CALL  SAVDC

C       OPEN OUTPUT DOC FILE NAME
        CALL OPENDOC(DOCNAM,.TRUE.,NLET,LUND,NICDOCO,.TRUE.,
     &         'CLASS MEMBERSHIP DOC FILE',
     &         .FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999

        CALL LUNDOCPUTCOM(NICDOCO,
     &                    'KEY, IMAGE/ELEMENT,    CLASS ',IRTFLG)

        DO I = 1, NOBJ
           DLIST(1) = NIMI(I)
           DLIST(2) = IC(I)
           CALL LUNDOCWRTDAT(NICDOCO,I, DLIST,2,IRTFLG)
        ENDDO

9999    CLOSE(LUND)
        CLOSE(LUNF)

        END
