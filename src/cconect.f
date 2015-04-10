
C **********************************************************************
C
C   CCONECT.FOR  -- CREATED OCT 90
C **********************************************************************
C *  AUTHOR: ArDean Leith 
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
C      CCONECT(NSAM,NROW,LUNOUT,SLICE1,SLICE2,LASTSLI,
C              IEQUIV,NEQUIV,NEQMAX,LASTCLUS,IRTFLG)
C
C      PURPOSE:  DETERMINES 3-D CONNECTIVITY USING 2 SLICES AT A TIME    
C
C      PARAMETERS:  
C
C      CALLED BY:  CONINT  
C
C      CALLS:       
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--********************************************************************

       SUBROUTINE CCONECT(NSAM,NROW,LUNOUT,SLICE1,SLICE2,LASTSLI,
     &                   IEQUIV,NEQUIV,NEQMAX,LASTCLUS,IRTFLG)

       INTEGER * 2  :: SLICE1(*),SLICE2(*)
       INTEGER      :: IEQUIV(2,NEQMAX)
       LOGICAL      :: INCLUS,LASTSLI

       IRTFLG = 0

       INCLUS = .FALSE.
       DO  IROW = 1, NROW
          IPTR0 = (IROW-1) * NSAM

       DO  ICOL = 1, NSAM
          IPTR1 = IPTR0 + ICOL
          NOW   = SLICE1(IPTR1)

          IF (NOW .EQ. 0) THEN
C            EMPTY VOXEL
             INCLUS = .FALSE.

          ELSEIF (NOW .GT. 0) THEN
C            HAVE ALREADY VISITED THIS VOXEL FROM LAST ROW OR SLICE

             IF (.NOT. INCLUS) THEN
C               PREVIOUSLY IN BACKGROUND REGION ON THIS ROW
                NCLUS  = NOW
                INCLUS = .TRUE.
                LASCON = 0

             ELSEIF (INCLUS .AND. NOW .NE. NCLUS .AND. 
     &                            NOW .NE. LASCON) THEN
C               Y SHAPED CLUSTER BRANCH FROM PREVIOUS SLICE
                CALL GOTBRANCH(NCLUS,NOW,IEQUIV,NEQUIV,NEQMAX,IRTFLG)
                IF (IRTFLG .NE. 0) RETURN
             ENDIF
             SLICE1(IPTR1) = NCLUS

          ELSEIF (NOW .LT. 0) THEN
C            VOXEL OCCUPIED AND NOT VISITED YET

             IF (.NOT. INCLUS) THEN
C               PREVIOUSLY IN A BACKGROUND AREA ON THIS LINE
                NCLUS    = LASTCLUS + 1
                LASTCLUS = NCLUS
                LASCON   = 0
                INCLUS   = .TRUE.
             ENDIF

C            OCCUPY THE VOXEL WITH THE CLUSTER NUMBER
             SLICE1(IPTR1) = NCLUS

          ENDIF

          IF (INCLUS) THEN
C            HAVE AN OCCUPIED VOXEL ON THIS ROW

             IF (IROW .LT. NROW) THEN

C              CHECK THE NEXT ROW DOWN FROM THIS OCCUPIED VOXEL
               IPTR2 = IPTR1 + NSAM
               NOW2  = SLICE1(IPTR2)
               
               IF (NOW2 .NE. 0) THEN
C                 VOXEL ON NEXT ROW DOWN IS ALSO IN THIS CLUSTER

                 IF (NOW2 .GT. 0 .AND. NOW2 .NE. NCLUS) THEN
C                  HAVE A VISITED VOXEL FROM A BRANCH CONNECTION
                   CALL GOTBRANCH(NCLUS,NOW2,IEQUIV,NEQUIV,
     &                            NEQMAX,IRTFLG)
                   IF (IRTFLG .NE. 0) RETURN
                 ENDIF
 
C                OCCUPY THE VOXEL ON THIS SECOND ROW
                 SLICE1(IPTR2) = NCLUS

                 IF (ICOL .GT. 1) THEN
C                  GO BACKWARDS ALONG SECOND ROW OCCUPYING THIS CLUSTER
                   DO  IPTRB = IPTR2-1,IPTR0+NSAM+1,-1
                     NOWB = SLICE1(IPTRB)

                     IF (NOWB .LT. 0) THEN
C                      THIS IS A CLUSTER VOXEL
                       SLICE1(IPTRB) = NCLUS

                     ELSE
C                      END OF CLUSTER ON SECOND ROW LEFT
                       GOTO 12
                     ENDIF
		   ENDDO
c10                 CONTINUE 
                 ENDIF
               ENDIF

12             CONTINUE
             ENDIF

             IF (.NOT. LASTSLI) THEN
               NOWD = SLICE2(IPTR1)
               IF (NOWD .LT. 0) THEN
C                 UNVISITED VOXEL FROM THIS CLUSTER ON NEXT SLICE DOWN
                  SLICE2(IPTR1) = NCLUS

                 IF (ICOL .GT. 1) THEN
C                  GO BACKWARDS ALONG 2ND SLICE ROW OCCUPYING THIS CLUSTER
                   DO  IPTRDB = IPTR1-1,IPTR0+1,-1
                     NOWDB = SLICE2(IPTRDB)

                     IF (NOWDB .LT. 0) THEN
C                      THIS IS A CLUSTER VOXEL
                       SLICE2(IPTRDB) = NCLUS
                     ELSE
C                      END OF CLUSTER ON THIS ROW OF SECOND SLICE
                       GOTO 32
                     ENDIF
                   END DO
32                 CONTINUE
                 ENDIF

               ENDIF
             ENDIF
          ENDIF

21        CONTINUE

	ENDDO
	ENDDO
       END

