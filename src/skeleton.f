
C ++********************************************************************
C                                                                      *
C SKELETON                             CREATED APR 23 2001 ARDEAN LEITH                  * 
C                                                                      *
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
C  SKELETON(LUN1,LUN2,NSAM,NROW,NSLICE)
C
C  PARAMETERS:
C
C  PURPOSE: ERODE (SHRINK) AN OBJECT IN AN IMAGE OR VOLUME UNTIL A
C           SKELETON IS REACHED USING A CROSS (4 FOLD CONNECTIVITY)
C           CONVENTION
C                                                                     *
C **********************************************************************

	SUBROUTINE SKELETON(LUN1,LUN2,NSAM,NROW,NSLICE)

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	REAL, ALLOCATABLE, DIMENSION(:,: )  :: VIN,VOUT

	CHARACTER(LEN=1) ::   MODE 

        IF (NSLICE .GT. 1) THEN
           CALL ERRT(101,'VOLUME SKELETONS NOT IMPLEMENTED',IER)
           RETURN
        ENDIF

cc 	CALL RDPRMC(MODE,NA,.TRUE.,'BOX OR CROSS NEIGHBORHOOD (B/C)',
cc     &               CHAR(0),IRTFLG)
cc        IF (IRTFLG .NE. 0) RETURN

        MODE = 'C'

        LENGTH = 3
        IF (MODE .EQ. 'B') THEN
           NEIGH = LENGTH * LENGTH - 1
           CALL ERRT(101,
     &     'BOX NEIGHBORHOOD (8 FOLD CONNECTIVITY) NOT IMPLEMENTED YET'
     &     ,IER)
           RETURN
        ELSE
           NEIGH = 2 * LENGTH - 2 
        ENDIF

        WRITE(NOUT,91) NEIGH
91      FORMAT(' NUMBER OF NEIGHBORS:',I5)

        ALLOCATE(VIN(NSAM,NROW),VOUT(NSAM,NROW),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'ER SK --VOUT, VIN',IER)
           RETURN
        ENDIF

C       LOAD INPUT IMAGE
        DO J = 1,NROW
           CALL REDLIN(LUN1,VIN(1,J),NSAM,J)
           DO I = 1,NSAM
              VOUT(I,J) = VIN(I,J)
           ENDDO
        ENDDO
   
        CALL SKELETON2(VIN,VOUT,NSAM,NROW,NSLICE,NEIGH, MODE,LUN2,
     &                   NPIXER)

        WRITE(NOUT,*) ' PIXELS ERODED: ',NPIXER

        DEALLOCATE(VIN,VOUT)

        END

C       ------------------------- SKELETON2 -----------------------------


	SUBROUTINE SKELETON2(VIN,VOUT,NSAM,NROW,NSLICE,NEIGH,MODE,
     &                      LUN2,NPIXERTOT)

	REAL, DIMENSION(NSAM,NROW) :: VIN,VOUT
        INTEGER, DIMENSION(8) ::      IXSQ,IYSQ
        INTEGER, DIMENSION(4) ::      IXCR,IYCR
	CHARACTER(LEN=1) ::           MODE 
         
        DATA IXSQ/-1, 0, 1,  1, 1,0,-1, -1/
        DATA IYSQ/-1,-1,-1,  0, 1,1,1,   0/

        DATA IXCR/ 0,  1,  0, -1/
        DATA IYCR/-1,  0,  1,  0/

        NPIXERTOT = 0

10      NPIXER    = 0

        DO J=1,NROW                      
           DO I=1,NSAM
              VOUT(I,J) = VIN(I,J)

              IF (VIN(I,J) .GT. 0.0)  THEN
C                CURRENT PIXEL IS OCCUPIED (NON-ZERO)

                 IF (MODE .EQ. 'C') THEN
C                   FIND NUMBER OF ORIGINAL OCCUPIED NEIGHBORS 

                    IHITS = 0
                    DO II = 1,4
                       IX = MOD(I+IXCR(II)+NSAM,NSAM) 
                       IY = MOD(J+IYCR(II)+NROW,NROW)
                       IF (VIN(IX,IY) .GT. 0.0) IHITS = IHITS + 1
                    ENDDO

                    IF (IHITS .GT. 1 .AND. IHITS .LT. 4) THEN
C                      THIS PIXEL COULD BE ERODED ON THIS CYCLE
C                      FIND NUMBER OF CURRENT OCCUPIED NEIGHBORS 

                       IHITS = 0
                       DO II = 1,4
                          IX = MOD(I+IXCR(II)+NSAM,NSAM) 
                          IY = MOD(J+IYCR(II)+NROW,NROW)
                          IF (VOUT(IX,IY) .GT. 0.0) IHITS = IHITS + 1
                       ENDDO

C                      INITIALIZE CLUSTER BREAKING TEST
                       ITRANS = 0
                       IX     = MOD(I+IXSQ(8)+NSAM,NSAM) 
                       IY     = MOD(J+IYSQ(8)+NROW,NROW)
                       IOLD   = VOUT(IX,IY)
                       IF (IOLD .GT. 0) IOLD = 1

C                      FIND TRANSITIONS (USED TO PREVENT CLUSTER BREAKS)
                       DO II = 1,8
                          IX    = MOD(I+IXSQ(II)+NSAM,NSAM) 
                          IY    = MOD(J+IYSQ(II)+NROW,NROW)
                          INEW  = VOUT(IX,IY)
                          IF (INEW .GT. 0) INEW = 1
                          IF (IOLD .NE. INEW) ITRANS = ITRANS + 1
                          IOLD = INEW
                       ENDDO

                       IF (IHITS .GT. 1 .AND. ITRANS .LE. 2) THEN
C                         ERODE CURRENT PIXEL (SET IT TO ZERO) IF
C                         IT DOESN'T DISCONNECT A CLUSTER.
C                         MORE THAN TWO TRANSISTIONS WOULD BREAK CLUSTER
                          NPIXER = NPIXER + 1
CC                        write(6,90) i,j,ihits,itrans
90                        format(' Eroded(',i3,',',i3,'): ',i2,
     &                           ' trans: ',i2)
                          VOUT(I,J) = 0.0
                       ENDIF
                    ENDIF

                 ELSE
C                   "BOX" OR 8 FOLD CONNECTIVITY (NOT IMPLEMENTED !!!!)
C                   FIND NUMBER OF ORIGINAL OCCUPIED NEIGHBORS 

                    IHITS = 0
                    DO II = 1,8
                       IX = MOD(I+IXSQ(II)+NSAM,NSAM) 
                       IY = MOD(J+IYSQ(II)+NROW,NROW)
                       IF (VIN(IX,IY) .GT. 0.0) IHITS = IHITS + 1
                    ENDDO

                    IF (IHITS .GT. 1 .AND. IHITS .LT. 8) THEN
C                      THIS PIXEL COULD BE ERODED ON THIS CYCLE
C                      FIND NUMBER OF CURRENT OCCUPIED NEIGHBORS 

                       IHITS = 0
                       DO II = 1,8
                          IX = MOD(I+IXSQ(II)+NSAM,NSAM) 
                          IY = MOD(J+IYSQ(II)+NROW,NROW)
                          IF (VOUT(IX,IY) .GT. 0.0) IHITS = IHITS + 1
                       ENDDO

C                      NEED FATE MAP IMPLEMENTED FOR 256 COMBINATIONS
C                      (I DO NOT HAVE TIME NOW FOR THIS) !!!!!!!

                       IF (IHITS .GT. 1 .AND. ITRANS .LE. 2) THEN
C                         ERODE CURRENT PIXEL (SET IT TO ZERO) IF
C                         IT DOESN'T DISCONNECT A CLUSTER.
C                         MORE THAN TWO TRANSISTIONS WOULD BREAK CLUSTER
                          NPIXER = NPIXER + 1
                          VOUT(I,J) = 0.0
                       ENDIF
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
        ENDDO

        IF (NPIXER .GT. 0) THEN
C          NEED MORE ERODE CYCLES
           NPIXERTOT = NPIXERTOT + NPIXER

           DO I = 1,NROW
              DO J = 1,NSAM
                 VIN(I,J) = VOUT(I,J)
              ENDDO
           ENDDO

           GOTO 10
        ENDIF
  
C       OUTPUT IMAGE
        DO J = 1,NROW
           CALL WRTLIN(LUN2,VOUT(1,J),NSAM,J)
        ENDDO
      END	

