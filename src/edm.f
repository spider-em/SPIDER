C ++********************************************************************
C                                                                      *
C EDM                     CREATED MAY 01 ARDEAN LEITH                  * 
C                         REMOVED AUTO ARRAYS (TOO BIG) APR 02 A. LEITH
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
C  EDM(LUN1,LUN2,NSAM,NROW,NSLICE,FMINT,FMAXT)
C
C  PARAMETERS:
C
C  PURPOSE: CREATE EUCLIDEAN DISTANCE MAP FROM AN IMAGE OR VOLUME 
C           SHOWING DISTANCE FROM A PIXEL/VOXEL TO NEAREST BACKGROUND
C           PIXEL/VOLUME
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C **********************************************************************

	SUBROUTINE EDM(LUN1,LUN2,NSAM,NROW,NSLICE,FMINT,FMAXT)

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	REAL, ALLOCATABLE, DIMENSION(:) ::  VIN
	REAL, ALLOCATABLE, DIMENSION(:) ::  VOUT

        THRESH = FMINT + (FMAXT - FMINT) / 2.0
10      CALL RDPRM1S(THRESH,NOT_USED,'BACKGROUND THRESHOLD',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (THRESH .GT. FMAXT) THEN
           CALL ERRT(101,'THRESHOLD MUST BE < THAN FMAX',IDUM) 
           GOTO 10
        ENDIF

        ALLOCATE(VIN(NSAM*NROW*NSLICE),
     &           VOUT(NSAM*NROW*NSLICE),STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'ER EDM, VIN,VOUT',IER)
            RETURN
        ENDIF

C       LOAD INPUT IMAGE/VOLUME
        CALL REDVOL(LUN1,NSAM,NROW,1,NSLICE,VIN,IRTFLG)
	IF (IRTFLG .NE. 0) GOTO 999

        IF (NSLICE .GT. 1) THEN
           CALL EDM3(VIN,VOUT,NSAM,NROW,NSLICE,THRESH,LUN2,
     &               NABOVE,LAYERS)

           WRITE(NOUT,*) ' VOXELS ABOVE THRESHOLD: ',NABOVE
           WRITE(NOUT,*) ' LAYERS: ',LAYERS

        ELSE
    
           CALL EDM2(VIN,VOUT,NSAM,NROW,NSLICE,THRESH,LUN2,
     &              NABOVE,LAYERS)
           WRITE(NOUT,*) ' PIXELS ABOVE THRESHOLD: ',NABOVE
           WRITE(NOUT,*) ' LAYERS: ',LAYERS
        ENDIF

999     IF (ALLOCATED(VIN))  DEALLOCATE(VIN)
	IF (ALLOCATED(VOUT)) DEALLOCATE(VOUT)

        END

C       ------------------------- EDM2 --------------------------------

	SUBROUTINE EDM2(X,X2,NSAM,NROW,NSLICE,THRESH,LUN2,NABOVE,LAYERS)

	REAL, DIMENSION(NSAM,NROW) :: X
	REAL, DIMENSION(NSAM,NROW) :: X2
	REAL, DIMENSION(-1:1,-1:1) :: DISTA

C       AUTOMATIC ARRAYS MAY BE TOO BIG, DO NOT USE

        DATA DISTA / 1.41, 1.0,  1.41,   
     &               1.00, 0.0,  1.00,  
     &               1.41, 1.0,  1.41/

C          INITIALIZE EDM
           NABOVE = 0
           LAYERS = 0

           DO J=1,NROW                      
              DO I=1,NSAM
                 IF (X(I,J) .LT. THRESH) THEN
                    X2(I,J) = 0.0
                    X (I,J) = 0.0
                 ELSE
                    X2(I,J) = -1.0
                    X (I,J) = -1.0
                    NABOVE  = NABOVE + 1
                 ENDIF
              ENDDO
           ENDDO 

10         NPIXER = 0
           LAYERS = LAYERS + 1

           DO J=1,NROW                      
              DO I=1,NSAM
                 IF (X(I,J) .LT. 0.0)  THEN
C                   CENTRAL PIXEL IS ABOVE THRESHOLD, FIND DISTANCE
                    NPIXER = NPIXER + 1

                    DISMIN = HUGE(DISMIN)
                    DO MJ=-1,1
                       MJM = MOD(J+MJ+NROW-1,NROW)+1
                       DO MI=-1,1
                          VALT = X(MOD(I+MI+NSAM-1,NSAM)+1,MJM)
                          IF (VALT .GE. 0.0) THEN
                             DISMIN = MIN(DISMIN,VALT + DISTA(MI,MJ))
                          ENDIF
                       ENDDO
C                      END OF: DO MI=-1,1

                    ENDDO
C                   END OF: DO MJ=-1,1
                    
                    IF (DISMIN .LT. HUGE(DISMIN)) X2(I,J) = DISMIN
	         ENDIF
C                END OF CLAUSE: IF (X(I,J) .LT. 0.0) 

              ENDDO
C             END OF: DO I=1,NSAM 

           ENDDO
C          END OF: DO J=1,NROW 

C          SEE IF WE NEED MORE DISTANCE CYCLES
           IF (NPIXER .GT. 0) THEN
              X = X2 
              GOTO 10
           ENDIF

C       OUTPUT IMAGE
        CALL WRTVOL(LUN2,NSAM,NROW,1,NSLICE,X2,IRTFLG)

        END	

C       ------------------------- EDM3 --------------------------------

	SUBROUTINE EDM3(X,X2,NSAM,NROW,NSLICE,THRESH,LUN2,NABOVE,LAYERS)

	REAL, DIMENSION(NSAM,NROW,NSLICE) :: X
	REAL, DIMENSION(NSAM,NROW,NSLICE) :: X2
	REAL, DIMENSION(-1:1,-1:1,-1:1) :: DISTA

C       AUTOMATIC ARRAYS MAY BE TOO BIG, DO NOT USE

        DATA DISTA /
     &               1.73, 1.41,  1.73,   
     &               1.00, 1.00,  1.00,  
     &               1.73, 1.41,  1.73,

     &               1.41, 1.00,  1.41,   
     &               1.00, 0.00,  1.00,  
     &               1.41, 1.00,  1.41,

     &               1.73, 1.41,  1.73,   
     &               1.00, 1.00,  1.00,  
     &               1.73, 1.41,  1.73 /

C       INITIALIZE EDM
        NABOVE = 0
        LAYERS = 0

        DO K=1,NSLICE
           DO J=1,NROW                      
              DO I=1,NSAM
                 IF (X(I,J,K) .LT. THRESH) THEN
                    X(I,J,K)  = 0.0
                    X2(I,J,K) = 0.0
                 ELSE
                    X(I,J,K)  = -1.0
                    X2(I,J,K) = -1.0
                    NABOVE   = NABOVE + 1
                 ENDIF
              ENDDO
           ENDDO
        ENDDO 

10      NPIXER = 0
        LAYERS = LAYERS + 1

        DO K=1,NSLICE
           DO J=1,NROW                      
              DO I=1,NSAM
                 IF (X(I,J,K) .LT. 0.0)  THEN
C                   CENTRAL VOXEL IS ABOVE THRESHOLD, FIND DISTANCE

                    DISMIN = HUGE(DISMIN)
                    NPIXER = NPIXER + 1

                    DO MK=-1,1
                    MKM = MOD(K+MK+NSLICE-1,NSLICE)+1

                    DO MJ=-1,1
                       MJM = MOD(J+MJ+NROW-1,NROW)+1
                       DO MI=-1,1
                          VALT = X(MOD(I+MI+NSAM-1,NSAM)+1,MJM,MKM)
                          IF (VALT .GE. 0.0) THEN
                             DISMIN = MIN(DISMIN,VALT + DISTA(MI,MJ,MK))
                          ENDIF
                       ENDDO
C                      END OF: DO MI=-1,1

                    ENDDO
C                   END OF: DO MJ=-1,1
                    ENDDO
C                   END OF: DO MK=-1,1
                    
                    IF (DISMIN .LT. HUGE(DISMIN)) X2(I,J,K) = DISMIN
	         ENDIF
C                END OF CLAUSE: IF (X(I,J,K) .LT. 0.0) 

              ENDDO
C             END OF: DO I=1,NSAM 

           ENDDO
C          END OF: DO J=1,NROW 

        ENDDO
C       END OF: DO K=1,NSLICE 


C       SEE IF WE NEED MORE DISTANCE CYCLES
        IF (NPIXER .GT. 0) THEN
           X = X2
           GOTO 10
        ENDIF

C       OUTPUT VOLUME
        CALL WRTVOL(LUN2,NSAM,NROW,1,NSLICE,X2,IRTFLG)

        END	
