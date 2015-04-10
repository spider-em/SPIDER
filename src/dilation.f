C ++********************************************************************
C                                                                      *
C DILATION                      ADDED 'NOFUSE" MAR 01 ARDEAN LEITH                                       *
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
C  DILATION
C  
C  PURPOSE:  DILATES  A BINARY IMAGE/VOLUME.  WHENEVER THE NUMBER OF
C            OCCUPIED NEIGHBORING PIXELS EXCEEDS A SPECIFIED LIMIT
C            A CENTRAL PIXEL IS ASSIGNED                                                                     *
C **********************************************************************

	SUBROUTINE DILATION(LUN1,LUN2,NSAM,NROW,NSLICE)

	INCLUDE 'CMBLOCK.INC'
	REAL, ALLOCATABLE, DIMENSION(:,:)   ::  VIN2
	REAL, ALLOCATABLE, DIMENSION(:,:,:) ::  VIN3
	CHARACTER (LEN=1) ::                    MODE,NULL
        LOGICAL ::                              NOFUSE

	NULL = CHAR(0)

C       FCHAR OPTION 'NF' DOESNT FUSE CLUSTERS
        NOFUSE = (FCHAR(4:5) .EQ. 'NF')

	CALL RDPRMC(MODE,NA,.TRUE.,'BOX OR CROSS NEIGHBORHOOD (B/C)',
     &               NULL,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

10      CALL RDPRI1S(LENGTH,NOT_USED,'LENGTH OF NEIGHBORHOOD',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (LENGTH .LE. 1) THEN
           CALL ERRT(102,'LENGTH MUST BE GREATER THAN',2) 
           GOTO 10
        ENDIF

        IF (MOD(LENGTH,2) .EQ. 0) THEN
           LENGTH = LENGTH + 1
           WRITE(NOUT,90) LENGTH 
90         FORMAT(' EFFECTIVE LENGTH OF NEIGHBORHOOD:',I5)
        ENDIF
        
        IF (NSLICE .LE. 1) THEN
	   IF (MODE .EQ. 'B') THEN
              LAT = LENGTH * LENGTH - 1
	   ELSE
              LAT = 2 * LENGTH - 2 
	   ENDIF
        ELSE
	   IF (MODE .EQ. 'B') THEN
              LAT = LENGTH * LENGTH * LENGTH - 1
	   ELSE
              LAT = 3 * LENGTH - 3 
	   ENDIF
	ENDIF

        WRITE(NOUT,91) LAT
91      FORMAT(' NUMBER OF NEIGHBORS:',I5)

	CALL RDPRI1S(LAT,NOT_USED,
     &     'DILATE IF NUMBER OF OCCUPIED NEIGHBORS IS > THAN',IRTFLG)

        IF (NSLICE .GT. 1) THEN
	   ALLOCATE(VIN3(NSAM,NROW,NSLICE),STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
               CALL ERRT(46,'DI, VIN3',IER)
               RETURN
           ENDIF

C          LOAD INPUT VOLUME
	   INDEX = 0
           DO K = 1,NSLICE
              DO J = 1,NROW
                 INDEX = INDEX + 1
                 CALL REDLIN(LUN1,VIN3(1,J,K),NSAM,INDEX)
	      ENDDO
	   ENDDO
     
           LH = LENGTH / 2 
           CALL DILATION3(VIN3,NSAM,NROW,NSLICE,LH,LAT,MODE,NOFUSE,LUN2)

	   DEALLOCATE(VIN3)
        ELSE
	   ALLOCATE(VIN2(NSAM,NROW),STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
               CALL ERRT(46,'DI, VIN2',IER)
               RETURN
           ENDIF

C          LOAD INPUT VOLUME
	   INDEX = 0
           DO J = 1,NROW
              INDEX = INDEX + 1
              CALL REDLIN(LUN1,VIN2(1,J),NSAM,INDEX)
	   ENDDO
     
           LH = LENGTH / 2 
           CALL DILATION2(VIN2,NSAM,NROW,NSLICE,LH,LAT,MODE,NOFUSE,LUN2)

	   DEALLOCATE(VIN2)
        ENDIF
        END
        

C       -------------------- DILATION2-------------------------------------

	SUBROUTINE DILATION2(X,NSAM,NROW,NSLICE,LH,LAT,MODE,NOFUSE,LUN2)

C       2D DILATION

	DIMENSION    X(NSAM,NROW),Y(NSAM)
	CHARACTER(LEN=1) ::   MODE
        LOGICAL ::            NOFUSE

           DO J=1,NROW                      
              IF (MODE .EQ. 'C')  THEN
C                "CROSS"
                 DO I=1,NSAM
C                   COPY UNDILATED VOXEL VALUE
                    Y(I) = X(I,J)

                    IF (X(I,J) .LE. 0.0) THEN
C                      CENTRAL VOXEL IS UNOCCUPIED (ZERO)
                       LB     = 0
                       VALGOT = 1.0
                       DO M=-LH,LH
                          IF (M. NE. 0)  THEN
C                            COUNT NUMBER OF OCCUPIED NEIGHBORS
                             CALL DILCHK(X(I,MOD(J+M+NROW-1,NROW)+1),
     &                                   VALGOT,LB,NOFUSE)

                             CALL DILCHK(X(MOD(I+M+NSAM-1,NSAM)+1,J),
     &                                   VALGOT,LB,NOFUSE)
                             IF (LB .LT. 0) EXIT
                          ENDIF
                       ENDDO

                       IF (LB .GT. LAT) THEN
C                         SET VOXEL TO OCCUPIED (DILATE) 
                          Y(I) = VALGOT
                      ENDIF
                   ENDIF
                 ENDDO
              ELSE
C                "BOX"
                 DO I=1,NSAM
C                   COPY UNDILATED VOXEL VALUE
                    Y(I) = X(I,J)

                    IF (X(I,J) .LE. 0.0)  THEN
C                      CENTRAL VOXEL IS NOT OCCUPIED (ZERO)
                       LB     = 0
                       VALGOT = 1.0

C                      COUNT NUMBER OF OCCUPIED NEIGHBORS
                       DO MJ=-LH,LH
                          MJM = MOD(J+MJ+NROW-1,NROW)+1
                          DO MI=-LH,LH
                             CALL DILCHK(X(MOD(I+MI+NSAM-1,NSAM)+1,MJM),
     &                                  VALGOT,LB,NOFUSE)
                             IF (LB .LT. 0) EXIT
                          ENDDO
                       ENDDO

                       IF (LB .GT. LAT) THEN
C                         SET VOXEL TO FILLED 
                          Y(I) = VALGOT
                       ENDIF
                     ENDIF
                 ENDDO
	      ENDIF

C             OUTPUT VOLUME
              CALL WRTLIN(LUN2,Y,NSAM,J)
           ENDDO
      END	


C       -------------------- DILCHK-------------------------------------

	SUBROUTINE DILCHK(VAL,VALGOT,LB,NOFUSE)

C       CHECKS FOR FUSION

        LOGICAL :: NOFUSE
        REAL    :: VALGOT

        IF (VAL .GT. 0.0) THEN
C          OCCUPIED NEIGHBOR PRESENT HERE
 
           IF (NOFUSE) THEN
C             DO NOT WANT TO FUSE CLUSTERS HAVING DIFFERENT VALUES

              IF (LB .GT. 0 .AND. VAL .NE. VALGOT) THEN
C                DILATION WOULD FUSE DIFFERENT CLUSTERS
C                MAKE LB SO NEGATIVE THAT IT WILL NEVER EXCEED LAT
                 LB = -100000
                 RETURN
              ENDIF
              VALGOT = VAL
           ENDIF

C          INCREMENT NUMBER OF OCCUPIED NEIGHBORS
           LB = LB + 1
        ENDIF

        RETURN
        END

C       -------------------- DILATION3-------------------------------------

	SUBROUTINE DILATION3(X,NSAM,NROW,NSLICE,LH,LAT,MODE,NOFUSE,LUN2)

C       3D DILATION

	DIMENSION    X(NSAM,NROW,NSLICE),Y(NSAM)
	CHARACTER(LEN=1) :: MODE
        LOGICAL ::          NOFUSE

        DO N=1,NSLICE                    
           DO J=1,NROW                      
              IF (MODE .EQ. 'C')  THEN
C                "CROSS"
                 DO I=1,NSAM
C                   COPY UNDILATED VOXEL VALUE
                    Y(I) = X(I,J,N)

                    IF (X(I,J,N) .LE. 0.0)  THEN
C                      CENTRAL VOXEL IS UNOCCUPIED (ZERO) MT
                       LB     = 0
                       VALGOT = 1.0

                       DO M=-LH,LH
                          IF (M. NE. 0)  THEN
C                            COUNT NUMBER OF OCCUPIED NEIGHBORS 
                             CALL DILCHK(X(I,MOD(J+M+NROW-1,NROW)+1,N),
     &                                  VALGOT,LB,NOFUSE)

                             CALL DILCHK(X(MOD(I+M+NSAM-1,NSAM)+1,J,N),
     &                                  VALGOT,LB,NOFUSE)

                             CALL DILCHK
     &                               (X(I,J,MOD(N+M+NSLICE-1,NSLICE)+1),
     &                               VALGOT,LB,NOFUSE)
                          ENDIF
                       ENDDO
                       IF (LB .GT. LAT) THEN
C                         SET VOXEL TO OCCUPIED (DILATE) 
                          Y(I) = VALGOT
                      ENDIF
                   ENDIF
                 ENDDO
              ELSE
C                "BAR"
                 DO I=1,NSAM
C                   COPY UNDILATED VOXEL VALUE
                    Y(I) = X(I,J,N)

                    IF (X(I,J,N) .LE. 0.0)  THEN
C                      VOXEL IS NOT OCCUPIED (ZERO)
                       LB     = 0
                       VALGOT = 1.0

C                      COUNT NUMBER OF OCCUPIED NEIGHBORS
                       DO MN=-LH,LH
                          MNM = MOD(N+MN+NSLICE-1,NSLICE)+1
                          DO MJ=-LH,LH
                             MJM = MOD(J+MJ+NROW-1,NROW)+1
                             DO MI=-LH,LH
                             CALL DILCHK(
     &                            X(MOD(I+MI+NSAM-1,NSAM)+1,MJM,MNM),
     &                            VALGOT,LB,NOFUSE)
                            ENDDO
                          ENDDO
                       ENDDO

                       IF (LB .GT. LAT) THEN
C                         SET CENTRAL VOXEL TO FILLED 
                          Y(I) = VALGOT
                       ENDIF
                     ENDIF
                 ENDDO
	      ENDIF

C             OUTPUT VOLUME
              CALL WRTLIN(LUN2,Y,NSAM,J+(N-1)*NROW)
           ENDDO
        ENDDO
      END	
