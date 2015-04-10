C ++********************************************************************
C                                                                      *
C EROSION       CREATED                            FEB 01 ARDEAN LEITH  
C               ADDED 'ER DOC'                     MAR 01 ARDEAN LEITH 
C               ADDED NPIXER                       MAR 01 ARDEAN LEITH 
C               INCORE LUNDOC                      JUL 03 ARDEAN LEITH 
C               'ER DOC' REMOVED  (BUGGY)          FEB 14 ARDEAN LEITH
C  
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C  EROSION(LUN1,LUN2,NX,NY,NZ)
C
C  PARAMETERS:
C
C  PURPOSE:      ERODE (SHRINK) AN OBJECT IN AN IMAGE OR VOLUME 
C    
C  NOTE:         'ER DOC' is buggy in concept --> REMOVED
C                                                                 
C **********************************************************************

        SUBROUTINE EROSION(LUN1,LUN2,NX,NY,NZ)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL, ALLOCATABLE      :: VIN2(:,: )
        REAL, ALLOCATABLE      :: VIN3(:,:,:)

        CHARACTER(LEN=1)       :: MODE 
        LOGICAL                :: NEWFILE

        CALL RDPRMC(MODE,NA,.TRUE.,'BOX OR CROSS NEIGHBORHOOD (B/C)',
     &               CHAR(0),IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        LENGTH = 3
        CALL RDPRI1S(LENGTH,NOT_USED,
     &              'NEIGHBORHOOD LENGTH (ODD)',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (LENGTH <= 1) THEN
           CALL ERRT(102,'LENGTH MUST BE GREATER THAN',2) 
           RETURN

        ELSEIF (MOD(LENGTH,2) == 0) THEN
           LENGTH = LENGTH + 1
           WRITE(NOUT,90) LENGTH 
90         FORMAT(' EFFECTIVE LENGTH OF NEIGHBORHOOD: ',I0)
        ENDIF
        LH = LENGTH / 2 

        IF (NZ <= 1) THEN
           IF (MODE == 'B') THEN
              LAT = LENGTH * LENGTH - 1
           ELSE
              LAT = 2 * LENGTH - 2 
           ENDIF
        ELSE
           IF (MODE == 'B') THEN
              LAT = LENGTH * LENGTH * LENGTH - 1
           ELSE
              LAT = 3 * LENGTH - 3 
           ENDIF
        ENDIF

        WRITE(NOUT,91) LAT
91      FORMAT('  NUMBER OF NEIGHBORS:',I0)

        CALL RDPRI1S(LAT,NOT_USED,
     &        'ERODE IF NUMBER OF OCCUPIED NEIGHBORS IS < THAN',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN


        IF (NZ > 1) THEN
           ALLOCATE(VIN3(NX,NY,NZ),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
               CALL ERRT(46,'ER, VIN3',NX*NY*NZ)
               RETURN
           ENDIF

C          LOAD INPUT VOLUME
           INDEX = 0
           DO K = 1,NZ
              DO J = 1,NY
                 INDEX = INDEX + 1
                 CALL REDLIN(LUN1,VIN3(1,J,K),NX,INDEX)
              ENDDO
           ENDDO
     
           CALL EROSION3(VIN3,NX,NY,NZ,LH,LAT, MODE,LUN2,
     &                   NPIXER)

           DEALLOCATE(VIN3)
           WRITE(NOUT,'(A,I0)') '  VOXELS ERODED: ',NPIXER

        ELSE
           ALLOCATE(VIN2(NX,NY),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
               CALL ERRT(46,'ER, VIN2',NX*NY)
               RETURN
           ENDIF

C          LOAD INPUT VOLUME
           INDEX = 0
           DO J = 1,NY
              INDEX = INDEX + 1
              CALL REDLIN(LUN1,VIN2(1,J),NX,INDEX)
           ENDDO
     
           CALL EROSION2(VIN2,NX,NY,NZ,LH,LAT, MODE,LUN2,
     &                   NPIXER)
           WRITE(NOUT,'(A,I0)') '  PIXELS ERODED: ',NPIXER
           DEALLOCATE(VIN2)
        ENDIF

        END

C       ------------------------- EROSION2 -----------------------------

        SUBROUTINE EROSION2(X,NX,NY,NZ,LH,LAT,MODE,LUN2,
     &                      NPIXER)

        REAL             :: X(NX,NY)

        REAL             :: Y(NX)

        CHARACTER(LEN=1) ::   MODE 

           IKEY   = 0
           NPIXER = 0

            DO J=1,NY                      
              IF (MODE == 'C')  THEN
C                "CROSS"

                 DO I=1,NX
C                   COPY UNERODED PIXEL VALUE
                    Y(I) = X(I,J)

                    IF (X(I,J) > 0.0)  THEN
C                      CENTRAL PIXEL IS OCCUPIED (NON-ZERO)
                       LB = 0

                       DO M=-LH,LH
                          IF (M. NE. 0) THEN
C                            COUNT NUMBER OF OCCUPIED NEIGHBORS 
                             IF (X(I,MOD(J+M+NY-1,NY)+1) > 0.0) 
     &                          LB = LB + 1 
                             IF (X(MOD(I+M+NX-1,NX)+1,J) > 0.0)
     &                          LB = LB + 1
                          ENDIF
                       ENDDO
                   ENDIF
                 ENDDO
              ELSE
C                "BOX" CONNECTIVITY
                 DO I=1,NX
                    Y(I) = X(I,J)

                    IF (X(I,J) > 0.0)  THEN
C                      CENTRAL PIXEL IS OCCUPIED (NON-ZERO)

C                      COUNT NUMBER OF OCCUPIED NEIGHBORS 
                       LB = 0
                       DO MJ=-LH,LH
                          MJM = MOD(J+MJ+NY-1,NY)+1
                          DO MI=-LH,LH
                             IF (X(MOD(I+MI+NX-1,NX)+1,MJM) > 0.0)
     &                       THEN 
                               LB = LB + 1
                             ENDIF
                          ENDDO
                       ENDDO

C                      ADJUST FOR THE CENTRAL ELEMENT
                       LB = LB - 1
       
                       IF (LB < LAT) THEN
C                         ERODE CENTRAL PIXEL (SET TO ZERO) 
                          Y(I)   = 0.0
                          NPIXER = NPIXER + 1
                       ENDIF
                    ENDIF
                 ENDDO
              ENDIF

C             OUTPUT IMAGE
              CALL WRTLIN(LUN2,Y,NX,J)
           ENDDO
        END       

C       ------------------------- EROSION3 -----------------------------

        SUBROUTINE EROSION3(X,NX,NY,NZ,LH,LAT,MODE,LUN2,
     &                      NPIXER)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        REAL             :: X(NX,NY,NZ)

        REAL             :: Y(NX)

        CHARACTER(LEN=1) :: MODE 
        INTEGER          :: NX,NY,NZ
        INTEGER          :: LH,LAT,LUN2,NPIXER,IKEY,N,J,I
        INTEGER          :: LB,M,BE,IRT,MN,MNM,MJ,MJM,MI,LBE

        IKEY   = 0
        NPIXER = 0

        DO N=1,NZ                    
           DO J=1,NY 
                     
              IF (MODE == 'C')  THEN
C                "CROSS"

                 DO I=1,NX
C                   COPY UNERODED VOXEL VALUE
                    Y(I) = X(I,J,N)

                    IF (X(I,J,N) > 0.0)  THEN
C                      CENTRAL VOXEL IS OCCUPIED (NON-ZERO)

C                      COUNT NUMBER OF OCCUPIED NEIGHBORS 
                       LB = 0
                       DO M=-LH,LH
                          IF (M. NE. 0) THEN
                             IF (X(I,MOD(J+M+NY-1,NY)+1,N) > 0.0) 
     &                          LB = LB + 1
                             IF (X(MOD(I+M+NX-1,NX)+1,J,N) > 0.0)
     &                          LB = LB + 1
                             IF (X(I,J,MOD(N+M+NZ-1,NZ)+1) > 0.0)
     &                          LB = LB + 1
                          ENDIF
                       ENDDO

                       IF (LB < LAT) THEN
C                         ERODE CENTRAL VOXEL (SET TO ZERO) 
                          Y(I)   = 0.0
                          NPIXER = NPIXER + 1
                       ENDIF
                    ENDIF
                 ENDDO
              ELSE
C                "BOX" CONNECTIVITY
                 DO I=1,NX

C                   COPY UNERODED VOXEL VALUE
                    Y(I) = X(I,J,N)

                    IF (X(I,J,N) > 0.0)  THEN
C                      CENTRAL VOXEL IS OCCUPIED (NON-ZERO)

C                      COUNT NUMBER OF OCCUPIED NEIGHBORS 
                       LB = 0
                       DO MN=-LH,LH
                          MNM = MOD(N+MN+NZ-1,NZ)+1
                          DO MJ=-LH,LH
                             MJM = MOD(J+MJ+NY-1,NY)+1
                             DO MI=-LH,LH
                                IF (X(MOD(I+MI+NX-1,NX)+1,MJM,MNM)
     &                             > 0.0) LB = LB + 1
                             ENDDO
                          ENDDO
                       ENDDO

C                      ADJUST COUNT FOR THE CENTRAL ELEMENT 
                       LB = LB - 1

                       IF (LB < LAT) THEN
C                         ERODE CENTRAL VOXEL (SET TO ZERO) 
                          Y(I)   = 0.0
                          NPIXER = NPIXER + 1
                       ENDIF
                    ENDIF
                 ENDDO
              ENDIF

C             OUTPUT VOLUME LINE
              CALL WRTLIN(LUN2,Y,NX,J+(N-1)*NY)
           ENDDO
        ENDDO
      END       

