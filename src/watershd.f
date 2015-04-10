
C ++********************************************************************
C                                                                      *
C WATERSHD                             CREATED APR 24 2002 ARDEAN LEITH                  * 
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
C  WATERSHD(LUN1,LUN2,NSAM,NROW,NSLICE)
C
C  PARAMETERS:
C
C  PURPOSE: WATERSHED AN IMAGE SKELETON IS REACHED USING A 
C           BOX (8 FOLD)  CONNECTIVITY CONVENTION
C                                                                     *
C **********************************************************************

	SUBROUTINE WATERSHD(LUN1,LUN2,NSAM,NROW,NSLICE)

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

	REAL, ALLOCATABLE, DIMENSION(:,: )  :: VIN

        LOGICAL          ::   THREED

        THREED = (NSLICE .GT. 1)
 
        IF (THREED) THEN
           WRITE(NOUT,*) 'THIS FILTER IS NOT IMPLEMENTED IN 3D!' 
           WRITE(NOUT,*) 'VOLUME WILL BE PROCESSED SLICE-BY-SLICE' 
           THREED = .FALSE.
        ENDIF

        LENGTH = 3
        LXD2   = LENGTH / 2
        LYD2   = LENGTH / 2

        ALLOCATE(VIN(NSAM,NROW),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'VIN',IER)
           RETURN
        ENDIF

        DO ISLICE = 1,NSLICE
C          LOAD INPUT IMAGE
           CALL REDVOL(LUN1,NSAM,NROW,ISLICE,ISLICE,VIN,IRTFLG)

           CALL WATERSHD2(VIN,NSAM,NROW,NSLICE,LXD2,LYD2,
     &                    ISLICE,LUN2,NPIX)
           IF (THREED) THEN
              WRITE(NOUT,*)'ISLICE: ',ISLICE,' WATERSHED PIXELS: ',NPIX 
           ELSE
              WRITE(NOUT,*)' WATERSHED PIXELS: ',NPIX 
           ENDIF
        ENDDO

        DEALLOCATE(VIN)

        END

C       ------------------------- WATERSHD2 -----------------------------


	SUBROUTINE WATERSHD2(VIN,NSAM,NROW,NSLICE,LXD2,LYD2,
     &                       ISLICE,LUN2,NPIX)

	REAL, DIMENSION(NSAM,NROW) :: VIN
        INTEGER, DIMENSION(8) ::      ILT

C       AUTOMATIC ARRAY
	REAL, DIMENSION(NSAM)      :: VOUT

        DO IY=1,NROW  
           VOUT = 0.0

           DO IX=1,NSAM

              NUMLT = 0
              VOC   = VIN(IX,IY)    
              ILOC  = 1

C             APPLY KERNAL 
              DO MY=-LYD2,LYD2
                 IYT = MOD(IY+MY+NROW-1,NROW)+1

                 DO MX=-LXD2,LXD2
C                   VALUE FOR IMAGE UNDER CURRENT KERNAL ELEMENT
                    VOK = VIN(MOD(IX+MX+NSAM-1,NSAM)+1,IYT)
              
                    ILT(ILOC) = 0
                    IF (VOC .LT. VOK) THEN
                       ILT(ILOC) = 1
                       NUMLT     = NUMLT + 1
                    ENDIF
                    ILOC      = ILOC + 1
                 ENDDO
              ENDDO

             IF (NUMLT .LE. 2) THEN
C                THIS IS A "CREEK" PIXEL 
                 VOUT(IX) = 1.0 
                 NPIX     = NPIX + 1

             ELSEIF (NUMLT .LE. 5) THEN    
C                MAYBE THIS IS A "CREEK" OR LAKE 

C                FIND TRANSITIONS 
                 ITRANS = 0
                 DO II = 1,7
                    IF (ILT(II) .NE. ILT(II+1)) ITRANS = ITRANS + 1
                 ENDDO
                 IF (ILT(1) .NE. ILT(8)) ITRANS = ITRANS + 1

                 IF ((NUMLT .EQ. 3 .AND. ITRANS .GE. 2) .OR.
     &               (NUMLT .EQ. 4 .AND. ITRANS .GE. 2) .OR.
     &               (NUMLT .EQ. 5 .AND. ITRANS .GE. 2)) 
     &              THEN
C                   THIS IS A "CREEK" PIXEL 
                    VOUT(IX) = 1.0 
                    NPIX     = NPIX + 1
                 ENDIF

              ENDIF
           ENDDO

C          OUTPUT ROW
           CALL WRTLIN(LUN2,VOUT,NSAM,(ISLICE-1)*NROW+IY)

        ENDDO
        END	

