
C++*********************************************************************
C
C REFORM0.F                CREATED        JAN 91
C                          OPFILEC        FEB 03   -- ArDean Leith
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
C    REFORM0(LUNIN,LUNOUT,NSAM,NROW,NSLICE,MAXDIM,IRTFLG)
C
C    PURPOSE:         REFORM AN IMAGE STAGE BY ROTATION OF 90, 180 OR 270 
C                     DEGREES AROUND  THE X,Y, OR Z AXIS. DOES NOT!!
C                     USE SAME CENTERS AS STANDARD SPIDER ROTATIONS!!!
C                      
C    PARAMETERS:      MAXDIM     COMMON BUFFER SIZE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE REFORM0(LUNIN,LUNOUT,NSAM,NROW,NSLICE,MAXDIM,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON /COMMUN/           FILNAM
        CHARACTER(LEN=MAXNAM) ::  FILNAM

        CHARACTER(LEN=1) ::       NULL,AXIS
        LOGICAL ::                ERRI2

        NULL       = CHAR(0)

C       INPUT FILE OPENED IN UTIL3, FIND OUTPUT FILE NAME
 2      CALL FILERD(FILNAM,NLET,NULL,'OUTPUT',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        WRITE(NOUT,900)
900     FORMAT(' FOLLOWING THREE QUESTIONS REFER TO ORIGINAL ',
     &         'VOLUME DIMENSIONS',/)
3       NSAM1 = 1
        NSAM2 = NSAM
        CALL RDPRIS(NSAM1,NSAM2,NOT_USED,
     &  'FIRST AND LAST X COLUMN NUMBER (OR <CR> FOR ALL)',IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 2
        IF (ERRI2(NSAM1,NSAM2,2,1,NSAM,NSAM1,NSAM)) GOTO 3

4       NROW1 = 1
        NROW2 = NROW
        CALL RDPRIS(NROW1,NROW2,NOT_USED,
     &    'FIRST AND LAST Y ROW NUMBER (OR <CR> FOR ALL)',IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 3
        IF (ERRI2(NROW1,NROW2,2,1,NROW,NROW1,NROW)) GOTO 4

5       NSLICE1 = 1
        NSLICE2 = NSLICE
        CALL RDPRIS(NSLICE1,NSLICE2,NOT_USED,
     &     'FIRST AND LAST Z SLICE NUMBER (OR <CR> FOR ALL)',
     &      IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 4
        IF (ERRI2(NSLICE1,NSLICE2,2,1,NSLICE,NSLICE1,NSLICE)) GOTO 5

6       CALL RDPRMC(AXIS,NLET,.TRUE.,
     &     'ROTATION AXIS (X,Y,Z) (<CR> IS Z)',NULL,IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 5
        IF (NLET .EQ. 0) AXIS = 'Z'

7       WRITE(NOUT,*) ' ROTATIONS ARE CLOCKWISE WHEN FACING ALONG AXIS.'
        WRITE(NOUT,*) ' X AXIS POINTS TO RIGHT, Y AXIS POINTS DOWN THE '
        WRITE(NOUT,*) ' SCREEN , Z AXIS POINTS OUT OF THE SCREEN.'
        IANG = 90
        CALL RDPRI1S(IANG,NOT_USED,
     &     'CLOCKWISE ROTATION ANGLE (90,180, OR 270)',IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 6
        IF (IANG .LT. 0) IANG = -IANG
        IF (IANG .NE. 90 .AND. IANG .NE. 180 .AND. IANG .NE. 270) THEN
           CALL ERRT(31,'REFORM',NE)
           GOTO 7
        ENDIF

C       SET IMAGE SIZE AND WINDOW
        NSAM3   = NSAM2   - NSAM1   + 1
        NROW3   = NROW2   - NROW1   + 1
        NSLICE3 = NSLICE2 - NSLICE1 + 1

        IF (AXIS .EQ. 'Z') THEN
C           3-D FILE WITH Z SLICE
            IF (IANG .EQ. 0 .OR. IANG .EQ. 180) THEN
               NSAMS   = NSAM3
               NROWS   = NROW3
            ELSE
               NSAMS   = NROW3
               NROWS   = NSAM3
            ENDIF
            NSLICES = NSLICE3
            NUMVOX  = MAX(NSAM,NSAMS)

        ELSEIF (AXIS .EQ. 'Y') THEN
C           3-D FILE WITH Y SLICE 

            IF (IANG .EQ. 0 .OR. IANG .EQ. 180) THEN
               NSAMS   = NSAM3
               NSLICES = NSLICE3
            ELSE
               NSAMS   = NSLICE3
               NSLICES = NSAM3
            ENDIF

            NROWS   = NROW3
            NUMVOX  = MAX(NSAM,NSAMS)

        ELSEIF (AXIS .EQ. 'X') THEN
C           3-D FILE WITH X SLICE

            IF (IANG .EQ. 0 .OR. IANG .EQ. 180) THEN
        
               NROWS   = NROW3
               NSLICES = NSLICE3
            ELSE
               NROWS   = NSLICE3
               NSLICES = NROW3
            ENDIF
            NSAMS   = NSAM3
            NUMVOX  = NSAMS * NROWS

        ENDIF

C       MAX. NO. OF VOXELS ALLOWED 
        MAXVOX  = MAXDIM 

        IF (NUMVOX .GT. MAXVOX) THEN
           CALL ERRT(9,'REFORM0',NE)
           GOTO 99
        ENDIF
     
C       OPEN OUTPUT FILE (HAD TO WAIT TO CALCULATE SIZE FROM AXIS INFO)
        MAXIM  = 0
        IRTFLG = 0
        CALL OPFILEC(LUNIN,.FALSE.,FILNAM,LUNOUT,'U',IFORM,NSAMS,
     &                 NROWS,NSLICES,
     &                 MAXIM,' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 99

C       ROTATE THE VOLUME
        CALL REFORM(LUNIN,LUNOUT,NSAM,NSAM1,NSAM2,NSAMS,
     &              NROW,NROW1,NROW2,NROWS,
     &       NSLICE,NSLICE1,NSLICE2,NSLICES,AXIS,IANG,IRTFLG)

99      CLOSE(LUNIN)
        CLOSE(LUNOUT)

	RETURN
	END

