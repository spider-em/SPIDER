C***********************************************************************
C
C FILDAT.F                     CREATED             DEC 87 ArDean Leith
C                              USED REG_           AUG 00 ArDean Leith
C                              IF (VERBOSE)        JUL 01 ArDean Leith
C                              MPI USE             NOV 05 ArDean Leith
C                              [] REGISTERS        NOV 05 ArDean Leith
C                              MRC SUPPORT         OCT 19 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2019  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C    FILDAT(LUN,NX)
C
C    PARAMETERS:     LUN       LOGICAL UNIT FOR FILE (ALREADY OPENED)
C                    NX        SAMPLES PER LINE
C
C    CALLED BY:      FILGEN    UTIL1
C
C    NOTE:   THE HEADER RECORD(S) OF THE FILE CONTAINS THE FOLLOWING 
C            BUFFER POSITIONS WHICH CAN BE RETRIEVED (AMONG OTHERS)
C                  7  FMAXD = IMAGE MAXIMUM
C                  8  FMIND = IMAGE MINIMUM
C		   9  AVD   = IMAGE AVERAGE
C                 10  SIG   = STANDARD DEVIATION (SQ. ROOT OF VARIANCE)
C                 11  IHIST = UNUSED
C                 14  IANGLE= FLAG INDICATING THAT TILT ANGLES STORED
C                 15  PHI   = TILT ANGLE
C                 16  THETA = TILT ANGLE
C                 17  PSI   = TILT ANGLE
C                 18  XOFF  = X OFFSET
C                 19  YOFF  = Y OFFSET
C                 20  ZOFF  = Z OFFSET
C                 21  SCALE = SCALE
C                 30  KANGLE= FLAG INDICATING THAT MORE ANGLES STORED
C                 31  PHI2  = PHI OF SECOND EULER ROTATION
C                 32  THETA2= THETA OF SECOND EULER ROTATION
C                 33  PSI2  = PSI OF SECOND EULER ROTATION
C                 34  PHI1  = PHI OF FIRST EULER ROTATION
C                 35  THETA1= THETA OF FIRST EULER ROTATION
C                 36  PSI1  = PSI OF FIRST EULER ROTATION
C
C--*********************************************************************

	   SUBROUTINE FILDAT(LUN,NX)

           IMPLICIT NONE

           INCLUDE 'CMBLOCK.INC'
           INCLUDE 'CMLIMIT.INC'
           INCLUDE 'LABLOCK.INC'

           INTEGER           :: LUN,NX

           REAL              :: BUF,VALUES
           INTEGER           :: NUMBER,ILIST
           COMMON /IOBUF/ BUF(NBUFSIZ)
           COMMON         NUMBER(20),VALUES(20),ILIST(20)

           INTEGER           :: NREG,NUMT,ILOW,IHI
           INTEGER           :: IRTFLG,LENNAME,I,NGOT

           CHARACTER(LEN=1)  :: CDUM
           CHARACTER(LEN=80) :: NAME

C          INCLUSION FOR OPTIONAL MPI INITIALIZATION. MAY ALTER MYPID  
           INTEGER           :: MYPID = -1
#include "MPI_INIT.INC"


           IF (MYPID <= 0) THEN
              WRITE(NOUT,446) FMIN,FMAX,AV,SIG
446	      FORMAT(6X,'FMIN= ',G10.3,3X,'FMAX= ',G10.3,3X,
     &                  'AV=   ',G10.3,3X,'SIG=  ',G10.3)

              IF (IANGLE > 0) WRITE(NOUT,447)PHI,THETA,PSI
447           FORMAT(6X,'PHI=  ',G10.3,3X,'THETA=',G10.3,3X,'PSI=',
     &                           G10.3) 

              WRITE(NOUT,448) XOFF,YOFF,ZOFF
448           FORMAT(6X,'XOFF= ',G10.3,3X,'YOFF= ',G10.3,3X,'ZOFF= ',
     &                           G10.3/)

              IF (SCALE .NE. 0) WRITE(NOUT,449) SCALE
449           FORMAT(6X,'SCALE=  ',G10.3)

              IF (KANGLE > 0) WRITE(NOUT,450) PHI2,THETA2,PSI2
450           FORMAT(6X,'PHI2= ',G10.3,3X,'THETA2=',G10.3,3X,'PSI2=',
     &                           G10.3)
 
              IF (KANGLE == 2) WRITE(NOUT,451) PHI1,THETA1,PSI1
451           FORMAT(6X,'PHI1= ',G10.3,3X,'THETA1=',G10.3,3X,'PSI1=',
     &                           G10.3)
        ENDIF
 

        CALL REG_GET_USED(NREG)

        IF (NREG > 0) THEN
C         REGISTERS SPECIFIED

          NUMT = MIN(NREG,20)
	  ILOW = 1
	  IHI  = 256
          CALL RDPRAI(NUMBER,20,NUMT,ILOW,IHI,
     &        'NUMBER(S) OF HEADER LOCATION TO BE RETRIEVED',
     &        CDUM,IRTFLG)

          IF (NUMBER(1) == 0) THEN
             IF (MYPID <= 0) WRITE(NOUT,900)
900          FORMAT(/,' SOME HEADER POSITIONS WHICH CAN BE RETRIEVED:'
     &           ,/,'      7  FMAX   = IMAGE MAXIMUM',/,
     &              '      8  FMIN   = IMAGE MINIMUM',/,
     &              '      9  AVD    = IMAGE AVERAGE',/,
     &              '     10  SIG    = IMAGE STANDARD DEVIATION',/,
     &              '     14  IANGLE = TILT ANGLES STORED FLAG',/,
     &              '     15  PHI    = TILT ANGLE',/,
     &              '     16  THETA  = TILT ANGLE',/,
     &              '     17  PSI    = TILT ANGLE',/,
     &              '     18  XOFF   = X OFFSET',/,
     &              '     19  YOFF   = Y OFFSET',/,
     &              '     20  ZOFF   = Z OFFSET',/,
     &              '     21  SCALE  = SCALE',/)
             RETURN
           ENDIF

C          GET HEADER VALUES FOR LOCATIONS CONTAINED IN: NUMBER ARRAY
           DO I = 1,NUMT
              CALL GETLAB(LUN,NUMBER(I),1,VALUES(I),IRTFLG)
           ENDDO

C          SET REGISTERS TO HEADER VALUES
           CALL REG_SET_NSELA(NUMT,VALUES,.FALSE.,IRTFLG)

           IF (VERBOSE) THEN
C             ECHO SETTINGS
              CALL REG_GET_SELS(ILIST,29,NGOT,IRTFLG)
              DO I = 1,NUMT
                 IF (MYPID <= 0) THEN
                    CALL REG_GET_NAME(ILIST(I),NAME,LENNAME,IRTFLG)
                    WRITE(NOUT,*) ' ',NAME(1:LENNAME),
     &                            ' SET TO:',VALUES(I)
                 ENDIF
              ENDDO
           ENDIF
        ENDIF

	END


C      ---------------------- FILDAT_MRC ------------------------- MRC

	SUBROUTINE FILDAT_MRC(LUN,NX)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'LABLOCK.INC'

        INTEGER           :: LUN,NX

        REAL              :: BUF,VALUES
        INTEGER           :: NUMBER,ILIST
        COMMON /IOBUF/ BUF(NBUFSIZ)
        COMMON         NUMBER(20),VALUES(20),ILIST(20)

        INTEGER           :: NREG,NUMT,ILOW,IHI
        INTEGER           :: IRTFLG,LENNAME,I,NGOT

        CHARACTER(LEN=1)  :: CDUM
        CHARACTER(LEN=80) :: NAME

C       INCLUSION FOR OPTIONAL MPI INITIALIZATION. MAY ALTER MYPID  
        INTEGER         :: MYPID = -1
#include "MPI_INIT.INC"
	   
        IF (MYPID <= 0) THEN
           WRITE(NOUT,446) FMIN,FMAX,AV,SIG
446	   FORMAT(6X,'FMIN= ',G10.3,3X,'FMAX= ',G10.3,3X,
     &                  'AV=   ',G10.3,3X,'SIG=  ',G10.3)

           IF (IANGLE > 0) WRITE(NOUT,447) PHI,THETA,PSI
447        FORMAT(6X,'PHI=  ',G10.3,3X,'THETA=',G10.3,3X,'PSI=',
     &                           G10.3) 

           WRITE(NOUT,448) XOFF,YOFF,ZOFF
448        FORMAT(6X,'XOFF= ',G10.3,3X,'YOFF= ',G10.3,3X,'ZOFF= ',
     &                           G10.3/)

           IF (SCALE .NE. 0) WRITE(NOUT,449) SCALE
449        FORMAT(6X,'SCALE=  ',G10.3)

           IF (KANGLE > 0)   WRITE(NOUT,451) PHI1,THETA1,PSI1
451          FORMAT(6X,'PHI1= ',G10.3,3X,'THETA1=',G10.3,3X,'PSI1=',
     &                           G10.3)
        ENDIF
 

        CALL REG_GET_USED(NREG)

        IF (NREG > 0) THEN
C         REGISTERS SPECIFIED

          NUMT = MIN(NREG,20)
	  ILOW = 1
	  IHI  = 57
          CALL RDPRAI(NUMBER,20,NUMT,ILOW,IHI,
     &        'NUMBER(S) OF HEADER LOCATION TO BE RETRIEVED',
     &        CDUM,IRTFLG)

          IF (NUMBER(1) == 0) THEN
             IF (MYPID <= 0) WRITE(NOUT,900)
900          FORMAT(/,' SOME HEADER POSITIONS WHICH CAN BE RETRIEVED:'
     &              '     20  FMIN   = IMAGE MINIMUM',/,
     &              '     21  FMAX   = IMAGE MAXIMUM',/,
     &              '     22  AVD    = IMAGE AVERAGE',/,
     &              '     55  SIG    = IMAGE STANDARD DEVIATION',/,
     &              '     43  PHI    = TILT ANGLE',/,
     &              '     45  THETA  = TILT ANGLE',/,
     &              '     33  PSI    = TILT ANGLE',/,
     &              '     11  CELLAX = X CELL DIMENSIONS IN ANG.' ,/,
     &              '     12  CELLAY = Y CELL DIMENSIONS IN ANG.' ,/,
     &              '     13  CELLAZ = Z CELL DIMENSIONS IN ANG.' ,/)

c    &              '     19  YOFF   = Y OFFSET',/,
c    &              '     20  ZOFF   = Z OFFSET',/,
c    &              '     21  SCALE  = SCALE',/)
             RETURN
           ENDIF

C          GET HEADER VALUES FOR LOCATIONS CONTAINED IN: NUMBER ARRAY
C          ALL INTERPRETED AS FLOATS !!!!!!!!!!!!!!

           DO I = 1,NUMT
              CALL GETLAB_R_MRC(LUN,NUMBER(I),1,VALUES(I),IRTFLG)
           ENDDO

C          SET REGISTERS TO HEADER VALUES
           CALL REG_SET_NSELA(NUMT,VALUES,.FALSE.,IRTFLG)

           IF (VERBOSE) THEN
C             ECHO SETTINGS
              CALL REG_GET_SELS(ILIST,29,NGOT,IRTFLG)
              DO I = 1,NUMT
                 IF (MYPID <= 0) THEN
                    CALL REG_GET_NAME(ILIST(I),NAME,LENNAME,IRTFLG)
                    WRITE(NOUT,*) ' ',NAME(1:LENNAME),
     &                            ' SET TO:',VALUES(I)
                 ENDIF
              ENDDO
           ENDIF
        ENDIF

	END


