
C++*********************************************************************
C
C  TRAF.F         LONG FILE NAMES JAN 89 ARDEAN LEITH
C                 FRAME INTRODUCED   ml 
C                 OPTC NSAM BUG FIXED & USED OPFILE  JUL 99 ARDEAN LEITH
C                 OPFILEC                            FEB 03 ARDEAN LEITH
C                 SETPRMB PARAMETERS                 MAY 09 ARDEAN LEITH
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
C   TRAF(LUN)
C
C   LUN    LOGICAL UNIT NUMBER OF FILE
C
C***********************************************************************

        SUBROUTINE TRAF(LUN)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        CHARACTER(LEN=MAXNAM) :: FILNAM
 
        REAL          B(512)
        COMMON        B

        CHARACTER     NULL,ANS,E,S,OPTC
        REAL          LAMBDA,KM

        NULL = CHAR(0)

        CALL FILERD(FILNAM,NLET,NULL,'OUTPUT',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL RDPRM(CS,NOT_USED,'CS [MM]')

        CALL RDPRM(LAMBDA,NOT_USED,'LAMBDA [A]')

        CALL RDPRM(DZ1,NOT_USED,'LOWER DEFOCUS LIMIT [A]')

        CALL RDPRM(DZ2,NOT_USED,'UPPER DEFOCUS LIMIT [A]')

        CALL RDPRMI(NSAM,NROW,NOT_USED,
     &       'NUMBER OF SPATIAL FREQ. POINTS and DEFOCUS GRID POINTS')

        CALL RDPRM(KM,NOT_USED,'MAXIMUM SPATIAL FREQUENCY [1/A]')

        CALL RDPRM(Q,NOT_USED,'SOURCE SIZE [1/A]')

        CALL RDPRM(DS,NOT_USED,'DEFOCUS SPREAD [A]')

        CALL RDPRM2(WGH,ENV,NOT_USED,
     &   'AMPLITUDE CONTRAST RATIO [0-1], GAUSSIAN ENV. HALFW [1/A]')
        IF (WGH .LT. 0.0 .OR. WGH .GT. 1.0) THEN
           CALL ERRT(31,'TRAF',NE)
           RETURN
        ENDIF

        ENV = 1.0 / ENV**2
        CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &     'DIFFRACTOGRAM, ENVELOPE or STRAIGHT (D/E/S)',NULL,IRTFLG)
        IF (ANS == 'E') IE = 1

        CALL RDPRMC(OPTC,NCHAR,.TRUE.,'FRAME WANTED? (Y/N)',NULL,IRTFLG)

        DZ  = DZ1
        DDZ = (DZ2-DZ1) / FLOAT(NROW)

C       FRAME OPTION
        IF (OPTC == 'Y') THEN
c          copied next two lines from above ml 2/2/95
           NSAMT    = NSAM + 2
           NROWT    = NROW + 2
           IOFF     = 1
           B(1)     = 1
           B(NSAMT) = 1
           IFRAME   = 1
        ELSE
           IOFF     = 2
           IFRAME   = 0
           NSAMT    = NSAM
           NROWT    = NROW
        ENDIF

C       OPEN CONVERTED TO OPEN3 JUNE 88 al
        MAXIM  = 0
        IFORM  = 1
        NSLICE = 1
        CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'U',IFORM,NSAMT,NROWT,NSLICE,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IDONE = 0
        DO  I=1,NROW

          CALL TF(B(2),CS,DZ,LAMBDA,KM,NSAM,Q,DS,IE,WGH,ENV)

          IF (OPTC == 'Y') B(NSAMT) = 1

C         ZERO DEFOCUS LINE AS PART OF FRAME
          IF (OPTC == 'Y'    .AND. 
     &        ABS(DZ) <  DDZ .AND.
     &        IDONE   == 0) THEN

            DO  K=1,NSAMT
               B(K) = 1
            END DO
            IDONE = 1
          ENDIF
C
          IF (ANS .NE. 'S') THEN
              DO  IA=2,NSAM+1
                 B(IA) = B(IA)*B(IA)
              ENDDO
           ENDIF

           CALL WRTLIN(LUN,B(IOFF),NSAMT,I+IFRAME)
           DZ = DZ + DDZ
        ENDDO

        IF (OPTC == 'Y') THEN
           DO  K=1,NSAMT
             B(K) = 1
           ENDDO

           CALL WRTLIN(LUN,B(IOFF),NSAMT,1)
           CALL WRTLIN(LUN,B(IOFF),NSAMT,NROWT)
        ENDIF

        CALL SETPRMB(LUN, 0.,0. ,0.,0.)

        END
