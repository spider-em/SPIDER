C++*********************************************************************
C
C DEFOCUS.F   
C              MAXNAM                              JUL 14 ARDEAN LEITH
C              PROBABLY NEVER WORKED               NOV 14 ARDEAN LEITH
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
C   DEFOCUS(IRTFLG)
C
C   PURPOSE:  CALCULATE DEFOCUS AND AMPLITUDE CONTRAST COMPONENT 
C             USING LEAST SQUARES
C
C   VARIABLES:
C      RKFR   ARRAY OF MINIMUM LOCATIONS FOR ONE IMAGE
C      A      AMPLITUDE OF ALL MINIMA      
C      NP     NUMBER OF MINIMUS CHOSEN FOR EACH IMAGE
C      KP     ARRAY OF SP. FREQ. POINTS OF MINIMUS
C      ABB    ARRAY OF ABBERATIONS COORESPONDING TO EACH MINIMUS
C      NUM    NUMBER OF IMAGE IN THE SERIES
C
C--*******************************************************************

        SUBROUTINE DEFOCUS(IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL   :: BUF
        COMMON /IOBUF/ BUF(NBUFSIZ)   ! FROM CMLIMITS.INC

C       NOTE: DEFO003 PASSES BACK DATA IN BUF BUT
C             IT APPEARS UNUSED HERE, ALTHOUGH NEEDED BY NOISE.F al

        INTEGER               :: IRTFLG
        INTEGER, PARAMETER    :: MAXMINS = 20
        INTEGER, PARAMETER    :: MAXIMGS = 20
        REAL                  :: RKFR(MAXMINS),SFR(MAXMINS),AMP(MAXMINS)
        REAL                  :: KP(MAXIMGS,MAXMINS)
        REAL                  :: ABB(MAXIMGS,MAXMINS)
        INTEGER               :: NP(MAXIMGS)
        INTEGER               :: NUMLIS(MAXIMGS)
        INTEGER               :: IDUM,NUM,NOT_USED,I,MAXIM,ITYPE
        INTEGER               :: NX,NY,NZ,NMIN,NCONSTRAIN,NC,J,NE,NDUM
        REAL                  :: FMINS,SPMAX
        LOGICAL               :: UNDERFOCUS
        CHARACTER(LEN=1)      :: USEUNDER
        CHARACTER(LEN=1)      :: NULL = CHAR(0)

        CHARACTER(LEN=MAXNAM) :: IMFILE
        INTEGER, PARAMETER    :: LUN1 = 18

        IRTFLG = 0

        CALL RDPRI1S(NUM,NOT_USED,
     &              'NUMBER OF IMAGES IN THE SERIES',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        IF (NUM < 1) THEN
           CALL ERRT(102,'INVALID NUMBER OF IMAGES',NUM)
           RETURN
        ENDIF

        IF (NUM == 1) THEN
           CALL RDPRMC(USEUNDER,NC,.TRUE.,'UNDERFOCUS? (Y/N)',
     &                 NULL,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
           UNDERFOCUS =  (USEUNDER .NE. 'N') 
        ELSE
C          USING IMAGE SERIES
           WRITE(NOUT,*)' INPUT IMAGE IN SEQUENCE'
        ENDIF

        DO  I=1,NUM   ! LOOP OVER ALL IMAGES IN SERIES
           IF (NUM > 1) WRITE(NOUT,'(A,I0,A)')'# ',I,' IMAGE'

           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,IMFILE,LUN1,'O',ITYPE,NX,NY,NZ,
     &             MAXIM,'IMAGE', .FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           WRITE(NOUT,90) NX,NY
90         FORMAT('  IMAGE DIMENSIONS: ', I0,' x ',I0)

C          FIND UP TO 20 MINIMA ALONG CURVE
           CALL DEFO003(I,NMIN,RKFR,SFR,AMP,NX,SPMAX,LUN1,
     &                  BUF,MAXMINS,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
           CLOSE(LUN1)

           CALL RDPRI1S(NP(I),NOT_USED,
     &                 'NUMBER OF MINIMA USED FOR CTF',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
           IF (NP(I) < 1 .OR. NP(I) > NMIN) THEN
              CALL ERRT(102,'MINIMA OUTSIDE VALID RANGE',NMIN)
              GOTO 9999
           ENDIF

           DO  J=1,NP(I)           ! LOOP OVER ALL DESIRED MINIMA

              KP(I,J)  = RKFR(J)    ! SPATIAL FREQ OF MINIMA
              ABB(I,J) = J - 1     ! ABBERATIONS IN PI UNITS
  
              IF (UNDERFOCUS) ABB(I,J) = -J   
           ENDDO
        ENDDO

        !write(6,*) np(1),  kp(1,1),kp(1,2),kp(1,3),kp(1,4)
        !write(6,*) abb(1,1),abb(1,2),abb(1,3),abb(1,4)

        CALL FLUSHRESULTS
        IF (NUM > 1) THEN
           WRITE(NOUT,*) ' ENTER CONSTRAINTS: '
           WRITE(NOUT,*) '           (1) SAME AMPLITUDE,'
           WRITE(NOUT,*)
     &     '           (2) SAME AMPLITUDE AND DEFINE DEFOCUS INTERVAL'

           CALL RDPRI1S(NCONSTRAIN, NOT_USED, 
     &                  'CONSTRAINTS: (1 or 2)',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           IF (NCONSTRAIN < 1 .OR. NCONSTRAIN > 2) THEN
              CALL ERRT(102,'UNDEFINED CONSTRAINT',NCONSTRAIN)
              RETURN

           ELSEIF (NCONSTRAIN == 2) THEN
              CALL DEFO001(NUM,NP,KP,ABB,NX,SPMAX)

           ELSE
              IF (NUM > 11) THEN
                 CALL DEFO2001(NUM,NP,KP,ABB,NX,SPMAX)
              ELSE
                 CALL DEFO1001(NUM,NP,KP,ABB,NX,SPMAX)
              ENDIF
           ENDIF
        ELSE
C          SINGLE INPUT FILE CASE
           CALL DEFO1001(NUM,NP,KP,ABB,NX,SPMAX)
        ENDIF

9999    CLOSE(LUN1)

        END
