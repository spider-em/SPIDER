
C   KNOWN BUG: IF YOU ASK FOR AXES IT FOULS UP SCALING !!!!1

C++*********************************************************************
C
C   D3DPLT.F
C
C **********************************************************************
C *  AUTHOR:  ArDean Leith                                                 *
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
C    D3DPLT(FILNM,LUN,NSAM,NROW,MAXDIM)
C
C    PURPOSE:       THIS SUBROUTINE MAKES A PERSPECTIVE THREE-
C                   DIMENSIONAL PLOT WITH HIDDEN LINES REMOVED.
C
C    PARAMETERS:    IMFILE   IMAGE FILE NAME
C                   EXTEN    EXTENSION FOR CNT FILE
C                   LUN      LOGICAL UNIT NUMBER OF IMAGE FILE
C                   NSAM     NUMBER OF SAMPLES PER ROW
C                   NROW     NUMBER OF ROWS
C                   MAXDIM   MAXIMUM BUFFER SPACE
C
C--*******************************************************************

      SUBROUTINE D3DPLT(LUN,NSAM,NROW,NSLICE,MAXDIMT)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

C     WORK BUFFER FOLLOWS DATA, IN DHIDE
      PARAMETER (NSIZE = 2000)
      COMMON DATA(3,NSIZE),BUF(1)     

      CHARACTER(LEN=MAXNAM) :: FILPOS
      CHARACTER(LEN=4)      :: ANSW
      CHARACTER(LEN=1)      :: NULL
      LOGICAL               :: PLOTNG
       
      DATA           XMINT/0./,XLNTH/-6./,YLNTH/-2.0/

      NULL   = CHAR(0)

      MAXDIM = MAXDIMT - 10000
      IF (MAXDIM .GT. 8000000) MAXDIM = 8000000

      LUNPOS  = 80

C     GET NAME OF POSTSCRIPT FILE AND OPEN AS SEQUENTIAL FORMATTED
10    CALL OPAUXFILE(.TRUE.,FILPOS,'ps',LUNPOS,0,'N',
     &               'POSTSCRIPT OUTPUT',.TRUE.,IRTFLGT)
      IF (IRTFLGT .NE. 0) RETURN
      NLETP = LNBLNKN(FILPOS)
 
      IERR = 0
      DO  I = 1,MAXDIM - 7000
        BUF(I) = .99
      ENDDO

      NSAM2  = NSAM + 2
      NGOP   = 0
      NG1OP  = -3
      N1     = -NSAM
      SCLMAX = ABS(YLNTH)
      XMAXT  = NSAM

      NSKIP = 1
1000  CALL RDPRIS(NSKIP,IDUM,NOT_USED,
     &   'SKIPPING FACTOR OR <CR> FOR NO SKIP',IRTFLG)
      IF (IRTFLG .EQ. -1) RETURN

      SCALFC = 1.0
1002  CALL RDPRM1S(SCALFC,NOT_USED,
     &   'SCALE FACTOR OR <CR> FOR NO SCALING',IRTFLG)
      IF (IRTFLG .EQ. -1) GOTO 1000

      SCLMAX = SCLMAX * SCALFC

C     SUBTRACT FROM MAXDIM, BUFFER SPACE FOR THE X AND Y ARRAYS.
C     THEN, DIVIDE REMAINING CORE INTO SIX PARTS.
C     THE ADDITIONAL '80' IS DUE TO AN EXTRA 40 BUFFER SPACES
C     FOUND AT THE END OF ARRAYS XG AND G.

200   NRES2 = (MAXDIM + 80 - 2 * NSAM2) / 6
      NRES  = NRES2 - 2
      NX    = 1     + NSAM2
      NXG   = NX    + NSAM2
      NG    = NXG   + NRES2 - 40
      NXH   = NG    + NRES2 - 40
      NH    = NXH   + NRES2
      NXG1  = NH    + NRES2
      NG1   = NXG1  + NRES2

C     BUF(1)    CORRESPONDS TO Y.
C     BUF(NX)   CORRESPONDS TO X.
C     BUF(NXG)  CORRESPONDS TO XG.
C     BUF(NG)   CORRESPONDS TO G.
C     BUF(NXH)  CORRESPONDS TO XH.
C     BUF(NH)   CORRESPONDS TO H.
C     BUF(NXG1) CORRESPONDS TO XG1.
C     BUF(NG1)  CORRESPONDS TO G1.

C     INSERT X-COORDINATES.
      DO  I = NX,NX+NSAM
        BUF(I) =I-NX
      ENDDO

1022  CALL RDPRMC(ANSW,NC,.TRUE.,'PLOT MINIMA ALSO? (N/Y)',NULL,IRTFLG)
      IF (IRTFLG .EQ. -1) GOTO 1002
      PLOTNG = (ANSW(1:1) .EQ. 'Y')

C     AXIS DISABLED AS POSTRCRIPT SCALIN IS WRONG, MAR 99
C      ANSW(:1) = 'N'
C      CALL RDPRMC(ANSW,NC,.TRUE.,
C     &   'DO YOU WANT A COORDINATE SYSTEM? (N/Y)',NULL,IRTFLG)
C      IF (IRTFLG .EQ. -1) GOTO 1022
C      IF (ANSW(1:1) .EQ. 'Y') THEN
C         WRITE(NOUT,*) ' SORRY, COORDINATE SYTEM BUGGY?? '
C         XLNTH = ABS(XLNTH)
C         YLNTH = ABS(YLNTH)
C      ENDIF

      PXLNTH = ABS(XLNTH)
      PYLNTH = ABS(YLNTH)

C     SEARCH FOR MAX AND MIN IF NOT YET AVAILABLE.
      IF (IMAMI.EQ.0) CALL NORM3(LUN,NSAM,NROW,NSLICE,FMAX,FMIN,AV)

      DELTAX = FLOAT(NSAM)   / ABS(XLNTH)
      DELTAY = (FMAX - FMIN) / SCLMAX
      IF (DELTAY .EQ. 0.0) THEN
        CALL ERRT(5,'D3DPLT',NE)
        GOTO 9999
      ENDIF

C     INITIALIZE & SET SCALING FOR POSTSCRIPT
      CALL POSTRT(-LUNPOS)

      XLL = -PXLNTH  
      YLL =  0.0
      XUR = PXLNTH + PYLNTH
      YUR = PYLNTH + (FMAX - FMIN) / DELTAY

        dtemp = (FMAX - FMIN) / DELTAY  
      ! write(6,*) '(FMAX - FMIN) / DELTAY:   ',dtemp
      !  write(6,*) 'DELTAX, DELTAY:           ',DELTAX, DELTAY
      !  write(6,*) 'XLL,YLL, XUR,YUR:         ',XLL, YLL, XUR,YUR

      CALL POSCALE(LUNPOS, 1.0,1.0,  XLL,YLL,  XUR,YUR)

C     READ IN AND PLOT EVERY NSKIPTH ROW.

      DO  I = 1,NROW,NSKIP

        CALL REDLIN(LUN,BUF,NSAM,I)
        DO  K=1,NSAM
          BUF(K) = BUF(K) - FMIN
	ENDDO

        CALL DHIDE(BUF(NX),BUF(1),BUF(NXG),BUF(NG),BUF(NXH),
     &    BUF(NH),NGOP,NRES,NSAM,NROW/NSKIP,XLNTH,YLNTH,XMINT,DELTAX,
     &    FMIN,DELTAY,XMAXT,LUNPOS)

        IF (NRES .LE. 0) THEN
           CALL ERRT(102,'D3DPLT; INSUFFICIENT BUFFER SPACE',MAXDIM)
           GO TO 9999
        ENDIF

        IF (PLOTNG) CALL DHIDE(BUF(NX),BUF(1),BUF(NXG1),BUF(NG1),
     &     BUF(NXH), BUF(NH),NG1OP,NRES,N1,0,XLNTH,YLNTH,XMINT,
     &     DELTAX,FMIN,DELTAY,XMAXT,LUNCS)

        IF (NRES .LE. 0) THEN
           CALL ERRT(102,'D3DPLT; INSUFFICIENT BUFFER SPACE2',MAXDIM)
           GO TO 9999
        ENDIF

      ENDDO

C     CLOSE THE POSTSCRIPT-FILE 
20    CALL POEND(LUNPOS)

      WRITE(NOUT,*) ' GRAPH PLACED IN: ',FILPOS(1:NLETP)

9999  CLOSE(LUNPOS)

      END

