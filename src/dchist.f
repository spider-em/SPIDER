
C++*********************************************************************
C
C   DCHIST                        I4 UPDATED        FEB 88 ArDean Leith
C                                 AXIS ALTERED     JULY 88 ArDean Leith
C                                 LONG FILE NAMES   FEB 89 ArDean Leith
C                                 POSTSCRIPT OUTPUT JAN 99 ArDean Leith
C                                 MAXNAM            JUL 14 ArDean Leith
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
C   DCHIST(DOCF1,LUN)
C
C   PURPOSE: CREATES HISTOGRAM IN A POSTSCRIPT FILE FROM A DOC FILE
C
C--*********************************************************************

	SUBROUTINE DCHIST(LUNDOC,LUNPOS)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        PARAMETER     (NSIZE  = 2000)
	COMMON        IFREQ(660),DATA(3,NSIZE)

        CHARACTER(LEN=MAXNAM) :: FILDOC,FILPOS
        LOGICAL       ERRI2

        INCLUDE 'F90ALLOC.INC'
        REAL, DIMENSION(:,:), POINTER :: PBUF

C       READ DOCUMENT FILE, ALL KEYS AND REGISTERS 
        MAXX   = 0
        MAXY   = 0
        CALL GETDOCDAT('HISTOGRAM DOCUMENT',.TRUE.,FILDOC,LUNDOC,
     &                 .TRUE.,MAXX, MAXY,PBUF,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        KEY1 = 1
        KEY2 = MAXY
8       CALL RDPRIS(KEY1,KEY2,NOT_USED,'FIRST & LAST KEY NUMBERS',
     &              IRTFLG)
        IF (IRTFLG .NE. 0 .AND. IRTFLG .NE. -3) GOTO 9999
        IF (ERRI2(KEY1,KEY2,2,1,MAXY,KEY1,MAXY)) GOTO 8

9       CALL RDPRI1S(ICOL,NOT_USED,
     &       'REGISTER (COLUMN) IN DOC. FILE USED FOR HISTOGRAM',IRTFLG)
        IF (IRTFLG .NE. 0 .AND. IRTFLG .NE. -3) GOTO 9999
        IF (ERRI2(ICOL,IDUM,1,0,MAXX-1,IDUM,IDUM)) GOTO 9

C       FIND RANGE FOR THIS COLUMN BETWEEN DESIRED KEYS
        SXMIN = PBUF(ICOL + 1,KEY1)
        SXMAX = SXMIN

        DO IKEY=KEY1,KEY2
           ICOUNT = PBUF( 1,IKEY)
           IF (ICOUNT .GT. 0) THEN
              VALUE = PBUF( ICOL + 1,IKEY)
              SXMIN = MIN(SXMIN,VALUE)
              SXMAX = MAX(SXMAX,VALUE)
           ENDIF
        ENDDO

        WRITE(NOUT,101) SXMIN, SXMAX 
101     FORMAT(' RANGE OF X VALUES IN THIS COLUMN:',G11.3,'...',G11.3,/)

11      XINT = (SXMAX - SXMIN) / 10
        CALL RDPRM2S(SXMIN,XINT,NOT_USED,
     &               'STARTING X VALUE, INTERVAL',IRTFLG)
        IF (IRTFLG .EQ. -1 .OR. IRTFLG .GT. 0) GOTO 9
        IF (IRTFLG .GT. 0) GOTO 9999

        CALL RDPRM1S(SXMAX,NOT_USED, 'ENDING X VALUE',IRTFLG)
        IF (IRTFLG .EQ. -1 .OR. IRTFLG .GT. 0) GOTO 11

C       CALCULATE HISTOGRAM -----------------------------

C       FIND NUMBER OF BINS
        FNBIN = (SXMAX - SXMIN) / XINT
        NBIN  = FNBIN
        IF (FLOAT(NBIN) .LT. FNBIN) NBIN = NBIN + 1
        NBIN  = MIN(660,NBIN)

        DO I = 1,NBIN
           IFREQ(I) = 0
        ENDDO

        MAXBIN = 0
        ICN    = 0
        ICM    = 0
        ICC    = 0
        DO IKEY=KEY1,KEY2
           ICOUNT = PBUF( 1, IKEY)
           IF (ICOUNT .GT. 0) THEN
              VALUE  = PBUF( ICOL + 1,IKEY)
              IDIV   = INT((VALUE - SXMIN) / XINT) + 1
              IF (IDIV .LT. 0) THEN
                 ICN = ICN + 1
              ELSEIF (IDIV .GT. NBIN) THEN
                 ICM = ICM + 1
              ELSE
                 ICC   = ICC+1
                 NFREQ = IFREQ(IDIV) + 1
                 IF (NFREQ .GT. MAXBIN) MAXBIN = NFREQ
                 IFREQ(IDIV) = NFREQ
              ENDIF
           ENDIF
        ENDDO

        WRITE(NOUT,111) ICC,ICN,ICM
111     FORMAT(' VALUES USED IN HISTOGRAM: ',I6,//,
     &         ' VALUES BELOW LOWER END: ',I6,' ABOVE UPPER END: ',I6)

C       HEADING FOR BIN DATA OUTPUT
        WRITE(NOUT,113) 
113     FORMAT(//,' INTERVAL:',14X,'NUMBER IN INTERVAL')

        DO L=1,NBIN
           F1 = (L - 1) * XINT + SXMIN
           F2 = L * XINT + SXMIN
           WRITE (NOUT,911) F1,F2,IFREQ(L)
911        FORMAT(1X,F12.4,' TO ',F12.4,2X,I4)
        ENDDO

C       GET NAME OF POSTSCRIPT FILE AND OPEN AS SEQUENTIAL FORMATTED
        LENREC = 0
        CALL OPAUXFILE(.TRUE.,FILPOS,'ps',LUNPOS,LENREC,'N',
     &                       'POSTSCRIPT OUTPUT',.TRUE.,IRTFLGT)
        IF (IRTFLGT .NE. 0) GOTO 9999

        CALL POSTRT(-LUNPOS)
        CALL POSCALE(LUNPOS,1.0,1.0,  -12.0,-7.0,  125.0,102.0)

C       MAKE AXIS FOR GRAPH
        SYMIN = 0.0
        SYMAX = MAXBIN
        CALL POSAXIS('Y',SYMIN,SYMAX, 0.0,0.0, 120.0,100.0, YFACTR,
     &              LUNPOS,IRTFLG)

        CALL POSAXIS('X',SXMIN,SXMAX, 0.0,0.0, 120.0,100.0, XFACTR,
     &              LUNPOS,IRTFLG)

C       TRANSFER DATA FOR PLOT:

        X     = 0.0
        XPLUS = XINT  * XFACTR
        DO  L=1,NBIN
          KK           = (L-1) * 3 + 1
          Y            = IFREQ(L) * YFACTR
          DATA(1,KK)   = X 
          DATA(2,KK)   = Y

          X            = X + XPLUS
          DATA(1,KK+1) = X
          DATA(2,KK+1) = Y

          DATA(1,KK+2) = X
          DATA(2,KK+2) = 0.0
        ENDDO

        NDATA = 3 * NBIN
        CALL POARAYF(LUNPOS,DATA,NDATA,.FALSE.,.FALSE.)

        NLET = LNBLNKN(FILPOS)
        WRITE(NOUT,*) ' GRAPH PLACED IN: ',FILPOS(1:NLET)

9998    CLOSE(LUNPOS)

C       DEALLOCATE DOC. FILE MEMORY
9999    DEALLOCATE(PBUF)

        END
