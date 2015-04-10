
C++*********************************************************************
C
C    PDATES.F                            FILNAME LENGTHENED   DEC  88 al
C                                        OUTPUT FIXED         SEPT 97 al
C                                        PAGING REMOVED       NOV  00 al
C                                        SPIREOUT             JUN  05 al  **********************************************************************
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
C    PDATES(STRING,IFMAT)
C
C    PURPOSE:  WRITES STRING WITH DATE AND TIME APPENDED
C
C    PARAMETERS:
C         STRING    CHAR. ARRAY, CONTAINING FILE NAME OR STRING (SENT)
C         IFMAT     SKIP SOME LINES BEFORE  OUTPUT              (SENT)
C
C--*********************************************************************

      SUBROUTINE PDATES(STRING,IFMAT)

      COMMON /UNITS/   LUN,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT

      CHARACTER(LEN=*)  :: STRING 
      CHARACTER(LEN=12) :: CDATT
      CHARACTER(LEN=8)  :: CTIMT
      CHARACTER(LEN=80) :: MESG

      IF (IFMAT .GE. 1) THEN
C        SKIP IFMAT LINES
         DO I = 1,IFMAT
            WRITE(NDAT,*) ' '
         ENDDO

      ELSEIF (IFMAT .NE. -1) THEN
C        SKIP ONE LINE
         WRITE(NDAT,*) ' '
      ENDIF

      CALL DATE_2K(CDATT)
      CALL MYTIME(CTIMT)

      NLET = lnblnk(STRING)
      IF (NLET .GT. 0) THEN
         WRITE(NDAT,91,ERR=223) STRING(1:NLET),CDATT(1:11),CTIMT
91       FORMAT(' ',A  ,2X,A,' at ',A,/)
         IF (IFMAT .EQ. -1) THEN
            WRITE(MESG,93) STRING(1:NLET),CDATT(1:11),CTIMT
93          FORMAT(' ',A  ,2X,A,' at ',A)
            CALL SPIREOUT(MESG,IRTFLG)
         ENDIF
      ELSE
         WRITE(NDAT,92,ERR=223) CDATT(1:11),CTIMT
92       FORMAT('  TIME: ',A,' at ',A,/)
      ENDIF

      RETURN


223   WRITE(NOUT,*) '*** PGM. ERROR IN PDATES FORMAT:'

      END

