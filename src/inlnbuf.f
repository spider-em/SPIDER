C++*********************************************************************
C
C    INLNBUF.FOR     AUG 14, 1995 ARDEAN LEITH
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
C    INLNBUF(FILNAM,NLET,NUMBUF,IRTFLG)
C
C    PARAMETERS:     FILNAM    CHAR. VARIABLE FOR FILENAME  (SENT)
C                    NLET      LENGTH OF INLINE FILE NAME   (RETURNED)
C                    NUMBUF    INLINE BUFFER NUMBER         (RETURNED)
C                    IRTFLG
C
C    PURPOSE:        RETRIEVE INLINE BUFFER NUMBER FROM FILENAME
C
C **********************************************************************

        SUBROUTINE INLNBUF(FILNAM,NLET,NUMBUF,IRTFLG)

        COMMON /UNITS/LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT

        CHARACTER *(*) FILNAM
        LOGICAL        ISDIGI

C       SET ERROR AND NON-INLINE BUFFER RETURN VALUES
        IRTFLG = 1
        NUMBUF = 0

        IF (FILNAM(1:1) .EQ. '_') THEN
C          INLINED FILE WANTED, EXTRACT BUFFER NUMBER
           LENF = LEN(FILNAM)
           IGO  = 0
           DO I =2,LENF
             IF (ISDIGI(FILNAM(I:I))) THEN
                 IF (IGO .EQ. 0) IGO = I
                 IEND = I
             ELSE
                 IF (IGO .GT. 0) GOTO 10
             ENDIF
           ENDDO

10         CONTINUE

           IF (IGO .EQ. 0) THEN 
            WRITE(NOUT,94) FILNAM
94           FORMAT('*** INVALID INLINE FILE: ',A)
             CALL ERRT(100,'INLNBUF',NE)
             RETURN
           ENDIF

           IDIG = IEND - IGO + 1
           IF (IDIG .LE. 0 .OR. IDIG .GT. 9) THEN
C            ZERO, OR MORE THAN 9 DIGITS AFTER _!!!
             WRITE(NOUT,94) FILNAM
             CALL ERRT(100,'INLNBUF',NE)
             RETURN
           ENDIF

           READ(FILNAM(IGO:IEND),'(I9)',IOSTAT=IER) NUMBUF

           IF (IER .NE. 0) THEN
              WRITE(NOUT,94) FILNAM
              CALL ERRT(100,'INLNBUF',NE)
              RETURN
           ENDIF
           IRTFLG = 0
           NLET   = IEND
        ELSE
C          NOT AN INLINE FILE
           NUMBUF = 0
           IRTFLG = 0
        ENDIF

        RETURN
        END

