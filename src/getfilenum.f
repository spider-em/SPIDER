
C++*********************************************************************
C
C GETFILENUM.F  -- NEW JAN 1999                   AUTHOR: ARDEAN LEITH
C                  EXTRACTED FROM LUNSETHDR AUG 02 ARDEAN LEITH
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
C *******************************************************************C **********************************************************************
C
C
C    GETFILENUM(FILNAM,IMGNUM,CALLERRT,IRTFLG)  
C
C    PURPOSE:    FINDS FILE NUMBER AT END OF FILENAME
C    
C    PARAMETERS:     FILNAM    CHAR. VARIABLE FILE NAME         (SENT)
C                    IMGNUM    NUMBER IN FILE NAME               (RET.)
C                    NDIGITS   NUMBER OF DIGITS                  (RET.)
C                    CALLERRT  CALL ERRT IF ERROR               (SENT)
C                    IRTFLG    ERROR FLAG                       (RET.)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************
 
       SUBROUTINE GETFILENUM(FILNAM,IMGNUM,NDIGITS,CALLERRT,IRTFLG)

       INCLUDE 'CMBLOCK.INC'

       CHARACTER *(*) FILNAM
       CHARACTER *1   CHARI
       LOGICAL        CALLERRT
    
C      FIND NUMBER OF CHAR. IN FILNAM
       NLET   = LNBLNKN(FILNAM)
       IGO    = NLET - 1
       IMGNUM = 0

C      EXTRACT IMGNUM FROM FILENAME
       DO I = NLET,1,-1
          CHARI = FILNAM(I:I)
          IF (CHARI .LT. '0' .OR. CHARI .GT. '9') GOTO 10
          IGO = I 
       ENDDO

10     NDIGITS = NLET - IGO + 1
       IF (NDIGITS .LE. 0 .OR. NDIGITS .GT. 10) THEN
C         NO NUMBER AT END OF FILNAM OR > 10 DIGITS
          IRTFLG = -1
          RETURN
       ENDIF

       READ(FILNAM(IGO:NLET),'(I10)',ERR=999) IMGNUM
       IRTFLG = 0
       RETURN
     
999    WRITE(NOUT,*) '*** CAN NOT GET FILE NUMBER FROM: ',FILNAM(:NLET)
       IF (CALLERRT) THEN
          CALL ERRT(100,'GETFILENUM',NE)
       ENDIF
       IRTFLG = 1

       RETURN
       END



