
C++*********************************************************************
C
C   WRTLIN8.F                            
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
C     WRTLIN8(LUN,LBUF,NB,NREC)
C
C     PURPOSE:    WRITE A LINE OF 8 BIT INTEGERS INTO A FILE
C
C     PARAMETERS:
C        LUN    LOGICAL UNIT NUMBER OF FILE TO BE WRITTEN TO
C        LBUF   LOGICAL * 1 BUFFER THAT RECORD IS READ FROM
C        NB     NUMBER OF VALUES IN RECORD TO BE WRITTEN
C        NREC   RECORD TO BE WRITTEN INTO
C
C        0         2         3         4         5         6         7     
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE   WRTLIN8(LUN,LBUF,NB,NREC)

 
      LOGICAL * 1     LBUF(NB)

      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100)

C     ADD LUNARA(LUN)(FOR LABEL REC) TO NREC TO GET THE CORRECT RECORD NUMBER
      I = NREC + LUNARA(LUN)
      WRITE(LUN,ERR=100,REC=I) LBUF

      RETURN

100   WRITE(6,*)' WRTLIN8 ERROR, UNIT:',LUN,' REC:',NREC,' LENGTH:',NB
      RETURN
      END

