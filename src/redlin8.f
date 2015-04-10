
C++*********************************************************************
C REDLIN8.F   
C                            RETURNED IOSTAT JAN 03 ArDean Leith
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
C    REDLIN8(LUN,LBUF,NB,NREC,IRTFLG)
C
C    PARAMETERS:
C        LUN    LOGICAL UNIT NUMBER OF FILE BEING READ
C        LBUF   BUFFER WHERE RECORD IS READ IN                    (RET.)
C        NB     NUMBER OF VALUES IN RECORD TO BE READ
C        NREC   RECORD TO BE READ
C        IRTFLG ERROR CODE >0 IS RETURNED IN CASE OF ERROR        (RET.)
C                           0 is NORMAL
C 
C--*******************************************************************

      SUBROUTINE REDLIN8(LUN,LBUF,NB,NREC,IRTFLG)

      LOGICAL * 1   LBUF(NB)
      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100)

C     ADD LUNARA(LUN)(FOR LABEL REC) TO NREC TO GET CORRECT RECORD NUMBER

      I = NREC + LUNARA(LUN)

      READ(LUN,REC=I,IOSTAT=IRTFLG) LBUF

      RETURN
      END

