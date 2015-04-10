
C++*********************************************************************
C
C WRTSEQ8.F                                    NEW JULY 97 ArDean Leith                 
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
C     WRTSEQ8:    WRITE A LINE OF 8 BIT INTEGERS INTO A SEQUENTIAL FILE
C
C     CALL WRTSEQ8(LUN,LBUF,NB,IRTFLG)
C
C        LUN      LOGICAL UNIT NUMBER OF FILE TO BE WRITTEN TO (INPUT)
C        BBUF     BYTE BUFFER CONTAINING DATA                  (INPUT)   
C        NB       NUMBER OF BYTES TO BE WRITTEN                (INPUT)
C        IRTFLG   ERROR FLAG (ZERO IS NORMAL)                  (OUTPUT)
C
C        0         2         3         4         5         6         7     
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE   WRTSEQ8(LUN,BBUF,NB,IRTFLG)


      COMMON /UNITS/LUND,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT

      INTEGER * 1    BBUF(NB)

      WRITE(LUN,ERR=999) BBUF

      IRTFLG = 0
C     WRITE(NOUT,*)' BYTES WRITTEN:',NB
      RETURN

999   WRITE(NOUT,*)' WRTSEQ8 ERROR, UNIT:',LUN,' LENGTH:',NB
      IRTFLG = 1
      RETURN
      END

