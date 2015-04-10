C++*********************************************************************
C   WRITEV.FOR                                 
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
C     WRITEV(LUN,BUF, NSAM1,NROW1, NSAM,NROW,NSLICE)
C
C     PURPOSE:  WRITE AN IMAGE FROM BUFFER TO A FILE USING WRTLIN. 
C		CAN WRITE FROM A ROW-LENGTH PADDED IMAGE E.G. FOURIER
C
C     PARAMETERS:
C        LUN               LOGICAL UNIT NUMBER FOR FILE BEING     (SENT)
C        BUF               BUFFER WHERE RECORD IS READ FROM       (SENT)
C	 NSAM1,NROW1       DIM. OF BUFFER CONTAINING IMAGE        (SENT)
C        NSAM,NROW,NSLICE  DIM. OF IMAGE WRITTEN TO FILE          (SENT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--****************************************************************

        SUBROUTINE WRITEV(LUN,BUF,NSAM1,NROW1,NSAM,NROW,NSLICE)
          
        REAL :: BUF(NSAM1,NROW1,1) 
	   
        DO K=1,NSLICE
           DO J=1,NROW
              CALL WRTLIN(LUN,BUF(1,J,K),NSAM,J+(K-1)*NROW)
           ENDDO
        ENDDO

        END
