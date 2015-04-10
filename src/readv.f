C++*********************************************************************
C   READV.F                                 
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
C     READV(LUN,BUF,NSAM1,NROW1,NSAM,NROW,NSLICE)
C
C     PURPOSE:  READ AN IMAGE FROM A FILE USING REDLIN AND STORE
C		IN A 3D ARRAY. THE IMAGE IS STORED AT THE UPPER LEFT
C		CORNER OF THE LARGER MATRIX.
C
C     PARAMETERS:
C        LUN              I/O UNIT # OF FILE BEING READ           (SENT)
C        BUF              STORAGE AREA                       (SENT/RET.)
C	 NSAM1,NROW1      DIMENSION OF LARGER MATRIX              (SENT)
C        NSAM,NROW,NSLICE DIMENSION OF SMALLER MATRIX             (SENT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--****************************************************************

          SUBROUTINE READV(LUN,BUF,NSAM1,NROW1,NSAM,NROW,NSLICE)
         
          REAL       :: BUF(NSAM1,NROW1,1)  ! ACTUALLY NSLICE1 
	   
          DO K=1,NSLICE
             DO J=1,NROW
                CALL  REDLIN(LUN,BUF(1,J,K),NSAM,J+(K-1)*NROW)
             ENDDO
          ENDDO

          END

#ifdef USE_MPI
C         THE SAME AS READV, BUT NO MPI_BCAST IN REDLIN1P

          SUBROUTINE READV1P(LUN,BUF,NSAM1,NROW1,NSAM,NROW,NSLICE)

          REAL    BUF(NSAM1,NROW1,1)
 
          DO    K=1,NSLICE
             DO J=1,NROW
                CALL  REDLIN1P(LUN,BUF(1,J,K),NSAM,J+(K-1)*NROW)
             ENDDO
          ENDDO

          END
#endif
