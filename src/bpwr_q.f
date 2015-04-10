C++*********************************************************************
C
C BPWR_Q.F
C                  OPFILEC                         FEB 03 ARDEAN LEITH
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
C   BPWR_Q(BUF,B,W,FM,ILIST,NANG,N2S,N2R,NSAM,NROW)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE BPWR_Q(BUF,B,W,FM,ILIST,NANG,N2S,N2R,
     &                     NSAM,NROW,FINPAT,NLET,FOUT,NLOT)

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC'

         CHARACTER (LEN=*)      :: FINPAT,FOUT

         CHARACTER (LEN=MAXNAM) :: FINPIC

         REAL                   :: B(N2S+2,N2R),W(N2S+2,N2R)
         DOUBLE PRECISION       :: AVE
         INTEGER                :: ILIST(NANG)
         REAL                   :: BUF(*)
 
         INTEGER, PARAMETER     :: LUNI = 21

         NS2 = N2S/2
         NR2 = N2R/2
         X1  = FLOAT(NS2)**2
         Y1  = FLOAT(NR2)**2

C        PREPARE WEIGHTING FUNCTION R**2
         DO J=1,N2R
            IY = J-1
            IF (IY .GT. NR2)  IY=IY-N2R
            DO    I=1,N2S+2,2
               IX = (I-1)/2
               IF (FM .GE. 0.0)  THEN
                  FQ = 0.25*(FLOAT(IX*IX)/X1+FLOAT(IY*IY)/Y1)
               ELSE
                  FQ = 0.5*SQRT(FLOAT(IX*IX)/X1+FLOAT(IY*IY)/Y1)
               ENDIF
               W(I,J)   = FQ
               W(I+1,J) = FQ
            ENDDO
	 ENDDO
         W(1,1) = 1.0
         W(2,1) = 1.0

C       PARZEN FILTER

	FM = ABS(FM)
	IF (.NOT.(FM.EQ.0.0.OR.FM.EQ.1.0))  THEN
           DO    J=1,N2R
              IY=J-1
              IF (IY .GT. NR2)  IY=IY-N2R
              DO    I=1,N2S+2,2
                 IX   = (I-1)/2
                 FQ   = 0.5*SQRT(FLOAT(IX*IX)/X1+FLOAT(IY*IY)/Y1)
                 PARZ = 0.0
                 IF (FQ .LE. FM)  THEN
                    IF (FQ.LE.FM/2.0)  THEN
                       PARZ = 1.0-6.0*(FQ/FM)**2*(1.0-FQ/FM)
                    ELSE
                       PARZ = 2.0*(1.0-FQ/FM)**3
                    ENDIF
                 ENDIF
                 W(I,J)   = W(I,J)*PARZ
                 W(I+1,J) = W(I+1,J)*PARZ
              ENDDO
           ENDDO
	ENDIF

         DO  K=1,NANG
C           READ ONE PROJECTION
            CALL  FILGET(FINPAT,FINPIC,NLET,ILIST(K),INTFLAG)
            MAXIM = 0
            CALL OPFILEC(0,.FALSE.,FINPIC,LUNI,'O',IFORM,LSAM,LROW,NSL,
     &               MAXIM,' ',.FALSE.,IRTFLG)
            IF (IRTFLG .NE. 0) THEN
               WRITE(NOUT,2032)  FINPIC
2032           FORMAT(' FILE SKIPPED: ',A)
	    ELSE

C              READ  IMAGE
               DO    I=1,NROW
                  CALL  REDLIN(LUNI,B(1,I),NSAM,I)
               ENDDO
               CLOSE(LUNI)

C              PADDING WITH BORDER AVERAGE
               IF (N2S.NE.NSAM .AND. N2R.NE.NROW)  THEN
                  AVE=0.0
                  DO    J=1,NROW
                     AVE=AVE+B(NSAM,J)
                     AVE=AVE+B(1,J)
                  ENDDO

                  DO    I=1,NSAM
                     AVE = AVE+B(I,1)
                     AVE = AVE+B(I,NROW)
                  ENDDO

                  AVE = AVE/(2*FLOAT(NSAM)+2*FLOAT(NROW))

                  DO    J=1,N2R
                     DO I=NSAM+1,N2S
                        B(I,J) = AVE
                     ENDDO
                  ENDDO

                  DO J=NROW+1,N2R
                     DO I=1,NSAM
                        B(I,J) = AVE
                     ENDDO
                  ENDDO
               ENDIF

               INV = +1
               CALL FMRS_2(B,N2S,N2R,INV)

C              APPLY FILTER

c$omp          parallel do private(i,j)
               DO J=1,N2R
                  DO I=1,N2S+2
                     B(I,J) = B(I,J)*W(I,J)
                  ENDDO
               ENDDO

               INV=-1
               CALL  FMRS_2(B,N2S,N2R,INV)

C              WRITE IMAGE
               CALL FILGET(FOUT,FINPIC,NLOT,ILIST(K),INTFLAG)

               MAXIM = 0
               CALL OPFILEC(0,.FALSE.,FINPIC,LUNI,'U',IFORM,
     &               LSAM,LROW,NSL,
     &               MAXIM,' ',.FALSE.,IRTFLG)

               DO I=1,NROW
                  CALL WRTLIN(LUNI,B(1,I),NSAM,I)
               ENDDO
               CLOSE(LUNI)
            ENDIF

C        ENDIF COMES FROM SKIPPING A FILE
	 ENDDO

         END
