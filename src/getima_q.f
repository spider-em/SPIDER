C++*********************************************************************
C
C  GETIMA_Q.F      ILIST, NIMA ADDED TO ARGUMENT SEP 2000 BIMAL RATH
C                  REPLACED 'INCORE'             SEP 2001 ARDEAN LEITH
C                  LEN=MAXNAM                    JUL 2014 ARDEAN LEITH
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
C GETIMA_Q
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE GETIMA_Q(X,LSD,NSAM,NROW,IMI,ILIST,NIMA,
     &                      Q,INPIC,INCORE)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM) :: FINPAT,FINPIC,DOCFIL,OUTIMA,OUTDOC
        COMMON /FISPEC/
     &     FINPAT,FINPIC,DOCFIL,OUTIMA,OUTDOC,NLET,NLETI,NLIMA,NLDOC

        DIMENSION       X(LSD,NROW), ILIST(NIMA)
        DIMENSION       Q(*)
        LOGICAL         INCORE

        IF (INCORE)  THEN
C          CAN RECOVER IMAGE STORED IN Q ARRAY
           LBASE = NSAM * NROW * (IMI - 1)
c$omp      parallel do private(i,j)
           DO  J=1,NROW
              DO  I=1,NSAM
                 X(I,J) = Q(I + (J - 1) * NSAM + LBASE)
              ENDDO
           ENDDO

        ELSE
C          MUST READ IN IMAGE FROM FILE
           L = ILIST(IMI)
           CALL FILGET(FINPAT,FINPIC,NLET,L,INTFLAG)

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,
     &                 NSAMT,NROWT,NSLICE, MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           CALL READV(INPIC,X,LSD,NROW,NSAM,NROW,NSLICE)
           CLOSE(INPIC)
        ENDIF

        END
