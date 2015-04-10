
C++*********************************************************************
C
C    SEEDFILL.F                    ArDean Leith   Mar 1995
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
C   SEEDFILL()
C
C   PURPOSE:   A SEED FILL PROGRAM
C
C   PARAMETERS: MAXDIM     MAX LENGTH OF COMMON BUFFER
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE SEEDFILL(LUNIN,LUNOUT,NSAM,NROW,FMINT,MAXDIM,IRTFLG)

	INCLUDE 'CMBLOCK.INC' 
 

        CHARACTER *1   NULL

	COMMON   BUF(1)
        INTEGER  LIST(1)
        EQUIVALENCE(BUF,LIST)

        IX = 1
        IY = 1
10      CALL RDPRIS(IX,IY,NOT_USED,'SEED LOCATION',IRTFLG)
        IF (IRTFLG .EQ. -1) RETURN
        IF (IX .LT. 0 .OR. IY .LT. 0 .OR. 
     &      IX .GT. NSAM .OR. IY .GT. NROW) THEN
            CALL ERRT(101,'LOCATION NOT WITHIN IMAGE',NE)
            GOTO 10
        ENDIF
 
        CALL RDPRM2S(THRESH,FILL,NOT_USED, 
     &       'THRESHOLD & FILL VALUE',IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 10

C       LOAD INPUT IMAGE IN INPUT BUFFER
        J = 1
        DO I = 1,NROW
          CALL REDLIN(LUNIN,BUF(J),NSAM,I)
          J = J + NSAM
        ENDDO

C       COPY INPUT IMAGE TO OUTPUT BUFFER
        ISIZE = NSAM * NROW
        DO I = 1,ISIZE
          BUF(ISIZE + I) = BUF(I) 
        ENDDO

C       DETERMINE A UNIQUE VALUE FOR VISTED 
        VISITED = FMINT -1

C       PUT SEED LOCATION IN TO-BE-VISITED LIST
        NEWSTART   = NSAM * NROW
        LSTART     = 2 * NSAM * NROW + 1
        LAST       = LSTART

        ILOC       = (IY -1) * NSAM + IX
        LIST(LAST) = ILOC 

        DO WHILE (LAST .GE. LSTART) 

C          GET CURRENT LOCATION FROM TO-BE-VISITED LIST
           ILOC = LIST(LAST)

C          DECREMENT THE TO-BE-VISITED  LIST POINTER        
           LAST = LAST - 1

C          CHECK IF CURRENT LOCATION HAS BEEN VISITED
           IF (BUF(ILOC) .GT. VISITED) THEN
C             NOT VISITED YET, UPDATE OUTPUT IMAGE

              IF (BUF(ILOC) .GT. THRESH) THEN
C                DO NOT FILL, COPY THIS PIXEL TO OUTPUT
                 BUF(ILOC+NEWSTART) = BUF(ILOC)

              ELSE
C                SET THIS PIXEL TO FILL VALUE
                 BUF(ILOC+NEWSTART) = FILL

C                PUT NEIGHBORS OF CURRENT LOCATION
C                IN TO-BE-VISITED LIST, IF NOT VISITED YET
                 CALL FILLPUT(BUF,ILOC,LIST,LAST,NSAM,NROW,
     &                  VISITED,MAXDIM)
              ENDIF
           ENDIF
 
C          MARK CURRENT LOCATION AS VISITED
           BUF(ILOC) = VISITED             

        ENDDO

C       REMOVE OUTPUT IMAGE FROM OUTPUT BUFFER
        J = NSAM * NROW + 1
        DO I = 1,NROW
          CALL WRTLIN(LUNOUT,BUF(J),NSAM,I)
          J = J + NSAM
        ENDDO

	RETURN
        END



C       ***************************************************************

        SUBROUTINE FILLPUT(BUF,ILOC,LIST,LAST,NSAM,NROW,VISITED,MAXDIM)

        DIMENSION BUF(1)
        INTEGER   LIST(1)

        IROW = ((ILOC - 1) / NSAM) + 1
        ICOL = ILOC - (IROW - 1) * NSAM

        NEXT = ILOC + 1
        IF (ICOL .LT. NSAM .AND. BUF(NEXT) .GT. VISITED) THEN
C          ADD LOCATION FROM NEXT COL
           LAST = LAST + 1
           IF (LAST .GT. MAXDIM) THEN
              CALL ERRT(6,'SEEDFILL',NE)
              LAST = 0
              RETURN
           ENDIF
           LIST(LAST) = NEXT
        ENDIF

        NEXT = ILOC + NSAM
        IF (IROW .LT. NROW .AND. BUF(NEXT) .GT. VISITED) THEN
C          ADD LOCATION FROM NEXT LINE DOWN
           LAST = LAST + 1
           IF (LAST .GT. MAXDIM) THEN
              CALL ERRT(6,'SEEDFILL',NE)
              LAST = 0
              RETURN
           ENDIF
           LIST(LAST) = NEXT
        ENDIF

        NEXT = ILOC - 1
        IF (ICOL .GT. 1 .AND. BUF(NEXT) .GT. VISITED) THEN
C          ADD LOCATION FROM PREVIOUS COL
           LAST = LAST + 1
           IF (LAST .GT. MAXDIM) THEN
              CALL ERRT(6,'SEEDFILL',NE)
              LAST = 0
              RETURN
           ENDIF
           LIST(LAST) = NEXT
        ENDIF

        NEXT = ILOC - NSAM
        IF (IROW .GT. 1 .AND. BUF(NEXT) .GT. VISITED) THEN
C          ADD LOCATION FROM PREVIOUS LINE
           LAST = LAST + 1
           IF (LAST .GT. MAXDIM) THEN
              CALL ERRT(6,'SEEDFILL',NE)
              LAST = 0
              RETURN
           ENDIF
           LIST(LAST) = NEXT
        ENDIF

        RETURN
        END
