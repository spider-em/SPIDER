C++*********************************************************************
C
C  DEFO003.F     FIX OF UNWORKING CODE           NOV 2014 ARDEAN LEITH
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
C  DEFO003(INUM,N,KFR,A,NX,SPMAX,LUN1,BUF,MAXMINS,IRTFLG)
C
C  PURPOSE: SEARCH FOR LOCAL MINIMA ALONG: BUF
C
C  PARAMETERS:
C     INUM     POSITION OF THE IMAGE IN THE SERIES               (SENT)
C     NMIN     NUMBER OF MINIMA                                  (RET.)
C     KFR      ARRAY OF RADII     POINTS OF MINIMA               (RET.)
C     B        ARRAY OF SP. FREQ. POINTS OF MINIMA               (RET.)
C     A        ARRAY OF AMPLITUDES OF MINIMA                     (RET.)
C     NX       DIMENSION OF IMAGE                                (SENT)
C     SPMAX    MAX SPATIAL FREQUENCY                             (RET.)
C     LUN1     INPUT UNIT                                        (SENT)
C     BUF      WORKING SPACE                                     (RET.)       
C     MAXMINS  MAX NUMBER OF MINIMA                              (SENT)
C     IRTFLG   ERROR FLAG                                        (RET.)
C 
C  NOTE:  FIXED ADDRESSING OUTSIDE BUFFER BY OMITTING SEARCH FOR
C         MIN AT LEFT.  al nov 2014
C      
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE DEFO003(INUM,NMIN,KFR,B,A,NX,SPMAX,LUN1,
     &                     BUF,MAXMINS,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC' 

        INTEGER            :: INUM,NMIN,NX
        REAL               :: SPMAX 
        REAL               :: KFR(MAXMINS),B(MAXMINS),A(MAXMINS)
        REAL               :: BUF(NX)
        INTEGER            :: LUN1 
        INTEGER            :: MAXMINS,IRTFLG 

        INTEGER            :: NEIB,NOT_USED,I,J,NC
        REAL               :: SC,BX,X1,X2,X3,Y1,Y2,Y3,A1,A2,A3
        CHARACTER          :: CHO
        LOGICAL            :: FLAG
        CHARACTER          :: NULL = CHAR(0)

        FLAG   = .TRUE.
        IRTFLG = 1

        IF (INUM == 1) THEN
           CALL RDPRM1S(SPMAX,NOT_USED,
     &                  'MAX SPATIAL FREQUENCY [A-1]',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF

20      WRITE(NOUT,25)
25      FORMAT('  SEARCHING FOR MINIMA')

        CALL RDPRI1S(NEIB,NOT_USED,
     &              'SEARCH NEIGHBORHOOD DISTANCE',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       READ FIRST LINE FROM IMAGE
        CALL REDLIN(LUN1,BUF,NX,1)

        SC     = SPMAX / FLOAT(NX)
        NMIN = 0

        DO I=2,NX-1
           BX = BUF(I)

           IF (BX < BUF(I+1)) THEN
C             THIS IS A LOCAL MINIMUM

C             CHECK IF NOT A LOCAL MINIMUM OVER CURRENT NEIGHBORHOOD
              DO J=I,I-NEIB,-1
                 IF (J > 0   .AND. BX > BUF(J)) FLAG = .FALSE.
              ENDDO

              DO J=I,I+NEIB
                 IF (J <= NX .AND. BX > BUF(J))  FLAG = .FALSE.
              ENDDO

              IF (FLAG) THEN
C               A LOCAL MINIMUM OVER CURRENT NEIGHBORHOOD

                IF (NMIN >= MAXMINS) THEN
                   CALL ERRT(102,'MINIMA OVERFLOW',NMIN)
                   IRTFLG = 1
                   RETURN
                ENDIF
                NMIN = NMIN + 1

C               FIT INTO PARABOLIC
C               Y = A1+A2*X+A3*X**2
C               A3 =[(Y1-Y3)(X2-X3)-(Y2-Y3)(X1-X3)]/[(X1^2-X3^2)(X2-X3)-
C                 (X2^2-X3^2)(X1-X3)]
C               A2 = [(Y1-Y3)(X2^2-X3^2)-(Y2-Y3)(X1^2-X3^2)]/[(X1-X3)
C                (X2^2-X3^2)-(X2-X3)
C               *(X1^2-X3^2)]
C               A1 = [(Y1X3^2-Y3X1^2)(X2X3^2-X3X2^2)-(Y2X3^2-Y3X2^2)
C                  (X1X3^2-X3X1^2)]/
C               [(X3^2-X1^2)(X2X3^2-X3X2^2)-(X3^2-X2^2)(X1X3^2-X3X1^2)]
C               XMIN = -0.5A2/A3
C               YMIN = A1-0.25A2^2/A3

                X1 = FLOAT(I-1)
                X2 = FLOAT(I)
                X3 = FLOAT(I+1)
                Y1 = BUF(I-1)
                Y2 = BUF(I)
                Y3 = BUF(I+1)
                A1 = ((Y1*X3**2-Y3*X1**2) * (X2*X3**2-X3*X2**2) -
     &                (Y2*X3**2-Y3*X2**2) * (X1*X3**2-X3*X1**2))
                A1 = A1 / 
     &                 ((X3**2-X1**2)*(X2*X3**2-X3*X2**2)-(X3**2-X2**2)*
     &                 (X1*X3**2-X3* X1**2))
                A2 = ((Y1-Y3)*(X2**2-X3**2)-(Y2-Y3)*(X1**2-X3**2))
                A2 = A2 / ((X1-X3)*(X2**2-X3**2)-(X2-X3)*(X1**2-X3**2))
                A3 = ((Y1-Y3)*(X2-X3) - (Y2-Y3)*(X1-X3))
                A3 = A3 / ((X1**2-X3**2)*(X2-X3)-(X2**2-X3**2)*(X1-X3))

                A(NMIN)   = A1 - 0.25 * A2**2 / A3
                KFR(NMIN) = -0.5 * A2 / A3
             ELSE
                FLAG = .TRUE.
             ENDIF      
           ENDIF
        ENDDO

        WRITE(NOUT,*) ' ' 
        WRITE(NOUT,90) NMIN
90      FORMAT('  CURVE HAS: ',I0,'  MINIMA:')

        WRITE(NOUT,'(A)') '     #  LOCATION  LOCATION       AMPLITUDE'
        WRITE(NOUT,'(A)') '       (PIXELS)    (A-1)'

        DO I=1,NMIN
           B(I) = ABS(KFR(I)) * SC

           WRITE(NOUT,93) ,I ,KFR(I),  B(I),           A(I)
93         FORMAT('  ',  I4, F8.2, 3X,F7.4,'         ',F8.4)

c          WRITE(NOUT,*) '#',I ,KFR(I), B(I), '(A-1)   A=',A(I)
        ENDDO
        WRITE(NOUT,*) ' ' 

        CALL RDPRMC(CHO,NC,.TRUE.,'CHANGE SEARCH NEIGHBORHOOD? (Y/N)',
     &              NULL,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (CHO == 'Y') GOTO 20

        IRTFLG = 0

        END
