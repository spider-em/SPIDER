C++*********************************************************************
C
C   ARITHL.F   MODIFIED FOR LOGICAL OPERATIONS         MAY   98 al
C              POLISH PARAMETERS                       DEC 2005 AL
C              I REG BUG                               JUL 2006 AL
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
C   ARITHL(LUN1,LUN2,NSAM,NROW,NSLICE)
C
C   PURPOSE:  CARRIES OUT ARITHMATIC OPERATION ON IMAGE PIXEL BY PIXEL
C             LOGICAL EXPRESSIONS SUPPORTED
C
C   PARAMETERS:
C        LUN1         LOGICAL UNIT NUMBER OF FILE 1
C        LUN2         LOGICAL UNIT NUMBER OF FILE 2
C        NSAM,NROW    X & Y DIMENSIONS OF FILES
C        NSLICE       Z DIMENSION2 OF FILES
C
C--*******************************************************************

      SUBROUTINE ARITHL(LUN1,LUN2,NSAM,NROW,NSLICE)

      INCLUDE 'CMBLOCK.INC'
 

      PARAMETER      (IVALEN = 40)
      PARAMETER      (IRPNLEN = 80)

      CHARACTER *80  EXPR1,EXPR2,EXPR3
      CHARACTER*1    NULL

      COMMON         IRPNIF1(IRPNLEN),IRPNIF2(IRPNLEN),IRPNIF3(IRPNLEN),
     &               IRPNEI1(IRPNLEN),IRPNEI2(IRPNLEN),IRPNEI3(IRPNLEN),
     &               IRPNEL1(IRPNLEN),
     &               VALIF1(IVALEN),  VALIF2(IVALEN),  VALIF3(IVALEN),
     &               VALEI1(IVALEN),  VALEI2(IVALEN),  VALEI3(IVALEN),
     &               VALEL1(IVALEN),  BUF(1) 

      COMMON /COMMUN/ ILIST(200)

      NULL = CHAR(0)
      
10    CALL RDPRMC(EXPR1,NLET1,.TRUE.,
     &    'IF (...) THEN P2=f(P1)',NULL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      IF (NLET1 .LE. 0 .OR. INDEX(EXPR1,'THEN') .LE. 0) THEN
          CALL ERRT(101,'NO -IF THEN- CLAUSE ENTERED',NE)
          RETURN
      ENDIF

C     CONVERT IF-THEN EXPRESSION FORMULAS TO RPN NOTATION 
      IF (NLET1 .GT. 0) THEN      
         CALL IFTORPN(EXPR1,ICOMPIF,
     &      IRPNIF1,NRPNIF1,VALIF1,
     &      IRPNIF2,NRPNIF2,VALIF2, 
     &      IRPNIF3,NRPNIF3,VALIF3,
     &      IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
      ENDIF

20    CALL RDPRMC(EXPR2,NLET2,.TRUE.,
     &    'ELSE IF (...) THEN P2=f(P1)',NULL,IRTFLG)
      IF (IRTFLG .EQ. -1) GOTO 10
      IF (IRTFLG .NE. 0) RETURN
      IF (NLET2 .GT. 0 .AND. INDEX(EXPR1,'THEN') .LE. 0) THEN
          CALL ERRT(101,'ILLEGAL -ELSE IF THEN- CLAUSE ENTERED',NE)
          RETURN
      ENDIF

C     CONVERT ELSE-IF-THEN EXPRESSION FORMULAS TO RPN NOTATION 
      IF (NLET2 .GT. 0) THEN      
         CALL IFTORPN(EXPR2,ICOMPEI,
     &      IRPNEI1,NRPNEI1,VALEI1,
     &      IRPNEI2,NRPNEI2,VALEI2, 
     &      IRPNEI3,NRPNEI3,VALEI3,
     &      IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
      ENDIF

30    CALL RDPRMC(EXPR3,NLET3,.TRUE., 'ELSE FORMULA',NULL,IRTFLG)
      IF (IRTFLG .EQ. -1) GOTO 20
      IF (IRTFLG .NE. 0) RETURN

C     CONVERT ELSE EXPRESSION FORMULA TO RPN NOTATION 
      IGO = INDEX(EXPR3,'=') + 1
      CALL POLISH(0,EXPR3(IGO:),NLET3-IGO+1,IRPNEL1,NRPNEL1,VALEL1,
     &            NVALEL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      DO  ISL=1,NSLICE
        IOFF = (ISL-1) * NROW

        DO  I = 1,NROW
           IROW = IOFF + I
           CALL REDLIN(LUN1,BUF,NSAM,IROW)

           DO  K = 1,NSAM
              CALL REG_SET_NSEL(1,3,FLOAT(K),FLOAT(I),FLOAT(ISL),
     &                          0.0,0.0,IRTFLG)
              IF (NLET1 .GT. 0) THEN
	         CALL CALC(IRPNIF1,NRPNIF1,VALIF1,BUF(K),A,IRTFLG)
	         CALL CALC(IRPNIF2,NRPNIF2,VALIF2,BUF(K),B,IRTFLG)

                 IF (A .LT. B) THEN
                    IF (ICOMPIF.EQ.3 .OR. ICOMPIF.EQ.4 .OR. 
     &                  ICOMPIF.EQ.6) THEN
                       CALL CALC(IRPNIF3,NRPNIF3,VALIF3,BUF(K),C,IRTFLG)                    
                       BUF(K) = C
                       CYCLE
                    ENDIF
                 ELSEIF (A .EQ. B) THEN
                    IF (ICOMPIF.EQ.1 .OR. ICOMPIF.EQ.2) THEN
                       CALL CALC(IRPNIF3,NRPNIF3,VALIF3,BUF(K),C,IRTFLG)                    
                       BUF(K) = C
                       CYCLE
                    ENDIF
                 ELSE
                    IF (ICOMPIF.EQ.5) THEN
                       CALL CALC(IRPNIF3,NRPNIF3,VALIF3,BUF(K),C,IRTFLG)                    
                       BUF(K) = C
                       CYCLE
                    ENDIF
                 ENDIF
              ENDIF

              IF (NLET2 .GT. 0) THEN
	         CALL CALC(IRPNEI1,NRPNEI1,VALEI1,BUF(K),A,IRTFLG)
	         CALL CALC(IRPNEI2,NRPNEI2,VALEI2,BUF(K),B,IRTFLG)

                 IF (A .LT. B) THEN
                    IF (ICOMPEI.EQ.3 .OR. ICOMPEI.EQ.4 .OR. 
     &                  ICOMPEI.EQ.6) THEN
                       CALL CALC(IRPNEI3,NRPNEI3,VALEI3,BUF(K),C,IRTFLG)                    
                       BUF(K) = C
                       CYCLE
                    ENDIF
                 ELSEIF (A .EQ. B) THEN
                    IF (ICOMPEI.EQ.1 .OR. ICOMPEI.EQ.2) THEN
                       CALL CALC(IRPNEI3,NRPNEI3,VALEI3,BUF(K),C,IRTFLG)                    
                       BUF(K) = C
                       CYCLE
                    ENDIF
                 ELSE
                    IF (ICOMPEI.EQ.5) THEN
                       CALL CALC(IRPNEI3,NRPNEI3,VALEI3,BUF(K),C,IRTFLG)                    
                       BUF(K) = C
                       CYCLE
                    ENDIF
                 ENDIF
              ENDIF

              IF (NLET3 .GT. 0) THEN
	         CALL CALC(IRPNEL1,NRPNEL1,VALEL1,BUF(K),C,IRTFLG)
                 BUF(K) = C
                 CYCLE
              ENDIF
	   ENDDO

           CALL WRTLIN(LUN2,BUF,NSAM,IROW)
        ENDDO
      ENDDO

      RETURN
      END



