C++*********************************************************************
C
C ARITHL.F   MODIFIED FOR LOGICAL OPERATIONS    MAY   98  ArDean Leith
C            POLISH PARAMETERS                  DEC 2005  ArDean Leith
C            I REG BUG                          JUL 2006  ArDean Leith
C            EXTRACTED PART FROM UTIL2          JAN 2020  ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2020  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C   ARITHL()
C
C   PURPOSE:  CARRIES OUT ARITHMATIC OPERATION ON IMAGE PIXEL BY PIXEL
C             LOGICAL EXPRESSIONS SUPPORTED
C
C--********************************************************************

      SUBROUTINE ARITHL()

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'
 
      REAL        :: BUF
      COMMON /IOBUF/ BUF(NBUFSIZ)         ! NBUFSIZ FROM CMLIMITS.INC
      INTEGER             :: ILIST(NIMAX) ! NIMAX   FROM CMLIMITS.INC

      INTEGER, PARAMETER  :: IVALEN  = 40
      INTEGER, PARAMETER  :: IRPNLEN = 80

      INTEGER  :: IRPNIF1(IRPNLEN),IRPNIF2(IRPNLEN),IRPNIF3(IRPNLEN),
     &            IRPNEI1(IRPNLEN),IRPNEI2(IRPNLEN),IRPNEI3(IRPNLEN),
     &            IRPNEL1(IRPNLEN)

      REAL     :: VALIF1(IVALEN),  VALIF2(IVALEN),  VALIF3(IVALEN),
     &            VALEI1(IVALEN),  VALEI2(IVALEN),  VALEI3(IVALEN),
     &            VALEL1(IVALEN) 

      CHARACTER(LEN=MAXNAM) :: EXPR1,EXPR2,EXPR3
      CHARACTER(LEN=MAXNAM) :: FILNAM1,FILNAM2
      CHARACTER(LEN=1)      :: NULL = CHAR(0)

      INTEGER               :: NX,NY,NZ,NDUM,ICOMPIF,ICOMPEI
      INTEGER               :: NRPNIF1,NRPNIF2,NRPNIF3
      INTEGER               :: NRPNEI1,NRPNEI2,NRPNEI3
      INTEGER               :: NRPNEL1,NRPNEL2,NRPNEL3
      INTEGER               :: IRTFLG,ITYPE,MAXIM1,MAXIM2
      INTEGER               :: NLET1,NLET2,NLET3,IROW,IGO
      INTEGER               :: ISL,I,K,NVALEL,IOFF
        
      INTEGER, PARAMETER    :: LUN1 = 21
      INTEGER, PARAMETER    :: LUN2 = 22

      REAL                  :: A,B,C

   
C     OPEN INPUT FILE,  BARE STACK NOT OK 
14    CALL OPFILEC(0,.TRUE.,FILNAM1,LUN1,'O',ITYPE,
     &               NX,NY,NZ,
     &               MAXIM1,'INPUT',.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     OPEN OUTPUT FILE
      CALL OPFILEC(LUN1,.TRUE.,FILNAM2,LUN2,'U',ITYPE,
     &             NX,NY,NZ,
     &             MAXIM2,'OUTPUT',.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9000

10    CALL RDPRMC(EXPR1,NLET1,.TRUE.,
     &    'IF (...) THEN P2=f(P1)',NULL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      IF (NLET1 <= 0 .OR. INDEX(EXPR1,'THEN') <= 0) THEN
          CALL ERRT(101,'NO -IF THEN- CLAUSE ENTERED',NDUM)
          RETURN
      ENDIF

C     CONVERT IF-THEN EXPRESSION FORMULAS TO RPN NOTATION 
      IF (NLET1 > 0) THEN      
         CALL IFTORPN(EXPR1,ICOMPIF,
     &      IRPNIF1,NRPNIF1,VALIF1,
     &      IRPNIF2,NRPNIF2,VALIF2, 
     &      IRPNIF3,NRPNIF3,VALIF3,
     &      IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
      ENDIF

20    CALL RDPRMC(EXPR2,NLET2,.TRUE.,
     &    'ELSE IF (...) THEN P2=f(P1)',NULL,IRTFLG)
      IF (IRTFLG  ==  -1) GOTO 10
      IF (IRTFLG .NE. 0) RETURN
      IF (NLET2 > 0 .AND. INDEX(EXPR1,'THEN') <= 0) THEN
          CALL ERRT(101,'ILLEGAL -ELSE IF THEN- CLAUSE ENTERED',NDUM)
          RETURN
      ENDIF

C     CONVERT ELSE-IF-THEN EXPRESSION FORMULAS TO RPN NOTATION 
      IF (NLET2 > 0) THEN      
         CALL IFTORPN(EXPR2,ICOMPEI,
     &      IRPNEI1,NRPNEI1,VALEI1,
     &      IRPNEI2,NRPNEI2,VALEI2, 
     &      IRPNEI3,NRPNEI3,VALEI3,
     &      IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
      ENDIF

30    CALL RDPRMC(EXPR3,NLET3,.TRUE., 'ELSE FORMULA',NULL,IRTFLG)
      IF (IRTFLG  ==  -1) GOTO 20
      IF (IRTFLG .NE. 0) RETURN

C     CONVERT ELSE EXPRESSION FORMULA TO RPN NOTATION 
      IGO = INDEX(EXPR3,'=') + 1
      CALL POLISH(0,EXPR3(IGO:),NLET3-IGO+1,IRPNEL1,NRPNEL1,VALEL1,
     &            NVALEL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      DO  ISL=1,NZ
        IOFF = (ISL-1) * NY

        DO  I = 1,NY
           IROW = IOFF + I
           CALL REDLIN(LUN1,BUF,NX,IROW)

           DO  K = 1,NX
              CALL REG_SET_NSEL(1,3,FLOAT(K),FLOAT(I),FLOAT(ISL),
     &                          0.0,0.0,IRTFLG)
              IF (NLET1 > 0) THEN
	         CALL CALC(IRPNIF1,NRPNIF1,VALIF1,BUF(K),A,IRTFLG)
	         CALL CALC(IRPNIF2,NRPNIF2,VALIF2,BUF(K),B,IRTFLG)

                 IF (A < B) THEN
                    IF (ICOMPIF == 3 .OR. ICOMPIF == 4 .OR. 
     &                  ICOMPIF == 6) THEN
                       CALL CALC(IRPNIF3,NRPNIF3,VALIF3,BUF(K),C,IRTFLG)                    
                       BUF(K) = C
                       CYCLE
                    ENDIF
                 ELSEIF (A  ==  B) THEN
                    IF (ICOMPIF == 1 .OR. ICOMPIF == 2) THEN
                       CALL CALC(IRPNIF3,NRPNIF3,VALIF3,BUF(K),C,IRTFLG)                    
                       BUF(K) = C
                       CYCLE
                    ENDIF
                 ELSE
                    IF (ICOMPIF == 5) THEN
                       CALL CALC(IRPNIF3,NRPNIF3,VALIF3,BUF(K),C,IRTFLG)                    
                       BUF(K) = C
                       CYCLE
                    ENDIF
                 ENDIF
              ENDIF

              IF (NLET2 > 0) THEN
	         CALL CALC(IRPNEI1,NRPNEI1,VALEI1,BUF(K),A,IRTFLG)
	         CALL CALC(IRPNEI2,NRPNEI2,VALEI2,BUF(K),B,IRTFLG)

                 IF (A < B) THEN
                    IF (ICOMPEI == 3 .OR. 
     &                  ICOMPEI == 4 .OR. 
     &                  ICOMPEI == 6) THEN
                       CALL CALC(IRPNEI3,NRPNEI3,VALEI3,BUF(K),C,IRTFLG)                    
                       BUF(K) = C
                       CYCLE
                    ENDIF
                 ELSEIF (A  ==  B) THEN
                    IF (ICOMPEI == 1 .OR. ICOMPEI == 2) THEN
                       CALL CALC(IRPNEI3,NRPNEI3,VALEI3,BUF(K),C,IRTFLG)                    
                       BUF(K) = C
                       CYCLE
                    ENDIF
                 ELSE
                    IF (ICOMPEI == 5) THEN
                       CALL CALC(IRPNEI3,NRPNEI3,VALEI3,BUF(K),C,IRTFLG)                    
                       BUF(K) = C
                       CYCLE
                    ENDIF
                 ENDIF
              ENDIF

              IF (NLET3 > 0) THEN
	         CALL CALC(IRPNEL1,NRPNEL1,VALEL1,BUF(K),C,IRTFLG)
                 BUF(K) = C
                 CYCLE
              ENDIF
	   ENDDO

           CALL WRTLIN(LUN2,BUF,NX,IROW)
        ENDDO
      ENDDO

9000  CLOSE(LUN1)
      CLOSE(LUN2)

      END



