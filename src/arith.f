C++*********************************************************************
C
C   ARITH.F 
C              POLISH PARAMETERS                       DEC 2005 AL
C              IOFFUP = -32 BUG                        JUN 2006 AL
C              NLETO =  LEN(EXPR) bug on ifc           NOV 2007 AL
C              CALC(I...,BUF(K),BUF(K)                 MAR 2009 AL
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
C   ARITH(LUN1,LUN2,NSAM,NROW,NSLICE)
C
C   PURPOSE:  CARRIES OUT ARITHMATIC OPERATION ON IMAGE PIXEL BY PIXEL
C
C   PARAMETERS:
C        LUN1         LOGICAL UNIT NUMBER OF FILE 1
C        LUN2         LOGICAL UNIT NUMBER OF FILE 2
C        NSAM,NROW    X & Y DIMENSIONS OF FILES
C        NSLICE       Z DIMENSION OF FILES
C
C--*******************************************************************

      SUBROUTINE ARITH(LUN1,LUN2,NSAM,NROW,NSLICE,EXPR)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      PARAMETER      (IVALEN  = 40)
      PARAMETER      (IRPNLEN = 80)
      COMMON         IRPN(IRPNLEN),VAL(IVALEN)

      COMMON /IOBUF/ BUF(NBUFSIZ)

c     CHARACTER(LEN=*) ::  EXPR  ifc compiler bug reported
      CHARACTER *(*)  EXPR
      LOGICAL :: INVAR

      PARAMETER      (IOFFUP = -32)

C     SQUISH ALL BLANKS OUT OF FORMULA
      NLETO =  LEN(EXPR)

C     WRITE(6,*) 'Before SHRINKQ',NLETO,':',EXPR 
      CALL SHRINKQ(EXPR,NLETO,EXPR,NLET)
C     WRITE(6,*) 'After SHRINKQ',NLET,':',EXPR

      INVAR = .FALSE.
      DO I = 1,NLET
         IF (EXPR(I:I) .EQ. '[') THEN
            INVAR = .TRUE.
            CYCLE
         ELSEIF (EXPR(I:I) .EQ. ']') THEN
            INVAR = .FALSE.
            CYCLE
         ENDIF
         IF (.NOT. INVAR) THEN
            IF (EXPR(I:I) .GE. 'a' .AND. EXPR(I:I) .LE. 'z') THEN
C              CONVERT OPERATION TO UPPERCASE
               EXPR(I:I) = CHAR(ICHAR(EXPR(I:I)) + IOFFUP)
            ENDIF
         ENDIF
      ENDDO

C     CONVERT INPUT FORMULA TO RPN NOTATION
      CALL POLISH(0,EXPR,NLET,IRPN,NRPN,VAL,NVAL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      DO  ISL=1,NSLICE
        IOFF = (ISL-1) * NROW
        DO  I = 1,NROW
           IROW = IOFF + I
           CALL REDLIN(LUN1,BUF,NSAM,IROW)
CCCC  omp parallel do private(k) removed sept 00 due to bug, al & pp
C          may be able to return this after ftemp use  mar 09 al
           DO  K = 1,NSAM
C             CALL CALC(IRPN,NRPN,VAL,BUF(K),BUF(K),IRTFLG)
              CALL CALC(IRPN,NRPN,VAL,BUF(K),FTEMP,IRTFLG)
              BUF(K) = FTEMP
	   ENDDO
           CALL WRTLIN(LUN2,BUF,NSAM,IROW)
        ENDDO
      ENDDO
      RETURN
      END


