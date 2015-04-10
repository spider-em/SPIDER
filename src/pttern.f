C ++********************************************************************
C PTTERN
C          CHANGED RECTANGLE TO BOX             AUG 2005 ARDEAN LEITH
C          'S' BUG FOR NEW IMAGE FIXED          OCT 2008 ARDEAN LEITH
C          MENU                                 JAN 2012 ARDEAN LEITH
C          LINE INTENSITY                       JAN 2014 ARDEAN LEITH
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
C PTTERN(LUN,NX,NY,FMAXT,FMINT)
C
C PURPOSE:  TO CREATE PATTERNS FOR USE AS MASKS  0....1 IF NEW FILE
C
C PARAMETERS:  LUN           LOGICAL UNIT NUMBER
C	       NX,NY         FILE DIMENSIONS
C              FMAXT,FMINT   FILE MIN & MAX
C
C **********************************************************************

	SUBROUTINE PTTERN(LUN,NX,NY,NZ,FMAXT,FMINT)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'

        INTEGER             :: LUN,NX,NY,NZ
	REAL, INTENT(INOUT) :: FMAXT,FMINT

	REAL                :: BUF(NY)
	REAL                :: RP
        INTEGER             :: JMP,IY,IDIM,I,IRTFLG,NCHAR,NF

        CHARACTER(LEN=3)    :: CODE
        CHARACTER(LEN=1)    :: YN
        CHARACTER(LEN=1)    :: NULL = CHAR(0)

	IF (COPT == 'I' .OR. VERBOSE) WRITE(NOUT,100)
 100    FORMAT(
     &      ' .MENU: P   -- POINT'/
     &      '        L   -- LINE'/
     &      '        C   -- FILLED   CIRCLE',/,
     &      '        CL  -- OUTLINED CIRCLE',/,
     &      '        CJ  -- CIRCLE FROM 3 POINTS',/,
     &      '        T   -- FILLED   TRIANGLE',/,
     &      '        TL  -- OUTLINED TRIANGLE'/
     &      '        B   -- FILLED   BOX',/,
     &      '        BL  -- OUTLINED BOX',//)

1	CODE(1:3) = '   '
	CALL RDPRMC(CODE,NCHAR,.TRUE.,
     &              'PATTERN (P,L,C,T,B,etc)',
     &              NULL,IRTFLG)

	RP = FMAXT                     ! FILL VALUE
        IF (FMINT == FMAXT) RP = 1.0   ! NEW IMAGE FILL VALUE

c       write(6,*) 'fmin,fmax:',fmin,fmax

        IF (INDEX(CODE,'S') > 0) THEN
C          FILL DENSITY IS FMINT, BACKGROUND IS FMAXT

           IF (FMINT  ==   FMAXT) THEN
C             RECREATE NEW IMAGE WITH BACKGROUND DENSITY = 1.0 FIRST
              BUF = 1.0   ! INITIALIZE WHOLE ARRAY
              DO IY = 1,NY*NZ
	         CALL WRTLIN(LUN,BUF,NX,IY)
              ENDDO
              FMINT = 0.0
           ENDIF
           RP = FMINT      ! FILL VALUE
        ENDIF
c       write(6,*) ' fmin,fmax,rp:',fmin,fmax,rp

	JMP  = 0
	IDIM = 2

	DO I=1,3
          IF (CODE(I:I) ==  'P') JMP  = 1   ! POINT
          IF (CODE(I:I) ==  'C') JMP  = 3   ! CIRCLE
          IF (CODE(I:I) ==  'T') JMP  = 4   ! TRIANGLE
          IF (CODE(I:I) ==  'R') JMP  = 5   ! OLD BOX
          IF (CODE(I:I) ==  'B') JMP  = 6   ! BOX

          IF (CODE(I:I) ==  'L') IDIM = 1   ! LINE
          IF (CODE(I:I) ==  'J') IDIM = -1  ! UNDOCUMENTED 3 PT CIR.
        ENDDO

C       DEFAULT IS LINE!!
	IF (IDIM == 1 .AND. JMP == 0) JMP = 2  ! LINE

        SELECT CASE(JMP)

        CASE(1)
C          POINT
 	   CALL MPOINT(LUN,NX,NY,RP)

        CASE(2)
C          LINE
	   IF (NZ  ==   1)  THEN
              CALL MLINE2(LUN,NX,NY,RP,IRTFLG)
	   ELSE
	      CALL MLINE3(LUN,NX,NY,NZ,RP)
	   ENDIF

        CASE(3)
C          CIRCLE
           CALL MCIRCL(LUN,NX,NY,RP,IDIM)

        CASE(4)
C          TRIANGLE
           CALL MTRIAN(LUN,NX,NY,RP,IDIM)

        CASE(5)
C          RECTANGLE  (RELATIVE ADDRESSSING OBSOLETE AUG 2005
           CALL MRECTL(LUN,NX,NY,RP,IDIM)

        CASE(6)
C          BOX
           CALL MBOX(LUN,NX,NY,RP,IDIM)

        CASE  DEFAULT 
C          ERROR HANDLING,UNKNOWN OPTION
           CALL ERRT(101,'UNKNOWN PATTERN',NF)
        END SELECT

        CALL RDPRMC(YN,NCHAR,.TRUE.,'CONTINUE? (Y/N)',NULL,IRTFLG)
	IF (YN .NE. 'N') GOTO 1

	END


C       ------------------- MBOX -------------------------------------

	SUBROUTINE MBOX(LUN,NX,NY,RP,IDIM)

        REAL, DIMENSION(NX) :: BUF
	
	CALL RDPRMI(IX,IY,NOT_USED,
     &      'COORDINATES OF UPPER  LEFT CORNER')
	IF (IX <= 0 .OR. IY <= 0) THEN
           CALL ERRT(101,'INCONSISTENT INPUT PARAMETERS',NF)
	   RETURN
        ENDIF

	CALL RDPRMI(IXR,IYR,NOT_USED,
     &      'COORDINATES OF LOWER RIGHT CORNER')
	IF (IXR <= 0 .OR. IYR <= 0) THEN
           CALL ERRT(101,'INCONSISTENT INPUT PARAMETERS',NF)
	   RETURN
        ENDIF

	IYSTRT = MAX(1,IY)
	IYEND  = MIN(NY,IYR)
	IF (IYSTRT > IYEND .OR. IYEND <  IYSTRT) THEN
           CALL ERRT(101,'INCONSISTENT INPUT PARAMETERS',NF)
	   RETURN
        ENDIF

	IXSTRT = MAX(1,IX)
	IXEND  = MIN(NX,IXR)
	IF (IXSTRT > IXEND .OR. IXEND <  IXSTRT) THEN
           CALL ERRT(101,'INCONSISTENT INPUT PARAMETERS',NF)
	   RETURN
        ENDIF

	DO I=IYSTRT,IYEND
	   CALL REDLIN(LUN,BUF,NX,I)
	   IF ((IDIM  ==   2) .OR.
     &         (I   == IYSTRT .AND. IYSTRT == IY) .OR.
     &         (I   == IYEND  .AND. IYEND  == IYR)) THEN
              DO J=IXSTRT,IXEND
                 BUF(J) = RP
              ENDDO
           ELSE 
	      IF (IX   == IXSTRT) BUF(IX)  = RP
	      IF (IXR  == IXEND)  BUF(IXR) = RP
	   ENDIF

 	   CALL WRTLIN(LUN,BUF,NX,I)
        ENDDO

	END
