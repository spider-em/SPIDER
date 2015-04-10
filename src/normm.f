C++*******************************************************************
C
C $$ NORMM.FOR
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
C
C $$ NORMM(LUN,LUNM,NSAM,NROW,FMAX,FMIN,AV)
C
C    PARAMETERS:
C        LUN          LOGICAL UNIT NUMBER OF IMAGE
C	 LUNM	      LOGICAL UNIT NUMBER OF MASK
C        NSAM,NROW    DIMENSIONS OF IMAGE
C        NSLICE
C        FMAX         MAXIMUM OF IMAGE
C        FMIN         MINIMUM OF IMAGE
C        AV           AVERAGE OF IMAGE
C
C--*******************************************************************

      SUBROUTINE NORMM(LUN,LUNM,NSAM,NROW,NSLICE,FMAX,FMIN,AV,NPOINT)

      COMMON BUF(1)
      COMMON /MASTER/NSAMC,NROWC,IREC,NLABEL,IFORM,IMAMI,FMAXC,FMINC,
     1               AVC,SIG,IHIST
      COMMON /UNITS/LUNC,NIN,NOUT
      DOUBLE PRECISION DAV,DAV2
      LOGICAL  SETMIN

      DAV    = 0.
      DAV2   = 0.
      SETMIN=.TRUE.
      NPOINT = 0
      DO  I = 1,NROW*NSLICE
         CALL REDLIN(LUN,BUF,NSAM,I)
         CALL REDLIN(LUNM,BUF(NSAM+1),NSAM,I)
         DO 20 K = 1,NSAM
          IF(BUF(NSAM+K).GE.0.5)  THEN
           B      = BUF(K)
            IF(SETMIN)  THEN
              FMAX   = B
              FMIN   = FMAX
              SETMIN=.FALSE.
            ENDIF
           NPOINT = NPOINT+1
           FMAX   = AMAX1(B,FMAX)
           FMIN   = AMIN1(B,FMIN)
           DAV    = DAV+B
           DAV2   = DAV2+B*DBLE(B)
	  ENDIF
20       CONTINUE
      ENDDO
      FNALL = NPOINT
      AV    = DAV/FNALL
      AVC   = AV
      FMAXC = FMAX
      FMINC = FMIN
      SIG   = DSQRT((DAV2-DAV*DAV/FNALL)/DBLE(FNALL-1.0))
      END
