C ++********************************************************************
C                                                                      *
C                                                                      *
C                                                                      *
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
C                                                                      *
C                                                                      *
C                                                                      *
C  PURPOSE:                                                            *
C
C  SHIFT 2-D IN FOURIER SPACE ROTATING 180 DEGS WITH COPY;
C  NO-POWER-OF-TWO DIMENSIONS
C  IF SX AND SY EQUAL ZERO THEN NO SHIFT
C
C  PARAMETERS:                                                         
C
C IMAGE_PROCESSING_ROUTINE                                                                    
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE  SH180_2(X,Y,LSD2,NSAM,NROW,SX,SY)

        COMPLEX  X(LSD2,NROW),Y(LSD2,NROW)
        COMPLEX  R
        DOUBLE PRECISION  PI2

        IF (SX.EQ.0.0  .AND.  SY.EQ.0.0)  THEN
           CALL  COP(X,Y,2*LSD2*NROW)
           RETURN
        ENDIF

C       INS=1
C       CALL  FMRS_2(X,NSAM,NROW,INS)
C       IF(INS.EQ.0)  THEN
C          sx=0.0
C          sy=0.0
C          RETURN
C       ENDIF

        PI2=8.0*DATAN(1.0D0)
        PX=PI2*SX/FLOAT(NSAM)
        PY=PI2*SY/FLOAT(NROW)
c$omp parallel do private(i,j,arg,argy,iy,ix)
        DO    J=1,NROW
           IY=J-1
           IF(IY.GT.NROW/2)  IY=IY-NROW
           ARGY=PY*IY
           DO    I=1,LSD2
              IX=I-1
              ARG=PX*IX+ARGY
              Y(I,J)=CONJG(X(I,J))*CMPLX(COS(ARG),SIN(ARG))
           ENDDO
        ENDDO

        END
