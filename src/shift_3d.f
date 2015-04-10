C++*********************************************************************
C
C $$ SHIFT_3D.FOR
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
C
C--*********************************************************************
C
C $$ SHIFT_3D.FOR
C
         SUBROUTINE  SHIFT_3D(LUN1,LUN2,BD,OUT,NSAM,NROW,NSLICE,
     $   SAMSI,ROWSI,SLISI)
         DIMENSION   BD(NSAM,4),OUT(NSAM)

         SAMS=-SAMSI
         ROWS=-ROWSI
         SLIS=-SLISI
1        IF(SAMS.LT.0.0)  THEN
         SAMS=SAMS+NSAM
         GOTO 1
         ENDIF
2        IF(ROWS.LT.0.0)  THEN
         ROWS=ROWS+NROW
         GOTO 2
         ENDIF
3        IF(SLIS.LT.0.0)  THEN
         SLIS=SLIS+NSLICE
         GOTO 3
         ENDIF
C
         NSAMS=SAMS
         DX=SAMS-NSAMS
         NROWS=ROWS
         DY=ROWS-NROWS
         NSLICES=SLIS
         DZ=SLIS-NSLICES
C
         C1=(1-DX)*(1-DY)*(1-DZ)
         C2=   DX *(1-DY)*(1-DZ)
         C3=(1-DX)*   DY *(1-DZ)
         C4=(1-DX)*(1-DY)*   DZ
         C5=   DX *   DY *(1-DZ)
         C6=   DX *(1-DY)*   DZ
         C7=(1-DX)*   DY *   DZ
         C8=   DX *   DY *   DZ
C
         NSR=NROW*NSLICE
         DO    K=1,NSLICE
         KS=MOD(K+NSLICES-1,NSLICE)+1
         DO    J=1,NROW
         IF(J.EQ.1)  THEN
         JS1=(KS-1)*NROW+MOD(J+NROWS-1,NROW)+1
         JS2=(KS-1)*NROW+MOD(J+NROWS,NROW)+1
         JS3=MOD(KS*NROW,NSR)+MOD(J+NROWS-1,NROW)+1
         JS4=MOD(KS*NROW,NSR)+MOD(J+NROWS,NROW)+1
         J1=1
         J2=2
         J3=3
         J4=4
         CALL  REDLIN(LUN1,BD(1,J1),NSAM,JS1)
         CALL  REDLIN(LUN1,BD(1,J2),NSAM,JS2)
         CALL  REDLIN(LUN1,BD(1,J3),NSAM,JS3)
         CALL  REDLIN(LUN1,BD(1,J4),NSAM,JS4)
         ELSE
         JS2=(KS-1)*NROW+MOD(J+NROWS,NROW)+1
         JS4=MOD(KS*NROW,NSR)+MOD(J+NROWS,NROW)+1
         JT1=J1
         JT3=J3
         J1=J2
         J3=J4
         J2=JT1
         J4=JT3
         CALL  REDLIN(LUN1,BD(1,J2),NSAM,JS2)
         CALL  REDLIN(LUN1,BD(1,J4),NSAM,JS4)
         ENDIF
         DO    I=1,NSAM
         IS=MOD(I+NSAMS-1,NSAM)+1
         IS1=MOD(I+NSAMS,NSAM)+1
C
         OUT(I)=
     &   +C1*BD(IS,J1)+C2*BD(IS1,J1)+C3*BD(IS,J2)+C4*BD(IS,J3)
     &   +C5*BD(IS1,J2)+C6*BD(IS1,J3)+C7*BD(IS,J4)+C8*BD(IS1,J4)
         ENDDO
         CALL  WRTLIN(LUN2,OUT,NSAM,(K-1)*NROW+J)
	 ENDDO
	 ENDDO
         END
