C++*********************************************************************
C
C ALROSI_Q.F    REFACTOR WHILE TRACING BUG       NOV  2017 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2017  Health Research Inc.,                         *
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
C ALROSI_Q
C
C VARIABLES:
C          A      : IMAGE DATA
C          ATMP=B : SHIFTED IMAGE A
C          B=D    : MT ,  OUT: ROTATED IMAGE ATMP --> FFT 
C          C      : CROSS CORRELATED REFER WITH B
C          REFER  : FFT OF SHIFTED IMAGE
C
C CALLING:  CALL ALROSI_Q(A,B,   D,C,REFER,LSD,NX,NY,NSI,.........
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE ALROSI_Q(A,ATMP,B,C,REFER,LSD,NX,NY,NSI,
     &                      PARA,NOCHANGE,
     &                      A_CIRC,REFER_CIRC,LCIRC,JACUP,
     &                      NUMR,NRING,MAXRIN,TEMP,MODE,KTN)

        PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
        PARAMETER (DGR_TO_RAD = (QUADPI/180))

        DIMENSION         A(LSD,NY),B(LSD,NY),C(LSD,NY)
        DIMENSION         REFER(LSD,NY),ATMP(LSD,NY),PARA(3)
        INTEGER           NUMR(3,NRING),MAXRIN
        DOUBLE PRECISION  TEMP(MAXRIN+2,2),TOTMIN
        DIMENSION         A_CIRC(LCIRC),REFER_CIRC(LCIRC)
        LOGICAL           NOCHANGE
        CHARACTER*1       MODE

        NSNR = LSD * NY

C       ATMP HAS TO BE UPDATED BEFORE THE CALL

        NOCHANGE = .FALSE.

        DO ITER = 1, 25

C         INTERPOLATE ATMP INTO POLAR COORDINATES A_CIRC
          CALL ALRQ_Q(ATMP,LSD,NX,NY,NUMR,A_CIRC,LCIRC,NRING,MODE,IPIC)

C         FFT OF POLAR RINGS IN A_CIRC
          CALL FOURING_Q(A_CIRC,LCIRC,NUMR,NRING,TEMP,MODE)

C         FIND MAX LOC OF CROSS-CORRELATED FFT RINGS
          CALL CROSRNG_Q
     &     (REFER_CIRC,A_CIRC,LCIRC,NRING,TEMP,TEMP(1,2),
     &      MAXRIN,JACUP,NUMR,TOTMIN,TOT,MODE)

C         ROTATE ATMP BY ANGLE ROTMP INTO B
          ROTMP = ANG(TOT,MODE)
          CALL RTQ_Q(ATMP,B,LSD,NX,NY,ROTMP)

C         FORWARD FFT ON: B
          INS = 1
          CALL FMRS_2(B,NX,NY,INS)

C         CROSS CORRELATE  REFER WITH B (BOTH FFT'S),  OUTPUT IS: C
          CALL CCRS_2(REFER,B,C, LSD,NX,NY)

C         FIND PEAK IN: C 
          CALL FINDMX_Q(C,LSD,NX,NY,NSI,CMX1,SX1,SY1)

          !write(6,*) ' In alrosi_q,ix:', ix,iy, ic,jc, sx,sy, nsi
          !write(6,*) ' In alrosi_q,D:', d(ix,iy)

          DD = ((COS(ROTMP*DGR_TO_RAD)-1.0) * NUMR(1,NRING)+SX1)**2
     &       + (SIN(ROTMP*DGR_TO_RAD) * NUMR(1,NRING)+SY1)**2

          IF (DD < 0.25 )  THEN
             IF (ITER == 1)  NOCHANGE = .TRUE.

C            NEXT LINE WOULD BE USED FOR NON-ZERO TOT ONLY
C            CALL UPDTC(REFER_CIRC,A_CIRC,LCIRC,NRING,NUMR,TOT,MAXRIN,KTN)

C            WEIGHTED ADD A_CIRC TO REFER_CIRC
C            REFER_CIRC(I)=(REFER_CIRC(I)*(IMI-1)+A_CIRC(I))/FLOAT(IMI)
             CALL UPDTF(REFER_CIRC,A_CIRC,LCIRC,KTN)

C            WEIGHTED ADD SHIFTED IMAGE B TO THE SHIFTED REFERENCE REFER
             CALL UPDTF(REFER,B,NSNR,KTN)
             RETURN
          ENDIF

          CALL SUMAP(PARA(1),PARA(2),PARA(3),ROTMP,SX1,SY1,P1,P2,P3)

C         ROTATE IMAGE A INTO ATMP
          PARA(1) = P1
          PARA(2) = P2
          PARA(3) = P3
          CALL RTQS_Q(A,ATMP,LSD,NX,NY,PARA(1),PARA(2),PARA(3))

        ENDDO

        END

