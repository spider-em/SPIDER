C++*********************************************************************
C
C ALROSI_Q.F
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
C IMAGE_PROCESSING_ROUTINE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE  ALROSI_Q(A,ATMP,B,C,REFER,LSD,NSAM,NROW,NSI,
     &     PARA,NOCHANGE,
     &    A_CIRC,REFER_CIRC,LCIRC,JACUP,NUMR,NRING,MAXRIN,TEMP,MODE,KTN)

        PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
        PARAMETER (DGR_TO_RAD = (QUADPI/180))

        DIMENSION   A(LSD,NROW),B(LSD,NROW),C(LSD,NROW)
        DIMENSION   REFER(LSD,NROW),ATMP(LSD,NROW),PARA(3)
        INTEGER  NUMR(3,NRING),MAXRIN
        DOUBLE PRECISION  TEMP(MAXRIN+2,2),TOTMIN
        DIMENSION  A_CIRC(LCIRC),REFER_CIRC(LCIRC)
        LOGICAL  NOCHANGE
        CHARACTER*1 MODE

        NSNR=LSD*NROW

C       ATMP HAS TO BE UPDATED BEFORE THE CALL

        ITER=0
        NOCHANGE=.FALSE.
101     ITER=ITER+1
        CALL  ALRQ_Q
     &     (ATMP,LSD,NSAM,NROW,NUMR,A_CIRC,LCIRC,NRING,MODE,IPIC)
        CALL  FOURING_Q
     &     (A_CIRC,LCIRC,NUMR,NRING,TEMP,MODE)

        CALL  CROSRNG_Q
     &     (REFER_CIRC,A_CIRC,LCIRC,NRING,TEMP,TEMP(1,2),
     &     MAXRIN,JACUP,NUMR,TOTMIN,TOT,MODE)

        ROTMP=ANG(TOT,MODE)
        CALL  RTQ_Q(ATMP,B,LSD,NSAM,NROW,ROTMP)

        INS=1
        CALL  FMRS_2(B,NSAM,NROW,INS)

        LSC = NSAM+2-MOD(NSAM,2)
        CALL CCRS_2(REFER,B,C, LSC,NSAM,NROW)

        CALL  FINDMX_Q(C,LSD,NSAM,NROW,NSI,CMX1,SX1,SY1)

        DD = ((COS(ROTMP*DGR_TO_RAD)-1.0)*NUMR(1,NRING)+SX1)**2
     &     +(SIN(ROTMP*DGR_TO_RAD)*NUMR(1,NRING)+SY1)**2
        IF (DD .LT. 0.25 )  THEN
           IF (ITER .EQ. 1)  NOCHANGE=.TRUE.

C          NEXT LINE WOULD BE USED FOR NON-ZERO TOT ONLY
C          CALL  UPDTC(REFER_CIRC,A_CIRC,LCIRC,NRING,NUMR,TOT,MAXRIN,KTN)

           CALL  UPDTF(REFER_CIRC,A_CIRC,LCIRC,KTN)
           CALL  UPDTF(REFER,B,NSNR,KTN)
           RETURN
        ENDIF

        CALL  SUMAP(PARA(1),PARA(2),PARA(3),ROTMP,SX1,SY1,P1,P2,P3)
        PARA(1 )= P1
        PARA(2) = P2
        PARA(3) = P3
        CALL  RTQS_Q(A,ATMP,LSD,NSAM,NROW,PARA(1),PARA(2),PARA(3))
        IF(ITER .GT. 25)  RETURN
        GOTO 101
        END
