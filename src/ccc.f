
C ++********************************************************************
C                                                                      *
C CCC.F                                                                *
C                             USED REG_SET      AUG 00 ARDEAN LEITH
C                             PULLED OPFILEC    FEB 03 ARDEAN LEITH
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
C  CCC(LUN1,FILNAM1,NSAM1,NROW1,NSLICE1,LUN2,FILNAM2,
C      NSAM2,NROW2,NSLICE2,LUNM,FILNAMM)
C  
C  PURPOSE: COMPUTES THE PEARSON CORRELATION COEFFICIENT AND 
C           EUCLIDEAN DISTANCE BETWEEN TWO IMAGES/VOLUMES BY 
C           DIRECT COMPUTATION.
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE CCC(LUN1,FILNAM1,NSAM1,NROW1,NSLICE1,
     &                 LUN2,FILNAM2,NSAM2,NROW2,NSLICE2,
     &                 LUNM,FILNAMM)

        INCLUDE 'CMBLOCK.INC'

        CHARACTER(LEN=*) :: FILNAM1,FILNAM2,FILNAMM 

        REAL             :: AIMG(NSAM1), BIMG(NSAM1), CIMG(NSAM1)

        DOUBLE PRECISION :: SUM,DAV1,DVAR1,DAV2,DVAR2,DCC

        DAV1  = 0.D0
        DVAR1 = 0.D0
        DAV2  = 0.D0
        DVAR2 = 0.D0
        DCC   = 0.D0
        N     = 0
        SUM   = 0.0

        DO  I=1,NROW1*NSLICE1

           CALL REDLIN(LUN1,AIMG,NSAM1,I)
           CALL REDLIN(LUN2,BIMG,NSAM1,I)
           CALL REDLIN(LUNM,CIMG,NSAM1,I)

           DO  K=1,NSAM1
              IF (CIMG(K) .GT. 0.5) THEN
                 DAV1  = DAV1  + AIMG(K)
                 DVAR1 = DVAR1 + AIMG(K) * DBLE(AIMG(K))
                 DAV2  = DAV2  + BIMG(K)
                 DVAR2 = DVAR2 + BIMG(K) * DBLE(BIMG(K))
                 DCC   = DCC   + AIMG(K) * DBLE(BIMG(K))
                 N     = N     + 1
                 SUM   = SUM + (AIMG(K) - BIMG(K))*
     &                     DBLE(AIMG(K) - BIMG(K))
              ENDIF
           ENDDO
        ENDDO
        DAV1  = DAV1 / FLOAT(N)
        DVAR1 = DVAR1/ FLOAT(N)-DAV1**2
        DAV2  = DAV2 / FLOAT(N)
        DVAR2 = DVAR2/ FLOAT(N)-DAV2**2
        DCC   = DCC  / FLOAT(N)-DAV1*DAV2
        DCC   = DCC  / SQRT(DVAR1*DVAR2)

        DVAR1 = DSQRT(DVAR1*N / FLOAT(N-1))
        DVAR2 = DSQRT(DVAR2*N / FLOAT(N-1))

        FDCC   = DCC
        FDAV1  = DAV1
        FDVAR1 = DVAR1
        FDAV2  = DAV2
        FDVAR2 = DVAR2
        FSUM   = SUM

        CALL REG_SET_NSEL(1,5,FDCC,FDAV1,FDVAR1,
     &                      FDAV2,FDVAR2,IRTFLG)
        CALL REG_SET_NSEL(6,1,FSUM,0.0, 0.0, 0.0, 0.0,IRTFLG)

        IF (VERBOSE) THEN
C          SKIP FILE INFO IF NOT-VERBOSE

           NLET = LNBLNKN(FILNAM1)
           WRITE(NOUT,501) FILNAM1(1:NLET)
501        FORMAT('  IMAGE NO. 1: ',A)

           WRITE(NOUT,502) DAV1,DVAR1
502        FORMAT('       AV = ',G12.5,'  S.D. = ',G12.5,/)

           NLET = LNBLNKN(FILNAM2)
           WRITE(NOUT,503) FILNAM2(1:NLET)
503        FORMAT('  IMAGE NO. 2: ',A)

           WRITE(NOUT,502) DAV2,DVAR2

           NLET = LNBLNKN(FILNAMM)
           WRITE(NOUT,504) FILNAMM(1:NLET)
504        FORMAT('  MASK USED:  ',A)
        ENDIF

        WRITE(NOUT,505) DCC
505     FORMAT('       CROSS-CORRELATION COEFFICIENT:',F10.5)

        WRITE(NOUT,508) SUM
508     FORMAT('       EUCLIDEAN DISTANCE:           ',G10.3/)

        RETURN
        END

