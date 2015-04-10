C++*********************************************************************
C
C  SAQB   NON-POWER-OF-TWO DIMENSIONS             JUL 91
C         OPFILEC                                 FEB 03  ARDEAN LEITH
C         SALB_P INSERTED                         MAR 12  ARDEAN LEITH
C         RDPRI1S BUG,COMMON REMOVED,COSMETIC     SEP 12  ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C  SAQB
C
C  PURPOSE:
C     SHIFT ALIGNMENT, NON-POWER-OF-TWO DIMENSIONS
C     SUBTRACTION OF AN IMAGE FROM THE AVERAGE.
C     QUADRATIC INTERPOLATION AS AN OPTION.
C     SCRATCH FILE ON THE DISK
C
C  PROCEDURES/FUNCTIONS CALLED:
C       SAQB
C       SAQB_P
C       SAQB_F
C       UPDTF(C,A,N,IMI)
C       COP(A,B,N)
C       CRSM_2(X,Y,O,NSAM,NROW,WRK)
C       MLC(X,Y,O,N)
C       SHFC_2(X,Y,NSAM,NROW,WRK,SX,SY)
C       SH180_2(X,Y,NSAM,NROW,WRK,SX,SY)
C       SHFM_2(X,NSAM,NROW,WRK,SX,SY)
C       CR180_2(X,Y,O,NSAM,NROW,WRK)
C       MJC(X,Y,O,N)
C       FMR_2(X,NSAM,NROW,WORK,INV)
C       FMR_1(X,N,WORK,INV)
C       FFTMCF (A,B,NTOT,N,NSPAN,ISN)
C       FINDMX(D,NSAM,NROW,CMX,SX,SY,JACUP)
C       ENFR_2(A,NSAM,NROW)
C       RTQ(X,OUT,NSAM,NROW,THETA)
C       QUADRI(XX, YY, NXDATA, NYDATA, FDATA)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE SAQB

        INCLUDE 'CMLIMIT.INC'   
        INCLUDE 'CMBLOCK.INC' 
 
        COMMON /IOBUF/ BUF(NBUFSIZ)

        CHARACTER (LEN=MAXNAM) :: FINPAT,DOCFIL,FINPIC

        REAL, ALLOCATABLE      :: CNEW(:,:)
        CHARACTER(LEN=1)       :: FLIP 
        CHARACTER(LEN=1)       :: NULL = CHAR(0) 

        INTEGER, PARAMETER     :: INPIC = 77
        INTEGER, PARAMETER     :: NDOC  = 55 

        CALL FILERD(FINPAT,NLET,NULL,
     &           'INPUT FILE TEMPLATE (E.G. PIC****)~',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

         CALL FILERD(DOCFIL,NLETI,NULL,
     &        'DOCUMENT FILE CONTAINING GROUP ASSIGNMENT~',IRTFLG)
        IF (IRTFLG .NE. 0)  RETURN

        CALL RDPRI1S(NGRP,NOT_USED,'GROUP NUMBER',IRTFLG)
        IF (IRTFLG .NE. 0)  RETURN

        CALL RDPRMC(FLIP,NA,.TRUE.,'CHECK 180 DEGREE ROTATION (Y/N)',
     &      NULL,IRTFLG)
        IF (IRTFLG .NE. 0)  RETURN

C       CALL RDPRI1S(JACUP,NOT_USED,PRECISION OF PEAK LOCATION (0..100)')
        JACUP = 0

        K     = 0
        K2    = 1
        NIMA  = 0
778     LERR  = -1
        KP1   = K+1
        CALL  UNSAV(DOCFIL,K,NDOC,KP1,BUF,4,LERR,K2)
        IF (LERR == 0)  THEN
           IF (IFIX(BUF(4)) == NGRP)  NIMA = NIMA + 1
           K = K + 1

C          PICK UP ONE OF THE IMAGES
           IMAGE = IFIX(BUF(1))
           GOTO  778
        ENDIF

        IF (NIMA == 0)  THEN
           WRITE(NOUT,*) ' *** DESIRED GROUP NOT FOUND !'
           CLOSE(NDOC)
           RETURN
        ENDIF

        write(nout,*) ' after unsav -------------------------------'
C       OPEN FIRST IMAGE FILE TO DETERMINE NSAM, NROW, NSL

        CALL FILGET(FINPAT,FINPIC,NLET,IMAGE,INTFLG)

        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NSAM,NROW,NSL,
     &             MAXIM,'DUMMY',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        CLOSE(INPIC)

        LSD = NSAM + 2 - MOD(NSAM,2)

        ALLOCATE (CNEW(LSD,NROW), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'SAQB; CNEW',LSD*NROW)
           GOTO 9999
        ENDIF

        IF (FLIP .EQ. 'Y')  THEN
           CALL SAQB_P(BUF,LSD,NSAM,NROW,NIMA,NGRP,JACUP,
     &                 CNEW,FINPAT,NLET,DOCFIL,FINPIC)           
        ELSE
           CALL SAQB_F(BUF,LSD,NSAM,NROW,NIMA,NGRP,JACUP,
     &                 CNEW,FINPAT,NLET,DOCFIL,FINPIC)     
        ENDIF

        MAXIM  = 0
        NSLICE = 1
        IFORM  = 1
        CALL OPFILEC(0,.TRUE.,FINPAT,INPIC,'U',IFORM,NSAM,NROW,NSLICE,
     &              MAXIM,'OUTPUT AVERAGE ',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999

C       INVERSE FOURIER
        INS = -1
        CALL  FMRS_2(CNEW,NSAM,NROW,INS)

        CALL WRITEV(INPIC,CNEW,LSD,NROW, NSAM,NROW,NSLICE)

5       CLOSE(INPIC)
        CLOSE(NDOC)

9999    IF (ALLOCATED(CNEW)) DEALLOCATE(CNEW)

        END

C++*********************************************************************
C
C SAQB_P.F
C                  OPFILEC                         FEB  03 ARDEAN LEITH
C
C **********************************************************************
C
C--*********************************************************************

         SUBROUTINE  SAQB_P(BUF,LSD,NSAM,NROW,NIMA,NGRP,JACUP,
     &                      CNEW,FINPAT,NLET,DOCFIL,FINPIC)

         INCLUDE 'CMBLOCK.INC'

         CHARACTER (LEN=*) :: FINPAT,DOCFIL,FINPIC

         REAL, ALLOCATABLE :: A(:,:),B(:,:),C(:,:),D(:,:)
         REAL, ALLOCATABLE :: X(:,:,:)
         REAL              :: CNEW(LSD,NROW)
         DOUBLE PRECISION  :: SOLD,SNEW,EAV,ENE(NIMA),ENFR_2
         REAL              :: DLIST(7),BUF(1024)
         REAL              :: DIST(NIMA),SHIFT(2,NIMA)
         LOGICAL           :: CKOT(NIMA)
         LOGICAL*1         :: CH_ANG

         DATA              INPIC/77/,NDOC/55/,NDOUT/56/

         ALLOCATE (A(LSD,NROW),
     &             B(LSD,NROW),
     &             C(LSD,NROW), 
     &             D(LSD,NROW), 
     &             X(LSD,NROW,NIMA), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            MWANT =  4*LSD*NROW + LSD*NROW*NIMA 
            CALL ERRT(46,'SAQB_P; A,...',MWANT)
            RETURN
         ENDIF

         NSNR = LSD * NROW
         SOLD = 8.0D0 * DATAN(1.0D0) / 360.0D0
         REWIND NDOC

         K1  = 1
         KT1 = 1
         KT2 = 1
         IMI = 0

778      LERR = -1
         CALL  UNSAV(DOCFIL,KT1,NDOC,K1,DLIST,4,LERR,KT2)
         IF (LERR .NE. 0)  GOTO  788

         K1 = K1+1
         IF (IFIX(DLIST(4)) .NE. NGRP)  GOTO  778
         IMI = IMI+1

C        AN(IMI)=DLIST(2)
         L = DLIST(1)

         CALL FILGET(FINPAT,FINPIC,NLET,L,INTFLAG)

         MAXIM = 0
         CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,NSAMT,NROWT,NSL,
     &                   MAXIM,'DUMMY',.FALSE.,IRTFLG)
          
         IF (IRTFLG .NE. 0)  THEN
            DEALLOCATE (A,B,C,D,X)
            RETURN
         ENDIF

         DO K2=1,NROW
            CALL  REDLIN(INPIC,A(1,K2),NSAM,K2)
         ENDDO
         CLOSE(INPIC)
         CALL  RTQ_Q(A,X(1,1,IMI),LSD,NSAM,NROW,DLIST(2))
         INS = 1

         CALL  FMRS_2(X(1,1,IMI),NSAM,NROW,INS)
         IF (INS .EQ. 0)  THEN
            CALL  ERRT(38,'AP SA',NE)
            DEALLOCATE (A,B,C,D,X)
            RETURN
         ENDIF
         GOTO  778
788      CONTINUE

c$omp parallel do private(imi)
         DO    IMI=1,NIMA
            ENE(IMI)=ENFR_2(X(1,1,IMI),LSD,NSAM,NROW)
         ENDDO

C        BUILD FIRST AVERAGE

C        TWO ESTIMATION OF INITIAL AVERAGE ARE USED
C        ONLY ONE !!  11/06/91

         DO    ICR=1,1

C           DIST  IS USED HERE FOR THE RANDOM CHOOSING OF IMAGES

c$omp       parallel do private(imi)
            DO    IMI=1,NIMA
                DIST(IMI) = 0.0
            ENDDO

            CALL RANDOM_NUMBER(CIID)
            IMI = MIN0(NIMA,MAX0(1,INT(CIID*NIMA+0.5)))

            DIST(IMI)    = 1.0
            CKOT(IMI)    = .FALSE.
            SHIFT(1,IMI) = 0.0
            SHIFT(2,IMI) = 0.0
            CALL  COP(X(1,1,IMI),C,NSNR)

            DO    KTN=2,NIMA

804            CALL RANDOM_NUMBER(CIID) 

               M = MIN0(NIMA,MAX0(1,INT(CIID*(NIMA-KTN+1)+0.5)))

               IMI = 0
               DO  I=1,NIMA
                  IF (DIST(I) .NE. 1.0)  THEN
                     IMI = IMI+1
                     IF (IMI .EQ. M)  GOTO  810
                  ENDIF
               ENDDO 
               GOTO  804

810            IMI       = I
               DIST(IMI) = 1.0

               LSC = NSAM+2-MOD(NSAM,2)
               CALL CCRS_2(C,X(1,1,IMI),D, LSC,NSAM,NROW)

               CALL FINDMX(D,LSD,NSAM,NROW,CMX1,SX1,SY1,JACUP)

               CALL CR180_2(C,X(1,1,IMI),D,LSD/2,NSAM,NROW)
               CALL FINDMX(D,LSD,NSAM,NROW,CMX2,SX2,SY2,JACUP)

               SX2 = -SX2
               SY2 = -SY2

               IF(CMX1.GE.CMX2)  THEN
                  CKOT(IMI)    = .FALSE.
                  SHIFT(1,IMI) = SX1
                  SHIFT(2,IMI) = SY1
                  CALL SHFC_2(X(1,1,IMI),A,LSD/2,NSAM,NROW,SX1,SY1)
                  CALL UPDTF(C,A,NSNR,KTN)
               ELSE
                  CKOT(IMI)    = .TRUE.
                  SHIFT(1,IMI) = SX2
                  SHIFT(2,IMI) = SY2
                  CALL  SH180_2(X(1,1,IMI),A,LSD/2,NSAM,NROW,SX2,SY2)
                  CALL  UPDTF(C,A,NSNR,KTN)
               ENDIF
            ENDDO
            SOLD=ENFR_2(C,LSD,NSAM,NROW)
            WRITE(NOUT,2055)  ICR,SOLD
2055        FORMAT('  Initial estimation #',I1,
     &             '   Squared sum=',1PE12.5)

C           WRITE(NOUT,2001) (CKOT(J),SHIFT(1,J),SHIFT(2,J),J=1,NIMA)
2001        FORMAT(4(1X,L1,2(1X,F8.3)))

            IF (ICR.EQ.1)  THEN
               CALL COP(C,CNEW,NSNR)
            ELSE
               LSC = NSAM+2-MOD(NSAM,2)
               CALL CCRS_2(C,CNEW,D, LSC,NSAM,NROW)

               CALL FINDMX(D,LSD,NSAM,NROW,CMX1,SX1,SY1,JACUP)

               CALL CR180_2(C,CNEW,D,LSD/2,NSAM,NROW)

               CALL FINDMX(D,LSD,NSAM,NROW,CMX2,SX2,SY2,JACUP)

               SX2 = -SX2
               SY2 = -SY2

               IF (CMX1 .GE. CMX2)  THEN
                  CALL SHFM_2(CNEW,LSD/2,NSAM,NROW,SX1,SY1)
                  CALL UPDTF(C,CNEW,NSNR,2)
               ELSE
                  CALL SH180_2(CNEW,A,LSD/2,NSAM,NROW,SX2,SY2)
                  CALL UPDTF(C,A,NSNR,2)
               ENDIF

               SOLD = ENFR_2(C,LSD,NSAM,NROW)
               WRITE(NOUT,2065)  SOLD
2065           FORMAT('  Initial average.  Squared sum=',1PE12.5)

c$omp          parallel do private(imi)
               DO IMI=1,NIMA
                  SHIFT(1,IMI) = 0.0
                  SHIFT(2,IMI) = 0.0
                  CKOT(IMI)    = .TRUE.
               ENDDO
            ENDIF
         ENDDO

C       ITERATIONS TO GET BETTER APPROXIMATION

         ITER = 0
901      CONTINUE

         ITER   = ITER+1
         CH_ANG = .FALSE.

         DO    IMI=1,NIMA
C           Subtract from the average ...
            IF (CKOT(IMI))  THEN
               CALL SH180_2
     &         (X(1,1,IMI),B,LSD/2,NSAM,NROW,SHIFT(1,IMI),SHIFT(2,IMI))

            ELSE
               CALL SHFC_2
     &         (X(1,1,IMI),B,LSD/2,NSAM,NROW,SHIFT(1,IMI),SHIFT(2,IMI))
            ENDIF

            CALL SUBAF(C,B,A,NSNR,NIMA)

            LSC = NSAM+2-MOD(NSAM,2)
            CALL CCRS_2(A,X(1,1,IMI),D, LSC,NSAM,NROW)

            CALL FINDMX(D,LSD,NSAM,NROW,CMX1,SX1,SY1,JACUP)

            CALL CR180_2(A,X(1,1,IMI),D,LSD/2,NSAM,NROW)
            CALL FINDMX(D,LSD,NSAM,NROW,CMX2,SX2,SY2,JACUP)

            SX2 = -SX2
            SY2 = -SY2

            IF (CMX1 .GE. CMX2)  THEN
               IF (CKOT(IMI)
     &            .OR.SHIFT(1,IMI).NE.SX1.OR.SHIFT(2,IMI).NE.SY1) THEN
                  CH_ANG       = .TRUE.
                  CKOT(IMI)    = .FALSE.
                  SHIFT(1,IMI) = SX1
                  SHIFT(2,IMI) = SY1

                  CALL  SHFC_2(X(1,1,IMI),B,LSD/2,NSAM,NROW,SX1,SY1)
               ENDIF

               CALL  UPDTF(CNEW,B,NSNR,IMI)
            ELSE
               IF (.NOT. CKOT(IMI)
     &         .OR.SHIFT(1,IMI).NE.SX2.OR.SHIFT(2,IMI).NE.SY2) THEN
                  CH_ANG       = .TRUE.
                  CKOT(IMI)    = .TRUE.
                  SHIFT(1,IMI) = SX2
                  SHIFT(2,IMI) = SY2

                  CALL  SH180_2(X(1,1,IMI),B,LSD/2,NSAM,NROW,SX2,SY2)
               ENDIF

               CALL UPDTF(CNEW,B,NSNR,IMI)
            ENDIF
         ENDDO
         SNEW = ENFR_2(CNEW,LSD,NSAM,NROW)

         WRITE(NOUT,2066)  ITER,SNEW
2066     FORMAT('  Iteration: ',I4,'   Squared sum:',1PE12.5)

         IF (SNEW .GE. SOLD .AND. CH_ANG)  THEN
            CALL  COP(CNEW,C,NSNR)
            SOLD = SNEW
            GOTO  901
         ENDIF

C sep 12   WRITE(NOUT,2001) (CKOT(J),SHIFT(1,J),SHIFT(2,J),J=1,NIMA)

         REWIND NDOC
         K1   = 1
         KT1  = 1
         KT2  = 1
         NIM  = 0
978      LERR = -1
         CALL  UNSAV(DOCFIL,KT1,NDOC,K1,DLIST,4,LERR,KT2)
         IF (LERR.NE.0)  GOTO  988

         K1 = K1+1
         IF(IFIX(DLIST(4)).NE.NGRP)  GOTO  978

         NIM=NIM+1
         DO    I=5,2,-1
            DLIST(I)=DLIST(I-1)
         ENDDO

C        CALCULATE DISTANCES ...
         IF (CKOT(NIM))  THEN
            CALL SH180_2
     &      (X(1,1,NIM),B,LSD/2,NSAM,NROW,SHIFT(1,NIM),SHIFT(2,NIM))

            CALL SUBAF(C,B,A,NSNR,NIMA)
            CALL CR180_2(A,X(1,1,NIM),D,LSD/2,NSAM,NROW)
            CALL FINDMX(D,LSD,NSAM,NROW,CMX1,SX2,SY2,JACUP)
         ELSE
            CALL SHFC_2
     &      (X(1,1,NIM),B,LSD/2,NSAM,NROW,SHIFT(1,NIM),SHIFT(2,NIM))

            CALL SUBAF(C,B,A,NSNR,NIMA)

            LSC = NSAM+2-MOD(NSAM,2)
            CALL CCRS_2(A,X(1,1,NIM),D, LSC,NSAM,NROW)

            CALL FINDMX(D,LSD,NSAM,NROW,CMX1,SX1,SY1,JACUP)
         ENDIF

         EAV  = ENFR_2(A,LSD,NSAM,NROW)
         SNEW = EAV+ENE(NIM)-2*CMX1

         DLIST(1) = NIM
         IF(CKOT(NIM))  DLIST(3) = DLIST(3)+180.0
         DLIST(5) = SHIFT(1,NIM)
         DLIST(6) = SHIFT(2,NIM)
         DLIST(7) = SNEW
         CALL  SAVD(NDOUT,DLIST,7,IRTFLG)
         GOTO  978

988      CALL  SAVDC
         CLOSE(NDOUT)
         DEALLOCATE(A,B,C,D,X)
         END


C++*********************************************************************
C
C    SAQB_F.F 
C              OPFILEC                             FEB  03 ARDEAN LEITH
C
C **********************************************************************
C
C--*********************************************************************

         SUBROUTINE SAQB_F(BUF,LSD,NSAM,NROW,NIMA,NGRP,JACUP,
     &                      CNEW,FINPAT,NLET,DOCFIL,FINPIC)

         INCLUDE 'CMBLOCK.INC'

         CHARACTER (LEN=*) :: FINPAT,DOCFIL,FINPIC

         REAL, ALLOCATABLE :: A(:,:),B(:,:),C(:,:),D(:,:)
         REAL, ALLOCATABLE ::X(:,:,:)

         DIMENSION         CNEW(LSD,NROW)
         DOUBLE PRECISION  SOLD,SNEW,EAV,ENE(NIMA),ENFR_2
         DIMENSION         DLIST(7),BUF(1024)
         DIMENSION         DIST(NIMA),SHIFT(2,NIMA)
         LOGICAL           CKOT(NIMA)
         LOGICAL*1         CH_POS

         DATA              IOUTMP/77/
         DATA              INPIC/77/,NDOC/55/,NDOUT/56/

         ALLOCATE (A(LSD,NROW),
     &             B(LSD,NROW),
     &             C(LSD,NROW), 
     &             D(LSD,NROW), 
     &             X(LSD,NROW,NIMA), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            MWANT =  4*LSD*NROW + LSD*NROW*NIMA 
            CALL ERRT(46,'SAQB_P; A,...',MWANT)
            RETURN
         ENDIF

         NSNR = LSD*NROW
         SOLD = 8.0D0*DATAN(1.0D0)/360.0D0
         REWIND NDOC
         K1   = 1
         KT1  = 1
         KT2  = 1
         IMI  = 0
778      LERR = -1
         CALL  UNSAV(DOCFIL,KT1,NDOC,K1,DLIST,4,LERR,KT2)
         IF (LERR .NE. 0)  GOTO  788

         K1 = K1+1
         IF (IFIX(DLIST(4)) .NE. NGRP)  GOTO  778
         IMI=IMI+1
C        AN(IMI)=DLIST(2)
         L=DLIST(1)

         CALL FILGET(FINPAT,FINPIC,NLET,L,INTFLAG)

         MAXIM = 0
         CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',ITYPE,NSAMT,NROWT,NSL,
     &                   MAXIM,'DUMMY',.FALSE.,IRTFLG)
          
         IF (IRTFLG .NE. 0)  THEN
            CALL ERRT(4,'SAQB_F',NE)
            DEALLOCATE (A,B,C,D,X)
            RETURN
         ENDIF

         DO K2=1,NROW
            CALL  REDLIN(INPIC,A(1,K2),NSAM,K2)
         ENDDO
         CLOSE(INPIC)

         CALL  RTQ_Q(A,X(1,1,IMI),LSD,NSAM,NROW,DLIST(2))

C        FORWARD FFT
         INS = 1
         CALL FMRS_2(X(1,1,IMI),NSAM,NROW,INS)
         IF (INS .EQ. 0)  THEN
            CALL  ERRT(38,'AP SA',NE)
            DEALLOCATE (A,B,C,D,X)
            RETURN
         ENDIF

         GOTO  778

788      CONTINUE

c$omp    parallel do private(imi)
         DO IMI=1,NIMA
            ENE(IMI)=ENFR_2(X(1,1,IMI),LSD,NSAM,NROW)
         ENDDO

C        BUILD FIRST AVERAGE
C        TWO ESTIMATION OF INITIAL AVERAGE ARE USED
C        ONLY ONE !!!  11/06/91

         DO ICR=1,1

C        DIST  IS USED HERE FOR THE RANDOM CHOOSING OF IMAGES

c$omp       parallel do private(imi)
            DO IMI=1,NIMA
               DIST(IMI) = 0.0
            ENDDO

            CALL RANDOM_NUMBER(CIID)

            IMI = MIN0(NIMA,MAX0(1,INT(CIID*NIMA+0.5)))

            DIST(IMI)    = 1.0
            CKOT(IMI)    = .FALSE.
            SHIFT(1,IMI) = 0.0
            SHIFT(2,IMI) = 0.0

            CALL COP(X(1,1,IMI),C,NSNR)

            DO KTN=2,NIMA

804            CALL RANDOM_NUMBER(CIID)
 
               M   = MIN0(NIMA,MAX0(1,INT(CIID*(NIMA-KTN+1)+0.5)))

               IMI = 0
               DO I=1,NIMA
                  IF (DIST(I) .NE .1.0)  THEN
                     IMI = IMI+1
                     IF (IMI .EQ. M)  GOTO  810
                  ENDIF
               ENDDO
               GOTO  804

810            IMI       = I
               DIST(IMI) = 1.0

               LSC = NSAM+2-MOD(NSAM,2)
               CALL CCRS_2(C,X(1,1,IMI),D, LSC,NSAM,NROW)

               CALL FINDMX(D,LSD,NSAM,NROW,CMX1,SX1,SY1,JACUP)

               SHIFT(1,IMI) = SX1
               SHIFT(2,IMI) = SY1

               CALL SHFC_2(X(1,1,IMI),A,LSD/2,NSAM,NROW,SX1,SY1)
               CALL UPDTF(C,A,NSNR,KTN)
            ENDDO

            SOLD = ENFR_2(C,LSD,NSAM,NROW)

            WRITE(NOUT,2055)  ICR,SOLD
2055        FORMAT('  INITIAL ESTIMATION #',I1, 
     &             '   SQUARED SUM:',1PE12.5)

C           WRITE(NOUT,2001)  (SHIFT(1,J),SHIFT(2,J),J=1,NIMA)
2001        FORMAT(4(2X,2(1X,F8.3)))

            IF (ICR.EQ.1)  THEN
               CALL  COP(C,CNEW,NSNR)
            ELSE
               LSC = NSAM+2-MOD(NSAM,2)
               CALL CCRS_2(C,CNEW,D, LSC,NSAM,NROW)

               CALL FINDMX(D,LSD,NSAM,NROW,CMX1,SX1,SY1,JACUP)
               CALL SHFM_2(CNEW,LSD/2,NSAM,NROW,SX1,SY1)
               CALL UPDTF(C,CNEW,NSNR,2)

               SOLD = ENFR_2(C,LSD,NSAM,NROW)
               WRITE(NOUT,2065)  SOLD
2065           FORMAT('  INITIAL AVERAGE.  SQUARED SUM=',1PE12.5)

c$omp          parallel do private(imi)
               DO IMI=1,NIMA
                  SHIFT(1,IMI)=0.0
                  SHIFT(2,IMI)=0.0
               ENDDO
            ENDIF
         ENDDO

C        ITERATIONS TO GET BETTER APPROXIMATION

         ITER=0
901      CONTINUE

         ITER=ITER+1
         CH_POS=.FALSE.


         DO    IMI=1,NIMA
C           SUBTRACT FROM THE AVERAGE ...
            CALL  SHFC_2
     &      (X(1,1,IMI),B,LSD/2,NSAM,NROW,SHIFT(1,IMI),SHIFT(2,IMI))
            CALL  SUBAF(C,B,A,NSNR,NIMA)

            LSC = NSAM+2-MOD(NSAM,2)
            CALL  CCRS_2(A,X(1,1,IMI),D, LSC,NSAM,NROW)

            CALL  FINDMX(D,LSD,NSAM,NROW,CMX1,SX1,SY1,JACUP)

            IF (SHIFT(1,IMI).NE.SX1.OR.SHIFT(2,IMI).NE.SY1) THEN
               CH_POS       = .TRUE.
               SHIFT(1,IMI) = SX1
               SHIFT(2,IMI) = SY1
               CALL  SHFC_2(X(1,1,IMI),B,LSD/2,NSAM,NROW,SX1,SY1)
            ENDIF
            CALL  UPDTF(CNEW,B,NSNR,IMI)
         ENDDO

         SNEW = ENFR_2(CNEW,LSD,NSAM,NROW)

         WRITE(NOUT,2066)  ITER,SNEW
2066     FORMAT('  ITERATION #',I4,'   SQUARED SUM=',1PE12.5)
         IF (SNEW .GE. SOLD .AND. CH_POS)  THEN
            CALL  COP(CNEW,C,NSNR)
            SOLD = SNEW
            GOTO  901
         ENDIF

c sep12al         WRITE(NOUT,2001)  (SHIFT(1,J),SHIFT(2,J),J=1,NIMA)

         REWIND NDOC
         K1   = 1
         KT1  = 1
         KT2  = 1
         NIM  = 0
978      LERR = -1
         CALL  UNSAV(DOCFIL,KT1,NDOC,K1,DLIST,4,LERR,KT2)
         IF (LERR.NE.0)  GOTO  988

         K1 = K1+1
         IF (IFIX(DLIST(4)) .NE. NGRP)  GOTO  978
         NIM = NIM+1

         DO I=5,2,-1
            DLIST(I) = DLIST(I-1)
         ENDDO

C        CALCULATE DISTANCES ...
         CALL SHFC_2
     &     (X(1,1,NIM),B,LSD/2,NSAM,NROW,SHIFT(1,NIM),SHIFT(2,NIM))

         CALL SUBAF(C,B,A,NSNR,NIMA)

         LSC = NSAM+2-MOD(NSAM,2)
         CALL CCRS_2(A,X(1,1,NIM),D, LSC,NSAM,NROW)

         CALL FINDMX(D,LSD,NSAM,NROW,CMX1,SX1,SY1,JACUP)

         EAV  = ENFR_2(A,LSD,NSAM,NROW)
         SNEW = EAV+ENE(NIM)-2*CMX1

         DLIST(1) = NIM
         DLIST(5) = SHIFT(1,NIM)
         DLIST(6) = SHIFT(2,NIM)
         DLIST(7) = SNEW
         CALL SAVD(NDOUT,DLIST,7,IRTFLG)
         GOTO  978

988      CALL SAVDC
         CLOSE(NDOUT)

         DEALLOCATE(A,B,C,D,X)

         END
