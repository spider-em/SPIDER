C++*********************************************************************
C
C HISTE.F                               USED OPFILE NOV 00 ARDEAN LEITH
C                        REWORKED MEMORY ALLOCATION OCT 05 ARDEAN LEITH
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
C  HISTE(UNUSED)
C
C  PURPOSE: Finds the linear transformation (applied to pixels) 
C           which fits the histogram of the image file to the 
C           histogram of the reference file. 
C
C--*********************************************************************

         SUBROUTINE HISTE(UNUSED)

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC' 
 
         CHARACTER(LEN=MAXNAM)   ::  FILNAM
         COMMON /COMMUN/FILNAM

         REAL, ALLOCATABLE, DIMENSION(:)    :: QK1,QK2
         LOGICAL, ALLOCATABLE, DIMENSION(:) :: QK6

         DATA  LUN1,LUN2,LUN3,LUN4/8,9,10,11/

         MAXIM = 0
         CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &                   MAXIM,'REFERENCE IMAGE',.FALSE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         MAXIM = 0
         CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',IFORM,
     &                NSAM1,NROW1,NSLICE1,
     &                MAXIM,'IMAGE TO BE CORRECTED~',.FALSE.,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 999

         MAXIM = 0
         CALL OPFILEC(0,.TRUE.,FILNAM,LUN4,'O',IFORM,
     &          NSAM2,NROW2,NSLICE2,
     &          MAXIM,'MASK FOR THE IMAGE TO CORRECTED~',.FALSE.,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 999

         IF (NSAM1   .NE. NSAM2 .OR.
     &       NROW1   .NE. NROW2 .OR. 
     &       NSLICE1 .NE. NSLICE2) THEN
            CALL ERRT(2,'CE FIT  ',NE)
            GOTO 999
         ENDIF

         MAXIM = 0
         CALL OPFILEC(LUN2,.TRUE.,FILNAM,LUN3,'U',IFORM,
     &             NSAM1,NROW1,NSLICE1,
     &             MAXIM,'OUTPUT',.FALSE.,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 999
 
         NSR    = NSAM  * NROW * NSLICE
         NSR1   = NSAM1 * NROW1 * NSLICE1
         LENH   = MIN0(NSR/16, NSR1/16, 256)

C        OLD MEMORY CALCULATION FOR COMMON BLOCK
         K1=1              ! NSAM  * NROW * NSLICE
         K2=K1+NSR         ! NSAM1 * NROW1 * NSLICE1
         K4=K2+NSR1        ! LENH * 3
         K5=K4+LENH*3      ! LENH * 3
         K6=K5+LENH*3      ! (NSR1 + 3) / 4

         ALLOCATE(QK1(NSR), QK2(NSR1), QK6(NSR1+3),STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            NEEDED = NSR + NSR1+ NSR1+3
            CALL ERRT(101,'UNABLE TO ALLOCATE MEMORY IN HISTE',NEEDED)
            GOTO 999
         ENDIF

C        READ REF. IMAGE INTO QK1
         CALL REDVOL(LUN1,NSAM,NROW,1,NSLICE,QK1,IRTFLG)

C        READ MASK FOR IMAGE TO CORRECT INTO QK2
         CALL REDVOL(LUN4,NSAM1,NROW1,1,NSLICE1,QK2,IRTFLG)

C        MAKE LOGICAL MASK
         QK6 = (QK2 .NE. 0.0)

C        READ IMAGE TO CORRECT INTO QK2
         CALL REDVOL(LUN2,NSAM1,NROW1,1,NSLICE1,QK2,IRTFLG)

         ITRMAX = 100
         CALL HISTC2(QK1,QK2,QK6,NSR,NSR1,LENH,ITRMAX,NOUT)

C        WRITE OUTPUT IMAGE 
         CALL WRTVOL(LUN3,NSAM1,NROW1,1,NSLICE1,QK2,IRTFLG)

999      CLOSE(LUN1)
         CLOSE(LUN2)
         CLOSE(LUN3)
         CLOSE(LUN4)

         IF (ALLOCATED(QK1)) DEALLOCATE(QK1)
         IF (ALLOCATED(QK2)) DEALLOCATE(QK2)
         IF (ALLOCATED(QK6)) DEALLOCATE(QK6)

         RETURN
         END



C--*********************************************************************

         SUBROUTINE  HISTC2(QK1,QK2,QK6,N,  NSR1,LENH,ITRMAX,NOUT)

         REAL      ::  QK1(N),QK2(NSR1)
         LOGICAL   ::  QK6(NSR1)

         REAL      ::  QK4(3*LENH),QK5(3*LENH)

         REAL      ::  AK(2),P(3,2),Y(3)
         REAL      ::  PR(2),PRR(2),PBAR(2)

         EXTERNAL  FHT2

         XRMI = QK1(1)
         XRMA = XRMI
         AVR  = XRMI
         SR   = XRMI**2
         DO I=2,N
            AVR  = AVR + QK1(I)
            SR   = SR  + QK1(I) ** 2
            XRMI = AMIN1(XRMI,QK1(I))
            XRMA = AMAX1(XRMA,QK1(I))
	 ENDDO

C        ximi=amin1(ximi,xi(i))
C1       xima=amax1(xima,xi(i))

         NT1  = 0
         XIMI = 1.E23
         XIMA = -XIMI
         AVI  = 0.0
         SI   = 0.0

         DO I=1,NSR1
            IF (QK6(I))  THEN
               NT1 = NT1 + 1
               AVI = AVI + QK2(I)
               SI  = SI  + QK2(I)**2
            ENDIF
	 ENDDO

         RXR = XRMA - XRMI
C        rxi = xima - ximi
         AVR = AVR / N
         AVI = AVI / NT1
         SR  = SQRT((SR - N   * AVR * AVR) / (N-1))
         SI  = SQRT((SI - NT1 * AVI * AVI) / (NT1-1))

         DO I=1,3*LENH
            QK4(I) = 0
	 ENDDO

         DO I=1,N
            L = INT((QK1(I) - XRMI) / RXR * (LENH-1) + LENH+1)
            QK4(L) = QK4(L) + 1
	 ENDDO

         DO I=1,3*LENH
            QK4(I) = QK4(I) * FLOAT(NT1) / FLOAT(N)
	 ENDDO

         A      = SR/SI
         P(1,1) = 0.9*A
         P(2,1) = A
         P(3,1) = 1.1*A
         B      =  AVR-A*AVI
         P(1,2) = -B
         P(2,2) = B
         P(3,2) = 3*B

         DO I=1,3
            AK(1) = P(I,1)
            AK(2) = P(I,2)
            Y(I)  = FHT2(AK,NSR1,LENH,QK2,QK4,QK5,QK6,RXR,XRMI)
	 ENDDO

         WRITE(NOUT,205)  A,B,Y(2)
205      FORMAT(' The transformation is  A*x + B',/,
     &   ' Initial parameters   A =',1pe12.5,'   B =',1pe12.5,/,
     &   ' Initial chi-square     =',1pe12.5)

         N2  = 2
         EPS = 0.0001

         CALL AMOEBA2(P,Y,N2,EPS,FHT2,ITER,ITRMAX,PR,PRR,PBAR,
     &               NSR1,LENH,QK2,QK4,QK5,QK6,RXR,XRMI)

         WRITE(NOUT,206)  ITER,P(2,1),P(2,2),Y(2)
206      FORMAT(' Minimum was found in ',i3,' iterations.',/,
     &          ' Parameters found     A =',1pe12.5,'   B =',1pe12.5,/,
     &          ' Final   chi-square     =',1pe12.5)

C        do  6  i=1,4
C6          print  203,(p(i,j),j=1,3)
C203     format(3(3x,e12.5))

         DO I=1,NSR1
            QK2(I) = QK2(I) * P(2,1) + P(2,2)
	 ENDDO

         END


C--*********************************************************************

         FUNCTION FHT2(AK,N,LENH,XI,H1,H2,IFP,RXR,XRMI)

         REAL    ::  AK(2),XI(N),H1(3*LENH),H2(3*LENH)
         LOGICAL ::  IFP(N)

         DO I=1,3*LENH
            H2(I) = 0
         ENDDO

         DO  I=1,N
            IF (IFP(I))  THEN
               XN    = XI(I) * AK(1) + AK(2)
               L     = INT((XN - XRMI) / RXR * (LENH - 1) + LENH + 1)
               IF (L.GE.1 .AND. L .LE. 3*LENH)   H2(L) = H2(L) + 1
            ENDIF
         ENDDO

         FHT2 = 0.0

         DO I=1,3*LENH
            FHT2 = FHT2 + (H1(I) - H2(I)) ** 2
         ENDDO

         END

C--*********************************************************************

         SUBROUTINE AMOEBA2(P,Y,NDIM,FTOL,FUNK,ITER,ITMAX,PR,PRR,PBAR,
     &                      NSR1,LENH,QK2,QK4,QK5,QK6,RXR,XRMI)

         INCLUDE 'CMBLOCK.INC'
 
         PARAMETER (ALPHA=1.0,BETA=0.5,GAMMA=2.0)
         DIMENSION P(NDIM+1,NDIM),Y(NDIM+1),PR(NDIM)
         DIMENSION PRR(NDIM),PBAR(NDIM)

         EXTERNAL  FUNK

         MPTS = NDIM+1
         ITER = 0

1        ILO = 1
         IF (Y(1) .GT. Y(2)) THEN
            IHI  = 1
            INHI = 2
         ELSE
            IHI  = 2
            INHI = 1
         ENDIF

         DO I=1,MPTS
            IF (Y(I) .LT. Y(ILO)) ILO=I
            IF (Y(I) .GT. Y(IHI)) THEN
               INHI=IHI
               IHI=I
            ELSE IF (Y(I) .GT. Y(INHI)) THEN
               IF ( I.NE. IHI) INHI=I
            ENDIF
         ENDDO

         RTOL = 2.*ABS(Y(IHI)-Y(ILO))/(ABS(Y(IHI))+ABS(Y(ILO)))
         IF (RTOL .LT. FTOL) RETURN
         IF (ITER .EQ. ITMAX) THEN
            WRITE(NOUT,*) ' Amoeba exceeding maximum iterations.'
            RETURN
         ENDIF
         ITER = ITER + 1
         DO J=1,NDIM
            PBAR(J) = 0.
         ENDDO

         DO I=1,MPTS
            IF (I .NE. IHI) THEN
               DO  J=1,NDIM
                  PBAR(J) = PBAR(J)+P(I,J)
               ENDDO
            ENDIF
         ENDDO

         DO J=1,NDIM
            PBAR(J) = PBAR(J) / NDIM
            PR(J)   = (1. + ALPHA) * PBAR(J) - ALPHA * P(IHI,J)
         ENDDO

         YPR = FUNK(PR,NSR1,LENH,QK2,QK4,QK5,QK6,RXR,XRMI)
C              FHT2(AK,NSR1,LENH,QK2,QK4,QK5,QK6,RXR,XRMI)

         IF (YPR .LE. Y(ILO)) THEN
            DO  J=1,NDIM
               PRR(J) = GAMMA*PR(J)+(1. - GAMMA) * PBAR(J)
            ENDDO

            YPRR = FUNK(PRR,NSR1,LENH,QK2,QK4,QK5,QK6,RXR,XRMI)
C                  FHT2(AK, NSR1,LENH,QK2,QK4,QK5,QK6,RXR,XRMI)

            IF (YPRR .LT. Y(ILO)) THEN
               DO J=1,NDIM
                  P(IHI,J) = PRR(J)
               ENDDO

               Y(IHI)=YPRR
            ELSE
               DO  J=1,NDIM
                  P(IHI,J) = PR(J)
               ENDDO

               Y(IHI) = YPR
            ENDIF

         ELSE IF (YPR .GE. Y(INHI)) THEN
            IF (YPR .LT. Y(IHI)) THEN
               DO  J=1,NDIM
                  P(IHI,J) = PR(J)
               ENDDO

               Y(IHI)=YPR
            ENDIF
            DO J=1,NDIM
               PRR(J) = BETA*P(IHI,J) + (1.-BETA)*PBAR(J)
            ENDDO

            YPRR = FUNK(PRR,NSR1,LENH,QK2,QK4,QK5,QK6,RXR,XRMI)
C                   FHT2(AK,NSR1,LENH,QK2,QK4,QK5,QK6,RXR,XRMI)

            IF (YPRR .LT. Y(IHI)) THEN
               DO J=1,NDIM
                  P(IHI,J) = PRR(J)
               ENDDO

               Y(IHI)=YPRR
            ELSE
               DO  I=1,MPTS
                  IF (I.NE.ILO) THEN
                     DO  J=1,NDIM
                        PR(J)  = 0.5*(P(I,J)+P(ILO,J))
                        P(I,J) = PR(J)
                     ENDDO

                     Y(I) = FUNK(PR,NSR1,LENH,QK2,QK4,QK5,QK6,RXR,XRMI)
                  ENDIF
               ENDDO
            ENDIF
         ELSE
            DO  J=1,NDIM
               P(IHI,J) = PR(J)
            ENDDO

            Y(IHI) = YPR
         ENDIF

         GO TO 1
         END

