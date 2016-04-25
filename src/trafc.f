C++*********************************************************************
C
C TRAFC.F
C          INCLUDE COMPLETE SIGN REVERSAL AND INTELLIGIBLE
C          GAUSSIAN PARAMETER                     9/5/94  JF
C          USED OPFILE                            NOV  00 ARDEAN LEITH
C          RECTANGULAR OUTPUTS                    OCT  01 BILL BAXTER
C          OPFILEC                                FEB  03 ARDEAN LEITH
C          WANT_FLIP                              MAR  14 ARDEAN LEITH
C          REFACTORED                             NOV  15 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2015  Health Research Inc.,                         *
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
C TRAFC(LUN,WANT_FLIP)
C
C OPERATIONS:   'TF C'  = STRAIGHT
C               'TF CT' = PHASE FLIPPING
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE TRAFC(LUN,WANT_FLIP)

         IMPLICIT NONE
         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC' 
 
         INTEGER               :: LUN   
         LOGICAL               :: WANT_FLIP

         COMPLEX               :: B(NBUFSIZ)  ! from cmlimit.inc

         CHARACTER(LEN=MAXNAM) :: FILNAM
         CHARACTER             :: NULL = CHAR(0)

         REAL                  :: CS,DZ,LAMBDA,FMAXSPFREQ, Q,DS
         REAL                  :: DZA,AZZ,ACR,GEH,FDUM
         REAL                  :: SIGN,SCX,SCY,AK,AZ,AZR,DZZ,TF

         INTEGER               :: I,KX
         INTEGER               :: NLET,LSM,IXC,IYC,NZ,MAXIM,IE,K,KY
   
         INTEGER               :: NDIM,NX,NY,IRTFLG   
         LOGICAL               :: WANT_AST,WANT_GEH,WANT_SIGN
         LOGICAL               :: WANT_SPFREQ,WANT_PIXSIZ

         DOUBLE PRECISION, PARAMETER :: QUADPI=3.1415926535897932384626

         IF (WANT_FLIP) THEN
            CALL FILERD(FILNAM,NLET,NULL,'PHASE FLIPPING CTF OUTPUT',
     &                  IRTFLG)
         ELSE
            CALL FILERD(FILNAM,NLET,NULL,'CTF OUTPUT',IRTFLG)
         ENDIF
         IF (IRTFLG .NE. 0) RETURN

C        GET COMMON TF INPUTS
         NDIM        =  2
         WANT_AST    = .TRUE.
         WANT_GEH    = .NOT. WANT_FLIP   
         WANT_SIGN   = .TRUE.
         WANT_SPFREQ = .TRUE.     ! ASK FOR SPFREQ
         WANT_PIXSIZ = .FALSE.    ! DO NOT ASK FOR PIXEL SIZE

         CALL GET_TF_INPUT(CS,DZ,LAMBDA,
     &                NDIM, NX, NY,
     &                WANT_SPFREQ,FMAXSPFREQ,
     &                WANT_PIXSIZ,FDUM,
     &                Q, DS,
     &                WANT_AST, DZA, AZZ,
     &                WANT_GEH, ACR, GEH,
     &                WANT_SIGN, SIGN,
     &                IRTFLG) 
         IF (IRTFLG .NE. 0) RETURN

         IF (MOD(NX,2) == 0)  THEN
            IFORM = -12
            LSM   = NX + 2
         ELSE
            IFORM = -11
            LSM   = NX + 1
         ENDIF
         IXC  = NX / 2 + 1

         IF (NY == 0)  THEN
            NY   = NX
            IYC  = IXC
         ELSE
            IYC  = NY / 2 + 1
         ENDIF

         NZ     = 1
         MAXIM  = 0
         IRTFLG = 0
         CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'U',IFORM,LSM,NY,NZ,
     &                MAXIM,' ',.TRUE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

C        IE  = 0 SELECTS STRAIGHT TRANSFER FUNCTION IN SUBROUTINE TFD
         IE  = 0

         IF (WANT_FLIP) THEN
C           WANT PHASE FLIPPING TRANSFER
 
C           IGNORE THE GAUSSIAN ENVELOPE FUNCTION
C           this is almost = 1/10000**2 setting? al
            GEH = 0.0

         ELSE
C           WANT FULL CTF TRANSFER FUNCTION

            IF (GEH .NE. 0.0) THEN
               GEH = 1. / GEH**2
            ENDIF
         ENDIF


C        sc  = fmaxspfreq / float(nx / 2) ! for ndim=1
         SCX = 2.0 / NX
         SCY = 2.0 / NY

         ACR = ATAN(ACR / (1.0 - ACR))
         CS  = CS * 1.E7

         DO  K=1,NY
            KY = K-1
            IF (K > IYC) KY = KY-NY

            DO  I=1,LSM,2
               KX = (I-1) / 2

C              Changed AK to handle rectangular images
C              AK = SQRT(FLOAT(KY)**2 + FLOAT(KX)**2)*SC

               AK = FMAXSPFREQ * SQRT((KX*SCX)**2 + (KY*SCY)**2)

C              AZ = QUADPI / 2.
               IF (KX .NE. 0) THEN
                  AZ = ATAN2(FLOAT(KY),FLOAT(KX)) + QUADPI / 2.
               ELSE
                  AZ = QUADPI / 2.
               ENDIF
               
               AZR = AZZ * (QUADPI/180.)
               DZZ = DZ + DZA / 2 * SIN(2 *(AZ - AZR))

               CALL TFD(TF,CS,DZZ,LAMBDA,Q,DS,IE,AK,ACR,GEH)

               IF (WANT_FLIP) THEN
C                 FOR PHASE FLIPPING
                  IF (TF >= 0.0) THEN
                     B(KX+1) = CMPLX( 1.0,0.0) * SIGN
                  ELSE
                     B(KX+1) = CMPLX(-1.0,0.0) * SIGN
                  ENDIF

               ELSE
C                 FOR FULL CTF CORRECTION
                  B(KX+1) = CMPLX(TF*SIGN, 0.0)

               ENDIF
            ENDDO

            CALL WRTLIN(LUN,B,LSM,K)
         ENDDO

         END

