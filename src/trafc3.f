C++*********************************************************************
C
C TRAFC3.F
C          INCLUDE COMPLETE SIGN REVERSAL AND INTELLIGIBLE
C          GAUSSIAN PARAMETER                     9/5/94  JF
C          FIXED                                  3/2/95  PP 
C          OPFILEC                                FEB 03  ARDEAN LEITH
C          CMLIMIT & PROMPTS                      NOV 10  ARDEAN LEITH
C          REFACTORED                             NOV 15  ARDEAN LEITH
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
C TRAFC3(LUN,WANT_FLIP)
C 
C
C OPERATION: 'TF C3'  = STRAIGHT
C PURPOSE:   CREATE A CUBIC COMPLEX CTF  FILE THAT CAN BE USED FOR
C            CTF CORRECTION.
C
C OPERATION: 'TF C3'  = STRAIGHT
C PURPOSE:   CREATE A CUBIC BINARY FILE (-1,1) THAT CAN BE USED FOR
C            PHASE FLIPPING CTF CORRECTION OF A VOLUME.
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

         SUBROUTINE TRAFC3(LUN,WANT_FLIP)

         IMPLICIT NONE
         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC'

         INTEGER               :: LUN   
         LOGICAL               :: WANT_FLIP

         COMPLEX               :: B(NBUFSIZ)  ! from cmlimit.inc

         CHARACTER(LEN=MAXNAM) :: FILNAM
         REAL                  :: LAMBDA,FMAXSPFREQ,SIGN
         CHARACTER(LEN=1)      :: NULL = CHAR(0)  

         REAL                  :: CS,DZ,Q,DS,DZA,AZZ,ACR,GEH,SC
         REAL                  :: AK,AZ,AZR,DZZ,TF,FDUM 
         INTEGER               :: NLET,IRTFLG,NDIM,NX,NY,NZ,LSM,MAXIM
         INTEGER               :: IE,NS2,NR2,NL2,K,IZ,J,IY,I,IX
         LOGICAL               :: WANT_AST,WANT_GEH,WANT_SIGN
         LOGICAL               :: WANT_SPFREQ,WANT_PIXSIZ
 
         DOUBLE PRECISION, PARAMETER :: PI = 3.1415926

         CALL FILERD(FILNAM,NLET,NULL,'OUTPUT',IRTFLG)
         IF (IRTFLG == -1) RETURN

C        GET COMMON TF INPUTS
         NDIM        = 1   
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
     &                WANT_AST, DZA,AZZ,
     &                WANT_GEH, ACR, GEH,
     &                WANT_SIGN, SIGN,
     &                IRTFLG) 
         IF (IRTFLG .NE. 0) RETURN

         NY = NX
         NZ = NY
         IF (MOD(NX,2) == 0)  THEN
            IFORM = -22
            LSM   = NX+2
         ELSE
            IFORM = -21
            LSM   = NX+1
         ENDIF

         MAXIM = 0
         CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'U',IFORM,LSM,NY,NZ,
     &               MAXIM,' ',.FALSE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         IF (WANT_FLIP) THEN
C           WANT PHASE FLIPPING TRANSFER
 
C           IGNORE SOURCE SIZE AND DEFOCUS SPREAD
            Q  = 0.0
            DS = 0.0

C           IGNORE THE GAUSSIAN ENVELOPE FUNCTION
            GEH = 0.0

         ELSE
C           WANT FULL CTF TRANSFER FUNCTION

            GEH = 1. / GEH**2
         ENDIF

         SC = FMAXSPFREQ / FLOAT(NX/2)

C        IE  = 0 SELECTS STRAIGHT TRANSFER FUNCTION IN SUBROUTINE TFD
         IE  = 0

         ACR = ATAN(ACR / (1.0 - ACR))
         CS  = CS * 1.E7

         NS2 = NX/2
         NR2 = NY/2
         NL2 = NZ/2
C
         DO K=1,NZ
           IZ = K-1
           IF (IZ > NL2) IZ = IZ-NZ

           DO J=1,NY
              IY = J-1
              IF (IY > NR2) IY = IY-NY

              DO I=1,LSM
                 IX = (I-1)/2

                 AK = SQRT(FLOAT(IX*IX)+FLOAT(IY*IY) + FLOAT(IZ*IZ))*SC
                 IF (AK .NE. 0.0) THEN
                    AZ = 0.0
                 ELSE
                    AZ = PI / 2.
                 ENDIF

                 AZR = AZZ * (PI / 180.)
                 DZZ = DZ + DZA / 2 * SIN(2 * (AZ-AZR))

                 CALL TFD(TF,CS,DZZ,LAMBDA,Q,DS,IE,AK,ACR,GEH)

                  IF (WANT_FLIP) THEN
C                    WANT PHASE FLIPPING TRANSFER FILE
                     IF (TF >= 0.0) THEN
                        B(IX+1) = CMPLX( 1.0,0.0) * SIGN
                     ELSE
                        B(IX+1) = CMPLX(-1.0,0.0) * SIGN
                     ENDIF
                  ELSE
C                    WANT STRAIGHT TRANSFER  FILE
                     B(IX+1) = CMPLX(TF*SIGN,0.0)
                  ENDIF
              ENDDO

              CALL WRTLIN(LUN,B,LSM,J+(K-1)*NY)

           ENDDO
         ENDDO

         END
