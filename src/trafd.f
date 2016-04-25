C++*********************************************************************
C
C TRAFD.F   USED OPFILE                           NOV 00 ArDean Leith
C           OPFILEC                               FEB 03 ArDean Leith
C           REWORKED 'TF L FLIP'                  NOV 15 ArDean Leith                           
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
C TRAFD(LUN)
C
C  PURPOSE: GENERATE THE PHASE CONTRAST TRANSFER FUNCTION  FOR
C           BRIGHT-FIELD ELECTRON MICROSCOPY. THIS OPERATION SHOWS
C            CTF IN REAL, DISPLAYABLE FORM TO AN IMAGE FILE.
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE TRAFD(LUN)

         IMPLICIT NONE
         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC' 
 
         INTEGER               :: LUN

         REAL                  :: B(NBUFSIZ)  ! From: cmlimit.inc

         CHARACTER(LEN=MAXNAM) :: FILNAM
         CHARACTER             :: NULL = CHAR(0)
         CHARACTER             :: ANS

         INTEGER               :: NDIM,NX,NY,NZ,IRTFLG   
         LOGICAL               :: WANT_AST,WANT_GEH,WANT_SIGN
         LOGICAL               :: WANT_SPFREQ,WANT_PIXSIZ

         INTEGER               :: NLET,MAXIM,NCHAR,IE,NS1,I,K
   
         REAL                  :: LAMBDA,FMAXSPFREQ,CS,DZ,Q,DS,DZA,AZZ
         REAL                  :: ACR,GEH,SIGN,SC,AK,AZ,PI,AZR,DZZ
         REAL                  :: FDUM

         CALL FILERD(FILNAM,NLET,NULL,'CTF OUTPUT',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

C        GET COMMON TF INPUTS
         NDIM        =  1
         WANT_AST    = .TRUE.
         WANT_GEH    = .TRUE.   
         WANT_SIGN   = .FALSE.
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

         IFORM = 1
         NY    = NX
         NZ    = 1
         MAXIM = 0
         CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'U',IFORM,NX,NY,NZ,
     &                   MAXIM,' ',.TRUE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &       'DIFFRACTOGRAM / ENVELOPE / STRAIGHT (D/E/S)',NULL,IRTFLG)

         IE = 0
         IF (ANS == 'E') IE = 1

         IF (GEH .NE. 0.0) GEH = 1. / GEH**2

         SC  = FMAXSPFREQ / FLOAT(NX / 2)
         ACR = ATAN(ACR / (1.0 - ACR))
         CS  = CS * 1.E7

         NS1 = NX / 2 + 1

         DO  I=1,NY

            DO  K=1,NX

               AK = SQRT(FLOAT(K-NS1)**2 + FLOAT(I-NS1)**2) * SC
               AZ = PI / 2.0

               AZ  = ATAN2(FLOAT(I-NS1), FLOAT(K-NS1)) + PI / 2.0
               AZR = AZZ * (PI / 180.)
               DZZ = DZ + DZA / 2 * SIN(2 * (AZ-AZR))

               CALL TFD(B(K),CS,DZZ,LAMBDA,Q,DS,IE,AK,ACR,GEH)

               IF (ANS .NE. 'S') B(K) = B(K) * B(K)

            ENDDO

            CALL WRTLIN(LUN,B,NX,I)

         ENDDO

         END
