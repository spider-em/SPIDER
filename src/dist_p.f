C++*********************************************************************
C
C DIST_P.F
C                MDIM_8                            MAY 13 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C PARAMETERS:
C        NFAC   - NUMBER OF PARAMETERS IN ONE RECORD
C        MAXFA  - ACTUAL NUMBER OF PARAMETERS USED
C        IDK    - IMAGE OR PIXEL IDS                             (RET.)
C        D      -                                                (RET.)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE DIST_P(D,MDIM_8,IDK,NKLA,W,INUM,
     &                     MAXFA,NFAC,LUNF,ITYPE)

         IMPLICIT NONE

         INTEGER * 8      :: MDIM_8
         REAL             :: D(MDIM_8)
         INTEGER          :: NKLA
         REAL             :: W(MAXFA)
         INTEGER          :: IDK(NKLA) 
         INTEGER          :: INUM(MAXFA)
         INTEGER          :: MAXFA,NFAC,LUNF,ITYPE
         
         DOUBLE PRECISION :: Q

         INTEGER          :: mono
         INTEGER          :: K1,K2
         REAL             :: ADUM,BDUM,CDUM,FDUM,FIM

         INTEGER          :: M,I,IDUM,J,KK,L,LL,NDUM

C        AUTOMATIC ARRAYS
         REAL             :: COO(NFAC),COB(NFAC)


C        INTERNAL FUNCTION
         MONO(K1,K2) = MIN(K1,K2)+((MAX(K1,K2)-1)*(MAX(K1,K2)-2)/2)

C        ZERO ARRAY
         D  = 0.0

         REWIND(LUNF)

         IF (ITYPE == 1) THEN
C           UNFORMATTED FILE FOR _SEQ IMAGE DATA
            READ(LUNF) NDUM,NDUM

            DO I=1,NKLA
               READ(LUNF) (COB(J),J=1,NFAC),FIM
               IDK(I) = FIM
            ENDDO
         ELSE
C           FORMATTED FILE FOR IMAGE OR PIXEL COOR.
            READ(LUNF,*) IDUM, IDUM, IDUM, IDUM, IDUM

            DO I=1,NKLA
                READ(LUNF,*) (COB(J),J=1,NFAC),ADUM,BDUM,FIM,CDUM
                IDK(I) = FIM
            ENDDO
         ENDIF

         DO K1=1,NKLA-1
            KK = NKLA-K1+1
            REWIND(LUNF)

            IF (ITYPE == 1) THEN
               READ(LUNF)      ! TO SET POSITION
            ELSE
               READ(LUNF,*) IDUM, IDUM, IDUM, IDUM, IDUM, IDUM
            ENDIF

            DO K2=1,NKLA-K1
               IF (ITYPE == 1) THEN
                  READ(LUNF)   (COO(J),J=1,NFAC),FIM
               ELSE
C                 IMC OR PIX FILE
                  READ(LUNF,*) (COO(J),J=1,NFAC),FDUM,FDUM,FDUM,FDUM
               ENDIF

               Q = 0.0D0
               DO L=1,MAXFA
                  LL = INUM(L)
                  Q  = Q + (DBLE(W(L)) * (COO(LL) - COB(LL)))**2
               ENDDO

               M    = MONO(K2,KK)
               D(M) = Q            ! Possible change of value
            ENDDO
            COB   = COO
         ENDDO

         CLOSE(LUNF)

         END
