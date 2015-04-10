 
C++********************************************************************
C
C  MASK.F        ADDED INNER BI                   SEP 99 ARDEAN LEITH
C                NXR = -9999999 RDRPIS FLOAT BUG  MAR 12 ARDEAN LEITH          
C                NZR BUG                          MAY 12 ARDEAN LEITH          
C                PROMTS FIXED                     FEB 14 ARDEAN LEITH          
C                ADDED MIN BACKGROUND             DEC 14 ARDEAN LEITH
C                    
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C   MASK(LUNI,LUNO,NX,NY,B,FMIN1)
C
C   PARAMETERS:
C        LUNI               LOGICAL UNIT NUMBER OF INPUT FILE    (SENT)
C        LUNO               LOGICAL UNIT NUMBER OF OUTPUT FILE   (SENT)
C        NX,NY,NZ           DIMENSIONS OF FILE                   (SENT)
C        B                  AVERAGE OF INPUT FILE      (SENT / ALTREED)
C        FMIN1              MIN OF INPUT FILE          (SENT)
C
C--********************************************************************

      SUBROUTINE MASK(LUNI,LUNO,NX,NY,NZ,B,FMIN1)

      INCLUDE 'CMBLOCK.INC' 

      INTEGER              :: LUNI,LUNO,NX,NY,NZ 
      REAL                 :: B,FMIN1

      REAL                 :: BUF(NX)

      INTEGER              :: ILIST(3) 
      CHARACTER (LEN=1)    :: MODE,NULL,ANS
      DOUBLE PRECISION     :: DAV,AVC,AVCI

      PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)

      IF (FCHAR(4:4) == 'L') THEN
C        MASK A LINE ACROSS IMAGE
         IF (NZ > 1) THEN
C           3-D FILE
            CALL RDPRI1S(NZR,NOT_USED,'SLICE NUMBER',IRTFLG)
            IF (IRTFLG .NE. 0) RETURN
         ELSE
            NZR = 1
         ENDIF

         CALL RDPRI1S(NYR,NOT_USED,
     &        'LINE NUMBER TO BE MASKED',IRTFLG)
            IF (IRTFLG .NE. 0) RETURN

         CALL RDPRM1S(B,NOT_USED,'BACKGROUND',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         IRECT = (NZR - 1) * NY + NYR

         DO I = 1,NY*NZ
            CALL REDLIN(LUNI,BUF,NX,I)

            IF (I == IRECT) THEN
C              MASK THIS LINE 
               BUF = B
            ENDIF
            CALL WRTLIN(LUNO,BUF,NX,I)
         ENDDO
         RETURN          

      ELSE
         RAD  = 0.0
         RADI = 0.0

         IF (FCHAR(4:4) == 'X' .OR.
     &       FCHAR(4:4) == 'Y' .OR.
     &       FCHAR(4:4) == 'Z') THEN
           CALL RDPRM2S(RAD,RADI,NOT_USED,
     &                  'OUTER & INNER HALFWIDTH EXTENTS',IRTFLG)
         ELSE
            CALL RDPRM2S(RAD,RADI,NOT_USED,
     &                   'OUTER & INNER RADII',IRTFLG) 
         ENDIF
         IF (IRTFLG .NE. 0) RETURN

         IF (RAD == 0.0) THEN
C           FOR INNER ONLY MASKING, SET OUTER RADIUS TO HUGE VALUE
            RAD = HUGE(RAD)
         ENDIF

         IF (RAD < 0.0 .OR. RADI < 0.0 .OR. RAD < RADI)  THEN
             CALL ERRT(101,'INCONSISTENT INPUT PARAMETERS',IER)
             RETURN
          ENDIF
       ENDIF

       IF (FCHAR(4:4) == 'X' .OR. 
     &     FCHAR(4:4) == 'Y' .OR.
     &     FCHAR(4:4) == 'Z') THEN

          CALL RDPRMC(MODE,NCHAR,.TRUE.,
     &     'SHARP, COSINE, GAUSSIAN EDGE, OR TRUE GAUSSIAN (S/C/G/T)',
     &      NULL,IRTFLG)
          IF (MODE == 'S') MODE = 'D'
       ELSE
          CALL RDPRMC(MODE,NCHAR,.TRUE.,
     &     'DISK, COSINE, GAUSSIAN EDGE, OR TRUE GAUSSIAN (D/C/G/T)',
     &     NULL,IRTFLG)
       ENDIF
       IF (IRTFLG .NE. 0) RETURN

C       BACKGROUND CHOICES:
C       DAV     BACKGROUND IS SET EQUAL TO THE AVERAGE OF THE
C               IMAGE BEFORE MASKING
C       PREV AV BACKGROUND IS SET EQUAL TO THE AVERAGE OF THE
C               IMAGE AREA PASSED BY THE MASK
C       CIRCUMF BACKGROUND IS SET EQUAL TO THE AVERAGE OF THE 
C               PIXELS ALONG THE MASK'S CIRCUMFERENCE
C       MIN     BACKGROUND IS SET EQUAL TO PREVIOUS FILE MINIMUM
C       EXTERNAL  BACKGROUND IS SET TO A VALUE SUPPLIED EXTERNALLY

        CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &'AVG, PREV AVG, CIRCUMF, MIN, OR EXTERNAL BACKGROUND (A/P/C/M/E)',
     &      NULL,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (ANS == 'E') THEN
           B = 0.0
           CALL RDPRM1S(B,NOT_USED,'BACKGROUND',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

        ELSEIF (ANS == 'M') THEN
           B = FMIN1
           WRITE(NOUT,'(A,1PG12.2)') '  Mask background: ',B
        ENDIF

        NXR = (NX/2) + 1
        NYR = (NY/2) + 1
        NZR = (NZ/2) + 1

        IF (FCHAR(4:4) == 'X') THEN
        
           NYR = 1
           NZR = 1
           CALL RDPRI1S(NXR,NOT_USED,
     &        'X CENTER OF MASK OR <CR> FOR IMAGE CENTER',
     &        IRTFLG)

         ELSEIF ( FCHAR(4:4) == 'Y') THEN
            
           NXR = 1
           NZR = 1
           NVAL = 2
           CALL RDPRAI(ILIST,2,NVAL,ILOW,IHI,
     &         'Y CENTER OF MASK OR <CR> FOR IMAGE CENTER',
     &         .FALSE.,IRTFLG)

           !write(6,*) nval, ilist

           IF     (NVAL == 0) THEN
               NYR = (NY/2) + 1
           ELSEIF (NVAL == 1) THEN
               NYR = ILIST(1)
           ELSE 
               NYR = ILIST(2)
           ENDIF

         ELSEIF ( FCHAR(4:4) == 'Z') THEN
            
           NXR = 1
           NYR = 1
           NVAL = 3
           CALL RDPRAI(ILIST,3,NVAL,ILOW,IHI,
     &         'Z CENTER OF MASK OR <CR> FOR VOLUME CENTER',
     &         .FALSE.,IRTFLG)

           IF     (NVAL == 0) THEN
               NZR = (NZ/2) + 1
           ELSEIF (NVAL == 1) THEN
               NZR = ILIST(1)
           ELSEIF (NVAL == 2) THEN
               NZR = (NZ/2) + 1
           ELSE 
               NZR = ILIST(3)
           ENDIF
        
         ELSEIF (NZ == 1) THEN

            CALL RDPRI2S(NXR,NYR,NOT_USED,
     &         'MASK CENTER LOCATION (X,Y) OR <CR> FOR IMAGE CENTER',
     &         IRTFLG)

         ELSE

            CALL RDPRI3S(NXR,NYR,NZR,NOT_USED,
     &        'MASK CENTER LOCATION (X,Y,Z) OR <CR> FOR VOLUME CENTER',
     &        IRTFLG)
         ENDIF

        IF (IRTFLG .NE. 0) RETURN


        !IF (NXR == HUGE(NXR)) THEN  ! FAILS AS HUGE(NXR) NOW < 0!!
        IF (NXR == -9999999) THEN
           NXR = (NX/2) + 1
           NYR = (NY/2) + 1
           IF (NZ > 1) NZR = (NZ/2) + 1
        ENDIF

C       ALTHOUGH MASK CENTER COORDINATES MAY BE OUTSIDE THE IMAGE
C       ITSELF, CHECK TO BE SURE HE RADIUS IS LARGE ENOUGH SO
C       THAT THE MASK WILL AT LEAST HAVE AN EFFECT ON THE IMAGE.
C       (APPLIES ONLY TO EXCLUSIVE OUTSIDE MASKING)

        IF (RAD > 0) THEN
C          POSITIVE OUTER RADIUS GIVEN
           IF (((NXR + RAD)      <  0.0) .OR.
     &         ((NXR - RAD - NX) >= 0.0) .OR.
     &         ((NYR + RAD)      <  0.0) .OR.
     &         ((NYR - RAD - NY) >= 0.0) .OR.
     &         ((NZR + RAD)      <  0.0) .OR.
     &         ((NZR - RAD-NZ)   >= 0.0)) THEN
              WRITE(NOUT,61)
61            FORMAT(' *** MASK OUTSIDE IMAGE; NO ACTION TAKEN')
              RETURN
           ENDIF

        ELSEIF (RAD <= 0.0) THEN
C          INNER ONLY MASKING, SET OUTER RADIUS TO HUGE VALUE
           RAD = HUGE(RAD)
        ENDIF

        RAD2   = RAD     **2
        RADI2  = RADI    **2
        RAD2P  = (RAD+1) **2
        RADI2P = (RADI-1)**2

C       RAD2P IS USED TO DEFINE A NON-EMPTY SET OF POINTS
C       ALONG THE MASK'S CIRCUMFERENCE

        IF (ANS == 'P' .OR. ANS == 'C') THEN
           DAV   = 0.0
           NAV   = 0.0
           AVC   = 0.0
           NAVC  = 0.0
           AVCI  = 0.0
           NAVCI = 0.0

           DO  J = 1,NZ
             FI2 = FLOAT(J-NZR)**2

             DO  I = 1,NY
               FI1 = FLOAT(I-NYR)**2 + FI2
               CALL REDLIN(LUNI,BUF,NX,I+(J-1)*NY)

               DO  K = 1,NX
                 CRAD2 = FI1 + FLOAT(K-NXR)**2
                 IF (CRAD2 <= RAD2 .AND. CRAD2 >= RADI2) THEN
                    DAV = DAV + BUF(K)
                    NAV = NAV + 1
                 ENDIF

C                PROVISION FOR CIRCUMFERENCE OPTION
                 IF (ANS == 'C') THEN
                    IF (CRAD2 >= RAD2 .AND. CRAD2 <= RAD2P) THEN
C                      FIND OUTER CIRCUMFERENCE
                       AVC  = AVC + BUF(K)
                       NAVC = NAVC + 1
                    ENDIF

                    IF (RADI2 > 0.0 .AND.
     &                 CRAD2 >= RADI2P .AND. CRAD2 <= RADI2) THEN
C                      FIND INNER CIRCUMFERENCE
                       AVCI  = AVCI  + BUF(K)
                       NAVCI = NAVCI + 1
                    ENDIF
                 ENDIF
               ENDDO
            ENDDO
          ENDDO

           DAV = DAV / NAV
           B   = DAV
           IF (ANS == 'C') THEN
              IF (NAVC >  0) THEN    
                 B = AVC / NAVC
                 WRITE(NOUT,21) B
21               FORMAT('  AVERAGE ALONG OUTER CIRCUMFERENCE:',1PG12.4)
              ENDIF

              IF (NAVCI > 0) THEN
                 BI = AVCI / NAVCI
                 WRITE(NOUT,22) BI
22               FORMAT('  AVERAGE ALONG INNER CIRCUMFERENCE:',1PG12.4)
              ENDIF
           ENDIF
        ENDIF

C       COSINE, GAUSSIAN OR STRAIGHT CUTOFF (DISK) ...
C       MASKS IN ONE DIRECTION
        IF (FCHAR(4:4)     == 'X') THEN
           SWITCHZ = 0.0
           SWITCHY = 0.0
           SWITCHX = 1.0

        ELSEIF (FCHAR(4:4) == 'Y') THEN
           SWITCHZ = 0.0
           SWITCHY = 1.0
           SWITCHX = 0.0

        ELSEIF (FCHAR(4:4) == 'Z') THEN
           SWITCHZ = 1.0
           SWITCHY = 0.0
           SWITCHX = 0.0

        ELSE
           SWITCHZ = 1.0
           SWITCHY = 1.0
           SWITCHX = 1.0
        ENDIF

        IF (MODE == 'C' ) THEN
C          COSINE EDGE MASKING

           CALL RDPRM1S(HW,NOT_USED,'FALLOFF WIDTH',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           DO  J = 1,NZ
            FI2 = FLOAT(J-NZR)**2*SWITCHZ
            DO  I = 1,NY
              FI1 = FLOAT(I-NYR)**2 *SWITCHY+FI2
              CALL REDLIN(LUNI,BUF,NX,I+(J-1)*NY)

              DO  K = 1,NX
                 CRAD2 = FI1 + SWITCHX * FLOAT(K-NXR)**2
                 IF (CRAD2 > RAD2) THEN
                    SRAD2  = SQRT(CRAD2)
                    WGH    =
     &              (1.0+COS(QUADPI*AMIN1(1.0,ABS(SRAD2-RAD)/HW)))*0.5
                    BUF(K) = B + (BUF(K)-B) * WGH

                 ELSEIF (CRAD2 < RADI2) THEN
                    SRAD2 = SQRT(CRAD2)
                    WGH   =
     &              (1.0+COS(QUADPI*AMIN1(1.0,ABS(SRAD2-RADI)/HW)))*0.5
                 BUF(K) = B + (BUF(K)-B) * WGH
                 ENDIF
              ENDDO
             CALL WRTLIN(LUNO,BUF,NX,I+(J-1)*NY)
            ENDDO
           ENDDO

        ELSE IF (MODE == 'G' ) THEN
C          GAUSSIAN EDGE MASKING

           CALL RDPRM1S(HW,NOT_USED,'FALLOFF HALFWIDTH',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           HW = -1./(HW**2)
           DO  J = 1,NZ
            FI2 = FLOAT(J-NZR)**2 * SWITCHZ

            DO  I = 1,NY
              FI1 = FLOAT(I-NYR)**2 * SWITCHY + FI2
              CALL REDLIN(LUNI,BUF,NX,I+(J-1)*NY)

              DO  K = 1,NX
                 CRAD2 = FI1 + SWITCHX*FLOAT(K-NXR)**2
                 IF (CRAD2 > RAD2) THEN
                    SRAD2  = SQRT(CRAD2)
                    WGH    = HW*(SRAD2-RAD)**2

                    IF (WGH < -50.0)  THEN
                       WGH = 0.0
                    ELSE
                       WGH = EXP(WGH)
                    ENDIF
                    BUF(K) = B + (BUF(K)-B) * WGH

                 ELSEIF (CRAD2 < RADI2) THEN
                    SRAD2  = SQRT(CRAD2)
                    WGH = HW * (SRAD2-RADI)**2
                    IF (WGH < -50.0)  THEN
                       WGH = 0.0
                    ELSE
                       WGH = EXP(WGH)
                    ENDIF
                    BUF(K) = B+(BUF(K)-B) * WGH
                 ENDIF
              ENDDO
             CALL WRTLIN(LUNO,BUF,NX,I+(J-1)*NY)
            ENDDO
          ENDDO

        ELSE IF (MODE == 'T' ) THEN
C          TRUE GAUSSIAN MASKING
           IF (RADI .NE. 0.0)  THEN
              CALL  ERRT(101,
     &           ' NO INNER MASKING FOR TRUE GAUSSIAN MASK',IER)
              RETURN
           ENDIF

           HW = -1.0 / (RAD**2)
           DO  J = 1,NZ
             FI2 = SWITCHZ*FLOAT(J-NZR)**2

             DO  I = 1,NY
               FI1 = SWITCHY*FLOAT(I-NYR)**2+FI2
               CALL REDLIN(LUNI,BUF,NX,I+(J-1)*NY)

               DO K=1,NX
                 CRAD2  = (FI1+SWITCHX*FLOAT(K-NXR)**2)*HW
                 IF (CRAD2 < -50.0)  THEN
                    BUF(K) = B
                 ELSE
                    BUF(K) = B + (BUF(K)-B) * EXP(CRAD2)
                 ENDIF
               ENDDO

               CALL WRTLIN(LUNO,BUF,NX,I+(J-1)*NY)
             ENDDO
           ENDDO

        ELSE
C          DISK MASKING
           DO  J = 1,NZ
             FI2 = SWITCHZ * FLOAT(J-NZR)**2
             DO  I = 1,NY
               FI1 = SWITCHY * FLOAT(I-NYR)**2 + FI2
               CALL REDLIN(LUNI,BUF,NX,I+(J-1)*NY)

               DO  K = 1,NX
                  CRAD2 = FI1 + SWITCHX * FLOAT(K-NXR)**2
                  IF (CRAD2 > RAD2 .OR. CRAD2 < RADI2) BUF(K) = B
                  IF (ANS == 'C'  .AND. CRAD2 < RADI2) BUF(K) = BI
               ENDDO

               CALL WRTLIN(LUNO,BUF,NX,I+(J-1)*NY)
             ENDDO
           ENDDO
        ENDIF

        END

