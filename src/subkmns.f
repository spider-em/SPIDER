C++*********************************************************************
C
C SUBKMNS.F                      USED REG_GET AUG 00 ARDEAN LEITH
C                                USED ALLOCATE JAN 01 ARDEAN LEITH
C                                _PIX NOT _EIG MAY 04 ARDEAN LEITH
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
C PURPOSE:
C       READS THE NAME OF THE FILE (IMC***, PIX***, EIG***)
C       TO GET NUMBER (AND THEIR COORDINATES) OF OBJECTS TO BE PARTITIONED
C
C--*********************************************************************

         SUBROUTINE SUBKMNS(LUNF,LUND)

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC'

         REAL, ALLOCATABLE, DIMENSION(:) :: Q
         CHARACTER(LEN=MAXNAM)           :: FILNAM
         CHARACTER(LEN=1)                :: NULL
         LOGICAL                         :: USE_DISK
 
        NULL = CHAR(0)

        WRITE(NOUT,*) ' YOU MAY USE A _SEQ, _PIX, or _IMC FILE'
        WRITE(NOUT,*) ' '

        CALL FILERD(FILNAM,NLET,NULL,
     &              'CORAN/PCA FILE (e.g. CORAN_IMC)~',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (INDEX(FILNAM,'_SEQ') .GT. 0) THEN
            CALL OPAUXFILE(.FALSE.,FILNAM,DATEXC,-LUNF,0,
     &                 'O',' ',.TRUE.,IRTFLG)
        ELSE
            CALL OPAUXFILE(.FALSE.,FILNAM,DATEXC,LUNF,0,
     &                 'O',' ',.TRUE.,IRTFLG)
        ENDIF
        IF (IRTFLG .NE. 0) RETURN

        IF (INDEX(FILNAM,'_SEQ') .GT. 0) THEN
            ITYPE = 1
            WRITE(NOUT,*)' In sequential coordinates file ---'
            READ(LUNF)  NOBJ, NFTOT

        ELSEIF(INDEX(FILNAM,'_IMC') .GT. 0) THEN
            ITYPE = 2
            WRITE(NOUT,*)' In image coordinates file ---'
            READ(LUNF,*)  NOBJ, NFTOT

        ELSE
            ITYPE = 3
            WRITE(NOUT,*)' In pixel coordinates file ---'
            READ(LUNF,*)  NOBJ, NFTOT
        ENDIF

        WRITE(NOUT,*)' Number of objects:', NOBJ
        WRITE(NOUT,*)' Number of factors:', NFTOT
        WRITE(NOUT,*)'  '

        CALL RDPRI1S(KCLASS,NOT_USED,'NUMBER OF CLASSES',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       CALCULATE NUMBER OF FACTORS : NFAC (RETURNED BY RDPRAI)
        NFAC = NIMAX
        CALL RDPRAI(INUMBR, NIMAX, NFAC, 0,9999999,
     &               'FACTOR NUMBERS', NULL, IRTFLG)

C	SIZE OF W() = Q(K_W) ARRAY IS: NFAC
	K_W    = 1
        K_CM   = K_W    + NFAC
        K_IC   = K_CM   + KCLASS * NFAC
        K_NC   = K_IC   + NOBJ
	K_NIMI = K_NC   + KCLASS
        K_WSS  = K_NIMI + NOBJ 
        K_AT   = K_WSS  + KCLASS
        K_COO  = K_AT   + NFAC
        K_X    = K_COO  + NFTOT
        MEMTOT = K_X    + NFAC * NOBJ

        ALLOCATE(Q(MEMTOT),STAT=IRTFLG)

        IF (IRTFLG .EQ. 0) THEN
C         GOT ENOUGH SPACE
          USE_DISK = .FALSE.
          IX_DIM   = NOBJ
        ELSE
	  MEMTOT = K_X + NFAC
          ALLOCATE(Q(MEMTOT),STAT=IRTFLG)
          IF (IRTFLG .NE. 0)  THEN
             CALL  ERRT(46,'Q',MEMTOT)
             RETURN
	  ENDIF

C	  NOT EVERYTHING FITS IN MEMORY. USE DISK SPACE.
C	  COPY ONLY ONE LINE IN X(1)
	  USE_DISK = .TRUE.
          IX_DIM    = 1
          WRITE(NOUT, *)'  WARNING - SLOW ON-DISK VERSION USED'
        ENDIF

        DO  I=0,NFAC-1
          Q(K_W + I) = 1.0
        ENDDO       

        DO I=0,NFAC-1
           W1 = 1.0
           CALL  RDPRM1S(W1,NOT_USED,'FACTOR WEIGHT',IRTFLG)
           IF (W1 == 0.0)  EXIT   ! WANT ALL SAME = 1.0
           Q(K_W + I) = W1
        ENDDO

        IF (W1 == 0.0) THEN
           WRITE(NOUT, 90)
90         FORMAT('  ALL FACTOR WEIGHTS: 1.0')

        ELSE
           WRITE(NOUT, 91) (Q(K_W + I), I = 0, NFAC-1)
91         FORMAT('  FACTOR WEIGHTS USED:', 10(F10.4, 1X))
        ENDIF

        ISEED = 1457
        CALL  RDPRI1S(ISEED,NOT_USED,
     &       'FOR RANDOM SEEDS GIVE NON-ZERO STARTING NUMBER',IRTFLG)
	ISEED = MOD(IABS(ISEED),340189)

        CALL SUBKMEANS(Q(K_X), Q(K_CM), Q(K_NIMI),
     &                  Q(K_IC), Q(K_NC), Q(K_WSS), Q(K_AT),
     &                  NOBJ, NFAC, KCLASS, Q(K_W), NFTOT, Q(K_COO),
     &                  INUMBR, BB, WW, DB, IX_DIM, USE_DISK,ISEED,
     &                  LUNF,LUND,ITYPE)

        C = BB * WW
        H = (BB / (KCLASS - 1)) / (WW / (NOBJ - KCLASS))

        WRITE (NDAT, 41)  NOBJ, NFAC, KCLASS
41      FORMAT(/,
     &         /,'   K - M E A N S   C L U S T E R   A N A L Y S I S',/,
     &         /,'   NUMBER OF OBJECTS: ',I8
     &         /,'   NUMBER OF FACTORS: ',I8
     &         /,'   NUMBER OF GROUPS : ',I8)

        CALL PRNTXX(Q(K_NC), Q(K_WSS), KCLASS, NDAT)

        WRITE(NDAT,42)  BB, WW, C, H, DB
42      FORMAT(/,
     &     '   Tr(B)                                              :',
     &            G12.5,/,
     &     '   Tr(W)                                              :',
     &            G12.5,/,
     &     '   Coleman criterion   B*W                            :',
     &            G12.5,/,
     &     '   Harabasz criterion (B/(kclass-1))/(W/(nobj-kclass)):',
     &            G12.5,/,
     &     '   Davies-Bouldin criterion                           :',
     &            G12.5,/)

        WRITE(NDAT, *)
        CALL REG_SET_NSEL(1,5,BB,WW,C,H,DB,IRTFLG)

9999    IF (ALLOCATED(Q)) DEALLOCATE(Q)

        END

C--*************************** PRNTXX ******************************

         SUBROUTINE PRNTXX(NC,WSS,K,NDAT)

         DIMENSION   NC(K),WSS(K)

         DATA  L/8/

         K1 = 1
         K2 = MIN(K,L)

         WRITE(NDAT,5)
5        FORMAT( /,'   NUMBER OF OBJECTS IN EACH CLUSTER AND',
     &             ' WITHIN SUM OF SQUARES.')

101      WRITE(NDAT,1)  (I,I=K1,K2)
1        FORMAT(/8(I7,5X))

         WRITE(NDAT,2)  (NC(I),I=K1,K2)
2        FORMAT(8(I7,5X))

         WRITE(NDAT,3)  (WSS(I),I=K1,K2)
3        FORMAT(8(2X,G10.3))

         K1 = K1 + L
         IF (K1 .GT. K)  RETURN

         K2 = MIN(K2+L,K)
         GOTO  101

         END
