
C ++********************************************************************
C
C COPYBRIX
C                   MAXNAM                         JUL 14 ARDEAN LEITH
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
C COPYBRIX(LUN1,LUN2,NX,NY,NZ)
C
C PURPOSE: CONVERTS 3D SPIDER IMAGES FOR INPUT DIRECTLY TO 'O' FORMAT
C          COMPATIBLE WITH 'O' VERSION 5.10.  JF 3/14/95
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C **********************************************************************

        SUBROUTINE COPYBRIX(LUN1,LUN2,NX,NY,NZ)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON /IOBUF/  BUF(NBUFSIZ)

        INTEGER               :: EXTENT(3),ORIGIN(3),GRID(3),PLUS,ERRCOD
        INTEGER               :: XYZBRIX(3)
        REAL                  :: CELL(6), PROD, SIGMA, GINT

        CHARACTER(LEN=MAXNAM) :: FILNAM
        CHARACTER(LEN=1)      :: NULL = CHAR(0)
        CHARACTER(LEN=1)      :: ANS 

	IF (IMAMI .NE. 1) CALL NORM3(LUN1,NX,NY,NZ,FMAX,FMIN,AV)
        NPT = NX*NY*NZ

C       OPEN BRIX FILE
        CALL OPAUXFILE(.TRUE.,FILNAM,NULL,LUN2,512,'N',
     &                       'BRIX OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        CALL RDPRM(GINT,NOT_USED,'SAMPLING DISTANCE')

        CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &     'DEFAULT UNIT CELL IS NX,NY,NZ. OVERWRITE (Y/N)?',
     &     NULL,IRTFLG)

C       PREPARE PARAMETERS FOR HEADER

        ORIGIN(1) = -NX/2
        ORIGIN(2) = -NY/2
	IF (MOD(NZ,2) .EQ. 0) THEN
	   ORIGIN(3) =  -NZ/2+1
	ELSE
           ORIGIN(3) = -NZ/2
	ENDIF

        EXTENT(1) = NX 
        EXTENT(2) = NY 
        EXTENT(3) = NZ

        CELL(1)   = EXTENT(1) * GINT
        CELL(2)   = EXTENT(2) * GINT
        CELL(3)   = EXTENT(3) * GINT
        DO I =4,6
          CELL(I) = 90.0
        ENDDO

        GRID(1)   = NX
        GRID(2)   = NY
        GRID(3)   = NZ

C       ALL DATA WILL BE RESCALED TO 1000, BECAUSE OF THE
C       FAILURE OF O TO READ IN CONTOUR LEVEL SPECIFICATIONS INVOLVING
C       DIGITS BEHIND THE DECIMAL POINT.

        FMAX  = FMAX*1000.
        FMIN  = FMIN*1000.
        PROD  = 255./(FMAX-FMIN)
        PLUS  = -FMIN*PROD
        SIGMA = SIG* 1000.

        IF (ANS.NE.'N' .AND. ANS.NE.'n') THEN
C          CALL RDPRMI(NX,NY,NOT_USED,'UNIT CELL DIMS (X,Y)')
C          CALL RDPRMI(NZ,NDUM,NOT_USED,'UNIT CELL DIMS (Z)')
           CALL RDPRI3S(NX,NY,NZ,NOT_USED,
     &                 'UNIT CELL DIMS (X,Y,Z)',IRTFLG)
           IF (IRTFLG .EQ. -1) THEN
              CLOSE (LUN1)
              RETURN
           ENDIF

          CELL(1) = NX* GINT
          CELL(2) = NY* GINT
          CELL(3) = NZ* GINT
        ENDIF 

        DO I=1,3
           XYZBRIX(I) = EXTENT(I)/8
           IF(MOD(EXTENT(I),8) .GE. 1) XYZBRIX(I) = XYZBRIX(I) +1
        ENDDO

C       SET FORMAT SWITCH TO "NEW"
        CALL PAGED_FORMAT('NEW')

C       NOW WRITE HEADER
        CALL PAGEDHDR(LUN2, ORIGIN, EXTENT, GRID, CELL, 
     &               SIGMA, PROD, PLUS, ERRCOD)
        IF (ERRCOD .NE. 0) THEN
            WRITE (NOUT,*)  '*** ERROR WRITING BRIX HEADER'
            CALL ERRT(100,'COPYBRIX',NE)
        ENDIF

C       NOW READ SPIDER VOLUME AND RECORD IN BRIX FILE
C       READ UPSIDE DOWN TO CORRECT HANDEDNESS FLIP
        INUM = 1
        DO ISEC = 1, NZ
           IRECO = (NZ-ISEC)*NY
	   DO J = 1,NY
              ISTART = (J-1) * NX + 1
              CALL REDLIN(LUN1, BUF(ISTART), NX, IRECO+J)
              DO I = 1,NX
                 BUF(ISTART+I-1) = (BUF(ISTART+I-1)*1000.-FMIN) * PROD
              ENDDO
           ENDDO

           CALL PAGED (LUN2, BUF(1), BUF(NX*NY+1),
     &        NX, NY, XYZBRIX(1), XYZBRIX(2), ISEC, INUM, ERRCOD)

           IF (ERRCOD .NE. 0) THEN
              WRITE (NOUT,*)  '*** ERROR WRITING BRICK'
              CALL ERRT(100,'COPYBRIX',NE)
           ENDIF
        ENDDO
C       END ISEC LOOP OVER SLICES

        ISEC = 0
        CALL PAGED (LUN2, BUF(1), BUF(NX*NY+1),
     &     NX, NY, XYZBRIX(1), XYZBRIX(2), ISEC, INUM, ERRCOD)

999     CLOSE (LUN1)
        CLOSE (LUN2)

        END


C=======================================================================

C WRITTEN BY MORTEN KJELDGAARD, NOV 1993.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE PAGEDHDR (OLUN, ORIGIN, EXTENT, GRID, CELL, 
     $                    SIGMA, PROD, PLUS, ERRCOD)

#ifndef SP_SUN4
      IMPLICIT     NONE
#endif

C     I DO NOT KNOW IF SAVE IS NEEDED FEB 99 al
      SAVE

      INTEGER  OLUN, ORIGIN(3), EXTENT(3), GRID(3), PLUS, ERRCOD
      REAL     CELL(6), PROD, SIGMA
      COMMON /PAGEDC/ FMT
      CHARACTER*4  FMT
      CHARACTER    STR*512

C     FOR OLD FORMAT
      INTEGER*2 FIRST(256)               
      INTEGER   I1, I2, I                  

      IF (FMT .EQ. 'OLD') THEN
C ---   FILL UP FIRST RECORD
         I1 = 80
         I2 = 100
         DO  I=1,3
            FIRST(I  ) = ORIGIN(I)
            FIRST(I+3) = EXTENT(I)
            FIRST(I+6) = GRID(I)
            FIRST(I+9) = I1*CELL(I)
            FIRST(I+12)= I1*CELL(I+3)
	 ENDDO
         FIRST(16) = I2*PROD
         FIRST(17) = PLUS
         FIRST(18) = I1
         FIRST(19) = I2
         DO  I=20,256
            FIRST(I) = 0
	 ENDDO

         WRITE (OLUN, REC=1, IOSTAT=ERRCOD) FIRST

      ELSE
         WRITE (STR, 10) ORIGIN, EXTENT, GRID, CELL, PROD, PLUS, SIGMA

 10      FORMAT (':-) ORIGIN', 3I5,' EXTENT', 3I5, ' GRID', 3I5,
     &        ' CELL ', 6F10.3, ' PROD', F12.5, ' PLUS',I8, 
     &        ' SIGMA ', F12.5)
      
         WRITE (OLUN, REC=1, IOSTAT=ERRCOD) STR

      ENDIF
      RETURN
      END



C=======================================================================
        SUBROUTINE PAGED_FORMAT (STR)
#ifndef SP_SUN4
      IMPLICIT     NONE
#endif
        CHARACTER*(*) STR
        COMMON /PAGEDC/ FMT
        CHARACTER*4 FMT

        IF (STR == 'OLD') THEN
           FMT = 'OLD'
        ELSE
           FMT = 'NEW'
        ENDIF
        END

C=======================================================================
        SUBROUTINE BYTSWP (REC)
#ifndef SP_SUN4
      IMPLICIT     NONE
#endif
        LOGICAL * 1 REC(2,256)
        LOGICAL * 1 ONE
        INTEGER     I

        DO  I=1,256
          ONE      = REC (1,I)
          REC(1,I) = REC(2,I)
          REC(2,I) = ONE
	ENDDO
        RETURN
        END


C A new level of density, store it and if necessary, write it out as 
C 3-d non-overlapping boxes of 8*8*8 values
C Original logic by Alwyn Jones, a long time ago.
C Modified as library routine, new brick format, Morten Kjeldgaard, Nov 1993
c=======================================================================

      SUBROUTINE PAGED (OLUN, RHOSEC, SLICES, IX, IY, NX, NY, 
     &                  ILEV, INUM, ERR)

#ifndef SP_SUN4
      IMPLICIT     NONE
#endif

C     I DO NOT KNOW IF SAVE IS NEEDED FEB 99 al
      SAVE

      INTEGER      OLUN, IX, IY, NX, NY, ERR
      REAL         RHOSEC(IY, IX), SLICES(8,IX*IY)

      COMMON      /PAGEDC/ FMT
      CHARACTER*4 FMT

      INTEGER     I, ICT, I1, ILEV, J, JCT, J1, J2, J3, K, K1, K2, K3 
      INTEGER     VALUE, INUM
      CHARACTER*1 RECORD(512)
      INTEGER*2   IREC(256)
      EQUIVALENCE (RECORD, IREC) 

      ERR = 0

      I1  = MOD (ILEV,8)
      IF (I1 .EQ. 0) I1 = 8
      ICT = 0

      IF (ILEV .NE. 0) THEN
         DO I=1,IY
            DO J=1,IX
               ICT = ICT+1
               SLICES(I1,ICT) = RHOSEC(I,J)
            ENDDO
         ENDDO
      ENDIF

C     PICK OUT OUR NON-OVERLAPPING BRICKS?
      IF (I1 .NE. 8 .AND. ILEV .NE. 0) RETURN

      VALUE = 0
C     LOOP OVER POSSIBLE Y-PAGES
      DO J=1,NY
         J1 = (J-1)*8+1
         J2 =   J  *8

C        LOOP OVER POSSIBLE X-PAGES
         DO K=1,NX
            K1 = (K-1)*8+1
            K2 =   K  *8
            ICT = 0

C           LOOP OVER Z-LEVELS
            DO I=1,8

C              LOOP OVER Y-INDECES OF CURRENT PAGE
               DO J3=J1,J2
                  JCT = (J3-1)*IX+K1-1

C                 LOOP OVER X-INDECES OF CURRENT PAGE
                  DO K3=K1,K2
                     ICT = ICT+1
                     JCT = JCT+1

C                    IF EITHER DIRECTION OVER EDGE, PACK RECORD
                     IF (J3 .GT. IY  .OR.  K3 .GT. IX 
     $                    .OR. I .GT. I1)  THEN
                        RECORD(ICT) = CHAR(0)
                     ELSE
                        RECORD(ICT) = CHAR(INT(SLICES(I,JCT)))
                     END IF
                  ENDDO
               ENDDO
            ENDDO

            INUM = INUM+1

            IF (FMT .EQ. 'OLD') THEN
               CALL BYTSWP (IREC)
               WRITE (OLUN, REC=INUM, IOSTAT=ERR) IREC
            ELSE
               WRITE (OLUN, REC=INUM, IOSTAT=ERR) RECORD
            ENDIF

         ENDDO
      ENDDO

      RETURN
      END



