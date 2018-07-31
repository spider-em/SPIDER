
C++*********************************************************************
C
C MRCP.F                    FILENAMES LENGTHENED   JAN 89 ARDEAN LEITH
C                           USED OPFILE            NOV 00 ARDEAN LEITH
C                           YEARS OLD BUG FIXED    MAR 16 ARDEAN LEITH
C                           INTEGER BUG            JUN 18 ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2018  Health Research Inc.,                         *
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
C  MRCP(LUN1,LUN2,LUN3)
C
C  VARIABLES:      NX          COLUMNS
C                  NY          ROWS
C                  NZ          SLICES
C                  LUN1        INPUT FILE
C                  LUN2        OUTPUT FILE
C                  LUN3        FOR DOC FILE
C                  A           ARRAY FOR SLICE
C                  BUF         IO ARRAY
C
C  NOTE:  CYLINDER X-SECTION IS IN XZ PLANE
C
C--*********************************************************************

        SUBROUTINE MRCP(LUN1,LUN2,LUN3)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

        INTEGER               :: LUN1,LUN2,LUN3

        INTEGER, PARAMETER    :: MAXREG = 7 
        INTEGER, PARAMETER    :: MAXKEY = 4000 
        REAL                  :: DBUF(MAXREG,MAXKEY,2)

        CHARACTER(LEN=MAXNAM) :: FLN1,DOCF1,FILNAM

        CHARACTER *1          :: NULL = CHAR(0)
        REAL                  :: PLIST(6)
    
        LOGICAL               :: MITO,DIST

        INTEGER               :: NX,NY,NZ,MAXIM,IRTFLG,MWANT,NLET
        INTEGER               :: PHI0
        REAL                  :: RADI,RAD,POINTS,AINC,WINK
        INTEGER               :: IPOINTS,NOT_USED,IRAD
        INTEGER               :: IRADI,NXP,IRX,IRY,NOPEN,NKEY,NREG,I

        INTEGER               :: NLETD,LERR,LAUF,N0,K,IX,IY
        REAL                  :: SUM,PHI,CP,SP,X,Y,XCENTER,YCENTER
        REAL                  :: DX,DY,DDX,DDY,B1,B2,W

        REAL, ALLOCATABLE     :: A(:,:), BUF(:)

        REAL, PARAMETER       :: PI = 3.14159


C       CYLINDRICAL PROJECTION
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,
     &                  NX,NY,NZ,
     &                  MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        MWANT = NX * NZ + NBUFSIZ
        ALLOCATE (BUF(NBUFSIZ), A(NX,NZ), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'MRCP, BUF...',MWANT)
           RETURN
        ENDIF

        CALL FILERD(FLN1,NLET,NULL,'OUTPUT',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 998

        PHI0 = 0
        CALL RDPRM1S(PHI0,NOT_USED,
     &              'STARTING ANGLE (0 = 3 OCLOCK)',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 998

        PHI0 = PI * (PHI0 / 180.0)

        RADI = 5
        RAD  = (NX / 2) -5

        CALL RDPRM2S(RADI,RAD,NOT_USED,'INNER, OUTER RADIUS',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 998

        MITO = .FALSE.
        DIST = .FALSE.
        IF (RADI  <   0.0) THEN
           RADI = - RADI
           MITO = .TRUE.
           WRITE(NOUT,*) 
     &         ' NEGATIVE INNER RADIUS  --> MODIFED MAX. PROJ.'
        ENDIF
        IF (RAD  <   0.0) THEN
           RAD  = - RAD
           DIST = .TRUE.
           WRITE(NOUT,*) 
     &         ' NEGATIVE OUTER RADIUS --> DISTANCE FROM CENTER.'
        ENDIF

        POINTS  = 2 * PI * RAD
        IPOINTS = POINTS
        AINC    = (360.0 / POINTS)                      
        WRITE(NOUT,100) AINC,IPOINTS
100     FORMAT('  ANGULAR INCREMENT: ',F7.2,
     &         ' DEGREES,      X DIMENSION:',I5)
        
88      CALL RDPRM(WINK,NOT_USED,
     &       'NEW ANGULAR INCREMENT OR <CR>')

        IF (WINK .NE. 0) THEN
          AINC    = WINK
          POINTS  = 360.0 / AINC
          IPOINTS = POINTS
          WRITE(NOUT,100) AINC,IPOINTS
          GOTO 88
        ENDIF           

        IRAD   = RAD
        IRADI  = RADI
        AINC   = PI * (AINC / 180.0)
        NXP    = POINTS                 

        IF (NXP > NBUFSIZ) THEN
C          TOO MANY POINTS FOR BUF ARRAY SIZE
           WRITE(NOUT,*) ' *** ONLY: ',NBUFSIZ,' POINTS ALLOWED'
           GOTO 88
        ENDIF

        IFORM = 1
        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FLN1,LUN2,'U',IFORM,NXP,NY,1,
     &                 MAXIM,' ',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

10      CALL RDPRMI(IRX,IRY,NOT_USED,
     &              'COLUMNS FOR XCENTER & ZCENTER')   
        NOPEN = 0       
        NKEY  = 1
        NREG  = 1
        I     = 1
        CALL FILERD(DOCF1,NLETD,DATEXC,'DOCUMENT',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL UNSDAL(DOCF1,NOPEN,LUN3,I,PLIST,NREG,
     &     DBUF,MAXKEY,MAXREG,NKEY,LERR)
        IF (LERR .NE. 0) GOTO 998

C     GO THROUGH THE VOLUME
      DO  LAUF=1,NY
        XCENTER = DBUF(IRX+1,LAUF,1)
        YCENTER = DBUF(IRY+1,LAUF,1)
        WRITE(NOUT,102) LAUF,XCENTER,YCENTER
102     FORMAT('  SLICE #:',I4,'  CENTER AT: (',F7.2,',',F7.2,')')

C       READ IN SLICE (PERPENDICULAR TO Y)
        DO I = 1,NZ
          N0 = (I-1) * NY + LAUF
          CALL REDLIN(LUN1,BUF,NX,N0)
          DO  K=1,NX
            A(K,I) = BUF(K)
          ENDDO
        ENDDO

C       SLICE IS THERE, NOW  DO THE REAL WORK          
C       LOOP ALONG PHI
        DO  I=1,IPOINTS
          SUM = 0.0
          IF (MITO .OR. DIST) SUM = FMIN

          PHI = PHI0 + I * AINC
          CP  = COS(PHI)
          SP  = SIN(PHI)

          DO  K=IRADI,IRAD
            RAD = K
            X   = RAD * CP + XCENTER
            Y   = RAD * SP + YCENTER

            IF (X >= 1 .AND. X <  NX .AND.
     &          Y >= 1 .AND. Y <  NZ) THEN
C             PIXEL IS WITHIN SLICE
              IX  = INT(X)
              IY  = INT(Y)


              DX  = X - IX
              DY  = Y - IY
              DDX = 1 - DX
              DDY = 1 - DY
              IF (MITO) THEN
C                MODIFED MAXIMUM PROJECTION 
                 SUM = MAX(A(IX,IY), SUM)
 
              ELSEIF (DIST) THEN
C                DISTANCE FROM CENTER OF PROJECTION 
                 IF (A(IX,IY) > 0.0) THEN
C                   HAVE A POSITIVE PIXEL, FIND DISTANCE FROM CENTER
                    SUM = RAD 
                 ENDIF
 
              ELSE
                 B1  = A(IX,IY)   * DDX + A(IX+1,IY)   * DX
                 B2  = A(IX,IY+1) * DDX + A(IX+1,IY+1) * DX
                 W   = B2 * DY + B1 * DDY
                 SUM = SUM + W
              ENDIF
            ENDIF
	  ENDDO

          BUF(I) = SUM
	ENDDO
     
        CALL WRTLIN(LUN2,BUF,NXP,LAUF)
      ENDDO

998   CLOSE(LUN2)
999   CLOSE(LUN1)
      IF (ALLOCATED(BUF))  DEALLOCATE(BUF)
      IF (ALLOCATED(A))    DEALLOCATE(A)

      END

