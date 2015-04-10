
C++*********************************************************************
C
C MRCP.F                    FILENAMES LENGTHENED JAN 89 ARDEAN LEITH
C                           USED OPFILE NOV 00 ARDEAN LEITH
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
C  MRCP(NSAM,NROW,NSLICE,LUN1,LUN2,LUN3,A,BUF,MAXSAM)
C
C  PARAMETERS:     NSAM        COLUMNS
C                  NROW        ROWS
C                  NSLICE      SLICES
C                  LUN1        INPUT FILE
C                  LUN2        OUTPUT FILE
C                  LUN3        FOR DOC FILE
C                  A           ARRAY FOR SLICE
C                  BUF         IO ARRAY
C                  MAXSAM      MAX LENGTH OF BUF ARRAY (SINCE A FOLLOWS IT)
C
C  NOTE:  CYLINDER X-SECTION IS IN XZ PLANE
C
C--*********************************************************************

        SUBROUTINE MRCP(NSAM,NROW,NSLICE,LUN1,LUN2,LUN3,BUF,A,MAXSAM)

        PARAMETER (MAXREG=7)
        PARAMETER (MAXKEY=1000)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
       
C       WARNING BUF AND A ARE IN UNLABELED COMMON IN VTIL2

        CHARACTER(LEN=MAXNAM)   ::   FLN1,DOCF1
        COMMON /COMMUN/  FLN1,DOCF1
        COMMON /DOC_BUF/ DBUF(MAXREG,MAXKEY,2)

        CHARACTER *1  NULL
        DIMENSION     PLIST(6),A(NSAM,NSLICE),BUF(MAXSAM)
        LOGICAL       MITO,DIST

        DATA PI/3.14159/

        NULL = CHAR(0)

        CALL FILERD(FLN1,NLET,NULL,'OUTPUT',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        CALL RDPRM(PHI0,NOT_USED,'STARTING ANGLE (0 = 3 OCLOCK)')
        PHI0    = PI * (PHI0 / 180.0)

        CALL RDPRM2(RADI,RAD,NOT_USED,'INNER, OUTER RADIUS')

        MITO = .FALSE.
        DIST = .FALSE.
        IF (RADI .LT. 0.0) THEN
           RADI = - RADI
           MITO = .TRUE.
           WRITE(NOUT,*) 
     &         ' NEGATIVE INNER RADIUS  --> MODIFED MAX. PROJ.'
        ENDIF
        IF (RAD .LT. 0.0) THEN
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

88      CALL RDPRM(WINK,NOT_USED,'NEW ANGULAR INCREMENT OR <RET>')

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
        NSAMP  = POINTS                 

        IF (NSAMP .GT. MAXSAM) THEN
C          TOO MANY POINTS FOR BUF ARRAY SIZE
           WRITE(NOUT,*) ' *** ONLY: ',MAXSAM,' POINTS ALLOWED'
           GOTO 88
        ENDIF

        IFORM = 1
        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FLN1,LUN2,'U',IFORM,NSAMP,NROW,1,
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
      DO  LAUF=1,NROW
        XCENTER = DBUF(IRX+1,LAUF,1)
        YCENTER = DBUF(IRY+1,LAUF,1)
        WRITE(NOUT,102) LAUF,XCENTER,YCENTER
102     FORMAT('  SLICE # ',I3,'  CENTER AT: (',F7.2,',',F7.2,')')

C       READ IN SLICE (PERPENDICULAR TO Y)
        DO I = 1,NSLICE
          N0 = (I-1)*NROW+LAUF
          CALL REDLIN(LUN1,BUF,NSAM,N0)
          DO  K=1,NSAM
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
            X   = RAD*CP+XCENTER
            Y   = RAD*SP+YCENTER

            IF (X.GE.1 .AND. X.LT.NSAM .AND.
     &          Y.GE.1 .AND. Y.LT.NSLICE) THEN
C             PIXEL IS WITHIN SLICE
              IX  = INT(X)
              IY  = INT(Y)


              DX  = X-IX
              DY  = Y-IY
              DDX = 1-DX
              DDY = 1-DY
              IF (MITO) THEN
C                MODIFED MAXIMUM PROJECTION 
                 SUM = MAX(A(IX,IY), SUM)
 
              ELSEIF (DIST) THEN
C                DISTANCE FROM CENTER OF PROJECTION 
                 IF (A(IX,IY) .GT. 0.0) THEN
C                   HAVE A POSITIVE PIXEL, FIND DISTANCE FROM CENTER
                    SUM = RAD 
                 ENDIF
 
              ELSE
                 B1  = A(IX,IY)  *DDX+A(IX+1,IY)  *DX
                 B2  = A(IX,IY+1)*DDX+A(IX+1,IY+1)*DX
                 W   = B2*DY+B1*DDY
                 SUM = SUM+W
              ENDIF
            ENDIF
	  ENDDO

          BUF(I) = SUM
	ENDDO
     
        CALL WRTLIN(LUN2,BUF,NSAMP,LAUF)
      ENDDO

998   CLOSE(LUN2)
999   CLOSE(LUN1)

      RETURN
      END

