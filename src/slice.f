
C++*********************************************************************
C
C  SLICE.F                  LONGFILENAMES           JAN 89 ARDEAN LEITH
C                           USED OPFILE             NOV 00 ARDEAN LEITH
C                           USED OPFILEC            FEB 03 ARDEAN LEITH
C                           SETPRMB PARAMETERS      MAY 09 ARDEAN LEITH
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
C      SLICE(MAXDIM,LUN1,LUN2,LUNO)
C
C      THIS SUBROUTINE MAKES A PERPENDICULAR SLICE THROUGH NUMT
C      FILES AT THE SPECIFIED ROW, NREC.
C 
C      PARAMETERS:
C        MAXDIM    MAXIMUM BUFFER SPACE
C        LUN1      LOGICAL UNIT NUMBER OF THE OLDER FILE
C        LUN2      LOGICAL UNIT NUMBER OF THE FILE MOST RECENTLY OPENED
C        LUNO      LOGICAL UNIT NUMBER OF OUTPUT FILE
C
C      CALLED BY: UTIL3
C
C--*******************************************************************

	SUBROUTINE SLICE(MAXDIM,LUN1,LUN2,LUNO)

        INCLUDE 'CMLIMIT.INC'

	COMMON BUF(1)

        COMMON/UNITS/LUNC,NIN,NOUT
	COMMON/MASTER/NSAM1,NROW1,IREC,NLABEL,IFORM,IMAMI,FMAX,FMIN,
     &	              AV,SIG,IHIST

        PARAMETER     (NUMTOT = 256)
	DIMENSION     LUNA(2),NREC1(NUMTOT)

        CHARACTER(LEN=MAXNAM) :: FILNAM(NUMTOT),FILOUT,FILPAT
        CHARACTER(LEN=1)      :: NULL

        NULL = CHAR(0)

	LUNA(1) = LUN1
	LUNA(2) = LUN2
	NMOD    = 2
	AVG     = 0.0
	NS      = 1
	NUMT    = 0

        WRITE(NOUT,*) ' SLICING A STACK OF 2D FILES'

        CALL FILSEQP(FILPAT,NLET,NREC1,NUMTOT,NUMT, 
     &      'FILE PREFIX OR TEMPLATE (EG. PIC****)', IRTFLG)
        IF (IRTFLG .EQ. -1) RETURN

200	DO  IFIL =1,NUMT
           CALL FILGET(FILPAT,FILNAM(IFIL),NLET,NREC1(IFIL),IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(19,'SLICE',NE)
              RETURN
           ENDIF
	ENDDO

        CALL FILERD(FILOUT,NLETO,NULL,'OUTPUT',IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(19,'SLICE',NE)
           RETURN
        ENDIF
	
	CALL RDPRMI(NREC,NDUM,NOT_USED,'ROW NUMBER OF SLICE')

        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FILNAM(1),LUNA(1),'O',IFORM,NSAM,NROW,
     &               NSL1, MAXIM,' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (IMAMI.EQ.0)
     &      CALL NORM3(LUNA(1),NSAM,NROW,NSL1,FMAX,FMIN,AV1)
	AVREF = AV1

C       THE EXTRA NSAM COMES FROM LABEL INFORMATION IN OPEN.
	NOFFS1 = 1+ NSAM
	NOFFS2 = NOFFS1 + 3 * NSAM
	NSAM2 = NSAM*2
	NSAM3 = NSAM*3

12  	IF (NREC.EQ.1 .OR. NREC+1.GT.NROW .OR. NSAM*8.GT.MAXDIM) THEN
           CALL ERRT(7,'SLICE',NE)
           RETURN
        ENDIF

C       NA1 CORRESPONDS TO THE 1ST OF 3 LINES OF ONE FILE.
C       NB1 CORRESPONDS TO THE 1ST OF 3 LINES OF THE SECOND FILE.
C       NAOUT CORRESPONDS TO 1 LINE CONTAINING INTERPOLATION BETWEEN THE
C       TWO FILES DESCRIBED ABOVE.

30	NA1 = NOFFS1
	NB1 = NOFFS2
	NAOUT = NB1+3*NSAM

        MAXIM = 0
        NROWT = 2*NUMT-1
        CALL OPFILEC(0,.FALSE.,FILOUT,LUNO,'U',IFORM,NSAM,NROWT,1,
     &                 MAXIM,' ',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
	   CLOSE(LUNA(1))
	   RETURN
        ENDIF
	IRECSV = IREC

C       READ IN THE THREE LINES TO BE USED FOR INTERPOLATION.
   	CALL REDLIN(LUNA(1),BUF(NA1),NSAM,NREC-1)
	CALL REDLIN(LUNA(1),BUF(NA1+NSAM),NSAM,NREC)

	B    = BUF(NA1+NSAM)
	AMAX = B
	AMIN = B
	AVG  = 0.0
        AVS  = 0.0

	DO  KL1 = 1,NSAM
          B   = BUF(NA1+NSAM+KL1-1)
          AVG = AVG + B
          AVS = AVS + B * B
          IF (B .GT. AMAX) AMAX = B
          IF (B .LT. AMIN) AMIN = B
	ENDDO

	CALL WRTLIN(LUNO,BUF(NA1+NSAM),NSAM,1)
	CALL REDLIN(LUNA(1),BUF(NA1+2*NSAM),NSAM,NREC+1)

C       FOR EFFICIENCY, DO THE FOLLOWING CALCULATIONS OUTSIDE LOOP.

	CONS1 = 1. / SQRT(2.)
	CONS2 = 2. + 6./ SQRT(2.)
	CONS3 = 2. + 8./ SQRT(2.)

C       LOOP TO FORM AN INTERPOLATED LINE BETWEEN EACH OF THE NUMT LINES.
	DO  J = 2,NUMT
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FILNAM(J),LUNA(NMOD),'O',IFORM,
     &                 NSAM,NROW,NSL2,MAXIM,' ',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0)  THEN
	      IF (NMOD .NE. 1) CLOSE(LUNA(1))
  	      CLOSE(LUNA(2))
	      CLOSE(LUNO)
	      RETURN
           ENDIF
           IF (IMAMI.EQ.0)
     &         CALL NORM3(LUNA(NMOD),NSAM,NROW,NSL2,FMAX,FMIN,AV)
	   AV2 = AV

C          READ IN THE THREE LINES TO BE USED FOR INTERPOLATION.

    	   CALL REDLIN(LUNA(NMOD),BUF(NOFFS2),NSAM,NREC-1)
	   CALL REDLIN(LUNA(NMOD),BUF(NOFFS2+NSAM),NSAM,NREC)
	   CALL REDLIN(LUNA(NMOD),BUF(NOFFS2+NSAM2),NSAM,NREC+1)

C          INTERPOLATE BETWEEN EACH OF NSAM POINTS.

C          NORMALIZE THE FILES
	   DIFF1 = AVREF - AV1
	   DIFF2 = AVREF - AV2
	   DIFF  = (2*DIFF1)+(2*DIFF2)

C          NORMALIZE LINE STARTING AT NOFFS2+NSAM SO THAT WHEN WE COPY
C          IT INTO THE OUTPUT FILE, IT IS NORMALIZED ALREADY.

	   DO  K = 1,NSAM
              KN = NOFFS2+NSAM+K-1
              BUF(KN) = BUF(KN)+DIFF2
	   ENDDO

C          FOR FIRST POINT,K=1, INTERPOLATE USING ONLY 3 POINTS.
C          NORMALIZE THE VALUES OF NREC-1 AND NREC+1 ONLY.

	   BUF(NAOUT) = (BUF(NOFFS1+NSAM)+BUF(NOFFS2+NSAM)+
     &        CONS1*(BUF(NOFFS1)+BUF(NOFFS1+NSAM+1)+BUF(NOFFS1+
     &        NSAM2)+BUF(NOFFS2)+BUF(NOFFS2+NSAM+1)+
     &        BUF(NOFFS2+NSAM2)+DIFF))/CONS2
	   B   = BUF(NAOUT)
	   AVG = AVG + B
           AVS = AVS + B * B
	   IF (B .GT. AMAX) AMAX = B
	   IF (B .LT. AMIN) AMIN = B

C          LOOP TO INTERPOLATE FOR ALL BUT THE ENDPOINTS.

	   DO  K = 2,NSAM-1
              BUF(NAOUT+K-1) = (BUF(NOFFS1+NSAM+K-1)+BUF(NOFFS2+NSAM+
     &        K-1)+CONS1*(BUF(NOFFS1+K-1)+BUF(NOFFS1+NSAM+K-2)+
     &        BUF(NOFFS1+NSAM+K)+BUF(NOFFS1+NSAM2+K-1)+
     &        BUF(NOFFS2+K-1)+BUF(NOFFS2+NSAM+K-2)+BUF(NOFFS2+
     &        NSAM+K)+BUF(NOFFS2+NSAM2+K-1)+DIFF))/CONS3
	      B   = BUF(NAOUT+K-1)
              AVG = AVG + B
              AVS = AVS + B * B
              IF (B .GT. AMAX) AMAX = B
              IF (B .LT. AMIN) AMIN = B
	   ENDDO

C          INTERPOLATE THE ENDPOINT.
	   BUF(NAOUT+NSAM-1) = (BUF(NOFFS1+NSAM2-1)+BUF(NOFFS2+
     &        NSAM2-1)+CONS1*(BUF(NOFFS1+NSAM-1)+BUF(NOFFS1+
     &        NSAM2-2)+BUF(NOFFS1+NSAM3-1)+BUF(NOFFS2+NSAM-1)+
     &        BUF(NOFFS2+NSAM2-2)+BUF(NOFFS2+NSAM3-1)+DIFF))/CONS2
	   B   = BUF(NAOUT+NSAM-1)
	   AVG = AVG + B
           AVS = AVS + B * B
	   IF (B .GT. AMAX) AMAX = B
	   IF (B .LT. AMIN) AMIN = B

	   CALL WRTLIN(LUNO,BUF(NAOUT),NSAM,J*2-2)
 	   CALL WRTLIN(LUNO,BUF(NOFFS2+NSAM),NSAM,J*2-1)
	   DO  KL2 = 1,NSAM
              B   = BUF(NOFFS2+NSAM+KL2-1)
              AVG = AVG + B
              AVS = AVS + B * B
              IF (B .GT. AMAX) AMAX = B
              IF (B .LT. AMIN) AMIN = B
	   ENDDO

C          CHANGE POINTER, THEN CLOSE OLDEST FILE.
	   IF (NMOD .EQ. 1) GO TO 90
	   NOFFS1 = NB1
	   NOFFS2 = NA1
	   AV1    = AV2
	   NMOD   = 1
	   GO TO 95

90	   NOFFS1 = NA1
	   NOFFS2 = NB1
	   AV1    = AV2
	   NMOD   = 2
95	   CLOSE(LUNA(NMOD))
	ENDDO

	FMIN  = AMIN
	FMAX  = AMAX
	AVE   = AVG / (FLOAT(NSAM) * FLOAT(NUMT*2-1))
        FNALL = FLOAT(NSAM) * FLOAT(NUMT*2-1)
        DTOP  = AVS - AVE * AVE / FNALL
        IF (DTOP .LT. 0.0D0 .OR. FNALL .EQ. 1.0) THEN
C          SQRT OF NEGATIVE NUMBER OR 1 PIXEL IMAGE
           SIG = 0.0
        ELSE
           SIG = SQRT( DTOP / (FNALL - 1.0))
        ENDIF

 	IREC  = IRECSV
	CALL SETPRMB(LUNO, FMAX,FMIN, AVE,SIG)

	RETURN
	END

