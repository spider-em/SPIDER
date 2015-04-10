
C ++********************************************************************
C                                                                      *
C                                                                      *
C                                                                      *
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
C                                                                      *
C	SPOTWT(MODE,LUN1,IXPOS,IYPOS,IXDIA,IYDIA,PA,PB,PC
C	1,NEGC,KOX,KOY,CAVG,CMX,BKG,AAVG,AMX)
C  
C  SPOTWT FINDS THE CORRECTED CENTER AND AVERAGE WEIGHT OF THE  SPOT
C   WITHIN THE WINDOW AREA. THE ACTIVE WINDOW  HAS AN AREA OF IXDIA X
C   IYDIA AND IS LOCATED WITH ITS CENTER AT IXPOS,IYPOS
C   
C   SPOTWT WORKS IN TWO STEPS. (1) FIND THE BEST FIT PLANE TO THE 
C   EDGES OF THE ACTIVE WINDOW. (2) SUBTRACT THIS BACKGROUND PLANE
C   FROM THE INNER WINDOW DATA, FIND THE LOCATION OF THE CENTER
C   POINT IN THE CORRECTED AREA,FIND THE AVERAGE VALUE.
C
C	MODE: IF MODE IS A NEGATIVE # SKIP BACKGD CORRECTION
C	     AND SET PA=PB=PC=0.
C	     IF( IABS(MODE).EQ.1)  CENTER PTS KOX,KOY IS MAX POINT
C	     IF(IABS(MODE).NEQ.1) CENTER IS CENTER OF DENSITY
C 	LUN1:DIFFRACTION DATA FILE
C	IXPOS,IYPOS: CENTER OF WINDOW
C	IXDIA,IYDIA: SIZE OF WINDOW
C	PA,PB,PC: XSLOPE YSLOPE AND OFFSET OF PLANE FIT TO EDGE 
C	NEGC: NUMBER OF CORRECTED WINDOW POINTS WITHIN EDGES < 0
C	KOX,KOY: POSITION OF CENTER IN CORRECTED WINDOW
C	CAVG: AVERAGE OF CORRECTED WINDOW
C       CMX: MAXIMUM OF DENSITY WITHIN BOXED AREA
C	BKG: AVERAGE VALUE OF EDGE POINTS OF BOXED AREA
C	AAVG: AVERAGE VALUE OF ALL POINTS IN AREA UNCORRECTED FOR BACKGD
C	AMX: ABSOLUTE MAXIMUM VALUE OF POINTS WITHIN BOXED AREA
C
C  CALL TREE: 'SP' --> VTIL3 -->   DIFF10  -->  SPOTWT  -->  SOLVE
C
C  *********************************************************************

	SUBROUTINE SPOTWT(MODE,LUN1,IXPOS,IYPOS,IXDIA,IYDIA,PA,
     &       PB,PC,NEGC,KOX,KOY,CAVG,CMX,BKG,AAVG,AMX,NSAM,NROW,CTOT)

C       I DO NOT KNOW IF SAVE IS NEEDED FEB 99 al
        SAVE

	COMMON BUFF(1)

	COMMON/EDGSUM/PN(3,4)
	COMMON/UNITS/LUN,NIN,NOUT

	IXCORN=IXPOS-IXDIA/2
	IYCORN=IYPOS-IYDIA/2

C       AVOID HITTING EDGES OF THE INPUT FILE
	IF(IXCORN.LT.1)IXCORN=1
	IF(IYCORN.LT.1)IYCORN=1
	IF(IXCORN+IXDIA-1.GT.NSAM)IXCORN=NSAM+1-IXDIA
	IF(IYCORN+IYDIA-1.GT.NROW)IYCORN=NROW+1-IYDIA

	BKG=0.
	DO  J=1,4
          DO  I=1,3
	    PN(I,J)=0.
	  ENDDO
	ENDDO

	CALL REDLIN(LUN1,BUFF,NSAM,IYCORN)
	Y=FLOAT(IYCORN)
	DO I=1,IXDIA
           X=FLOAT(IXCORN-1+I)
           B=BUFF(IXCORN+I-1)
	   BKG=B+BKG
	   CALL PLNEDG(X,Y,B)
        ENDDO

	CALL REDLIN(LUN1,BUFF,NSAM,IYCORN-1+IYDIA)
	Y=FLOAT(IYCORN-1+IYDIA)
	DO I=1,IXDIA
          X=FLOAT(IXCORN-1+I)
          B=BUFF(IXCORN+I-1)
          BKG=B+BKG
          CALL PLNEDG(X,Y,B)
        ENDDO

	DO J=2,IYDIA-1
           CALL REDLIN(LUN1,BUFF,NSAM,IYCORN-1+J)
           Y=FLOAT(IYCORN-1+J)
           X=FLOAT(IXCORN)
           B=BUFF(IXCORN)
           BKG=B+BKG
           CALL PLNEDG(X,Y,B)
           X=FLOAT(IXCORN-1+IXDIA)
           B=BUFF(IXCORN-1+IXDIA)
           BKG=B+BKG
           CALL PLNEDG(X,Y,B)
        ENDDO

	IF (MODE .LE. 0) THEN
           PA=0.
           PB=0.
           PC=0.
        ELSE
           CALL SOLVE(PN,3,3)
           PA=PN(1,4)
           PB=PN(2,4)
           PC=PN(3,4)
        ENDIF

C       FIND MAX AND CENTER OF DENSITY POSITIONS OF 
C       THE WINDOW, ALSO COMPUTE AVERAGE DENSITY

	CMX=0.
	AMX=0.
	AAVG=0.
	CAVG=0.
	XBAR=0.
	YBAR=0.
	NEGC=0
	XCNTR=FLOAT(IXCORN+IXDIA/2)
	YCNTR=FLOAT(IYCORN+IYDIA/2)
	DO  J=1,IYDIA
           CALL REDLIN(LUN1,BUFF,NSAM,IYCORN-1+J)

           Y=FLOAT(IYCORN-1+J)

           DO  I=1,IXDIA
             X=FLOAT(IXCORN-1+I)
             BC=BUFF(IXCORN-1+I)
             ABC=BC
             IF (ABC .GE. AMX)AMX=ABC

             CBC=BC-PA*X-PB*Y-PC
             AAVG=BC+AAVG
             CAVG=CBC+CAVG

             XBAR=XBAR+X*CBC
             YBAR=YBAR+Y*CBC

             IF(CBC.LT.0.)NEGC=NEGC+1
             IF(CBC .GE. CMX) THEN
                CMX=CBC
                MKOX=IXCORN-1+I
                MKOY=IYCORN-1+J
             ENDIF
	  ENDDO
	ENDDO

          IF(ABS(CAVG).LT..1E-10) THEN
            WRITE(68,1001) CMX,AMX,CAVG,AAVG,NEGC,IXDIA,IYDIA,
     &                     ABC,CBC,PA,PB,PC,XBAR,YBAR
1001        FORMAT(4F12.6,3I5/7F12.6)
          ENDIF

	XBAR=XBAR/CAVG
	YBAR=YBAR/CAVG
        CTOT = CAVG
	AREACT=FLOAT(IXDIA*IYDIA)
	CAVG=CAVG/AREACT
	AAVG=AAVG/AREACT
	BKG=BKG/FLOAT(2*IXDIA+2*IYDIA-4)

	IF (IABS(MODE) .EQ. 1) THEN
           KOX=MKOX
           KOY=MKOY
        ELSE
           KOX=IFIX(XBAR+.5)
           KOY=IFIX(YBAR+.5)
        ENDIF

	RETURN
	END

