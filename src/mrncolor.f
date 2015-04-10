C ++********************************************************************
C  MRNCOLOR                                                                    *
C                  OPFILEC                         FEB 03 ARDEAN LEITH
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
C   MRNCOLOR                                                            *
C***********************************************************************

	SUBROUTINE MRNCOLOR

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

        COMMON ADUM(80),OUT(2048),BUF(30000)

        CHARACTER(LEN=MAXNAM)   ::   FIL1,FIL2,FIL3,FIL4,FIL5
        COMMON /COMMUN/ FIL1,FIL2,FIL3,FIL4,FIL5

        DIMENSION     ILUNZ(10),ILUNS(10),IB(20)
        DIMENSION     FMAXZ(10),FMINZ(10),FMAXS(10),FMINS(10)

        CHARACTER *1  NULL,YN

        NULL=CHAR(0)

        LUN5=15

        CALL RDPRMC(YN,NLET,.TRUE.,
     &     'CREATE (I)MAGES, (C)OLORTABLE, (B)OTH (=DEFAULT)'
     &     ,NULL,IRTFLG)

        CALL RDPRMI(IMANY,IDUM,NOT_USED,
     &     'NUMBER OF COLORS OR IMAGES IN IMAGE (MAX 10)')

        IF (YN.EQ.'C') GOTO 1000

        CALL FILERD(FIL5,NLET,NULL,'OUTPUT',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        DO  ICOL=1,IMANY
           ILUNZ(ICOL)=ICOL+20
           ILUNS(ICOL)=ICOL+30

           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FIL1,ILUNZ(ICOL),'O',IFORM,
     &                  NSAM,NROW,NSL,MAXIM,'Z-BUF',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

	   IF (IMAMI.NE.1)
     &        CALL NORM3(ILUNZ(ICOL),NSAM,NROW,NSL,FMAX,FMIN,AV)
           FMAXZ(ICOL)=FMAX
           FMINZ(ICOL)=FMIN

           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FIL2,ILUNS(ICOL),'O',IFORM,NSAM1,NROW1,
     &                 NSL,MAXIM,'Z-BUF',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 1


           IF (NSAM1.NE.NSAM .OR. NROW1.NE.NROW) THEN
             CALL ERRT(1,'MRNCOLOR',NDUM)
             GOTO 1
           ENDIF

	   IF (IMAMI.NE.1)
     &        CALL NORM3(ILUNS(ICOL),NSAM1,NROW1,NSL,FMAX,FMIN,AV)
           FMAXS(ICOL)=FMAX
           FMINS(ICOL)=FMIN
	ENDDO
        RANGE = FMAXS(1)-FMINS(1)
        IB(1) = 0
        IB(2) = NSAM
        DO  ICOL=2,IMANY

C          SET BUFFER BOUNDARIES:
           IG      = (ICOL-1)*2+1
           IG2     = IG+1
           IB(IG)  = IB(IG-1)+NSAM
           IB(IG2) = IB(IG)+NSAM

C          FIND LARGEST DENSITY RANGE:
           R = FMAXS(ICOL)-FMINS(ICOL)
           IF (R.GT.RANGE) RANGE=R
	ENDDO

        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FIL5,LUN5,'U',IFORM,NSAM,NROW,NSL,
     &                   MAXIM,' ',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 1

C       SET RANGE OF IMAGES TO LARGEST OF THE 2.
        KSTEP=117./FLOAT(IMANY)-2.

c       PERC=3*RANGE/STEP 
c       D3=IMANY*RANGE+IMANY*PERC

C       PUT THE 2 RANGE PICTURE TOGETHER:

        DO  I=1,NROW

C          READ Z-BUFFERS FOR DECISION OF HIDDEN SURFACES:
           DO  II=1,IMANY

              IBZ=(II-1)*2+1
              IBS=IBZ+1
              CALL REDLIN(ILUNZ(II),BUF(IB(IBZ)+1),NSAM,I)

C             READ S-BUFFERS:
              CALL REDLIN(ILUNS(II),BUF(IB(IBS)+1),NSAM,I)

	   ENDDO

C          NOW MERGE LINE I:

           DO  K=1,NSAM
             IBS=1
             IPERC=1
             IWHICH=1
             SHORT=NSAM

             DO J=1,IMANY

             IBZ=(J-1)*2+1
             IBS=IBZ+1
             VAL=BUF(IB(IBZ)+K)
             IF(VAL.NE.0..AND.VAL.LT.SHORT) THEN
                SHORT=VAL
                IWHICH=IBS
                IPERC=J
             ENDIF

          ENDDO

          OUT(K)=BUF(IB(IWHICH)+K)/RANGE*KSTEP+(IPERC-1)*(KSTEP+2)

	ENDDO

        IF(I.EQ.1) OUT(1)=117

        CALL WRTLIN(LUN5,OUT,NSAM,I)

	ENDDO

1       CONTINUE

C       CLOSE ALL FILES:
        DO  I=1,IMANY
           CLOSE(ILUNZ(I))
           CLOSE(ILUNS(I))
	ENDDO
        CLOSE(LUN5)
        IF (YN .EQ. 'I') GOTO 200

1000    CONTINUE

C       OPEN(UNIT=LUN5,FILE=FIL5,STATUS='NEW',FORM='FORMATTED',ERR=300)
        CALL OPAUXFILE(.TRUE.,FIL5,NULL,LUN5,0,
     &                 'U','COLOR TABLE',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

c       3 BECAUSE OF PERC BEING 3 STEPS LARGE
        COLRANGE=117./FLOAT(IMANY)-2.
        ICOLRANGE=COLRANGE
        COLSTEP=256./(iCOLRANGE)
        ISTART=1

        WRITE(NOUT,*) 
     &    'COLOUR INTENSITIES IN ARBITRARY UNITS.'
        WRITE(NOUT,*) 
     &     'THE HIGHEST INTENSITY WILL BE SET TO 1 AND THE OTHERS'
        WRITE(NOUT,*) 
     &     'NORMALIZED ACCORDINGLY'

        KOUNT=0
        DO  I=1,IMANY
           WRITE(NOUT,110) I
110        FORMAT(1X,I3,'. COLOUR,:')
           CALL RDPRM2(RED,GREEN,NOT_USED,'RED,GREEN$')
           CALL RDPRM(BLUE,NOT_USED,'BLUE$') 
           CMAX=RED
           IF (GREEN.GT.CMAX) CMAX=GREEN
           IF(BLUE.GT.CMAX) CMAX=BLUE 
           RED=RED/CMAX
           GREEN=GREEN/CMAX
           BLUE=BLUE/CMAX
           REDS=COLSTEP*RED
           GREENS=COLSTEP*GREEN
           BLUES=COLSTEP*BLUE
           DO  IC=1,ICOLRANGE 
              IBLUE=(IC-1)*BLUES
              IRED=(IC-1)*REDS
              IGREEN=(IC-1)*GREENS
              KOUNT=KOUNT+1
              WRITE(LUN5,103) IRED,IGREEN,IBLUE,KOUNT+9
103           FORMAT(4I5)
	   ENDDO
C          WRITE 2 LINES MORE OF THIS:
           DO  IC=1,2
              KOUNT=KOUNT+1
              WRITE(LUN5,103) IRED,IGREEN,IBLUE,KOUNT+9
	   ENDDO
	ENDDO
        CLOSE(LUN5)
        GOTO 200

200     CONTINUE

        RETURN
        END
