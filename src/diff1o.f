
C ++********************************************************************
C                                                                      
C    DIFF1O                                                          
C                   USED OPFILE                    NOV 00 ARDEAN LEITH
C                   OPFILEC                        FEB 03 ARDEAN LEITH
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
C  DIFF1O(IRTFLG)                                                      
C                                                                      
C  PURPOSE:                                                            
C           OPTION: LATTICE :: EICONIX OR PE DENSITOMETER AND LEXIDATA
C           INTERACTIVE SELECTION OF DIFFRACTION SPOTS WITH LATTICE
C           REFINEMENT AND CORRECTION FOR LOCAL BACKGROUND. 
C
C           OPTION: SINGLE :: : EICONIX OR PE DENSITOMETER AND LEXIDATA
C           SELECTION OF A SINGLE SPOT; INFORMATION OUTPUT TO TERMINAL.
C
C           OPTION RING :: EICONIX OR PE DENSITOMETER AND LEXIDATA
C           SELECTION OF UP TO 20 POINTS ON ONE RING OF A POWDER
C           PATTERN AND LEAST-SQUARES FIT TO CENTER AND RADIUS.
C
C           THIS IS REALLY ANCIENT!! al 2014
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C **********************************************************************

        SUBROUTINE DIFF1O(IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON /FRAME/IFRAME
        COMMON /POINTR/IPOINT
        COMMON /LATTRX/XM(3,4),YM(3,4)
        COMMON /EDGSUM/PN(3,4)
        COMMON /SIZE/ISIZ

        DIMENSION     PLIST(7)

        CHARACTER*8   CTIME
        CHARACTER*12  CDATE

        CHARACTER(LEN=MAXNAM) :: FILNAM,DOCNAM

        CHARACTER *3  ANS,ANSP
        CHARACTER     NULL
 
        REAL VDUM
        DATA  VDUM/0.0/

        DATA JKEY/0/
 
        NULL = CHAR(0)
        LUN1   = 10
        LUNDOC = 13
        IDUM=0 

1000    CALL RDPRMC(ANS,NC,.TRUE.,'OPTION L/S/R ',NULL,IRTFLG)
        IF (IRTFLG.NE.0) RETURN
        
        CALL RDPRMI(ISIZ,ID,NOT_USED,'SIZE FACTOR')
        IF (ISIZ.EQ.0) RETURN

        IF (ANS(1:1).EQ.'R') THEN
           CALL WFTCIRC(X1,Y,RR)
           WRITE(NOUT,2001) X1,Y,RR
2001       FORMAT(' THE CENTER OF THIS RING IS AT ',E12.5,1X,E12.5,/,
     &           ', AND THE RADIUS IS ',E12.5)
           WRITE(NDAT,2001) X1,Y,RR
           RETURN

        ELSEIF (ANS(1:1).EQ.'S') THEN
           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &                   MAXIM,'INPUT',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           CALL SINGLE(NSAM,NROW)
           CLOSE(LUN1)
           RETURN

        ELSE IF(ANS(1:1).NE.'L') THEN
           RETURN
        ENDIF

1       MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &                   MAXIM,'INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

910     CONTINUE
        CALL RDPRMC(ANS,NC,.TRUE.,
     &     'DO YOU WANT A WINDOW POSITION DOCUMENT? (Y/N)',
     &     NULL,IRTFLG)

        CALL RDPRMC(ANSP,NC,.TRUE.,
     &     'DO YOU WANT A PATTERSON FUNCTION DOCUMENT? (Y/N)',
     &     NULL,IRTFLG)

C       INITIALIZE LATTIC SUBROUTINE
        CALL LATICE(0,IDUM,IDUM,VDUM,VDUM,VDUM,VDUM,VDUM,
     &              VDUM,VDUM,VDUM,VDUM,VDUM)
 
        CALL RDPRMI(IXDIA,IYDIA,NOT_USED,'WINDOW SIZE')
        IF (IXDIA .EQ. 0) GO TO 1000
        IF (IYDIA .EQ. 0) IYDIA = IXDIA

        CALL RDPRMI(MODE,ID,NOT_USED,
     &   'MODE: MAX (1)/ CNTR OF DEN (2) // NO BK CORR (-)')
333     CONTINUE

C       FIND INITIAL VALUES FOR LATTICE VECTORS

20      CALL FILERD(DOCNAM,NLET,NULL,'REFLECTION DOC.',IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 333

        NOPEN = 0
        DO IKEY = 1,3
          CALL UNSAV(DOCNAM,NOPEN,LUNDOC,IKEY,PLIST,5,LERR,1)
          NOPEN = 1
          IF (LERR .NE. 0) GOTO 20

          IH    = PLIST(2)
          IK    = PLIST(3)
          IXPOS = PLIST(4)*ISIZ
          IYPOS = PLIST(5)*ISIZ

          CALL SPOTWT(MODE,LUN1,IXPOS,IYPOS,IXDIA,IYDIA,PA,PB,PC,
     &         NEGC,KOX,KOY,CAVG,CMX,BKG,AAVG,AMX,NSAM,NROW,CTOT)

          WRITE(NOUT,4029)KOX,KOY,CAVG,CMX
4029      FORMAT(' CENTER AT:',2(I5,2X),' WINDOW AVG,MAX',
     &             F10.5,', ',F10.5)
          SX = FLOAT(KOX)
          SY = FLOAT(KOY)
C          CALL LATICE(1,IH,IK,SX,SY,CAVG,,,,,,,)
          CALL LATICE(1,IH,IK,SX,SY,CAVG,VDUM,VDUM,
     &              VDUM,VDUM,VDUM,VDUM,VDUM)
 

        ENDDO

C       ALL REFLECTIONS OBTAINED
C        CALL LATICE(2,,,,,,CX,CY,AX,AY,BX,BY,ANGLE)
        CALL LATICE(2,IDUM,IDUM,VDUM,VDUM,VDUM,
     &              CX,CY,AX,AY,BX,BY,ANGLE)
 
        WRITE(NOUT,4035)CX,CY,AX,AY,BX,BY,ANGLE
4035    FORMAT(' LAT. CNTR CX,CY:',2(F10.4,', '),' AX,AY:',
     &     2(F10.4,', '),/,' BX,BY:',2(F10.4,', '),' ANGLE:',F7.2)

400     CONTINUE     
C       OPTIMIZATION STEP
        CALL RDPRMI(IRAD,IMN,NOT_USED,'MAX RADIUS,MIN RADIUS')
        IF (IRAD .EQ. 0) GO TO 910

        RADMX  = FLOAT(IRAD)
        RADMIN = FLOAT(IMN)
        IF (RADMX.LE.RADMIN) THEN
           WRITE(NOUT,*) ' MAX RADIUS LESS THAN MIN RADIUS; TRY AGAIN!'
        GO TO 400
        ENDIF
C       LET REFERENCE REFLECTION BE ' S '
        CALL RDPRMI(ISH,ISK,NOT_USED,'REFLN. INDEX FOR REL. INT.')

        DO  IPASS=1,2
          ARAD = SQRT(AX*AX+AY*AY)
          BRAD = SQRT(BX*BX+BY*BY)
          NHR  = RADMX/ARAD
          NKR  = RADMX/BRAD
C          CALL LATICE(0,,,,,,,,,,,,)
          CALL LATICE(0,IDUM,IDUM,VDUM,VDUM,VDUM,VDUM,VDUM,
     &              VDUM,VDUM,VDUM,VDUM,VDUM)
 
          CALL DATE_2K(CDATE)
          CALL MYTIME(CTIME)

          WRITE(NDAT,432)CDATE,CTIME
432       FORMAT(/'***********   ON  ',A11,'   AT  ',A8,'********',/)

          WRITE(NDAT,4036) IPASS,MODE,IXDIA,IYDIA,IRAD,IMN,CX,CY,
     &                     AX,AY,BX,BY,ANGLE,ISH,ISK
4036      FORMAT('  PASS #:',I2,
     &    ' MODE OF SPOT MEASURING:',I2,' IXDIA,IYDIA:',2I4,
     &    ' MAX,MIN RADIUS ', 2I6,/
     &    '  LATTICE CENTER X,Y:',2(F8.1,' '),'A VECTOR:',2(G10.3,' '),
     &    ' B VECTOR:',2(G10.3,' '),/ ' ANGLE: ',F7.2,
     &    '  INTENSITIES REL. TO  ',2I4,/)

          WRITE(NDAT,4037)
4037      FORMAT('    H   K     R.I.      AVER(COR)    MAX(COR) ',
     &     '   BKG      AVER(UNCOR)  MAX(UNCOR)  CEN X CEN Y',
     &     '  IN X, IN Y,   R,  CTOT')
          SX  = AX*ISH+BX*ISK
          SY  = AY*ISH+BY*ISK
          ISX = IFIX(SX+CX+0.5)
          ISY = IFIX(SY+CY+0.5)
          CALL SPOTWT(MODE,LUN1,ISX,ISY,IXDIA,IYDIA,PA,PB,PC
     &     ,NEGC,KOX,KOY,CAVG,CMX,BKG,AAVG,AMX,NSAM,NROW,CTOT)
          SINT = CAVG
C
          DO  KP=1,NKR+1
            DO 4050 IHP=1,NHR+1

C             COUNTER FOR R.I. AVERAGING
              ICOUNT = 0
              RISUM  = 0.

C             GROUP THE INDICIES BY ABSOLUTE VALUE

              DO 4048 ISC=1,4
                IF (ISC .NE. 1) GO TO 4041
                IH = IHP-1
                IK = KP-1

C               CHECK THAT THE INDEXED SPOT IS IN ALLOWED LIMITS, IF NO SKIP IT

                XPOS = AX*IH+BX*IK
                YPOS = AY*IH+BY*IK
                R    = SQRT(XPOS*XPOS+YPOS*YPOS)
                IF (R .LT. RADMX .AND. R .GT. RADMIN) GO TO 4051
C               SPOT OUT OF MIN MAX RADIUS BOUNDS
                GO TO 4050

4041            IF(ISC.NE.2)GO TO 4043
                IH = -IHP+1
                IK = -KP+1
                GO TO 4051

4043            IF (ISC .NE. 3) GO TO 4045
                IF (IH .EQ. 0 .OR. IK .EQ. 0) GO TO 4049
                IH = -IHP+1
                IK = KP-1
                GO TO 4051

4045            IH = IHP-1
                IK = -KP+1

4051            XPOS  = AX*IH+BX*IK
                YPOS  = AY*IH+BY*IK
                IXPOS = IFIX(XPOS+CX+0.5)
                IYPOS = IFIX(YPOS+CY+0.5)

                IF (IXPOS .GT. 0 .AND. IXPOS .LE. NSAM) GO TO 4052

C               INDEXED SPOT OUT OF RANGE OF FILE
                GO TO 4048

4052            IF (IYPOS .LE.0 .OR. IYPOS .GT. NROW) GO TO 4048

C               MEASURE, ENTER INTO LATTICE OPTIMIZATION AND RECORD SPOT DATA
        CALL SPOTWT(MODE,LUN1,IXPOS,IYPOS,IXDIA,IYDIA,PA,PB,PC
     &           ,NEGC,KOX,KOY,CAVG,CMX,BKG,AAVG,AMX,NSAM,NROW,CTOT)
        IF (ANS(1:1).EQ.'Y'.AND.IPASS.EQ.2) THEN
           JKEY = JKEY+1
           PLIST(1) = JKEY
           PLIST(2) = IH
           PLIST(3) = IK
           PLIST(4) = IXPOS
           PLIST(5) = IYPOS
           PLIST(6) = IXDIA
           PLIST(7) = IYDIA
           CALL SAVD(30,PLIST,7,IRTFLG)
           IF(ANSP(1:1).EQ.'Y') THEN
              PLIST(4) = CTOT
              CALL SAVD(31,PLIST,4,IRTFLG)
           ENDIF
         ENDIF

         RELATI = CAVG/SINT
         ICOUNT = ICOUNT+1
         RISUM  = RISUM+RELATI
         RIAVG  = RISUM/ICOUNT
         WRITE(NDAT,4150)IH,IK,RELATI,CAVG,CMX,BKG,AAVG,AMX,KOX,
     &                 KOY,IXPOS, IYPOS,R,CTOT
4150     FORMAT(1X,2I4,2X,6(G10.4,', '),2I5,',',2I5,',',F7.1
     &              ,', ',F8.2)

         CALL LATICE(1,IH,IK,FLOAT(KOX),FLOAT(KOY),CAVG*R,VDUM,VDUM,
     &              VDUM,VDUM,VDUM,VDUM,VDUM)

C        DISTANCE WEIGHTED POSITIONING -----INSTALLED
4048     CONTINUE

C        ******** CHECK COUNTER IF USED
4049     IF (ICOUNT .EQ. 0) GO TO 4050

         WRITE(NDAT,4151) RIAVG
4151     FORMAT('  AVERAGE R. I. =  ',G10.4,/)
4050     CONTINUE
         ENDDO


          WRITE(NDAT,4152) JKEY
4152      FORMAT('  TOTAL NUMBER OF WINDOWS = ',I6)

C         PARAM(35) = JKEY
          CALL REG_SET(34,FLOAT(JKEY),NULL,IRTFLG)
C         THIS PUTS JKEY INTO SPIDER REGISTER X34
  
C         COMPUTE OPTIMUM LATTICE VECTORS FOR NEXT ITERATION
          IF (IPASS .NE. 2)
     &         CALL LATICE(2,IDUM,IDUM,VDUM,VDUM,VDUM,
     &                  CX,CY,AX,AY,BX,BY,ANGLE)
        ENDDO
        CLOSE(30)
        CLOSE(31)

        CLOSE(LUN1)
        GO TO 910

        END





