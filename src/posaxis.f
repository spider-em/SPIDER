
C++*********************************************************************
C
C POSAXIS.FOR -- CSAXIS CREATED JAN 87                 ARDEAN LEITH
C                CONVERTED FROM METAFILE TO PS  FEB 99 ARDEAN LEITH 
C                USED RDPRI1S                   APR 01 ARDEAN LEITH
C                NO 'ENTER'                     MAY 13 ARDEAN LEITH
C
C **********************************************************************
C * AUTHOR: ARDEAN LEITH
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
C    POSAXIS(TYPET,AMIN,AMAX,XORG,YORG,XEND,YEND,FACTR,LUNPOS,IRTFLG
C
C    PURPOSE:       FINDS SCALING FACTOR FOR X OR Y AXIS ON A PLOT.
C                   UNLESS TYPE IS LOWER CASE THE AXIS IS ENTERED IN
C                   A POSTSCRIPT FILE WHICH BE OPENED ON LUNPOS BEFORE
C                   CALLING THIS ROUTINE.
C
C    PARAMETERS:    TYPET   CHAR. VARIABLE FOR AXIS X OR Y        (SENT)
C                           IF LOWERCASE DO NOT OUTPUT THE AXIS TO FILE
C                   AMIN    MIN VALUE TO BE PLOTTED               (SENT)
C                   AMAX    MAX VALUE TO BE PLOTTED               (SENT)
C                   XORG    LOCATION OF X ORIGIN FOR PLOT         (SENT)
C                   YORG    LOCATION OF Y ORIGIN FOR PLOT         (SENT)
C                   XEND    LOCATION OF END OF X AXIS FOR PLOT    (SENT)
C                   YEND    LOCATION OF TOP OF Y AXIS FOR PLOT    (SENT)
C                   FACTR   SCALE FACTOR FOR THIS AXIS DIMENSION  (RET.)
C                   LUNPOS  IO UNIT                               (SENT)
C                   IRTFLG  ERROR RETURN FLAG (0=NORMAL)          (RET.)
C                           IRTFLG  = -1 ON INPUT MEANS SKIP AXIS QUERY
C                           IRTFLG  = -9 ON INPUT MEANS SKIP ALL QUERY
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE POSAXIS(TYPET,AMIN,AMAX,XORG,YORG,XEND,YEND,FACTR,
     &                  LUNPOS,IRTFLG)

      INCLUDE 'CMBLOCK.INC'

      PARAMETER      (NSIZE=20)
      PARAMETER      (LABNO=50)
      DIMENSION      DATA(3,NSIZE),APOSA(LABNO)
      CHARACTER *(*) TYPET
      CHARACTER *20  PROMPT
      CHARACTER *10  LABEL(LABNO),LABELT
     
      CHARACTER * 6  FMT
      CHARACTER * 1  NULL,TYPE
      LOGICAL        XAXIS,AXPRNT,QUERY,OFFIT
      LOGICAL        ERRF2,ERRI2,USEINT,ASK
      INTEGER        LENLA(LABNO)
       
      DATA FLTZER/10E-30/,FLTMAX/10E30/
      DATA PROMPT/'.AXIS OFFSET'/
      DATA FMT/'(    )'/

      NULL       = CHAR(0)

      QUERY      = .TRUE.
C     AVOID NAN ON FDUM
      FDUM       = 0.0

      ASK =  (IRTFLG .NE. -9)

      IF (IRTFLG == -1) THEN
         QUERY = .FALSE.
         IRTFLG = 0

      ELSE
         ITSIZE = 12
         ITANGL = 0
         JUST   = 0
         IRTFLG = 0
      ENDIF

      TYPE = TYPET     
C     PREVENTS CONSTANT CHANGE ERROR

C     DO NOT WORRY ABOUT CONTOUR DIRECTION
      OFFIT  = .FALSE.
      XAXIS  = .TRUE.
      AXPRNT = .TRUE.

C     LOWER CASE MEANS DO NOT PRINT THE AXIS FIND SCALE ONLY
      IF (TYPE(1:1) == 'y') THEN
         XAXIS  = .FALSE.
         AXPRNT = .FALSE.
         TYPE   = 'Y'
      ELSE IF (TYPE(1:1) == 'Y') THEN
         XAXIS  = .FALSE.
      ELSE IF (TYPE(1:1) == 'x')THEN
         AXPRNT = .FALSE.
         TYPE   = 'X'
      ENDIF

      AXOFF = 0.0
      IF (LEN(TYPET) > 1 .AND. TYPET(2:2) == 'O' .AND. ASK) THEN
C         WANT TO OFFSET AXIS FROM ORIGIN
          PROMPT(8:8) =  TYPE(1:1)
104       CALL RDPRM1S(AXOFF,NOT_USED,PROMPT,IRTFLGT)
          IF (IRTFLGT == -1) RETURN
          OFFIT = .TRUE.
      ENDIF                      

10    USEINT = .FALSE.
      IF (MOD(AMIN,1.0) == 0 .AND. MOD(AMAX,1.0) == 0 .AND.
     &    AMIN > -10E5 .AND. AMIN < 10E5 .AND.
     &    AMAX > -10E5 .AND. AMAX < 10E5) THEN
          USEINT = .TRUE.
          IMIN = AMIN
          IMAX = AMAX
          WRITE(NOUT,909)TYPE,IMIN,IMAX
  909     FORMAT(/,2X,A1,' AXIS OF PLOT:',I10,'....',I10)

      ELSEIF(((AMIN .GE. -10E3  .AND. AMIN .LE.-10E-4) .OR.
     &        (AMIN .GE.  10E-4 .AND. AMIN .LE. 10E4 )) .AND.
     &       ((AMAX .GE. -10E3  .AND. AMAX .LE.-10E-4) .OR.
     &        (AMAX .GE.  10E-4 .AND. AMAX .LE. 10E4 ))) THEN
           
    
         WRITE(NOUT,902)TYPE,AMIN,AMAX
  902    FORMAT(/,2X,A1,' AXIS OF PLOT:',F11.5,'....',F11.5)
      ELSE
         WRITE(NOUT,90)TYPE,AMIN,AMAX
   90    FORMAT(/,2X,A1,' AXIS OF PLOT:',1PG11.4,'....',1PG11.4)
      ENDIF

1115  IF (QUERY .AND. ASK) THEN
  105    CALL RDPRM2S(AMIN,AMAX,NOT_USED,
     &     'NEW LOWER, UPPER AXIS BOUNDS OR <CR>',IRTFLGT)

         IF (IRTFLGT == -1) THEN
            IF (OFFIT) THEN
               CALL RDPRM1S(AXOFF,NOT_USED,PROMPT,IRTFLGT)
               IF (IRTFLGT == -1) RETURN
               OFFIT = .TRUE.
	       GOTO 10
	    END IF 
            RETURN
         ENDIF
         IF (ERRF2(AMAX,FDUM,1,AMIN,FLTMAX,FDUM,FDUM)) GOTO 105
      ENDIF           

      IF (AMIN == AMAX) THEN
         FACTR = 1.0
      ELSE IF (XAXIS) THEN
         FACTR = (XEND-XORG) / ABS(AMAX - AMIN)
      ELSE
         FACTR = (YEND-YORG) / ABS(AMAX - AMIN)
      ENDIF

C     RETURN IF AXIS IS NOT TO BE PLOTTED
      IF (.NOT. AXPRNT) RETURN

      FMAGNF  = 1.0
      NTICK   = 5

C     FIND DEFAULT AXIS LABEL UNITS (6 LABELS)
      AUNIT = ABS(AMAX - AMIN) / 5.0

      AUNITT = AUNIT

      IF (MOD(AUNIT,1.0) .NE. 0.0) THEN
         AUNITT = ABS(AMAX - AMIN) / 6.0
         IF (MOD(AUNITT,1.0) .NE. 0.0) THEN
           AUNITT = ABS(AMAX - AMIN) / 4.0
           IF (MOD(AUNITT,1.0) .NE. 0.0) THEN
             AUNITT = ABS(AMAX - AMIN) / 7.0
             IF (MOD(AUNITT,1.0) .NE. 0.0) GOTO 12
           ENDIF
         ENDIF
      ENDIF
      AUNIT = AUNITT
 
12    IF (MOD(AUNIT,1.0) == 0 .AND. 
     &    AUNIT > -10E5 .AND. AUNIT < 10E5) THEN
C        USE INTEGER FORMAT FOR AXIS UNIT PROMPTING
         IUNIT = AUNIT

         WRITE(NOUT,913)TYPE,IUNIT,NTICK
  913    FORMAT(/,'  ',A1,' AXIS LABELED EVERY:',I9,
     &         ' UNITS WITH',I3,' TICKS / LABEL'/)

      ELSEIF ((AUNIT > -10E3  .AND. AUNIT <-10E-4) .OR.
     &    (AUNIT >  10E-4 .AND. AUNIT < 10E4 )) THEN
C        USE FLOATING POINT FORMAT FOR AXIS UNIT PROMPTING

         WRITE(NOUT,912)TYPE,AUNIT,NTICK
  912    FORMAT(/,'  ',A1,' AXIS LABELED EVERY:',F11.5,
     &         ' UNITS WITH',I3,' TICKS / LABEL'/)
      ELSE
C        USE G FORMAT FOR AXIS UNIT PROMPTING

         WRITE(NOUT,91)TYPE,AUNIT,NTICK
   91    FORMAT(/,'  ',A1,' AXIS LABELED EVERY:',1PG11.4,
     &         ' UNITS WITH',I3,' TICKS / LABEL'/)
      ENDIF

      FTICK = NTICK
        
   14   IF (ASK) THEN
          CALL RDPRM2S(AUNIT,FTICK,NOT_USED,
     &    'NEW AXIS LABEL UNIT AND TICKS / LABEL OR <CR>',IRTFLGT)

          IF (IRTFLGT == -1) THEN
             IF (QUERY) GOTO 1115
             IF (OFFIT) THEN
               CALL RDPRM1S(AXOFF,NOT_USED,PROMPT,IRTFLGT)
               IF (IRTFLGT == -1) RETURN
               OFFIT = .TRUE.
	       GOTO 10
	     ENDIF
             RETURN
          ENDIF
        ENDIF

        IF (ERRF2(AUNIT,FTICK,2,FLTZER,FLTMAX,FLTZER,1000.0)) GOTO 14

        NTICK = FTICK
        IF (NTICK .LE. 0) NTICK = 1

C************************ FUTURE ADDITION
C        FTEMP = 1
C        IMAG  = 0
C        IF (FMAGNF > 1.0) THEN
C   18      IF (MOD(FMAGNF,FTEMP) .NE. 0.0) THEN
C              IMAG = IMAG + 1
C              FTEMP = FTEMP * 10
C              IF (FTEMP < FLTMAX) GOTO 18
C           ENDIF
C           WRITE(NOUT,*) ' *** MAG. FACTOR MUST BE POWER OF TEN'
CC           GOTO 12
C        ELSE IF (FMAGNF < 1.0) THEN
C   18      IF (MOD(FMAGNF,FTEMP) .NE. 0.0) THEN
C              IMAG = IMAG + 1
C              FTEMP = FTEMP + 1
CC              IF (FTEMP < FLTMAX) GOTO 18
C           ENDIF
C           WRITE(NOUT,*) ' *** MAG. FACTOR MUST BE POWER OF TEN'
C           GOTO 12
C************************************************

   20   CONTINUE
     
C       ADD  THIS AXIS
        IF (XAXIS) THEN
           DATA(1,1) = XORG
           DATA(2,1) = YORG + AXOFF
           DATA(1,2) = XEND
           DATA(2,2) = DATA(2,1)
        ELSE
           DATA(1,1) = XORG + AXOFF
           DATA(2,1) = YORG 
           DATA(1,2) = DATA(1,1)
           DATA(2,2) = YEND 
        ENDIF
        CALL POSEG(LUNPOS,DATA(1,1),DATA(2,1), DATA(1,2),DATA(2,2))

C       ADD TICK MARKS
        TUNIT = AUNIT * FACTR / NTICK

        IF (XAXIS) THEN
           TPOS      = XORG
           DATA(2,2) = YORG + AXOFF
           TEND      = XEND + 0.001
        ELSE
           TPOS      = YORG
           DATA(1,2) = XORG + AXOFF
           TEND      = YEND + 0.001
        ENDIF

        ITICK = -1
        IPUSH = 0

   22   ITICK = ITICK + 1
          TLEN  = -1.0
          IF (MOD(ITICK,NTICK) == 0) TLEN = -2.5
          IF (XAXIS) THEN
             DATA(1,1) = TPOS
             DATA(1,2) = TPOS
             DATA(2,1) = TLEN + YORG + AXOFF
          ELSE
             DATA(1,1) = TLEN + XORG + AXOFF
             DATA(2,1) = TPOS
             DATA(2,2) = TPOS
          ENDIF
          CALL POSEG(LUNPOS,DATA(1,1),DATA(2,1), DATA(1,2),DATA(2,2))
          
          TPOS = TPOS + TUNIT
        IF (TPOS. LT. TEND) GOTO 22

        CALL POSEG(LUNPOS,DATA(1,1),DATA(2,1), DATA(1,2),DATA(2,2))
        
C       ADD LABELS TO AXIS

        IF (XAXIS) THEN
           APOS = XORG
           JUST = 0
           YPOS = -7.0 + YORG + AXOFF
        ELSE
           APOS = YORG
           JUST = 2
           XPOS = -5.0 + YORG + AXOFF
        ENDIF
        ALAB  = AMIN / FMAGNF
        FUNIT = AUNIT * FACTR

        IPUSH = 1

        NLABEL = 0
        NDEC   = 0
   30   CONTINUE
C         WRITE EACH LABEL
          NLABEL = NLABEL + 1
          IF (NLABEL > LABNO) THEN
             WRITE(NOUT,9045)
9045         FORMAT(' *** SORRY: TOO MANY LABELS ON THIS AXIS'/)
             GOTO 35
          ENDIF

          IF ((ALAB .LE. 1000000.0 .AND. ALAB .GE. -1000000) .AND. 
     &         MOD(ALAB,1.0) == 0.0 .AND. NDEC == 0) THEN
C            USE INTEGERS AS LABEL
             ILAB = ALAB
             WRITE(LABEL(NLABEL)(1:10),8004,ERR=70) ILAB
 8004        FORMAT(I10)
             
          ELSE
C            USE FLOATING POINT NUMBER FOR LABEL
             WRITE(LABELT,8001,ERR=70) ALAB
 8001        FORMAT(1PE10.2)

             READ(LABELT(8:10),8002,ERR=70) NEXP
 8002        FORMAT(I3)     
             LABEL(NLABEL)(1:10) = '          '

            IF (ALAB == 0.0 .AND. 
     &         (USEINT .OR. NDEC == 0)) THEN
C              ZERO IS SPECIAL CASE TO AVOID EXPONENT
               LABEL(NLABEL)(1:10) = '         0'

            ELSEIF (ALAB == 0.0) THEN
C              ZERO IS SPECIAL CASE TO AVOID EXPONENT
               LABEL(NLABEL)(1:10) = '      0.00'

             ELSEIF (NEXP .LE. 3 .AND. NEXP .GE. -2) THEN
C              USE FLOATING POINT FOR LABELS
               NA    = MAX(0,(3 - NEXP - 1))

               NS = 0
               IF (LABELT(6:6) == '0') THEN
                  NS = 1
                  IF (LABELT(5:5) == '0') THEN
                     NS = 2
                  ENDIF
               ENDIF

               NDEC = MAX(NDEC,NA-NS)
               WRITE(FMT(2:5),8005,ERR=70) NDEC
 8005          FORMAT('F5.',I1)
               WRITE(LABEL(NLABEL)(1:10),FMT,ERR=70) ALAB

            ELSE
C              USE G FORMAT LABEL
               WRITE(LABEL(NLABEL)(1:10),8001,ERR=70)ALAB
            ENDIF 
          ENDIF
     
C         REMOVE LEADING BLANKS
          DO MB = 1,10
             IF (LABEL(NLABEL)(MB:MB) .NE. ' ')GOTO 33
          ENDDO

   33     ME = 10
          IF (.NOT. XAXIS) THEN
C           REMOVE TRAILING BLANKS FROM Y AXIS LABELS
            DO  ME = 10,1,-1
              IF (LABEL(NLABEL)(ME:ME) .NE. ' ')GOTO 38
            ENDDO
          ENDIF
  38      LENL = ME - MB + 1
          LENLA(NLABEL) = LENL
          APOSA(NLABEL) = APOS
          LABEL (NLABEL)(1:LENL) = LABEL(NLABEL)(MB:ME)
          GOTO 71

   70     WRITE(NOUT,*) ' FORMATTING ERROR OCCURED:',ALAB
          WRITE(NOUT,*) ' FMT:    ',FMT
         
  71      APOS = APOS + FUNIT
          ALAB = ALAB + AUNIT / FMAGNF
        IF (APOS < TEND) GOTO 30

        NED = 0

35      IF (XAXIS) THEN
          WRITE(NOUT,9349)
9349      FORMAT(/,'  ',
     &    'X AXIS LABEL NUMBER : LABEL (LEFT JUSTIFIED ON PLOT)')
        ELSE
          WRITE(NOUT,9350)
9350      FORMAT(/,'  ',
     &    'Y AXIS LABEL NUMBER : LABEL (RIGHT JUSTIFIED ON PLOT)')
        ENDIF
        DO  I = 1,NLABEL
          LENL = LENLA(I)
          WRITE(NOUT,92) I, LABEL(I)(1:LENL)
92        FORMAT(16X,I4,' : ',A)
        ENDDO
        WRITE(NOUT,*) ' '

37      NED  = MIN((NED + 1),NLABEL)
        NEDT = NED
        WRITE(NOUT,900) NEDT
900     FORMAT(' .LABEL NO. TO EDIT.  USE  0  FOR LABEL:',I3)

        NED = -3
        IF (ASK) THEN
           CALL RDPRI1S(NED,NOT_USED,
     &      ' -1  TO RELIST ALL LABELS,  OR  <CR> TO CONTINUE',IRTFLGT)
        ENDIF

        IF (ASK .AND. IRTFLGT == -1) THEN
C          GOTO PREVIOUS QUESTION
           IF (NED .LE. 1) GOTO 12
           NED = NEDT - 2
           GOTO 37

        ELSEIF (NED < -2) THEN
C          A <CR> RESPONSE MEANS CONTINUE
           GOTO 683

        ELSEIF (NED == -1) THEN
C          LIST ALL LABELS
           NED = NEDT
           GOTO 35     
       
        ELSEIF (NED > 0) THEN
           IF (ERRI2(NED,IDUM,1,1,NLABEL,IDUM,IDUM)) GOTO 37
           LENL = LENLA(NED)
           WRITE(NOUT,93) LABEL(NED)(1:LENL)
93         FORMAT('  OLD LABEL: ',A)
           CALL RDPRMC(LABEL(NED),NCHAR,.FALSE.,'NEW LABEL',
     &                 NULL,IRTFLGT)
           LENLA(NED) = NCHAR
           GOTO 37

        ELSEIF (NED == 0) THEN
           NED  = NEDT
           LENL = LENLA(NED)
           WRITE(NOUT,93) LABEL(NED)(1:LENL)
           CALL RDPRMC(LABEL(NED),NCHAR,.FALSE.,'NEW LABEL',
     &                 NULL,IRTFLGT)
           LENLA(NED) = NCHAR
           GOTO 37
        ENDIF

683      IPUSH = 0
         DO I = 1,NLABEL
           IF (I .GE. NLABEL) IPUSH = 1
           IF (XAXIS) THEN
             XPOS = APOSA(I)
           ELSE
             YPOS = APOSA(I)
           ENDIF
           LENL = LENLA(I)
           CALL POTEX(LUNPOS,LABEL(I)(:LENL),LENL,XPOS,YPOS,
     &                   ITSIZE,ITANGL,JUST)
          
         ENDDO

        END



