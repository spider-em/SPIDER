
C***********************************************************************
C
C  POSTRT.FOR        CREATED                     AUG 88  ARDEAN LEITH
C                    UPDATED                     JUL 92 ARDEAN LEITH
C                    UPDATED                     MAR 99   ARDEAN LEITH
C                    LINEWIDTH BUG               SEP 11   ARDEAN LEITH
C                    MAXNAM                      JUL 14 ARDEAN LEITH
C
C***********************************************************************
C * AUTHOR:   ARDEAN LEITH
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
C  POSTRT
C
C       PURPOSE: CONTAINS PLOTTING COMMANDS FOR POSTSCRIPT OUTPUT
C
C       ENTRY POINTS ARE:
C
C             POSTRT(LUNPOS)
C                  START PLOT OUPTUT ON UNIT NO. LUNPOS
C    
C             POMOVE(LUNPOS,X,Y)
C                  MOVE PEN TO (X,Y) WITHOUT DRAWING
C
C             PODRAW(LUNPOS,X,Y)
C                  DRAW A LINE FROM CURRENT PEN POSITION TO (X,Y)
C
C             POEND(LUNPOS)
C                  TO END A PLOTTING SEQUENCE
C  
C             POFONT(LUNPOS,IFONT)
C                  TO CHANGE CHARACTER FONTS
C
C             POLINE(LUNPOS,INTEN,IPEN,LINTYP)
C                  TO CHANGE LINE TYPE, LINTYP IS LINE TYPE (1-10)
C                  USED FOR BOTH TEXT AND LINES
C
C             POSEG(LUNPOS,X1,Y1,X2,Y2)
C                  DRAW A LINE. 
C
C             POARAY(LUNPOS,IDATA,NDATA,CLOSD,FILLED)
C                  DRAW A POLYGON. (USE POPATH TO START IT)
C
C             POARAYF(LUNPOS,DATA,NDATA,CLOSD,FILLED)
C                  DRAW A POLYGON. (NO POPATH NEEDED)
C
C             POTEX(LUNPOS,TXLINE,NCHAR,X,Y,ITSIZE,ITANGL,JUST)
C                  PRINTS TEXT AT (X,Y) WITH SIZE=ITSIZE, ANGLE=ITANGL,
C                  AND JUSTIFICATION CODE=JUST
C
C             POTEXT(LUNPOS,TEXT)
C                  PRINTS TEXT
C
C             POSIZE(LUNPOS,ITSIZE)
C                  SET TEXT SIZE,  DEFAULT SIZE IS 12 POINT
C
C             POWIND(LUNPOS,WXMIN,WYMIN,WXMAX,WYMAX)
C                   SETS OUTPUT CLIPPING WINDOW
C
C             POPATH(LUNPOS,X,Y)
C                   START NEW PATH. MOVE PEN TO (X,Y) 
C
C      VARIABLES:   LASSIZ IS LAST TEXT SIZE IN POINTS
C                   
C      NOTE:        THIS DRIVER USES THE SAME SCALING AS USED FOR THE 
C                   RASTER TECH. TERMINAL. ORIGIN (0,0) IS CENTER OF SCREEN.
C                   X INCREASES TO RIGHT,  Y INCREASES UP THE SCREEN
C                   SCREEN IS 1280 X 1024 (-640..639,-512..511)
C               
C                   (-640,511)......................(639,511)
C                        :                              :
C                        :                              :
C                        :                              :
C                        :                              :
C                   (-640,-512)......................(639,-512)
C
C
C      NOTE:        SETGREY WILL NOT WORK IN POSTSCRIPT. IT MUST BE
C                   SPELLED SETGRAY!
C
C--*********************************************************************

        SUBROUTINE POSTRT(LUNPOS)

C       SUBROUTINE TO INITIALIZE POSTSCRIPT DISPLAY

C       SAVE IS NEEDED FEB 99 al
        SAVE

        CHARACTER * 80  TXLINE
        CHARACTER * 18  FONT(0:9)
        CHARACTER * 7   LINT(0:9)
        CHARACTER * 1   NULL,YN
        LOGICAL         CLOSD,FILLED

        PARAMETER       (NSIZE=2000)
        INTEGER         IDATA(2*NSIZE)
        DIMENSION       DATA(3,NSIZE)
        INTEGER         LINWID(0:9)
         
        DATA LINT/'[] 0   ','[2] 0  ','[4] 0  ','[4 1] 0','[6 2] 0',
     &            '[] 0   ','[2] 0  ','[] 0   ','[2] 0  ','[] 0   '/

        DATA LINWID/1,1,1,1,1,2,2,3,3,4/

        DATA FONT(0) /'/Times-Roman      '/
        DATA FONT(1) /'/Times-Bold       '/
        DATA FONT(2) /'/Times-Italic     '/
        DATA FONT(3) /'/Helvetica        '/
        DATA FONT(4) /'/Helvetica-Bold   '/
        DATA FONT(5) /'/Helvetica-Oblique'/
        DATA FONT(6) /'/Courier          '/
        DATA FONT(8) /'/Courier-Bold     '/
        DATA FONT(9) /'/Times-Roman      '/

        DATA SCALEK / 0.333333 /

C         CAN HANDLE NEW NON-RASTER TECH SCALING IF LUNPOST IS < 0
          LUNPOST = LUNPOS
          IF (LUNPOS .LE. 0) LUNPOST = -LUNPOS

          LASANG      = 0
          LASTYP      = 0
          LASTEN      = -1
C         SETS DEFAULT CHARACTER FONTS TO SIZE 12 POINTS
          LASSIZ      = 12
C         SETS DEFAULT FONT TO HELVETICA
          LASFON      = 3

          SCALED      = 1.0

          WRITE(LUNPOST,333) 
333       FORMAT(T1,'%!PS-Adobe-1.0')

C         CONVERT RASTER TECH UNITS TO POSTSCRIPT UNITS (1/72 INCH)
C         XCENT = 4.25 * 72  , YCENT = 5.5 * 72,  SCALE = 72.0 / 128.0

C         SET CURRENT FONT AND FONT SIZE
          FLASSIZ = LASSIZ
          WRITE(LUNPOST,98) FONT(LASFON),FLASSIZ

          IF (LUNPOS .GT. 0) THEN 
C            SET ORIGIN AT CENTER, LANDSCAPE ORIENTATION
             WRITE(LUNPOST,*) ' 306 396 translate 90 rotate'
C            SCALE UNITS SO SAME AS RASTER TECH
             WRITE(LUNPOST,*) ' .5625 .5625 scale'
             SCALED = SCALE / 0.5625

C            RESCALE TEXT
             WRITE(LUNPOS,89) SCALED
          ENDIF

          RETURN

        ENTRY POGETSCALE(LUNPOS,SCALET,SCALEDT)
C         GETS SCALING ------------------------------------- POGETSCALE

          SCALET  = SCALE
          SCALEDT = SCALED

          RETURN

 
        ENTRY POSCALE(LUNPOS,XMARGIN,YMARGIN, XMINS,YMINS, XMAXS,YMAXS)
C         SETS SCALING ---------------------------------------- POSCALE

          XSCALE  = (( 8.5 - (2 * XMARGIN)) * 72)  / (XMAXS - XMINS)
          YSCALE  = ((11.0 - (2 * YMARGIN)) * 72)  / (YMAXS - YMINS)

          SCALE   = MIN(XSCALE,YSCALE)

          WRITE(LUNPOST,334) 
334       FORMAT(/,T1,'% -- scale and center output -(unit=1/72 in.)--')
        
C         SET ORIGIN AT LOWER LEFT, PORTRAIT ORIENTATION
C         TRANSLATE SO OUTPUT IS CENTERED

          XT = 0.5 * ((8.5 * 72) - SCALE * (XMAXS - XMINS)) -
     &         0.5 * SCALE * XMINS

          YT = 0.5 * ((11 * 72)  - SCALE * (YMAXS - YMINS)) - 
     &         0.5 * SCALE * YMINS

          WRITE(LUNPOS,394) XT,YT
394       FORMAT(1X,F7.2,1X,F7.2,' translate ')

C         SCALE UNITS TO FIT MARGIN
          WRITE(LUNPOS,399) SCALE,SCALE
399       FORMAT(F8.2,1X,F8.2,' scale')

          WRITE(LUNPOST,335) 
335       FORMAT(T1,'% ------------------------',/)

C         RESCALE TEXT
          WRITE(LUNPOS,89) 1.0 / SCALE
          SCALED = SCALED / SCALE

C         RESCALE LINE WIDTH
          SLINWID = LINWID(LASTYP) * SCALED
          WRITE(LUNPOS,393) SLINWID
393       FORMAT(2X,F7.3,' setlinewidth ')

          RETURN


        ENTRY POWIND(LUNPOS,WXMIN,WYMIN,WXMAX,WYMAX)
C         SETS OUTPUT WINDOW ----------------------------------- POWIND
          WRITE(LUNPOS,90) WXMIN,WYMIN,WXMIN,WYMAX,
     &                    WXMAX,WYMAX,WXMAX,WYMIN
90        FORMAT( '  initclip newpath ',
     &               F8.2,1X,F8.2,' moveto ',
     &               F8.2,1X,F8.2,' lineto ',/,
     &            2X,F8.2,1X,F8.2,' lineto ',
     &               F8.2,1X,F8.2,' lineto closepath clip')
        RETURN


        ENTRY POMOVE(LUNPOS,X,Y)
C..       MOVE PEN TO (X,Y) ------------------------------------ POMOVE
          WRITE(LUNPOS,92)X,Y
92        FORMAT(2X,F8.2,1X,F8.2,' moveto')
        RETURN


        ENTRY POPATH(LUNPOS,X,Y)
C         START NEW PATH. MOVE PEN TO (X,Y) -------------------- POPATH
          WRITE(LUNPOS,97) X,Y
97        FORMAT('  newpath ',F8.2,1X,F8.2,' moveto')
        RETURN


        ENTRY PODRAW(LUNPOS,X,Y)
C         DRAW LINE TO (X,Y) ----------------------------------- PODRAW
          WRITE(LUNPOS,91)X,Y
91        FORMAT(2X,F8.2,1X,F8.2,' lineto')
        RETURN


        ENTRY POSEG(LUNPOS,X1,Y1,X2,Y2)
C         DRAW LINE --------------------------------------------- POSEG
          WRITE(LUNPOS,997) X1,Y1,X2,Y2
997       FORMAT(' newpath ',F8.2,1X,F8.2,' moveto ',F8.2,1X,F8.2,
     &           ' lineto stroke')
        RETURN


        ENTRY POSHOW(LUNPOS,CLOSD,FILLED)
C         FINISH OFF A LINE AND DISPLAY IT --------------------- POSHOW
          IF (CLOSD .OR. FILLED) THEN
             WRITE(LUNPOS,*) ' closepath'
             IF (FILLED) THEN
                 WRITE(LUNPOS,*) ' gsave fill grestore 1.0 setgray'
                 LASTEN = -1
             ENDIF
          ENDIF
          WRITE(LUNPOS,*) ' stroke'
        RETURN


        ENTRY POLINE(LUNPOS,INTEN,IPEN,LINTYP)
C..       SELECT LINE TYPE AND WIDTH --------------------------- POLINE
          IF (LASTYP .NE. LINTYP) THEN

C******************888DEBUG
            IF (LINTYP .LT. 0 .OR. LINTYP .GT. 9) THEN
               WRITE(6,*) 'LINTYP:',LINTYP
               STOP
            ENDIF
C**************************

C           RESCALE LINE WIDTH
            !SLINWID = LINWID(LASTYP) * SCALE sept 2011 al
            SLINWID = LINWID(LASTYP) 

            !WRITE(LUNPOS,93) LINT(LINTYP)(1:7),SLINWIDsept 2011 al
            WRITE(LUNPOS,93) LINT(LINTYP)(1:7)
93          FORMAT(2X,A7,' setdash ')
!93          FORMAT(2X,A7,' setdash ',F8.3,' setlinewidth ')
            LASTYP = LINTYP
          ENDIF
          IF (LASTEN .NE. INTEN) THEN
            IF (INTEN .GE. 0) THEN
              GRAY = 1.0 - ((INTEN + 1) / 10.0)
            ELSE
              GRAY = 1.0 - (-INTEN / 256.0)
            ENDIF
            GRAY = MIN(1.0,GRAY)
            WRITE(LUNPOS,100) GRAY
100         FORMAT(2X,F5.2,' setgray')
            LASTEN = INTEN
          ENDIF
          RETURN


        ENTRY POARAY(LUNPOS,IDATA,NDATA,CLOSD,FILLED)
C..       DRAW A POLYGON --------------------------------------- POARAY

          NDATA2 = NDATA * 2
          ISTOP = 1

C         USE BATCHES OF 400 COORDINATES TO AVOID PS OVERFLOW
7334      IGO   = ISTOP + 399 
          IF (IGO .GT. NDATA2) IGO = NDATA2
          NOWD2 = (IGO - ISTOP + 1) / 2
          WRITE(LUNPOS,95) (IDATA(I),I= IGO,ISTOP,-1)
95        FORMAT(7(I5,1X,I5,1X))
          WRITE(LUNPOS,94) NOWD2
94        FORMAT(2X,I4,' {lineto} repeat')
          IF (IGO .LT. NDATA2) THEN
C            ANOTHER BATCH NEEDED
             ISTOP = IGO + 1
             GOTO 7334
          ENDIF

          IF (CLOSD .OR. FILLED) THEN
             WRITE(LUNPOS,*) ' closepath'
             IF (FILLED) THEN
                 WRITE(LUNPOS,*) ' gsave fill grestore 1.0 setgray'
                 LASTEN = -1
             ENDIF
          ENDIF

          WRITE(LUNPOS,*) ' stroke'
          RETURN

        ENTRY POARAYF(LUNPOS,DATA,NDATA,CLOSD,FILLED)
C         DRAW A POLYGON -------------------------------------- POARAYF

C         START NEW PATH. MOVE PEN TO (X,Y) 
          WRITE(LUNPOS,697) DATA(1,1),DATA(2,1)
697       FORMAT('  newpath ',F8.2,1X,F8.2,' moveto')
c          CALL SETMINMAX(DATA(1,I),DATA(2,I),XMIN,XMAX,YMIN,YMAX)

C         USE BATCHES OF 200 COORDINATES TO AVOID PS OVERFLOW
          ISTOP = 2

6334      IGO   = ISTOP + 200 - 1 
          IF (IGO .GT. NDATA) IGO = NDATA

          WRITE(LUNPOS,695) (DATA(1,I),DATA(2,I),I= IGO,ISTOP,-1)
695       FORMAT(4(F7.2,1X,F7.2,1X))

C         NOWD2 IS NUMBER OF VALUES IN EACH BATCH
          NOWD2 = (IGO - ISTOP + 1) 
          WRITE(LUNPOS,694) NOWD2
694       FORMAT(2X,I4,' {lineto} repeat')

          IF (IGO .LT. NDATA) THEN
C            ANOTHER BATCH NEEDED
             ISTOP = IGO + 1
             GOTO 6334
          ENDIF

          IF (CLOSD .OR. FILLED) THEN
             WRITE(LUNPOS,*) ' closepath'
             IF (FILLED) THEN
                 WRITE(LUNPOS,*) ' gsave fill grestore 1.0 setgray'
                 LASTEN = -1
             ENDIF
          ENDIF

          WRITE(LUNPOS,*) ' stroke'
          RETURN

#ifdef NEVER
          DO I = 2,NDATA
             WRITE(LUNPOS,695) DATA(1,I),DATA(2,I)
695          FORMAT(2X,F8.2,1X,F8.2,' lineto')
c             CALL SETMINMAX(DATA(1,I),DATA(2,I),XMIN,XMAX,YMIN,YMAX)
          ENDDO
#endif

        ENTRY POTEX(LUNPOS,TXLINE,NCHAR,X,Y,ITSIZE,ITANGL,JUST)
C         PRINT TEXT AT (X,Y) ----------------------------------- POTEX

          IPOINT = ITSIZE

          IF (LASSIZ .EQ. 0) THEN
C            ADDING TO OLD FILE
C            SETS DEFAULT CHARACTER FONTS TO SIZE 12 POINTS
             LASSIZ  = 12
             LASFON  = 3
             FITSIZE = ITSIZE
             WRITE(LUNPOS,98) FONT(LASFON),FITSIZE

          ELSEIF (IPOINT .NE. LASSIZ) THEN
            RELSIZ = FLOAT(ITSIZE) / LASSIZ
            LASSIZ = ITSIZE
            WRITE(LUNPOS,89) RELSIZ
89          FORMAT('  currentfont ',F7.3,' scalefont setfont')
          ENDIF

C         MOVE TO UNJUSTIFIED LOCATION 
          WRITE(LUNPOS,96) X,Y
96        FORMAT(2X,F8.2,1X,F8.2,' moveto')

          IF (ITANGL .NE. LASANG) THEN
C            SET TEXT ANGLE
             WRITE(LUNPOS,88) ITANGL
88           FORMAT('  gsave ',I4,' rotate')
          ENDIF

C         PUT TEXT ON STACK
          WRITE(LUNPOS,80) TXLINE(1:NCHAR)
80        FORMAT(1X,'(',A,')')

C         FIND RELATIVE LOCATION FOR JUSTIFICATION
          YR =  -ITSIZE * SCALEK / SCALE 

          IF (JUST .EQ. 0) THEN
C            LEFT JUSTIFIED TEXT CENTERED ON Y
             WRITE(LUNPOS,103) YR
103          FORMAT(' 0 ',F8.2,' rmoveto show')

          ELSE IF (JUST .EQ. 1) THEN
C            CENTERED IN X AND Y TEXT
             WRITE(LUNPOS,101) YR
101         FORMAT(' dup stringwidth pop 2 div neg ',F8.2,
     &             ' rmoveto show')

          ELSEIF (JUST .EQ. 2) THEN
C            RIGHT JUSTIFIED TEXT CENTERED ON Y
             WRITE(LUNPOS,102) YR
102          FORMAT(' dup stringwidth pop neg ',F8.2,' rmoveto show')

          ELSE
C            ??? JUSTIFIED TEXT 
             WRITE(LUNPOS,104)
104          FORMAT(' show')

          ENDIF

#ifdef NEVER
c*******************new stuff
          WRITE(LUNPOS,*) '  currentpoint exch pop dup dup '
          WRITE(LUNPOS,*) '  ymin lt {dup /ymin exch def} if'
          WRITE(LUNPOS,*) '  ymax gt {    /ymax exch def} {pop} ifelse'
          WRITE(LUNPOS,*) '  currentpoint pop dup dup '
          WRITE(LUNPOS,*) '  xmin lt {dup /xmin exch def} if'
          WRITE(LUNPOS,*) '  xmax gt {    /xmax exch def} {pop} ifelse'
C*************************************
#endif

          IF (ITANGL .NE. LASANG) THEN
             WRITE(LUNPOS,*) '  grestore'
             LASANG = ITANGL
          ENDIF

1011      RETURN



        ENTRY POTEXT(LUNPOS,TXLINE,NCHAR)
C         PUT TEXT ON STACK ------------------------------------ POTEXT
          WRITE(LUNPOS,81) TXLINE(1:NCHAR)
81        FORMAT(1X,'(',A,')  show')
          RETURN


        ENTRY POSIZE(LUNPOS,ITSIZE)
C         SET TEXT SIZE,  DEFAULT SIZE IS 12 POINT ------------- POSIZE

          IF (LASSIZ .EQ. 0) THEN
C            ADDING TO OLD FILE
C            SETS DEFAULT CHARACTER FONTS TO SIZE 12 POINTS
             LASSIZ  = 12
             LASFON  = 3
             FITSIZE = ITSIZE
             WRITE(LUNPOS,98) FONT(LASFON),FITSIZE

          ELSEIF (ITSIZE .NE. LASSIZ) THEN
C           CAHNGE SIZE
            RELSIZ = SCALED * FLOAT(ITSIZE) / FLOAT(LASSIZ)
            LASSIZ = ITSIZE
            WRITE(LUNPOS,89) RELSIZ
          ENDIF
          RETURN



        ENTRY POFONT(LUNPOS,IFONT)
C         SETS CHARACTER FONTS --------------------------------- POFONT
          IF (LASSIZ .EQ. 0) THEN
C            ADDING TO OLD FILE,SETS DEFAULT FONT TO SIZE 12 POINTS
             LASSIZ = 12
          ENDIF

          LASFON     = IFONT
          SCALEDSIZE = LASSIZ * SCALED
          WRITE(LUNPOS,98) FONT(LASFON),SCALEDSIZE
98        FORMAT(2X,A13,' findfont ',F8.3,' scalefont setfont')
          RETURN


        ENTRY POEND(LUNPOS)
C         END PLOTTER USE --------------------------------------- POEND
          WRITE(LUNPOS,99) 
99        FORMAT('  showpage')
          RETURN

        ENTRY POZERO(LUNPOS)
C         INITIALIZE THE TEXT EXTENT RECORDING VARIABLES ------- POZERO
          WRITE(LUNPOS,*) '  /xmax -20000 def'
          WRITE(LUNPOS,*) '  /ymax -20000 def'
          WRITE(LUNPOS,*) '  /xmin  20000 def'
          WRITE(LUNPOS,*) '  /ymin  20000 def'
          RETURN

        ENTRY POUL(LUNPOS)
C         GOTO UPPER LEFT OF TEXT EXTENTS ------------------------ POUL
          WRITE(LUNPOS,*) '  xmin ymax moveto'
          RETURN

        ENTRY POLL(LUNPOS) 
C         GOTO LOWER LEFT OF TEXT EXTENTS ------------------------ POLL
          WRITE(LUNPOS,*) '  xmin ymin moveto'
          RETURN

        END

        SUBROUTINE SETMINMAX(X,Y,XMIN,XMAX,YMIN,YMAX)

        XMIN = MIN(XMIN,X)
        XMAX = MAX(XMAX,X)
        YMIN = MIN(YMIN,Y)
        YMAX = MAX(YMAX,Y)
        RETURN
        END

C       PRINT THE POSTSCRIPT FILE  ----------------------------- POPRINT

        SUBROUTINE POPRINT(FILENAME) 

        INCLUDE 'CMBLOCK.INC'

        CHARACTER(LEN=*):: FILENAME
        CHARACTER * 90     LINE
        CHARACTER *1       YN,NULL

        NULL = CHAR(0)

        CALL RDPRMC(YN,NA,.TRUE.,'PRINT NOW? (Y/N)',NULL,IRTFLG)

        IF (YN .NE. 'N' .AND. YN .NE. 'n') THEN
C          THIS WILL HAVE TO BE ALTERED AT DIFFERENT SITES!!!! 

           WRITE(NOUT,*) ' WARNING: SITE SPECIFIC COMMAND IN POSTRT'
           LINE = 'lp ' // FILENAME // NULL
           CALL CSVMS(LINE,.TRUE.,IERR)
           WRITE(NOUT,*) ' '
        ENDIF

        END

#ifdef NEVER

  /Plot {
.
.
.
  } def

  Plotit
  showpage
#endif
#ifdef NEVER
C         INITIALIZE THE TEXT EXTENT RECORDING VARIABLES
          WRITE(LUNPOST,*) '  /xmax -20000 def'
          WRITE(LUNPOST,*) '  /ymax -20000 def'
          WRITE(LUNPOST,*) '  /xmin  20000 def'
          WRITE(LUNPOST,*) '  /ymin  20000 def'

C         SET NEW SCALING VARIABLES
          XMIN = 10E10
          YMIN = XMIN
          XMAX = -XMIN
          YMAX = -YMIN
#endif
