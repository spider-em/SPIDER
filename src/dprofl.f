
C++*********************************************************************
C
C DPROFL_G.F   -- CREATED FROM DPROFL            DEC 2014 ArDean Leith
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
C    DPROFL_G(IMFILE,NLETI,GPLFILE,NLETC,LUNI,LUNPOS,NX,NY)
C
C    PURPOSE:       PLOTS INTENSITIES IN A LINE ACROSS IMAGE INTO
C                   GNUPLOT FILE
C
C    PARAMETERS:
C         IMFILE    CHAR. VARIABLE CONTAINING IMAGE FILE NAME
C         NLETI     LENGTH OF IMAGE FILE NAME
C         GPLFILE   CHAR. VARIABLE CONTAINING GNUPLOT FILE NAME  (NEW)
C         NLETC     LENGTH OF GNUPLOT FILE NAME
C         LUNI      LOGICAL UNIT NUMBER OF IMAGE FILE
C         NX        NUMBER OF SAMPLES ON AN IMAGE LINE
C         NY        NUMBER OF ROWS IN IMAGE
C
C    CALLED BY:     PLOT1
C
C--*********************************************************************

      SUBROUTINE DPROFL_G(IMFILE,NLETI,GPLFILE,NLETC,LUNI,LUNGPL,
     &                  NX,NY)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      
      CHARACTER (LEN=*)  :: IMFILE,GPLFILE
      INTEGER            :: NLETI,NLETC,LUNI,LUNGPL,NX,NY

      REAL               :: BUF(NBUFSIZ),BUFSUM(NBUFSIZ) ! FROM CMBLOCK

      INTEGER, PARAMETER :: NTOTL = 17008 
      INTEGER            :: ILIST(NTOTL)

      CHARACTER (LEN=80) :: LINE
      CHARACTER (LEN=1)  :: NULL = CHAR(0)
      CHARACTER (LEN=1)  :: YN
      LOGICAL            :: SUMALL,GLOBAL
      INTEGER            :: NDUM,NTOTLN,NC,NLIST,IRTFLG,IT,I,L,LENL
      INTEGER            :: NNUM,NA


        NTOTLN = NTOTL

        IF (NX <= 0 .OR. NY <= 0) THEN
C         NO. COL OR ROW NOT SPECIFIED SO IT IS A BLANK IMAGE
          CALL ERRT(101, 'UNSPECIFIED NX, NY',NDUM)
          RETURN
        ENDIF

 1      CALL RDPRMC(YN,NC,.TRUE.,
     &    'INDIVIDUAL, GLOBAL, OR SUM SCALE PLOT? (I/G/S)',
     &    NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        GLOBAL = (YN == 'G')
        SUMALL = (YN == 'S')

	WRITE(NOUT,92) NY
   92   FORMAT('  FILE HAS: ',I0,'  ROWS')

        NLIST = NTOTL
        CALL RDPRAI(ILIST,NTOTL,NLIST,1,NY,'ROW NUMBER(S)',
     &              NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999


        WRITE(LUNGPL,'(A)') 'set ylabel "INTENSITY"' 
        WRITE(LUNGPL,'(A)') 'set xlabel "COLUMN"' 

        NNUM = 10
C       PUT UP TO NNUM LINE NUMBERS ON TOP OF PLOT -------------------

        NC = 0
        IF (SUMALL) THEN
           LINE(NC+1:NC+7) = 'SUM OF '
           NC = NC + 7
        ENDIF
        IF (NLIST > 1) THEN
           LINE(NC+1:NC+6 ) = 'ROWS: '
           NC = NC + 6
        ELSE
           LINE(NC+1:NC+5)  = 'ROW: '
           NC = NC + 5
        ENDIF

        DO I = 1,NLIST
           IF (I > NNUM) THEN
              WRITE(NOUT,89)NNUM
   89         FORMAT(/,' SORRY: TITLE CAN ONLY HOLD ',I0,' LINES')
              EXIT
           ENDIF

           L = ILIST(I)
           CALL INTTOCHAR(L,LINE(NC+1:),LENL,1)

           NC = NC + LENL 
           IF (I == NLIST) EXIT

           LINE(NC+1:NC+2) = ', '
           NC = NC + 2
        ENDDO

        WRITE(LUNGPL,'(3A)') 'set title "INTENSITY vs COLUMN \\n',
     &                       LINE(1:NC),'"'

        IF (SUMALL) THEN

           BUFSUM(1:NX) = 0.0      ! ARRAY OP

           DO IT = 1,NLIST         ! LOOP OVER ALL LINES IN LIST
              I = ILIST(IT)        ! CURRENT LINE

              CALL REDLIN(LUNI,BUF,NX,I)

              BUFSUM(1:NX) = BUFSUM(1:NX) + BUF(1:NX)      ! ARRAY OP
           ENDDO

           WRITE(LUNGPL,'(3A)') 'plot \\'
           WRITE(LUNGPL,'(A,I0,A,I0,A)') 
     &           '"-" using 1:2 with line notitle'

           WRITE(LUNGPL,94) (I, BUFSUM(I),I = 1,NX)
94         FORMAT(I10,ES12.4)
        
        ELSE
           IF (GLOBAL) THEN
C             SCALE PLOT BASED ON MAX AND MIN Y VALUES IN WHOLE FILE (FMIN, FMAX)

c              WRITE(LUNGPL,'(A,ES12.4,A,ES12.4,A)') 
              WRITE(LUNGPL,'(A,ES12.4,A,1PG12.4,A)') 
     &              'set yrange [',FMIN,':',FMAX,']'
           ENDIF

           WRITE(LUNGPL,'(3A)') 'plot \\'

           DO IT = 1,NLIST         ! LOOP OVER ALL LINES IN LIST
              I = ILIST(IT)        ! CURRENT LINE

              IF (IT < NLIST) THEN
                 WRITE(LUNGPL,'(A,I5,A)') 
     &             '"-" using 1:2 with line title "Line:',I,',", \\'
              ELSE
                 WRITE(LUNGPL,'(A,I5,A)') 
     &             '"-" using 1:2 with line title "Line:',I,'"'
              ENDIF
           ENDDO

           DO IT = 1,NLIST         ! LOOP OVER ALL LINES IN LIST
              I = ILIST(IT)        ! CURRENT LINE

              CALL REDLIN(LUNI,BUFSUM,NX,I)
              WRITE(LUNGPL,94) (I, BUFSUM(I), I = 1,NX)
              WRITE(LUNGPL,'(A,I0,A)') 'end' 

           ENDDO
        ENDIF

        WRITE(NOUT,*)' GRAPH PLACED IN: ',GPLFILE(1:NLETC)

        CALL FLUSHFILE(LUNGPL)

C       QUERY ABOUT PRINT
        CALL RDPRMC(YN,NA,.TRUE.,'DISPLAY PLOT NOW? (Y/N)',NULL,IRTFLG)

        IF (YN .NE. 'N') THEN
c          WRITE(NOUT,*)' WARNING: EXTERNAL SOFTWARE DEPENDENCY'

           LINE = 'gnuplot -persist ' // GPLFILE(1:NLETC)
           CALL CSVMS(LINE,.TRUE.,IRTFLG)
           WRITE(NOUT,*) ' '
        ENDIF

9999    CLOSE(LUNI)
        CLOSE(LUNGPL)

        END


C++*********************************************************************
C
C DPROFL.F   -- CREATED                           JAN 87 ArDean Leith
C               EXTENSIVELY REWRITTEN             JUL 87 ArDean Leith
C               POSTSCRIPT OUTPUT                 JAN 99 ArDean Leith
C
C **********************************************************************
C *  AUTHOR:  ArDean Leith                                                 *
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
C    DPROFL(IMFILE,NLETI,POSFILE,NLETC,LUNI,LUNPOS,NSAM,NROW)
C
C    PURPOSE:       PLOTS INTENSITIES IN A LINE ACROSS IMAGE INTO
C                   POSTSCRIPT FILE
C
C    PARAMETERS:
C         IMFILE    CHAR. VARIABLE CONTAINING IMAGE FILE NAME
C         NLETI     LENGTH OF IMAGE FILE NAME
C         POSFILE    CHAR. VARIABLE CONTAINING CONTOUR FILE NAME (NEW)
C         NLETC     LENGTH OF CONTOUR FILE NAME
C         LUNI      LOGICAL UNIT NUMBER OF IMAGE FILE
C         NSAM      NUMBER OF SAMPLES ON AN IMAGE LINE
C         NROW      NUMBER OF ROWS IN IMAGE
C
C    CALLED BY:     PLOT1
C
C--*********************************************************************

      SUBROUTINE DPROFL(IMFILE,NLETI,POSFILE,NLETC,LUNI,LUNPOS,
     &                  NSAM,NROW)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      PARAMETER      (NSIZE = 2000)
      PARAMETER      (NTOTL = 17008)
      COMMON         FDATA(3,NSIZE),ILIST(NTOTL)

      COMMON /IOBUF/ BUFF(NBUFSIZ)

      CHARACTER *(*) IMFILE,POSFILE
      CHARACTER * 80 LINE
      CHARACTER * 1  NULL
      LOGICAL        SUMALL,GLOBAL
       
      DATA FLTMIN/-10E30/,FLTMAX/10E30/

        NULL   = CHAR(0)

        NTOTLN = NTOTL

        IF (NSAM .LE. 0 .OR. NROW .LE. 0) THEN
C         NO. COL OR ROW NOT SPECIFIED SO IT IS A BLANK IMAGE
          WRITE(NOUT,*) ' *** ERROR: UNSPECIFIED NSAM, NROW'
          RETURN
        ENDIF

 1      CALL RDPRMC(LINE,NC,.TRUE.,
     &    'INDIVIDUAL,  GLOBAL,  OR SUM  SCALE PLOT? (I/G/S)',
     &    NULL,IRTFLG)
        IF (IRTFLG .EQ. -1) THEN
           CLOSE(LUNPOS)
           CLOSE(LUNI)
           RETURN
        ENDIF

        GLOBAL = .FALSE.
        SUMALL = .FALSE.
        IF (LINE(1:1) .EQ. 'G')THEN
C         SCALE PLOT BASED ON MAX AND MIN Y VALUES IN WHOLE FILE (FMIN, FMAX)
          GLOBAL = .TRUE.
          SYMIN = FMIN
          SYMAX = FMAX

        ELSE IF (LINE(1:1) .EQ. 'S') THEN
C         SUM ALL LINES AND SCALE PLOT BASED ON MAX AND MIN Y OF THE SUM
          SUMALL = .TRUE.
          SYMIN = FLTMAX
          SYMAX = FLTMIN
          DO  I=1,NSAM
             BUFF(NSAM+I) = 0.0
          ENDDO

        ELSE
C         SCALE THE PLOT BASED ON MAX AND MIN Y VALUES OF PLOTTED LINES ONLY
          SYMIN = FLTMAX
          SYMAX = FLTMIN
        ENDIF

	WRITE(NOUT,92) NROW
   92   FORMAT('  FILE HAS:',I4,' ROWS')

   18   NLIST = NTOTL
        CALL RDPRAI(ILIST,NTOTL,NLIST,1,NROW,'ROW NUMBER(S)',
     &              NULL,IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 1

        IF (.NOT. GLOBAL) THEN
C          GO THRU ALL REQUESTED LINES TO FIND MIN AND MAX
           DO  IT = 1,NLIST
              I = ILIST(IT)
              CALL REDLIN(LUNI,BUFF,NSAM,I)
              DO  J = 1,NSAM
                 VAL = BUFF(J)
                 IF (SUMALL) THEN
                    BUFF(NSAM+J) = BUFF(NSAM + J) + VAL
                 ELSE
                    IF (VAL .LT. SYMIN) SYMIN = VAL
                    IF (VAL .GT. SYMAX) SYMAX = VAL
                 ENDIF
              ENDDO
	   ENDDO
        ENDIF

        IF (SUMALL) THEN
C          FIND MAX AND MIN Y VALUES OF THE SUMS
           DO  J = 1,NSAM
              VAL = BUFF(J + NSAM)
              IF (VAL .LT. SYMIN) SYMIN = VAL
              IF (VAL .GT. SYMAX) SYMAX = VAL
          ENDDO
        ENDIF

C       INITIALIZE & SET SCALING FOR POSTSCRIPT
        CALL POSTRT(-LUNPOS)
        CALL POSCALE(LUNPOS,1.0,1.0,  -12.0,-7.0,  125.0,102.0)

C       SET POSTSCRIPT LINE TYPES 
        INTEN  = 9
        IPEN   = 0
        LINTYP = 0

C       PLOT AXES
        SXMIN = 0.0
        SXMAX = NSAM
 
        CALL POSAXIS('X',SXMIN,SXMAX,0.0,0.0,120.0,100.0,XFACTR,LUNPOS,
     &  IRTFLG)

        CALL POSAXIS('Y',SYMIN,SYMAX,0.0,0.0,120.0,100.0,YFACTR,LUNPOS,
     &  IRTFLG)

        IF (.NOT. SUMALL) THEN
C          GET THE VALUES AGAIN
           DO  L = 1,NLIST
              I = ILIST(L)
              CALL REDLIN(LUNI,BUFF,NSAM,I)
              LINTYP = MOD(L-1,10)
              CALL POLINE(LUNPOS,INTEN,IPEN,LINTYP)

C****************
C**              IPEN  = MOD(L,10)  DISABLED FOR WORKSTATION
C*****************

              NDATA = 0
              DO J = 1,NSAM
                 NDATA          = NDATA + 1
                 FDATA(1,NDATA) = J * XFACTR
                 FDATA(2,NDATA) = (BUFF(J) - SYMIN) * YFACTR
                 IF (NDATA .GE. NSIZE) THEN
C                   ARRAY FULL, PUSH INTO FILE            
                    CALL POARAYF(LUNPOS,FDATA,NDATA,.FALSE.,.FALSE.)

                    FDATA(1,1) = FDATA(1,NDATA)
                    FDATA(2,1) = FDATA(2,NDATA)
                    NDATA = 1
                 ENDIF
              ENDDO
              IF (NDATA .GT. 0) 
     &           CALL POARAYF(LUNPOS,FDATA,NDATA,.FALSE.,.FALSE.)
	  ENDDO

        ELSE IF (SUMALL) THEN
           LINTYP = 0
           CALL POLINE(LUNPOS,INTEN,IPEN,LINTYP)
           NDATA = 0
           DO  J = 1,NSAM
              NDATA = NDATA + 1
              FDATA(1,NDATA) = J * XFACTR
              FDATA(2,NDATA) = (BUFF(J+NSAM) - SYMIN) * YFACTR
              IF (NDATA .GE. NSIZE) THEN
C                ARRAY FULL, PUSH INTO FILE            
                 CALL POARAYF(LUNPOS,FDATA,NDATA,.FALSE.,.FALSE.)

                 FDATA(1,1) = FDATA(1,NDATA)
                 FDATA(2,1) = FDATA(2,NDATA)
                 NDATA = 1
              ENDIF
	   ENDDO
           IF (NDATA .GT. 0) 
     &        CALL POARAYF(LUNPOS,FDATA,NDATA,.FALSE.,.FALSE.)
        ENDIF          

C     PUT POSTSCRIPT FILENAME AT TOP,  PUT IMAGE FILENAME AT TOP
      JUST = 0
      LINE = POSFILE(1:NLETC) // '      ' // IMFILE(1:NLETI) // 
     &       '.' // DATEXC(1:3) // '            ' // CHAR(0)

      NC = NLETC + NLETI + 12 + 4

      NNUM = 4
C     PUT UP TO NNUM LINE NUMBERS ON TOP OF PLOT
      XPOS = 66.0

      IF (SUMALL) THEN
         LINE(NC+1:NC+7) = 'SUM OF '
         NC = NC + 7
      ENDIF
      IF (NLIST .GT. 1) THEN
         LINE(NC+1:NC+6 ) = 'ROWS: '
         NC = NC + 6
      ELSE
         LINE(NC+1:NC+5)  = 'ROW: '
         NC = NC + 5
      ENDIF

      DO I = 1,NLIST
         IF (I .GT. NNUM) THEN
            WRITE(NOUT,89)NNUM
   89       FORMAT(/,' SORRY: LABEL LIMITED TO ',I3,' LINES')
            GOTO 61
         ENDIF
         L = ILIST(I)
         CALL INTTOCHAR(L,LINE(NC+1:),LENL,1)

         NC = NC + LENL + 2
         IF (I .EQ. NLIST .OR. I .EQ. 10) GOTO 61
         LINE(NC+1:NC+2) = ', '
         NC = NC + 2
      END DO

   61 CONTINUE

C     SET TEXT CHARACTARISTICS FOR LABEL
      ITSIZE = 12
      ITANGL = 0
      JUST   = 0

C     SET TEXT POSITION FOR LABEL
      YPOS = 114.0
      XPOS = 0.0
      NC   = NC - 2
      CALL POTEX(LUNPOS,LINE,NC,XPOS,YPOS, ITSIZE,ITANGL,JUST)

      WRITE(NOUT,*)' GRAPH PLACED IN: ',POSFILE(1:NLETC)

C     CLOSE THE POSTSCRIPT-FILE 
      CALL POEND(LUNPOS)
      CLOSE(LUNPOS)

      RETURN
      END





