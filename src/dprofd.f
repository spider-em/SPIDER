
C++*********************************************************************
C
C  DPROFD.FOR -- CREATED                          JAN  87 ArDean Leith
C                MODIFIED                         JUNE 87 ArDean Leith
C                POSTSCRIPT OUTPUT                JAN  99 ArDean Leith
C                MISC                             MAY  13 ArDean Leith
c
C **********************************************************************
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
C    DPROFD(LUNPOST,LUNDOC)
C
C    PARAMETERS:    LUNPOST,LUNDOC    IO UNITS                   (SENT)

C    PURPOSE:       PREPARES PLOT FROM SELECTED COLUMNS OF A
C                   DOCUMENT FILE.
C
C    CALLED BY:     PLOT1
C
C--*******************************************************************

      SUBROUTINE DPROFD(LUNDOC,LUNPOS)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

C     COMMON /IMGMAX/ INUMBR(NIMAX)

      INTEGER, PARAMETER    :: NSIZE = 2000
      REAL                  :: DATA(3,NSIZE)
      
      CHARACTER(LEN=MAXNAM) :: FILDOC,FILPOS,LINE
      CHARACTER(LEN=1)      :: NULL = CHAR(0)
      LOGICAL               :: VSKEY,ERRI2
     
      INCLUDE 'F90ALLOC.INC'
      REAL, POINTER         :: PBUF(:,:)

C     READ DOCUMENT FILE, ALL KEYS AND REGISTERS 
      MAXX   = 0
      MAXY   = 0
      CALL GETDOCDAT('DOCUMENT',.TRUE.,FILDOC,LUNDOC,
     &                 .TRUE.,MAXX, MAXY,PBUF,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

C     GET NAME OF POSTSCRIPT FILE AND OPEN AS SEQUENTIAL FORMATTED
      CALL OPAUXFILE(.TRUE.,FILPOS,'ps',LUNPOS,0,'N',
     &                       'POSTSCRIPT OUTPUT',.TRUE.,IRTFLGT)
      IF (IRTFLGT .NE. 0) GOTO 9999

        WRITE(NOUT,99)
   99   FORMAT(
     &    ' ASSIGN COLUMN(S) OF REGISTER(S) TO BE PLOTTED ON Y AXIS',
     &    ' AND THE COLUMN OF ',/,
     &    ' REGISTERS TO BE PLOTTED ON THE',
     &    ' X AXIS.  IF THE COLUMN FOR THE X AXIS IS ZERO,',/,
     &    ' THEN THE COLUMN OF REGISTERS IS PLOTTED VERSUS REGISTER ',
     &    'KEY.'/)

  13    NLIST  = 6
        NTOTLN = NIMAX
        CALL RDPRAI(INUMBR,NTOTLN,NLIST,0,MAXX-1,
     &     'REGISTER NUMBER(S) PLOTTED ON Y AXIS',NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9998

        IXREG = 1
   14   CALL RDPRIS(IXREG,IDUM,NOT_USED,
     &     'REGISTER NUMBER PLOTTED ON X AXIS',IRTFLG)
        IF (IRTFLG == -1) GOTO 13
        IF (ERRI2(IXREG,IDUM,1,0,MAXX-1,IDUM,IDUM)) GOTO 14
        VSKEY = (IXREG == 0)

        IKEY1 = 1
        IKEY2 = MAXY
   15   CALL RDPRIS(IKEY1,IKEY2,NOT_USED,
     &        'INITIAL AND FINAL KEY NUMBERS',IRTFLG)
        IF (IRTFLG == -1) GOTO 14
        IF (ERRI2(IKEY1,IKEY2,2,1,MAXY,IKEY1,MAXY)) GOTO 15
  
      
C       GO THRU ALL Y REGISTERS TO FIND OVERALL SYMIN AND SYMAX
        ICOL  = INUMBR(1)  
        SYMIN = PBUF( ICOL + 1, IKEY1)
        SYMAX = SYMIN
        DO  N = 1,NLIST
           ICOL = INUMBR(N) 
           DO IKEY=IKEY1,IKEY2
              ICOUNT = PBUF(1,IKEY)
              IF (ICOUNT > 0) THEN
                 VALUE = PBUF(ICOL + 1,IKEY)
                 SYMIN = MIN(SYMIN,VALUE)
                 SYMAX = MAX(SYMAX,VALUE)
              ENDIF
           END DO
        END DO

        IF (VSKEY) THEN
C          X AXIS IS KEY 
           SXMIN = IKEY1
           SXMAX = IKEY2
        ELSE
C          GO THRU X REGISTER TO FIND SXMIN AND SXMAX
           ICOL  = IXREG 
           SXMIN = PBUF(ICOL + 1, IKEY1)
           SXMAX = SXMIN
           DO IKEY=IKEY1,IKEY2
              ICOUNT = PBUF( 1, IKEY)
              IF (ICOUNT > 0) THEN
C                KEY EXISTS
                 VALUE = PBUF( ICOL + 1, IKEY)
                 SXMIN = MIN(SXMIN,VALUE)
                 SXMAX = MAX(SXMAX,VALUE)
              ENDIF
           END DO
        ENDIF

C       INITIALIZE & SET SCALING FOR POSTSCRIPT
        CALL POSTRT(-LUNPOS)
        CALL POSCALE(LUNPOS,1.0,1.0, -12.0,-7.0,  120.0,102.0)

C       PUT AXES IN POSTSCRIPT FILE
        CALL POSAXIS('X',SXMIN,SXMAX,0.0,0.0, 120.0,100.0,XFACTR,LUNPOS,
     &               IRTFLG)
        CALL POSAXIS('Y',SYMIN,SYMAX,0.0,0.0, 120.0,100.0,YFACTR,LUNPOS,
     &               IRTFLG)

C       SET POSTSCRIPT LINE TYPES 
        INTEN  = 9
        IPEN   = 0
        LINTYP = 0

C       GO THRU  REGISTERS AND PLOT CONTENTS

        DO  N = 1,NLIST
           LINTP = MOD(N-1,10)
           CALL POLINE(LUNPOS,INTEN,IPEN,LINTYP)

           ICOL  = INUMBR(N) 
           NDATA = 0

           DO IKEY = IKEY1,IKEY2
              ICOUNT = PBUF( 1, IKEY)
              IF (ICOUNT > 0) THEN
C                KEY EXISTS, GRAPH IT

                 NDATA = NDATA + 1

C                FIND X VALUE
                 IF (VSKEY) THEN
                    DATA(1,NDATA) = (IKEY - SXMIN) * XFACTR
                 ELSE
                    ICOLX  = IXREG
                    VALUEX = PBUF( ICOLX + 1,IKEY)
                    DATA(1,NDATA) = (VALUEX - SXMIN) * XFACTR
                 ENDIF

C                FIND Y VALUE
                 VALUE         = PBUF( ICOL + 1, IKEY)
                 DATA(2,NDATA) = (VALUE - SYMIN) * YFACTR

                 IF (NDATA >= NSIZE) THEN
C                   ARRAY FULL
                    CALL POARAYF(LUNPOS,DATA,NDATA,.FALSE.,.FALSE.)
                    DATA(1,1) = DATA(1,NDATA)
                    DATA(2,1) = DATA(2,NDATA)
                    NDATA     = 1
                 ENDIF
              ENDIF
           ENDDO
           IF (NDATA > 0)
     &        CALL POARAYF(LUNPOS,DATA,NDATA,.FALSE.,.FALSE.)
   	ENDDO

C       PUT POSTSCRIPT AND DOC FILENAME AT TOP
        NLETP = LNBLNKN(FILPOS)
        NLETD = LNBLNKN(FILDOC)
        LINE  = FILPOS(:NLETP) // '      '// FILDOC(:NLETD) // '      ' 
        NC    = NLETP + NLETD + 12

C       PUT UP TO 6 REG NUMBERS ON TOP OF PLOT
        IF (NLIST == 1) THEN
          LINE(NC+1:NC+5) = 'REG: '
          NC = NC + 5
        ELSE
          LINE(NC+1:NC+6) = 'REGS: '
          NC = NC + 6
        ENDIF

        DO  I = 1,NLIST
          CALL INTTOCHAR(INUMBR(I),LINE(NC+1:),LENL,1)
          NC = NC + LENL 
          IF (I .LT. NLIST) THEN
             LINE(NC+1:NC+2) = ', '
             NC = NC + 2
          ENDIF
       ENDDO
           
       LINE(NC+1:NC+6) = '  VS  '
       NC = NC + 6
       IF (VSKEY) THEN
C         WRITE 'KEY'
          LINE(NC+1:NC+3) = 'KEY'
          NC = NC + 3
       ELSE
C         WRITE X REG NUMBER
          CALL INTTOCHAR(IXREG,LINE(NC+1:),LENL,1)
          NC   = NC + LENL
       ENDIF

C      SET TEXT CHARACTARISTICS FOR LABEL
       ITANGL = 0
       ITSIZE = 20
       JUST   = 0

C      SET POSITION FOR LABEL
       YPOS   = 114.0
       XPOS   = -2.0
       CALL POTEX(LUNPOS,LINE,NC,XPOS,YPOS, ITSIZE,ITANGL,JUST)

       WRITE(NOUT,*)' GRAPH PLACED IN: ',FILPOS(1:NLETP)

C      CLOSE THE POSTSCRIPT-FILE 
       CALL POEND(LUNPOS)
9998   CLOSE(LUNPOS)

C      DEALLOCATE DOC. FILE MEMORY
9999   DEALLOCATE(PBUF)

C      QUERY ABOUT PRINT
       CALL POPRINT(FILPOS(1:NLETP))
 
       END


C++*********************************************************************
C
C DPROFD_GPL.F   MODIFIED FROM GNUPLOT OUTPUT   DEC 14 ArDean Leith
c
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
C    DPROFD_GPL(LUNGPL,LUNDOC)
C
C    PARAMETERS:    LUNGPL,LUNDOC    IO UNITS                   (SENT)
C
C    PURPOSE:       PREPARES GNUPLOT SCRIPT FILE FROM SELECTED 
C                   COLUMNS OF A DOCUMENT FILE.
C
C    CALLED BY:     PLOT1
C
C--*********************************************************************

        SUBROUTINE DPROFD_GPL(LUNDOC,LUNGPL)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER               :: LUNDOC,LUNGPL

C       COMMON /IMGMAX/ INUMBR(NIMAX)

      
        CHARACTER(LEN=MAXNAM) :: FILDOC,FILGPL,LINE
        CHARACTER(LEN=MAXNAM) :: YLABEL,XLABEL,TITLE

        CHARACTER(LEN=1)      :: NULL = CHAR(0)
        CHARACTER             :: YN
        INTEGER               :: MAXX,MAXY,IRTFLG,NLETG,LNBLNKN,NLIST
        INTEGER               :: NTOTLN,IXREG,IYREG,IDUM,NOT_USED
        INTEGER               :: IKEY1,IKEY2,NCHARY,NCHARX,NCHART
        INTEGER               :: I,IKEY,ICOUNT,J,NA

        LOGICAL               :: VSKEY,ERRI2
     
        INCLUDE 'F90ALLOC.INC'
        REAL, POINTER         :: PBUF(:,:)

C       READ DOCUMENT FILE, ALL KEYS AND REGISTERS 
        MAXX   = 0
        MAXY   = 0
        CALL GETDOCDAT('DOCUMENT',.TRUE.,FILDOC,LUNDOC,
     &                 .TRUE.,MAXX, MAXY,PBUF,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       GET NAME OF GNUPLOT FILE AND OPEN AS SEQUENTIAL FORMATTED
C       OPEN FORMATTED, SEQUENTIAL FILE FOR GNUPLOT COMMANDS
        CALL OPAUXFILE(.TRUE.,FILGPL,'gpl',LUNGPL,0,'N',
     &                 'GNUPLOT OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        NLETG = lnblnkn(FILGPL)

        WRITE(NOUT,99)
   99   FORMAT(
     &    '  ASSIGN REGISTER(S) TO BE PLOTTED ON Y AND X AXES. ',
     &    'SPECIFY REGISTER ZERO TO PLOT KEYS ON X AXIS',/)

  13    NLIST  = 6
        NTOTLN = NIMAX
        CALL RDPRAI(INUMBR,NTOTLN,NLIST,0,MAXX-1,
     &     'REGISTER NUMBER(S) PLOTTED ON Y AXIS',NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9998

        IXREG = 1
   14   CALL RDPRIS(IXREG,IDUM,NOT_USED,
     &     'REGISTER NUMBER PLOTTED ON X AXIS',IRTFLG)
        IF (IRTFLG == -1) GOTO 13
        IF (ERRI2(IXREG,IDUM,1,0,MAXX-1,IDUM,IDUM)) GOTO 14
        VSKEY = (IXREG == 0)
        IXREG = IXREG + 1

        IKEY1 = 1
        IKEY2 = MAXY
   15   CALL RDPRIS(IKEY1,IKEY2,NOT_USED,
     &        'INITIAL AND FINAL KEY NUMBERS',IRTFLG)
        IF (IRTFLG == -1) GOTO 14
        IF (ERRI2(IKEY1,IKEY2,2,1,MAXY,IKEY1,MAXY)) GOTO 15
  
        IRTFLG = -999     ! NOT UPPERCASE
        CALL RDPRMC(YLABEL,NCHARY,.FALSE.,'Y AXIS LABEL',NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9998

        IRTFLG = -999     ! NOT UPPERCASE
        CALL RDPRMC(XLABEL,NCHARX,.FALSE.,'X AXIS LABEL',NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9998

        IRTFLG = -999     ! NOT UPPERCASE
        CALL RDPRMC(TITLE,NCHART,.FALSE.,'TITLE',NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9998

        IF (NCHARY > 0) 
     &     WRITE(LUNGPL,'(3A)') 'set ylabel "',YLABEL(:NCHARY),'"' 
        IF (NCHARX > 0) 
     &     WRITE(LUNGPL,'(3A)') 'set xlabel "',XLABEL(:NCHARX),'"' 
        IF (NCHART > 0) 
     &     WRITE(LUNGPL,'(3A)') 'set title "' ,TITLE(:NCHART) ,'"'

        WRITE(LUNGPL,'(3A)') 'plot \\'
 
        IF (NLIST == 1) THEN
           WRITE(LUNGPL,'(A,I0,A,I0,A)') 
     &       '"-" using ', IXREG, ':',INUMBR(1)+1, ' with line notitle'
        ELSE
           DO I = 1, NLIST    
              IF (I < NLIST) THEN
                WRITE(LUNGPL,'(A,I0, A,I0, A,I0, A)') 
     &            '"-" using ',IXREG, 
     &            ':', INUMBR(I)+1,
     &            ' with line title "Reg: ',INUMBR(I),'", \\'
              ELSE
                WRITE(LUNGPL,'(A,I0, A,I0, A,I0, A)') 
     &            '"-" using ',IXREG,
     &            ':',INUMBR(I)+1,
     &            ' with line title "Reg: ',INUMBR(I),'"'
              ENDIF
           ENDDO
        ENDIF

C       KLUDGE FOR MULTIPLE LINES, EMBED DOC FILE FOR EACH LINE
        DO I = 1,NLIST

           DO IKEY = IKEY1,IKEY2
              ICOUNT = PBUF( 1, IKEY)

              IF (ICOUNT > 0) THEN
C               KEY EXISTS, GRAPH IT

                 !write(6,*) ' key:',ikey,icount,IKEY1,IKEY2
                 WRITE(LUNGPL,'(I8," ",30(ES13.6," "))') 
     &                 IKEY,(PBUF(J,IKEY),J=2,ICOUNT+1)

              ENDIF
           ENDDO 
           
           WRITE(LUNGPL,'(A,A,A)') 'end' 
           WRITE(LUNGPL,'(A,A,A)') ' ' 
     &                 
        ENDDO            

        WRITE(NOUT,*)' GRAPH PLACED IN: ',FILGPL(1:NLETG)

C       CLOSE THE GNUPLOT-FILE 
9998    CLOSE(LUNGPL)
        CALL FLUSHFILE(LUNGPL)

C       QUERY ABOUT PRINT
        CALL RDPRMC(YN,NA,.TRUE.,'DISPLAY PLOT NOW? (Y/N)',NULL,IRTFLG)

        IF (YN .NE. 'N') THEN
c          WRITE(NOUT,*)' WARNING: EXTERNAL SOFTWARE DEPENDENCY'

           LINE = 'gnuplot -persist ' // FILGPL(1:NLETG)
           CALL CSVMS(LINE,.TRUE.,IRTFLG)
           WRITE(NOUT,*) ' '
        ENDIF

C       DEALLOCATE DOC. FILE MEMORY
9999    IF (ASSOCIATED(PBUF)) DEALLOCATE(PBUF)

        END


