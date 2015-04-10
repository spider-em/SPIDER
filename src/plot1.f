
C++**************************************************** 6/23/80  2/5/81  VAX
C
C    PLOT1.FOR           CNINT.FOR CALLS ADDED    JAN 86 ARDEAN LEITH
C                        CS DRIVERS ALTERED       NOV 86 ARDEAN LEITH
C                        PROFILE PGMS ALTERED     JAN 86 ARDEAN LEITH
C                        VERSATEC CALLS REMOVED   MAR 90 ARDEAN LEITH
C                        CONVERTED TO SUBROUTINE  JUL 87 ARDEAN LEITH
C                        FILENAME CONVERSION      NOV 88 ARDEAN LEITH
C                        CHANGED 4/30/93           MAHIEDDINE LADJADJ     
C                        REMOVED CNTERM          2/05/99 ARDEAN LEITH
C                        MOST OUTPUT NOW POSTSCRIPT 2/99 ARDEAN LEITH
C                        OPFILEC                  FEB 03 ARDEAN LEITH
C                        MAXNAM                   JUL 14 ARDEAN LEITH
C                        DPROFF_GPL ADDED, CASE   DEC 14 ARDEAN LEITH
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
C       PLOT1(MAXDIM)
C
C       COMMANDS:
C       CO              CONTOURS A IMAGE TO POSTSCRIPT 
C       CO S            CONTOURS A IMAGE FILE INTO STERECON INPUT FILES
C       PF              PROFILES A IMAGE TO POSTSCRIPT 
C       PL              PLOTS CONTENTS OF DOC FILE TO POSTSCRIPT
C       PL G            PLOTS CONTENTS OF DOC FILE TO GNUPLOT SCRIPT
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE PLOT1(MAXDIM)

      PARAMETER (MAXREG = 7)
      PARAMETER (MAXKEY = 5000)
      PARAMETER (NSIZE  = 2000)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

C     DATA IS IN BLANK COMMON FOR USE IN DPROFD ??? al 2014
      COMMON DATA(3,NSIZE),DBUF(MAXREG,MAXKEY)

C     COMMON /COMMUN/ IS USED IN CNINT3!

      REAL                  :: PLIST(7)
      CHARACTER(LEN=MAXNAM) :: IMFILE,DOFILE,POSFILE
      CHARACTER * 1         :: ANSW
      CHARACTER * 1         :: NULL = CHAR(0)
      LOGICAL               :: DOCIT

      INTEGER, PARAMETER    :: LUNIM   = 18
      INTEGER, PARAMETER    :: LUNDOC  = 81
      INTEGER, PARAMETER    :: LUNPOST = 70

C     IFUNCS:   CO  PF   PL 

      SELECT CASE(FCHAR(1:2))

      CASE ('CO')    !  --------- IMAGE CONTOURING --------------- 'CO'

C       TRACE CONTOURS INTO POSTSCRIPT FILE
C       OPEN SPIDER IMAGE TYPE FILE AS INPUT

1       MAXIM = 0
        CALL OPFILEC(0,.TRUE.,IMFILE,LUNIM,'O',ITYPE,NX,NY,NZ,
     &             MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       ESTABLISH FMIN AND FMAX IF NOT IN COMMON YET
        IF (IMAMI .NE. 1) 
     &       CALL NORM3(LUNIM,NX,NY,NZ,FMAX,FMIN,AV)

        CALL CNINT3(LUNIM,LUNPOST,NX,NY,NZ,FMIN,FMAX,MAXDIM)


      CASE ('PF')    !  --------- PROFILE ACROSS IMAGE------------'PF'

C       OPEN A SPIDER IMAGE TYPE FILE AS INPUT
2       MAXIM = 0
        CALL OPFILEC(0,.TRUE.,IMFILE,LUNIM,'O',IFORM,
     &               NX,NY,NZ,
     &               MAXIM,'IMAGE',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9000

        SELECT CASE(FCHAR(4:4))

        CASE ('G')  

C          GET NAME FOR GNUPLOT FILE & OPEN AS SEQUENTIAL FORMATTED
           CALL OPAUXFILE(.TRUE.,POSFILE,'gpl',LUNPOST,0,'N',
     &                 'GNUPLOT OUTPUT',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           NLETC = LNBLNKN(POSFILE)
           NLETI = LNBLNKN(IMFILE)

           CALL DPROFL_G(IMFILE,NLETI,POSFILE,NLETC,LUNIM,LUNPOST,NX,NY)

         CASE DEFAULT

C          GET NAME OF POSTSCRIPT FILE & OPEN AS SEQUENTIAL FORMATTED
           CALL OPAUXFILE(.TRUE.,POSFILE,'ps',LUNPOST,0,'N',
     &                    'POSTSCRIPT OUTPUT',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           NLETC = LNBLNKN(POSFILE)
           NLETI = LNBLNKN(IMFILE)

           CALL DPROFL(IMFILE,NLETI,POSFILE,NLETC,LUNIM,LUNPOST,NX,NY)

           CALL POPRINT(POSFILE)
        END SELECT


      CASE ('PL')    !  ---------  PL OPTIONS ------------------'PL --'

         SELECT CASE(FCHAR(4:5))

         CASE ('HI')    !  HISTOGRAM OF DOC FILE -------------- 'PL HI'

            CALL DCHIST(LUNDOC,LUNIM)


         CASE ('F')    !   FIT CURVE TO DOC FILE (NO PLOT!) --- 'PL FIT'

C          OPEN DOC FILE FOR INPUT
           CALL FILERD(DOFILE,NLETD,DATEXC,'DOCUMENT',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000
           IOP = 0
           CALL UNSDAL(DOFILE,IOP,LUNDOC,1,PLIST,1,DBUF,
     &                 MAXKEY,MAXREG, NKEY,LER)
           CLOSE(LUNDOC)

           CALL POLLSQ(DBUF,MAXREG,MAXKEY,LUNDOC)


         CASE ('G')       ! GNUPLOT REGISTER CONTENTS  ---------- 'PL G'

           CALL DPROFD_GPL(LUNDOC,LUNPOST)


         CASE DEFAULT     ! PS PLOT REGISTER CONTENTS  ---------- 'PL'

           CALL DPROFD(LUNDOC,LUNPOST)

         END SELECT 

       CASE DEFAULT

       END SELECT 


9000   CLOSE(LUNIM)
       CLOSE(LUNPOST)

       END
