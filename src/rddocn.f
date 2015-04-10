
C++*********************************************************************
C
C RDDOCN.F     ADAPTED FROM RDDOCA2Q.F          FEB 1997 ARDEAN LEITH
C              OPENDOC PARAMETERS               DEC 2000 ARDEAN LEITH
C              INCORE OPENDOC                   JUL 2003 ARDEAN LEITH
C              ADDED 'UD M'                     MAY 2009 ARDEAN LEITH
C              ADDED KEYUSED RETURN             AUG 2009 ARDEAN LEITH
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
C   RDDOCN(DOCNAM,NDOC,MAXMIN)
C
C   PURPOSE:  FINDS NUMBER OF COLUMNS, HIGHEST KEY, & NUMBER OF KEYS 
C             IN DOC A FILE  OR
C             FINDS MAX/MIN VALUE IN SPECIFED COLUMN OF DOC FILE. CAN
C             PLACE THESE VALUES IN REGISTERS GIVEN ON COMMAND LINE
C C 
C   PARAMETERS:  DOCNAM    DOC FILE                            (SENT)
C                NDOC      LUN FOR DOC FILE                    (SENT)
C   PARAMETERS:  NDOC      LUN FOR DOC FILE                    (SENT)
C
C   USAGE:       UD N   [maxkeys],[maxcols],[nkeys]
C                UD MAX [maxval],[minval]
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

	SUBROUTINE RDDOCN(DOCNAM,NDOC,MAXMIN)

	INCLUDE 'CMBLOCK.INC' 

        CHARACTER(LEN=*)                :: DOCNAM
        REAL, ALLOCATABLE, DIMENSION(:) :: PLIST
        LOGICAL                         :: NEWFILE,MAXMIN,EMPTY

        CALL SET_MPI(ICOMM,MYPID,IRTFLG) ! SETS ICOMM AND MYMPID

C       GET DOC FILE NAME
        CALL OPENDOC(DOCNAM,.FALSE.,NLET,NDOC,NICDOC,.FALSE.,'',
     &               .TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
  
        IF (MAXMIN) THEN
C          FIND MAX/MIN VALUE IN SPECIFED COLUMN OF DOC FILE.

           ICOL  = 0
           CALL RDPRI1S(ICOL,NOT_USED,
     &            'REGISTER COLUMN TO BE SEARCHED (0 = KEY)',IRTFLG)
           IF (IRTFLG .NE. 0)  GOTO 9998
           IF (ICOL .LT. 0) THEN
              CALL ERRT(102,'ILLEGAL COLUMN',ICOL)
              GOTO 9999
           ENDIF

           IW = MAX(ICOL,1)
           ALLOCATE(PLIST(IW),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'PLIST',IW)
              GOTO 9999
           ENDIF

           VMIN   = HUGE(VMIN)
           VMAX   = -VMIN
           EMPTY  = .TRUE.

           !write(6,*) ' ikeygot,icount,plist(1),irtflg '
 
           DO        ! ENDLESS LOOP
              CALL LUNDOCREDNXT(NICDOC,IKEYGOT,PLIST,IW,
     &                          IGO,ICOUNT,IRTFLG)
              !write(6,*) ikeygot,icount,plist(1),irtflg
 
              IF (IRTFLG .EQ. 1)  THEN
                 CALL ERRT(101,'ERROR READING DOC. FILE',NE)
                 GOTO 9998

              ELSEIF (IRTFLG .EQ. 2 .AND. EMPTY) THEN
C                END OF DOC FILE AND NO KEYS FOUND
                 CALL ERRT(101,'DID NOT FIND ANY KEYS IN DOC. FILE',NE)
                 GOTO 9998

              ELSEIF (IRTFLG .NE. 0) THEN
C                HAVE FINISHED ALL KEYS IN DOC FILE
                 EXIT

              ELSEIF (ICOUNT .LT. ICOL) THEN
                 CALL ERRT(102,'REGISTER MISSING IN DOC. FILE',ICOL)
                 GOTO 9998
              ENDIF
              EMPTY  = .FALSE.

              IF (ICOL .GT. 0) THEN
C                WANT MAX/MIN VALUE IN REGISTER
                 VMAX = MAX(PLIST(ICOL),VMAX)
                 VMIN = MIN(PLIST(ICOL),VMIN)
              ELSE
C                WANT MAX/MIN KEY
                 VMAX = MAX(FLOAT(IKEYGOT),VMAX)
                 VMIN = MIN(FLOAT(IKEYGOT),VMIN)
              ENDIF
           ENDDO

C          SET OPERATION LINE REGISTERS TO VMIN & VMAX
           CALL REG_SET_NSEL(1,2, VMAX,VMIN, 0.0,0.0,0.0, IRTFLG)

           IF (MYPID .LE. 0) THEN
              WRITE(NOUT,*) ' '
              WRITE(NOUT,90) ICOL,VMAX,VMIN
90            FORMAT('  IN COLUMN: ',I7,
     &                '  MAXIMUM=', G12.5,'  MINIMUM=', G12.5,/)
           ENDIF

9998       IF (ALLOCATED(PLIST)) DEALLOCATE (PLIST)
        ELSE
C          FIND NUMBER OF COLUMNS AND HIGHEST KEY IN DOC A FILE.
           CALL LUNDOCINFO(NICDOC,MAXKEYS,MAXCOLS,KEYUSED,.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

C          SET OPERATION LINE REGISTERS TO MAXKEYS, MAXCOLS, & KEYUSED
           CALL REG_SET_NSEL(1,3,FLOAT(MAXKEYS),FLOAT(MAXCOLS),
     &                        FLOAT(KEYUSED),0.0,0.0,IRTFLG)

           IF (MYPID .LE. 0) WRITE(NOUT,*) ' '
        ENDIF

9999    IF (MYPID .LE. 0) CLOSE(NDOC)
	END

