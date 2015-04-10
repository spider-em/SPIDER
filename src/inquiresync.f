C++*********************************************************************
C                                                                      
C  INQUIRESYNC.F              NEW ROUTINE:         JUNE 99 ARDEAN LEITH
C                             ADDED DELAY INPUT    AUG  00 ARDEAN LEITH
C                             ADDED VERBOSE        MAR  02 ARDEAN LEITH
C                             ADDED VERBOSE        MAR  02 ARDEAN LEITH
C                             ISECMAXT < 0         JAN  08 ARDEAN LEITH
C                             REMOVE               AUG  10 ARDEAN LEITH
C                             TYPET                MAY  12 ARDEAN LEITH
C                                                                      
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C   INQUIRESYNC    
C
C   PURPOSE:   DETERMINES IF A FILE EXISTS OR IS GONE.  
C              WAITS TILL IT  EXISTS OR IS GONE
C              CAN ALSO DELETE THE FILE AFTERWARDS
C
C--*********************************************************************

        SUBROUTINE INQUIRESYNC(GONE,REMOVE)

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 
 
        CHARACTER(LEN=MAXNAM) :: FILNAM
        LOGICAL               :: EX,ISOPEN,GONE,INCORE,REMOVE
        CHARACTER(LEN=3)      :: TYPET 

        INTEGER, PARAMETER    :: LUN = 80   ! IO UNIT FOR DELETE OPEN

        IWAITED  = 0

        IF (GONE) THEN
          CALL FILERD(FILNAM,NLET,DATEXC,'WAIT FOR DISAPPEARANCE OF~9',
     &                IRTFLG)
        ELSE
          CALL FILERD(FILNAM,NLET,DATEXC,
     &                'WAIT FOR EXISTANCE OF~9',IRTFLG)
        ENDIF
        IF (IRTFLG .NE. 0) RETURN

5       IDELAY   = 1
        ISECMAX  = 9999999
        CALL RDPRIS(IDELAY, ISECMAXT, NOT_USED,
     &      'DELAY INTERVAL AND MAXIMUM WAIT (<0 --> NO ERROR)',IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 5
        IF (IRTFLG .NE. 0) RETURN
        ISECMAX = IABS(ISECMAXT)

C       FIND IF FILE EXISTS, (CAN EVEN BE A STACKED FILE)
10      TYPET = 'FI'
        CALL INQUIREIF1(33,FILNAM,TYPET,EX,ISOPEN,LUNOP,
     &                  INLNED,IMGNUM,IRTFLG)

        IF ((GONE .AND.  EX) .OR. ((.NOT. GONE) .AND. (.NOT. EX))) THEN
C          WAIT AWHILE AND CHECK AGAIN
           IF (IWAITED .GT. ISECMAX) THEN
              IF (ISECMAXT .GE. 0) THEN
                 CALL ERRT(102,'IQ SYNC -- WAIT TIME EXCEEDED',ISECMAX)
              ELSE
                 WRITE(NOUT,*)'IQ SYNC -- WAIT TIME EXCEEDED:',ISECMAX
                 IF (NDAT .NE. NOUT)
     &             WRITE(NOUT,*)'IQ SYNC -- WAIT TIME EXCEEDED',ISECMAX
                 CALL REG_SET_NSEL(1,1,FLOAT(IWAITED),
     &                            0.0, 0.0, 0.0, 0.0,IRTFLG)
              ENDIF
              RETURN
           ENDIF
#if defined(SP_IBMSP3)
           CALL sleep_(IDELAY)
#else
           CALL SLEEP(IDELAY)
#endif
           IWAITED = IWAITED + IDELAY
           GOTO 10
        ENDIF

        NLET = lnblnk(FILNAM)
        IF (GONE) THEN
           WRITE(NOUT,90) IWAITED,FILNAM(1:NLET)
90         FORMAT('  WAITED: ',I9,' SECONDS FOR DISAPPEARANCE OF: ',A)
        ELSE
           WRITE(NOUT,91) IWAITED,FILNAM(1:NLET)
91         FORMAT('  WAITED: ',I9,' SECONDS FOR CREATION OF: ',A)
        ENDIF

         IF (REMOVE) THEN
C          REMOVE THE SIMPLE FILE
           ILOCAT  = INDEX(FILNAM,'@')
           INCORE  = (FILNAM(1:1) .EQ. '_')

           IF (ILOCAT .LE. 0 .AND. .NOT. INCORE ) THEN
C             CAN NOT REMOVE STACK, STACKED FILE, OR INCORE FILE
              IF (VERBOSE) WRITE(NOUT,*) ' REMOVING: ',FILNAM(1:NLET)
              OPEN(LUN,FILE=FILNAM(1:NLET),STATUS='OLD',IOSTAT=IER)
              IF (IER .EQ. 0) CLOSE(LUN,STATUS='DELETE',IOSTAT=IER)
           ENDIF
        ENDIF              
        IF (VERBOSE) WRITE(NOUT,*) ' '

        CALL REG_SET_NSEL(1,1,FLOAT(IWAITED),0.0, 0.0, 0.0, 0.0,IRTFLG)

        RETURN
        END





