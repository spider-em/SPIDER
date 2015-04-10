C ++********************************************************************
C                                                                      *
C COPYFROMXPLO    REMARKS BLANK BUG             JAN 2014  Ardean Leith *                                                      *
C                                                                      *
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
C                                                                      *
C COPYFROMXPLO(LUNIN,LUNOUT,UNUSED)                                                                     *
C 
C PURPOSE: CONVERT EXPLOR (CCP4) VOLUME TO SPIDER VOLUME
C                                                                     *
C PARAMETERS:                                                                             *
C                                                                      *
C **********************************************************************

        SUBROUTINE  COPYFROMXPLO(LUNIN,LUNOUT,UNUSED)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        CHARACTER(LEN=MAXNAM)   :: FILNAM

        REAL, ALLOCATABLE       ::  Q(:)

        CHARACTER(LEN=1)        :: NULL = CHAR(0)
        CHARACTER(LEN=8)        :: CTEXT

        CALL OPAUXFILE(.TRUE.,FILNAM,NULL,LUNIN,0,
     &                 'O','X-PLOR ASCII',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
       
C       FIND FIRST 'REMARKS'
        DO
           READ(LUNIN,821,END=9999,ERR=9999) CTEXT
           !write(6,*) ' text:',ctext

821        FORMAT(A8)
           IF ( INDEX(CTEXT,'REMARKS') > 0 ) EXIT
        ENDDO

C       FIND LAST 'REMARKS'
        DO
           READ(LUNIN,821,END=9999,ERR=9999) CTEXT
           IF ( INDEX(CTEXT,'REMARKS') <= 0 ) EXIT
        ENDDO
        BACKSPACE LUNIN

        READ(LUNIN,822,ERR=9999) I1,I2,I3, J1,J2,J3, K1,K2,K3
822     FORMAT(9(I8))
        NX   = I3-I2+1
        NY   = J3-J2+1
        NZ   = K3-K2+1

        ALLOCATE (Q(NX*NY), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'Q',NX*NY)
           GOTO 9999
        ENDIF

        IFORM = 3
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNOUT,'U',IFORM,NX,NY,NZ,
     &              MAXIM,'SPIDER OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO9999

        READ(LUNIN,206,IOSTAT=IERR)  DX,DY,DZ
206     FORMAT(6(1PE12.5))
           IF (IERR .NE. 0) THEN
              CALL ERRT(101,'READING X-PLOR ASCII FILE',NDUM)
              GOTO 9999
           ENDIF

        READ(LUNIN,*)

        DO   K=1,NZ
           READ(LUNIN,823,IOSTAT=IERR)  ISL
823        FORMAT(I8)
           IF (IERR .NE. 0) THEN
              CALL ERRT(101,'READING X-PLOR ASCII FILE',NDUM)
              GOTO 9999
           ENDIF

           IF ( (K-1) .NE. ISL)  THEN
              CALL ERRT(102,'WRONG SLICE',K)
              GOTO 9999
           ENDIF

           READ(LUNIN,824,IOSTAT=IERR) (Q(I),I=1,NX*NY)
824        FORMAT(6E12.5)
           IF (IERR .NE. 0) THEN
              CALL ERRT(101,'READING X-PLOR ASCII FILE',NDUM)
              GOTO 9999
           ENDIF
              
           DO  J=1,NY
              NREC = J+(K-1)*NY
              CALL WRTLIN(LUNOUT,Q(1+(J-1)*NX),NX,NREC)
           ENDDO
        ENDDO

9999    CLOSE(LUNIN)
        CLOSE(LUNOUT)

        IF (ALLOCATED(Q)) DEALLOCATE (Q)

        END
