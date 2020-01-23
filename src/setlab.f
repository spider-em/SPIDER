
C++*********************************************************************
C
C SETLAB.F      CREATED                           NOV 87  ArDean Leith
C               LUNRED                            FEB 03  ArDean Leith
C               1PG FORMAT                        NOV 10  ArDean Leith
C               LINE FORMAT, MPISET               NOV 13  ArDean Leith
C               MRC                               JUL 19  ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2019  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C    SETLAB_R    (LUN,IGO,NBUF,VALUES,PRNT,IRTFLG)
C    SETLAB_R_MRC(LUN,IGO,NBUF,VALUES,PRNT,IRTFLG)
C
C    PURPOSE: SETS HEADER PARAMETERS BY BUFFER NUMBER
C             WRITES THE HEADER INTO THE FILE.
C
C    PARAMETERS:
C           LUN          LOGICAL UNIT NUMBER OF FILE 
C           IGO          FIRST BUFFER POSITION TO BE SET
C           NBUF         NUMBER OF BUFFER POSITIONS TO BE SET
C           VALUES       ARRAY FOR BUFFER VALUES TO BE SET
C           PRNT         PRINT-OUT WANTED
C           IRTFLG       ERROR FLAG (-1 ON ENTRY SUPRESSES PRINT-OUT)
C
C--*********************************************************************
  
        SUBROUTINE SETLAB_R(LUN,IGO,NBUF,VALUES,PRNT,IRTFLG)
  
        IMPLICIT   NONE

        INCLUDE 'CMBLOCK.INC' 

        INTEGER           :: LUN,IGO,NBUF,IRTFLG
        REAL              :: VALUES(*)
        LOGICAL           :: PRNT

C       AUTOMATIC ARRAYS    
        REAL              :: OLDVALUES(NBUF)

        INTEGER           :: ISTOP,NVAL,J,I
        INTEGER           :: ICOMM,MYPID,MPIERR

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID
                
C       UPDATE THE HEADER BUFFER

        ISTOP = MIN(256,IGO+NBUF-1)
        NVAL  = ISTOP - IGO + 1

        IF (PRNT .AND. VERBOSE) THEN
C          GET RELEVANT CURRENT  HEADER VALUES (ONLY USED FOR PRNT)
           CALL LUNGETVALS(LUN,IGO,NVAL,OLDVALUES,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF

C       ALTER RELEVANT  HEADER VALUES
        CALL LUNSETVALS(LUN,IGO,NVAL,VALUES,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (PRNT .AND. VERBOSE) THEN
           J = 0
           DO I = IGO,ISTOP
              J = J + 1
              IF (MYPID <= 0) THEN
                 WRITE(NOUT,9999) I,OLDVALUES(J),VALUES(J)
9999             FORMAT('  HEADER LOCATION: ',I3,' CHANGED FROM: ',
     &                     1PG10.3,' TO: ',1PG10.3)
              ENDIF
           ENDDO
        ENDIF

C       WRITE ALTERED HEADER BACK IN THE FILE
        CALL LUNWRTCURHED(LUN,IRTFLG) 

        END


C       ----------- SETLAB_R_MRC -------------------------------------

        SUBROUTINE SETLAB_R_MRC(LUN,IGO,NBUF,VALUES,PRNT,IRTFLG)
  
        IMPLICIT   NONE

        INCLUDE 'CMBLOCK.INC' 

        INTEGER         :: LUN,IGO,NBUF,IRTFLG
        REAL            :: VALUES(*)

C       AUTOMATIC ARRAYS    
        REAL            :: OLDVALUES(NBUF)
        REAL            :: VALUES_ASI(NBUF)

        LOGICAL         :: PRNT
        INTEGER         :: ISTOP,NVAL,J,I
        INTEGER         :: ICOMM,MYPID,MPIERR

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C       UPDATE THE HEADER BUFFER

        ISTOP = MIN(256,IGO+NBUF-1)
        NVAL  = ISTOP - IGO + 1

        IF (PRNT .AND. VERBOSE) THEN
C          GET RELEVANT CURRENT  HEADER VALUES (NEEDED ONLY FOR PRNT)
           CALL LUNGETVALS_R_MRC(LUN,IGO,NVAL,OLDVALUES,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF

C       ALTER RELEVANT  HEADER VALUES
        CALL LUNSETVALS_R_MRC(LUN,IGO,NVAL,VALUES,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (PRNT .AND. VERBOSE) THEN
           J = 0
           DO I = IGO,ISTOP
              J = J + 1
              IF (MYPID <= 0) THEN
                 WRITE(NOUT,9999) I,OLDVALUES(J),VALUES(J)
9999             FORMAT('  HEADER LOCATION: ',I3,' CHANGED FROM: ',
     &                     1PG10.3,' TO: ',1PG10.3)
              ENDIF
           ENDDO
        ENDIF

C       PUSH ALTERED HEADER BACK IN THE FILE
        CALL LUNWRTHED_MRC(LUN,IRTFLG) 

        END


