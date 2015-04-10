
C **********************************************************************
C
C   MAKTAB.FOR  -- CREATED OCT 90
C **********************************************************************
C *  AUTHOR: ArDean Leith 
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
C      MAKTAB(IEQUIV,NEQUIV,TABLE,LASTCLUS,NTAB,IRTFLG)
C
C      PURPOSE:      MAKES A CONNECTIVITY LOOK-UP-TABLE FROM IEQUIV  
C
C      PARAMETERS:   IEQUIV    EQUIVALENCE TABLE            INPUT
C                    NEQUIV    NUMBER OF EQUIV. IN IEQUIV   INPUT
C                    TABLE     NEW CLUSTER MAPPINGS         OUTPUT
C                    LASTCLUS  NUMBER OF CLUSTERS           INPUT
C                    NTAB      FINAL NUMBER OF CLUSTERS     OUPUT
C                    IRTFLG    ERROR FLAG                   OUTPUT
C
C      CALLED BY:    CONINT    MAPDIST
C
C      CALLS:           
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--********************************************************************

       SUBROUTINE MAKTAB(IEQUIV,NEQUIV,TABLE,LASTCLUS,NTAB,IRTFLG)

       INCLUDE 'CMBLOCK.INC'

       REAL      :: TABLE(LASTCLUS)
       INTEGER   :: IEQUIV(2,*)

       LOGICAL   :: DEBUGING,DEBUGING2
       LOGICAL   :: NOTDONE(LASTCLUS),NOTUSED(NEQUIV)
       INTEGER   :: ISTACK(LASTCLUS)      !???DIM???

C      INITIALIZE THE NOTDONE ARRAY FOR EACH POSSIBLE CLUSTER
       DO IT = 1,LASTCLUS
         NOTDONE(IT) = .TRUE.
       ENDDO

C      INITIALIZE ALREADY ARRAY FOR EACH CONNECTION
       DO IT = 1,NEQUIV
         NOTUSED(IT) = .TRUE.
       END DO

C      NTAB IS NEW CONSECUTIVE NUMBER FOR THIS CLUSTER
       NTAB = 0

C*************DEBUGING
C        IT1 = 1048
C        IT2 = 1085
C        IT3 = 1224
C        IT4 = 1469
C        IT5 = 1358
C*************

C      CHECK EQUIVALENCES FOR EACH OF THE CLUSTERS
       DO  ICLUS = 1,LASTCLUS

c*******************
C         DEBUGING = .FALSE.
C            
C         IF (ICLUS .EQ. IT1 .OR. ICLUS .EQ. IT2 .OR.
C     &       ICLUS .EQ. IT3 .OR. ICLUS .EQ. IT4 .OR.
C     &       ICLUS .EQ. IT5) DEBUGING = .TRUE.
C         IF (DEBUGING) THEN
C             WRITE(10,*)'ICLUS,NOTDONE(ICLUS): ',ICLUS,NOTDONE(ICLUS)
C         ENDIF
C*******************

         IF (NOTDONE(ICLUS)) THEN
C          THIS CLUSTER HAS NOT BEEN CONNECTED ANYWHERE YET

C          INCREMENT NTAB FOR NEW CLUSTER NUMBER
           NTAB = NTAB + 1

         call limchk('TABLE',ICLUS,LASTCLUS)
           TABLE(ICLUS) = NTAB
           ISTACK(1)     = ICLUS
           NUMSTK        = 1
           NOWSTK        = 1

 20        IF (NOWSTK .LE. NUMSTK) THEN
C            CHECK FOR CONNECTIONS INVOLVING THIS STACK MEMBER (JCLUS)
         call limchk('ISTACK',NOWSTK,LASTCLUS)
             ISTKNOW = ISTACK(NOWSTK)

             DO 25 JCLUS = 1,NEQUIV
               IF (NOTUSED(JCLUS)) THEN
                 MAXT = IEQUIV(1,JCLUS)
                 MINT = IEQUIV(2,JCLUS)

                 IF (MINT .EQ. ISTKNOW .OR. MAXT .EQ. ISTKNOW) THEN
C                  ISTKNOW HAS A CONNECTION INVOLVING ICLUS CLUSTER  
                   IT = MINT
                   IF (MAXT .NE. ISTKNOW) IT = MAXT

C********************* DEBUGING
C         DEBUGING2 = .FALSE.
C         IF (MINT .EQ. IT1 .OR. MAXT .EQ. IT1 .OR.
C     &       MINT .EQ. IT2 .OR. MAXT .EQ. IT2 .OR.
C     &       MINT .EQ. IT3 .OR. MAXT .EQ. IT3 .OR.
C     &       MINT .EQ. IT4 .OR. MAXT .EQ. IT4 .OR.
C     &       MINT .EQ. IT5 .OR. MAXT .EQ. IT5) DEBUGING2 = .TRUE.
C          IF (DEBUGING2) THEN
C            WRITE(10,*) ' ICLUS,JCLUS,ISTKNOW,IEQUIV:',ICLUS,JCLUS,
C     &                    ISTKNOW,
C     &                    IEQUIV(1,JCLUS),IEQUIV(2,JCLUS)
C         ENDIF
C************************

C                  MAKE SURE THIS CONNECTION IS NOT ALREADY IN STACK
                   DO ISTK = 1,NUMSTK
                      IF (ISTACK(ISTK) .EQ. IT) GOTO 25
                   END DO
 
C                  PUT THE NEW CONNECTED CLUSTER IN THE STACK                 
                   NUMSTK = NUMSTK + 1
         call limchk('ISTACK',NUMSTK,LASTCLUS)
                   ISTACK(NUMSTK) = IT

C                  ISTKNOW IS CONNECTED TO ICLUS
         call limchk('TABLE',IT,LASTCLUS)
                   TABLE(IT) = NTAB

C                  NO NEED TO CHECK IT IN OUTSIDE LOOP
         call limchk('NOTDONE',IT,LASTCLUS)
                   NOTDONE(IT) = .FALSE.

C                  NO NEED TO CHECK THIS CONNECTION IN INSIDE LOOP
         call limchk('NOTUSED',JCLUS,NEQUIV)
                   NOTUSED(JCLUS) = .FALSE.
C****************************
C           IF (DEBUGING2) THEN
C              WRITE(10,*) 'NUMSTK,IT,JCLUS,TABLE(IT):',NUMSTK,IT,JCLUS,
C     &                    TABLE(IT),ISTACK(NUMSTK)
C           ENDIF
C************************
                 ENDIF
               ENDIF
  25 	     CONTINUE

             NOWSTK = NOWSTK + 1
             GOTO 20
           ENDIF 
C          PUT IN FLAG THAT THIS CLUSTER HAS BEEN CONNECTED             
         call limchk('NOTDONE',ICLUS,LASTCLUS)
           NOTDONE(ICLUS) = .FALSE.                  
         ENDIF
       ENDDO

       WRITE(NOUT,90) LASTCLUS,NTAB
90     FORMAT('  Original:',I6, ' clusters mapped to:',I6/)

       RETURN
       END


C      **************************** LIMCHK *************

       SUBROUTINE LIMCHK(CARRAY,I,LIMIT)

       CHARACTER(LEN=*) :: CARRAY
         
       IF (I > LIMIT) THEN
          WRITE(6,90) CARRAY,LIMIT,I
90        FORMAT(' ARRAY: ',A,'  OVERFLOWS WITH: ',I6,'  >',I6)
          STOP
       ENDIF
       END

