C++*********************************************************************
C                                                                      *
C  SETHEAD.F     NEW                            NOV  2010 ArDean Leith *                       
C                                                                      *
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
C***********************************************************************
C
C   SETHEAD(LUN,NX,NY,NZ.IRTFLG)    
C
C   PURPOSE: SETS VALUES OF HEADER VARIABLE(S) BY NAME
C            I/O FILE ALREADY OPEN ON: LUN
C
C   PARAMETERS:  LUN       I/O UNIT                           (SENT)
C                NX..      FILE SIZE                          (SENT)
C                IRTFLG    ERROR FLAG                         (RET.)
C                          0  IS NORMAL
C                          1 INQUIRY WAS NOT AS EXPECTED
C
C***********************************************************************

      SUBROUTINE SETHEAD(LUN,NX,NY,NZ,IRTFLG)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC' 
      INCLUDE 'CMLIMIT.INC' 

      INTEGER                 :: LUN,NX,NY,NZ,IRTFLG 

      CHARACTER(LEN=160)      :: RESPONSE,ARGNOW,MSG
      CHARACTER(LEN=1)        :: NULL = CHAR(0)

      INTEGER, PARAMETER      :: MAXREGNAM = 10
      CHARACTER(LEN=MAXREGNAM):: REGNAME

C     MAXNSEL IS CURRENTLY SAME AS IN REG_SET.F !!!
      INTEGER, PARAMETER      :: MAXNSEL = 24  ! SEARCH & REGISTER LIST
      INTEGER                 :: LOCATION(MAXNSEL)
      INTEGER                 :: IREGSELS(MAXNSEL)
      REAL                    :: VALUES(MAXNSEL)
      CHARACTER(LEN=8)        :: NAMES(MAXNSEL)
      LOGICAL                 :: IS_MRC
      INTEGER                 :: I,IVALS,NGOT,IHEDLOC
      INTEGER                 :: ICOMM,MYPID,MPIERR 
        
      CALL SET_MPI(ICOMM,MYPID,IRTFLG)   ! RETURNS MPI PROCESS ID OR -1

      CALL LUNGETIS_MRC(LUN,IS_MRC,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET LIST OF ALL VARIABLE NUMBERS TO BE SET IN: LOCATION 
      IF (IS_MRC) THEN
         CALL INQUIREHEAD_LOC_MRC(LOCATION,NAMES,IVALS,IRTFLG)
      ELSE
         CALL INQUIREHEAD_LOC(    LOCATION,NAMES,IVALS,IRTFLG)
      ENDIF     
      IF (IRTFLG .NE. 0) RETURN
      ! write(6,*)' Matching values; ',ivals

C     GET LIST OF VARIABLE VALUES IN: VALUES 
      CALL RDPRA('VALUE(S)',IVALS,0,.FALSE.,VALUES,NGOT,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (IVALS .NE. NGOT) THEN
         CALL ERRT(102,'INCONSISTENT # OF VARIABLES',NGOT)
         RETURN
      ENDIF


      DO I=1,IVALS    ! LOOP OVER ALL  WANTED HEADER VALUES
        IHEDLOC = LOCATION(I)     ! HEADER LOCATION WANTED
        !write(6,*)' In sethead - i,values(i): ',i,values(i)

C       SET HEADER VALUE FOR THIS HEADER LOCATION  
        IF (IS_MRC) THEN
C         MRC FILE
          CALL LUNSETVALS_R_MRC(LUN,IHEDLOC,1,VALUES(I),IRTFLG)
        ELSE
C         SPIDER FILE
          CALL SETLAB_R(LUN,IHEDLOC,1,VALUES(I),.FALSE.,IRTFLG)
        ENDIF
        ! write(6,*) ' IHEDLOC: ',IHEDLOC,VALUES(I)
        IF (IRTFLG .NE. 0) RETURN


C       KLUDGE FOR NZ < 0
        IF (LOCATION(I) == 1 .AND. VALUES(I) < 0) 
     &      VALUES(I) = ABS(VALUES(I))

        IF (VERBOSE .AND. (MYPID <= 0) ) THEN
C          ECHO HEADER VALUE & REGISTER SETTING
           ! WRITE(NOUT,91) IHEDLOC,NAME(IHEDLOC),VALUES(I)
91         FORMAT('  HEADER VARIABLE: ',I3,'  NAME: ',A,' = ',1PG11.3)
        ENDIF  ! END OF: IF (VERBOSE .AND. (MYPID <= 0) ) THEN
      ENDDO    ! END OF: DO I=1,IVALS  
  
      IF (IS_MRC) THEN
C       WRITING TO MRC FILE (SPIDER DONE IN SETLAB_R)
        CALL LUNWRTHED_MRC(LUN,IRTFLG)
      ENDIF

      END

C     ---------------- SETHEADEM2 -----------------------------------

C     PURPOSE: ALLOWS FIXING BAD FILES FROM: EM2EM O

      SUBROUTINE SETHEADEM2(LUN,NY,IRTFLG)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC' 
 
      INTEGER :: LUNARA,LUNSTK,LUNARB,LUNFLIP
      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

      INTEGER :: LUN,NY,IRTFLG

      INTEGER :: NRECS,IMGNUM,NOT_USED,LABREC,NDUM,LABBYT,LENBYT
      INTEGER :: IRECOLD
      REAL    :: OLDIMGNUM

      INTEGER :: ICOMM,MYPID

      CALL SET_MPI(ICOMM,MYPID,IRTFLG)   ! RETURNS MPI PROCESS ID OR -1

C     RETRIEVE IMAGE NUMBER
      CALL LUNGETVALS(LUN,27,1,OLDIMGNUM,    IRTFLG)
      IMGNUM = IFIX(OLDIMGNUM)

C     RETRIEVE OLD IREC  AND LABREC FROM HEADER
      CALL LUNGETLAB(LUN,LABREC,NDUM,IRECOLD,LABBYT,LENBYT,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     FIND CORRECT IREC (NUMBER OF RECORDS IN IMAGE + HEADER)
      IREC = NY + LABREC

C     SET CORRECT IREC IN STACKED IMAGE HEADER LOCATION: 3
      CALL LUNSETVALS(LUN,3, 1,FLOAT(IREC),     IRTFLG)

C     SET CORRECT ISTACK IN STACKED IMAGE HEADER LOCATION: 24
      CALL LUNSETVALS(LUN,24,1, 0.0,IRTFLG)

C     SET UNUSED VALUE IN STACKED IMAGE HEADER LOCATION: 25
      CALL LUNSETVALS(LUN,25,1, 0.0 ,IRTFLG)

      IF (VERBOSE .AND. (MYPID <= 0) ) THEN
C         ECHO HEADER VALUE CHANGE

         WRITE(NOUT,'(/,A,I0,A,I0,A,I0,A,I0,/)')
     &              '  Image: ',        IMGNUM,
     &              '   Stack offset:', LUNSTK(LUN),
     &              '   Header irec: ', IRECOLD,
     &              ' --> ',            IREC        

      ENDIF 
   
      IRTFLG = 0
      END




