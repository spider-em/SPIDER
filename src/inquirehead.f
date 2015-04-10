C++*********************************************************************
C                                                                      *
C  INQUIREHEAD.F     NEW                         MAY 2009 ARDEAN LEITH *
C                    NORM CALL                   NOV 2010 ARDEAN LEITH *                          
C                    PIXSIZ LOC.                 NOV 2010 ARDEAN LEITH * 
C                    1PG FORMAT                  NOV 2010 ARDEAN LEITH *                       
C                    NX.., PROJ...               JUN 2011 ARDEAN LEITH *                       
C                    GLONUM                      OCT 2013 ARDEAN LEITH *                       
C                                                                      *
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
C   INQUIREHEAD(LUN,NX,NY,NZ.IRTFLG)    
C
C   PURPOSE: RETRIEVES VALUES OF HEADER VARIABLE(S) INTO REGISTERS
C
C   PARAMETERS:  LUN       I/O UNIT                           (SENT)
C                NX..      FILE SIZE                          (SENT)
C                IRTFLG    ERROR FLAG                         (RET.)
C                          0  IS NORMAL
C                          1 INQUIRY WAS NOT AS EXPECTED
C
C--*********************************************************************

      SUBROUTINE INQUIREHEAD(LUN,NX,NY,NZ,IRTFLG)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC' 
      INCLUDE 'CMLIMIT.INC' 
 
      INTEGER                   :: LUN,NX,NY,NZ,IRTFLG

      INTEGER, PARAMETER        :: MAXREGNAM = 10
      CHARACTER(LEN=MAXREGNAM)  :: REGNAME

C     MAXNSEL IS CURRENTLY SAME AS IN REG_SET.F !!!
      INTEGER, PARAMETER        :: MAXNSEL = 24  ! SEARCH & REGISTER LIST
      INTEGER                   :: LOCATION(MAXNSEL)
      INTEGER                   :: IREGSELS(MAXNSEL)
      REAL                      :: VALUES(MAXNSEL)
      CHARACTER(LEN=8)          :: NAMES(MAXNSEL)

      INTEGER                   :: ICOMM,MYPID,IVALS,I,NREG,IHEDLOC
      INTEGER                   :: UNUSED,LEN,LENT
      INTEGER                   :: lnblnkn

      CALL SET_MPI(ICOMM,MYPID,IRTFLG)   ! RETURNS MPI PROCESS ID OR -1

      CALL INQUIREHEAD_LOC(LOCATION,NAMES,IVALS,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      DO I=1,IVALS     ! LOOP OVER ALL  WANTED HEADER VALUES
         !write(6,*) 'locs:', location(I),names(i)
         IF (LOCATION(I) >= 7 .AND. LOCATION(I) <= 10) THEN
            IF (IMAMI .NE. 1) THEN
C              STATISTICS NOT DETERMINED YET!
               CALL NORM3(LUN,NX,NY,NZ, FMAX,FMIN,AV)
            ENDIF
            EXIT
         ENDIF
      ENDDO

C     RETRIEVE REGISTER NUMBER(S) (IF ANY)  FROM OPERATION LISTING
      CALL REG_GET_SELS(IREGSELS,MAXNSEL,NREG,IRTFLG)

      DO I=1,IVALS    ! LOOP OVER ALL  WANTED HEADER VALUES

        IHEDLOC = LOCATION(I)     ! HEADER LOCATION WANTED

C       GET HEADER VALUE FROM  THIS HEADER LOCATION  
        CALL GETLAB(LUN,NX,UNUSED,IHEDLOC,1,VALUES(I),IRTFLG)

        !write(6,*) ' IHEDLOC: ',IHEDLOC,names(i), VALUES(I)

C       KLUDGE FOR NZ < 0
        IF (LOCATION(I) == 1 .AND. VALUES(I) .LT. 0) 
     &      VALUES(I) = ABS(VALUES(I))

        IF (VERBOSE .AND. (MYPID <= 0) ) THEN
C          ECHO VALUE

           IF (NREG <= 0) THEN
C             NOT SETTING ANY REGISTER, ECHO VALUE ONLY
              !WRITE(NOUT,90) NAMES(I)(1:LENT),VALUES(I),LENT
              WRITE(NOUT,90) NAMES(I),VALUES(I)
90            FORMAT('  ',A,' = ',1PG11.3)
           ELSE
C             GET REGISTER NAME
              CALL REG_GET_NAME(IREGSELS(I),REGNAME,LEN,IRTFLG)

C             ECHO HEADER VALUE & REGISTER SETTING
              WRITE(NOUT,91) REGNAME(1:LEN),NAMES(I),VALUES(I)
91            FORMAT('  REGISTER VARIABLE: ',A,'  HOLDS: ',A,
     &                  ' = ',1PG11.3)
           ENDIF
        ENDIF  ! END OF: IF (VERBOSE .AND. (MYPID <= 0) ) THEN
      ENDDO    ! END OF: DO I=1,IVALS  
  
C     SET REGISTER VARIABLES TO HEADER VALUES
      CALL REG_SET_NSELA(IVALS,VALUES,.FALSE.,IRTFLG)
 
      END

C     -------------------- INQUIREHEAD_LOC ----------------------------


      SUBROUTINE INQUIREHEAD_LOC(LOCATION,NAMEGOT,IVALS,IRTFLG)

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC' 
      INCLUDE 'CMLIMIT.INC' 
 
C     MAXNSEL IS CURRENTLY SAME AS IN REG_SET.F !!!
      INTEGER, PARAMETER        :: MAXNSEL = 24  ! SEARCH & REGISTER LIST
      INTEGER                   :: LOCATION(MAXNSEL)
      CHARACTER(LEN=8)          :: NAMEGOT(MAXNSEL)
      INTEGER                   :: IVALS
      INTEGER                   :: IRTFLG

      LOGICAL                   :: STRIP
      CHARACTER(LEN=80)         :: RESPONSE,ARGNOW,MSG
      CHARACTER(LEN=1)          :: NULL = CHAR(0)
      INTEGER                   :: ICOMM,MYPID,NT,ILOC,LENVAR,LENT,NE
      INTEGER                   :: I,IFIRST, IGO, IEND, NLETA, IVAR
      INTEGER                   :: lnblnk

      INTEGER, PARAMETER        :: MAXHEDNAM = 50
      CHARACTER(LEN=8)          :: NAME(MAXHEDNAM)

      CALL SET_MPI(ICOMM,MYPID,IRTFLG)   ! RETURNS MPI PROCESS ID OR -1

      DO I=1,MAXHEDNAM
        NAME(I) = ''
      ENDDO

      NAME(1)  =  'NZ'
      NAME(2)  =  'NY'
      NAME(3)  =  'IREC'
      NAME(4)  =  'UNUSED'
      NAME(5)  =  'IFORM'
      NAME(6)  =  'IMAMI'
      NAME(7)  =  'FMAX'
      NAME(8)  =  'FMIN'
      NAME(9)  =  'AV'
      NAME(10) =  'SIG'
      NAME(11) =  'UNUSED'
      NAME(12) =  'NX'
      NAME(13) =  'LABREC'
      NAME(14) =  'IANGLE'
      NAME(15) =  'PHI'
      NAME(16) =  'THETA'
      NAME(17) =  'PSI'
      NAME(18) =  'XOFF'
      NAME(19) =  'YOFF'
      NAME(20) =  'ZOFF'
      NAME(21) =  'SCALE'
      NAME(22) =  'LABBYT'
      NAME(23) =  'LENBYT'
      NAME(24) =  'ISTACK'
      NAME(25) =  'UNUSED'
      NAME(26) =  'MAXIM'
      NAME(27) =  'IMGNUM'
      NAME(28) =  'LASTINDX'
      NAME(29) =  'UNUSED'
      NAME(30) =  'UNUSED'
      NAME(31) =  'KANGLE'
      NAME(32) =  'PHI1'
      NAME(33) =  'THETA1'
      NAME(34) =  'PSI1'
      NAME(35) =  'PHI2'
      NAME(36) =  'THETA2'
      NAME(37) =  'PSI2'
      NAME(38) =  'PIXSIZ'
      NAME(39) =  'EV'
      NAME(40) =  'PROJ'
      NAME(41) =  'MIC'
      NAME(42) =  'NUM'
      NAME(43) =  'GLONUM'
      NAME(44) =  'UNUSED'
      NAME(45) =  'UNUSED'
      NAME(46) =  'UNUSED'
      NAME(47) =  'UNUSED'
      NAME(48) =  'UNUSED'
      NAME(49) =  'UNUSED'
      NAME(50) =  'UNUSED'

      STRIP     = .TRUE.
      CALL RDPRMC(RESPONSE,NLETA,STRIP,'HEADER VARIABLE(S)',
     &            NULL,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IFIRST = 1
      IVALS  = 0
      DO WHILE (IFIRST .LT. NLETA) 
         CALL GETNEXTTOKEN(RESPONSE,IFIRST,IGO,IEND)
         IF (IGO <= 0) EXIT
         ARGNOW = RESPONSE(IGO:IEND)
         IFIRST = IEND + 1         ! START OF NEXT TOKEN
         NT     = IEND - IGO + 1   ! CHAR. IN ARGNOW

         ! write(6,*) ' token: ',argnow(1:nt)

         IVAR = 0

         IF     (ARGNOW(1:NT) == 'MAX') THEN
C           SPECIAL ALLOWS 'MAX' FOR 'FMAX'
            IVAR = 7
         ELSEIF (ARGNOW(1:NT) == 'AVG') THEN
            IVAR = 9               ! ALTERNATE LABEL

         ELSEIF (ARGNOW(1:NT) == 'NSAM') THEN
            IVAR = 12              ! ALTERNATE LABEL

         ELSEIF (ARGNOW(1:NT) == 'NROW') THEN
            IVAR = 2               ! ALTERNATE LABEL

         ELSEIF (ARGNOW(1:NT) == 'NSLICE') THEN
            IVAR = 1               ! ALTERNATE LABEL

         ELSE
            DO I = 1,MAXHEDNAM
               ILOC = INDEX(NAME(I),ARGNOW(1:NT))
               IF (ILOC > 0) THEN
C                 POSSIBLE MATCH

                  LENVAR = lnblnk(NAME(I))
                  !write(6,*) ' matches: ',name(i),nt,lenvar

                  IF (NT == LENVAR) THEN
                     IVAR = I
                     EXIT       ! COMPLETE MATCH

                  ELSEIF (IVAR > 0) THEN
C                    DOUBLE MATCH
                     MSG = 'HEADER VARIABLE NAME: "'// ARGNOW(1:NT) // 
     &                     '" IS NOT UNIQUE, USE LONGER NAME'
                     LENT = lnblnk(MSG)
                     CALL ERRT(101,MSG(1:LENT),NE)
                     RETURN
                  ENDIF            
                  IVAR = I
               ENDIF
            ENDDO
         ENDIF

         IF (IVAR <= 0) THEN
            MSG = 'UNKNOWN HEADER VARIABLE: "' // ARGNOW(1:NT) //'"' 
            CALL ERRT(101,MSG,NE)
            IRTFLG = 1
            RETURN

         ELSEIF (IVAR > 0) THEN
C           GOT A MATCHING HEADER VARIABLE
            IVALS = IVALS + 1
            IF (IVALS > MAXNSEL) THEN
                CALL ERRT(102,
     &                 'TOO MANY VARIABLES REQUESTED, LIMIT',MAXNSEL)
                IRTFLG = 1
                RETURN
            ENDIF
            LOCATION(IVALS)      = IVAR
            NAMEGOT(IVALS)(1:8)  = NAME(IVAR)(1:8)
            !write(6,*) ' matched',ivar,namegot(ivals)
         ENDIF
      ENDDO  ! END OF: DO WHILE (IFIRST .LT. NLETA) 
C     NOW HAVE LISTED ALL MATCHING VARIABLE NUMBERS IN: LOCATION 

      ! write(6,*) ' matching values',ivals
 
      END



