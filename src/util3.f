 
C++*********************************************************************
C 
C  UTIL3.F        LONG FILENAMES           ArDean Leith           10/88  
C                                          MAHIEDDINE LADJADJ   4/23/93 
C                                          JING SU              8/31/93 
C                 REWRITTEN                ArDean Leith         1/15/98 
C                 ROT CALL CHANGED         ArDean Leith         8/02/00  
C                 ADDED 'ER'               ArDean Leith         2/16/01  
C                 ADDED 'ER SK'            ArDean Leith         4/23/01 
C                 ADDED 'ER EDM'           ArDean Leith         5/16/01 
C                 ADDED 'RT 3DQ'           ArDean Leith         4/10/01 
C                 ADDED 'ER WA '           ArDean Leith         4/25/02  
C                 OPFILEC                  ArDean Leith         2/18/03  
C                 'RT 3D' LUN closed       ArDean Leith         4/30/03 
C                 MRQLI1 -> APMASTER       ArDean Leith         9/ 5/03 
C                 MPI                      CHAO YANG           10/30/03 
C                 REMOVED RCONV            ArDean Leith        11/24/03 
C                 'RT QS' SELECTED FILES   ArDean Leith        12/15/06 
C                 'RTD **'                 ArDean Leith         1/15/07 
C                 ORACFMSKM CALL           ArDean Leith         3/19/08
C                 'OR S' CALL              ArDean Leith         6/06/08
C                 'OR R' CALL              ArDean Leith         6/06/08
C                 SETPRMB PARAMETERS       ArDean Leith         5/19/09
C                 LUNXM                    ArDean Leith        12/16/10
C                 ROTQSS REPLACED          ArDean Leith         1/12/11
C                 ROTQSS RTKSQ             ArDean Leith         5/24/11
C                 ROTQSS MIRROR            ArDean Leith         3/06/12
C                 TF FIND                  ArDean Leith         5/07/12
C                 OR MAP                   ArDean Leith        12/12/12
C                 OBSOLETE OR'S REMOVED    ArDean Leith        13/08/13
C                 ROT32 IOFF               ArDean Leith         3/12/13
C                 'AF' CALL                ArDean Leith         1/10/20
C
C ********************************************************************** 
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2020  Health Research Inc.,                         *
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
C  UTIL3(MAXDIM) 
C 
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12 
C--********************************************************************* 
 
        SUBROUTINE UTIL3(MAXDIM) 
 
	INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 
 
        INTEGER                :: MAXDIM

        INTEGER                :: ICOMM,MYPID,MPIERR
 
        CHARACTER(LEN=MAXNAM)  :: FILNAM,FILNM1,FILNAMO 
        CHARACTER (LEN=1)      :: CDUM 
	LOGICAL                :: MIRROR,CHKMIRROR 
        LOGICAL                :: ASKPEAKS
        CHARACTER (LEN=1)      :: NULL = CHAR(0)
 
	INTEGER,PARAMETER      :: LUN1   = 7 
	INTEGER,PARAMETER      :: LUN2   = 8 
	INTEGER,PARAMETER      :: LUN3   = 9 
	INTEGER,PARAMETER      :: LUNDOC = 87
	INTEGER,PARAMETER      :: LUNXM1 = 88
	INTEGER,PARAMETER      :: LUNXM2 = 89

        
C       DATA FUNC/'ED', 'RC', 'RT', 'BC', 'CT', 
c    &            'OR', 'FC', 'SL', 'RO', 'OD', 
c    &            'MK', 'AF', 'OP', 'DI', 'ER',
c    &            '13'/ 
C         13 is for 'RTD' operation


      CALL SET_MPI(ICOMM,MYPID,MPIERR)  ! SET MYPID

      IRTFLG = 0 
      MAXIM1 = 0 
      MAXIM2 = 0 

      SELECT CASE(FCHAR(1:2))
 
 
      CASE ('ED') !  -------  EDGE ------------------------------ 'ED' 

 	CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &              NX,NY,NZ, 
     & 		   MAXIM1,'INPUT',.FALSE.,IRTFLG) 
	IF (IRTFLG .NE. 0) GOTO 9000 
 
        IF (NZ .GT. 1) THEN 
C          DOES NOT WORK ON 3-D FILES 
           CALL ERRT(101,'SORRY DOES NOT WORK ON 3-D FILES',NE) 
           GOTO 9000 
        ENDIF 
 
C       OPEN AN OUTPUT FILE WITH DIMENSIONS BASED ON FIRST INPUT FILE 
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &             NX,NY,NZ, 
     &  	   MAXIM2,'OUTPUT',.TRUE.,IRTFLG) 
	IF (IRTFLG .NE. 0) GOTO 9000 
 
        CALL EDGE(LUN1,LUN2,LUN3,NX,NY) 
        GOTO 9000 
 
      CASE ('RC') !  ------- REAL CONVOLUTION ------------------ 'RC' 
 
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NX,NY,NZ, 
     &               MAXIM1,'INPUT',.FALSE.,IRTFLG) 
        IF (IRTFLG .NE. 0) GOTO 9000 
 
C       OPEN OUTPUT FILE WITH SAME DIMENSIONS AS INPUT FILE 
        CALL OPFILEC(LUN1,.TRUE.,FILNAMO,LUN2,'U',ITYPE,
     &               NX,NY, 
     &               NZ,MAXIM2,'OUTPUT',.TRUE.,IRTFLG) 
        IF (IRTFLG .NE. 0) GOTO 9000 
 
        CALL RCONV(LUN1,LUN2,LUN3,NX,NY,NZ,1,MAXDIM) 
        GOTO 9000 
 
      CASE ('13') !  ------- ROTATE & SCALE ------------------ 'RTD S*' 

        CHKMIRROR = (FCHAR(7:7) == 'M')

        IF (FCHAR(4:5) == 'SQ') THEN
C          OPERATION -- ROTATE & SCALE  ---------------------- 'RTD SQ' 
           CALL ROTQSS(LUN1,LUN2,LUNDOC,LUNXM1,LUNXM2,
     &                 .FALSE.,.FALSE.,CHKMIRROR,IRTFLG)

        ELSEIF ( FCHAR(4:5) == 'SF') THEN                   

C          OPERATION -- ROTATE & SCALE  ---------------------- 'RTD SF' 
           CALL ROTQSS(LUN1,LUN2,LUNDOC,LUNXM1,LUNXM2,
     &                 .FALSE.,.TRUE.,CHKMIRROR,IRTFLG)
        ENDIF

      CASE ('14') !  ------- ROTATE  ------------------------- 'ROT *' 
           CALL ROTATES()


      CASE ('RT') !  ------- ROTATE -------------------------- 'RT **' 
                    
        CHKMIRROR = (FCHAR(7:7) == 'M')

        IF (FCHAR(4:5) == 'SQ') THEN
C          OPERATION -- ROTATE & SCALE  ---------------------- 'RT SQ' 
           CALL ROTQSS(LUN1,LUN2,LUNDOC,LUNXM1,LUNXM2,
     &                   .TRUE.,.FALSE.,CHKMIRROR,IRTFLG)
           GOTO 9000 

        ELSEIF (FCHAR(4:5) == 'SF') THEN
C          OPERATION -- ROTATE & SCALE  ---------------------- 'RT SF' 
           CALL ROTQSS(LUN1,LUN2,LUNDOC,LUNXM1,LUNXM2,
     &                   .TRUE.,.TRUE.,CHKMIRROR,IRTFLG)
           GOTO 9000 
        ENDIF

        IF (FCHAR(4:5) == '3 ')  THEN 
            CALL ERRT(101," OBSOLETE OPERATION, USE 'ROT A' NOW",NE)
            GOTO 9000  
        ENDIF 

C       OPEN INPUT FILE 
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &                  NX,NY,NZ, 
     & 		        MAXIM1,'INPUT',.FALSE.,IRTFLG) 
        IF (IRTFLG .NE. 0) GOTO 9000 
 
        IF (FCHAR(4:4) == '3' .AND. ITYPE .NE. 3)  THEN 
            CALL ERRT(101,'NOT A VOLUME',NE) 
            GOTO 9000  
        ENDIF 
 
C       RECORD INPUT IMAGE  AVERAGE IN CASE NEEDED FOR BACKGROUND 
        AV1    = AV 
        FMAX1  = FMAX 
        FMIN1  = FMIN 
        IMAMI1 = IMAMI 
 
C       OPERATION -------- ROTATE 90, 180, 270 ---------------- 'RT 90' 
C       SPECIAL ARRANGEMENT FOR 'RT 90' (DO NOT KNOW DIMENSIONS YET) 
        IF (FCHAR(4:5) == '90') THEN 
C          90, 180, OR 270 DEGREE ROTATION OF VOLUME 
           CALL REFORM0(LUN1,LUN2,NX,NY,NZ,MAXDIM,IRTFLG) 
           GOTO 9000 
        ENDIF 
                      
C       OPEN OUTPUT FILE WITH SAME DIMENSIONS AS INPUT FILE 
        MAXIM2 = 0 
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &              NX,NY,NZ, 
     &              MAXIM2,'OUTPUT',.FALSE.,IRTFLG) 
        IF (IRTFLG .NE. 0) RETURN 
 
        IF (FCHAR(4:5) == '3D' .OR. FCHAR(4:5) == '3A')THEN 
C          OPERATION -------- ROTATE 3D 3A -------------------- 'RT 3??' 
           WRITE(NOUT,*) " OBSOLETE OPERATION, USE 'ROT C' NOW"
           CALL ROTAS3(LUN1,LUN2,NX,NY,NZ,FCHAR(4:6)) 
           GOTO 9000 
 
	ELSEIF (FCHAR(4:5) == '3L') THEN 
C          OPERATION -------- ROTATE 3D AROUND LINE ------------ 'RT 3L' 
           WRITE(NOUT,*) " OBSOLETE OPERATION, USE 'ROT L' NOW"
           CALL ROTAL3(LUN1,LUN2,NX,NY,NZ,FCHAR(4:6)) 
           GOTO 9000 
        ENDIF 
 
C       OPERATION -------- ROTATE ------------------------------ 'RT' 
C       OPERATION -------- ROTATE USE MIN AS BACK. ------------- 'RT M' 
C       OPERATION -------- ROTATE SPECIFY BACK ----------------- 'RT B' 
C       OPERATION -------- ROTATE AROUND ARBITRARY CENTER  ----- 'RT C' 
 
        IF (FCHAR(4:4) .NE. 'C') THEN 
           WRITE(NOUT,*) " OBSOLETE OPERATION, USE 'ROT' NOW"
        ENDIF

        CALL RDPRM(THETA,NOT_USED,'ROTATION ANGLE') 

        TH = THETA * DATAN(1.0D0) / 45.0D0 

        IF (FCHAR(4:4) .NE. 'B') THEN 
C          BACKGROUND IS INPUT IMAGE AVERAGE, 
C          EXCEPT FOR "RT M" WHERE BACKGROUND IS MINIMUM 
           IF (IMAMI1 .NE. 1) 
     &        CALL NORM3(LUN1,NX,NY,NZ,FMAX1,FMIN1,AV1) 
           IF (FCHAR(4:4) == 'M') AV1 = FMIN1 
        ELSE 
C          FOR "RT B" MUST SUPPLY BACKGROUND 
           CALL RDPRM(AV1,NOT_USED,'BACKGROUND') 
        ENDIF 
 
C       SET ROTATION CENTER 
        SHX = 0.0 
        SHY = 0.0 

        IF (FCHAR(4:4) == 'C') THEN 
C          ROTATE AROUND AN ARBITRARY CENTER  
           WRITE(NOUT,*) " OBSOLETE OPERATION, USE 'ROT C' NOW"
           IF (NX*(1+NY) > MAXDIM) THEN 
              CALL ERRT (101, 
     &           'VARIABLE CENTER NOT AVAILABLE FOR THIS SIZE IMAGE',NE) 
              GOTO 9000 
           ENDIF 

           SHX = 0
           CALL RDPRM1S(SHX,NOT_USED,'X-SHIFT',IRTFLG) 
           SHY = 0
           CALL RDPRM1S(SHY,NOT_USED,'Y-SHIFT',IRTFLG) 
        ENDIF 
 
C       CAN DO IMAGE OR VOLUME 
        DO  ISLICE = 1, NZ 

C           ROTATE SLICE IN-CORE 
            NYS  = (ISLICE-1) * NY + 1 
            NYE  = NYS + NY - 1 
            IOFF = NYS - 1
            CALL ROT32(LUN1,LUN2,NX,NYS,NYE,1, TH,AV1,SHX,SHY,IOFF) 
        ENDDO 
 
             
      CASE ('BC') !  ------- BOX CONVOLUTION --------------------- 'BC' 
 
 	CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &             NX1,NY1,NZ1, 
     & 		   MAXIM1,'INPUT',.FALSE.,IRTFLG) 
	IF (IRTFLG .NE. 0) GOTO 9000 
 
        IMAMI1 = IMAMI 
        AV1    = AV 
 
C       OPEN AN OUTPUT FILE WITH DIMENSIONS SAME AS INPUT FILE 
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE, 
     &               NX1,NY1,NZ1, MAXIM2,'OUTPUT',.TRUE.,IRTFLG) 
	IF (IRTFLG .NE. 0) GOTO 9000 
 
C       BOX CONVOLUTION 
        IMAMI = IMAMI1 
        AV    = AV1 
        CALL BOXX(LUN1,LUN2,NX1,NY1,NZ1,MAXDIM) 
 
C       SET HEADER FOR ALTERATIONS IN IMAGE DUE TO OPERATION 
        CALL  SETPRMB(LUN2, 0.0,0.0, 0.0,0.0) 
                                 
      CASE ('CT') !  -------- CT --------------------------------'CT' 

        CALL ERRT (101, 'OBSOLETE OPERATION REMOVED',NE) 
        RETURN 
 
      CASE ('OR') !  -------- ORIENTATIONAL SEARCH ------------- 'OR' 

       IF (FCHAR(4:4) == '2') THEN 
C          OPERATION -------- ORIENTATIONAL SEARCH ------------- 'OR 2' 
C          OPERATION -------- ORIENTATIONAL SEARCH ------------- 'OR 2M' 
           IF (mypid <= 0) WRITE(NOUT,888) 
888        FORMAT(" OBSOLETE OPERATION, USE: 'OR R' INSTEAD",/)
           RETURN 

        ELSEIF (FCHAR(4:4) == 'R') THEN 
C          OPERATION -------- ORIENTATIONAL SEARCH -------------- 'OR R' 
C          REPLACES 'OR 2' & 'OR 2M' 
           CALL ORMD(.TRUE.) 
           RETURN 
 
        ELSEIF (FCHAR(4:4) == 'A')  THEN 
C          OPERATION -------------------------------------------- 'OR A' 
	   IF (FCHAR(5:5) == 'M')  THEN 
              CALL ERRT(101,"OPERATION IS NOW: 'OR A'",NE)
              RETURN
           ENDIF
           CALL ORACFMSK() 
           RETURN 

        ELSEIF (FCHAR(4:4) == 'Q')  THEN 
C          OPERATION -------------------------------------------- 'OR Q' 
	   IF (FCHAR(5:5) == 'M')  THEN 
              CALL ERRT(101,"OPERATION IS NOW: 'OR Q'",NE)
              RETURN
           ENDIF
           CALL ORACFMSKM() 
           RETURN 
 
        ELSEIF (FCHAR(4:4) == 'S')  THEN 
C          OPERATION REPLACES 'MQ' & 'NQ' ---------------------- 'OR SH' 
           CALL APMASTER('F','ORS') 
           RETURN 
 
        ELSEIF (FCHAR(4:5) == 'MA')  THEN 
C          OPERATION ------------------------------------------ 'OR MAP' 
           CALL ORMAP() 
           RETURN 

        ELSEIF (FCHAR(4:5) == 'MQ')  THEN 
C          OPERATION ------------------------------------------- 'OR MQ' 
           IF (mypid <= 0) WRITE(NOUT,887) 
887        FORMAT(" OBSOLETE OPERATION, USE: 'OR SH' INSTEAD",/)
           RETURN 
 
        ELSEIF (FCHAR(4:5) == 'NQ')  THEN 
C          OPERATION  ------------------------------------------ 'OR NQ' 
           IF (mypid <= 0) WRITE(NOUT,887) 
           RETURN 
 
        ELSEIF (FCHAR(4:5) == '3Q')  THEN 
C          OPERATION ------------------------------------------- 'OR 3Q' 
           CALL QALI('Q') 
           GOTO 9911 
 
        ELSEIF (FCHAR(4:5) == '3A')  THEN 
C          OPERATION ------------------------------------------- 'OR 3A' 
           CALL QALI('A') 
           GOTO 9911 
        END IF 
        GOTO 9911 

      CASE ('FC') !  ------ ------------------------------------- 'FC' 
 
	CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,NX1,
     &               NY1,NZ1, 
     &               MAXIM1,'INPUT',.FALSE.,IRTFLG) 
	IF (IRTFLG .NE. 0) GOTO 9000 
 
        IF (NZ1 .NE. 1) THEN 
C          DOES NOT WORK ON 3-D FILES 
           CALL ERRT(101,'DOES NOT WORK ON VOLUMES',NE) 
           GOTO 9000 
        END IF 
 
        FMAX1  = FMAX 
        FMIN1  = FMIN 
        IMAMI1 = IMAMI 
 
C       OPEN AN OUTPUT FILE WITH SAME DIMENSIONS AS INPUT FILE 
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE, 
     &               NX1,NY1,NZ1, 
     &               MAXIM2,'OUTPUT',.TRUE.,IRTFLG) 
	IF (IRTFLG .NE. 0) GOTO 9000 
 
C       FILE CONTOUR 
        FMAX  = FMAX1 
        FMIN  = FMIN1 
        IMAMI = IMAMI1 
        CALL CNTRFL(LUN1,LUN2,NX1,NY1,NZ1,MAXDIM) 
 
C       SOME OPTIONS REQUIRE RESCALING (AL 28 JAN 92) 
        CALL SETPRMB(LUN2, 0.0,0.0, 0.0,0.0) 
 
      CASE ('SL') !  ------ -------------  SLICE ------------------'SL' 
 
        CALL SLICE(MAXDIM,LUN1,LUN2,LUN3) 
      
      CASE ('RO') !  ------ ROTATIONAL AVERAGE --------------------'RO' 
  
        IF (FCHAR(4:4) == 'I') THEN 
	   CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE, 
     &                  NX1,NY1,NZ1, 
     &                  MAXIM1,'INPUT',.FALSE.,IRTFLG) 
	   IF (IRTFLG .NE. 0) GOTO 9000 
 
           IF (NZ1 .NE. 1) THEN 
C             OPERATION DOES NOT WORK ON 3-D FILES 
              CALL ERRT(101,'DOES NOT WORK ON VOLUMES',NE) 
              GOTO 9000 
	   ENDIF 
 
C          OPEN AN OUTPUT FILE  
	   IFORM = 1 
           CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE, 
     &                  NX1,NY1,NZ1,
     &                  MAXIM2,'OUTPUT',.TRUE.,IRTFLG) 
	   IF (IRTFLG .NE. 0) GOTO 9000 
 
C          ROTATIONAL AVERAGING 
           CALL RADAVI(LUN1,LUN2,NX1,NY1,MAXDIM) 
        ELSE 
           CALL RADAV(LUN1,LUN2) 
        ENDIF 
 
      CASE ('OD') !  ------  OPTICAL DENSITY -------------------- 'OD' 
 
	CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NX1,NY1,NZ1,
     &               MAXIM1,'INPUT',.FALSE.,IRTFLG) 
	IF (IRTFLG .NE. 0) GOTO 9000 
 
        IF (NZ1 .NE. 1) THEN 
C          DOES NOT WORK ON 3-D FILES 
           CALL ERRT(101,'DOES NOT WORK ON VOLUMES',NE) 
           GOTO 9000 
        END IF 
 
C       OPEN AN OUTPUT FILE WITH DIMENSIONS SAME AS INPUT FILE 
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE, 
     &               NX1,NY1,NZ1, 
     &               MAXIM2,'OUTPUT',.TRUE.,IRTFLG) 
	IF (IRTFLG .NE. 0) GOTO 9000 
 
        CALL OD(LUN1,LUN2,LUN3,NX1,NY1,MAXDIM) 
 

      CASE ('MK') !  ------------ MARKER ------------------------ 'MK' 

        IF (FCHAR(4:5) == '3')  THEN      !-------------------- 'MK 3'
           CALL MRK3(MAXDIM) 
        ELSEIF (FCHAR(4:6) == 'RT') THEN  !-------------------- 'MK RT'
           CALL MRRT(MAXDIM) 	
        ELSE 
	   CALL MRK(MAXDIM)                 !-------------------- 'MK'
	ENDIF 
 

      CASE ('AF') !  ------- TRANSFORMATION --------------------- 'AF' 
C       UNDOCUMENTED OPERATION !!
        MAXIM1 = -3
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE, 
     &               NX1,NY1,NZ1, MAXIM1,'INPUT',
     &              .FALSE.,IRTFLG) 
        IF (IRTFLG .NE. 0) GOTO 9000 
 
        IF (MAXIM1  >=  0) THEN   ! FEB 2020 al
           CALL ERRT(101,
     &        'THIS OPERATION NO LONGER ACCEPTS BARE STACK INPUT',NDUM)
           GOTO 9000
        ENDIF
  
C       OPEN AN OUTPUT FILE WITH SAME DIMENSIONS AS INPUT FILE 
        CALL OPFILEC(LUN1,.TRUE.,FILNAMO,LUN2,'U',ITYPE, 
     &               NX1,NY1,NZ1,MAXIM2,
     &               'OUTPUT',.TRUE.,IRTFLG) 
        IF (IRTFLG .NE. 0) GOTO 9000 
 
	CALL AF(MAXDIM, LUN1,LUN2, NX1,NY1,NZ1, IRTFLG) 
 
      CASE ('OP') !  ----- ORIENTATION OF PROJECTIONS ------------ 'OP' 
  	CALL POLQS(MAXDIM) 
 
      CASE ('DI') !  ----- DILATION ------------------------------ 'DI' 
 	CALL  OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE, 
     &        NX,NY,NZ, MAXIM1,'INPUT',.FALSE.,IRTFLG) 
        IF (IRTFLG .NE. 0) RETURN 
 
        CALL  OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE, 
     &                NX,NY,NZ,
     &                MAXIM2,'OUTPUT',.FALSE.,IRTFLG) 
        IF (IRTFLG .NE. 0) GOTO 9000 
 
	CALL  DILATION(LUN1,LUN2,NX,NY,NZ) 
	GOTO 9000 
 
      CASE ('ER') !  -----  EROSION ----------------------------- 'ER' 
 	CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &       NX,NY,NZ, 
     &       MAXIM1,'INPUT',.FALSE.,IRTFLG) 
        IF (IRTFLG .NE. 0) RETURN 
        FMIN1 = FMIN 
        FMAX1 = FMAX 
 
        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &              NX,NY,NZ, 
     &              MAXIM2,'OUTPUT',.FALSE.,IRTFLG) 
        IF (IRTFLG .NE. 0) GOTO 9000 
 
        IF (FCHAR(4:5) == 'SK') THEN 
C          BINARY SKELETON 
	   CALL SKELETON(LUN1,LUN2,NX,NY,NZ) 
 
        ELSEIF (FCHAR(4:4) == 'E') THEN 
C          EUCLIDEAN DISTANCE MAP 
	   CALL EDM(LUN1,LUN2,NX,NY,NZ,FMIN1,FMAX1) 
 
        ELSEIF (FCHAR(4:4) == 'W') THEN 
C          WATERSHED  MAP 
	   CALL WATERSHD(LUN1,LUN2,NX,NY,NZ) 
 
        ELSE 
	   CALL EROSION(LUN1,LUN2,NX,NY,NZ) 
        ENDIF 
        
      END SELECT
 

9000  CLOSE(LUN1) 
      CLOSE(LUN2) 
      CLOSE(LUN3) 
 
9911  RETURN 
      END 
 









