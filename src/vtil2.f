
C++************************************************ 6/23/80 2/20/81 VAX
C
C VTIL2.F          FILENAMES LENGTHENED               1/89 ARDEAN LEITH
C                  CHANGED                            4/93 M. LADJADJ
C                  CHANGED                            8/93 JING SU
C                  PJ CYL NEEDED OPFILE               9/01 ARDEAN LEITH
C                  OPFILEC                         3/18/03 ARDEAN LEITH 
C                  IRTFLG = 0                     10/28/03 ARDEAN LEITH
C                  PJ  CASE                       04/13/05 ARDEAN LEITH
C                  'PJ RG' REMOVED                10/18/05 ARDEAN LEITH
C                  CASE                           12/20/06 ARDEAN LEITH
C                  'RB' ADDED                      1/02/07 ARDEAN LEITH
C                  'BPD' ADDED                     1/23/07 ARDEAN LEITH
C                  'PJ RG' REMOVED                10/18/05 ARDEAN LEITH
C                  'PJ 3G' ADDED                   3/28/07 ARDEAN LEITH
C                  'BP 3G' ADDED                   3/28/07 ARDEAN LEITH
C                  'BPD --> BP, BP --> OLD'        6/08/08 ARDEAN LEITH
C                  'BP NF --> BP 3N'               6/16/08 ARDEAN LEITH
C                  WIW3D_OLD                      10/17/08 ARDEAN LEITH
C                  WIW3* RENAMED BP*               1/10/11 ARDEAN LEITH
C                  REPCG RENAMED BPCG              1/12/11 ARDEAN LEITH
C                  REPS RENAMED BPRP               5/16/11 ARDEAN LEITH
C                  REMOVED 'PJ 3O' & 'PJ 3QO'     11/07/11 ARDEAN LEITH
C                  FBS SUPPORT IN PJ3Q            11/07/11 ARDEAN LEITH
C                  ADDED 'BP RP 3'                04/09/12 ARDEAN LEITH
C                  REMOVED 'BP 3F O'               1/04/13 ARDEAN LEITH
C                  'BP 3F P','BP 3F M'             2/14/13 ARDEAN LEITH
C                  'SK' PARAMETERS                 5/20/13 ARDEAN LEITH
C                  'DR ERR' REMOVED                2/10/14 ARDEAN LEITH
C                  'PJ 2D' ADDED                   5/24/18 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2018  Health Research Inc.,                         *
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
C  VTIL2(MAXDIM)
C
C  HANDLES: 'PS SK CS PJ BP DC DR MF RB BP '
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE VTIL2(MAXDIM)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER               :: MAXDIM

        LOGICAL               :: REFRINGS,REPLACE

        CHARACTER(LEN=MAXNAM) :: FILNAM,FIL

        INTEGER               :: MAXIM  
        INTEGER               :: MAXIM2 
        INTEGER               :: IRTFLG 

        INTEGER, PARAMETER    :: LUN1   = 8            
        INTEGER, PARAMETER    :: LUN2   = 9
        INTEGER, PARAMETER    :: LUN3   = 10
        INTEGER, PARAMETER    :: LUN4   = 11

        MAXIM  = 0
        MAXIM2 = 0
        IRTFLG = 0

        SELECT CASE(FCHAR(1:2))

        CASE ('RB') !  -------- ROTATE & BACK PROJECT  -----------   RB
           SELECT CASE(FCHAR(4:5))

           CASE ('32') !   
              CALL BP32F(.TRUE.)
           CASE ('3F') !   
              CALL BP3F(.TRUE.,.FALSE.)
           CASE DEFAULT
              CALL ERRT(101,'UNKNOWN OPERATION',NDUM)
           END SELECT

        CASE ('DC') !  -------------- DECIMATE -------------------   DC
           CALL  DECIMATE

        CASE ('DR') !  --------------- ERROR ---------------------   DR

           SELECT CASE(FCHAR(4:6))
       
           CASE ('ERR') !   
C             CALCULATE ERROR MEASURES BETWEEN 2 VOLUMES
              CALL ERRT(101,
     &            'OBSOLETE, USE <DR DIFF> INSTEAD',NDUM)
              !CALL COMP3D(LUN1,LUN2)

           CASE ('DIF') !   
C             CALCULATE ERROR MEASURES BETWEEN 2 VOLUMES WITHIN 
C             BOUNDARIES OF A MASK, SCALE VOLS AND CALCULATE 
C             DIFFERENCE VOL
              CALL COMP3DMAD(LUN1,LUN2,LUN3,LUN4)

           CASE DEFAULT
              CALL ERRT(101,'UNKNOWN/OBSOLETE OPERATION',NDUM)
           END SELECT


        CASE ('MF') !  --------------------------------------------  MF
C          FIT SPHERE MODEL TO A 3-D RECONSTRUCTION
           CALL ERRT(101,
     &         'OBSOLETE SUBROUTINE (LUNA OR MATVEC) CALLED',NE)

        CASE ('SK') !  -------------- STACK SLICES ----------------  SK

          !WRITE(NOUT,'(A,/)') '  OBSOLETE OPERATION, USE: <CP TO VOL>'
          IF (FCHAR(4:) == 'R') THEN
             CALL STACK_REPLACE()
          ELSE
             CALL STACK()
          ENDIF

        CASE ('CS') ! ---- ARBITRARY  SLICING (SAME AS "PS A") ----  CS
C         ARBITRARY DIRECTION OF SLICING (SAME AS "PS A")
          CALL CSLICE

        CASE ('PS') ! -----------------PICK SLICE -----------------  PS

          IF (FCHAR(4:4) == 'A') THEN
C           ARBITRARY DIRECTION OF SLICING
            CALL CSLICE

          ELSE
            CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &               NX1,NY1,NZ1, MAXIM,
     &               'INPUT',.FALSE.,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999
            IF (IFORM .NE. 3) THEN
               CALL ERRT(2,'VTIL2',NE)
               GOTO 9999
            ENDIF

            FMIN1 = FMIN
            FMAX1 = FMAX
            AV1   = AV
            SIG1  = SIG
            ITYPE = 1
            IF (FCHAR(4:4) == 'Z') THEN

C             WANT Z SLICE
              NZ = 1
              CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &                     NX1,NY1,NZ,MAXIM2,
     &                    'OUTPUT',.FALSE.,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999
              CALL PICKSL(LUN1,LUN2,NX1,NY1,NZ1)
              NX = NX1
            ELSE

C             WANT X OR Y SLICE
              NX = NX1
              IF (FCHAR(4:4) == 'X') NX = NY1
              NZ = 1
              !write(6,*) 'opening maxim2:',maxim2

              CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',ITYPE,
     &                   NX,NZ1,NZ,MAXIM2,
     &                   'OUTPUT',.FALSE.,IRTFLG)
              CALL PICKSV(LUN1,LUN2,NX1,NY1,NZ1)
            ENDIF
            IF (FCHAR(5:5) == 'N') THEN
C              KEEP FMIN AND FMAX SAME FOR ALL SLICES
               IF (MYPID <= 0) WRITE(NOUT,*) 
     &            ' SETTING FMIN & FMAX:',FMIN1,FMAX1
               SIG = SIG1
               CALL SETPRM(LUN2,NX,IDUM,FMAX1,FMIN1,AV1,'U')
            ENDIF
          ENDIF


        CASE ('PJ') !  -------------- PROJECTION -----------------   PJ
C         MOST "PJ" ROUTINES OPEN THEIR OWN FILES

          NCT = lnblnkn(FCHAR)
          SELECT CASE(FCHAR(4:NCT))

          CASE ('2')
              CALL PJ2D(LUN1,LUN2,LUN3)

          CASE ('3')
              IF (USE_FBS_INTERP) THEN 
                 CALL PJ3_FBS()
              ELSE
                 CALL PJ3_N()
              ENDIF

          CASE ('3O')
              CALL ERRT(101,'OBSOLETE OPERATION, USE <PJ 3>',NE)
 
          CASE ('3Q RR')  ! NOT FBS, MAKE REFRINGS 
              REFRINGS = .TRUE.
              CALL PJ3Q_N(.FALSE., REFRINGS)

          CASE ('3Q')  ! NOT FBS 
              REFRINGS = .FALSE.
              CALL PJ3Q_N(.FALSE., REFRINGS)

          CASE ('3F RR')  ! FBS, MAKE REFRINGS
              REFRINGS = .TRUE.
              CALL PJ3Q_N(.TRUE., REFRINGS)

          CASE ('3F')  ! FBS   
               REFRINGS = .FALSE.
             CALL PJ3Q_N(.TRUE., REFRINGS)

          CASE ('3G')
              CALL PJ3G()     ! GRIDDED PROJECTION? MAR 07

          CASE ('3Q O','3QO')
              CALL ERRT(101,'OBSOLETE OPERATION, USE <PJ 3Q>',NE)

          CASE ('ST')
              CALL MRRSURF
 
          CASE ('SU')
              CALL MRSURF
 
          CASE ('SHAD')
              CALL MRREFL
 
          CASE ('COL')
              CALL MRNCOLOR

          CASE ('A','AT')
C            PROJECT VOLUME USING EXPONENTIAL ATTENUATION

C            OPEN FIRST INPUT FILE 
             CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,
     &               NX1,NY1,NZ1,
     &             MAXIM,'INPUT',.FALSE.,IRTFLG)
             IF (IRTFLG .NE. 0) GOTO 9999

             CALL MRPRREP(LUN1,LUN2,MAXDIM,IER)

          CASE ('C','CY','CYL')

C            CYLINDRICAL PROJECTION
             CALL MRCP(LUN1,LUN2,LUN3)

          CASE DEFAULT
              CALL ERRT(101,'UNKNOWN/OBSOLETE OPERATION',NDUM)

          END SELECT
          CALL FLUSHRESULTS  


        CASE ('BP') !  -------- BACK PROJECTION ------------------   BP

C          ALL 'BP' ROUTINES OPEN THEIR OWN FILES

           IF (FCHAR(7:7) == 'O' .OR. FCHAR(8:8) == 'O' ) THEN !  'BP 3F O'
              CALL ERRT(101,'OBSOLETE OPERATION REMOVED',NDUM)
              GOTO 9999
           ENDIF

           SELECT CASE(FCHAR(4:5))
       
           CASE ('W2')                      !  'BP W2'  
              CALL WGBP2(MAXDIM)

           CASE ('RP')                      !  'BP RP' 
              IF (FCHAR(7:7) == '3') THEN
                 CALL BPRP3()               !  'BP RP 3'
              ELSE
                 CALL BPRP()                !  'BP RP'
              ENDIF

           CASE ('R2')                      !  'BP R2'  
              CALL BPWR(MAXDIM)

           CASE ('S2')                      !  'BP S2'   
              CALL BPS2(MAXDIM)

           CASE ('3D')                      !  'BP 3D'   
              CALL BCQ(MAXDIM)

           CASE ('3F')                      !  'BP 3F'  

              IF (FCHAR(7:7) == 'M') THEN
                 CALL BP3F_MERGE()          !  'BP 3F M'

              ELSEIF (FCHAR(7:7) == 'P' ) THEN
                 CALL BP3F(.FALSE.,.TRUE.)  !  'BP 3F P' 

              ELSE
                 CALL BP3F(.FALSE.,.FALSE.) !  'BP 3F' 
              ENDIF

           CASE ('3N')                      ! 'BP 3N'      
              CALL NN4        

           CASE ('32')                      !  'BP 32'     
              IF (FCHAR(6:6) == 'N') THEN
                ! UNDOC. OPERATION (HIGH MEMORY)
                CALL NN24                   ! 'BP 32N' 
              ELSE
                CALL BP32F(.FALSE.)
              ENDIF

           CASE ('CG')                      !  'BP CG'   
             CALL BPCG

           CASE ('3G')                      !  'BP CG'  BUGGY??
                CALL WIW3G

           CASE DEFAULT
              CALL ERRT(101,'UNKNOWN/OBSOLETE OPERATION',NDUM)
           END SELECT

           RETURN
        END SELECT

C----------------------------------------------------------------------  

9999    CLOSE(LUN1)
        CLOSE(LUN2)

        END
