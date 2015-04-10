C++*********************************************************************
C
C  COPY1.F                   REMOVED FROM UTIL2    JUL 96 ArDean Leith
C                            ADDED NT INTERFACES   OCT 98 ArDean Leith
C                            REMOVED 'CP FROM LUM' FEB 99 ArDean Leith
C                            ADDED 'CP FROM XP'    OCT 00 Pawel Penczek
C                            'CP TO NT' ON NT      JUN 01 ArDean Leith
C                            'CP FROM SG' ON NT    SEP 01 ArDean Leith
C                            'CP TO/FROM CCP4'     FEB 02 ArDean Leith
C                            'CP FROM EMI'         FEB 02 ArDean Leith
C                            INDEXED STACK         JAN 02 ArDean Leith
C                            'CP TO OPEND'         FEB 03 ArDean Leith
C                            'CP TO' NORM          JUL 03 ArDean Leith
C                            'CP FIX' REMOVED      OCT 03 ArDean Leith
C                            MPI                   OCT 03 Chao Yang
C                            'CP TO SF3' GONE      DEC 04 ArDean Leith
C                            'CP FROM NIK'         JAN 05 ArDean Leith
C                            'CP TO PDS' GONE      MAR 09 ArDean Leith
C                            'CP TO SGI' GONE      MAY 09 ArDean Leith
C                            'CP TO VV'  GONE      DEC 10 ArDean Leith
C                            'CCP4' NOW JUST MRC   FEB 12 ArDean Leith
C                            'CP TO JPG'           APR 13 ArDean Leith
C                            'CP TO VOL'           MAY 13 ArDean Leith
C                            'CP FROM TIF'         MAR 14 ArDean Leith
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
C  COPY1(MAXDIM)
C
C  PARAMETERS:      MAXDIM         MAX LENGTH FOR UNLABELED COMMON
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE COPY1(MAXDIM)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM) :: FILOLD,FILNEW
 
        CHARACTER(LEN=1)      :: NULL = CHAR(0)
        LOGICAL               :: INDXD,FLIPOUT

        INTEGER,PARAMETER     :: IDELAY  = 3

        INTEGER, PARAMETER    :: LUN1    = 14 
        INTEGER, PARAMETER    :: LUN2    = 15 
        INTEGER, PARAMETER    :: LUNT    = 16 
        INTEGER, PARAMETER    :: LUNDOC  = 81 
        INTEGER, PARAMETER    :: LUNDOC2 = 82 
        INTEGER, PARAMETER    :: LUNXM1  = 83 
        INTEGER, PARAMETER    :: LUNXM2  = 84 

C       FROM/'AS','MR','PD','RA','NI','TE','VA','EM','NT','XP'/
  
C       TO/'AS','BR','MR','PO','RA','TI','XP','OP','JPG'/  

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

        IF (FCHAR(4:7) == 'TO O') THEN !----------------- 'CP TO OPEND'
C          STANDARD COPY WITH FLIPPED ENDEDNESS
           INDXD   = .FALSE.
           FLIPOUT = .TRUE.
           CALL COPYD(LUN1,LUN2,LUNDOC,LUNXM1,LUNXM2,INDXD,FLIPOUT)
           GOTO 9000

        ELSEIF (FCHAR(4:7) == 'TO V') THEN !--------------- 'CP TO VOL'
C          SPIDER IMAGE(S) FILE INTO VOLUME FILES '
           CALL COPYTOVOL()
           GOTO 9000
        ENDIF

C                                    12345678     1234567890
        SELECT CASE(FCHAR(4:5))   ! 'CP TO XX' , 'CP FROM XX'

#if defined(SP_IFC) || defined(__INTEL_COMPILER) 
        CASE ('I')      ! --------------------------------------- 'CP I'
#else
#if defined(__GFORTRAN__) 
        CASE ('I')      ! --------------------------------------- 'CP I'
#else
        CASE ('I','I ') ! --------------------------------------- 'CP I'
#endif
#endif

C          STANDARD COPY OR STANDARD COPY TO INDEXED STACK
           INDXD   = .TRUE.
           FLIPOUT = .FALSE.
           CALL COPYD(LUN1,LUN2,LUNDOC,LUNXM1,LUNXM2,INDXD,FLIPOUT)

        CASE ('')    ! -------------------------------------------- 'CP'
C          STANDARD COPY 
           INDXD   = .FALSE.
           FLIPOUT = .FALSE.
           CALL COPYD(LUN1,LUN2,LUNDOC,LUNXM1,LUNXM2,INDXD,FLIPOUT)



        CASE ('TO')  ! -------------------------------- ----- 'CP TO **'

C          OPEN INPUT FILE, WHOLE STACK NOT ALLOWED
           MAXIM = 0
           CALL OPFILEC(0,.TRUE.,FILOLD,LUN1,'O',ITYPE,
     &                   NX1,NY1,NZ1,MAXIM,'SPIDER INPUT',
     &                   .FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9000

           IF (IMAMI .NE. 1) 
     &         CALL NORM3(LUN1,NX1,NY1,NZ1,FMAX,FMIN,AV)

           SELECT CASE(FCHAR(7:8)) 

           CASE ('AS') 
C             FROM SPIDER IMAGE INTO EDITABLE IMAGE ------ 'CP TO ASCII'
              CALL COPYE(LUN1,LUN2,NX1,NY1,NZ1)

           CASE ('BR')
C             FROM 3D SPIDER FILE TO BRIX FORMAT----------- 'CP TO BRIX'
              CALL COPYBRIX(LUN1,LUN2,NX1,NY1,NZ1)

           CASE ('MR','CC')
C             FROM SPIDER IMAGE FILE INTO MRC FORMAT ------- 'CP TO MRC'
              CALL COPYCCP4(LUN1,LUN2,NX1,NY1,NZ1)
              ! OBSOLETE CALL COPYMRC(LUN1,LUN2,NX1,NY1,NZ1)

           CASE ('PO')
C             FROM SPIDER TO POSTSCRIPT IMAGE (8 BIT)------ 'CP TO POST'
              CALL COPYPOS(FILOLD,LUN1,LUN2,NX1, NY1, NZ1)

           CASE ('RA')
C             FROM SPIDER IMAGE FILE INTO RAW IMAGE FILE --- 'CP TO RAW'
              CALL COPYU(LUN1,LUN2,NX1,NY1,NZ1)
 
           CASE ('TI')
C             FROM SPIDER FORMAT TO TIFF FORMAT ----------- 'CP TO TIFF'
              CALL SPDTOTIFF(LUN1,LUN2,NX1,NY1,NZ1,IRTFLG)
 
           CASE ('XP')
C             SPIDER IMAGE FILE INTO XPLOR FILES --------- 'CP TO XPLOR'
              CALL COPYTOXPLOR(LUN1,FILOLD,LUN2,NX1,NY1,NZ1)

           CASE ('JP')
C             SPIDER IMAGE FILE INTO JPG FILES ------------- 'CP TO JPG'
              VERBOSET = VERBOSE  ! FROM CMBLOCK
              FILNEW   = NULL
              CALL COPYTOJPG(LUN1,LUNT,FILNEW,NX1,NY1,NZ1,
     &                       VERBOSET,IDELAY)

           CASE DEFAULT
C             NO SUCH COPY FUNCTION
              CALL ERRT(101,'NO SUCH CP OPERATION, CHECK MENU',IDUM)
          
           END SELECT



        CASE ('FR') !  CP FROM OPERATIONS ----------------------------

           SELECT CASE(FCHAR(9:10))

           CASE ('AS')
C             EDITABLE IMAGE FILE TO SPIDER IMAGE FILE - 'CP FROM ASCII'
              CALL COPYF(LUN1,LUN2)

           CASE ('MR','CC')
C             FROM MRC FORMAT TO SPIDER ------------------ 'CP FROM MRC'
              CALL COPYCCP4(LUN1,LUN2,NSAM,NROW,NSLICE)
              ! NO LONGER SUPPORTED: COPYMRC(LUN1,LUN2,NX,NY,NZ)

           CASE ('PD')
C             FROM PDB FILE TO SPIDER VOLUME FILE -------- 'CP FROM PDB'
              CALL READPDB

           CASE ('RA')
C             COPY RAW IMAGE FILE INTO SPIDER IMAGE FILE - 'CP FROM RAW'
              CALL RAWTOSPIDER(LUN1,LUN2,IRTFLG)

           CASE ('NI')
C             NIKON TIFF IMAGE FILE INTO SPIDER FILE ----- 'CP FROM NIK'
              CALL COPYFROMNIK(LUN1,LUN2,IRTFLG)

           CASE ('TE')
C             FROM TERMINAL INTO SPIDER IMAGE FILE ------ 'CP FROM TERM'
              CALL COPYR(LUN2)
          
           CASE ('VA')
C             VAX SPIDER TO UNIX SPIDER ------------------ 'CP FROM VAX'
              FILOLD = NULL
              FILNEW = NULL
              CALL VAXTOUNIX(FILOLD,FILNEW,LUN1,LUN2,IRTFLG)
 
           CASE ('EM')
C             EMI FORMAT TO SPIDER ----------------------- 'CP FROM EMI'
              CALL COPYEMI(LUN1,LUN2)

           CASE ('NT')
C             FROM NT FORMAT TO SPIDER -------------------- 'CP FROM NT'
              CALL ERRT(101,'USE OPERATION: <CP TO OPEND> ',IDUM)

           CASE ('XP')
C             FROM XPLO FORMAT TO SPIDER ------------------ 'CP FROM XP'
              CALL COPYFROMXPLO(LUN1,LUN2,MAXDIM)

           CASE ('SG')
C             FROM TO SGI BYTE FORMAT --------------------- 'CP FROM SG'
              CALL ERRT(101,"USE OPERATION: <CP TO OPEND> ",IDUM)

           CASE ('TI')
C             FROM TIFF TO SPIDER FORMAT ----------------- 'CP FROM TIF'
               CALL COPYFROMTIF(LUN1,LUN2,IRTFLG)


           CASE DEFAULT

C             NO SUCH COPY FUNCTION
              CALL ERRT(101,"NO SUCH 'CP' OPERATION, CHECK MENU",IDUM)

           END SELECT

        CASE DEFAULT
           CALL ERRT(101,'NO SUCH CP OPERATION, CHECK MENU',IDUM)

        END SELECT

9000    CLOSE(LUN1)
        CLOSE(LUN2)

        END
