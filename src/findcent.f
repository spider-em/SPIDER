C++*********************************************************************
C
C FINDCENT.F
C                 ADDED CENT_SYM           JAN 2012 GREGORY KISHCHENKO *
C                 RENAMED FROM CENT        MAR 2012 ARDEAN LEITH       *
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
C  FINDCENT()
C
C  PURPOSE: DETERMINATION OF CENTER OF OBJECT USING PHASE OR SYMMETRY
C           INFO                         
C
C--*********************************************************************

         SUBROUTINE FINDCENT()

         IMPLICIT NONE

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC' 
 
         CHARACTER(LEN=MAXNAM)  :: FILNAM
         REAL, ALLOCATABLE      :: AIMG(:,:),F1(:,:)

         INTEGER, PARAMETER     :: LUN1 = 21

         INTEGER                :: ICOMM,MYPID,MPIERR
         INTEGER                :: MAXIM,ITYPE,NX,NY,NZ
         INTEGER                :: IRTFLG,M,NS,NR,NC,NXLD
         REAL                   :: SNS,SNR,SNC
         LOGICAL                :: SHIFTIT = .FALSE.

         CALL SET_MPI(ICOMM,MYPID,MPIERR)  ! SETS ICOMM AND MYPID

         MAXIM = 0
         CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,NX,NY,NZ,
     &                MAXIM,'INPUT',.FALSE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
         
         IF (ITYPE == 1) THEN
C           IMAGE ----------------------------------------- IMAGE

            IF (FCHAR(4:4) == 'S') THEN
C              FIND QUASI-SYMMETRY CENTER ----------------- IMAGE SYM

               NXLD = NX + 2 - MOD(NX,2)

               ALLOCATE (AIMG(NXLD,NY), 
     &                   F1  (NXLD,NY), STAT=IRTFLG)
               IF (IRTFLG .NE. 0) THEN 
                  CALL ERRT(46,'CENT, AIMG,F1',2*NXLD*NY)
                  GOTO 505
               ENDIF

               CALL READV(LUN1,AIMG,NXLD,NY, NX,NY,1) 
          
               CALL CENT_SYM(AIMG,F1,SHIFTIT,
     &                       NXLD, NX,NY,SNS,SNR,IRTFLG)
               IF (IRTFLG .NE. 0) GOTO 505

            ELSE
C              FIND PHASE CENTER --------------------------- IMAGE PH

               ALLOCATE (AIMG(NX,NY), STAT=IRTFLG)
               IF (IRTFLG .NE. 0) THEN 
                  CALL ERRT(46,'CENT, AIMG',NX*NY)
                  GOTO 505
               ENDIF

               CALL READV(LUN1,AIMG,NX,NY, NX,NY,1)
           
               CALL CENT_PH(AIMG, NX,NY,SNS,SNR)

               M   = NX/2+1
               SNS = SNS - M
               M   = NY/2 + 1
               SNR = SNR - M
            ENDIF

            NS  = INT(SNS)
            NR  = INT(SNR)

            IF (VERBOSE .AND. MYPID <= 0 ) THEN
               WRITE(NOUT,*) ' '
               WRITE(NOUT,2222) NS,NR, SNS,SNR
2222           FORMAT('  Approximate center offset:',/,
     &         '             X            Y',/,
     &         '         ', I5,'        ',I5,/,
     &         '          ',F8.2,'     ',F8.2,/)

            ELSEIF (MYPID <= 0) THEN
               WRITE(NOUT,2223) NS,NR, SNS,SNR
2223           FORMAT('  Approximate center offset: ',
     &                '(',I5,','I5,    ')  == ',
     &                '(',F8.2,','F8.2,')')
            ENDIF

            CALL REG_SET_NSEL(1,4,FLOAT(NS),FLOAT(NR),SNS,SNR,
     &                        0.0,IRTFLG)


         ELSEIF(ITYPE .NE. 3 )  THEN
C           FOURIER FILE? ------------------------------
            CALL ERRT(102,
     &      'OPERATION NOT IMPLEMENTED FOR THIS DATA FORMAT',ITYPE)
            GOTO 505

         ELSEIF(ITYPE == 3 .AND. FCHAR(4:4) == 'S' )  THEN
C           VOLUME ------------------------------------
            CALL ERRT(101,
     &          'OPERATION NOT IMPLEMENTED FOR VOLUMES',IRTFLG)
            GOTO 505

         ELSEIF(ITYPE == 3)  THEN
C           VOLUME ------------------------------------------- VOLUME

            CALL CENT_3PH(LUN1,NX,NY,NZ, SNS,SNR,SNC)

            M   = NX/2+1
            SNS = SNS-M
            M   = NY/2+1
            SNR = SNR-M
            M   = NZ/2+1
            SNC = SNC-M
            NS  = INT(SNS)
            NR  = INT(SNR)
            NC  = INT(SNC)

            IF (VERBOSE .AND. MYPID <= 0) WRITE(NOUT,*) ' '
            IF (VERBOSE .AND. MYPID <= 0 ) THEN
               WRITE(NOUT,3333) NS,NR,NC, SNS,SNR,SNC
3333        FORMAT('  Approximate center offset:',/,
     &      '         X             Y             Z',/,
     &      '     ', I5,'         ',I5,'         ',I5,/,
     &      '        ',F8.2,'      ',F8.2,'      ',F8.2,/)
            ELSEIF (MYPID <= 0) THEN
               WRITE(NOUT,3334) NS,NR,NC, SNS,SNR,SNC
3334           FORMAT('  Approximate center offset: ',
     &                '(',I5,','I5,','I5,    ')  == ', 
     &                '(',F8.2,','F8.2,','F8.2,')')
            ENDIF

            CALL REG_SET_NSEL(1,5,FLOAT(NS), FLOAT(NR), FLOAT(NC),
     &                            SNS, SNR, IRTFLG)
            CALL REG_SET_NSEL(6,1,SNC, 0.0, 0.0, 0.0, 0.0,IRTFLG)
         ENDIF

505      CLOSE(LUN1)
         IF (ALLOCATED(AIMG)) DEALLOCATE(AIMG)
         IF (ALLOCATED(F1))   DEALLOCATE(F1)

         END
