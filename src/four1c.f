C++*********************************************************************
C
C FOUR1C.F
C                                   OPFILEC          FEB 03 ARDEAN LEITH    
C                                   CLEANED          JUL 08 ARDEAN LEITH    
C **********************************************************************
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
C  FOUR1C
C
C  PURPOSE:  FORWARD OR REVERSE FOURIER TRANSFORMS AN IMAGE OR VOLUME.
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE FOUR1C

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 
        
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: QA

        CHARACTER(LEN=MAXNAM)               :: FILNAM

        DATA  LUN1,LUN2/21,22/

        MAXIM  = 0      ! NO BARE STACK SUPPORT
        MAXIM2 = 0      ! NO BARE STACK SUPPORT
        IRTFLG = 0

C       OPEN INPUT FILE.  CAN BE FOURIER OR NOT
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &             MAXIM,'INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0)  RETURN

        INV = ISIGN(1,IFORM)                ! FORWARD OR REVERSE

        IF (IFORM.EQ.1)  THEN
C          NON-FOURIER IMAGE
           IF (MOD(NSAM,2) .EQ. 0)  THEN
              IFORM = -12
              LS    = NSAM+2
           ELSE
              IFORM = -11
              LS    = NSAM+1
           ENDIF
           
        ELSEIF (IFORM.EQ.3)  THEN
C          NON-FOURIER VOLUME
           IF (MOD(NSAM,2) .EQ. 0)  THEN
              IFORM = -22
              LS    = NSAM+2
           ELSE
              IFORM = -21
              LS    = NSAM+1
           ENDIF
           
        ELSEIF (IFORM .EQ. -11)  THEN
C          FOURIER IMAGE
           IFORM = 1
           LS    = NSAM-1
           
        ELSEIF (IFORM .EQ. -12)  THEN
C          FOURIER IMAGE
           IFORM = 1
           LS    =NSAM-2
           
        ELSEIF (IFORM .EQ. -21)  THEN
C          FOURIER VOLUME
           IFORM = 3
           LS    = NSAM-1
 
        ELSEIF (IFORM .EQ. -22)  THEN
C          FOURIER VOLUME
           IFORM = 3
           LS    = NSAM-2
                                               
        ELSE
C          UNKNOWN FORMAT
           CALL ERRT(2,'FT',NE)
           CLOSE(LUN1)
           RETURN
        ENDIF
        
        ISPACE = MAX(NSAM,LS)
        ISAM   = MIN(NSAM,LS)

        CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'N',IFORM,LS,NROW,NSLICE,
     &             MAXIM2,'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CLOSE(LUN1)
           RETURN
        ENDIF

 	ALLOCATE (QA(ISPACE,NROW,NSLICE), STAT=IRTFLG)
 	
        IF (IRTFLG .EQ. 0) THEN
C          ADEQUATE MEMORY, READ INPUT ARRAY
           CALL READV(LUN1,QA,ISPACE,NROW,NSAM,NROW,NSLICE)
	   CLOSE(LUN1)

           IF (IFORM .EQ. 3 .OR. IFORM .LT. -20)  THEN
C             VOLUME
              CALL FMRS_3(QA,ISAM,NROW,NSLICE,INV)
           ELSE
C             UNAGE
              CALL FMRS_2(QA,ISAM,NROW,INV)
           ENDIF

           IF (INV .EQ. 0) THEN
              CALL ERRT(38,'FT',NE)
           ELSE
C             SAVE OUTPUT ARRAY                             
              CALL WRITEV(LUN2,QA,ISPACE,NROW,LS,NROW,NSLICE)	             
           ENDIF

        ELSE
C          INADEQUATE MEMORY
               
           NNNN = ISAM+2-MOD(ISAM,2)
           WRITE(NOUT,*)' WARNING: USING SLOW ON-DISK VERSION OF FFT.'

           IF (NSLICE .LE. 1)  THEN
C             2D FFT 
              NC   = 2
              LR   = ISPACE / NC
              CALL FMRS_2DR(LUN1,LUN2,LR,NNNN,ISAM,NROW,INV)
             
           ELSE
C             3D FFT  
              CALL FMRS_3DR(LUN1,LUN2,NNNN,ISAM,NROW,NSLICE,INV)      
           ENDIF
        ENDIF
        
9999    CLOSE(LUN2)
        IF (ALLOCATED(QA)) DEALLOCATE (QA) 
       
        END
