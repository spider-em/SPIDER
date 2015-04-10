C++*********************************************************************
C
C    FOUR1B.  
C                                   OPFILEC          FEB 03 ARDEAN LEITH    
C
C	
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
C  PURPOSE: CALCULATES POWER SPECTRUM
C
C IMAGE_PROCESSING_ROUTINE
C
C    ---------------- POWER SPECTRUM ------------------------ 'PW'
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

	SUBROUTINE FOUR1B

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: QA

        COMMON /COMMUN/ FILNAM
        CHARACTER (LEN=MAXNAM) ::  FILNAM

	DATA LUN1,LUN2/21,22/

C       FOURIER MODULI FROM COMPLEX FOURIER TRANSFORM 
C       ON DISK FOR 2-D OR 3-D PICTURES.
        MAXIM = 0
	CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &		   MAXIM,'INPUT',.TRUE.,IRTFLG)
	IF (IRTFLG .NE. 0) RETURN

C       FOURIER INPUT 2D & 3D
        NSLICE = MAX0(1,NSLICE)
	IFI=IFORM
        IF (IFORM.EQ.1)  THEN
           IF(MOD(NSAM,2).EQ.0)  THEN
              LS=NSAM+2
           ELSE
              LS=NSAM+1
           ENDIF
           IRL=NSAM
          
        ELSEIF (IFORM.EQ.3)  THEN
           IF(MOD(NSAM,2).EQ.0)  THEN
              LS=NSAM+2
           ELSE
              LS=NSAM+1
           ENDIF
           IRL=NSAM
          
        ELSEIF (IFORM.EQ.-11)  THEN
           IFORM=1
           LS=NSAM
           IRL=NSAM
           NSAM=NSAM-1
          
        ELSEIF (IFORM.EQ.-12)  THEN
           IFORM=1
           LS=NSAM
           IRL=NSAM
           NSAM=NSAM-2
           
        ELSEIF (IFORM.EQ.-21)  THEN
           IFORM=3
           LS=NSAM
           IRL=NSAM
           NSAM=NSAM-1
           
        ELSEIF(IFORM.EQ.-22)  THEN
           IFORM=3
           LS=NSAM
           IRL=NSAM
           NSAM=NSAM-2
           
        ELSE
           CALL ERRT(2,'FT',NE)
           CLOSE(LUN1)
           RETURN
        ENDIF

        MAXIM = 0
	CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'N',IFORM,NSAM,NROW,NSLICE,
     &		       MAXIM,'OUTPUT',.FALSE.,IRTFLG)
	IF (IRTFLG .NE. 0)  RETURN

 	ALLOCATE (QA(LS,NROW,NSLICE), STAT=IRTFLG)
        IF (IRTFLG .EQ. 0) THEN 
          
	   CALL READV(LUN1,QA,LS,NROW,IRL,NROW,NSLICE)
	   CLOSE(LUN1)

	   IF (IFI .GT. 0) THEN
C             REAL INPUT 2D & 3D
	      INV=+1
              IF (NSLICE.EQ.1)  THEN
	         CALL  FMRS_2(QA,NSAM,NROW,INV)
	      ELSE
	         CALL  FMRS_3(QA,NSAM,NROW,NSLICE,INV)
	      ENDIF
	      IF (INV.EQ.0)  THEN
                 CALL  ERRT(38,'PW',NE)
                 CLOSE(LUN1)
                 CLOSE(LUN2)
                 DEALLOCATE (QA)                   
                 RETURN
	      ENDIF
	   ENDIF
           
           IF (NSLICE.EQ.1)  THEN
              CALL PW2SR(QA,NSAM,NROW,FCHAR(4:4))
	   ELSE
              CALL PW3SR(QA,NSAM,NROW,NSLICE,FCHAR(4:4))
           ENDIF
	   
           CALL WRITEV(LUN2,QA,LS,NROW,NSAM,NROW,NSLICE)
	   CLOSE(LUN2)

        ELSEIF(IFI.GT.0)  THEN
C          REAL INPUT FOR ON DISK VERSION NOT SUPPORTED
           CALL ERRT(6,'PW',NE)
           CLOSE(LUN1)
           CLOSE(LUN2)
           DEALLOCATE (QA)   
           RETURN

	ELSE
C          ON DISK VERSION
           WRITE(NOUT,*)' ** WARNING: SLOW ON-DISK VERSION USED.'

	   IF (IFORM.EQ.1)  THEN
	      CALL  PW2SDR(LUN1,LUN2,NSAM,NROW,FCHAR(4:4))
	   ELSE
	      CALL  PW3SDR(LUN1,LUN2,NSAM,NROW,NSLICE,FCHAR(4:4))
	   ENDIF

	   CLOSE(LUN1)
	   CLOSE(LUN2)
	ENDIF

        DEALLOCATE (QA)        
        END
