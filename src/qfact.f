C ++********************************************************************
C                                                                      *
C                                                                      *
C QFACT                 USED FILELIST             APR 2004 ARDEAN LEITH    
C                
C                                                                      *
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
C                                                                      *
C  QFACT
C                                                                    *                                                                     *
C  PURPOSE: ADDS A SERIES OF IMAGES REPRESENTING SINGLE,
C         ALIGNED PARTICLES OR PATCHES OF A CRYSTAL LATTICE, AND
C         COMPUTES A Q FACTOR MAP. IMAGES HAVE TO HAVE EVEN
C         DIMENSIONS (INTERNALLY FFT IS USED)                                                           *
C                                                                      *
C  PARAMETERS:                                                         *
C
C IMAGE_PROCESSING_ROUTINE
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE QFACT(LUN1,LUN2,LUNQ)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
       
        COMMON /COMMUN/          FILNAM,FILA,FILQ,FILPAT
        CHARACTER(LEN=MAXNAM) :: FILNAM,FILA,FILQ,FILPAT        
        INTEGER               :: HREF,KREF
      
C       INPUT IMAGE SERIES
	
        CALL FILELIST(.TRUE.,LUNQ,FILPAT,NLET,INUMBR,NIMAX,NUMT,
     &                 'INPUT FILE TEMPLATE (E.G. PIC****)',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        CALL FILGET(FILPAT,FILNAM,NLET,INUMBR,IRTFLG)

        MAXIM = 0
        CALL OPFILEC(0,.FALSE.,FILNAM,LUN2,'O',ITYPE,NSAM,NROW,NSLICE,
     &               MAXIM,' ',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 980
        CLOSE(LUN2)

        IF( MOD(NSAM,2).NE.0 .OR. MOD(NROW,2).NE.0) GOTO  950
        IF (IFORM.NE.1)  GO TO 970
        IF (NSLICE.GT.1) GO TO 970
        MAXIM  = 0
        ITYPE  = 1
        CALL OPFILEC(0,.TRUE.,FILA,LUN1,'U',ITYPE,
     &            NSAM,NROW,NSLICE,MAXIM,'AVERAGE',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 980


        MAXIM = 0
        ITYPE = 1
        CALL OPFILEC(0,.TRUE.,FILQ,LUNQ,'U',ITYPE,
     &             NSAM,NROW,NSLICE,MAXIM,'Q FACTOR',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 980

C       ENTER FOURIER INDEX FOR VECTOR LISTING (0,0=NO LISTING)

        CALL RDPRMI(HREF,KREF,NOT_USED,'INDICES FOR VECTOR LISTING')

        CALL  QFACT_P(LUN1,LUN2,LUNQ,INUMBR,NUMT,NSAM,NROW,HREF,KREF)

        CLOSE (LUNQ)
        CLOSE (LUN1)
        RETURN


950     IER=10
        GOTO  990

970     IER=2
        GO TO 990

980     IER=4
990     CALL ERRT(IER,'AS F',NE)
        END





