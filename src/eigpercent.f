C **********************************************************************
C
C  EIGPERCENT   NEW                               FEB 13 ARDEAN LEITH
C
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
C   EIGPERCENT()
C
C   PURPOSE:  COPY EIGENVECTORS AND % TO DOC FILE 
C
C   OPERATIONS SUPPORTED : CA EIGDOC  -  
C
C   VARIABLES:        NFAC     NUMBER OF EIGENVECTORS             
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C-----------------------------------------------------------------------------

        SUBROUTINE EIGPERCENT()
        
        IMPLICIT NONE       
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
 
        INTEGER               :: lnblnkn

        CHARACTER(LEN=MAXNAM) :: FILPRE,FILNME,DOCFIL
        CHARACTER(LEN=160)    :: COMMEN
        
        CHARACTER(LEN=1)      :: NULL = CHAR(0)
         
        LOGICAL               :: USE_PCA,NEWFILE
        INTEGER               :: NLET,IRTFLG,NFAC,KIND_PCA,NVAL,NC
        INTEGER               :: NOUTDOC,IKEY
        REAL                  :: SUMP,TRACE
        REAL                  :: DLIST(3)

        INTEGER, PARAMETER    :: LUNE   = 81
        INTEGER, PARAMETER    :: LUNDOC = 82

        CALL FILERD(FILPRE,NLET,NULL,
     &     'CORAN/PCA FILE PREFIX FOR EIGENVECTORS (E.G. CORAN)~',
     &     IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        FILNME = FILPRE(1:NLET) // '_EIG'//NULL

C       OPEN EIGENVALUE FILE
        CALL OPAUXFILE(.FALSE.,FILNME,DATEXC,LUNE,0,
     &                       'O', ' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       READ _EIG FILE HEADER
	READ(LUNE,*) NFAC, SUMP, TRACE, KIND_PCA, NVAL
80      FORMAT(I10,G12.5,1X,G12.5,I10,I10)

        WRITE(NOUT,*)' NUMBER OF FACTORS: ',NFAC

        USE_PCA = (KIND_PCA == 1)

        CALL OPENDOC(DOCFIL,.TRUE.,NLET,LUNDOC,NOUTDOC,.TRUE.,
     &       'EIGENVALUE DOC',.FALSE.,.FALSE.,.TRUE.,
     &       NEWFILE,IRTFLG)

        COMMEN = ' SPIDER doc file with % of variance accounted' //
     &           ' for by factor. Key = Factor'
        NC = lnblnkn(COMMEN)
        CALL LUNDOCPUTCOM(NOUTDOC,COMMEN(1:NC),IRTFLG)

        COMMEN = 'Key     Eigenvalue       %      Cumulative-%'
        NC = lnblnkn(COMMEN)
        CALL LUNDOCPUTCOM(NOUTDOC,COMMEN(1:NC),IRTFLG)

C       READ _EIG EIGEN VALUES & STORE IN DOC FILE

        DO IKEY = 1,NFAC
	   READ(LUNE,*) DLIST(1), DLIST(2), DLIST(3)
           CALL LUNDOCWRTDAT(NOUTDOC,IKEY,DLIST,3,IRTFLG)
        ENDDO

        CALL REG_SET_NSEL(1,1,FLOAT(NFAC),0.0,0.0,0.0,0.0,IRTFLG)

9999    CLOSE (LUNDOC)
        CLOSE (LUNE)

        END




