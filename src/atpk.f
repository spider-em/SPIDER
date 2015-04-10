C	
C ++********************************************************************
C
C ATPK    OPENDOC USED                           OCT  2013 ArDean Leith
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
C ATPK                                                                   
C
C PURPOSE: PEAK PICKING
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE ATPK(LUN1,NX,NY,NZ)
	
        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

        INTEGER               :: LUN1,NX,NY,NZ

        REAL, ALLOCATABLE     :: XYZ(:,:)
        CHARACTER(LEN=MAXNAM) :: DOCNAM
 
        REAL                  :: DLIST(3)
        LOGICAL               :: ADDEXT,GETNAME,ISOLD
        LOGICAL               :: APPEND,MESSAGE,NEWFILE
        CHARACTER(LEN=90)     :: COMMENT

        INTEGER               :: NGH,NOT_USED,NQ,ITMP,IRTFLG
        INTEGER               :: LUNDOCNO,ne,nnx,nny,l,nlet,nlist,key,i1
        REAL                  :: THRSH

        INTEGER, PARAMETER    :: LUNDOCN = 55

        CALL RDPRI1S(NGH,NOT_USED,
     &              'PIXEL NEIGHBOURHOOD FOR SEARCH',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
 
        IF (NGH <. 3) THEN
           CALL ERRT(101,'NEIGHBOURHOOD MUST BE >= 3',NE)
           RETURN
        ENDIF

        NQ = NGH / 2
       
        THRSH = 0.7
        CALL RDPRM1S(THRSH,NOT_USED,'CC THRESHOLD FOR PEAKS',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN 

        WRITE(NOUT,*) '   '
        WRITE(NOUT,*)
     &        ' PEAKS CLOSE TO EDGE OF MICROGRAPH CAN BE EXCLUDED'

        CALL RDPRI2S(NNX,NNY,NOT_USED,
     &               'MICROGRAPH EDGE DIMENSION',IRTFLG)
        WRITE(NOUT,*)'   '

C       MAXIMUM NUMBER OF POSSIBLE PEAKS 
        ITMP = (NX/NGH +1) * (NY/NGH +1)

       !write(6,*) 'ngh, itmp:',ngh,itmp

       CALL FLUSHRESULTS

        ALLOCATE (XYZ(3,ITMP), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'XYZ',3*ITMP)
           RETURN
        ENDIF

        CALL PKD(LUN1,NX,NY,NQ,XYZ,ITMP,THRSH,L, NNX,NNY)

        ADDEXT  = .TRUE.
        GETNAME = .TRUE.
        ISOLD   = .FALSE.
        APPEND  = .FALSE.
        MESSAGE = .TRUE.
        IRTFLG  = -8         ! NO IC USE

        CALL OPENDOC(DOCNAM,ADDEXT,NLET,LUNDOCN,LUNDOCNO,GETNAME,
     &           'PEAK LOCATION DOC',ISOLD,APPEND,MESSAGE,
     &            NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C                  123456789 123456789 123456789 123456789 123456789 123456789 
        COMMENT = '          NUMBER OF PEAKS     '
        CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:30),IRTFLG)

        NLIST    = 2
        KEY      = -1
        DLIST(1) = L
        CALL LUNDOCWRTDAT(LUNDOCNO,KEY,DLIST,1,IRTFLG)

C                  123456789 123456789 123456789 123456789 123456789 123456789 
        COMMENT = '  X          Y          Z     '
        CALL LUNDOCPUTCOM(LUNDOCNO,COMMENT(1:30),IRTFLG)

        DO I1 = 1,L
           KEY      = I1
           DLIST(1) = XYZ(1,I1)
           DLIST(2) = XYZ(2,I1)
           DLIST(3) = XYZ(3,I1)
           CALL LUNDOCWRTDAT(LUNDOCNO,KEY,DLIST,3,IRTFLG)
        ENDDO
	CLOSE(LUNDOCN)

        IF (ALLOCATED(XYZ)) DEALLOCATE(XYZ)

        END


