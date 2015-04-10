C++*********************************************************************
C
C  COPYTOXPLOR.F  --              CREATED    OCT 99     PAWEL PENCZEK                
C                                 REFINED    10/21/99   BIMAL RATH
C                                 KX... +1   03/08/06   ARDEAN LEITH
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
C COPYTOXPLOR(LUN1,FIOLD,LUN2,NSAM,NROW,NSLICE)
C 
C PURPOSE:     CONVERTS SPIDER IMAGE FILE TO EXPLORER FORMAT
C  
C PARAMETERS:
C        LUN1                LOGICAL UNIT NUMBER OF INPUT IMAGE
C        FIOLD               INPUT FILE NAME (WITHOUT EXTENSION)
C        LUN2                LOGICAL UNIT NUMBER OF OUTPUT IMAGE
C        NSAM,NROW,NSLICE    DIMENSIONS OF VOLUME
C        
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE COPYTOXPLOR(LUN1,FILOLD,LUN2,NSAM,NROW,NSLICE)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL, ALLOCATABLE, DIMENSION(:,:) :: AIMG
        CHARACTER(LEN=MAXNAM) :: FILOLD,FILNEW,FILENAME 
        
        CALL RDPRM(A,NOT_USED,'PIXEL SIZE')
	IF (IMAMI .NE. 1) THEN
           CALL NORM3(LUN1,NSAM,NROW,NSLICE,FMAX,FMIN,AV)
        ENDIF

        BMIN   = FMIN
        BMAX   = FMAX
        SUM    = AV
        SIGMA  = SIG
        LENREC = 0

        CALL OPAUXFILE(.TRUE.,FILNEW,DATEXC,LUN2,LENREC,'N',
     &	    'OUTPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       WRITE HEADER
        WRITE(LUN2,201)
201     FORMAT(/,'       2')
        CALL FILNAMANDEXT(FILOLD,DATEXC,FILENAME,NLET,.TRUE.,IRTFLG)
	IF (IRTFLG .NE. 0) RETURN 
	
        WRITE(LUN2,202) FILENAME
202     FORMAT('REMARKS FILENAME= ',A40)

        WRITE(LUN2,203)
cc203     FORMAT('REMARKS CREATED BY SPIDER     CP TO XPLOR')
203     FORMAT('REMARKS CREATED BY MAPMAN V. 960827/4.6.2')

C       AL Mar 06
        IF (MOD(NSAM,2) .EQ. 0) THEN               !EVEN
           KX  = -NSAM/2    + 0
        ELSE
           KX  = -NSAM/2    + MOD(NSAM+1,2)
        ENDIF

        IF (MOD(NROW,2) .EQ. 0) THEN               !EVEN
           KY  = -NROW/2    + 0
        ELSE
           KY  = -NROW/2    + MOD(NROW+1,2)
        ENDIF

        IF (MOD(KX,2) .EQ. 0) THEN               !EVEN
           KZ  = -NSLICE/2    + 0
        ELSE
           KZ  = -NSLICE/2    + MOD(NSLICE+1,2)
        ENDIF

C       BEFORE AL Mar 06
        IEX = NSAM/2
        IEY = NROW/2
        IEZ = NSLICE/2

C       AL Mar 06
        IEX = KX + NSAM   - 1
        IEY = KY + NROW   - 1
        IEZ = KZ + NSLICE - 1

	
        WRITE(LUN2,204) NSAM,KX,IEX,NROW,KY,IEY,NSLICE,KZ,IEZ
204     FORMAT(9I8)

        ANG = 90.0
        WRITE(LUN2,206)  A*NSAM,A*NROW,A*NSLICE,ANG,ANG,ANG
206     FORMAT(6(1PE12.5))

        WRITE(LUN2,205)
205     FORMAT('ZYX')
      
        ALLOCATE (AIMG(NSAM,NROW), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'COPYTOEXPLOR',IER)
           RETURN
        ENDIF

        DO  K=NSLICE,1,-1
           DO J=1,NROW
              CALL  REDLIN(LUN1,AIMG(1,J),NSAM,J+(K-1)*NROW)
           ENDDO

           WRITE(LUN2,504)  NSLICE-K
504        FORMAT(I8)
           WRITE(LUN2,505)  ((AIMG(I,J),J=1,NROW),I=1,NSAM)
505        FORMAT(6E12.5)
        ENDDO
	
        K     = -9999
        WRITE(LUN2,504)  K
        WRITE(LUN2,506)  SUM,SIGMA**2
506     FORMAT(1X,6E12.5)
        DEALLOCATE (AIMG)

        END
