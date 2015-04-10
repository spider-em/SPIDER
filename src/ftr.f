C++*********************************************************************
C
C    FTR.F                                         BIMAL RATH 8/14/2001         
C                  OPFILEC                         FEB 03 ARDEAN LEITH
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
C FTR
C
C IMAGE_PROCESSING_ROUTINE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C-************************************************************************

        SUBROUTINE FTR

        INCLUDE 'CMBLOCK.INC'      
        INCLUDE 'CMLIMIT.INC' 
 
        CHARACTER(LEN=MAXNAM)   ::   FILNAM1,FILNAM2,FILNAM3
  
        REAL, ALLOCATABLE, DIMENSION(:,:,:) ::    AIMG,BIMG,CIMG
 
        CHARACTER*1       NULL

        DATA  LUN1,LUN2,LUN3/21,22,23/

C       INPUT FIRST IMAGE
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM1,LUN1,'O',ITYPE1,NSAM,NROW,
     &          NSLICE,MAXIM,'FIRST INPUT',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (ITYPE1 .GT. 1)THEN
           CLOSE(LUN1)
           CALL ERRT(2,'RF',NE)
           RETURN
        ENDIF

C       INPUT SECOND IMAGE
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM2,LUN2,'O',ITYPE2,NSAM2,NROW2,
     &          NSLICE2,MAXIM,'SECOND INPUT',.TRUE.,IRTFLG)     

        IF (IRTFLG .NE. 0) THEN
           CLOSE(LUN1)
           RETURN

        ELSEIF (ITYPE1 .NE. ITYPE2) THEN
           CALL ERRT(2,'RF',NE)
           GOTO 9999

        ELSEIF (NSAM.NE.NSAM2 .OR. NROW.NE.NROW2) THEN        
           CALL ERRT(1,'RF',NE)
           GOTO 9999
        ENDIF
        
C       CHECK FOURIER INPUT FILES        
        IF ((ITYPE1 .GT. 0) .OR. (ITYPE2 .GT. 0)) THEN
           CALL ERRT(2,'RF',NE)
           GOTO 9999 
        ENDIF   
        
C       OUTPUT IMAGE

        MAXIM = 0
        CALL OPFILEC(LUN1,.TRUE.,FILNAM3,LUN3,'U',ITYPE1,NSAM,NROW,
     &          NSLICE,MAXIM,'OUTPUT',.TRUE.,IRTFLG)

        ALLOCATE (AIMG(NSAM,NROW,NSLICE),BIMG(NSAM,NROW,NSLICE), 
     &                       CIMG(NSAM,NROW,NSLICE),STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'RF, AIMG & BIMG & CIMG',IER)
           GOTO 9999
        ENDIF

        CALL READV(LUN1,AIMG,NSAM,NROW,NSAM,NROW,NSLICE)
        CALL READV(LUN2,BIMG,NSAM,NROW,NSAM,NROW,NSLICE)
 
        DO K=1,NSLICE
           DO J=1,NROW
              DO I=1,NSAM,2
                    IF(BIMG(I,J,K).NE.0.0)THEN
                       PHB=ATAN2(BIMG(I+1,J,K),BIMG(I,J,K))
                    ELSE
                       PHB=0.0
                    ENDIF 
                    QA=SQRT(AIMG(I,J,K)**2+AIMG(I+1,J,K)**2)
                    CIMG(I,J,K)   =  QA*COS(PHB)
                    CIMG(I+1,J,K) =  QA*SIN(PHB)
              ENDDO
           ENDDO
        ENDDO
        
        CALL WRITEV(LUN3,CIMG,NSAM,NROW,NSAM,NROW,NSLICE)         

9999    CLOSE(LUN1)
        CLOSE(LUN2)
        CLOSE(LUN3)
        IF (ALLOCATED(AIMG))  DEALLOCATE (AIMG)
        IF (ALLOCATED(BIMG))  DEALLOCATE (BIMG)
        IF (ALLOCATED(CIMG))  DEALLOCATE (CIMG)        
        RETURN
        END
        

        
