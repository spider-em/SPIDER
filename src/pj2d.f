
C++*********************************************************************
C
C PJ2D.F               NEW                         JUN 18 ARDEAN LEITH
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
C  PJ2D(LUN1,LUN2,LUN3)
C
C  PURPOSE:        2D PROJECTION BY ROWS OR COLUMNS
C
C  PARAMETERS:     LUN1        INPUT FILE
C                  LUN2        OUTPUT FILE
C                  LUN3        OUTPUT FILE
C
C--*********************************************************************

       SUBROUTINE PJ2D(LUN1,LUN2,LUN3)

       IMPLICIT NONE

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC' 

       INTEGER               :: lnblnkn

       INTEGER               :: LUN1,LUN2,LUN3

       CHARACTER(LEN=MAXNAM) :: FILNAM

       CHARACTER *1          :: NULL = CHAR(0)
       CHARACTER *1          :: PROJ_AXIS
    
       INTEGER               :: NX,NY,NZ,MAXIM,IRTFLG,MWANT,NLET
       INTEGER               :: NOT_USED,IX,IY,NE,NXP
 
       REAL, ALLOCATABLE     :: BUF(:,:), PROJ(:), PROJ_VIS(:,:)

       INTEGER, PARAMETER    :: NVISY = 32


C      GET INPUT IMAGE
       CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,
     &                  NX,NY,NZ,
     &                  MAXIM,'INPUT',.FALSE.,IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       IF (NZ > 1) THEN
           CALL ERRT(101,'SORRY, DOES NOT WORK ON VOLUMES',NE)
           RETURN
       ENDIF

       CALL RDPRMC(PROJ_AXIS,NLET,.TRUE.,
     &             'PROJECT ALONG X or Y AXIS? (X/Y)', NULL,IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

C      SET PROJECTION VISUALIZATION SIZE
       NXP   = NY
       IF (PROJ_AXIS == 'Y') THEN
          NXP = NX
       ENDIF

       MWANT = NX * NY + NXP + NXP * NVISY 
       ALLOCATE (BUF(NX,NY),PROJ(NXP),PROJ_VIS(NXP,NVISY),STAT=IRTFLG)
       IF (IRTFLG .NE. 0) THEN
          CALL ERRT(46,'PJ2D, BUF...',MWANT)
          RETURN
       ENDIF

C      LOAD INPUT IMAGE
       CALL REDVOL(LUN1,NX,NY, 1,1, BUF,IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 999

C      CREATE OUTPUT LINE PROJECTION FILE
       IFORM = 1
       MAXIM = 0
       CALL OPFILEC(LUN1,.TRUE.,FILNAM,LUN2,'U',IFORM,
     &                  NXP,1,1,
     &                  MAXIM,'OUTPUT',.FALSE.,IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 999

C      ZERO PROJECTION LINE
       PROJ = 0.0
 
C      PROJECT THROUGH THE IMAGE
       IF (PROJ_AXIS == 'Y') THEN

C         SUM ALONG COLUMNS
          DO IY =1,NY
             PROJ = PROJ + BUF(1:NX,IY) 
          ENDDO

       ELSE
C         SUM ALONG ROWS
          DO IX =1,NX
             PROJ = PROJ + BUF(IX,1:NY) 
          ENDDO
       ENDIF

C      OUTPUT PROJECTION FILE   
       CALL WRTLIN(LUN2,PROJ,NXP,1)

C      DUPLICATE THE PROJECTED LINE IMAGE    
       DO  IY =1,NVISY
         PROJ_VIS(1:NXP,IY) = PROJ 
       ENDDO

C      CREATE PROJECTED LINE VISUALIZATION FILE
       IFORM  = 1
       MAXIM  = 0
       NLET   = lnblnkn(FILNAM)
       FILNAM = FILNAM(1:NLET) // '_vis'
       CALL OPFILEC(0,.FALSE.,FILNAM,LUN3,'U',IFORM,
     &              NXP,NVISY,1,
     &              MAXIM,' ',.FALSE.,IRTFLG)

C      OUTPUT PROJECTION VISUALIZATION FILE   
       CALL WRTVOL(LUN3,NXP,NVISY,1,NXP,PROJ_VIS,IRTFLG)

999    CLOSE(LUN1)
       CLOSE(LUN2)
       CLOSE(LUN3)

       IF (ALLOCATED(BUF))      DEALLOCATE(BUF)
       IF (ALLOCATED(PROJ))     DEALLOCATE(PROJ)
       IF (ALLOCATED(PROJ_VIS)) DEALLOCATE(PROJ_VIS)

       END

