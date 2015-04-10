C++*********************************************************************
C
C  INQUIREALLOC.F                 NEW ROUTINE  JAN 2001 ArDean Leith
C                                 REWRITE      JAN 2006 ArDean Leith
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
C   INQUIREIF    DETERMINES AMOUNT OF ALLOCABLE MEMORY
C
C   PARAMETERS:  RSTART   STARTING GB AMOUNT
C                IGOT     ALLOCATABLE AMOUNT MB
C                SAYIT    LOGICAL FLAG TO LIST AMOUNT

C--*********************************************************************

        SUBROUTINE INQUIREALLOC(RSTART,IGOT,SAYIT,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
 
        LOGICAL    :: SAYIT
#ifndef SP_32
        INTEGER *8    IASK8,IOK
#else
        INTEGER    :: IASK8,IOK
#endif
        INTEGER    :: IASK4

        DOUBLE PRECISION                   :: RASK,RGOT,GASK
        INTEGER, ALLOCATABLE, DIMENSION(:) :: IBUF

        IRTFLG = 1

C       -n32 ALLOCATION ON SGI LIMITED TO:  2147483647 

C       CONVERT FROM GB TO BYTE UNITS
        RASK  = RSTART * DBLE(1e9) 
        IASK8 = RASK

        CALL BIGALLOC(IASK8,IOK,.TRUE.,.FALSE.,IRTFLG)
        IOK   = MIN(IOK,IASK8)

        GASK  = DBLE(IOK) / DBLE(1e9) 
        WRITE(NOUT,*) ' STARTING WITH REQUEST FOR GBYTES: ',GASK

C       CONVERT STARTING REQUEST TO 4 BYTE (INTEGER) UNITS
        IASK8 = DBLE(IOK) / DBLE(4) 

        DO
           GASK  = DBLE(4) * DBLE(IASK8) /  DBLE(1e9) 
           WRITE(NOUT,*) ' REQUESTING GBYTES: ',GASK
           IF (NOUT .NE. NDAT) WRITE(NDAT,*)' REQUESTING GBYTES: ',GASK 

           ALLOCATE(IBUF(IASK8),STAT=IRTFLG)

           IF (IRTFLG .EQ. 0) THEN
              WRITE(NOUT,*)' ALLOCATED GBYTES: ',GASK
              IF (NOUT .NE. NDAT) WRITE(NDAT,*)
     &                     ' ALLOCATED GBYTES: ',GASK 
              RGOT = DBLE(IASK8) * DBLE(4) /  DBLE(1e9)
              IGOT = RGOT
              EXIT
           ELSE
              RASK  = DBLE(IASK8) * DBLE(0.8) 
              IASK8 = RASK
              IGOT  = 0
           ENDIF
        ENDDO     
  
        IF (ALLOCATED(IBUF)) DEALLOCATE(IBUF)

        RETURN
        END






