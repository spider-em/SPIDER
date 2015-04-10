
C **********************************************************************
C
C   GOTBRANCH.FOR  -- CREATED OCT 90
C **********************************************************************
C *  AUTHOR: ArDean Leith 
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
C      GOTBRANCH(NCLUS,NOW,IEQUIV,NEQUIV,NEQMAX,IRTFLG)
C
C      PURPOSE:      SETS BRANCH CONNECTIVITY ARRAY  
C
C      PARAMETERS:  
C
C      CALLED BY:    CCONECT
C
C      CALLS:       
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--********************************************************************

      SUBROUTINE GOTBRANCH(NCLUS,NOW,IEQUIV,NEQUIV,NEQMAX,IRTFLG)

 

      COMMON /UNITS/LUN,NIN,NOUT

      INTEGER  IEQUIV(2,*)

      NOWMIN = MIN(NCLUS,NOW)
      NOWMAX = MAX(NCLUS,NOW)

      DO I = 1,NEQUIV
 
C       CHECK IF ALREADY HAVE RECORDED THIS Y IN NEQUIV
        IF (IEQUIV(2,I) .EQ. NOWMIN .AND.
     &      IEQUIV(1,I) .EQ. NOWMAX) RETURN
      END DO

c**********************************
C**      WRITE(NOUT,*) ' NOW,NCLUS: ',NOW,NCLUS
c**********************************

      NEQUIV = NEQUIV + 1
      IERR   = 0
      IF (NEQUIV .GT. NEQMAX) THEN
C        RAN OUT OF SPACE IN IEQUIV ARRAY
         IF (IERR .LE. 1) THEN
            WRITE(NOUT,*) '*** BRANCH LIMIT (NEQMAX) IS:', NEQMAX
            WRITE(NOUT,*) ' OUTPUT IS JUNK !!!!!!!!!!!!!!!'
            IERR = IERR + 1
         ENDIF             
         IRTFLG = 1
         RETURN
      ENDIF

      IEQUIV(1,NEQUIV) = NOWMAX
      IEQUIV(2,NEQUIV) = NOWMIN

      IRTFLG = 0

      RETURN
      END
