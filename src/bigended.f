
C ++********************************************************************
C                                                                      *
C  BIGENDED                                                          *            *
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
C **********************************************************************                                                                     *
C 
C  PURPOSE:  CHECK ENDEDNESS OF CURRENT ARCHITECTURE. 
C            RETURNS TRUE  IF BIG ENDIAN (IBM, SUN, SGI)
C            RETURNS FALSE IF LITTLE ENDIAN (VAX, COMPAQ, INTEL)   
C                                                                      *
C  PARAMETERS:  IDUM  UNUSED                                                       *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	
         LOGICAL FUNCTION BIGENDED(IDUM)

        INTEGER * 2    ITWO
        INTEGER * 1    I1ARRAY(2)
        EQUIVALENCE (ITWO,I1ARRAY(1))

C       GET CURRENT ARCHITECTURE ENDED-NESS
        I1ARRAY(1) = 0
        I1ARRAY(2) = 0

        ITWO       = 1

        BIGENDED   = (I1ARRAY(1) .EQ. 0) 
        END


 
