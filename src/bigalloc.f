
C++*********************************************************************
C
C BIGALLOC.F   -- NEW JAN 2001 ARDEAN LEITH
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
C    BIGALLOC(IASK8,IOK,SAYIT,CALLERRT,IRTFLG)
C
C    PURPOSE: COMPLAIN ABOUT POSSIBLE EXCESSIVE ALLOCATION OF 
C             MEMORY
C
C    PARAMETERS:   IASK8           MEMORY WANTED                 (SENT)
C                  IOK             MAX. MEMORY OK                (RET.)
C                  SAYIT           ECHO COMPILATION TYPE         (SENT)
C                  CALLERRT        CALL ERRT IF NOT OK           (SENT)
C                  IRTFLG          ZERO IS OK                    (RET.)
C
C--*********************************************************************

        SUBROUTINE BIGALLOC(IASK8,IOK,SAYIT,CALLERRT,IRTFLG)

#ifndef SP_32
        INTEGER *8       IASK8,IOK
#else
        INTEGER *4       IASK8,IOK
#endif
        INTEGER * 4 :: IASK4
        LOGICAL     :: SAYIT,CALLERRT

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        IRTFLG = 0
        IOK    = HUGE(IASK8)

#ifdef SP_32
        IF (SAYIT) WRITE(NOUT,*) ' COMPILED FOR NT 32'
        RETURN
#endif

#if _MIPS_SZPTR==64
C       COMPILED ON 64 BIT SGI MACHINE
        IF (SAYIT) WRITE(NOUT,*) ' COMPILED FOR -64'
        RETURN
#endif

#ifdef __osf
C       COMPILED ON 64 BIT COMPAQ??
        IF (SAYIT) WRITE(NOUT,*) ' COMPILED FOR Alpha 64'
        RETURN
#endif

#ifdef SP_IBMSP3
C       COMPILED ON 64 BIT IBM??
        IF (SAYIT) WRITE(NOUT,*) ' COMPILED FOR IBM 64'
        RETURN
#endif

        IASK4  = HUGE(IASK4)
        IOK    = HUGE(IASK4) 
        IRTFLG = 1

#if _MIPS_SZPTR==32 
C       COMPILED ON 32 BIT SGI MACHINES
        IF (SAYIT) WRITE(NOUT,*) ' COMPILED FOR -n32'

        IF (IASK8 .GT. IASK4) THEN
           WRITE(NOUT,90) IASK4,IASK8
90         FORMAT(' *** WARNING -n32 ALLOCATION ON SGI LIMITED TO: ',I11,
     &               ' YOU ARE USING >',I13)
           IF (CALLERRT) THEN
              CALL ERRT(101,'EXCESSIVE ALLOCATION REQUEST',NDUM)
           ENDIF
        ENDIF
#else
C       COMPILED ON UNKNOWN PROBABLY 32 BIT PLATFORM
        IF (IASK8 .GT. IASK4) THEN
           WRITE(NOUT,91) IASK4,IASK8
91         FORMAT('*** ALLOCATION MAY LIMITED TO: ',I11,
     &               ' YOU ARE USING > ',I13)
        ENDIF
#endif

        RETURN
	END
