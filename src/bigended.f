
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


C       UNUSED BELOW FOR LEGACY REFERENCE ONLY FROM MRC -------------

#ifdef NEVER

C ---------------------- SETSTAMP -----------------------------------
 
      SUBROUTINE SETSTAMP(MACHSTMP,ISSWABT)

C     PURPOSE: SETS MACHINE STAMP FOR THIS ARCHITECTURE
C  
C     NOTE: I HAVE EXTRACTED THIS FROM THE MRC 2000 CODE AND
C           CONVERTED TO FORTRAN. BUT I MAY HAVE BOTCHED IT? al
C  
C     Bytes 213 and 214 contain 4 `nibbles' (half-bytes) indicating 
C     the representation of float, complex, integer and character 
C     datatypes. Bytes 215 and 216 are unused. The CCP4 library contains 
C     a general representation of datatypes, but in practice it is 
C     safe to use 0x44 0x44 0x00 0x00 for little endian machines, and 
C     0x11 0x11 0x00 0x00 for big endian machines. The CCP4 library 
C     uses this information to automatically byte-swap data if 
C     appropriate, when tranferring data files between machines.  

      INTEGER * 4 :: MACHSTMP
      LOGICAL     :: ISSWABT

      INCLUDE 'CMBLOCK.INC'

c       Little-ended Intel and some MIPS
#if defined(MIPSEL) || defined(i386) || defined(i860)
#  define NATIVEIT 4
#  define NATIVEFT 4
#endif

C       AN attempt at machines using the powerPC chip.             
#if defined (SP_PPC)
#  define NATIVEIT 1
#  define NATIVEFT 1
#endif

C       Compaq alpha's running Unix.             
#ifdef __alpha
#  define NATIVEFT 4
#  define NATIVEIT 4
#endif

C        SGI Altix running GNU/Linux.             
#ifdef __ia64
#  define NATIVEFT 4
#  define NATIVEIT 4
#endif

C       Big-endian ieee includes SGI machines,       
C       HP (68k-based or RISC), RS/6000 and all        
C       Suns except obsolete i386-based ones.       

#if defined(MIPSEB) || defined(__hpux) || defined(_AIX) || defined(m68k) || defined(mc68000) || defined(sparc)
#  define NATIVEIT 1
#  define NATIVEFT 1
#endif

C       Little-endian AMD OPTERON,       
#if defined  __x86_64__ || defined __ORDER_LITTLE_ENDIAN__
#  define NATIVEIT 4
#  define NATIVEFT 4
#endif


#if defined(SP_IBMSP3) 
#  define NATIVEIT 1
#  define NATIVEFT 1
#endif


#ifndef NATIVEFT
#error "Can't determine machine number format"
#endif

      NATIVEFTT = NATIVEFT
      NATIVEITT = NATIVEIT

C     write(6,*) '  ISSWABT:   ',isswabt
C     write(6,*) '  NATIVEFTT: ',nativeftt
C     write(6,*) '  NATIVEITT: ',nativeitt
      write(6,*) '  MACHSTMP:  ',machstmp
      IF (ISSWABT) THEN
C        RUNNING WITH NON-NATIVE BYTE-SWAPPING
         IF (NATIVEFTT == 1) THEN
            MACHSTMP = 4369
         ELSE
            MACHSTMP = 286326784
         ENDIF
      ELSE
         IF (NATIVEFTT == 1) THEN
            MACHSTMP = 286326784
         ELSE
            MACHSTMP = 4369
         ENDIF
       ENDIF
#endif



 
