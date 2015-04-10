
C **********************************************************************
C  
C UTIL5.F  AUTHOR: M.RADERMACHER
C          REMOVED 'MR FILTER' AT M.R.'s REQUEST  OCT 2005 ArDean Leith
C          REMOVED CPUTRMN ON GNU COMPILER        APR 2012 ArDean Leith
C          REMOVED SUPPORT                        SEP 2014 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C OPERATIONS: 2DN, COMPLN, FTN, SUM3N, 3DN, MAKE3, FTINVN, INV2N, INV3N
C
C NOTE:
C   RADON transform routines are out-of-date, (no longer same as MR) and
C   some do not compile with gfortran compiler so I have moved them to
C   Attic and no longer support them in SPIDER.  Please contact 
C   M. Radermacher  at Univ. Vermont for improved versions.
C
C   mv rad2calln.f       RCS/rad2calln.f*     Attic  
C   mv radcompleten.f    RCS/radcompleten.f*  Attic
C   mv fouradcalln.f     RCS/fouradcalln.f*   Attic
C   mv cputrmn.f         RCS/cputrmn.f*       Attic
C   mv radon3dn.f        RCS/radon3dn.f*      Attic
C   mv mrrmake.f         RCS/mrrmake.f*       Attic
C   mv fouradinvcn.f     RCS/fouradinvcn.f*   Attic
C   mv callrminvn.f      RCS/callrminvn.f*    Attic
C   mv call3drminvn.f    RCS/call3drminvn.f*  Attic
C   mv rmreset.f         RCS/rmreset.f*       Attic
C   mv radon2en.f        RCS/radon2en.f*      Attic
C   mv ccopyangles_1.f   RCS/ccopyangles_1.f* Attic
C   mv c1gfilt.f         RCS/c1gfilt.f*       Attic
C   mv fouradn.f         RCS/fouradn.f*       Attic
C   mv symangs.f         RCS/symangs.f*       Attic
C   mv putrmrealn.f      RCS/putrmrealn.f*    Attic
C   mv putrmn.f          RCS/putrmn.f*        Attic
C   mv wrtarray.f        RCS/wrtarray.f*      Attic
C   mv rminvn.f          RCS/rminvn.f*        Attic
C   mv mrfour1.f         RCS/mrfour1.f*       Attic
C   mv fouradinvn.f      RCS/fouradinvn.f*    Attic
C   mv mrfft.f           RCS/mrfft.f*         Attic
C
C **********************************************************************

        SUBROUTINE UTIL5(IDUM)
      
        CALL ERRT(101,
     &      '<RM> OPERATIONS NOT SUPPORTED IN THIS SPIDER VERSION', NE)

        END
