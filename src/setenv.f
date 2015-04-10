
C++*********************************************************************
C
C SETENV   NEW APR 2001 ArDean Leith
C          NT HAS DIFFERENT FUNCTION NAME        Sep 2001 ArDean Leith
C          GNU COMPILER  HAS NO FUNCTION         Apr 2012 ArDean Leith
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
C   SETENV(VAR,VALUE,IRTFLG)
C
C   PURPOSE:   SETS AN ENVIRONMENTAL VARIABLE FROM SPIDER
C
C   PARAMETERS:   VAR       NAME OF ENVIRONMENTAL VARIABLE   (SENT)
C                 VALUE     VALUE OF ENVIRONMENTAL VARIABLE  (SENT)
C                 ITRFLG    ERROR FLAG   (0 IS NORMAL)       (RETURNED)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE SETENV(VAR,VALUE,IRTFLG)

	INCLUDE 'CMBLOCK.INC' 

        CHARACTER(LEN=*)   :: VAR,VALUE
        CHARACTER(LEN=160) :: STRING 

#if defined (SP_GFORTRAN)
        WRITE(NOUT,*) 'OPERATION NOT AVAILABLE WITH GNU COMPILER'
        IRTFLG = 1
        RETURN
#else
        IF ((LEN(VAR) + LEN(VALUE)) .GE. 160) THEN
           CALL ERRT(101,'STRING OVERFLOW IN SETENV',NDUM)
           IRTFLG = 1
           RETURN
        ENDIF
 
        STRING = VAR // '=' // VALUE
        
C       SET THE ENVIRONMENT VARIABLE
#if defined (SP_NT) || defined (SP_IFC) || defined(__INTEL_COMPILER) 
        CALL setenvqq(STRING)
#else 
        CALL putenv(STRING)
#endif

        IRTFLG = 0
#endif

        END

