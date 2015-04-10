
C++*********************************************************************
C
C REVERSEBYTES 
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
C    REVERSEBYTES(CLINE,NCHAR,IRTFLG) 
C
C    PURPOSE: REVERSE BYTES IN A CHARACTER STRING      
C
C    PARAMETERS:  CSTRING  STRING                  (SENT & RETURNED) 
C                 NSAM     STRING LENGTH *4        (SENT)
C                 IRTFLG   ERROR FLAG (O = OK)     (RETURNED)    
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

	SUBROUTINE  REVERSEBYTES(CLINE,NCHART4,IRTFLG)

        CHARACTER *(*) CLINE
        CHARACTER      TDAT

C       NOTE THAT THE ALPHA-NUMERICAL DATA (ABCD) WILL BE WRITTEN 
C       IN INVERTED ORDER (DCBA)
        DO K = 1, NCHART4, 4
           TDAT           = CLINE(K:K)
           CLINE(K:K)     = CLINE(K+3:K+3)
           CLINE(K+3:K+3) = TDAT
           TDAT           = CLINE(K+1:K+1)
           CLINE(K+1:K+1) = CLINE(K+2:K+2)
           CLINE(K+2:K+2) = TDAT
        ENDDO

        IRTFLG = 0

        RETURN
        END
