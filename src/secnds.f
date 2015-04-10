 
 

C++*********************************************************************
C
C SECNDS.FOR  -- CREATED JAN 85
C
C **********************************************************************
C *  AUTHOR: R. BANERJEE                                                   *
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
C    SECNDS(TIME)
C
C    PURPOSE:     TO FIND ELAPSED TIME SINCE LAST CALL TO SECNDS
C
C    PARAMETERS:  TIME      A DUMMY VALUE ON UNIX
C
C--********************************************************************

#ifdef SP_IBMSP3

	REAL FUNCTION SECNDS(FTIME)

C       CLOCK RESETS AT ICMAX!! THIS WILL FAIL TO BE USEFULL !!!!!
        CALL SYSTEM_CLOCK(ICOUNT,ICPSEC,ICMAX)

        SECNDS = (FLOAT(ICOUNT) / FLOAT(ICPSEC)) 

        RETURN
        END

#else
	FUNCTION SECNDS(TIME)

C       I DO NOT KNOW IF SAVE IS NEEDED FEB 99 al
        SAVE

	REAL  TIME,SECNDS,dtime,T2(2)

        SECNDS = dtime(T2)
 
	RETURN
	END
#endif
