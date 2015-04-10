
C ++********************************************************************
C                                                                      *
C GETTHREADS                                                           *
C               EXTRACTED FROM OTHER ROUTINES JAN 2000 ARDEAN LEITH
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
C **********************************************************************
C  GETTHREADS(NUMTH)
C
C  PURPOSE:  RETURNS NUMBER JOF CURRENT TRHEASD FOR OMP PROGRAMS
C
C  PARAMETERS:  NUMTH   NUMBER OF THREADS                      (RET.)
C
C***********************************************************************

       SUBROUTINE GETTHREADS(NUMTH)

       INTEGER  OMP_GET_NUM_THREADS

#ifdef SP_MP

c$omp parallel private(np)
	np = OMP_GET_NUM_THREADS()
c$omp single
	NUMTH = np
c$omp end single
c$omp end parallel

#else
	NUMTH = 1
#endif

       RETURN
       END
