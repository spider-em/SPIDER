
C ++********************************************************************
C                                                                      *
C SETTHREADS                                                           *
C               EXTRACTED FROM OTHER ROUTINES MAR 2000 ARDEAN LEITH
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
C  SETTHREADS(NUMTH)
C
C  PURPOSE:  SETS NUMBER JOF CURRENT TRHEADS FOR OMP PROGRAMS
C
C  PARAMETERS:  NUMTH   NUMBER OF THREADS                      (SENT.)
C
C***********************************************************************

       SUBROUTINE SETTHREADS(NUMTH)


#if defined (SP_MP) && defined (SP_IBMSP3)
C      IBM SP HAS PROBLEMS WITH THIS AT NERSC

       CALL omp_set_num_threads(NUMTH)

#else
#ifdef SP_MP

       CALL omp_set_num_threads(NUMTH)

#endif
#endif

       RETURN
       END
