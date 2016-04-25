
C ++********************************************************************
C                                                                      *
C FLUSHRESULTS.F    NEW                         NOV 2000 ARDEAN LEITH
C                   CHANGED VERBOSE             JUN 2004 ARDEAN LEITH                                                           *
C                   defined (SP_GFORTRAN)       DEC 2005 ARDEAN LEITH
C                   OUTPUT MSG                  NOV 2013 ARDEAN LEITH
C                                                                   *
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
C                                                                      *
C  FLUSHRESULTS
C
C  PURPOSE:  CLOSES AND REOPENS RESULTS FILE TO ENSURE IT IS UPDATED 
C            IN CASE OF CRASH
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE FLUSHRESULTS

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'

        INTEGER  :: ICOMM,MYPID,MPIERR,IRET

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

#if defined(SP_GFORTRAN) 
        IF (NOUT .NE. NDAT) CALL flush(NOUT)
        CALL flush(NDAT)
#else
#if defined (SP_IBMSP3) 
        IF (NOUT .NE. NDAT) CALL flush_(NOUT)
        CALL flush_(NDAT)
#else
        IF (NOUT .NE. NDAT) CALL flush(NOUT,IRET)
        CALL flush(NDAT,IRET)
#endif
#endif

        IF (MYPID <= 0) 
     &     CALL PDATES(' RESULTS FILE FLUSHED: ',-1)

        END



        SUBROUTINE FLUSHRESULTS_Q(SAYIT)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'

        LOGICAL  :: SAYIT

        INTEGER  :: ICOMM,MYPID,MPIERR,IRET

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

#if defined(SP_GFORTRAN) 
        IF (NOUT .NE. NDAT) CALL flush(NOUT)
        CALL flush(NDAT)
#else
#if defined (SP_IBMSP3) 
        IF (NOUT .NE. NDAT) CALL flush_(NOUT)
        CALL flush_(NDAT)
#else
        IF (NOUT .NE. NDAT) CALL flush(NOUT,IRET)
        CALL flush(NDAT,IRET)
#endif
#endif

        IF (SAYIT .AND. MYPID <= 0) 
     &     CALL PDATES(' Results file flushed: ',-1)

        END
