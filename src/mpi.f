
C ++********************************************************************
C
C MPI.F   MPI SUPPORT                               OCT 08 ARDEAN LEITH
C         ADDED SET_MPI                             MAY 09 ARDEAN LEITH
C         LONG_LONG _INT                            DEC 12 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012 Health Research Inc.,                         *
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
C  CONTAINS:
C    BCAST_MPI(OPER,VARNAME,VAR,ILEN,CMPITYPE,ICOMM)
C    SEND_MPI (OPER,VARNAME,VAR,ILEN,CMPITYPE,IDEST,ITAG,ICOMM)
C    RECV_MPI (OPER,VARNAME,VAR,ILEN,CMPITYPE,ISOURCE,ITAG,ICOMM)
C    ALLREDUCE_MPI(OPER,VARNAME,VARIN,VAROUT,ILEN,CMPITYPE,CMPIHOW,ICOMM)
C    SET_MPI(ICOMM,MYPID,IRTFLG)
C
C  PURPOSE:  MY MPI SUPPORT ROUTINES
C            
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************


C **********************************************************************
C
C  BCAST_MPI(OPER,VARNAME, VAR,ILEN, CMPITYPE,ICOMM)
C
C  PURPOSE:  BRROADCASTS MSG IN MPI
C            
C  PARAMETERS:  OPER        OPERATION                             (SENT)
C               VARNAME     VARIABLE NAME                         (SENT)
C               VAR         VARIABLE POITER                       (SENT)
C               ILEN        VARIABLE LENGTH                       (SENT)
C               MPITYPE     VARIABLE TYPE (I,F,C,L,D..)           (SENT)
C               ICOMM       MPI SERIES                            (SENT)
C                                                     
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE BCAST_MPI(OPER,VARNAME, VAR,ILEN, CMPITYPE,ICOMM)

#ifdef USE_MPI
        include 'mpif.h'

        CHARACTER(LEN=*), INTENT(IN) :: OPER
        CHARACTER(LEN=*), INTENT(IN) :: VARNAME
        REAL, INTENT(IN)             :: VAR      ! ACTUALLY A POINTER
        INTEGER, INTENT(IN)          :: ILEN,ICOMM
        CHARACTER(LEN=1), INTENT(IN) :: CMPITYPE
        INTEGER                      :: MPIERR,ILENT,NLETO,NLETV
        INTEGER                      :: MPITYPE

        CALL TYPE_MPI(CMPITYPE,MPITYPE)

        ILENT = MAX(ILEN,1)
        CALL MPI_BCAST(VAR,ILENT,MPITYPE, 0,ICOMM,MPIERR)

        IF (MPIERR .NE. 0) THEN
           NLETO = lnblnkn(OPER)
           NLETV = lnblnkn(VARNAME)
           WRITE(0,*)OPER(:NLETO),'; FAILED TO BCAST: ',VARNAME(:NLETV)
           STOP 'FAILED IN BCAST_MPI'
        ENDIF
#else
C       DUMMY ROUTINE, FOR NON-MPI USE
        RETURN
#endif
        END

C **********************************************************************
C
C  RECV_MPI(OPER,VARNAME, VAR, ILEN, CMPITYPE, ISOURCE,ITAG, ICOMM)
C
C  PURPOSE:  RECEIVES MSG IN MPI
C            
C  PARAMETERS:  OPER        OPERATION                             (SENT)
C               VARNAME     VARIABLE NAME                         (SENT)
C               VAR         VARIABLE POITER                       (SENT)
C               ILEN        VARIABLE LENGTH                       (SENT)
C               CMPITYPE    VARIABLE TYPE (I,F,C,L,D..)           (SENT)
C               ISOURCE     SOURCE OF MESSAGE                     (SENT)
C               ITAG        MSG. TAG                              (SENT)
C               ICOMM       MPI SERIES                            (SENT)
C                                                     
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE RECV_MPI(OPER,VARNAME, VAR, ILEN, CMPITYPE,
     &                      ISOURCE,ITAG, ICOMM)

#ifdef USE_MPI
        include 'mpif.h'

        CHARACTER(LEN=*), INTENT(IN) :: OPER
        CHARACTER(LEN=*), INTENT(IN) :: VARNAME
        REAL, INTENT(IN)             :: VAR       ! ACTUALLY A POINTER
        INTEGER, INTENT(IN)          :: ILEN,ISOURCE,ITAG,ICOMM
        CHARACTER(LEN=1), INTENT(IN) :: CMPITYPE

        INTEGER                      :: MPIERR,ILENT,NLETO,NLETV,ITAGT
        INTEGER                      :: MPITYPE

        CALL TYPE_MPI(CMPITYPE,MPITYPE)

        ILENT = MAX(ILEN,1)

        ITAGT = ITAG
        IF (ITAGT .LE. 0) ITAGT = MPI_ANY_TAG  

        CALL MPI_RECV(VAR,ILENT,MPITYPE, ISOURCE,ITAGT,
     &                ICOMM, MPI_STATUS_IGNORE, MPIERR)

        IF (MPIERR .NE. 0) THEN
           NLETO = lnblnkn(OPER)
           NLETV = lnblnkn(VARNAME)
           WRITE(0,*)OPER(:NLETO),'; FAILED TO RECV: ',VARNAME(:NLETV)
           STOP
        ENDIF
#else
C       DUMMY ROUTINE, FOR NON-MPI USE
        RETURN
#endif
        END

C **********************************************************************
C
C  SEND_MPI(OPER,VARNAME, VAR, ILEN, CMPITYPE, ISOURCE,ITAG, ICOMM)
C
C  PURPOSE:  SENDS MSG IN MPI
C            
C  PARAMETERS:  OPER        OPERATION                             (SENT)
C               VARNAME     VARIABLE NAME                         (SENT)
C               VAR         VARIABLE POITER                       (SENT)
C               ILEN        VARIABLE LENGTH                       (SENT)
C               CMPITYPE    VARIABLE TYPE (I,F,C,L,D..)           (SENT)
C               IDEST       DESTINATION FOR MESSAGE               (SENT)
C               ITAG        MSG. TAG                              (SENT)
C               ICOMM       MPI SERIES                            (SENT)
C                                                     
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE SEND_MPI(OPER,VARNAME, VAR, ILEN, CMPITYPE,
     &                      IDEST,ITAG, ICOMM)

#ifdef USE_MPI
        include 'mpif.h'

        CHARACTER(LEN=*), INTENT(IN) :: OPER
        CHARACTER(LEN=*), INTENT(IN) :: VARNAME
        REAL, INTENT(IN)             :: VAR    ! ACTUALLY A POINTER
        INTEGER, INTENT(IN)          :: ILEN,ICOMM,IDEST,ITAG
        CHARACTER(LEN=1), INTENT(IN) :: CMPITYPE

        INTEGER                      :: MPIERR,ILENT,NLETO,NLETV
        INTEGER                      :: MPITYPE

        CALL TYPE_MPI(CMPITYPE,MPITYPE)

        ILENT = MAX(ILEN,1)
         
        CALL MPI_SEND(VAR,ILENT,MPITYPE, IDEST,ITAG, ICOMM,MPIERR)
        IF (MPIERR .NE. 0) THEN
           NLETO = lnblnkn(OPER)
           NLETV = lnblnkn(VARNAME)
           WRITE(0,*)OPER(:NLETO),'; FAILED TO SEND: ',VARNAME(:NLETV)
           STOP
        ENDIF
#else
C       DUMMY ROUTINE, FOR NON-MPI USE
        RETURN
#endif
        END


C *************************** ALLREDUCE_MPI ****************************
C
C  ALLREDUCE_MPI(OPER,VARNAME, VARIN,VAROUT,ILEN,CMPITYPE,CMPIHOW,ICOMM)
C
C  PURPOSE:  REDUCES VARIN TO: VAROUT 
C            
C  PARAMETERS:  OPER        OPERATION                             (SENT)
C               VARNAME     VARIABLE NAME                         (SENT)
C               VARIN       VARIABLE POITER                       (SENT)
C               VAROUT      VARIABLE POITER                       (RET.)
C               ILEN        VARIABLE LENGTH                       (SENT)
C               CMPITYPE    VARIABLE TYPE (I,F,C,L,D..)           (SENT)
C               CMPITYPE    REDUCTION TYPE (S..)                  (SENT)
C               ITAG        MSG. TAG                              (SENT)
C               ICOMM       MPI SERIES                            (SENT)
C                                                     
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE ALLREDUCE_MPI(OPER,VARNAME,VARIN,VAROUT,
     &                           ILEN, CMPITYPE,CMPIHOW,ICOMM)

#ifdef USE_MPI
        include 'mpif.h'

        CHARACTER(LEN=*), INTENT(IN) :: OPER
        CHARACTER(LEN=*), INTENT(IN) :: VARNAME
        REAL, INTENT(IN)             :: VARIN     ! ACTUALLY A POINTER
        REAL, INTENT(OUT)            :: VAROUT    ! ACTUALLY A POINTER
        INTEGER, INTENT(IN)          :: ILEN,ICOMM
        CHARACTER(LEN=1), INTENT(IN) :: CMPITYPE
        CHARACTER(LEN=1), INTENT(IN) :: CMPIHOW
        INTEGER                      :: MPIERR,ILENT,NLETO,NLETV
        INTEGER                      :: MPITYPE

        CALL TYPE_MPI(CMPITYPE,MPITYPE)

        IF (CMPIHOW == 'S') THEN
            IHOW = MPI_SUM
        ELSE
            CALL ERRT(101,'UNKNOWN REDUCTION METHOD',NE)
            STOP
        ENDIF

        ILENT = MAX(ILEN,1)
        
        CALL MPI_ALLREDUCE(VARIN,VAROUT,ILENT, MPITYPE,IHOW,
     &                     ICOMM,MPIERR)
        IF (MPIERR .NE. 0) THEN
           NLETO = lnblnkn(OPER)
           NLETV = lnblnkn(VARNAME)
           WRITE(0,*)OPER(:NLETO),'; FAILED TO REDUCE: ',VARNAME(:NLETV)
           STOP
        ENDIF
#else
C       DUMMY ROUTINE, FOR NON-MPI USE
        RETURN
#endif
        END



C       ------------------------ TYPE_MPI -----------------------

        SUBROUTINE TYPE_MPI(CMPITYPE,MPITYPE)

C       CONVERTS CHAR. MPI TYPE TO MPI VARIABLE TYPE SPECIFICATION
#ifdef USE_MPI
        include 'mpif.h'

        CHARACTER(LEN=1), INTENT(IN) :: CMPITYPE
        INTEGER, INTENT(OUT)         ::  MPITYPE
 
        IF     (CMPITYPE == 'C') THEN
            MPITYPE = MPI_CHARACTER
        ELSEIF (CMPITYPE == 'I') THEN
            MPITYPE = MPI_INTEGER
        ELSEIF (CMPITYPE == 'R') THEN
            MPITYPE = MPI_REAL
        ELSEIF (CMPITYPE == 'D') THEN
            MPITYPE = MPI_DOUBLE_PRECISION
        ELSEIF (CMPITYPE == 'L') THEN
            MPITYPE = MPI_LOGICAL
        ELSEIF (CMPITYPE == 'X') THEN
            MPITYPE = MPI_COMPLEX
        ELSE
            CALL ERRT(101,'UNKNOWN MPI VARIABLE TYPE',NE)
            STOP
        ENDIF

#else
C       DUMMY ROUTINE, FOR NON-MPI USE
        RETURN
#endif
        END

C ********************************* SET_MPI **************************
C
C   SET_MPI(ICOMM,MYPID,IRTFLG)    
C
C   PURPOSE: SETS ICOMM AND MYMPID
C
C   PARAMETERS:  ICOMM    MPI ICOMM                             (RET.)
C                MYPID    MPI PROCESS ID                        (RET.)
C                           -1 IF NOT USING MPI
C                IRTFLG   UNUSED ERROR FLAG                     (RET.)
C                            0  IS NORMAL
C
C--*********************************************************************

      SUBROUTINE SET_MPI(ICOMM,MYPID,IRTFLG)

#ifdef USE_MPI
      INCLUDE 'mpif.h'

      ICOMM = MPI_COMM_WORLD
      CALL MPI_COMM_RANK(ICOMM, MYPID, IERR)
      IRTFLG = IERR
#else
      MYPID  = -1
      IRTFLG =  0
#endif

      END






