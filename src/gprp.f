C++*********************************************************************
C
C GPRP.F          REMOVED FROM UTIL@           ARDEAN LEITH  9/05/03    
C                 SETPRMB PARAMETERS           ARDEAN LEITH  5/19/09
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
C  GPRP(LUN,NSAM,NROW,NSLICE,FCHART)
C
C  PURPOSE:         GET OR RPLACE SINGLE PIXEL FROM IMAGE/VOLUME
C
C  PARAMETERS:      FCHART         FCHAR
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE GPRP(LUN,NSAM,NROW,NSLICE,FCHART)

C       USE INLINE BUFFER COMMON AREA
        INCLUDE 'INLN_INFO.INC'

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'

        COMMON /IOBUF/ BUF(NBUFSIZ)

        CHARACTER(LEN=*) :: FCHART
        CALL SET_MPI(ICOMM,MYPID,MPIERR)  ! SET MYPID

        IF (NSLICE .LE. 1) THEN
	   CALL RDPRIS(IX,IY,NOT_USED,'COLUMN & ROW',IRTFLG)
           IZ = 1
        ELSE
           IZ = 0
           CALL RDPRI3S(IX,IY,IZ,NOT_USED,
     &                    'COLUMN, ROW, & SLICE',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
           IF (IZ .LE. 0) THEN
               CALL RDPRI1S(IZ,NOT_USED,'SLICE',IRTFLG)
           ENDIF
        ENDIF
        IF (IRTFLG .NE. 0) RETURN

        IF (IX .LT. 1 .OR. IX .GT. NSAM .OR.
     &      IY .LT. 1 .OR. IY .GT. NROW .OR.
     &      IZ .LT. 1 .OR. IZ .GT. NSLICE) THEN
           CALL ERRT(101,'OUTSIDE IMAGE/VOLUME BOUNDS',NE)
           IRTFLG = 1
           RETURN
        ENDIF


        IF (FCHART(1:2) .EQ. 'GP') THEN 
C          OPERATION  GP -- (GET PIXEL VALUE)--------------------- 'GP' 

           IF (ISINLINE(LUN)) THEN
C             USE INLINED BUFFER FOR I/O (SEE OPENINLN.F)
              CALL INLN_REDVOX(LUN,NROW,BVAL,1,
     &                         IX,IY,IZ,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN
           ELSE
	      IREC = (IZ-1) * NROW + IY
 	      CALL REDLIN(LUN,BUF,NSAM,IREC)
              BVAL = BUF(IX)
           ENDIF

	   IF (MYPID .LE. 0) WRITE(NOUT,90) IX,IY,IZ,BVAL
90         FORMAT('  LOCATION: ',3I7,'  VALUE: ',1PE12.5)

           CALL REG_SET_NSEL(1,1,BVAL, 0.0, 0.0, 
     &                       0.0, 0.0, IRTFLG)

        ELSE
C          OPERATION  RP ---(REPLACE PIXEL) --------------------- 'RP' 

	   CALL RDPRM(BVAL,NOT_USED,'NEW PIXEL VALUE')

           IF (ISINLINE(LUN)) THEN
C             USE INLINED BUFFER FOR I/O (SEE OPENINLN.F)
              CALL INLN_WRTVOX(LUN,NROW,BVAL,1,
     &                         IX,IY,IZ,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

           ELSE
	      IREC = (IZ-1) * NROW + IY
	      CALL REDLIN(LUN,BUF,NSAM,IREC)
              BUF(IX) = BVAL
 	      CALL WRTLIN(LUN,BUF,NSAM,IREC)
           ENDIF

C          SET UNDETERMINED STATISTICS FLAG
           IF (IFORM.GT.0) CALL SETPRMB(LUN, 0.0,0.0, 0.0,0.0)
        ENDIF

        IRTFLG = 0
        RETURN
        END



