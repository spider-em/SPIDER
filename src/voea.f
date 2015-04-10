
C++*********************************************************************
C                                                                       
C VOEA.F                                         04/23/92                
C         MODIFIED TO GET EVEN SPACING IN PHI.   10/23/96                
C         USED REG_SET_NSEL                      AUG 00   ARDEAN LEITH   
C         USED LUNDOC                            JAN 14   ARDEAN LEITH   
C         FIO EQUIV REMOVED, CHECKED OUTPUT      MAR 14   ARDEAN LEITH   
C         LEGACY DEFAULT FOR T1,P1..             APR 14   ARDEAN LEITH                                                     
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
C  PURPOSE: PREPARE ANG DOC FILE FOR 'PJ 3Q'                              
C           WITH QUASI-EQUAL ANGULAR SPACING                            
C                                                                       
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE VOEA(IDUM)

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC'

         REAL                   :: FIO(3)

	 LOGICAL                :: SKIP
	 REAL, PARAMETER        :: QUADPI = 3.14159265358979323846
	 REAL, PARAMETER        :: DGR_TO_RAD = (QUADPI/180)
         CHARACTER (LEN=MAXNAM) :: DOCNAM
         LOGICAL                :: NEWFILE
         CHARACTER (LEN=80)     :: COMMENT

         INTEGER, PARAMETER     :: NDOC = 80

         CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

         LITER = 0

         PSI   = 0.0

         DELTA = 15.0
         CALL  RDPRM1S(DELTA,NOT_USED,'DELTA THETA',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         T1 =  0.0
         T2 = 90.0
         CALL RDPRM2S(T1,T2,NOT_USED,'RANGE OF THETA (0,90)',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         IF (T1 == 0.0 .AND. T2 == 0.0) THEN
C           FOR LEGACY DEFAULT
            T1 =  0.0
            T2 = 90.0
         ENDIF

         P1 =   0.0
         P2 = 359.9
         CALL  RDPRM2S(P1,P2,NOT_USED,'RANGE OF PHI (0,359.9)',IRTFLG)
         IF (P1 == 0.0 .AND. P2 == 0.0) THEN
C           FOR LEGACY DEFAULT
            P1 =   0.0
            P2 = 359.9
         ENDIF

         IF (IRTFLG .NE. 0) RETURN

         SKIP  = T1   <  90.0 .AND.
     &           T2  ==  90.0 .AND.
     &           P1  ==   0.0 .AND.
     &           P2   > 180.0

         CALL OPENDOC(DOCNAM,.TRUE.,NLET,NDOC,NDOCT,.TRUE.,
     &           'ANGULAR DOCUMENT',.FALSE.,.TRUE.,.TRUE.,
     &            NEWFILE,IRTFLG)
         IF (IRTFLG  ==  -1) GOTO 9999

         COMMENT = '          PSI,         THETA,          PHI'
         CALL LUNDOCPUTCOM(NDOCT,COMMENT,IRTFLG)

         LITER = 0

         DO THETA=T1,T2,DELTA
            IF (THETA == 0.0 .OR. THETA ==  180.0)  THEN
		DETPHI = 360.0
		LT     = 1
            ELSE
		DETPHI = DELTA / SIN(THETA*DGR_TO_RAD)
		LT     = MAX(INT((P2-P1)/DETPHI)-1,1)
		DETPHI = (P2-P1) / LT
            ENDIF

!           DO     PHI=P1,P2,DETPHI
 	    DO I=1,LT
	       PHI = P1+(I-1) * DETPHI
	       IF (SKIP .AND. THETA == 90.0 .AND. PHI >= 180.0) CYCLE
               LITER  = LITER + 1
               FIO(1) = PSI
               FIO(2) = THETA
               FIO(3) = PHI

               CALL LUNDOCWRTDAT(NDOCT,LITER,FIO,3,IRTFLG)

            ENDDO
         ENDDO

9999     IF (MYPID <= 0) THEN
            WRITE(NOUT,90)  LITER
90          FORMAT('  NUMBER OF DIRECTIONS: ',I0)

            CLOSE(NDOC)

            CALL FLUSHRESULTS

            CALL REG_SET_NSEL(1,1,FLOAT(LITER),0.0,0.0,0.0,0.0,IRTFLG)
         ENDIF

         END
