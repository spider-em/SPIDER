C++*********************************************************************
C
C  AF.F                               
C                 GETNEWSTACK PARAM.             FEB 03   ARDEAN LEITH
C                 SETPRMB PARAMETERS             MAY 09   ARDEAN LEITH 
C                 GETNEWSTACK PARAM.             OCT 10   ARDEAN LEITH
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
C   AF(MAXDIM,LUN1,LUN2,LUN3,FILNAM,FILNAMO,NSAM,NROW,NSLICE,
C         MAXIM,IRTFLG)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE AF(MAXDIM,LUN1,LUN2,LUN3,FILNAM,FILNAMO,
     &                 NSAM,NROW,NSLICE,MAXIM,IRTFLG)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      COMMON         BUF(1)

      INTEGER               :: ILIST(6)
      CHARACTER(LEN=*)      :: FILNAM,FILNAMO
      CHARACTER(LEN=MAXNAM) :: DOCNAM
      CHARACTER(LEN=1)      :: NULL

      INCLUDE 'F90ALLOC.INC'
      REAL, DIMENSION(:,:), POINTER :: IPQ
      REAL, DIMENSION(:,:), POINTER :: QBUF

      NULL   = CHAR(0)
      IRTFLG = 1

C     FIND MEMORY NEEED 
      K1     = 1 + NSAM * NROW
      MEMTOT = K1 + NSAM

      IF (MEMTOT .GT. MAXDIM)  THEN
C        OVERFLOW WILL OCCUR
         CALL  ERRT(6,'AF',NE)
         RETURN
      ENDIF

      IF (MAXIM .LT. 0) THEN
C        SINGLE IMAGE OPERATION
         AAA = 1.0
         BBB = 0.0
         CALL RDPRM2S(AAA,BBB,NOT_USED,
     &             'TRANSFORMATION PARAMETERS A,B',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         CCC = 0.0
         DDD = 1.0
         CALL RDPRM2S(CCC,DDD,NOT_USED,
     &             'TRANSFORMATION PARAMETERS C,D',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

	 SHXI = 0.0
	 SHYI = 0.0
         CALL  RDPRM2S(SHXI,SHYI,NOT_USED,
     &                'SHIFTS IN X AND Y',IRTLFG)
         IF (IRTFLG .NE. 0) RETURN

C        CHECK THAT DETERMINANT OF THE TRANSFORMATION IS NOT TOO SMALL
	 DET = AAA * DDD - BBB * CCC
	 IF (ABS(DET) .LT. 1.0E-4)  THEN
	    CALL  ERRT(34,'AF ',NE)
	    RETURN
	 ENDIF

C        LOAD AND TRANSFORM SLICE BY SLICE
         DO L=1,NSLICE
           LB = (L-1)*NROW
           DO J=1,NROW
              CALL  REDLIN(LUN1,BUF(1+(J-1)*NSAM),NSAM,LB+J)
           ENDDO

C          ROTATE THIS SLICE
           CALL AFS(BUF,BUF(K1),NSAM,NROW,
     &		     AAA,BBB,CCC,DDD,SHXI,SHYI,LUN2,LB)
         ENDDO

C        RESET FILE HEADER FOR ALTERATIONS IN STATISTICS
         CALL SETPRMB(LUN2, 0.0,0.0, 0.0,0.0)
         IRTFLG = 0

      ELSE
C        WHOLE STACK OPERATION
         NUMB = 6
         CALL RDPRAI(ILIST,6,NUMB,1,6,
     &      'REG. NUMBERS FOR A,B,C,D, X, & Y SHIFT',
     &      'T',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

C        MAXX IS 1 + NUM OF REGISTERS SINCE PBUB CONTAINS KEY ALSO
         MAXX   = MAX(ILIST(1),ILIST(2),ILIST(3),
     &                ILIST(4),ILIST(5),ILIST(6)) + 1
         MAXY   = 0
         CALL GETDOCDAT('TRANSFORMATION PARAMETERS DOCUMENT',
     &                 .TRUE.,DOCNAM,LUN3,.TRUE.,MAXX, MAXY,IPQ,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
         QBUF = > IPQ

         IMGNUM = 1
         DO WHILE (IMGNUM .LE. MAXIM)

C           GET INPUT IMAGE FROM STACK            
            CALL GETOLDSTACK(LUN1,NSAM,IMGNUM,
     &                       .TRUE.,.FALSE.,.TRUE.,IRTFLG)

            IF (IRTFLG .EQ. 0) THEN
C              CREATE OUTPUT IMAGE IN STACK
               CALL GETNEWSTACK(LUN1,LUN2,.FALSE.,NSAM,IMGNUM,IRTFLG)

               IF (IRTFLG .EQ. 0) THEN
C                 OUTPUT IMAGE CREATED OK
                  AAA  = QBUF( ILIST(1) + 1, IMGNUM )
                  BBB  = QBUF( ILIST(2) + 1, IMGNUM )
                  CCC  = QBUF( ILIST(3) + 1, IMGNUM )
                  DDD  = QBUF( ILIST(4) + 1, IMGNUM )
                  SHXI = QBUF( ILIST(5) + 1, IMGNUM )
                  SHYI = QBUF( ILIST(6) + 1, IMGNUM )

C                 CHECK THAT DETERMINANT OF TRANSFORMATION IS NOT TOO SMALL
	          DET = AAA*DDD-BBB*CCC
	          IF (ABS(DET) .LT. 1.0E-4)  THEN
	             CALL  ERRT(34,'AF ',NE)
	             RETURN
	          ENDIF
C                 LOAD AND ROTATE SLICE BY SLICE
                  DO L=1,NSLICE
                     LB = (L-1)*NROW
                     DO J=1,NROW
                        CALL  REDLIN(LUN1,BUF(1+(J-1)*NSAM),NSAM,LB+J)
                     ENDDO
              
                     CALL  AFS(BUF,BUF(K1),NSAM,NROW,
     &		        AAA,BBB,CCC,DDD,SHXI,SHYI,LUN2,LB)
                  ENDDO

                  WRITE(NOUT,90) AAA,BBB,CCC,DDD,SHXI,SHYI 
90                FORMAT(' A: ',G10.3,'  B: ',G10.3,'  C: ',G10.3,
     &                   '  D:',G10.3,/,
     &                   ' SHIFTS(',G10.3,',',G10.3,')',/)

C                 RESET HEADER FOR ALTERATIONS IN STATISTICS
                  CALL SETPRMB(LUN2, 0.0,0.0, 0.0,0.0)
               ENDIF
            ENDIF
            IMGNUM = IMGNUM + 1
         ENDDO

C        DEALLOCAT DOC. FILE MEMORY
         DEALLOCATE(IPQ)
      ENDIF

      RETURN
      END

