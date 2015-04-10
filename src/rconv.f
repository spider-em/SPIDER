C++*********************************************************************
C
C RCONV.F               LONG FILENAMES           JAN 89 al
C                       USED OPFILE              NOV 00 ARDEAN LEITH
C                       OPFILEC                  FEB 03 ARDEAN LEITH
C                       REFACTORED               JAN 14 ARDEAN LEITH
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
C RCONV(LUN1,LUN2,LUNP,NSAM,NROW,NSLICE,MODE,MBUF)
C        LUN1               LOGICAL UNIT NUMBER OF FILE
C        LUN2               LOGICAL UNIT NUMBER OF FILE
C        LUNP               LOGICAL UNIT NUMBER OF PSF FILE
C        NSAM,NROW,NSLICE   DIMENSIONS OF FILE
C        MODE
C        MBUF
C
C--*******************************************************************

      SUBROUTINE RCONV(LUN1,LUN2,LUNP,NSAM,NROW,NSLICE,MODE,MBUF)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC' 
 
      REAL,ALLOCATABLE       :: PSF(:,:), Q(:)  

      CHARACTER(LEN=MAXNAM)  :: FILEP
      DOUBLE PRECISION       :: A
      CHARACTER(LEN=1)       :: NULL = CHAR(0)
C                                          123456789  
      CHARACTER (LEN=MAXNAM) :: PROMPTL = 'ROW('
      LOGICAL                :: SAVEPSF

C     OPTIONS WITH MODE > 1 RESERVED FOR STANDARD PSFS
C     TO BE IMPLEMENTED LATER
      IF (MODE .NE. 1) THEN
	 CALL ERRT(102,' OPTION NOT YET IMPLEMENTED',MODE)
         RETURN
      ENDIF

      CALL FILERD(FILEP,NLET,NULL,'PSF INPUT',IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      SAVEPSF = (FILEP(1:1) == '#')

      IF (NSLICE == 1)  THEN
C        IMAGE INPUT

C        -----------------------------------------------------------
C        NEW PSF, INPUT FROM CONSOLE OR INPUT STREAM EXPECTED,
C        PSF NOT TO BE SAVED

         IF (FILEP(1:1) .NE. '*' .AND. 
     &       FILEP(1:1) .NE. '#') THEN
C            READ PSF FROM FILE --------------------------------

            MAXIM = 0
            CALL OPFILEC(0,.FALSE.,FILEP,LUNP,'O',
     &                   IFORM,NSPRD,NSPRD,NDUM,
     &                   MAXIM,' ',.TRUE.,IRTFLG)
            IF (IRTFLG .NE. 0) RETURN

C           MEMORY AVAILABLE ?
	    K_Q    = 1
	    K_B    = K_Q + NSPRD*NSAM
            MUBUF  = K_B + NSAM

            ALLOCATE (PSF(NSPRD,NSPRD),
     &                Q(MUBUF), STAT=IRTFLG)
            IF (IRTFLG .NE. 0) THEN
               CALL ERRT(46,'RCONV; PSF..',MUBUF*NSPRD*NSPRD)
               GOTO 999
            ENDIF

            CALL REDVOL(LUNP,NSPRD,NSPRD,1,1,PSF,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 999

         ELSE

C           PSF INPUT FROM CONSOLE OR INPUT STREAM 

            IF (SAVEPSF) THEN
C              PSF TO BE SAVED
               CALL FILERD(FILEP,NLETO,NULL,'PSF OUTPUT',IRTFLG)
               IF (IRTFLG .NE. 0) GOTO 999
            ENDIF

            CALL RDPRI1S(NSPRD,NOT_USED,'PSF WIDTH',IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 999

C           MEMORY AVAILABLE ?
	    K_Q    = 1
	    K_B    = K_Q + NSPRD*NSAM
            MUBUF  = K_B + NSAM

            ALLOCATE (PSF(NSPRD,NSPRD),
     &                Q(MUBUF), STAT=IRTFLG)
            IF (IRTFLG .NE. 0) THEN
               CALL ERRT(46,'RCONV; PSF..',MUBUF*NSPRD*NSPRD)
            ENDIF

            IF (SAVEPSF) THEN
C              OPEN FILE FOR PSF OUTPUT 
               IFORM = 1
               MAXIM = 0
               CALL OPFILEC(0,.FALSE.,FILEP,LUNP,'U',
     &                   IFORM,NSPRD,NSPRD,1,
     &                   MAXIM,' ',.TRUE.,IRTFLG)
               IF (IRTFLG .NE. 0) GOTO 999
            ENDIF

            NSPRD1 = NSPRD  / 2 + 1
            NSPRD2 = NSPRD1 + 1
            NSPRDS = NSPRD  * NSPRD

            WRITE(NOUT,43) NSPRD, NSPRD
            IF (NDAT .NE. 0 .AND. NDAT .NE. NOUT)
     &         WRITE(NDAT,43)NSPRD, NSPRD
43          FORMAT(/,'  PSF MATRIX(',I2,' x ',I2,')')

            DO  IROW = 1,NSPRD

               CALL INTTOCHAR(IROW,PROMPTL(5:),NC,1)
               PROMPTL(5+NC:5+NC) = ')'

               CALL RDPRA(PROMPTL(1:5+NC), NSPRD,0,.FALSE.,
     &                 PSF(1,IROW),NGOT,IRTFLG)

               IF (IRTFLG .NE. 0 ) THEN
                 CALL ERRT(101,'INPUT CAN NOT BE UNDERSTOOD',NE)
                 GOTO 999
               ELSEIF (NGOT < NSPRD) THEN
                 CALL ERRT(102,'INSUFFICIENT VALUES FOR THIS LINE',NGOT)
                 GOTO 999
               ENDIF
            ENDDO

C           NORMALIZE IF POSSIBLE
            A = SUM(PSF)
            IF (A == 0.0) THEN
               A = 1.0
               WRITE(NOUT,*) ' *** PSF CANNOT BE NORMALIZED'
            ELSE
               WRITE(NOUT,*) ' PSF NORMALIZED, DIVIDED BY:',A
            ENDIF
   
            PSF = PSF / A     ! ARRAY OPERATION

            IF (SAVEPSF) CALL WRTVOL(LUNP,NSPRD,NSPRD, 1,1,PSF,IRTFLG)
            CLOSE(LUNP)

         ENDIF

C        WRITE OUT PSF FOR VERIFICATION
         WRITE(NOUT,*) ' '
         WRITE(NOUT,*) ' PSF IN USE:'
         DO  IROW = 1,NSPRD
            WRITE(NOUT,90) (PSF(ICOL,IROW),ICOL=1,NSPRD)
 90	    FORMAT(10(1X,ES10.3))
         ENDDO
         WRITE(NOUT,*) ' '

C        REAL 2D, CONVOLUTION
         NQ = NSPRD / 2
         CALL RCNV2_P(LUN1,LUN2,Q(K_B),Q(K_Q),NSAM,NROW,PSF,NQ)

C------------------------------------------------------------

      ELSE
C        3-D REAL CONVOLUTION WITH PSF READ FROM FILE
         MAXIM = 0
         CALL OPFILEC(0,.FALSE.,FILEP,LUNP,'O',IFORM,
     &                NSPRD,NSPRD,NSPRD,
     &                MAXIM,' ',.TRUE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         MUBUF = NSAM*NROW*NSPRD+NSAM+NSPRD**3
         ALLOCATE (Q(MUBUF), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'RCONV; Q..',MUBUF*NSPRD*NSPRD)
            GOTO 999
         ENDIF
 
         DO K=1,NSPRD
            DO I = 1,NSPRD
	       NR =(K-1)*NSPRD+I
               I1 = (K-1)*NSPRD*NSPRD+(I-1)*NSPRD + 1
               CALL REDLIN(LUNP,Q(I1),NSPRD,NR)
            ENDDO
         ENDDO
         CLOSE(LUNP)

C        WRITE OUT PSF FOR VERIFICATION
         WRITE(NOUT,*) ' '
         WRITE(NOUT,*) ' PSF IN USE:'

	 K1 = 1
	 K2 = NSPRD
	 DO K=1,NSPRD
	    WRITE(NOUT,*) K

	    DO I=1,NSPRD

	       WRITE(NOUT,892) (Q(J),J=K1,K2)
892	       FORMAT(10(1X,ES10.3))

	       K1 = K1 + NSPRD
	       K2 = K2 + NSPRD
	    ENDDO
	 ENDDO

	 K1 = NSPRD*NSPRD*NSPRD

	 CALL RCNV3_P(LUN1,LUN2,
     &	    Q(NSAM*NROW*NSPRD+K1+1),Q(K1+1),NSAM,NROW,NSLICE,Q,NSPRD)
      ENDIF
	
999   IF (ALLOCATED(PSF)) DEALLOCATE(PSF)
      IF (ALLOCATED(Q))   DEALLOCATE(Q)

      END
