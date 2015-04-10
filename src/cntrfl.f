
C++*******************************************************************
C
C CNTRFL.F     TAKEN FROM [SPIDER.FORTRAN]CNTRFL2.FOR MAY 27 88 al
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
C    CNTRFL(LUN1,LUN2,NSAM,NROW,NSLICE,MAXDIM)
C
C    PURPOSE:   OUTPUT CONTUR FILE
C
C    PARAMETERS:    LUN1         LOGICAL UNIT NUMBER
C                   LUN2         LOGICAL UNIT NUMBER
C                   NSAM,NROW    DIMENSIONS OF FILE
C                   MAXDIM       MAXIMUM BUFFER SPACE AVAILABLE
C
C    CALLED BY: UTIL3
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

	SUBROUTINE CNTRFL(LUN1,LUN2,NSAM,NROW,NSLICE,MAXDIM)

 

        INCLUDE 'CMBLOCK.INC'
 
	COMMON        DARRAY(1024)

	INTEGER       PRESPT,PASTPT,NEXTPT
	LOGICAL       ONTOP,ADJ,LEVSET
        CHARACTER     IANS,IANS1,NULL

        DATA          FLTZER/ 10E-7/

	NULL = CHAR(0)

c**************************debug
C        FVAL = 2.0
C        WRITE(7,*) ' ------ FVAL:',FVAL, '  ----------------------'
C        DO IB = -2,2
C          BOTT = IB
C          CRIB = BOTT - FVAL
C          WRITE(7,*) ' '
C          WRITE(7,*) ' ------ BOTT:',BOTT, ' -- CRIT:',CRIB,' ----'
C          DO I = -8,6
C             DF   = I
C             DIA  = DF - AMOD(DF - BOTT, FVAL)
C             D2 = DIA
C             IF (DF .LT. BOTT) D2 = CRIB
C             WRITE(7,900)  DF,DIA,D2
C900          FORMAT(3F6.0)
C          ENDDO
C        ENDDO
C        FMIN = -6.0
C        FMAX = 7.0
C        NCON = 2
C        FVAL = (FMAX-FMIN) / FLOAT(NCON)
C        BOTT = FMIN 
C        FMM  = FMIN - FVAL
C        WRITE(7,*) ' FMIN,FMAX,FVAL:',FMIN,FMAX,FVAL, '  -------'
C
C        WRITE(7,*) ' '
C        WRITE(7,*) ' ------ BOTT:',BOTT, ' -- FMM:',FMM,' ----'
C
C        DO I = FMIN,FMAX
C             DF   = I
C             DIA  = DF - AMOD(DF - BOTT, FVAL)
C             D2   = DIA
C             IF (DF .EQ. FMIN) D2 = FMM
C             WRITE(7,900)  DF,DIA,D2
C900          FORMAT(3F6.0)
C        ENDDO
CC********************

	IF (5*NSAM .GT. MAXDIM)THEN
C          INSUFFICIENT BUFFER SPACE
           CALL ERRT(6,'CNTRFL',NE)
           RETURN
        ENDIF

        IF (IMAMI.NE.1) CALL NORM3(LUN1,NSAM,NROW,NSLICE,FMAX,FMIN,AV)

        WRITE(NOUT,90) FMIN,FMAX
90      FORMAT(' IMAGE RANGE: ', 1PG11.3,'...',1PG11.3)


8       IF (FCHAR(4:4) .EQ. 'E') THEN
C          SPECIFIED CONTOUR LEVELS
10         CALL RDPRM1S(BOTT,NOT_USED,'BOTTOM LEVEL',IRTFLG)
           IF (IRTFLG .EQ. -1) RETURN

           CALL RDPRM1S(FVAL,NOT_USED,'CONTOUR LEVEL SEPERATION',IRTFLG)
           IF (IRTFLG .EQ. -1) GOTO 10

           LEVSET = .TRUE.
           CRIT   = BOTT
           PVAL   = BOTT - FVAL

        ELSE
C          RANGE FROM FMIN...FMAX WILL HAVE NO-1 CONTOUR LEVELS
           CALL RDPRI1S(NCON,NOT_USED,'NO. OF CONTOUR LEVELS',IRTFLG)
           IF (IRTFLG .EQ. -1) RETURN

           IF (NCON .LE. 2) THEN
              CALL ERRT(30,'CNTRFL',NE)
              RETURN
           ENDIF

           NCON   = NCON - 1
           LEVSET = .FALSE.

           FVAL   = (FMAX-FMIN) / FLOAT(NCON)
           BOTT   = FMIN 
           FMM    = FMIN - FVAL
        ENDIF


	CALL RDPRMC(IANS,NCHAR,.TRUE.,'OVERWRITE? (Y/N)',NULL,IRT)
	IF (IRT .EQ. -1) GOTO 8
 
        ONTOP = .FALSE.
C       SET CONTOUR INTENSITY LEVELS FOR NOT OVERLAY
        FM1   = 0.0
        FM2   = 2.0
        ADJ   = .FALSE.

	IF (IANS .NE. 'N' .AND. IANS .NE. 'n') THEN
C          OVERLAY CONTOURS ON ORIGINAL IMAGE
           ONTOP = .TRUE.
C          SET INTENSITY LEVELS FOR OVERLAY CONTOUR
           FM1 = FMIN
           FM2 = FMAX
	   CALL RDPRMC(IANS1,NCHAR,.TRUE.,
     &       '(W)HITE, (B)LACK OR (A)ADJUSTED  CONTOURS',NULL,IRT)

           HF = FMIN + (FMAX - FMIN) / 2.0
           IF (IANS1 .EQ. 'A') ADJ = .TRUE.
           IF (IANS1 .NE. 'W') THEN
C             INVERT CONTOUR INTENSITY LEVELS
              FM1 = FMAX
              FM2 = FMIN
           ENDIF
        ENDIF


C       CALCULATES ROW 1 -------------------------------------------------

	NSET = NSAM + 1
	IP   = 3 * NSAM
	IQ   = IP + NSAM
	IF (.NOT. ONTOP) IQ = IP
	CALL REDLIN(LUN1,DARRAY,NSAM,1)
	CALL REDLIN(LUN1,DARRAY(NSET),NSAM,2)

	DO IA=1,NSAM
           NDUM = NSAM + IA
           DF   = DARRAY(IA)
           DIA  = DF - AMOD(DF - BOTT, FVAL)

           IF (ABS(DIA) .LT. FLTZER) DIA = 0.0

           IF (ADJ) THEN
              IF (DIA .GT. HF) THEN
                 FM2 = FMIN
              ELSE
                 FM2 = FMAX
              ENDIF
           ENDIF

           IF (LEVSET) THEN
              IF (DF .LT. CRIT)  DIA = PVAL
           ELSE
              IF (DF .EQ. FMIN)  DIA = FMM
           ENDIF
           DARRAY(IA) = DIA


           DF    = DARRAY(NDUM)
           DNDUM = DF - AMOD(DF - BOTT,FVAL)
           IF (ABS(DNDUM) .LT. FLTZER) DNDUM = 0.0
	   IF (LEVSET) THEN
              IF (DF .LT. CRIT) DNDUM = PVAL
           ELSE
              IF (DF .EQ. FMIN) DNDUM = FMM
           ENDIF
           DARRAY(NDUM) = DNDUM
        ENDDO

C       WRITES ROW 1 TO OUTPUT FILE

	IF (ONTOP)CALL REDLIN(LUN1,DARRAY(IQ+1),NSAM,1)
	DARRAY(IP+1) = FM1

	IF ((DARRAY(1) - DARRAY(2)      .GT. FLTZER) .OR.
     &      (DARRAY(1) - DARRAY(NSAM+1) .GT. FLTZER))
     &       DARRAY(IQ+1) = FM2

	DO  IA=2,NSAM-1
           DARRAY(IP+IA) = FM1
           DIA           = DARRAY(IA)
           IF ((DIA - DARRAY(IA+1)    .GT. FLTZER) .OR.
     &         (DIA - DARRAY(IA-1)    .GT. FLTZER) .OR. 
     &         (DIA - DARRAY(NSAM+IA) .GT. FLTZER)) 
     &          DARRAY(IQ+IA) = FM2
        ENDDO

	DARRAY(IP+NSAM) = FM1
	IF ((DARRAY(NSAM) - DARRAY(NSAM-1) .GT. FLTZER).OR.
     &      (DARRAY(NSAM) - DARRAY(2*NSAM) .GT. FLTZER)) 
     &       DARRAY(IQ+NSAM) = FM2
	CALL WRTLIN(LUN2,DARRAY(IQ+1),NSAM,1)

C       CALCULATES ROWS 2 TO NROW -1 ------------------------------------

	PRESPT = NSAM
	NEXTPT = 2*NSAM
	PASTPT = 0

	DO I=3,NROW
          NSET = NSET+NSAM
          IF (NSET .GT. 3*NSAM) NSET = 1
          CALL REDLIN(LUN1,DARRAY(NSET),NSAM,I)

          IF (ONTOP) CALL REDLIN(LUN1,DARRAY(IQ+1),NSAM,I-1)

          DF    = DARRAY(NSET)
          DNSET = DF - AMOD(DF - BOTT,FVAL)
          IF (ABS(DNSET) .LT. FLTZER) DNSET = 0.0

          IF (LEVSET) THEN
             IF (DF .LT. CRIT) DNSET = PVAL
          ELSE
             IF (DF .EQ. FMIN) DNSET = FMM
          ENDIF

          DARRAY(NSET) = DNSET

          NDUM = NSET+1
          DF    = DARRAY(NDUM)
          DNDUM = DF - AMOD(DF - BOTT,FVAL)
          IF (ABS(DNDUM) .LT. FLTZER) DNDUM = 0.0

          IF (LEVSET) THEN
            IF (DF .LT. CRIT) DNDUM = PVAL
          ELSE
            IF (DF .EQ. FMIN) DNDUM = FMM
          ENDIF
          DARRAY(NDUM) = DNDUM

          PASTPT = PASTPT + 1
          PRESPT = PRESPT + 1
          NEXTPT = NEXTPT + 1
          PASTPT = MOD(PASTPT, 3*NSAM)
          PRESPT = MOD(PRESPT, 3*NSAM)
          NEXTPT = MOD(NEXTPT, 3*NSAM)
          DARRAY(IP+1) = FM1

          DPRESPT = DARRAY(PRESPT)
          IF (ADJ) THEN
            IF (DPRESPT .GT. HF) THEN
               FM2 = FMIN
            ELSE
               FM2 = FMAX
            ENDIF
          ENDIF
          IF ((DPRESPT - DARRAY(PRESPT+1) .GT. FLTZER) .OR.
     &        (DPRESPT - DARRAY(PASTPT)   .GT. FLTZER) .OR.
     &        (DPRESPT - DARRAY(NEXTPT)   .GT. FLTZER))  
     &         DARRAY(IQ+1) = FM2

          DO IA=2,NSAM-1
            NDUM  = NSET+IA
            DF    = DARRAY(NDUM)
            DNDUM = DF - AMOD(DF - BOTT,FVAL)
            IF (ABS(DNDUM) .LT. FLTZER) DNDUM = 0.0

            IF (LEVSET) THEN
              IF (DF .LE. CRIT) DNDUM = PVAL
            ELSE
              IF (DF .EQ. FMIN) DNDUM = FMM
            ENDIF
            DARRAY(NDUM) = DNDUM

            PASTPT = PASTPT+1
            PRESPT = PRESPT+1
            NEXTPT = NEXTPT+1
            DARRAY(IP+IA) = FM1

            DPRESPT = DARRAY(PRESPT)
            IF (ADJ) THEN
              IF (DPRESPT .GT. HF) THEN
                 FM2 = FMIN
              ELSE
                 FM2 = FMAX
              ENDIF
            ENDIF

            IF ((DPRESPT - DARRAY(NEXTPT)   .GT. FLTZER) .OR.
     &          (DPRESPT - DARRAY(PASTPT)   .GT. FLTZER) .OR.
     &          (DPRESPT - DARRAY(PRESPT-1) .GT. FLTZER) .OR.
     &          (DPRESPT - DARRAY(PRESPT+1) .GT. FLTZER)) THEN
                    DARRAY(IQ+IA) = FM2
            ENDIF

          ENDDO

          DARRAY(IP+NSAM) = FM1
          PRESPT = PRESPT+1
          PASTPT = PASTPT+1
          NEXTPT = NEXTPT+1
          DPRESPT = DARRAY(PRESPT)

          IF (ADJ) THEN
            IF (DPRESPT .GT. HF) THEN
               FM2 = FMIN
            ELSE
               FM2 = FMAX
            ENDIF
          ENDIF

          IF ((DPRESPT - DARRAY(PASTPT)   .GT. FLTZER) .OR.
     &        (DPRESPT - DARRAY(NEXTPT)   .GT. FLTZER) .OR.
     &        (DPRESPT - DARRAY(PRESPT-1) .GT. FLTZER)) 
     &        DARRAY(IQ+NSAM) = FM2
          CALL WRTLIN(LUN2,DARRAY(IQ+1),NSAM,I-1)

        ENDDO


C       WRITES LAST ROW TO OUTPUT FILE ------------------------------------

	IF (ONTOP)CALL REDLIN(LUN1,DARRAY(IQ+1),NSAM,NROW)
	PRESPT = PRESPT+1
	PASTPT = PASTPT+1
	PASTPT = MOD(PASTPT,3*NSAM)
	PRESPT = MOD(PRESPT,3*NSAM)
	DARRAY(IP+1) = FM1

        DPRESPT = DARRAY(PRESPT)
	IF ((DPRESPT - DARRAY(PRESPT+1) .GT. FLTZER) .OR.
     &      (DPRESPT - DARRAY(PASTPT)   .GT. FLTZER))   
     &       DARRAY(IQ+1) = FM2

	DO IA=2,NSAM-1
          DARRAY(IP+IA) = FM1
          PRESPT  = PRESPT + 1
          PASTPT  = PASTPT + 1
          DPRESPT = DARRAY(PRESPT)
          IF (ADJ) THEN
            IF (DPRESPT .GT. HF) THEN
              FM2 = FMIN
            ELSE
              FM2 = FMAX
            ENDIF
          ENDIF

          IF ((DPRESPT - DARRAY(PRESPT+1) .GT. FLTZER) .OR.
     &        (DPRESPT - DARRAY(PRESPT-1) .GT. FLTZER) .OR.
     &        (DPRESPT - DARRAY(PASTPT)   .GT. FLTZER))   
     &         DARRAY(IQ+IA) = FM2
        ENDDO

	PRESPT = PRESPT + 1
	PASTPT = PASTPT + 1
        DPRESPT = DARRAY(PRESPT)
	IF (ADJ) THEN
          IF (DPRESPT .GT. HF) THEN
             FM2 = FMIN
          ELSE
             FM2 = FMAX
          ENDIF
        ENDIF

        DARRAY(IP+NSAM) = FM1
	IF ((DPRESPT - DARRAY(PASTPT)   .GT. FLTZER) .OR.
     &      (DPRESPT - DARRAY(PRESPT-1) .GT. FLTZER)) 
     &       DARRAY(IQ+NSAM) = FM2
	CALL WRTLIN(LUN2,DARRAY(IQ+1),NSAM,NROW)

	RETURN
	END
