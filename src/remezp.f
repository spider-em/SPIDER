
C ++********************************************************************
C                                                          
C  REMEZP                                                              
C                  USED OPFILE                NOV 00 ARDEAN LEITH
C                  USED ALLOCATE              JAN 01 ARDEAN LEITH 
C                  OPFILEC                    FEB 03 ARDEAN LEITH
C                  MAXNAM                     JUL 14 ARDEAN LEITH
C                                                              
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
C REMEZP(MAXDIM)                                                       
C                                          
C PURPOSE: 	PROGRAM TO DESIGN THE REMEZ FILTER
C                                                                   
C PARAMETERS:                                                         
C                                                                      
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE REMEZP(MAXDIM)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON A(80),BUF(1024),ROT(29,29,29),BQ(1)

        REAL, ALLOCATABLE, DIMENSION(:) :: Q

C       CHANGES IN COMMON LENGTH SHOULD BE CONSULTED WITH
C       SUBROUTINES CONTAINING ARRAY ROT:
C       ROT_P, ROT3_P, INF_P, INF3_P, INFD3_P, OPSF_P, OPSF3_P

	DOUBLE PRECISION        ROT
	DIMENSION               H(66)
        CHARACTER(LEN=MAXNAM):: FILNAM
	CHARACTER*1             NULL,WHAT

	DATA  NFREQ/128/,IDIM/33/,LUN/9/

        NULL = CHAR(0)

C       LENRO - MEMORY RESERVED IN COMMON BY ROT AND OTHERS
	LENRO=80+1024+2*29**3 

111     CALL RDPRMI(NFILT,IDUM,NOT_USED,'IMPULSE RESPONSE LENGTH')

	NFILT=2*(NFILT/2)+1
	NFCNS=NFILT/2+1
	IF (NFCNS .GT. 29)  GOTO 111
	CALL  RMZ_P(H,NFILT)

        CALL RDPRMC(WHAT,NCHAR,.TRUE.,
     &	      'PLOT FREQUENCY RESPONSE (Y/N)',NULL,IRTFLG)

	IF (WHAT .EQ. 'Y')  THEN
           DO    J=1,NFREQ
              BUF(J)=0.0
           ENDDO
           BUF(1)=H(NFCNS)
           KK=NFCNS-1
           DO    J=KK,1,-1
              K=NFREQ-KK+J
              JJ=KK-J+2
              BUF(JJ)=H(J)
              BUF(K)=H(J)
           ENDDO
           NV=LOG2(NFREQ)
           CALL  FFTR_Q(BUF,NV)
C          50 SCALING FACTOR FOR PICTURE
           DELF=50
           NF2=NFREQ/2+1
           W=BUF(2)
           K=1
           DO J=3,NFREQ,2
              K=K+1
              BUF(K)=BUF(J)
           ENDDO
           BUF(NF2)=W
           LR=NFREQ/((IDIM-1)*2)
	   IF (LR.GT.1)  THEN
	      K=1
	      DO J=2,NF2-1,LR
	         K=K+1
	         Z=0.0
	         DO M=1,LR
 	            Z=Z+BUF(J+M-1)
	         ENDDO
 	         BUF(K)=Z/FLOAT(LR)
	      ENDDO
	      BUF(IDIM)=BUF(NF2)
	   ENDIF
	   DO J=1,IDIM
              BUF(J)=BUF(J)*DELF
	   ENDDO
	   BMA=BUF(1)
	   BMI=BMA
	   DO J=2,IDIM
	      BMA=AMAX1(BMA,BUF(J))
	      BMI=AMIN1(BMI,BUF(J))
	   ENDDO
	   T = (BMA-BMI)/50.0
	   DO J=1,IDIM
	      BUF(J)=(BUF(J)-BMI)/T
	   ENDDO
           WRITE(NOUT,606)
 606       FORMAT(//15X,'***** FREQUENCY RESPONSE *****')
           CALL MRKUR3(BUF,IDIM,0.,0,60)
        ENDIF

 222    CALL RDPRMC(WHAT,NCHAR,.TRUE.,
     &	    'DO YOU LIKE YOUR FILTER (Y/N)',NULL,IRTFLG)

	IF (WHAT .EQ. 'N') THEN
	   GOTO  111
	ELSE
	  CALL  FILERD(FILNAM,NLET,NULL,'FILTER',IRTFLG)
          IF (IRTFLG .NE. 0)  RETURN

333	  WRITE(NOUT,788)
788	  FORMAT(' DIMENSIONS OF RESULTING FILTER IN FOURIER SPACE'/
     &	     ' IF 0 THEN POINT SPREAD FUNCTION IS STORED'/
     &       ' TO BE USED FOR REAL SPACE FILTRATION WITH RC COMMAND') 
	  CALL  RDPRMI(NSAM,NROW,NOT_USED,'(NSAM,NROW)')
	  IF (NSAM .EQ. 0)  THEN
789	     CALL  RDPRMI(IPSF,NDUM,NOT_USED,
     &		  'PRODUCE 2-D OR 3-D PSF? (2/3)')
	     IF (IPSF.NE.2 .AND. IPSF.NE.3)  GOTO  789
C            PSF TO BE USED IN RC
	     IF (IPSF.EQ.3)  THEN
	        IFORM = 3
	        NFIL3 = NFILT
	     ELSE
	        IFORM = 1
	        NFIL3 = 1
	     ENDIF

             MAXIM = 0
             CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'U',IFORM,
     &                   NFILT,NFILT,NFIL3,
     &                   MAXIM,' ',.TRUE.,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

	      IF (IPSF.EQ.2)  THEN
	         CALL  ROT_P(H,NFCNS,ROT)
	         CALL  OPSF_P(LUN,ROT,NFCNS,BQ,NFILT)
	      ELSE
	         CALL  ROT3_P(H,NFCNS,ROT)
	         CALL  OPSF3_P(LUN,ROT,NFCNS,BQ,NFILT)
	      ENDIF
	      GOTO 9999
	   ENDIF
	
	   CALL  RDPRMI(NSLICE,NDUM,NOT_USED,
     &	       'NUMBER OF SLICES (IF 0 THEN 2-D FILTER), (NSLICE)')
	   NSLICE = MAX0(NSLICE,1)

	   IF (MOD(NSAM,2).EQ.0)  THEN
              IF (NSLICE.EQ.1)  THEN
                 IFORM=-12
              ELSE
	         IFORM=-22
              ENDIF
           ELSE
              IF (NSLICE.EQ.1)  THEN
                 IFORM=-11
              ELSE
                 IFORM=-21
              ENDIF
           ENDIF
           LRCL = NSAM+2-MOD(NSAM,2)

           MEMWANT = LRCL*NROW*NSLICE
           ALLOCATE(Q(MEMWANT),STAT=IRTFLG)
           IF (IRTFLG .NE. 0)  THEN
              CALL  ERRT(102,'UNABLE TO ALLOCATE Q',MEMWANT)
              GOTO 9999
	   ENDIF

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'U',IFORM,LRCL,NROW,NSLICE,
     &                   MAXIM,' ',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9998

           INV = +1
           IF (NSLICE.EQ.1)  THEN
	      CALL  ROT_P(H,NFCNS,ROT)
	      CALL  INF_Q(Q(1),LRCL,NSAM,NROW,ROT,NFCNS)
	      CALL  FMRS_2(Q(1),NSAM,NROW,INV)
           ELSE
	      CALL  ROT3_P(H,NFCNS,ROT)
              CALL  INF3_P(Q(1),LRCL,NSAM,NROW,NSLICE,ROT,NFCNS)
	      CALL  FMRS_3(Q(1),NSAM,NROW,NSLICE,INV)
           ENDIF
           DO K=1,NROW*NSLICE
              CALL WRTLIN(LUN,Q(1+(K-1)*LRCL),LRCL,K)
           ENDDO

9998       CLOSE(LUN)  
9999       IF (ALLOCATED(Q)) DEALLOCATE(Q)
           RETURN

        ENDIF

        END
