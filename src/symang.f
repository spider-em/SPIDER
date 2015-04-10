C ++********************************************************************
C                                                                      *
C SYMANG.F          NEW                       APR 2004  JAMIE LEBARRON * 
C                   ROTATION INTEGER BUG      OCT 2007  ARDEAN LEITH   *
C                   ORTH. NAN BUG             MAY 2010  ARDEAN LEITH   *
C                   GROEL                     NOV 2011  ARDEAN LEITH   *
C                                                                      *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C  SYMANG(CSYM,IFOLD,LUNDOC,IRTFLG)
C
C  PURPOSE: PRINT ANGLES OF SYMMETRY FOR DIFFERENT SYMMETRY GROUPS.
C           FOR THE TETRAHEDRAL AND ICOSAHEDRAL GROUPS, OUTPUT IS 
C           CHANGED SO THAT LARGEST AXIS OF SYMMETRY IS ALIGNED WITH Z
C           AXIS. 
C
C  VARIABLES:
C       X_OR_Z             = MATRIX CREATION FLAG, X=0 Z =1
C       DEGREE             = DEGREES OF TILT FOR MATRIX CREATION
C       DELTA#             = 360 DEGREES / #
C       COUNT#             = VARIOUS LOOP COUNTERS
C       T_OR_O             = SYMMETRY DIFFERENCE BETWEEN TETRA & OCTAL
C       CSYM               = COUNTER FOR WHAT KIND OF SYMMETRY
C       IFOLD              = LEVEL OF SYMMETRY (FIRST 2 SYMMETRY GROUPS)
C       ALPHA, BETA, GAMMA = SYMMETRY ANGLE PARAMETERS
C       CSYM               = WHAT TYPE OF SYMMETRY NEEDED
C       DLIST              = HOLDER FOR PSI, THETA, PHI
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE SYMANG(CSYM,IFOLD,LUNDOC,IRTFLG)
    
        IMPLICIT NONE

        CHARACTER(LEN=2)          :: CSYM
        INTEGER                   :: IFOLD, LUNDOC, IRTFLG

        INTEGER                   :: X_OR_Z, COUNT1 
	INTEGER                   :: COUNT2, T_OR_O, COUNT3, COUNT4
        INTEGER                   :: IKEY, IC

	REAL, PARAMETER           :: DELTA2 = 180.0
	REAL, PARAMETER           :: DELTA3 = 120
	REAL, PARAMETER           :: DELTA5 = 72

	DOUBLE PRECISION          :: DEGREE, TESTER, DELTAN
	DOUBLE PRECISION          :: PSI, THETA, PHI 
	DOUBLE PRECISION          :: DT 
  
       
        DOUBLE PRECISION, PARAMETER           :: ALPHA = 58.282524
	DOUBLE PRECISION, PARAMETER           :: BETA  = 20.905157
	DOUBLE PRECISION, PARAMETER           :: GAMMA = 54.735611
	DOUBLE PRECISION, DIMENSION (1:3,1:3) :: A, B, C, G, AB, N, P 
	DOUBLE PRECISION, DIMENSION (1:3,1:3) :: BA, D, E, F, H, M
 
        CHARACTER(LEN=60)                     :: COMMENT
        REAL                                  :: DLIST(3)

	DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626
  
	IKEY = 0
  
C -------------------------- ROTATIONAL --------------------------------	
C       TRY TO GET Dn group RUNNING
        IF ((CSYM(1:1) .EQ. 'C') .OR. (CSYM(1:1) .EQ. 'c')) THEN
          IF ((CSYM(2:2) .EQ. 'I') .OR. (CSYM(2:2) .EQ. 'i')) THEN
             COMMENT = 
     &         'ROTATIONAL SYMMETRY WITH INVERSION ::  PSI, THETA, PHI'
          ELSE
             COMMENT = 'ROTATIONAL SYMMETRY  ::  PSI, THETA, PHI'
          ENDIF
          CALL LUNDOCPUTCOM(LUNDOC,COMMENT,IRTFLG)
	
	  DELTAN = 360.0 / IFOLD
	  
C         START X LOOP
          DO COUNT1=0, 1
	    
	    DEGREE = COUNT1 * DELTA2
	    X_OR_Z = 0

C           X(J*delta2)	    	  	 
            CALL MATCREATE(X_OR_Z, DEGREE, A)
	    
C           START Z LOOP
            DO COUNT2=0, (IFOLD -1)
	    
	      DEGREE = COUNT2 * DELTAN
	      X_OR_Z = 1

C	      Z(I* DELTA N)
	      CALL MATCREATE(X_OR_Z, DEGREE, B)
	      
	      G = MATMUL(A,B)
	      CALL MATEXTRACT(G,PSI, THETA, PHI)
C	      PRINT "(f9.4,f9.4,f9.4)",PSI, THETA, PHI

              DLIST(1) = PSI
              DLIST(2) = THETA
              DLIST(3) = PHI
              IKEY     = IKEY + 1
              CALL LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,3,IRTFLG)
	    ENDDO
	    
	    IF ((CSYM(2:2) .NE. 'I').AND.(CSYM(2:2) .NE. 'i'))GOTO 9999
          ENDDO
	RETURN
	ENDIF
	
C --------------------------  TETRAHEDRAL & OCTOHEDRAL ----------------    

C       FOR 'T' TYPE AND 'O' TYPE, 
C       CHANGED SO THAT RUNS ONLY FOR OCT TYPE, LT.5 = FOR BOTH

        IF ((CSYM(1:1) .EQ. 'O') .OR.(CSYM(1:1) .EQ. 'o')) THEN

          COMMENT = 'OCTAHEDRAL/CUBIC SYMMETRY  ::  PSI, THETA, PHI'
          CALL LUNDOCPUTCOM(LUNDOC,COMMENT,IRTFLG)
	  
C	  ASSUME ITS T
          T_OR_O = 1

C         IF NOT CHANGE
          IF ((CSYM(1:1) .EQ. 'O').OR.(CSYM(1:1) .EQ. 'o')) T_OR_O = 3
	  
C         CREATE THE X(GAMMA)Z(45) MATRIX
          X_OR_Z = 0
	  DEGREE = GAMMA
	  CALL MATCREATE(X_OR_Z, DEGREE, A)
	  
	  X_OR_Z = 1
	  DEGREE = 45.0
	  CALL MATCREATE(X_OR_Z, DEGREE, B)
	  
	  AB = MATMUL(A,B)
	  BA = TRANSPOSE(AB)
	  
C         START X*K LOOP
          DO COUNT1=0,1
	  
	    X_OR_Z = 0
	    DEGREE = COUNT1 * DELTA2
	    CALL MATCREATE(X_OR_Z, DEGREE, C)
	    
C           START X(2 OR 4)
            DO COUNT2 =0, T_OR_O
	      
	      DELTAN= 360.0 /(T_OR_O + 1.0)
	      DEGREE = COUNT2 * DELTAN
	      X_OR_Z = 1
	      CALL MATCREATE(X_OR_Z, DEGREE, D)
	      
C             START P*I LOOP
              DO COUNT3 =0,2
	      
	        DEGREE = COUNT3*DELTA3
		X_OR_Z = 1
		CALL MATCREATE(X_OR_Z, DEGREE, E)
		
C               CALC ZI*XGAMMA*Z45
                F = MATMUL(E,AB)

C               CALC XGAMMA'*Z45'*ZI*XGAMMA*Z45 =PI
                E = MATMUL(BA,F)

C               CALC CHANGEABLEZ * PI
                F = MATMUL(D,E)

C               CALC XDELTA2*CHANGEABLEZ * PI
                G = MATMUL(C,F)

		CALL MATEXTRACT(G,PSI, THETA, PHI)
C               PRINT "(f9.4,f9.4,f9.4)",PSI, THETA, PHI

                DLIST(1) = PSI
                DLIST(2) = THETA
                DLIST(3) = PHI
                IKEY     = IKEY + 1
                CALL LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,3,IRTFLG)

	      ENDDO
	    ENDDO
	  ENDDO

	RETURN
	ENDIF

C -------------------------------- TETRAHEDRAL ----------------------	
C       TETRAHEDRAL, WITH 3AXIS ALIGN W/Z AXIS, POINT ON +VE X AXIS

        IF ((CSYM(1:1) .EQ. 'T') .OR. (CSYM(1:1) .EQ. 't'))THEN
	
          COMMENT = 'TETRAHEDRAL SYMMETRY  ::  PSI, THETA, PHI'
          CALL LUNDOCPUTCOM(LUNDOC,COMMENT,IRTFLG)

C         VALUE FROM (90 -DIHEDRAL ANGLE) + 90 =? 109.47
          TESTER = 1.0 / 3.0

          !DEGREE = 180.0 - (ACOS(TESTER) * (180.0/PI))
          DT      = MAX(-1.0D0, MIN(1.0D0,TESTER))
          DEGREE  = 180.0 - (DACOS(DT) * (180.0/PI))
	  
	  DO COUNT1=0,2
	    
	    PSI   = 0.0
	    THETA = 0.0
	    PHI   = COUNT1 * DELTA3
C           PRINT "(f9.4,f9.4,f9.4)",PSI, THETA, PHI

            DLIST(1) = PSI
            DLIST(2) = THETA
            DLIST(3) = PHI
            IKEY     = IKEY + 1
            CALL LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,3,IRTFLG)
	    
	    PSI = PHI
	    DO COUNT2=0,2
	    	      
	      PHI = 60.0 + COUNT2 * DELTA3
	      THETA = DEGREE
	      
C             PRINT "(f9.4,f9.4,f9.4)",PSI, THETA, PHI

              DLIST(1) = PSI
              DLIST(2) = THETA
              DLIST(3) = PHI
              IKEY     = IKEY + 1
              CALL LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,3,IRTFLG)
	      
	    ENDDO
	  ENDDO
	RETURN
	ENDIF	  


C -------------------------------- GROEL ----------------------	
C 7 FOLD ROTATIONAL SYMMETRY WITH INVERSION ::  PSI, THETA, PHI

        IF ((CSYM(1:1) .EQ. 'G') .OR. (CSYM(1:1) .EQ. 'g'))THEN
	
          COMMENT = '7 FOLD ROT. SYM. WITH INVERSION: PSI,THETA,PHI'
          CALL LUNDOCPUTCOM(LUNDOC,COMMENT,IRTFLG)

          COMMENT = 'LONG AXIS OF GROEL IS IN Z'
          CALL LUNDOCPUTCOM(LUNDOC,COMMENT,IRTFLG)

          DT       = 51.429
	  
	  PSI      = 0.0
	  THETA    = 0.0
          DLIST(1) = PSI
          DLIST(2) = THETA

	  DO IKEY=1,7
	     PHI      = (IKEY -1) * DT
             DLIST(3) = PHI
             CALL LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,3,IRTFLG)
	  ENDDO
       
	  PSI      =   0.0
	  THETA    = 180.0
          DLIST(1) = PSI
          DLIST(2) = THETA

          IKEY     = 8
	  PHI      = 0
          DLIST(3) = PHI
          CALL LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,3,IRTFLG)

	  DO IKEY=9,14
	     IC       = 15 - IKEY
	     PHI      = IC * DT
             DLIST(3) = PHI
             CALL LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,3,IRTFLG)
	  ENDDO

	  RETURN
	ENDIF	  

	
C -------------------------- ICOSAHEDRAL --------------------------------	

C       FOR THE LAST CASE OF ICOSAHEDRAL GROUP
C       INCLUDES EXTRA ROTATIONS, TO MOVE THEIR AXES TO ALIGN WITH OURS
        IF ((CSYM(1:1) .EQ. 'I') .OR.(CSYM(1:1) .EQ. 'i'))THEN
	
          COMMENT = 'ICOSAHEDRAL SYMMETRY  ::  PSI, THETA, PHI'
          CALL LUNDOCPUTCOM(LUNDOC,COMMENT,IRTFLG)
	
C         -----THIS PRINTS OUT THE "STANDARD" ROTATION FOUND ONLINE
          DELTAN = 36.0
          DO COUNT1 =0,1
	    DO COUNT2 =0,4
	  
	      PSI   = 0.0
	      THETA = COUNT1 * DELTA2
	      PHI   = COUNT2 * DELTA5
C	      PRINT "(f9.4,f9.4,f9.4)",PSI, THETA, PHI

              DLIST(1) = PSI
              DLIST(2) = THETA
              DLIST(3) = PHI
              IKEY     = IKEY + 1
              CALL LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,3,IRTFLG)
	    ENDDO
	  ENDDO
	
	  THETA = 63.4349
	  DO COUNT1 = 0,4
	    DO COUNT2 = 0,4
	  
	      PSI = COUNT1 * DELTA5
	      PHI = COUNT2 * DELTA5 + DELTAN
C	      PRINT "(f9.4,f9.4,f9.4)",PSI, THETA, PHI

              DLIST(1) = PSI
              DLIST(2) = THETA
              DLIST(3) = PHI
              IKEY     = IKEY + 1
              CALL LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,3,IRTFLG)

	    ENDDO
	  ENDDO
	
	  THETA = 116.5651
	  DO COUNT1 = 0,4
	    DO COUNT2 = 0,4
	  
	      PSI = COUNT1 * DELTA5 + DELTAN
	      PHI = COUNT2 * DELTA5 
C	      PRINT "(f9.4,f9.4,f9.4)",PSI, THETA, PHI

              DLIST(1) = PSI
              DLIST(2) = THETA
              DLIST(3) = PHI
              IKEY     = IKEY + 1
              CALL LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,3,IRTFLG)
	    ENDDO
	  ENDDO
C         -----END THE "STANDARD" PRINTOUT
	  


C         create Xbeta =A and Xbeta' = B
C          X_or_Z = 0
C	  DEGREE = BETA
C	  CALL MATCREATE(X_or_Z, DEGREE, A)
C	  B = TRANSPOSE(A)
C         create Xalpha =C and Xalpha' =D
C         X_or_Z = 0
C	  DEGREE = ALPHA
C	  CALL MATCREATE(X_or_Z, DEGREE, C)
C	  D = TRANSPOSE(C)
C         start Q loop
C         DO COUNT1=0,4
C           using 72.0 because = 360/5	    
C	    DEGREE = COUNT1 * DELTA5
C	    X_or_Z = 1	    
C	    CALL MATCREATE(X_or_Z, DEGREE, E)
C           calc Zdelta5 * Xalpha	    
C	    F = MATMUL(E,C)
C           calc X'alpha *Zdelta5 * Xalpha
C           E = MATMUL(D,F)
C           start Xk loop
C           DO COUNT2 =0,1
C             create X(k*delta2)	    
C	      DEGREE = COUNT2 * DELTA2
C	      X_or_Z = 0
C	      CALL MATCREATE(X_or_Z, DEGREE, G)
C             start Zj loop
C             DO COUNT3 =0,1
C               create Z(j*delta2) =AB	        
C		DEGREE = COUNT3 * DELTA2
C		X_or_Z = 1
C	        CALL MATCREATE(X_or_Z, DEGREE, AB)
C               start P loop
C               DO COUNT4=0,2
C		  create Z(i*delta3)
C                 DEGREE = COUNT4 * DELTA3
C		  X_or_Z = 1
C	          CALL MATCREATE(X_or_Z, DEGREE, BA)
C		  calc Z*X(beta)
C                  H = MATMUL(BA, A)
C                 calc X'(beta)*Z*X(beta) =P
C                  BA = MATMUL(B,H)
C                 calc Z(delta2)*P
C                  H = MATMUL(AB,BA)
C                 calc X(delta2)*Z(delta2)*P
C                  BA = MATMUL(G, H)
C                 calc Q*X(delta2)*Z(delta2)*P
C                  H = MATMUL(E,BA)
		  
C                 rotate to align 5fold axis with Z axis	  
C		  X_or_Z = 0
C		  DEGREE = ALPHA
C		  CALL MATCREATE(X_or_Z, DEGREE, M)		  
C                 rotate to align their X, our Y
C                 X_or_Z = 1
C		  DEGREE = 54
C		  CALL MATCREATE(X_or_Z, DEGREE, N)
C		  P = MATMUL(M,N)		  		  
C		  N = MATMUL(P,H)
C                 H = N
C		  CALL MATEXTRACT(H,PSI, THETA, PHI)
C		  PRINT "(f9.4,f9.4,f9.4)",PSI, THETA, PHI
C		ENDDO
C              ENDDO
C	    ENDDO
C	  ENDDO

	RETURN
	ENDIF

C       print warning if ever gets here
C        PRINT *, "Enter a supported value"  
 
9999    END SUBROUTINE SYMANG
        




C       ------MATCREATE------------------------------------------------

C       PURPOSE:CREATES A ROTATION MATRIX AROUND EITHER X OR Z ASSUMES
C       ROTATION CCW LOOKING TOWARDS 0 FROM +AXIS ACCEPTS 2 ARGUMENTS, 
C       0=X OR 1=Z (FIRST REG.)AND ROT ANGLE IN DEG.(SECOND REG.)

        SUBROUTINE MATCREATE(INXORZ, INDEGR, NEWMAT)
	
	  DOUBLE PRECISION, DIMENSION (1:3,1:3) :: NEWMAT
	  DOUBLE PRECISION                      :: INDEGR
	  DOUBLE PRECISION, PARAMETER:: PI = 3.14159265358979323846264
	  INTEGER                    :: INXORZ

C         create blank rotation matrix
          NEWMAT = 0
	  
C         change input of degrees to radians for fortran
	  INDEGR = INDEGR * (PI / 180.0)
	  
C         for X rot matrix, place 1 and add cos&sin values
          IF (INXORZ .EQ. 0) THEN
	    
	    NEWMAT(1,1) = 1
	    
	    NEWMAT(2,2) = COS(INDEGR)
	    NEWMAT(3,3) = COS(INDEGR)
	    
	    NEWMAT(3,2) = SIN(INDEGR)
	    NEWMAT(2,3) = -SIN(INDEGR)
	    
C         FOR Z ROT MATRIX, DO AS ABOVE, BUT FOR Z
          ELSEIF(INXORZ.EQ.1)THEN
	  
	    NEWMAT(3,3) = 1
	    
	    NEWMAT(1,1) = COS(INDEGR)
	    NEWMAT(2,2) = COS(INDEGR)
	    
	    NEWMAT(2,1) = SIN(INDEGR)
	    NEWMAT(1,2) = -SIN(INDEGR)
	    
C         ERROR CATCH
          ELSE
	    PRINT *, "Did not specify X or Z matrix"
	  ENDIF
	END SUBROUTINE MATCREATE





C       -------MATEXTRACT----------------------------------------------

C	USED TO CALCULATE ANGLES SPIDER EXPECTS FROM ROT. MATRIX.
C       ASSUMES SIN(THETA) IS POSITIVE (0-180 DEG)

        SUBROUTINE MATEXTRACT(ROTMAT, PSI, THETA, PHI)
	
	  DOUBLE PRECISION, DIMENSION (1:3,1:3) :: ROTMAT, W
	  DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979323846264

	  DOUBLE PRECISION            :: PSI, THETA, PHI
	  DOUBLE PRECISION            :: RADTHA
	  DOUBLE PRECISION            :: SINTHA
	  DOUBLE PRECISION            :: RADONE
	  DOUBLE PRECISION            :: RADTWO
	  DOUBLE PRECISION            :: ROUNDER,DT

C         CALCULATE THETA FROM LOWER/RIGHT CORNER
          !RADTHA = ACOS(ROTMAT(3,3))
          DT     = MAX(-1.0D0, MIN(1.0D0,ROTMAT(3,3)))
          RADTHA = DACOS(DT)
	  THETA  = RADTHA * (180.0 / PI)
	  SINTHA = SIN(RADTHA)

C	  CLOSE ENOUGH TEST, SET CORNER TO -1
          IF ( ABS(1-(ROTMAT(3,3)/(-1.0))) .LT. (1.0e-6) )THEN
	    ROTMAT(3,3) = -1.0
	  ENDIF

C	  CLOSE ENOUGH TEST, SET CORNER TO 1
          IF ( ABS(1-(ROTMAT(3,3)/(1.0))).LT.(1.0E-6) ) THEN
	    ROTMAT(3,3) = 1.0
	  ENDIF
 
C         SPECIAL CASE OF THETA ROTATION/ Y ROTAION = 180 OR 0
C         IF WE DO THIS, NEED ONLY ONE Z ROTATION
          IF((ABS(ROTMAT(3,3))).EQ.1) THEN
	  
	    PSI = 0
	    
C           FIND PHI, IF SIN=-, SWITCH SIGN ALL IN RADIANS
            !RADONE = ACOS(ROTMAT(1,1))
            DT      = MAX(-1.0D0, MIN(1.0D0,ROTMAT(1,1)))
            RADONE  = DACOS(DT)
            RADTWO  = ROTMAT(2,1)
	    
	    IF(RADTWO.LT.0) THEN
	      PHI = 2.0 * PI - RADONE
	      PHI = PHI * (180.0 / PI)
	      
	    ELSE
	      PHI = RADONE * (180.0 / PI)
	    ENDIF
	    
          ELSE
C           NORMAL CASE OF ALL THREE ROTATIONS
 
C           FIND PSI, IF SIN(PSI) =- THEN SWITCH AROUND
            !RADONE = ACOS(-ROTMAT(3,1) / SINTHA)
            DT      = -ROTMAT(3,1) / SINTHA
            DT      = MAX(-1.0D0, MIN(1.0D0,DT))
            RADONE  = DACOS(DT)
            RADTWO  = (ROTMAT(3,2) / SINTHA)
	    IF (RADTWO.LT.0) THEN
	      PSI = 2.0 * PI - RADONE
	      PSI = PSI * (180.0/PI)

	    ELSE
	      PSI = RADONE * (180.0/PI)
	    ENDIF
	    
C           FIND PHI, SIMILAR TO BEFORE
            !RADONE = ACOS(ROTMAT(1,3) / SINTHA)
            DT      = ROTMAT(1,3) / SINTHA
            DT      = MAX(-1.0D0, MIN(1.0D0,DT))
            RADONE  = DACOS(DT)
	    RADTWO  = ROTMAT(2,3) / SINTHA
	    IF(RADTWO.LT.O) THEN
	      PHI = 2.0 * PI - RADONE
	      PHI = PHI * (180.0 / PI)
	    ELSE
	      PHI = RADONE * (180.0 / PI)
	    ENDIF
	  ENDIF
	  
C         CATCH TO CHANGE 360 PHI TO 0
          IF( ABS(1.0-(PHI/(360.0))) .LT. (1.0e-2) )THEN
	    PHI = 0.0
	  ENDIF
	  
C         CATCH TO CHANGE REALLY SMALL PSI TO 0, FOR OCT.
          IF( ABS(1.0- ( (PSI+1)/1.0)).LT.(1.0E-4) )THEN
	    PSI = 0.0
	  ENDIF
	  
C         CATCH TO CHANGE REALLY SMALL PHI TO 0, FOR OCT.
          IF( ABS(1.0- ( (PHI+1)/1.0)).LT.(1.0e-4) )THEN
	    PHI = 0.0
	  ENDIF

	END SUBROUTINE MATEXTRACT
