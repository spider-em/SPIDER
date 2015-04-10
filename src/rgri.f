
C **********************************************************************
C                                                                      
C  RGRI     NEW                                  2/6/86 J.F.  
C           COSMETIC OUTPUT CHANGES              DEC 2008 ARDEAN LEITH                                                *
C           IMC FILE FDUM                        JUN 09 ARDEAN LEITH
C           COSMETIC                             NOV 11 ARDEAN LEITH
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
C  RGRI(NUMIM,KFAC,NKLA,KLAS,IDI,PK,GT,IV,LUNK,LUNI,NFAC)
C
C  PURPOSE:                                                             
C     READ IMC_ CLUSTER FILE PRODUCED BY SCLASSI AND PRINT OUT
C     LISTS OF CLASS MEMBERS, CENTERS OF GRAVITIES FOR EACH CLASS, AND
C     RE-CLASSIFICATION LOOKUP TABLE
C                                                                       
C  PARAMETERS:                                                          
C
C  NOTE: IMPORTANT CHANGE: NKLA IS REDUCED BEFORE DEUCL IS CALLED. 
C        THE ORIGINAL NUMBER OF CLUSTERS IS STORED IN NKLAO 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE RGRI(NUMIM,KFAC,NKLA,KLAS,IDI,PK,GT,IV,
     &                LUNK,LUNI,NFAC,KV)

      INCLUDE 'CMBLOCK.INC' 

      INTEGER :: NUMIM,KFAC,NKLA
      INTEGER :: KLAS(NUMIM), IDI(NUMIM)
      REAL    :: PK(NKLA), GT(NKLA, KFAC)
      INTEGER :: IV(NKLA)
      INTEGER :: LUNK,LUNI,NFAC
      INTEGER :: KV(KFAC)

C     NUMBER OF MAJOR CLASSES TO BE ANALYSED IN TERMS OF DISPERSIONS
C     10/20/87 TEMPORARILY CHANGED TO 11 JF
      INTEGER, PARAMETER :: NMAJ=11

      INTEGER, PARAMETER :: MAXPRT=1000  ! MAXIMUM NO. PER CLASS TO BE PRINTED
	
      INTEGER            :: BELONG(MAXPRT)
      REAL               :: DIST(NMAJ,NMAJ)
      COMMON /COMMUN/ DIST,BELONG

      REAL               :: COO(NFAC)
      DIMENSION          :: DISP(NMAJ),DELTA(3),NEXT(3)
 
C     READ CLASSIFICATION DATA FROM CLUSTER FILE ON: LUNK                                                                              
      READ(LUNK) (KLAS(I), I=1,NUMIM), 
     &           (IDI(I), I=1,NUMIM),
     &           (PK(L),L=1,NKLA), 
     &           ((GT(L,J),L=1,NKLA),J=1,KFAC)

	NMAJOR = MIN(NMAJ,NKLA)

C       LIST IMAGES BY CLASSES
	WRITE(NDAT,13)
13	FORMAT(/,'  LIST OF CLASS MEMBERS',/,'  CLASS')

	DO  K=1,NKLA
	   IKLA = 0
           DO I=1, NUMIM
             IF (KLAS(I) .EQ. K) THEN
                IKLA         = IKLA+1
                BELONG(IKLA) = IDI(I)
             ENDIF
           ENDDO
	   IF (IKLA .EQ. 0) CYCLE

	   IF (IKLA > MAXPRT) THEN
              WRITE(6,12)K,MAXPRT
12            FORMAT(I4,' *** MORE THAN ',I4,' IMAGES')
              CYCLE
	   ENDIF

	   WRITE(NDAT,11)K,(BELONG(I),I=1,IKLA)
11	   FORMAT(I6,10I7,(6X,10I7))

        ENDDO

C       LIST CLASS CENTER COORDINATES

	WRITE(NDAT,1001)(KV(IFAC),IFAC=1,KFAC)
1001	FORMAT(/,'  LIST OF CLASS CENTER COORDINATES',/,
     &           '   CLASS  SIZE    ',8(5X,I2,3X)/)
        DO N=1, NKLA                               
           WRITE(NDAT, 1002) N,INT(PK(N)+0.5),(GT(N,IFAC),IFAC=1,KFAC)
 1002      FORMAT(2X,I4,3X,I4,3X,12(F9.4,1X))                 
        ENDDO

C       RE-CLASSIFICATION LOOKUP TABLE

	WRITE(NDAT,1004)((I),I=1,NKLA)
1004	FORMAT(/,'  RE-CLASSIFICATION LOOKUP TABLE'/,
     &           '  ORIGINAL CLASS',/,
     &         8X,(20I3))

C       READ RE-CLASSIFICATION DATA FROM CLUSTER FILE ON: LUNK                                                                              
        KPART = NKLA - 1                                            
        DO KPRO = 2,KPART                                   
           READ(LUNK) (IV(J), J=1,NKLA)
                
           WRITE(NDAT, 1003) KPRO, (IV(J), J=1,NKLA) 
 1003      FORMAT(I6,2X,(40I3)) 
        ENDDO

C       COMPUTE CLASS DISPERSIONS AND INTERCLASS DISTANCES FOR NMAJOR
C       MAIN CLASSES

C       CLEAR ARRAY DIST
	DO K=1,NMAJOR
           DISP(K) = 0.0
           DO K1=1,NMAJOR
              DIST(K,K1) = 0.0
  	   ENDDO
        ENDDO

	CALL REWF(LUNI,1)     ! REWIND _IMC FILE 
	DO I=1,NUMIM
C          READ _IMC FILE ON: LUNI                                                                              
           READ(LUNI,*) (COO(IFAC),IFAC=1,NFAC), FDUM,FDUM,FDUM,FDUM

           K = KLAS(I)
           IF (K .LE. NMAJOR) THEN
              DO IFAC=1,KFAC
                 DISP(K) = DISP(K) + (COO(IFAC) - GT(K,IFAC))**2
  	      ENDDO
           ENDIF
  	ENDDO

	DO K=1,NMAJOR
  	   DISP(K) = SQRT(DISP(K) / PK(K))
        ENDDO

C       COMPUTE INTERCLASS DISTANCES
	DO K=1,NMAJOR
           DO K1=1,NMAJOR
             IF (K1 .NE. K) THEN
                DO IFAC=1,KFAC
                  DIST(K,K1) = DIST(K,K1)+(GT(K,IFAC) - GT(K1,IFAC))**2
  	        ENDDO
             ENDIF
          ENDDO
        ENDDO

C       SCALE DISTANCES
	DO K=1,NMAJOR
           DO K1=1,NMAJOR
              DIST(K,K1 )= SQRT(DIST(K,K1))
           ENDDO
        ENDDO

C       WRITE HEADING
        WRITE(NDAT,*) ' '
	WRITE(NDAT,1202)(I,I=1,10)
1202	FORMAT('  DISPERSIONS AND INTER-CLASS DISTANCES OF 10 LARGEST',
     &         ' CLUSTERS',//,
     &         '   CLASS      DISP       NEIGHBORS   ',10I7/)

C       FOR EACH CLUSTER, DETERMINE THE 3 CLOSEST NEIGHBORS

	DO K=1,NMAJOR
C          CLEAR TABLES
           DO J=1,3
              DELTA(J) = 100000.
              NEXT(J)  = 0
           ENDDO

           DO K1=1,NMAJOR
             IF (K1 .EQ.K) CYCLE
             IF (DIST(K,K1) .LT. DELTA(1)) THEN
		DELTA(3) = DELTA(2)
		DELTA(2) = DELTA(1)
		DELTA(1) = DIST(K,K1)
		NEXT(3)  = NEXT(2)
		NEXT(2)  = NEXT(1)
		NEXT(1)  = K1
		CYCLE
             ENDIF

             IF (DIST(K,K1) .LT. DELTA(2)) THEN
		DELTA(3) = DELTA(2)
		DELTA(2) = DIST(K,K1)
		NEXT(3)  = NEXT(2)
		NEXT(2)  = K1
		CYCLE
             ENDIF

             IF (DIST(K,K1) .LT. DELTA(3)) THEN
		DELTA(3) = DIST(K,K1)
		NEXT(3)  = K1
             ENDIF
   	   ENDDO

	   WRITE(NDAT,1201)K,DISP(K),(NEXT(J),J=1,3),
     &                     (DIST(K,K1),K1=1,K-1)
1201	   FORMAT(5X,I4,3X,F7.4,5X,3I3,6X,10F7.4)

      ENDDO
      WRITE(NDAT,*)' '
      WRITE(NDAT,*)' '

      END                                                                       
                                              
