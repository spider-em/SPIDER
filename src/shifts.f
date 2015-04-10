C++*********************************************************************
C
C SHIFTS.F    EXTRACTED FROM UTIL2               ARDEAN LEITH 03/24/12
C             INTEGER SHIFT IN FOURIER NOW       ARDEAN LEITH 09/19/14
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
C  SHIFTS(FOURIER, LUN1,LUN2, BUF1,BUF2,
C         NX,NY,NZ, SHX,SHY,SHZ, IRTFLG)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE SHIFTS(FOURIER, LUN1,LUN2, BUF1,BUF2,
     &                  NX,NY,NZ, SHX,SHY,SHZ, IRTFLG)
    
      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      LOGICAL              :: FOURIER 
      INTEGER              :: LUN1,LUN2  
      REAL                 :: BUF1(*), BUF2(NX*6)
      INTEGER              :: NX,NY,NZ
      REAL                 :: SHX,SHY,SHZ
      INTEGER              :: IRTFLG

      INTEGER              :: IX,IY,IZ,IER,INS,I
      INTEGER              :: NXLD

         IRTFLG  = 1

         IF (SHX == 0 .AND.
     &       SHY == 0 .AND.
     &       SHZ == 0)  THEN
C           NO SHIFT, JUST COPY -------------------------------- COPY
            DO I = 1,NY*NZ
               CALL REDLIN(LUN1,BUF1,NX,I)
               CALL WRTLIN(LUN2,BUF1,NX,I)
            ENDDO

         ELSEIF (IFORM == 1  .AND. .NOT. FOURIER .AND. 
     &          (REAL(IFIX(SHX)) == SHX) .AND.
     &          (REAL(IFIX(SHY)) == SHY)) THEN
C              INTEGER SHIFT ------------------------------------- 2D I

               IX = SHX
               IY = SHY
               CALL SHIFT2_N(LUN1,LUN2, BUF1,BUF2, NX,NY, IX,IY)

         ELSEIF (IFORM == 3  .AND..NOT. FOURIER .AND.
     &          (SHX   == FLOAT(IFIX(SHX))) .AND.
     &          (SHY   == FLOAT(IFIX(SHY))) .AND.
     &          (SHZ   == FLOAT(IFIX(SHZ)))) THEN

C              INTEGER SHIFT  ----------------------------------- 3D I
               IX  = SHX
               IY  = SHY
               IZ  = SHZ

               CALL SHIFT3_N(LUN1,LUN2, BUF2, 
     &                       NX,NY,NZ, IX,IY,IZ)
 
         ELSEIF (FOURIER)  THEN
C           FOURIER SHIFT  ----------------------------------- FOURIER

            NXLD = NX+2-MOD(NX,2)

            CALL READV(LUN1,BUF1, NXLD,NY, NX,NY,NZ)

	    INS = +1          ! FORWARD FFT
	    CALL  FMRS_3(BUF1, NX,NY,NZ, INS)
	    IF (INS == 0) RETURN
               
	    CALL SHIFT_PF(BUF1,NXLD/2, NX,NY,NZ,
     &                    SHX,SHY,SHZ)

	    INS = -1          ! INVERSE FFT
	    CALL FMRS_3(BUF1, NX,NY,NZ, INS)
	    IF (INS == 0) RETURN

            CALL WRITEV(LUN2,BUF1, NXLD,NY, NX,NY,NZ)

         ELSEIF (IFORM == 1)  THEN
C           BILINEAR  REAL SHIFT ------------------------------ 2D REAL
            CALL SHIFTR_N(LUN1,LUN2,BUF2, NX,NY, 
     &                       1,NY,1, SHX,SHY)

         ELSEIF (IFORM == 3) THEN

C           BILINEAR REAL SHIFT  ------------------------------ 3D REAL

            CALL SHIFT_3D_N(LUN1,LUN2,BUF2,BUF1,
     &                       NX,NY,NZ, SHX,SHY,SHZ)
         ENDIF   ! END OF: IF (IFORM == 3) THEN

         IRTFLG = 0

         END



C ++********************************************************************
C                                                                      *
C SHIFT2                                                                     *
C                                                                      *
C **********************************************************************
C 
C SHIFT2(LUN1,LUN2,BUF1,BUF2,NX,NY,NXS,NYS)
C
C PURPOSE: INTEGER SHIFT A 2D ARRAY FROM LUN1 CIRCULARLY BY NXS,NYS.  
C          OUTPUT TO:  LUN2.
C
C PARAMETERS:
C         LUN1       LOGICAL UNIT NUMBER OF INPUT IMAGE
C         LUN2       LOGICAL UNIT NUMBER OF OUTPUT IMAGE
C         BUF1,BUF2  REAL BUFFERS OF SIZE NX
C         NX,NY      DIMENSIONS OF IMAGE
C         NXS,NYS    NUMBER OF SAMPLES AND SHY TO BE SHIFTED
C
C--*******************************************************************

      SUBROUTINE SHIFT2_N(LUN1,LUN2, BUF1,BUF2, NX,NY, NXS,NYS)

      IMPLICIT NONE

      REAL    :: BUF1(NX),BUF2(NX)
      INTEGER :: LUN1,LUN2, NX,NY, NXS,NYS

      INTEGER :: NS,NR,I,I1

      NS = MOD(NXS,NX)
      NR = MOD(NYS,NY)

      IF (NR < 0) NR = NR + NY

      IF (NS .NE. 0)  THEN
         DO  I = 1,NY                  ! X & Y SHIFT 
            CALL REDLIN(LUN1,BUF1,NX,I)

            CALL SHIFT1_N(BUF1,BUF2, NX,NXS)

            I1 = I + NR
            IF (I1 > NY) I1 = I1-NY

            CALL WRTLIN(LUN2,BUF2,NX,I1)
         ENDDO

      ELSE
         DO  I = 1,NY                   ! Y SHIFT ONLY
            CALL REDLIN(LUN1,BUF1,NX,I)

            I1 = I + NR
            IF (I1 > NY) I1 = I1 - NY   ! CIRCULAR

            CALL WRTLIN(LUN2,BUF1,NX,I1)
         ENDDO
      ENDIF


      END

C **********************************************************************
C
C SHIFT1(ARRAY1,ARRAY2,NX,IS)
C                                                                      
C PURPOSE: INTEGER SHIFT AN ARRAY, A, CIRCULARLY BY IS
C
C PARAMETERS: ARRAY*   REAL ARRAY OF DIMENSION NX 
C             B        OUTPUT REAL ARRAY
C             NX       DIMENSION
C             IS       INTEGER SHIFT
C
C--*******************************************************************

      SUBROUTINE SHIFT1_N(ARRAY1,ARRAY2, NX,IS)

      IMPLICIT NONE

      REAL    :: ARRAY1(NX), ARRAY2(NX)
      INTEGER :: NX,IS

      INTEGER :: ISH,NX1,I
 
      ISH = MOD(IS,NX)

      IF (ISH < 0) ISH = ISH + NX
       
      NX1 = -ISH - 1 + NX

      DO  I = 1,NX
	 ARRAY2(I) = ARRAY1(MOD(I+NX1,NX) + 1)
      ENDDO

      END

C ++********************************************************************
C                                                                      *
C SHIFT3                                                               *
C                                                                      *
C **********************************************************************
C                                                                      
C SHIFT3(LUN1,LUN2, BUF2, NX,NY,NZ, NXS,NYS,NZS)
C
C PURPOSE:  INTEGER SHIFT VOLUME STORED ON LUN1 
C           CIRCULARLY BY NXS,NYS,NZS.  OUTPUT TO  LUN2.
C
C PARAMETERS:
C         LUN1         LOGICAL UNIT NUMBER OF INPUT IMAGE
C         LUN2         LOGICAL UNIT NUMBER OF OUTPUT IMAGE
C         BUF2         REAL BUFFER OF SIZE NX*3
C         NX,NY,NZ     DIMENSIONS OF IMAGE
C         NXS,NYS,NZS  NUMBER OF SAMPLES AND SHY TO BE SHIFTED
C
C--*******************************************************************

	SUBROUTINE SHIFT3_N(LUN1,LUN2, BUF2,
     &                    NX,NY,NZ, NXS,NYS,NZS)

        IMPLICIT NONE

        INTEGER :: LUN1,LUN2 
        REAL    :: BUF2(3*NX)
        INTEGER :: NX,NY,NZ, NXS,NYS,NZS

        INTEGER :: NZT,L,JROFF,LEND,JWOFF

	NZT = MOD(NZS,NZ)
	IF (NZT < 0) NZT = NZT + NZ

	IF (NZT == 0) THEN
 	   DO L = 1,NZ
	      JROFF = (L-1) * NY
	      CALL SHIFTT_N(LUN1,LUN2, BUF2, 
     &                    NX,NY, NXS,NYS, JROFF,JROFF)
	   ENDDO

        ELSE

	   LEND = NZ - NZT
	   DO L = 1,LEND
	      JROFF = (L-1) * NY
	      JWOFF = (L+NZT-1) * NY
	      CALL SHIFTT_N(LUN1,LUN2, BUF2, 
     &                    NX,NY, NXS,NYS,JROFF,JWOFF)
	   ENDDO

	   DO L = LEND+1,NZ
	      JROFF = (L-1) * NY
	      JWOFF = (L-LEND-1) * NY
	      CALL SHIFTT_N(LUN1,LUN2, BUF2, 
     &                    NX,NY, NXS,NYS,JROFF,JWOFF)
	   ENDDO
   	ENDIF

	END

C++*******************************************************************
C                                                                     
C SHIFTT(LUN1,LUN2,BUF2(3*NX),NX,NY,NXS,NYS,JROFF,JWOFF)        
C                                                                     
C PURPOSE: SHIFT A TWO-DIMENSIONAL ARRAY THAT IS PART OF A VOLUME                              
C                                                                       
C       INTEGER SHIFT A TWO-DIMENSIONAL ARRAY OR PICTURE   
C       STORED ON LUN1 CIRCULARLY BY NXS,NYS;  OUTPUT TO    
C       ON LUN2. RECORD NUMBER OFFSETS ARE USED TO GET INPUT FROM       
C       APPROPRIATE SLICE AND PUT IT INTO POSITION WITH Z-OFFSET        
C                                                                       
C PARAMETERS:        
C         LUN1           LOGICAL UNIT NUMBER OF INPUT IMAGE             
C         LUN2           LOGICAL UNIT NUMBER OF OUTPUT IMAGE            
C         NX,NY          DIMENSIONS OF IMAGE                            
C         NXS,NYS        NUMBER OF SAMPLES AND SHY TO BE SHIFTED       
C         JROFF,JWOFF    RECORD NUMBER OFFSETS ON LUN1,LUN2 FOR READ    
C                        AND WRITE, RESPECTIVELY                                       
C                                                                       
C--*******************************************************************  
C                                                                       
      SUBROUTINE SHIFTT_N(LUN1,LUN2, BUF2,
     &                  NX,NY, NXS,NYS, JROFF,JWOFF)    

      IMPLICIT NONE

      INTEGER :: LUN1,LUN2 
      REAL    :: BUF2(3*NX)
      INTEGER :: NX,NY, NXS,NYS, JROFF,JWOFF

      INTEGER :: NS,NR,NA,NA2,I,NS1,K,I1

      NS  = MOD(NXS,NX)
      NR  = MOD(NYS,NY)

      NA  = NX + 1    ! OFFSET 1        
      NA2 = NX + NX   ! OFFSET 2 
    
      IF (NR < 0) NR = NR + NY

      IF (NR == 0) THEN
         DO I = 1,NY
            CALL REDLIN(LUN1,BUF2(NA),NX,I+JROFF)

            IF (NS < 0) THEN
               NS1 = -NS
               DO  K = 1,NS1
                 BUF2(NA2+K) = BUF2(NA+K-1)
               ENDDO

            ELSEIF (NS > 0) THEN
               NS1 = NS
               DO  K = 1,NS1
                 BUF2(NA-K) = BUF2(NA2-K+1)
               ENDDO
            ENDIF

            CALL WRTLIN(LUN2,BUF2(NA-NS),NX,I+JWOFF)
         ENDDO



      ELSE
         DO  I = 1,NY
            CALL REDLIN(LUN1,BUF2(NA),NX,I+JROFF)

            I1 = I + NR
            IF (I1 > NY) I1 = I1-NY

            IF (NS < 0) THEN
               NS1 = -NS
               DO  K = 1,NS1
                  BUF2(NA2+K) = BUF2(NA+K-1)
               ENDDO

            ELSEIF (NS > 0) THEN
               NS1 = NS
               DO  K = 1,NS1
                  BUF2(NA-K) = BUF2(NA2-K+1)
               ENDDO
            ENDIF

            CALL WRTLIN(LUN2,BUF2(NA-NS),NX,I1+JWOFF)
         ENDDO
      ENDIF

      END



C **********************************************************************
C
C  SHIFTR(LUN1,LUN2,BUF, NX,NY, NXGO,NXEND,NXSHIP,SHX,SHY)
C
C  PURPOSE: SHIFT A 2D IMAGE  9/26/81 NON-INTEGER SHIFT
C           STORED ON LUN1 CIRCULARLY BY SHX,SHY;  OUTPUT TO LUN2.
C
C  PARAMETERS
C         LUN1           LOGICAL UNIT NUMBER OF INPUT IMAGE
C         LUN2           LOGICAL UNIT NUMBER OF OUTPUT IMAGE
C         BUF            REAL BUFFER OF SIZE NX
C         NX,NY          DIMENSIONS OF IMAGE
C         NXGO,NXEND
C         NXSHIP         STARTING ROW, ENDING, AND SKIPING FACTOR
C         SHX,SHY      SHIFT VECTOR COMPONENTS IN SAMPLE AND ROW DIR.
C
C  JMC, 1/1/87, CHANGED TO ACCEPT
C               A GENERAL VALUE FOR THE STARTING AND ENDING SHY, AS
C               WELL AS A SKIPPING FACTOR. CHANGES ARE FOR USE IN 3-D    
C
C--*******************************************************************

      SUBROUTINE SHIFTR_N(LUN1,LUN2, BUF, NX,NY,
     &                  NXGO,NXEND,NXSHIP, SHX,SHY)

        REAL :: BUF(NX*6)

	NXS = SHX
	NYS = SHY
	DX  = SHX - NXS

	IF (SHX .LT. 0.) DX = 1 + DX
	DY = SHY - NYS
	IF (SHY .LT. 0.) DY = 1 + DY

	C1 = (1-DX) * (1-DY)
	C2 = DX * (1-DY)
	C3 = DY * (1-DX)
	C4 = DX * DY

C	WRITE(5,999)DX,DY,C1,C2,C3,C4
C999	FORMAT(6F10.4)
 
        NS = MOD(NXS,NX) + 1
        NR = MOD(NYS,NY) + 1
	IF (SHX .LT. 0.) NS = NS-1
	IF (SHY .LT. 0.) NR = NR-1

C       ADDRESSES USED
        NA  = NX + 1
C       NA  = INPUT BUFFER I OFFSET +1
	NA3 = 3 * NX + 1
C       NA3 = INPUT BUFFER II OFFSET +1
	NA5 = 5 * NX
C       NA5 = OUTPUT BUFFER OFFSET
	NA1 = 2 * NX + 1
        NA2 = 3 * NX
        IF (NR .LT. 0) NR = NR + NY

C       INITIALIZE FIRST BUFFER LINE
        CALL REDLIN(LUN1,BUF(NA),   NX,NXEND)
	CALL REDLIN(LUN1,BUF(1),    NX,NXEND)
	CALL REDLIN(LUN1,BUF(NA+NX),NX,NXEND)

	I     = 0
        NNNY  = NXGO - NXSHIP

80	I     = I + 1
        NNNY  = NNNY + NXSHIP
	NOFF1 = NA
	NOFF2 = NA3

	IF (MOD(I,2) .EQ. 0) NOFF2 = NA
	IF (MOD(I,2) .EQ. 0) NOFF1 = NA3

	NO1 = NOFF1 - NS
	NO2 = NOFF2 - NS

C       NOFF1  POSITION OF FIRST ELEMENT OF OLD LINE
C       NOFF2  POSITION OF FIRST ELEMENT OF NEW LINE

	CALL REDLIN(LUN1,BUF(NOFF2),NX,NNNY)

C       FOR NS=0, TWO NEIGHBORS ON EITHER END OF BUFFER ARE NEEDED
	BUF(NOFF2+NX) = BUF(NOFF2)
	BUF(NOFF2-1)  = BUF(NOFF2+NX-1)

        I1 = I + NR- 1
        IF (I1 .GT. NY) I1 = I1 - NY
        IF (I1 .LT. 1)  I1 = I1 + NY
        NNROI1 = NXGO + (I1-1) * NXSHIP

        IF (NS < 0) THEN
           NS1 = -NS
           DO K = 1,NS1+1
	      BUF(NOFF2+NX+K-1) = BUF(NOFF2+K-1)
	   ENDDO

        ELSEIF (NS > 0) THEN
           NS1 = NS
           DO K = 1,NS1+1
              BUF(NOFF2-K) = BUF(NOFF2+NX-K)
	   ENDDO
        ENDIF

   	DO K=1,NX
           BUF(NA5+K) = BUF(NO2+K-1)*C2 + BUF(NO2+K)*C1 + 
     &                  BUF(NO1+K-1)*C4 + BUF(NO1+K)*C3
        ENDDO

        CALL WRTLIN(LUN2,BUF(NA5+1),NX,NNROI1)

C	IF (I.LT.5) WRITE(3,201)(BUF(K),K=1,6*NX)
C201	FORMAT(16F7.3)

	IF (I .LT. NY) GOTO 80

	END


C++*********************************************************************
C
C SHIFT_3D.F 
C
C **********************************************************************
C  PURPOSE: SHIFT A 3D ARRAY OR PICTURE
C           STORED ON LUN1 CIRCULARLY BY SHX,SHY;  OUTPUT TO LUN2.
C
C  PARAMETERS
C         LUN1           LOGICAL UNIT NUMBER OF INPUT IMAGE
C         LUN2           LOGICAL UNIT NUMBER OF OUTPUT IMAGE
C         BD,OUT         REAL BUFFER OF SIZE NX+
C         NX,NY          DIMENSIONS OF IMAGE
C         NXGO,NXEND
C         NXSHIP         STARTING ROW, ENDING, AND SKIPPING FACTOR
C         SHX,SHY      SHIFT VECTOR COMPONENTS IN SAMPLE AND ROW DIR.
C
C **********************************************************************

         SUBROUTINE SHIFT_3D_N(LUN1,LUN2,BD,OUT,NX,NY,NZ,
     $                        SHXI,SHYI,SLISI)

         REAL :: BD(NX,4),OUT(NX)

         SHX = -SHXI
         SHY = -SHYI
         SLIS = -SLISI

1        IF (SHX .LT. 0.0)  THEN
            SHX = SHX + NX
            GOTO 1
         ENDIF

2        IF (SHY .LT. 0.0)  THEN
            SHY = SHY + NY
            GOTO 2
         ENDIF

3        IF (SLIS .LT. 0.0)  THEN
            SLIS = SLIS + NZ
            GOTO 3
         ENDIF
C
         NXS = SHX
         DX  = SHX-NXS
         NYS = SHY
         DY  = SHY-NYS
         NZS = SLIS
         DZ  = SLIS-NZS

         C1 = (1-DX)*(1-DY)*(1-DZ)
         C2 =    DX *(1-DY)*(1-DZ)
         C3 = (1-DX)*   DY *(1-DZ)
         C4 = (1-DX)*(1-DY)*   DZ
         C5 =    DX *   DY *(1-DZ)
         C6 =    DX *(1-DY)*   DZ
         C7 = (1-DX)*   DY *   DZ
         C8 =    DX *   DY *   DZ

         NSR = NY * NZ
         DO K=1,NZ
            KS = MOD(K+NZS-1,NZ) + 1
            DO J=1,NY
               IF (J .EQ. 1)  THEN
                  JS1 = (KS-1)*NY + MOD(J+NYS-1,NY)+1
                  JS2 = (KS-1)*NY + MOD(J+NYS,NY)+1
                  JS3 = MOD(KS*NY,NSR) + MOD(J+NYS-1,NY)+1
                  JS4 = MOD(KS*NY,NSR) + MOD(J+NYS,NY)+1
                  J1  = 1
                  J2  = 2
                  J3  = 3
                  J4  = 4
                  CALL  REDLIN(LUN1,BD(1,J1),NX,JS1)
                  CALL  REDLIN(LUN1,BD(1,J2),NX,JS2)
                  CALL  REDLIN(LUN1,BD(1,J3),NX,JS3)
                  CALL  REDLIN(LUN1,BD(1,J4),NX,JS4)
               ELSE
                  JS2 =      (KS-1)*NY + MOD(J+NYS,NY)+1
                  JS4 = MOD(KS*NY,NSR) + MOD(J+NYS,NY)+1
                  JT1 = J1
                  JT3 = J3
                  J1  = J2
                  J3  = J4
                  J2  = JT1
                  J4  = JT3
                  CALL REDLIN(LUN1,BD(1,J2),NX,JS2)
                  CALL REDLIN(LUN1,BD(1,J4),NX,JS4)
               ENDIF

               DO I=1,NX
                  IS  = MOD(I+NXS-1,NX)+1
                  IS1 = MOD(I+NXS,NX)+1

                  OUT(I)= C1*BD(IS,J1) + C2*BD(IS1,J1) + C3*BD(IS,J2) +
     &                    C4*BD(IS,J3) + C5*BD(IS1,J2) + C6*BD(IS1,J3) +
     &                    C7*BD(IS,J4)+C8*BD(IS1,J4)
               ENDDO
               CALL  WRTLIN(LUN2,OUT,NX,(K-1)*NY+J)
	    ENDDO
	 ENDDO

         END
