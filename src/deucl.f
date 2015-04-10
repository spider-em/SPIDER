
C ++********************************************************************
C                                                                      *
C DEUCL.F                                                   02.03.81   *         
C                                                           01.12.86   *         
C           IMC FILE FDUM                      JUN 2009 ARDEAN LEITH
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
C*
C*  DEUCL(NKDIM,NKLA,MCARD,NUMIM,KFAC,NFAC,                       
C*              KLAS, D, T, PK, CI, LUNI )                             
C*                                                                     
C*  CONSTRUCTION OF THE MONO-INDEXED TABLE OF DISTANCES BETWEEN                 
C*  THE NKLA CLASSES IN THE EUCLIDEAN SPACE CHARACTERIZED BY THE                
C*  KFAC FIRST FACTORIAL COORDINATES.                                           
C*                                                                              
C*  INPUT:                                                                      
C*      NKDIM        = MAJORIZER OF NKLA FOR T(*,KFAC)                          
C*      NKLA         = NUMBER OF CLUSTERS FOR THE GIVEN PARTITION               
C*      MCARD        = NKLA * (NKLA-1) / 2  = DIMENSION FOR D(*)                 
C*      JFIN         = 2 * NKLA - 1   = DIMENSION OF WORKING ARRAYS             
C*      NUMIM        = NUMBER OF OBJECTS TO BE CLUSTERED                        
C*      KFAC,NFAC    = FIRST KFAC FACTORIAL COORDINATES FROM A TOTAL            
C*                     OF NFAC                                                  
C*      KLAS(NUMIM)  = GIVEN PARTITION IN NKLA CLUSTERS.                        
C*      LUNI         = FILE CONTAINING THE FACTORIAL COORDINATES OF             
C*                     THE GRAVITY CENTERS OF THE CLUSTERS.                     
C*                                                                              
C*   OUTPUT:                                                                    
C*      T(NKLA,KFAC) = TABLE CONTAINING THE COORDINATES OF THE NKLA             
C*                     CENTERS                                                  
C*      PK(JFIN)     = WEIGHT OF THE NKLA CLUSTERS                              
C*      D(MCARD)     = MATRIX OF DISTANCES BETWEEN CLUSTERS                     
C*                                                                              
C*      WORKING ARRAY ... CI(KFAC)                                        
C*                                                                              
C*   INTERNAL FUNCTION  ... MONO.                                               
C*                                                                              
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE DEUCL(NKDIM,NKLA,MCARD,NUMIM,KFAC,NFAC,                       
     &                  KLAS, D, T, PK, CI, LUNI)                              

      DIMENSION KLAS(NUMIM),D(MCARD),PK(NKLA),CI(NFAC),T(NKDIM,KFAC)                                                   

C     MONO-INDEXING OF THE TABLE OF DISTANCES                              

      MONO(K1,K2) = MIN(K1,K2) + ((MAX(K1,K2)-1)*(MAX(K1,K2)-2)/2)            

C     GRAVITY CENTER OF THE NKLA CLASSES                                   

      CALL REWF(LUNI, 1)                                                         
      DO  J=1,NKLA                                                       
         PK(J) = -1.0                                                             
      ENDDO 

      DO I = 1,NUMIM                                                      
         READ(LUNI,*) (CI(KF),KF=1,NFAC), FDUM,FDUM,FDUM,FDUM

         J = KLAS(I)

C        FOLLOWING LINE ADDED FEB. 04 al TO PREVENT ACCESS BEYOND PK
         IF (J .GT. NKLA)  CYCLE
                                                           
         IF (PK(J) .LE. -0.99) THEN                                
            PK(J) = 1.0                                                               
            DO K = 1,KFAC                                                       
               T(J,K) = CI(K)
	    ENDDO                                                             
         ELSE                       
            PK(J) = PK(J) + 1.0                                                       
            PKJ   = 1.0 / PK(J)                                                       
            DO   K = 1,KFAC                                                       
               CI(K)  = CI(K)  - T(J,K)                                                    
               T(J,K) = T(J,K) + PKJ*CI(K)                                                
	    ENDDO
        ENDIF
      ENDDO                                                                

C       CALCULATION OF THE MATRIX OF DISTANCES BETWEEN OBJECTS               
        DO K1 = 2,NKLA                                                      
           K3 = K1 - 1
                                                            
           DO K2 = 1,K3                                                        
              K1K2    = MONO (K1, K2)                                                     
              D(K1K2) = 0.0 
                                                             
              DO K = 1,KFAC                                                       
                 AJ      = T(K1,K) - T(K2,K)                                                 
                 D(K1K2) = D(K1K2) + AJ*AJ                                                  
              ENDDO
           ENDDO
	ENDDO

        DO K1 = 2,NKLA                                                      
           K3 = K1 - 1                                                            
           DO K2 = 1,K3                                                        
              K1K2    = MONO (K1, K2)                                                     
              D(K1K2) = PK(K1)*PK(K2)*D(K1K2) / (PK(K1)+PK(K2))                          
           ENDDO
        ENDDO

        END  
