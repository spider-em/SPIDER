
C **********************************************************************
C                                                                      *         
C  STABK                                                     06.07.80  *         
C                                                            01.12.86  *         
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
C  STABK(NUMIM, NKLA, KLAS, U, IU, IV, KTOT, KFIN,NIND)
C 
C  PURPOSE:
C                                                                   
C  INPUT KLAS(NUMIM) CONTAINS THE CROSSED-PARTITION WITH KTOT                   
C       CATEGORIES.  PRODUCT OF THE NBASE INITIAL PARTITIONS.                   
C  OUTPUT KLAS(NUMIM) CONTAINS A PARTITION OF NKLA CLASSES.                     
C       THE (NKLA-1) FIRST CLASSES ARE COMPOSED OF STABLE OBJECTS IN            
C       THE NBASE INITIAL PARTITIONS.  NIND OBJECTS IN THE LAST ONE.           
C       WORKING ARRAYS U(*), IU(*), IV(*).                                      
C       DIMENSION KTOT = NCLAS ** NBASE.  OUTPUT, IU(1,...,KTOT)                
C       CONTAINS THE SIZE OF THE STABLE CLASSES U(1,...,KFIN)                   
C       CUMULATIVE PERCENTAGES OF SIZE.                                         
C 
C  CALLS:  SHELK                                                                
C                                                                               
C **********************************************************************

      SUBROUTINE STABK(NUMIM, NKLA, KLAS, U, IU, IV, KTOT, KFIN,NIND)          

      INCLUDE 'CMBLOCK.INC' 

      DIMENSION  U(KTOT), IU(KTOT), IV(KTOT), KLAS(NUMIM)                       

C     THE SIZES OF THE CROSSED-PARTITION ARE IN U(*)                        

      DO K = 1,KTOT                                                       
         U(K) = 0.0                                                                
      ENDDO

      DO I = 1,NUMIM                                                      
         K    = KLAS(I)                                                           
         U(K) = U(K) + 1.0                                                        
      ENDDO

C     THE MOST IMPORTANT (NKLA - 1) CLASSES ARE KEPT.  THE REMAINING            
C     ARE IN NKLA.
                                                              
C     SORT VECTOR U() IN ASCENDING ORDER, DESTROYS U, SORTED IN: IU
      CALL SHELK (KTOT, U, IU ) 
                                                
      DO K = 1,KTOT                                                       
         M     = IU(K)                                                             
         IV(M) = KTOT - K + 1                                                      
      ENDDO

      DO I = 1,NUMIM                                                       
         L     = KLAS(I)                                                           
         KO    = IV(L)                                                             
         IF (KO .GT. NKLA) KO = NKLA                                
         KLAS(I) = KO
      ENDDO 
                                                            
      DO K = 1,KTOT                                                       
         KS    = KTOT - K + 1                                                      
         IU(K) = U(KS) + 0.00001                                                   
      ENDDO

C     RESULTS OF THE CLASSIFICATION                                        

      KFIN  = 0                                                                 
      NIND  = 0                                                                 
      DO K = 1,KTOT                                                       
         IF (K .LT. NKLA) NIND = NIND + IU(K)                      
         IF (IU(K) .LE. 0) EXIT                                
         KFIN  = KFIN + 1                                                          
         U(K)  = 100.0 * FLOAT(IU(K)) / FLOAT(NUMIM)                                 
      ENDDO

      NIND  = NUMIM - NIND                                                      
      DO K = 2,KFIN                                                       
         U(K) = U(K) + U(K-1)                                                     
      ENDDO
           
      RETURN                                                                 
      END                                                                    
                                                       

