
C **********************************************************************
C  CLAST.F                                                             *
C           COSMETIC   CHANGES                  DEC 08 ARDEAN LEITH    *                                               *
C           PIXEL FILE FDUM                     JUN 09 ARDEAN LEITH    *
C           IMC FILE FDUM                       JUN 09 ARDEAN LEITH    *
C           KV                                  NOV 11 ARDEAN LEITH    *
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
C                                                                      *         
C   CLAST(N2DIM,NFAC,NUMIM,KFAC,NCLAS,G,KLAS,CI,PCLAS,LUNI)
C                                                                      *         
C   PURPOSE:
C       DETERMINATION OF CLASSES FROM COORDINATES.  THE DISTANCE               
C       BETWEEN OBJECT I AND THE NCLAS CENTERS IS CALCULATED.                  
C       THE COORDINATES OF THE NCLAS CENTERS ARE IN G(NCLAS,*).                
C       OBJECT I IS ASSIGNED TO THE NEAREST CLASS.  THE GRAVITY                
C       CENTERS OF THE NEWLY FORMED NCLAS CLASSES ARE DETERMINED               
C       AND USED IN THE NEXT ITERATION.                                        
C       FILLINF OF KLAS(NUMIM).                                                
C                                                                      *         
C  CALL TREE:   SCLASSI - NOYAU - PARST  - RETIR                       *
C                  |        |            - CLAST                       *
C                  |        |            - STABK - SHELK               *
C                  |        |                                          *
C                  |      NOYAU - DEUCL                                * 
C                  |      NOYAU - CHAVA                                *
C                  |      NOYAU - DENDRO - DENLST                      * 
C                  |      NOYAU - COUPE                                *
C                  |                                                   *
C               SCLASSI - RGRI                                         *    
C                                                                      *         
C **********************************************************************

      SUBROUTINE CLAST(N2DIM,NFAC,NUMIM,KFAC,NCLAS,G,KLAS,CI,PCLAS,            
     &                 KV,LUNI)
                                                                               
      INTEGER :: N2DIM,NFAC,NUMIM,KFAC,NCLAS
      REAL    :: G(N2DIM,KFAC)
      INTEGER :: KLAS(NUMIM)
      REAL    :: CI(NFAC),PCLAS(NCLAS)
      INTEGER :: KV(KFAC)
      INTEGER :: LUNI

      CALL REWF(LUNI, 1)      ! REWIND _IMC FILE TO RECORD 1
                                                         
      DO N = 1,NCLAS                                                      
         DO KF = 1,KFAC                                                      
            G(NCLAS+N,KF) = 0.0                                                        
	 ENDDO

         PCLAS(N) = 0.0                                                             
      ENDDO
                                                                         
C     ASSIGNMENT OF OBJECT TO THE NEAREST CENTER          
      DO I = 1,NUMIM   
         READ(LUNI,*) (CI(NF),NF=1,NFAC), FDUM,FDUM,FDUM,FDUM
         SUP   = 1.0E+15                                                          
         KAT   = 1 
                                                                
         DO  N=1,NCLAS                                                      
            DIST  = 0.0  

            DO KF = 1,KFAC
               KT   = KV(KF)
               DC   = CI(KT) - G(N,KF)                                                  
               DIST = DIST   + DC*DC 
            ENDDO

            IF (DIST > SUP) CYCLE                                
            KAT   = N                                                                 
            SUP   = DIST                                                              
         ENDDO 
                                                               
         KLAS(I)    = KAT                                                              
         PCLAS(KAT) =  PCLAS(KAT) + 1.0  
                                            
         DO KF = 1,KFAC
            KT   = KV(KF)
            G(NCLAS+KAT,KF) = G(NCLAS+KAT,KF) + CI(KT)                                 
         ENDDO
      ENDDO
                                                                              
C     PROTECTION AGAINST AN EMPTY CLASS (REJECTED FURTHER) 
                                                                              
      DO N=1,NCLAS                                                      
         IF (PCLAS(N) > 0.0)  CYCLE
                             
         DO KF = 1,KFAC
            G(NCLAS+N,KF) = 1.0E+8                                   
         ENDDO
      ENDDO 
                                                                               
C     UPDATE COORDINATES OF THE GRAVITY CENTERS 
      DO N=1,NCLAS                                                     
         IF (PCLAS(N) <= 0.0) PCLAS(N) = 1.0E-8 
          
         DO KF = 1,KFAC 
            G(N,KF) = G(NCLAS+N,KF) / PCLAS(N)
         ENDDO
      ENDDO

      END                                                                    
