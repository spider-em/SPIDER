C **********************************************************************
C
C  CHAVA.F                                                    04.11.81          
C                                                             01.12.86          
C               COSMETIC OUTPUT CHANGES            DEC 08 ARDEAN LEITH                                                *
C               MDIM_8                             MAY 13 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
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
C  CHAVA(NKLA, MDIM_8, JFIN, D, PK, VAL, LA, LB, NT, NO, LUNK, MODE)
C*  
C    PURPOSE:                                                                          
C*      HIERARCHICAL ASCENDENT CLASSIFICATION USING THE VARIANCE AS            
C*      AN AGGREGATION CRITERION.  THE DISTANCES ARE IN THE MONO-              
C*      INDEXED TABLE D(MCAR).                                                 
C*                                                                            
C*   INPUT:                                                                    
C*      NKLA            = NUMBER OF ELEMENTS TO BE CLASSIFIED                  
C*      MDIM_8          = NKLA * (NKLA-1) /2                             
C*      JFIN            = 2*NKLA - 1                                          
C*      D(MDIM_8)       = TABLE OF DISTANCES DESTROYED AFTER EXECUTION         
C*      PK(1,...,NKLA)  = WEIGHT OF THE NKLA ELEMENTS                          
C*                                                                             
C*   OUTPUT:                                                                   
C*      PK(NKLA+1,...)      = WEIGHT OF THE NODES OF THE HIERARCHY             
C*      LA(NKLA,LB(NKLA)    = SENIORS AND JUNIORS OF THE                       
C*                            CLASSIFICATION                                   
C*      VAL(NKLA+1,...)     = NODE INDICES                                     
C*      NT(NKLA+1,...)      = SIZE OF CLASSES                                  
C*      NO(NKLA+1,...)      = WORKING ARRAY                                    
C*                                                                             
C*      INTERNAL FUNCTION ... MONO                                             
C*                                                                             
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

       SUBROUTINE CHAVA(NKLA, MDIM_8, JFIN,  D, PK,                         
     &                  VAL, LA, LB, NT, NO, LUNK, MODE)

       INCLUDE 'CMBLOCK.INC'

       INTEGER            :: NKLA,JFIN,LUNK,MODE

       INTEGER*8          :: MDIM_8
       REAL               :: D(MDIM_8), PK(JFIN),VAL(JFIN)                 
       INTEGER            :: LA(NKLA), LB(NKLA),NT(JFIN),NO(JFIN)                                

       CHARACTER (LEN= 1) :: IAST,KLIGN(110)
       INTEGER            :: K1,K2,L,J,JDEB,II1,II2,II,I1,I3,I2

       REAL               :: VMAX,DINF

       INTEGER            :: mono

       DATA  IAST /'*'/

C      MONO-INDEXING OF THE TABLE OF DISTANCES                             
       MONO(K1,K2)= MIN0(K1,K2) + ((MAX0(K1,K2)-1)*(MAX0(K1,K2)-2)/2)            

C      GENERATION OF IDENTIFIERS, NUMBERS IN A4                            

        IF     (MODE == 1) THEN
           WRITE(NDAT,*) ' CLUSTERING CRITERION:  SINGLE LINKAGE'
        ELSEIF (MODE == 2) THEN
           WRITE(NDAT,*) ' CLUSTERING CRITERION:  COMPLETE LINKAGE' 
        ELSEIF (MODE == 3) THEN
           WRITE(NDAT,*) ' CLUSTERING CRITERION:  AVERAGE LINKAGE' 
        ELSEIF (MODE == 4) THEN
           WRITE(NDAT,*) ' CLUSTERING CRITERION:  CENTROID METHOD' 
        ELSEIF (MODE == 5) THEN
           WRITE(NDAT,*) ' CLUSTERING CRITERION:  WARDS METHOD' 
        ENDIF

C       ARRAY OPERATION
        IF (MODE <= 3)  D = SQRT(D)

        WRITE (NDAT,650)                                                           
  650   FORMAT (/,'  DESCRIPTION OF THE HIERARCHY NODES',//,                  
     &            '   NO. SENIOR JUNIOR NO.  WEIGHT  INDEX',/)

        DO L=1,100                                                          
            KLIGN(L) = IAST 
        ENDDO 
                                                          
        DO J = 1,JFIN                                                       
           NO(J)  = J                                                                 
           VAL(J) = 0.0                                                               
           NT(J)  = 1                                                                 
        ENDDO

C       CALCULATIONS AND LISTINGS                                           
C       DETERMINE PAIRWISE DISTANCES FOR ALL PAIRS

        JDEB  = NKLA + 1                                                          
        VMAX  = 0.0 
                                                              
        REWIND(LUNK)          ! REWIND CLUSTER FILE   
                                                         
        DO J=JDEB,JFIN                                                   
           DINF  = D(1)                                                              
           II1   = 1                                                                 
           II2   = 2                                                                 
           II    = 2*NKLA  - J + 1 
                                                  
           DO I1 = 2,II                                                       
              I3    = I1 - 1                                                            
              DO I2 = 1,I3                                                       
                 I1I2  = MONO (I1, I2)                                                     
                 IF (D(I1I2) .GE. DINF) CYCLE                               
                 DINF  = D(I1I2)                                                           
                 II1   = I1                                                                
                 II2   = I2                                                                
              ENDDO                                                                
           ENDDO  
                                                              
           IR      = J - JDEB + 1                                                      
           LA(IR)  = MIN0 (NO(II1), NO(II2))                                           
           LB(IR)  = MAX0 (NO(II1), NO(II2))                                           
           VAL(J)  = DINF                                                              
           LAI     = LA(IR)                                                            
           LBI     = LB(IR)                                                            
           NT(J)   = NT(LAI) + NT(LBI)                                                 
           VA      = VAL(LAI)                                                          
           VB      = VAL(LBI)                        
           P1      = PK(LAI)
           P2      = PK(LBI)
           PK(J)   = P1 + P2
           NO(II1) = J 
                                                               
           DO III = 1,II                                                      
              IF (III==II1 .OR. III == II2) CYCLE                               
              M1    = MONO (III, II1)                                                   
              M2    = MONO (III, II2)                                                   
              JS    = NO(III)                                                           
              VJS   = PK(JS)

              IF (MODE == 1) THEN
                    D(M1) = AMIN1(D(M1),D(M2))
              ELSEIF (MODE == 2) THEN
                    D(M1) = AMAX1(D(M1),D(M2))
              ELSEIF (MODE == 3) THEN
                    D(M1) =(P1*D(M1)+P2*D(M2))/(P1+P2)
              ELSEIF (MODE == 4) THEN
                    D(M1) = ((P1*D(M1)+P2*D(M2))/(P1+P2))-
     &                      (P1*P2*DINF/((P1+P2)**2))
              ELSEIF (MODE == 5) THEN
                    D(M1) = ((VA+VJS)*D(M1)+(VB+VJS)*D(M2)-
     &                   VJS*DINF)/(VA+VB+VJS)
              ENDIF 
           ENDDO 
                                        
           IF (II2 .NE. II) THEN                               
               NO(II2)= NO(II)                                                           
               IT    = II - 1 
                                                           
               DO III = 1,IT                                                      
                  IF (III == II2) CYCLE                               
                  M1    = MONO (III, II2)                                                   
                  M2    = MONO (III, II)                                                    
                  D(M1) = D(M2)                                                             
               ENDDO 
            ENDIF                                                  
            IF (VAL(J) > VMAX) VMAX = VAL(J)                            

C           OVERWRITE CLUSTER FILE ON: LUNK                                                       
            WRITE(LUNK) J,LA(IR),LB(IR),NT(J),PK(J),VAL(J)                         
        ENDDO
                                                                
      REWIND(LUNK)             ! REWIND CLUSTER FILE                                                          
      NKLA1 = NKLA - 1  
   
C     WRITES DESCRIPTION OF HIERARCHY NODES
                                                     
      DO  JJ = 1,NKLA1                                                    
         J    = JJ + NKLA                                                          
         IR   = J - JDEB + 1

C        READ FROM CLUSTER FILE ON: LUNK                                                       
         READ(LUNK) J,LA(IR),LB(IR),NT(J),PK(J),VAL(J)                         
         LIG   = 90.0 * VAL(J) / VMAX + 1.0                                            
         IF (LIG > 75) LIG = 75 
                                
         WRITE(NDAT,660) J,LA(IR),LB(IR),NT(J),PK(J),VAL(J),                       
     &                              (KLIGN(L),L=1,LIG)                          
  660    FORMAT(1X,I5,3I5,2(1PG10.2,1X),90A1)
      ENDDO

      WRITE(NDAT,*) ' '

      END                                                                    

