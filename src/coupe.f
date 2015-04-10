
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

C*--------------------------------------------------------------------*         
C*                                                                    *         
C*                                                                    *         
C*      INITIAL DENDROGRAM CHARACTERIZED BY NKLA SUMMITS IS TRUN-     *         
C*      CATED INTO KPART FINAL CLASSES.  OBJECTS ARE ASSIGNED         *         
C*      TO THE NEW CLASSES.                                           *         
C*                                                                    *         
C*   INPUT:                                                           *         
C*      NUMIM,NKLA,KPART,PK(NKLA),LA(NKLA),LB(NKLA) (CF,CHAVA)        *         
C*                                                                    *         
C*   OUTPUT:                                                          *         
C*      KLAS(NUMIM)      = NEW CLASSIFICATION CONTAINING KPART CLASSES*         
C*      IV(NKLA)         = ASSIGNMENT OF OLD CLASSES TO NEW CLASSES   *         
C*      NT(NKLA)         = SIZES OF THE KPART CLASSES                 *         
C*                                                                    *         
C*   WORKING ARRAY ... IW(NKLA)                                       *         
C*                                                                    *         
C*--------------------------------------------------------------------*         

      SUBROUTINE COUPE(NUMIM,NKLA,KPART, PK,LA,LB,IV,KLAS,NT, IW) 
             
      DIMENSION IV(NKLA),LA(NKLA),LB(NKLA),KLAS(NUMIM),PK(NKLA),                
     &          NT(NKLA),IW(NKLA)                                               

      JDEB  = NKLA + 1                                                          
      JFIN  = 2*NKLA - 1                                                        
                                                                              
C     INDICATOR IV(NKLA) OF THE AGGREGATION IN KPART CLASSES               
                                                                               
      DO   L = 1,NKLA                                                       
         IV(L) = 0 
      ENDDO 
                                                               
      DO   J = JDEB,JFIN                                                    
         NT(1) = J                                                                 
         KPT   = 0                                                                 
         JI    = 1 
                                                                
   20    IF (NT(JI) .LE. NKLA) THEN                                
            K       = NT(JI)                                                            
            KPT     = KPT + 1                                                           
            IW(KPT) = K                                                                
            JI      = JI - 1                                                            
         ELSE 
                               
            IJ      = JI + 1                                                            
            NI      = NT(JI) - JDEB + 1                                                 
            NT(IJ)  = LA(NI)                                                            
            NT(JI)  = LB(NI)                                                            
            JI      = JI + 1 
         ENDIF                                                           
         IF (JI .NE. 0) GO TO  20
                                
         IF (J .LE. JFIN-KPART+1) THEN                                
            DO   KK = 1,KPT                                                       
               JPP     = IW(KK)                                                            
               IV(JPP) = J                                                                
	    ENDDO
         ENDIF
         I1    = IW(1)                                                             
         I2    = IW(KPT)                                                           
      ENDDO

      KKK   = 1                                                                 
      NKLA1 = NKLA - 1                                                          
      DO  IL1 = 1,NKLA1                                                   
         IF(IV(IL1) .NE. 0) THEN                                
            IF (IV(IL1) .LT. KKK) CYCLE                               
            IV1     = IV(IL1)                                                           
            IV(IL1) = KKK                                                              
            IL3     = IL1 + 1
                                                           
            DO IL2 = IL3,NKLA                                                   
               IF (IV(IL2) .NE. IV1) CYCLE                                
               IV(IL2) = IV(IL1)                                                          
            ENDDO                                                                
         ELSE 
            IV(IL1 )= KKK 
         ENDIF                                                             
         KKK   = KKK + 1
      ENDDO
                                                               
      IF (IV(NKLA) .EQ. 0) IV(NKLA) = KKK 
                          
      END                                                                    
