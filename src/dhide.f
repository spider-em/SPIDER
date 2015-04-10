
C++*********************************************************************
C
C  DHIDE.FOR        ALTERED BY al JAN 87 TO PLACE PLOT IN PLOT METAFILE
C                   MAKES POSTSCRIPT OUTPUT NOW MAR 99 ARDEAN LEITH
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
C     DHIDE(X,Y,XG,G,XH,H,NG,MAXDIM,N1,NFNS,
C          XLNTH,YLNTH,XMINT,DELTAX,DELTAY)
C
C    PURPOSE:
C
C      THIS SUBROUTINE PRODUCES A 2-DIMENSIONAL REPRESENTATION OF       
C      A 3-DIMENSIONAL FIGURE OR SURFACE.                               
C
C      THE FIRST CALL TO HIDE IS FOR INITIALIZATION AND PLOTTING   
C      THE CURVE FARTHEST IN THE FOREGROUND.  ON EACH SUBSEQUENT        
C      CALL, A CURVE FARTHER IN THE BACKGROUND IS PLOTTED.             
C      X IS THE ABCISSA ARRAY FOR THE CURVE TO BE PLOTTED BY        
C      HIDE ON THIS CALL.  THE X VALUES MUST BE INCREASING.            
C      IF X(I) GE X(I+1) FOR SOME I, MAXDIM WILL BE SET TO ZERO,     
C      AND A RETURN WILL BE EXECUTED.                                
C      Y IS THE ORDINATE ARRAY.                                       
C      G VS. XG IS THE CURRENT VISUAL MAXIMUM FUNCTION ON EACH         
C      RETURN FROM HIDE.                                          
C      XH AND H ARE WORKING ARRAYS.                                   
C      ON EACH RETURN FROM HIDE, NG IS THE NUMBER OF POINTS IN      
C      THE CURRENT MAXIMUM FUNCTION.                                 
C      ON THE FIRST CALL, NG IS A NONPOSITIVE INTEGER WHICH         
C      SPECIFIES CERTAIN OPTIONS.                                    
C      -2 PLOT UNHIDDEN MINIMUM RATHER THAN MAXIMUM.  IN THIS         
C         CASE, G VS. XG WILL BE THE NEGATIVE OF THE VISUAL       
C         MINIMUM FUNCTION.                                       
C      0  PLOT MAXIMUM.                               
C      MAXDIM IS THE DIMENSION IN THE CALLING PROGRAM OF THE      
C      ARRAYS XG, G, XH, AND H.  IF ONE OF THESE ARRAYS WOULD     
C      HAVE BEEN OVERFLOWED, MAXDIM IS SET EQUAL TO ITS NEGATIVE, 
C      AND A RETURN IS EXECUTED.                                  
C      N1 IS THE NUMBER OF POINTS (X(I),Y(I)) TO BE PLOTTED IN    
C      A GIVEN CALL TO HIDE.                                      
C      IF N1 IS LESS THAN 0, Y VS. X WILL NOT BE PLOTTED, BUT ON  
C      SUBSEQUENT CALLS, PLOTTING WILL BE DONE AS IF              
C      ((X(I),Y(I)),I=1,-N1) HAD BEEN PLOTTED (WHERE UNHIDDEN).   
C      N1 WILL BE RETURNED AS ITS ABSOLUTE VALUE.                 
C      NFNS IS THE TOTAL NO. OF CURVES TO BE PLOTTED FOR THIS     
C      GRAPH IF TRANSLATING THE ARRAYS TO SIMULATE STEPPING IN    
C      THE DEPTH DIMENSION IS DESIRED.  IF NO TRANSLATION IS      
C      DESIRED, NFNS SHOULD BE NEGATIVE.  IF THE SAME TRANSLATION 
C      AS IN THE PREVIOUS CALL TO HIDE IS DESIRED, NFNS SHOULD BE 
C      ZERO.  THE NFNS=0 OPTION MAY BE SPECIFIED FOR INDIVIDUAL   
C      CURVES AFTER THE FIRST FOR A GIVEN GRAPH.  ALL             
C      TRANSLATIONS WHICH ARE PERFORMED WILL HAVE EQUAL STEP SIZE 
C      DETERMINED BY THE VALUES IN THE INITIAL CALL FOR XLNTH,    
C      YLNTH, AND NFNS.                                           
C
C      XLNTH IS THE LENGTH IN INCHES OF THE HORIZONTAL AXIS.      
C      IF XLNTH IS LESS THAN 0, THE X-AXIS AND THE DEPTH AXIS     
C      WILL NOT BE DRAWN.  IN ANY CASE, UNLESS THIS OPTION IS     
C      SUPPRESSED THROUGH NFNS, THE ITH CURVE WILL BE TRANSLATED  
C      (I-1)*(9.-ABS(XLNTH))/(NFNS-1) INCHES TO THE LEFT.  THIS   
C      PLUS A SIMILAR VERTICAL TRANSLATION IS DONE TO SIMULATE    
C      STEPPING IN THE DEPTH DIMENSION.                           
C
C      XMIN-(9.-ABS(XLNTH))*DELTAX WILL BE THE ABCISSA VALUE AT   
C      THE PLOTTING REFERENCE POINT (WHICH IS WHERE THE           
C      HORIZONTAL AND VERTICAL AXES WOULD INTERSECT IF DRAWN).    
C
C      YLNTH IS THE LENGTH OF THE VERTICAL AXIS IN INCHES.        
C      IF YLNTH IS LESS THAN 0, THE VERTICAL AND DEPTH AXES WILL  
C      NOT BE DRAWN.  BUT IN ANY CASE, UNLESS THIS OPTION IS      
C      SUPPRESSED THROUGH NFNS, THE ITH CURVE WILL BE TRANSLATED  
C      (I-1)*(6.-ABS(YLNTH))/(NFNS-1) INCHES UP TO SIMULATE       
C      STEPPING IN THE DEPTH DIMENSION.  YMINT-(6.-ABS(YLNTH))*   
C      DELTAY WILL BE THE ORDINATE VALUE AT THE PLOTTING          
C      REFERENCE POINT.                                           
C
C      IF TRANSLATIONS ARE PERFORMED, X AND Y WILL BE RESTORED TO 
C      THEIR ORIGINAL VALUES BEFORE THE RETURN TO THE CALLING     
C      PROGRAM.                                                   
C
C      NOTE THAT IF ABS(XLNTH)=9, AND ABS(YLNTH)=6, THERE WILL BE 
C      NO TRANSLATION, AND, IF AXES ARE NOT DRAWN, THE 
C      DIMENSIONS OF THE PLOT ARE UNSPECIFIED.                    
C
C      DELTAX IS THE X DATA INCREMENT PER INCH FOR THE PLOT. DX IS 
C      1 / DELTAX.  DX, XMIN AND DELTAX DETERMINE THE PLOTTING SCALE 
C      FOR X.(SEE ABOVE.)
C
C      DY, YMINT AND DELTAY, SIMILARLY, DETERMINE THE SCALE FOR Y.        
C      IF AN ERROR RETURN IS MADE FROM HIDE, ALL FURTHER CALLS        
C      WILL RESULT ONLY IN THE EXECUTION OF A RETURN UNLESS           
C      MAXDIM IS RESET TO A POSITIVE VALUE.                           
C
C--*******************************************************************

       SUBROUTINE DHIDE(X,Y,XG,G,XH,H,NG,MAXDIM,N1,NFNS,                 
     &        XLNTH,YLNTH,XMINT,DELTAX,FMIN,DELTAY,XMAXT,LUNPOS)

       PARAMETER (NSIZE = 2000)
       COMMON DATA(3,NSIZE)     
C      WORK BUFFER FOLLOWS DATA

       DIMENSION X(*),Y(*),XG(*),G(*),H(*),XH(*) 

C     THE ONLY PURPOSE OF THE FOLLOWING EQUIVALENCE STATEMENT IS
C     TO SAVE STORAGE.                                          
      EQUIVALENCE (K1,IWHICH),(K2,SLOPE),(FNSM1,Z1),        
     &            (IGGP1,K1),(K1,N2)                        

C     EPS1 IS THE RELATIVE ABCISSA INCREMENT USED TO SIMULATE   
C     DISCONTINUITIES IN THE MAXIMUM FUNCTION.                  
      DATA EPS1/ .00001/                                    
         
C     MAKE VARIABLES STATIC.                                                    
      SAVE

C     THE FOLLOWING STATEMENT FUNCTION COMPUTES THE ORDINATE ON 
C     THE LINE JOINING (XI,YI) AND (XIP1,YIP1) CORRESPONDING TO 
C     THE ABCISSA XX.    
                                         
      F(XX,XI,YI,XIP1,YIP1) = YI+(XX-XI)*(YIP1-YI)/(XIP1-XI)

      IF (MAXDIM .LE .0) RETURN                                       
      DO  I = 2,N1                            
         IF (X(I-1) .GE. X(I)) THEN
C          OVERFLOW??
           MAXDIM = 0
           RETURN
         ENDIF
      END DO                                            

      IFPLOT = 1                                          
      IF (N1 .LE. 0) THEN
         N1 = -N1
         IFPLOT = 0
      ENDIF

C     IF > FIRST PLOT GOTO 5000
      IF (NG .LE. 0) THEN
        IF (N1+4 .GT. MAXDIM) THEN
C          OVERFLOW ??
           MAXDIM = -MAXDIM        
           RETURN
        ENDIF

        YMINT = FMIN
        DX    = 1 / DELTAX
        DY    = 1 / DELTAY

C       WE WANT SIGN = 1 IF WE ARE PLOTTING MAXIMUM, = -1 IF MINIMUM.
        SIGN = 1.                                    
        IF (NG .LT. -1) SIGN = -1.
                                                      
C       THE KTH CURVE TO BE PLOTTED WILL (OPTIONALLY) BE        
C       TRANSLATED BY THE VECTOR (-DXIN,DYIN)*(K-1) TO SIMULATE 
C       STEPPING IN THE DEPTH DIMENSION.                      
        IF (NFNS .GT. 0) THEN
           FNSM1 = NFNS-1  
           DXIN  = 3.*DELTAX/FNSM1  
           DYIN  = 3.*(DELTAY/FNSM1)
        ENDIF

        IF (XLNTH .GE. 0.) THEN
C         AXIS ARE DISABLED!!!!
C         CALL ROUTINE TO DRAW THE HORIZONTAL AXIS.  THE 
C         LEFT END IS SPECIFIED IN INCHES RELATIVE TO THE REFERENCE    
C         POINT BY THE FIRST TWO ARGUMENTS.                            
          CALL POSAXIS('X',XMINT,XMAXT,0,0,XLNTH,YLNTH,DX,LUNPOS,IRTFLG)
        ENDIF

        IF (YLNTH .GE. 0.) THEN
C         AXIS ARE DISABLED!!!!
C          DRAW THE DEPTH AXIS.
           CALL POSEG(LUNPOS, 0.0,0.0, -XLNTH,YLNTH)

C          DRAW THE VERTICAL AXIS.  THE BOTTOM POINT IS SPECIFIED IN
C          INCHES RELATIVE TO THE REFERENCE POINT BY THE FIRST TWO  
C          ARGUMENTS.                                               
           CALL POSAXIS('Y',YMINT,YMAXT,0,0,XLNTH,YLNTH,DY,
     &                  LUNPOS,IRTFLG)
        ENDIF

C       CURVES SUCCESSIVELY FARTHER IN THE BACKGROUND WILL BE   
C       PLOTTED WHERE THEY ARE NOT HIDDEN BY G VS. XG.  G VS XG 
C       WILL BE UPDATED EACH TIME A NEW CURVE IS DRAWN AND WILL BE
C       THE VISUAL MAXIMUM (OR MINIMUM) FUNCTION OF THE CURVES    
C       ALREADY PLOTTED.

        INDEXT=3        
        DO  J = 1,N1   
           XG(INDEXT) = X(J) 
           G(INDEXT) = SIGN*Y(J)
           INDEXT = INDEXT+1
        END DO    

C       THE FOLLOWING PRECAUTIONARY STEP IS USED IN PLACE OF A        
C       TEST IN SUBROUTINE LOOKUP TO SEE IF THE VALUE FOR WHICH WE     
C       WANT AN INDEX IS OUTSIDE THE TABLE.                          
C       THE LAST XG VALUE WILL BE SET EQUAL TO THE LAST ABCISSA   
C       OF THE CURVE TO BE PLOTTED IN THE NEXT CALL TO HIDE.       

        EPS      = EPS1*(ABS(XMINT)+ABS(DELTAX))                          
        NG       = N1+4                                            
        XG(1)    = -FNSM1*DXIN+XMINT-ABS(XMINT)-ABS(XG(3))-1.         
        XG(2)    = XG(3)-EPS                                        
        XG(N1+3) = XG(N1+2)+EPS                                   
        ZZ       = YMINT                                               
        IF (SIGN.LT.0.) ZZ = -YMINT-50.*DELTAY                   
        G(1)     = ZZ                                                
        G(2)     = ZZ                                                
        G(N1+3)  = ZZ                                              
        G(NG)    = ZZ                                                

C       CALL ROUTINE TO PRODUCE A LINE PLOT OF                 
C       (X(I),Y(I),I=1,N1) - THE CURVE FARTHEST IN THE FOREGROUND.

C       XSTART IS THE X VALUE AT THE REFERENCE POINT.        
        XSTART = XMINT-(9.-ABS(XLNTH))*DELTAX             
    
C       WORKS BETTER FOR POSTSCRIPT, MAR 99 AL       
        XSTART = 0.0

        DO  I = 1,N1
           DATA(1,I) = X(I) * DX + XSTART
           DATA(2,I) = Y(I) * DY
        END DO
        IF (IFPLOT.EQ.1) CALL POARAYF(LUNPOS,DATA,N1,.FALSE.,.FALSE.)

        DXKK   = 0.                                               
        DYKK   = 0.                                               
        RELINC = DELTAX / DELTAY                                  
        XG(NG) = SIGN                                            
        RETURN                                                     
      END IF

C     STATEMENT 5000 IS REACHED IF ANY CURVE EXCEPT THE FARTHEST   
C     IN THE FOREGROUND IS TO BE PLOTTED.                              
C5000 CONTINUE
 
      SIGN = XG(NG)                                             
      XG(NG) = X(N1)                                             

C     TRANSLATE THE ARRAYS BEFORE PLOTTING TO SIMULATE STEPPING        
C     IN THE DEPTH DIMENSION.                                     
      IF(NFNS) 52,48,49     
              
   49 DXKK = DXKK+DXIN                       
      DYKK = DYKK+DYIN      
                
   48 DO J = 1,N1                  
         Y(J) = SIGN*(Y(J)+DYKK)     
         X(J) = X(J)-DXKK  
      END DO
                   
   52 CALL LOOKUP(X(1),XG(1),JJ)  
             
      IF(JJ.GE.MAXDIM) GO TO 700               
      DO J = 1,JJ                             
         XH(J) = XG(J)                           
         H(J) = G(J)     
      END DO                              
      IG = JJ+1                                 
      XH(IG) = X(1)                             
      H(IG) = F(X(1),XG(JJ),G(JJ),XG(IG),G(IG))              

C     WE WILL BE MAKING TABLE LOOKUPS FOR AN INCREASING SEQUENCE      
C     OF NUMBERS - THEREFORE, WE DO NOT HAVE TO SEARCH FROM 
C     FIRST OF THE (XG AND X) TABLES EACH TIME.  HENCE INDEXG   
C     AND INDEXT.   
                                             
      INDEXG = JJ                                        
      INDEXT = 1                               
      Z1 = X(1)                                        
      F1 = H(IG)-Y(1)                                   
      IT = 2                                             
      JJ = IG                                            
      IF (H(IG) .LT. Y(1)) THEN
        IF (JJ .GE. MAXDIM) GO TO 700                     
        JJ = IG+1                                          
        H(JJ) = Y(1)                                       
        XH(JJ) = Z1+EPS                                    
      ENDIF
      LAST = 0                                       
      X1 = Z1                                 
C                                      
C     FIND THE FIRST ZERO, Z2, OF THE FUNCTION G-Y TO THE RIGHT OF Z1.
C 1100 IF(XG(IG) .LT. X(IT)) GO TO 1001 
                   
1100  IF(XG(IG) .GE. X(IT)) THEN                                                              
C        LOOK FOR A ZERO BETWEEN X1 AND  X(I).
         IWHICH = 0                                   
         X2 = X(IT)                                    
         F2 = F( X2, XG(IG-1), G(IG-1), XG(IG), G(IG) ) - Y(IT)    
         IT = IT+1                                
C      GO TO 1002                                      
      ELSE
C        LOOK FOR A ZERO BETWEEN X1 AND XG(IG).
C 1001   CONTINUE
         X2 = XG(IG)                                           
         IWHICH = 1                                     
         F2 = G(IG) - F( X2, X(IT-1), Y(IT-1), X(IT), Y(IT) )    
         IG = IG+1 
      END IF                                     

C     THE FUNCTION (G-Y) HAS A ZERO Z2 SUCH THAT X1 LE Z2 LE X2     
C     IF AND ONLY IF (G-Y AT X1) * (G-Y AT X2) LE 0.            
C     (G-Y IS ASSUMED, FOR PLOTTING PURPOSES, TO BE LINEAR ON     
C     EACH INTERVAL (X1,X2).)                                  
C1002 CONTINUE
      IF (F1*F2.GT.0.) GO TO 1005                           
      IGG = IG-1-IWHICH                                     
      ITT = IT-2+IWHICH                        
      IF (X2 .EQ. X1) THEN
         Z2=X1
        GOTO 1006
      ENDIF

      SLOPE = (F2-F1)/(X2-X1)
      IF (ABS(SLOPE*RELINC) .GT. 1.E-6) THEN
C       OTHERWISE, COMPUTE THE ZERO Z2.                               
        Z2 = X1-F1/SLOPE                                          

      ELSE
C       IF G AND Y DIFFER IMPERCEPTIBLY (FOR PLOTTING PURPOSES) 
C       ON THE INTERVAL (X1,X2), SET Z2=X2.  THIS STEP PREVENTS        
C       DIVISION BY ZERO.                                             
        Z2 = X2                                                 
      ENDIF
      GO TO 1006                                        

C     IF NO ZERO WAS FOUND BETWEEN X1 AND X2, CONTINUE THE ZERO SEARCH
 1005 X1 = X2                                                  
      F1 = F2                                                   
      IF (IT .LE. N1) GO TO 1100                             

C     IF THE END OF THE X TABLE HAS BEEN REACHED, CONSIDER THE          
C     INTERVAL FROM THE LAST ZERO FOUND TO THE END OF THE X        
C     TABLE (PLOT, UPDATE MAXIMUM FUNCTION AS INDICATED).          

 1008 LAST = 1                                          
      Z2 = X(N1)                                           
      CALL LOOKUP(Z2,XG(INDEXG),IGG)                          
      IGG = INDEXG+IGG-1                                     
      ITT = N1-1                                          

C     IT IS NECESSARY TO PLOT Y VS. X ON THE INTERVAL (Z1,Z2)    
C     ONLY IF Y IS UNHIDDEN AT EACH ZZ SUCH THAT Z1 LT ZZ LT Z2. 
C     EFFICIENCY IN THE TABLE LOOKUP.                            
C     WE CHOOSE ZZ NEAR THE LEFT END OF THE INTERVAL FOR         
C     NOTE THAT IT IS MORE EFFICIENT TO CHOOSE THIS VALUE FOR ZZ 
C     THAN, SAY, .99*X(INDEXT)+.01*X(INDEXT+1), WHICH WOULD      
C     ELIMINATE ONE OF THE TWO TABLE LOOKUPS, BUT WOULD          
C     NECESSITATE A TEST TO DETERMINE IF ZZ WAS BETWEEN Z1 AND Z2. 

 1006 ZZ = .99*Z1+.01*Z2                                     
      CALL LOOKUP(ZZ,X(INDEXT),K1)                           
      CALL LOOKUP(ZZ,XG(INDEXG),K2)                          
      K1 = K1+INDEXT-1                                       
      K2 = K2+INDEXG-1                                       
      IF (F(ZZ,X(K1),Y(K1),X(K1+1),Y(K1+1)).GT.               
     &   F(ZZ,XG(K2),G(K2),XG(K2+1),G(K2+1))) GO TO 7        

C     IF Y IS HIDDEN BETWEEN Z1 AND Z2, UPDATE THE MAXIMUM FUNCTION. 
C     FOR GENERALITY, THE MAXIMUM FUNCTION IS UPDATED EVEN IF    
C     THIS IS THE (NFNS)TH CURVE.                                

      IF(JJ+IGG-INDEXG.GE.MAXDIM) GO TO 700                  
      IF (INDEXG .NE. IGG) THEN
        J1 = INDEXG+1                                          
        DO  I = J1,IGG                                       
           JJ = JJ+1                                           
           XH(JJ) = XG(I)                                      
           H(JJ) = G(I) 
        END DO                                       
      ENDIF

  712 JJ = JJ+1                                              
      XH(JJ) = Z2                                            
      H(JJ) = F(Z2,XG(IGG),G(IGG),XG(IGG+1),G(IGG+1))        
      INDEXG = IGG                                           
      INDEXT = ITT                                           
      GO TO 60                                               

C     IF T IS NOT HIDDEN BETWEEN Z1 AND Z2, UPDATE THE MAXIMUM   
C     FUNCTION AND PLOT.                                         
    7 NGRAPH = ITT-INDEXT+2                                  
      IF (JJ+NGRAPH-1 .GT. MAXDIM) GO TO 700                    
      N2 = JJ                                                
      IF (NGRAPH .NE. 2) THEN
        J1 = INDEXT+1                                          
        DO I = J1,ITT                                       
          JJ = JJ+1                                           
          XH(JJ) = X(I)                                       
          H(JJ) = Y(I)
        END DO                                        
      ENDIF
   9  JJ = JJ+1                                              
      XH(JJ) = Z2                                            
      H(JJ) = F(Z2,X(ITT),Y(ITT),X(ITT+1),Y(ITT+1))          

C     CALL ROUTINE TO PRODUCE PLOT OF (XH(I),H(I),I=N2,N2+NGRAPH-1)

      N = 0
      DO  I = N2,N2+NGRAPH-1
         N = N + 1
         DATA(1,N) = XH(I) * DX + XSTART
         DATA(2,N) =  H(I) * DY
      END DO
      IF (IFPLOT.EQ.1) CALL POARAYF(LUNPOS,DATA,N,.FALSE.,.FALSE.)

      INDEXT = ITT           
      INDEXG = IGG           

   60 IF (LAST .NE. 1) THEN
        X1 = X2                
        F1 = F2                
        Z1 = Z2                
C       AFTER PLOTTING AND/OR UPDATING THE MAXIMUM FUNCTION ON THE 
C       INTERVAL (Z1,Z2), SEARCH FOR THE NEXT ZERO IF THE END OF   
C       THE ABCISSA TABLE XT HAS NOT BEEN REACHED.                 
        IF (IT .LE. N1) GO TO 1100                                
        GO TO 1008                                             
      ENDIF

C     AFTER Y VS. X HAS BEEN PLOTTED, FINISH UPDATING AND STORE  
C     THE NEW MAXIMUM FUNCTION.                                  
C     ALLOW FOR THE POSSIBILITY THAT THE PREVIOUS MAXIMUM        
C     FUNCTION EXTENDS TO THE RIGHT OF THE FUNCTION JUST PLOTTED.

   61 IF (XG(NG) .LE. XG(NG-1)) NG = NG-1                       
      IF (XG(NG) .GT. X(N1)) THEN
        IF (JJ+3+NG-IGG .GT. MAXDIM) GO TO 700                    
        XH(JJ+1) = XH(JJ)+EPS                                  
        JJ = JJ+1                                              
        H(JJ) = F(X(N1),XG(IGG),G(IGG),XG(IGG+1),G(IGG+1))     
        IGGP1 = IGG+1                                          
        DO J = IGGP1,NG                                     
          JJ = JJ+1                                           
          XH(JJ) = XG(J)                                      
          H(JJ) = G(J)
        END DO                                        
      ENDIF

      NG = JJ+2                                              
      IF (NG.GT.MAXDIM) GO TO 700                             
      DO  I = 1,JJ                                         
         G(I) = H(I)                                         
         XG(I) = XH(I)     
      END DO                                  
                                                            
C     THE FOLLOWING PRECAUTIONARY STEP IS USED IN PLACE OF A     
C     TEST IN SUBROUTINE LOOKUP TO SEE IF THE VALUE FOR WHICH WE 
C     WANT AN INDEX IS OUTSIDE THE TABLE.                        
C     THE LAST XG VALUE WILL BE SET EQUAL TO THE LAST ABCISSA    
C     OF THE NEXT CURVE TO BE PLOTTED. 
                          
      XG(JJ+1) = XG(JJ)+EPS                                  
      G(JJ+1) = YMINT+DYKK                                   
      IF (SIGN .LT. 0.) G(JJ+1) = -YMINT-50.*DELTAY+DYKK        
      G(NG) = G(JJ+1)                                        

C     RESTORE ARRAYS X AND Y BEFORE RETURNING.                   
   66 CONTINUE
C      IF(NFNS.LT.0) GOTO 53
      IF(NFNS .GE. 0) THEN
         DO I = 1,N1                                         
           X(I) = X(I)+DXKK                                    
           Y(I) = SIGN*Y(I)-DYKK
         END DO
      END IF    
C53   CONTINUE
                           
      XG(NG) = SIGN                                          
      RETURN                                                 


C     IF STATEMENT 700 IS REACHED, DIMENSIONS WOULD HAVE BEEN    
C     EXCEEDED.  SEE COMMENTS ON CALLING SEQUENCE FOR HIDE.      
  700 MAXDIM = -MAXDIM                                       
      GO TO 66                                               
      END                                                    



C++*********************************************************************
C
C LOOKUP(X,XTBL,J)
C
C  PURPOSE:    THIS SUBROUTINE PERFORMS A TABLE LOOK UP
C
C  PARAMETERS:   X       X VALUE TO BE LOOKED UP
C                XTBL    ARRAY CONTAINING TABLE
C                J       INDEX FOR THE TABLE
C
C  CALLED BY:    DHIDE
C
C  NOTE:         BECAUSE OF PRECAUTIONS TAKEN IN DHIDE, A TEST TO
C                SEE IF X IS OUTSIDE THE TABLE IS UNNECESSARY.
C
C--********************************************************************

      SUBROUTINE LOOKUP(X,XTBL,J)

      DIMENSION XTBL(1)

      J = 2

    4 IF (XTBL(J) - X) 1,2,3

    1 J = J+1
      GO TO 4

    3 J = J - 1

    2 RETURN
      END

