
C ++********************************************************************
C                                                                      *
C HISMAP.F   ADAPTED FROM HPLAN.FOR  ON NOV 8 1986 BY ARDEAN LEITH
C            BORDERS CHANGED TO SOLID LINES FEB 88 ARDEAN LEITH
C            CHAR  ID PASSED AUG 99  ARDEAN LEITH
C                                                                      *
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
C   HISMAP(IDIM,NPTS,JX,JY,X,Y,ID,MOD,NLIGN,NPAGE,PEX,NDAT,NBAND) 
C   PREPARES PLOT FILE 
C
C   GRAPH OF NPTS POINTS WITH NLIGN-ROWS, AND NPAGE-PAGES                      
C   (IF NLIGN=0, AUTOMATIC DETERMINATION OF NLIGN. NPAGE = 1 TO 8)              
C   WARNING: NPAGE MAY NOT BE MORE THAN 3 WITH MODE 1.                          
C   COORDINATES X(*) FOR HORIZONTAL AXIS JX, Y(*) FOR VERTICAL AXIS JY          
C   LABELS ARE IN ID(*), FORMAT A1 IF MOD=1, FORMAT A4 IF MOD=4               
C   POINTS AT MORE THAN PEX STANDARD DEVIATIONS ARE POSITIONED ON THE           
C   EDGES OF THE GRAPH (SUBROUTINE EPUR4).                                      
C   IF NOR=1 THE ORIGIN IS ADDED TO THE POINTS TO BE DISPLAYED.                 
C   WARNING:  X(*), Y(*), ID(NPTS+1) ARE DESTROYED UPON RETURN.                
C   GRAPH IS ABORTED IF MORE THAN 2048 POINTS ARE ON THE EDGES,                 
C   THE FIRST NHID HIDDEN POINTS ARE PRINTED WITH THEIR COORDINATES.         
C
C   CALLED BY:  SGRAF1 (IN SGRAP.F)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C **********************************************************************

       SUBROUTINE HISMAP(IDIM,NPTS,JX,JY,X,Y,ID,MOD,NLIGN,NPAGE,PEX,
     &                   NDAT,NBAND)     

       INTEGER, PARAMETER :: NHID = 200
       INTEGER            :: LD1(NHID),LD2(NHID), KLIC(371)
       REAL               :: EX(48)
       CHARACTER(LEN= 1)  :: NULL = CHAR(0)

       CHARACTER * 4   KLAC(371)
       CHARACTER * 4   LA(8),NA(8)
       CHARACTER * 1   MA(8)
       INTEGER         FPAGE,FL,FC
       REAL            X(IDIM),Y(IDIM),KKO,LLO,IREM,EY
       CHARACTER * 7   ID(IDIM)
       CHARACTER * 4   LDA,LDA2,LDB,LDB2
       CHARACTER * 4   IKLAC

       DATA MA/'-',' ','|','|','.','+',':','^'/
       DATA NA/'----','    ','|   ','   |',' .  ',' +  ',':   ',' ^  '/ 
       DATA LDA2/'    '/,LDB2/'    '/

       IF (NPAGE == 0) NPAGE = 1                                
       REWIND NBAND
  
       NPTSP1 = NPTS + 1   ! altered jan 2014 al                
       NPTSP1  = NPTS                     

C      CHECK TO SEE THAT NUMBER OF POINTS IS NOT EXCESSIVE
       IF (NPTSP1 > IDIM) THEN
         CALL ERRT(102,'HISMAP; NPTS EXCEEDS ARRAY DIMENSION',IDIM)
         RETURN
       ENDIF         

       WRITE (NDAT,1110)  NPTS, JX, JY, JX, JY
 1110  FORMAT (///' ',25X,'MAP OF'                 ,I6, '  POINTS' ,
     &   ' ON AXES'       ,I2,'  AND'  , I2,/' ',130('-')//,' ',21X,
     &   'AXIS',I2,'  /HORIZONTAL' ,10X, 'AXIS',I2,'  /VERTICAL' //)

C      FIND POINTS ON BOUNDARY OF MAP
 
       CALL EPUR4(IDIM,NPTS, X,Y,ID, MOD, PEX, KP, IDUM,IRTFLG,NDAT)
       IF (IRTFLG == 1)  RETURN     

C      SET SPECIAL SYMBOLS ACCORDING TO MODE
       DO K = 1,8                                                          
          IF (MOD == 1) THEN
             LA(K) = MA(K)
          ELSE
             LA(K) = NA(K)                           
          ENDIF
       ENDDO

       !ID(NPTSP1) = LA(6) ! jan 2014 al 
       !X(NPTSP1)  = 0.0   ! jan 2014 al
       !Y(NPTSP1)  = 0.0   ! jan 2014 al
       !CALL BORNS (NPTSP1,X,XMINT,XMAXT)
       !CALL BORNS (NPTSP1,Y,YMINT,YMAXT) 
             
C      FIND MIN/MAX
       XMINT = MINVAL(X(1:NPTSP1))
       XMAXT = MAXVAL(X(1:NPTSP1))
       YMINT = MINVAL(Y(1:NPTSP1))
       YMAXT = MAXVAL(Y(1:NPTSP1))

C      FIND NUMBER OF POSITIONS ON A LINE OF A PAGE (LENP)
       IF (MOD == 1) THEN
         LENP = 123    
       ELSEIF (MOD == 4) THEN
         LENP =  30   
       END IF

C      FIND TOTAL NUMBER OF POSITIONS ON A LINE                              
       LENTOT = LENP * NPAGE         
       NLINES = NLIGN
       FC     = LENTOT
       FPAGE  = NPAGE

C      FIND DEFAULT NUMBER OF LINES
       IF (NLINES == 0) THEN
         NLINES = ((YMAXT - YMINT) / (XMAXT - XMINT)) * FPAGE * 74.0
       END IF       

   30  IF (NLINES <= 12) NLINES = 12 
       FL    = NLINES
       S     = (XMAXT - XMINT) / FC    
       T     = (YMAXT - YMINT) / FL         
       NINT  = 5 * NPAGE + 1         
       ESPX  = FC / (5.0 * FPAGE)

       DO  J = 1,NINT      
C        FIND X AXIS LABELS

         EX(J) = XMINT + S * ESPX * FLOAT(J-1) 
       ENDDO

       KKO   = 0.50001 - (XMINT / S)
       LLO   = 0.50001 + ABS(YMAXT / T)

       DO   I = 1,NPTSP1 
C        SCALE EACH POINT TO FIT PAGE

         K     = (X(I) - XMINT) / S + 0.500001                   
         L     = (YMAXT - Y(I)) / T + 0.500001

C        PUT SCALED POINT IN X,Y
         IF (K == 0) K = 1 
         IF (L == 0) L = 1                                    
         X(I)  = K + 0.0001                             
         Y(I)  = L + 0.0001 
       ENDDO                  

       NSTACKED = 0
       DO   LL = 1,NLINES
C        PREPARE EACH LINE, EY IS Y AXIS LABEL

         EY    = YMAXT - T * FLOAT(LL-1)           

         DO  KK = 1,LENTOT
C          FILL OUTPUT LINE WITH BLANKS EXCEPT FOR FIRST, LAST, ETC. LINE
           KLIC(KK) = 0             
           KLAC(KK) = LA(2)                                   
           IF (LL == 1   .OR. LL == NLINES) KLAC(KK) = LA(1)      
           IF (KK == KKO .OR. LL == LLO)    KLAC(KK) = LA(5)
         ENDDO 

         KLAC(1)          = LA(3)        
         KLAC(LENTOT)     = LA(4)      
         KLAC(LENTOT + 1) = LA(2) 
         KLIC(LENTOT + 1) = 0  
         KLIC(LENTOT+2)   = 0 
  
         DO  I = 1,NPTSP1
C          PUT EACH POINT ID ON THE PAGE MAP, KLIC(K) HOLDS POINT ID

           L     = Y(I)  
           IF (L == LL) THEN
             K     = X(I) 
             IF (KLIC(K) == 0) THEN
C               PUT NUMBER OF POINT IN KLIC
                KLIC(K) = I  
             ELSE
C               ALREADY HAVE POINT AT THIS LOCATION
                IK        = KLIC(K)   
                NSTACKED  = NSTACKED + 1

C	        for pixel map, obviously, nhid=200 is not enough.  <ML
                IF ((MOD .NE. 1 .OR. NSTACKED > 100) .AND. 
     &                               NSTACKED <= NHID) THEN
                  LD1(NSTACKED) = IK 
                  LD2(NSTACKED) = I
                ENDIF
             ENDIF
           ENDIF
         ENDDO

         DO KK = 1,LENTOT
C          PUT CORRECT SYMBOLS ON THE OUTPUT LINE

           IK = KLIC(KK)                                      
           IF (IK .NE. 0) THEN
C             HAVE POINT AT THIS LOCATION
              IF (MOD == 1) THEN
                 KLAC(KK) = ID(IK)(1:1)
              ELSE
                 IF (ID(IK)(3:3) == ' ') THEN
                    KLAC(KK) = '  ' // ID(IK)(1:2)
                 ELSEIF (ID(IK)(4:4) == ' ') THEN
                    KLAC(KK) = ' ' // ID(IK)(1:3)
                 ELSE
                    KLAC(KK) =        ID(IK)(1:3)
                 ENDIF
              ENDIF
           ENDIF
         ENDDO

         LENTOTP1 = LENTOT + 1        

C        WRITE OUTPUT LINE TO TEMP FILE       
         WRITE (NBAND)  EY, (KLAC(K), K=1,LENTOTP1)                   
       ENDDO

       KB    = 0 
       KBB   = 1 
       IPUSH = 0

C      PRINT EACH PAGE
       DO   IPAGE = 1,NPAGE  
         REWIND NBAND          

         KA    = KB + 1   
         KB    = KA + LENP - 1                              
         IF (KB == LENTOT .AND. IPAGE > 1)  KB = LENTOT + 1 
         KAA   = KBB + 1 
         KBB   = KAA + 4 

         DO   LL = 1,NLINES 
C          PRINT EACH LINE ON THIS PAGE,  GET LINE FROM TEMPORARY FILE
           READ (NBAND)  EY, (KLAC(K), K=1,LENTOTP1)          

           IREM = LL - 5 * (LL / 5)

           IF (MOD == 1) THEN
C            SYMBOL IS ONLY ONE CHAR. LONG

             IF (IPAGE == 1) THEN
                WRITE(NDAT,1010) EY, (KLAC(K),K=KA,KB)
 1010           FORMAT (' ',F8.2,1X,123A1)                  
             ELSE
                WRITE(NDAT,1030) (KLAC(K),K=KA,KB)    
 1030           FORMAT (' ',130A1)
             ENDIF

           ELSE
C            SYMBOL IS 4 CHAR. LONG
             IF (IPAGE == 1) THEN
                WRITE (NDAT,1020) EY, (KLAC(K),K=KA,KB) 
 1020           FORMAT (' ',F8.2,1X,31A4)             
             ELSE
                WRITE (NDAT,1040) (   KLAC(K),K=KA,KB)    
 1040           FORMAT (' ',32A4)
             ENDIF
           ENDIF
  	 ENDDO

         IF (IPAGE == 1) THEN
C          PAGE ONE
           WRITE(NDAT,1050) (EX(K),K=1,6)     
 1050      FORMAT (' ',2X,5(F10.4,14X),F10.4)         

         ELSE IF (IPAGE > 1) THEN
           WRITE (NDAT,1060) (EX(K),K=KAA,KBB)
 1060      FORMAT (' ',18X,4(F10.4,14X),F10.4) 
        ENDIF

         IF (IPAGE < NPAGE) THEN
C           NEW PAGE
            WRITE (NDAT,90) 
90          FORMAT(/////)
         ENDIF
       ENDDO

C      CHECK IF THERE ARE STACKED UP POINTS
       IF (NSTACKED == 0 .OR. MOD == 1) RETURN

C      FOR STACKED CASES ONLY ---------------------------------------                             

C      GIVE DATA ON HIDDEN (STACKED) POINTS
       WRITE (NDAT,1070) NHID
 1070  FORMAT(//,' ',50X,'STACKED POINTS   ( MAXIMUM OF ',I3,' )/',/,    
     &       ' ',2(' POINT',5X,'*',4X,'ABSCISSA',9X,'*',6X,'ORDINATE',
     &     7X,'*',10X)/,
     &     '   VISABLE',2X,'*',4X,'APPROX.',8X,'*',6X,'APPROX.',6X,
     &     '*',11X,'HIDDEN',4X,'*',4X,'APPROX.',8X,'*',6X,
     &     'APPROX.',6X,'*')

C      DRAW LINE ACROSS PAGE
       WRITE (NDAT,1080)
 1080  FORMAT (' ',2(28(' *'),10X) )

       LTOT = MIN0(NSTACKED,NHID)     
       DO  L = 1, LTOT
         IK    = LD1(L)      
         I     = LD2(L)      
         LDA   = ID(IK)    
         LDB   = ID(I)
         XD1   = X(IK) * S + XMINT                     
         YD1   = YMAXT - Y(IK) * T          
         XD2   = X(I) * S + XMINT           
         YD2   = YMAXT - Y(I) * T           
         WRITE(NDAT,1090) LDA,LDA2,XD1,YD1,LDB,LDB2,XD2,YD2   
 1090    FORMAT(' ',2(' *' ,2A4,' *',2(F13.3,8X,'*'),10X)  ) 
       ENDDO          

C      DRAW LINE ACROSS PAGE
       WRITE(NDAT,1080)

       WRITE(NDAT,1100) NSTACKED        
 1100  FORMAT(/ ' ',' NUMBER OF STACKED POINTS = ',I6, /)

       END                                                                     
