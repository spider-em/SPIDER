
C***********************************************************************
C
C  DENDRO.F -- ADAPTED FOR METAFILE             3 NOV 1986 ARDEAN LEITH
C              USED POSTSCRIPT OUTPUT             MAR 1999 ARDEAN LEITH
C              COSMETIC OUTPUT CHANGES            DEC 2008 ARDEAN LEITH                                                *
C              MERGED WITH DENDRO2 & ARBRE        JUN 2009 ARDEAN LEITH                                                *
C              SCALE CHANGED                      JUN 2009 ARDEAN LEITH                                                *
C              INPUT ORDER CHANGED                JUN 2009 ARDEAN LEITH                                                *
C              FILPOS * BUG                       NOV 2009 ARDEAN LEITH                                                *
C**********************************************************************
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
C PURPOSE: DRAWS DENDROGRAM AND FORMS A POSTSCRIPT PLOT FILE FOR IT
C
C ORIGIN:
C    BASED ON ARBRE.F,  A DENDROGRAM PGM BY - JEAN-PIERRE BRETAUDIERE 
C    THE UNIVERSITY OF TEXAS HEALTH SCIENCE CENTER AT HOUSTON                      
C    MEDICAL SCHOOL - DEPARTMENT OF PATHOLOGY AND LABORATORY MEDICINE              
C    P.O. BOX 20708, HOUSTON, TX 77225.                                            
C                                                                 
C    TREE DESCRIPTION 
C                                                                  
C       NKLA SUMMITS              JFIN = 2 * NKLA - 1                
C       ARRAYS PROVIDED BY CHAVA  VAL(JFIN), LA(NKLA), LB(NKLA)     
C                                 PK(JFIN)      
C                                                               
C       WORKING ARRAYS            NO(JFIN), HIT(NKLA), NT(NKLA),     
C                                 IW(NKLA), W(NKLA), IV(NKLA),    
C                                 NUM(NKLA)                          
C
C       CALLED BY:     NOYAU
C
C **********************************************************************

      SUBROUTINE DENDRO(NKLA, JFIN, VAL, LA, LB, PK, IDK, KLAS, NUMIM, 
     &                  IDI, ADDHIDDEN, NO,NUM,NT,IV,IW,V,W)

      INCLUDE 'CMLIMIT.INC' 
      INCLUDE 'CMBLOCK.INC' 

      INTEGER, PARAMETER    :: LUNPOS = 80
      INTEGER, PARAMETER    :: LUNDOC = 81
      INTEGER, PARAMETER    :: NSIZE  = 3
      REAL                  :: DATA(3,NSIZE)
      
      REAL                  :: VAL(JFIN),PK(JFIN),HIT(NKLA),W(NKLA)
      INTEGER               :: LA(NKLA),LB(NKLA),IDK(NKLA)
      INTEGER               :: NO(JFIN),NT(NKLA)
      INTEGER               :: IV(NKLA),NUM(NKLA),IW(NKLA)
      INTEGER               :: KLAS(NUMIM),IDI(NUMIM)

      REAL                  :: Y1(MAX(NKLA,NUMIM))
      CHARACTER(LEN=3)      :: CWGHT1(MAX(NKLA,NUMIM))
      INTEGER               :: IM1(MAX(NKLA,NUMIM))
      INTEGER               :: ICLAS1(MAX(NKLA,NUMIM))

      CHARACTER(LEN=1)      :: NULL
      CHARACTER(LEN=8)      :: LABEL
      CHARACTER(LEN=10)     :: CVMIN,CVMAX
      CHARACTER(LEN=MAXNAM) :: FILPOS,ANSW
      LOGICAL               :: ASKTRUNC,ADDHIDDEN,WANTPOS

      INTEGER, PARAMETER    :: LENDOT = 102
      CHARACTER(LEN=LENDOT) :: LDOT

       NULL = CHAR(0)

       CALL FILERD(FILPOS,NLET,NULL,'DENDROGRAM POSTSCRIPT OUTPUT~9',
     &             IRTFLG)

       WANTPOS  = (FILPOS(1:1) .NE. '*')
       IF (WANTPOS .AND. IRTFLG .NE. 0) RETURN

       ANSW = FILPOS(1:NLET)
       CALL SSUPCAS(ANSW)

       ASKTRUNC = .FALSE.
       PMIN     = 0.0

C      LONG KLUDGE FOR BACKWARDS INPUT COMPATIBILITY
       IF ((NLET .EQ. 1 .AND. ANSW(1:1) .EQ. 'Y') .OR.
     &     (NLET .EQ. 3 .AND. ANSW(1:3) .EQ. 'YES')) THEN
C         WANT UN-TRUNCATED OUTPUT FILE, ASK FOR FILE NAME AGAIN
          CALL FILERD(FILPOS,NLET,NULL,
     &                'DENDROGRAM POSTSCRIPT OUTPUT~9',IRTFLG)
          
       ELSEIF ((NLET .EQ. 1 .AND. ANSW(1:1) .EQ. 'N') .OR.
     &         (NLET .EQ. 2 .AND. ANSW(1:2) .EQ. 'NO')) THEN
C          DO NOT WANT POSTSCRIPT OUTPUT FILE, BUT NEED RESULTS PLOT
           WANTPOS  = .FALSE.

       ELSEIF ((NLET .EQ. 1 .AND. ANSW(1:1).EQ. 'T') ) THEN
C          OLD INPUT FORMAT, ASK FOR TRUNCATION  LEVEL HERE
           CALL RDPRM1S(PMIN,NOT_USED,'PLOT CUTOFF % (0 .. 100)',
     &                IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          ASK FOR FILE NAME AGAIN
           CALL FILERD(FILPOS,NLET,NULL,
     &                 'DENDROGRAM POSTSCRIPT OUTPUT~9',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
       ELSEIF (WANTPOS) THEN
            ASKTRUNC = .TRUE.
       ENDIF

       IF (WANTPOS) THEN
C         OPEN POSTSCRIPT FILE AS SEQUENTIAL FORMATTED
          CALL OPAUXFILE(.FALSE.,FILPOS,NULL,LUNPOS,0,'N',
     &                ' ',.TRUE.,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN
          NLETP = LNBLNKN(FILPOS)
       ENDIF

       JDEB   = NKLA + 1     
       VMIN   = HUGE(VAL)    
       VMAX   = -VMIN          

       DO J = JDEB,JFIN                                                   
          IR     = J - JDEB + 1           
          IA     = LA(IR)        
          IB     = LB(IR)         
          VMIN   = MIN(VAL(J),VMIN)
          VMAX   = MAX(VAL(J),VMAX)
          NO(IA) = J          
          NO(IB) = J         
       ENDDO
       NO(JFIN) = JFIN  

C      RESCALE THE INDEX (HEIGHT) LEVELS 0..SCLMAX
       SCLMAX = 100.0
       VSCAL  = SCLMAX / (VMAX - VMIN)

       DO I = JDEB,JFIN
          VAL(I) = (VAL(I) - VMIN) * VSCAL
       ENDDO

      IF (ASKTRUNC) THEN
C        NEW INPUT FORMAT, ASK FOR TRUNCATION  LEVEL HERE
         CALL RDPRM1S(PMIN,NOT_USED,'PLOT CUTOFF (0 .. 100)',
     &                IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
      ENDIF
C     INTEGERIZE CUT-OFFS
      IPMIN = PMIN   + 0.5       
      IVMAX = SCLMAX + 0.5

C     CREATE DESCRIPTION OF THE HIERARCHY CLASSES 
      WRITE(NDAT,90)     
   90 FORMAT (/,' NODE  INDEX  SENIOR  JUNIOR  SIZE'    ,6X,  
     &         'DESCRIPTION OF HIERARCHY CLASSES',/,
     &         ' ',64('. '),/)

      DO J = JDEB,JFIN         
        NT(1) = J 
        KPT   = 0 
        JI    = 1 

        DO WHILE (JI .NE. 0)
           IF (NT(JI) .LE. NKLA) THEN    
              K               = NT(JI)   
              KPT             = KPT + 1  
              IV(KPT)         = IDK(K)
              NUM(NKLA-KPT+1) = IDK(K)
              IW(KPT)         = K       
              JI              = JI - 1 

              !iit = NKLA-KPT+1
              !write(6,9900) iit,num(iit), j,kpt,ji,k,iv(kpt),iw(kpt) 
9900          format(' NUM(',i4,'):',i4,'  ----',6i6)
        
           ELSE 
              IJ              = JI + 1         
              NI              = NT(JI) - JDEB + 1
              NT(IJ)          = LA(NI)           
              NT(JI)          = LB(NI)           
              JI              = JI + 1
           ENDIF
        ENDDO
     
C       PRINT DESCRIPTION OF THE HIERARCHY CLASSES 
        IR = J - JDEB + 1     
        WRITE(NDAT,610) J, VAL(J),LA(IR),LB(IR),KPT,
     &		        (MOD(IV(KK),10000),KK=1,KPT)
  610   FORMAT (1X,I4,F9.3,I5,I6,I7,5X,
     &               18(1X,I4,:),/,:,
     &       ' ',36X,18(1X,I4,:),/,:,
     &       ' ',36X,18(1X,I4,:),/,:,
     &       ' ',36X,18(1X,I4,:),/,:,
     &       ' ',36X,18(1X,I4,:))

        I1    = IW(1)            
        I2    = IW(KPT)          
        W(I1) = NO(J)  + 0.0001   
        W(I2) = NO(I2) + 0.0001  
      ENDDO

C     CALCULATE TRUNCATED DENDROGRAM  --------------------------------

C     COMPUTE INDEX FOR EACH ORIGINAL CLASS
      DO J = 1,NKLA
         K    = W(J)      
         W(J) = VAL(K)
      ENDDO

      DO J = 1,NKLA   
         IWI           = IW(J)     
         HIT(NKLA-J+1) = W(IWI)     ! HEIGHT OR INDEX
      ENDDO

      NBR   = 0        ! NUMBER OF LEAVES OR NODES
 
      DO J = 1,NKLA
         !write(6,*) ' j,hit(j),pmin:',j,hit(j),pmin

C        DO NOT LABEL TRUNCATED LINES
         IF (HIT(J) .LT. PMIN)  CYCLE


C        GET CURRENT NODE WEIGHT 
         JJ      = IW(NKLA-J+1)
         IWGHT   = PK(JJ)

         IF (ADDHIDDEN) THEN
C           ADD UNLABELD CLASS LEAVES 

            DO I=1, NUMIM
               IF (KLAS(I) .EQ. NUM(j)) THEN
C                 IMAGE (OR ID) IS IN THIS CLASS
                  NBR     = NBR + 1
                   
C                 SET  LEAF HEIGHT (Y)
                  Y1(NBR)  = PMIN               ! HEIGHT 

C                 SET LEAF ID
                  IM1(NBR) = IDI(I)            ! ID  

C                 SET LEAF CLASS
                  ICLAS1(NBR) = NUM(J)         ! CLASS NUMBER

C                 SET LEAF WEIGHT
                  CWGHT1(NBR) = '  1'
               ENDIF
            ENDDO

         ELSE
C           NO NEED TO ADD UNLABELED CLASSES  
            NBR     = NBR + 1

C           SET LEAF WEIGHT
            WRITE(CWGHT1(NBR),FMT='(I3)') IWGHT

C           SET LEAF CLASS 
            ICLAS1(NBR) = NUM(J)            !  CLASS NUMBER

C           SET LEAF ID
            IM1(NBR) = NUM(J)               ! ID  
         ENDIF

C        SET (OR RESET) LEAF HEIGHT (Y) 
         Y1(NBR) = HIT(J)  
       
      ENDDO	

C     DRAW DENDROGRAM  ------------------------------------------

C     SET WINDOW SIZE (SCALING WILL BE DONE IN PLOT ROUTINES)
c     TO ALTER X SPACING OF NODES

      IF (NBR .LT. 20) THEN
         XWIN    = 100     
         YWIN    = 120
         ITSIZEX = 9
      ELSEIF (NBR .LT. 100) THEN
         XWIN    = 200    
         YWIN    = 220
         ITSIZEX = 5
      ELSE
         XWIN    = 300     
         YWIN    = 270
         ITSIZEX = 5
      ENDIF

C     FIND Y PLOT SCALE FOR TREE (NOT INCLUDING LABELS)
      YSCALE = YWIN / (SCLMAX - PMIN)
  
C     SET HORIZ DISTANCE BETWEEN LEAVES
      XDIS = XWIN / (NBR - 1)

C     ADDED LENGTH FOR LEAVES
      YLEAF = -5.0

C     INITIALIZE & SET SCALING FOR POSTSCRIPT
      IF (WANTPOS) CALL POSTRT(-LUNPOS)

      XLL = -66.0  
      YLL = -26.0
      XUR = XWIN
      YUR = YWIN + 5
      IF (WANTPOS) CALL POSCALE(LUNPOS, 1.0,1.0,  XLL,YLL,  XUR,YUR)

C     SET TEXT CHARACTARISTICS FOR Y AXIS LABELS
      ITSIZE = 9
      ITANGL = 0
      JUST   = 0

C     PUT POSTSCRIPT FILENAME AT TOP
      YPOS = YWIN
      XPOS = 20.0
      IF (WANTPOS)
     &    CALL POTEX(LUNPOS,FILPOS,NLETP,XPOS,YPOS,ITSIZE,ITANGL,JUST)

C     RIGHT JUSTIFIED Y LABELS
      JUST = 2

C     LABEL Y AXIS
      XPOS  = -16.0
      YPOS  = -30.0
      LABEL = 'ID'
      IF (WANTPOS .AND. ADDHIDDEN) 
     &   CALL POTEX(LUNPOS,LABEL,2,XPOS,YPOS, ITSIZE,ITANGL,JUST)

C     LABEL FOR WEIGHT
      XPOS  = -16.0
      YPOS  = -24.0
      LABEL = 'WEIGHT'
      IF (WANTPOS) 
     &   CALL POTEX(LUNPOS,LABEL,6,XPOS,YPOS, ITSIZE,ITANGL,JUST)

C     LABEL FOR INDEX
      YPOS = -18.0
      LABEL = 'INDEX '
      IF (WANTPOS) 
     &   CALL POTEX(LUNPOS,LABEL,6,XPOS,YPOS, ITSIZE,ITANGL,JUST)

C     LABEL FOR CLASS
      YPOS = -12.0
      LABEL = 'CLASS '
      IF (WANTPOS) 
     &   CALL POTEX(LUNPOS,LABEL,6,XPOS,YPOS, ITSIZE,ITANGL,JUST)

C     LABEL FOR SCALE RANGE
      YPOS = 0.0   
      IF (WANTPOS) THEN
         WRITE(CVMIN,FMT='(I4)') IPMIN
         WRITE(CVMAX,FMT='(I4)') IVMAX
         CALL POTEX(LUNPOS,CVMIN,4, XPOS,0.0,  ITSIZE,ITANGL,JUST)
         CALL POTEX(LUNPOS,CVMAX,4, XPOS,YWIN, ITSIZE,ITANGL,JUST)
      ENDIF

      XPOS   = -20.0
      YPOS   = 50.0
      ITANGL = 90
      LABEL   = 'SCALE '

      IF (WANTPOS) THEN
         CALL POTEX(LUNPOS,LABEL,6,XPOS,YPOS, ITSIZE,ITANGL,JUST)
C        ADD TICK MARKS AT Y = 0 AND Y = YMAX
         CALL POSEG(LUNPOS, -9.0, 0.0,  -4.0, 0.0)
         CALL POSEG(LUNPOS, -9.0,YWIN,  -4.0,YWIN)
      ENDIF

C     PRINT LABEL FOR DENDROGRAM IN RESULTS FILE

      WRITE(NDAT,*)'  '
      WRITE(NDAT,*)' xxxxxxxxxxxxxxxxxxxxxxxx DENDROGRAM xxxxxxxxxxx'//
     &             'xxxxxxxxxxxxxxxxxxxxxxxxx'
      WRITE(NDAT,*)'  '

      IF (ADDHIDDEN) THEN
         WRITE(NDAT,620) IPMIN, IVMAX
  620    FORMAT (///,'    ID   INDEX   CLASS           DENDROGRAM    ',    
     &           '(SCALE:  ',I4,'.....',I4,' )',//)                 
      ELSE
         WRITE(NDAT,621) IPMIN, IVMAX
  621    FORMAT (///,'  WEIGHT INDEX  CLASS          DENDROGRAM    ',    
     &               '(SCALE:  ',I4,'.....',I4,' )',//)                 
      ENDIF

C     SET TEXT CHARACTARISTICS FOR X AXIS LABELS
      ITANGL         = 0
      JUST           = 1    ! CENTERED

      LDOT(1:LENDOT) = ' '   ! BLANK THE RESULTS FILE DOTTED LINE
      DOTSCL         = FLOAT(LENDOT-1) / ( SCLMAX - PMIN )

C     DRAW  NBR TREE BRANCHES INTO POSTSCRIPT FILE 
      DO I=1,NBR

         XLOC = I * XDIS
         YLOC = (Y1(I) - PMIN) * YSCALE

C        DRAW VERTICAL LINE
         DATA(1,1) = XLOC
         DATA(2,1) = YLEAF

         DATA(1,2) = XLOC
         DATA(2,2) = YLOC 
         NDATA  = 2

         !write(6,*) ' i,im1,y1:',i,im1(i),y1(i)

         IF (I .LT. NBR) THEN
C           FIND LENGTH OF HORIZONTAL LINE
            DO J = I+1,NBR
               IF (((Y1(J) - PMIN) * YSCALE) .GE. YLOC) THEN
C                 ADD HORIZONTAL LINE TO THIS BRANCH
                  DATA(1,3) = J * XDIS
                  DATA(2,3) = YLOC   
                  NDATA     = 3
                  EXIT
               ENDIF
  	    ENDDO
            ! IF NO HIGHER BRANCH FOUND, POSSIBLE ERROR 
         ENDIF

C        PUSH  LINES INTO POSTSCRIPT FILE
         IF (WANTPOS) 
     &      CALL POARAYF(LUNPOS,DATA,NDATA,.FALSE.,.FALSE.)

         IF (WANTPOS .AND. ADDHIDDEN) THEN
C          LABEL LEAF WITH ID
           YPOS  = -30.0
           WRITE(LABEL,FMT='(I4)') IM1(I)
           CALL POTEX(LUNPOS,LABEL,4,XLOC,YPOS, ITSIZEX,ITANGL,JUST)
         ENDIF

C        LABEL LEAF WITH WEIGHT
         YPOS  = -24.0
         IF (WANTPOS) 
     &     CALL POTEX(LUNPOS,CWGHT1(I),3,XLOC,YPOS, ITSIZEX,ITANGL,JUST)

C        LABEL LEAF WITH INDEX 
         YPOS  = -18.0
         INDX  = Y1(I) + 0.5
         WRITE(LABEL,FMT='(I3)') INDX
         IF (WANTPOS) 
     &       CALL POTEX(LUNPOS,LABEL,3,XLOC,YPOS, ITSIZEX,ITANGL,JUST)

C        LABEL LEAF WITH CLASS NUMBER
         YPOS = -12.0
         WRITE(LABEL(1:3),FMT='(I3)')ICLAS1(I)
         IF (WANTPOS) 
     &     CALL POTEX(LUNPOS,LABEL(:3),3,XLOC,YPOS, ITSIZEX,ITANGL,JUST)

         !write(6,*) ' i: ',i,im1(i),indx,xloc,iclas1(i)

C        CREATE DOTTED LINE FOR RESULTS FILE
         FLX        = (Y1(I) - PMIN)  * DOTSCL + 1.0
         LX         = FLX 
         DO IT = 1,LX
            LDOT(IT:IT) = '.'
         ENDDO
         LT = lnblnkn(LDOT)

C        WRITE LEAF IN RESULTS FILE
         IF (ADDHIDDEN) THEN
            WRITE(NDAT,96) IM1(I),    INDX,  ICLAS1(I),  LDOT(1:LT)
  96        FORMAT (' ',   I4,4X,     I4,5X, I4,'  ..',  A)
         ELSE 
            WRITE(NDAT,95) CWGHT1(I), INDX,  ICLAS1(I),  LDOT(1:LT)
  95        FORMAT (' ',   A,4X,      I4,5X, I4,'  ..',  A) 
         ENDIF 
C        BLANK THE LOWER PART TO MAKE HORIZONTAL LINE ON NEXT WRITE 
         IF (LX .NE. 1)  LDOT(1:LX-1) = ' '

      ENDDO
      WRITE(NOUT,*) ' '

C     CLOSE THE POSTSCRIPT-FILE 
      IF (WANTPOS) THEN
          CALL POEND(LUNPOS)
          WRITE(NOUT,*) ' PLOT PLACED IN: ',FILPOS(1:NLETP)
          CLOSE(LUNPOS)
      ENDIF

C     CALCULATE UNTRUNCATED DENDROGRAM  ------------------------------

      NBR   = 0        ! NUMBER OF LEAVES OR NODES

      DO J = 1,NKLA
         !write(6,*) ' j,hit(j),pmin:',j,hit(j),pmin

         IF (ADDHIDDEN) THEN
C           ADD UNLABELD CLASS LEAVES 

            DO I=1, NUMIM
               IF (KLAS(I) .EQ. NUM(J)) THEN
C                 IMAGE (OR ID) IS IN THIS CLASS
                  NBR         = NBR + 1
                   
                  ICLAS1(NBR) = NUM(J)       ! CLASS NUMBER
                  Y1(NBR)     = 0.0          ! HEIGHT 
                  IM1(NBR)    = IDI(I)       ! ID  
               ENDIF
            ENDDO
            Y1(NBR) = HIT(J)                 ! HEIGHT 
         ELSE
C           NO NEED TO ADD UNLABELED CLASSES  
            NBR         = NBR + 1

            ICLAS1(NBR) = NUM(J)             ! CLASS NUMBER
            Y1(NBR)     = HIT(J)             ! HEIGHT 
            IM1(NBR)    = NUM(J)             ! ID  
         ENDIF
      ENDDO	

C     PLACE UNTRUNCATED DENDROGRAM SPECIFICATIONS IN DOC FILE
      CALL DENLST(LUNDOC,NBR, Y1,ICLAS1,IM1, IRTFLG)

      END

