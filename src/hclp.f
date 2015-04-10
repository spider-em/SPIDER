
C++*********************************************************************
C
C HCLP.F  
C                DOC FILE *                        MAY 09 ARDEAN LEITH
C                DENDRO REWRITE                    JUN 09 ARDEAN LEITH
C                MDIM_8                            MAY 13 ARDEAN LEITH
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
C   HCLP(NKLA,MDIM_8,JFIN,WD,INUM,COO,COB,MAXFAC,NFAC,
C        D,PK,VAL,IDK,LA,LB,NT,NO,NUM,IV,IW,V,W)
C        LUNF,LUNT,LUNDOC,MODE,ITYPE)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

         SUBROUTINE HCLP(NKLA,MDIM_8,JFIN,WD,INUM,MAXFAC,NFAC,
     &                    D,PK,VAL,LA,LB,NT,NO,IV,IW,V,W,
     &                    LUNF,LUNT,LUNDOC,MODE,ITYPE)

         INCLUDE 'CMLIMIT.INC' 
         INCLUDE 'CMBLOCK.INC'

         INTEGER *8      :: MDIM_8
         REAL            :: D(MDIM_8),PK(JFIN),VAL(JFIN)
         INTEGER         :: NT(JFIN),NO(JFIN),LA(NKLA),LB(NKLA)

         INTEGER         :: IV(NKLA),IW(NKLA)
         REAL            :: V(NKLA),W(NKLA)
         REAL            :: WD(MAXFAC)
         INTEGER         :: INUM(MAXFAC)

C        AUTOMATIC ARRAYS
         INTEGER         :: NUM(NKLA),IDK(NKLA)

C        INPUT DATA
         CALL DIST_P(D,MDIM_8,IDK,NKLA,WD,INUM,MAXFAC,NFAC,LUNF,ITYPE)

C        SET PK ARRAY
         PK = 1.0

         OPEN(LUNT,STATUS='SCRATCH',FORM='UNFORMATTED')

         CALL CHAVA(NKLA, MDIM_8, JFIN,
     &             D,  PK, VAL, LA, LB, NT, NO, LUNT, MODE)

C        CLASSIFICATION TREE OF THE NKLA CENTERS
         KDUM = 1
         CALL DENDRO(NKLA, JFIN, VAL, LA, LB, PK, IDK, KDUM,KDUM,
     &               KDUM, .FALSE., NO,NUM,NT,IV,IW,V,W)

        CLOSE(LUNT)

        END
