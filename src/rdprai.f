
C++*********************************************************************
C
C RDPRAI.F  -- CREATED JAN 87   AUTHOR:   ARDEAN LEITH
C              COMMENTS ALLOWED   JAN 98    ARDEAN LEITH
C              USED RDPRANC       AUG 99    ARDEAN LEITH
C              INTEGER RDPRANC    MAY 02    ARDEAN LEITH
C              RDPRANC PARAMETERS MAY 07    ARDEAN LEITH
C                                                                      *
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
C    RDPRAI(ILIST,NMAX,NUMBER,ILOW,IHI,PROMPT,LIMITEDC,IRTFLG)
C
C    PURPOSE:  READ AN ALPHANUMERIC STRING, CHECK FOR ANY SPECIAL OPERATION,
C              RETURNS ARRAY OF INTEGERS, THE NUMBER OF INTEGERS, OR A
C              FLAG  TO INDICATE THAT ONE SHOULD RETURN TO PREVIOUS QUESTION.
C              CHECKS TO SEE THAT INTEGERS ARE WITHIN SPECIFIED RANGE.
C              (CURRENTLY NO CHECKING IS DONE!!!!)
C
C              ALLOWABLE STRINGS           NUMBERS ENTERED IN ARRAY
C                 1,2,8,66                  1, 2, 8, 66
C                 1 2 8,66                  1, 2, 8, 66
C                 1,2-8 66                  1, 2, 3, 4, 5, 6, 7, 8, 66
C                 X11,X12....              CONTENTS OF X11 & X12...
C                 X11-X12                  CONTENTS OF X11....X12
C
C              (A - ONLY SIGNIFIES A RANGE IF IT IMMEDIATELY FOLLOWS
C               PRECEEDING VALUE)
C
C  PARAMETERS : ILIST     ARRAY FOR ANSWERS                      (RET.)
C               NMAX      MAX LENGTH OF ARRAY                    (SENT)
C               NUMBER    ON ENTRY IS MAX NUMBER OF ANSWERS (SENT/RET.) 
C                         TO BE RETURNED!!
C                         (<0 ON ENTRY IS FLAG TO ACCEPT NULL 
C                         RESPONSE)
C                         ON RETURN IS NUMBER OF ANSWERS ACTUALLY 
C                         RETURNED
C               ILOW      BOTTEM OF RANGE (CAN'T BE BELOW THIS)   (SENT)
C               HI        TOP OF RANGE (CAN'T BE ABOVE THIS)      (SENT)
C                         ILOW & IHI CURRENTLY UNUSED!!!!!!!!!!!!!!
C               PROMPT    SOLICITATION MESSAGE                    (SENT)
C               LIMITEDC  REQUIRES NMAX VALUES / CALL, SO THAT
C                           IT CAN BE USED IN BATCH DO-LOOP
C               IRTFLG    RETURN FLAG (0 IS NORMAL, -1 IS GOTO    (SENT)
C                              PREVIOUS QUESTION)
C      
C    NOTE: CAN PERFORM REGISTER SUBSTITUTION NOW al.
C          IF THE LAST CHAR. ON A LINE OF INPUT IS A COMMA IT WILL ASK 
C          FOR ANOTHER LINE OF INPUT.
C
C          BE SURE NUMBER IS SET TO NMAX BEFORE CALLING THIS ROUTINE
C
C--*********************************************************************

        SUBROUTINE RDPRAI(ILIST,NMAX,NUMBER,ILOW,IHI,PROMPT,
     &             LIMITEDC,IRTFLG)        

        INCLUDE 'CMBLOCK.INC'

        DIMENSION       ILIST(NMAX)
        CHARACTER *(*)  PROMPT,LIMITEDC
        LOGICAL         LIMITED

        LIMITED = (LIMITEDC .EQ. 'T')

        CALL RDPRANC(ILIST,NMAX,NUMBER,ILOW,IHI,PROMPT,
     &               LIMITED,IRTFLG)        

        RETURN
        END
       
