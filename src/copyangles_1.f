C++*********************************************************************
C
C COPYANGLES_1.F          
C                   
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
C  COPYANGLES_1(LUN1,LUN2,NSAM1,NSAM2)
C
C  PURPOSE:  COPIES ANGLES FROM HEADER OF INPUT FILE (HEADER POSITIONS,
C            (14-17) TO HEADER OF OUTPUT FILE 
C  LUN1      LOGICAL UNIT FOR INPUT FILE              (SENT)
C  LUN       LOGICAL UNIT FOR OUTPUT FILE             (SENT)
C  NSAM1     INPUT   IMAGE SIZE                       (SENT)
C  NSAM2     OUTPUT  IMAGE SIZE                       (SENT)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE COPYANGLES_1(LUN1,LUN2,NSAM1,NSAM2)

        DIMENSION      VALUE(60), BUF1(1000)
        
        ILOC1 = 14
        CALL GETLAB(LUN1,NSAM1,BUF1,ILOC1,4,VALUE,IRTFLG)
        ILOC1 = 14
        CALL SETLAB(LUN2,NSAM2,BUF1,ILOC1,4,VALUE,'U',IRTFLG)
        RETURN
        END
