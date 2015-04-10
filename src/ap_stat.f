 
C++*********************************************************************
C
C    AP_STAT.F    NEW                               MAR 05 ARDEAN LEITH
C                 CCDIF = CCROT-CCROTLAS            AUG 10 ARDEAN LEITH
C                 NBORDER,NSUBPIX                   OCT 10 ARDEAN LEITH
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
C  AP_STAT(NIDI,ANGDIFTHR,IBIGANGDIF, ANGDIFAVG, CCROTAVG,
C           IMPROVCCROT,CCROTIMPROV,IWORSECCROT, CCROTWORSE,LUNDOC)
C
C  PURPOSE: ACCUMULATE AND LIST REFINEMENT STATISTICS IN DOC FILE
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       SUBROUTINE AP_STAT_ADD(NGOTPAR,CCROT,ANGDIF,ANGDIFTHR,CCROTLAS,
     &                       CCROTAVG,IBIGANGDIF,ANGDIFAVG,IMPROVCCROT,
     &                       CCROTIMPROV,IWORSECCROT,CCROTWORSE)

        IF (NGOTPAR .LT. 0) THEN
C          INITIALIZE CCROT CHANGE STATISTICS

           CCROTAVG    = 0.0
           IBIGANGDIF  = 0
           ANGDIFAVG   = 0.0

           IMPROVCCROT = 0
           CCROTIMPROV = 0.0
           IWORSECCROT = 0
           CCROTWORSE  = 0.0
           RETURN
        ENDIF  

C       COMPILE CCROT CHANGE STATISTICS

        CCROTAVG = CCROTAVG + CCROT

        IF (NGOTPAR .GE. 8) THEN
           IF (ANGDIF .GT. ANGDIFTHR) IBIGANGDIF = IBIGANGDIF + 1
           ANGDIFAVG = ANGDIFAVG + ANGDIF
           CCDIF     = CCROT - CCROTLAS     ! DIFFERENCE

           IF (CCDIF .GE. 0.0) THEN
              IMPROVCCROT = IMPROVCCROT + 1
              CCROTIMPROV = CCROTIMPROV + CCDIF
           ELSE
              IWORSECCROT = IWORSECCROT + 1
              CCROTWORSE  = CCROTWORSE + ABS(CCDIF)
           ENDIF
        ENDIF   ! END OF: IF (NGOTPAR .GE. 8)

        END

C       -------------------- AP_STAT -------------------------------

       SUBROUTINE AP_STAT(NUMEXP,ANGDIFTHR,IBIGANGDIF,
     &                    ANGDIFAVG, CCROTAVG,
     &                    IMPROVCCROT,CCROTIMPROV,
     &                    IWORSECCROT, CCROTWORSE,
     &                    NBORDER,NSUBPIX,LUNDOC)
       
        REAL, DIMENSION(8)  :: DLIST
        CHARACTER (LEN=118) :: COMMENT 


C       SAVE CCROT & ANGULAR DISPLACEMENT STATISTICS
        FNUMEXP  = NUMEXP
        DLIST(1) = 0.0
        IF (ANGDIFTHR .GT. 0)DLIST(1) = 100 * FLOAT(IBIGANGDIF)/FNUMEXP

        DLIST(2) = ANGDIFAVG / FNUMEXP
        DLIST(3) = CCROTAVG  / FNUMEXP
        DLIST(4) = 100 * FLOAT(IWORSECCROT) / FNUMEXP

        DLIST(5) = 0.0
        IF (IWORSECCROT .GT. 0)DLIST(5) = CCROTWORSE/FLOAT(IWORSECCROT)

        DLIST(6) = 0.0
        IF (IMPROVCCROT .GT. 0)DLIST(6) = CCROTIMPROV/FLOAT(IMPROVCCROT)

        IF (NBORDER > 0 .OR. NSUBPIX > 0) THEN

           DLIST(7) = NBORDER
           DLIST(8) = NSUBPIX

C            123456789 123456789 123456789 123456789 123456789 123456789 1
           COMMENT='      '//
     &      '  %BIG-ANGDIF,   AVG-ANGDIF,     AVG-CCROT,     %WORSE, '//
     &      '    AVG-WORSE,    AVG-BETTER       #BORDER,    #SUB_PXL'
           CALL LUNDOCPUTCOM(LUNDOC,COMMENT,IRTFLG)
           CALL LUNDOCWRTDAT(LUNDOC,-2,DLIST,8,IRTFLG)
      
        ELSE
C                123456789 123456789 123456789 123456789 123456789 123456789 1
           COMMENT='      '// 
     &          '%BIG-ANGDIF,   AVG-ANGDIF,     AVG-CCROT,   %WORSE,'//
     &          '  AVG-WORSE,   AVG-BETTER'
           CALL LUNDOCPUTCOM(LUNDOC,COMMENT,IRTFLG)
           CALL LUNDOCWRTDAT(LUNDOC,-2,DLIST,6,IRTFLG)
        ENDIF

        CLOSE(LUNDOC)

        END
