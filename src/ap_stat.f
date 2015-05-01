 
C++*********************************************************************
C
C    AP_STAT.F    NEW                               MAR 05 ARDEAN LEITH
C                 CCDIF = CCROT-CCROTLAS            AUG 10 ARDEAN LEITH
C                 NBORDER,NSUBPIX                   OCT 10 ARDEAN LEITH
C                 FNUMPIX                           APR 15 ARDEAN LEITH
C                 AP_STAT_R ADDED                   APR 15 ARDEAN LEITH
C               
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2015  Health Research Inc.,                         *
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

       IF (NGOTPAR < 0) THEN
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

C      COMPILE CCROT CHANGE STATISTICS

       CCROTAVG = CCROTAVG + CCROT

       IF (NGOTPAR >= 8) THEN
           IF (ANGDIF > ANGDIFTHR) IBIGANGDIF = IBIGANGDIF + 1
           ANGDIFAVG = ANGDIFAVG + ANGDIF
           CCDIF     = CCROT - CCROTLAS     ! DIFFERENCE

           IF (CCDIF >= 0.0) THEN
              IMPROVCCROT = IMPROVCCROT + 1
              CCROTIMPROV = CCROTIMPROV + CCDIF
           ELSE
              IWORSECCROT = IWORSECCROT + 1
              CCROTWORSE  = CCROTWORSE + ABS(CCDIF)
           ENDIF
       ENDIF   ! END OF: IF (NGOTPAR >= 8)

       END

C       -------------------- AP_STAT_ADD_R -------------------------------

      SUBROUTINE AP_STAT_ADD_R(NGOTPAR,
     &              ANGDIF,ANGDIFTHR, IBIGANGDIF,ANGDIFAVG,
     &              CC,CCLAS,CCAVG,
     &              IMPROVCC,CCIMPROV,IWORSECC,CCWORSE,
     &              FIRSTCC)

       IMPLICIT  NONE

       INTEGER  :: NGOTPAR,IBIGANGDIF, IMPROVCC,IWORSECC

       REAL     :: ANGDIF,ANGDIFTHR,ANGDIFAVG,
     &             CC,CCLAS,CCAVG,CCIMPROV,CCWORSE 

       LOGICAL  :: FIRSTCC
       REAL     :: CCDIF

       IF (NGOTPAR < 0) THEN
C          INITIALIZE CCR CHANGE STATISTICS

           IBIGANGDIF = 0
           ANGDIFAVG  = 0.0

           CCAVG      = 0.0
           IMPROVCC   = 0
           CCIMPROV   = 0.0
           IWORSECC   = 0
           CCWORSE    = 0.0

           RETURN
        ENDIF  

C       COMPILE CC CHANGE STATISTICS

        CCAVG = CCAVG + CC

        IF (NGOTPAR >= 8) THEN
           IF (ANGDIF > ANGDIFTHR) IBIGANGDIF = IBIGANGDIF + 1

           ANGDIFAVG = ANGDIFAVG + ANGDIF

           IF (.NOT. FIRSTCC) THEN
              CCDIF     = CC - CCLAS              ! DIFFERENCE

              IF (CCDIF >= 0.0) THEN
                 IMPROVCC = IMPROVCC + 1
                 CCIMPROV = CCIMPROV  + CCDIF
              ELSE
                 IWORSECC = IWORSECC + 1
                 CCWORSE  = CCWORSE + ABS(CCDIF)
              ENDIF
           ENDIF  ! END OF: IF (.NOT. FIRSTCC)
        ENDIF     ! END OF: IF (NGOTPAR >= 8)

        END


C       -------------------- AP_STAT -------------------------------

       SUBROUTINE AP_STAT(NUMEXP,ANGDIFTHR,IBIGANGDIF, ANGDIFAVG, 
     &            CCRAVG, IMPROVCCR,CCRIMPROV, IWORSECCR, CCRWORSE,
     &            NBORDER,NSUBPIX,LUNDOC)

        IMPLICIT NONE 
          
        INTEGER         :: NUMEXP,IBIGANGDIF,IMPROVCCR,IWORSECCR 
        INTEGER         :: NBORDER,NSUBPIX,LUNDOC

        REAL            :: ANGDIFTHR,ANGDIFAVG,CCRAVG,CCRIMPROV,CCRWORSE
        REAL            :: FNUMEXP
      
        INTEGER         :: NLIST, NC, IRTFLG
        REAL            :: DLIST(9)

        CHARACTER(LEN=132) :: COMMENT 


C       SAVE CCR & ANGULAR DISPLACEMENT STATISTICS
        FNUMEXP  = NUMEXP

        DLIST(1) = 0.0
        IF (ANGDIFTHR > 0)DLIST(1) = 100 * FLOAT(IBIGANGDIF)/FNUMEXP

        DLIST(2) = ANGDIFAVG / FNUMEXP
        DLIST(3) = CCRAVG  / FNUMEXP
        DLIST(4) = 100 * FLOAT(IWORSECCR) / FNUMEXP

        DLIST(5) = 0.0
        IF (IWORSECCR > 0)DLIST(5) = CCRWORSE/FLOAT(IWORSECCR)

        DLIST(6) = 0.0
        IF (IMPROVCCR > 0)DLIST(6) = CCRIMPROV/FLOAT(IMPROVCCR)

        DLIST(7) = FNUMEXP

        NLIST    = 9
        NC       = 132

        IF (NBORDER > 0 .OR. NSUBPIX > 0) THEN

           DLIST(8) = NBORDER
           DLIST(9) = NSUBPIX

C           123456789 123456789 123456789 123456789 123456789 123456789 1
           COMMENT =  
     &     '        ' // 
     &     '%BIG-ANGDIF,   AVG-ANGDIF,    AVG-CCR,      %WORSE,    '//
     &     ' AVG-WORSE,   AVG-BETTER,    #PARTICLES,        #BORDER,'//
     &     '    #SUB_PXL'
      
        ELSE
           NLIST = 7
           
C           123456789 123456789 123456789 123456789 123456789 123456789 1
           COMMENT =      
     &     '        '// 
     &     '%BIG-ANGDIF,   AVG-ANGDIF,      AVG-CCR,      %WORSE,'//
     &     '     AVG-WORSE,   AVG-BETTER,    #PARTICLES'
        ENDIF

        CALL LUNDOCPUTCOM(LUNDOC,COMMENT(:NC),IRTFLG)
        CALL LUNDOCWRTDAT(LUNDOC,-2,DLIST,NLIST,IRTFLG)

        CLOSE(LUNDOC)

        END



C       -------------------- AP_STAT_SHC -------------------------------

       SUBROUTINE AP_STAT_SHC(NUMEXP,ANGDIFTHR,IBIGANGDIF, ANGDIFAVG, 
     &            CCRAVG, IMPROVCCR,CCRIMPROV, IWORSECCR, CCRWORSE,
     &            NBORDER,NSUBPIX,LUNDOC,FIRST)

        IMPLICIT NONE 
          
        INTEGER         :: NUMEXP,IBIGANGDIF,IMPROVCCR,IWORSECCR 
        INTEGER         :: NBORDER,NSUBPIX,LUNDOC
        LOGICAL         :: FIRST

        REAL            :: ANGDIFTHR,ANGDIFAVG,CCRAVG,CCRIMPROV,CCRWORSE
        REAL            :: FNUMEXP
      
        INTEGER         :: NLIST, NC, IRTFLG
        REAL            :: DLIST(9)

        CHARACTER(LEN=132) :: COMMENT 


C       SAVE CCR & ANGULAR DISPLACEMENT STATISTICS
        FNUMEXP  = NUMEXP

        DLIST    = 0.0

        IF (ANGDIFTHR > 0) DLIST(1) = 100 * FLOAT(IBIGANGDIF) / FNUMEXP

        DLIST(2) = ANGDIFAVG / FNUMEXP
        DLIST(3) = CCRAVG    / FNUMEXP
        DLIST(4) = 100 * FLOAT(IWORSECCR) / FNUMEXP

        IF (IWORSECCR > 0) DLIST(5) = CCRWORSE / FLOAT(IWORSECCR)
        IF (IMPROVCCR > 0) DLIST(6) = CCRIMPROV / FLOAT(IMPROVCCR)

        DLIST(7) = FNUMEXP
        DLIST(8) = NBORDER
        DLIST(9) = NSUBPIX

        NLIST    = 9
        NC       = 132

        IF (FIRST) THEN

           DLIST(1) = 0.0
           DLIST(2) = 0.0
           DLIST(4) = 0.0
           DLIST(5) = 0.0
           DLIST(6) = 0.0

C           123456789 123456789 123456789 123456789 123456789 123456789 1
           COMMENT =  
     &     '        ' // 
     &     '     UNUSED,       UNUSED,     AVG-CCR,       UNUSED,    '//
     &     '    UNUSED,      UNUSED,     #PARTICLES,     #BORDER,'    //
     &     '     #SUB_PXL'
      
        ELSE
           
C           123456789 123456789 123456789 123456789 123456789 123456789 1
           COMMENT =      
     &     '        '// 
     &     '%BIG-ANGDIF,   AVG-ANGDIF,     AVG-CCR,       %WORSE,    '//
     &     ' AVG-WORSE,   AVG-BETTER,    #PARTICLES,     #BORDER, '   //
     &     '    #SUB_PXL'
        ENDIF

        CALL LUNDOCPUTCOM(LUNDOC,COMMENT(:NC),IRTFLG)
        CALL LUNDOCWRTDAT(LUNDOC,-2,DLIST,NLIST,IRTFLG)

        CLOSE(LUNDOC)

        END





C       -------------------- AP_STAT_R -------------------------------

        SUBROUTINE AP_STAT_R(NUMEXP,ANGDIFTHR,IBIGANGDIF,ANGDIFAVG, 
     &                      CCAVG,IMPROVCC,CCIMPROV, IWORSECC,CCWORSE,
     &                      FIRST,FIRSTCC, LUNDOC)
 
        IMPLICIT NONE

        INTEGER             :: NUMEXP,IBIGANGDIF,IMPROVCC,IWORSECC
        LOGICAL             :: FIRST,FIRSTCC
        INTEGER             :: LUNDOC
      
        REAL                :: ANGDIFTHR, ANGDIFAVG, CCAVG
        REAL                :: CCIMPROV, CCWORSE
 
        INTEGER             :: NLIST,NC,IRTFLG
        REAL                :: FNUMEXP
        REAL                :: DLIST(7)
        CHARACTER (LEN=132) :: COMMENT 

C       SAVE CC & ANGULAR DISPLACEMENT STATISTICS IN DOC FILE

        FNUMEXP  = NUMEXP
        DLIST    = 0.0
        DLIST(3) = CCAVG     / FNUMEXP
        DLIST(7) = FNUMEXP

        NLIST    = 7
        NC       = 120

        IF (FIRST) THEN

C               123456789 123456789 123456789 123456789 123456789 123456789 12
           COMMENT =      
     &         '        ' // 
     &         '     UNUSED,        UNUSED,         AVG-CC,    '  //
     &         'UNUSED,         UNUSED,       UNUSED,   #PARTICLES'

        ELSE

           IF (ANGDIFTHR > 0)DLIST(1) = 100 * FLOAT(IBIGANGDIF)/FNUMEXP

           DLIST(2) = ANGDIFAVG / FNUMEXP
           DLIST(3) = CCAVG     / FNUMEXP
           DLIST(4) = 100 * FLOAT(IWORSECC) / FNUMEXP

           IF (IWORSECC > 0) DLIST(5) = CCWORSE / FLOAT(IWORSECC)
           IF (IMPROVCC > 0) DLIST(6) = CCIMPROV / FLOAT(IMPROVCC)

           IF ( FIRSTCC) THEN

C               123456789 123456789 123456789 123456789 123456789 123456789 1
              COMMENT =      
     &         '        ' // 
     &         '%BIG-ANGDIF,    AVG-ANGDIF,         AVG-CC,    '  //
     &         'UNUSED,      UNUSED,       UNUSED,     #PARTICLES'
           ELSE
C               123456789 123456789 123456789 123456789 123456789 123456789 1
              COMMENT =      
     &         '        ' // 
     &         '%BIG-ANGDIF,    AVG-ANGDIF,         AVG-CC,    '  //
     &         '%WORSE,      AVG-WORSE,   AVG-BETTER,   #PARTICLES'
           ENDIF

        ENDIF

C       WRITE TO DOC FILE
        CALL LUNDOCPUTCOM(LUNDOC,COMMENT(:NC),IRTFLG)
        CALL LUNDOCWRTDAT(LUNDOC,-2,DLIST,NLIST,IRTFLG)

        CLOSE(LUNDOC)

        END



