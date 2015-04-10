C++*********************************************************************
C
C    MX.F                                          BIMAL RATH 3/14/2003         
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
C MX
C 
C  READS 2 REAL IMAGE FILES AND COMPARES CORRESPONDING PIXEL VALUES. PICKS
C  THE HIGHER PIXEL VALUE AND WRITES IT AT THE CORRESPONDING PIXEL POSITION
C  OF THE OUTPUT FILE 
C
C IMAGE_PROCESSING_ROUTINE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C-************************************************************************

        SUBROUTINE MX(LUN1,LUN2,LUN3,NSAM,NROW,NSLICE)

        INCLUDE 'CMBLOCK.INC'      
        INCLUDE 'CMLIMIT.INC' 
        
        REAL, DIMENSION(NSAM) :: BUF1
        REAL, DIMENSION(NSAM) :: BUF2
        REAL, DIMENSION(NSAM) :: BUF3             
                
        INTEGER  LUN1,LUN2,LUN3,NSAM,NROW,NSLICE     
            
        NREC  = NROW * NSLICE
        
        DO IREC=1,NREC
           CALL REDLIN(LUN1,BUF1,NSAM,IREC)
           CALL REDLIN(LUN2,BUF2,NSAM,IREC) 
                               
           DO  ISAM=1,NSAM
              PIXEL = BUF1(ISAM)
              IF (BUF1(ISAM) .LT. BUF2(ISAM)) PIXEL  = BUF2(ISAM)
              BUF3(ISAM) = PIXEL
           ENDDO
           
           CALL  WRTLIN(LUN3,BUF3,NSAM,IREC)

           
        ENDDO
         
        RETURN
        END
        
