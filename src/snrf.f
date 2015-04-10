C **********************************************************************
C *  SNRF.F
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2002, L. Joyeux & P. A. Penczek                       *
C=*                                                                    *
C=* University of Texas - Houston Medical School                       *
C=*                                                                    *
C=* Email:  pawel.a.penczek@uth.tmc.edu                                *
C=*                                                                    *
C=* This program is free software; you can redistribute it and/or      *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* This program is distributed in the hope that it will be useful,    *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
C=* General Public License for more details.                           *
C=*                                                                    *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program; if not, write to the                      *
C=* Free Software Foundation, Inc.,                                    *
C=* 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.      *
C=*                                                                    *
C=**********************************************************************
C
C **********************************************************************

      SUBROUTINE SNRF

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'F90ALLOC.INC'

      REAL      FP, FS
      INTEGER   NB,NOT_USED
       REAL   EPS, AA
      PARAMETER(EPS=0.882)
      PARAMETER(AA=10.624)      
      INTEGER       I, IRTFLG
      REAL         ORD, RAD, FREQ, BUT, SNRFACTOR
        INTEGER      NLIST
      PARAMETER  (NLIST=3)
      REAL         DLIST(NLIST)
      INTEGER      NDOC

      DATA        NDOC/88/
      
      CALL RDPRM2(FP,FS,NOT_USED,'FREQUENCIES PASS AND STOP')

      CALL RDPRMI(NB,I,NOT_USED, 'HOW MANY POINTS')

      CALL RDPRM2(FMIN,FMAX,NOT_USED, 
     &   'RANGE [0,1] OF FSC IS MAPPED TO [FSCMIN:FSCMAX]')

      CALL RDPRM(SNRFACTOR, NOT_USED, 
     &     'FACTOR APPLIED ON FSC/(1-FSC) ')
      
      ORD = 2*ALOG10(EPS/SQRT(AA**2-1))/ALOG10(FP/FS)
      RAD = FP/EPS**(2/ORD)
      
      DO I=1, NB
         FREQ     = (I-1.)/(2*(NB-1))
         BUT      = 1/SQRT(1+(FREQ/RAD)**ORD)
         DLIST(1) = I
         DLIST(3) = BUT
         BUT      = (FMAX-FMIN)*BUT+FMIN
         DLIST(2) = SNRFACTOR*BUT/(1-BUT)
         CALL SAVD(NDOC,DLIST,NLIST,IRTFLG)         
      END DO

      CALL SAVDC
      CLOSE(NDOC)
      
      RETURN
      END
   
