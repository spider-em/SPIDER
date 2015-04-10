

C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2000  M. Radermacher                                  *
C=*                                                                    *
C=* Email:  spider@wadsworth.org                                       *
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


        SUBROUTINE GETSYM(IKIND,ISYM,IER)
       
        INCLUDE 'CMBLOCK.INC'
        
        IQ = 0

1       WRITE(NOUT,100)
100     FORMAT(' CHOICE OF SYMMETRIES:',/,
     &         ' (1) SYMMETRY AROUND Z-AXIS')

        CALL RDPRMI(IKIND,ISYM,NOT_USED,'ENTER SYMMETRY, ORDER OF SYM.')

        IF (IKIND .EQ. 0) ISYM=1
        IF (IKIND .GT. 1) THEN
           WRITE(NOUT,*)  '***** UNKNOWN SYMMETRY, REENTER'
           IQ = IQ+1
           IF (IQ. GT. 3) THEN
              WRITE(NOUT,101) 
101           FORMAT(' *** WELL, THREE TIMES IS ENOUGH ****',/,
     &        ' YOU HAVE THE CHOICE OF:'/1X,
     &        ' 1) TALK TO THE PROGRAMMER WHO WROTE THIS',/,
     &        '    AND TELL HIM THERE IS A MISTAKE.',/,
     &        ' 2) START FROM THE START',/,
     &        ' 3) READ THE MANUAL AND FIND OUT WHAT YOU DID WRONG',/,
     &        ' 4) THINK OF A NEW CAREER',/,
     &        ' 5) GO HOME AND GET SOME SLEEP')
              RETURN
           ENDIF  
           GOTO 1
        ENDIF

        RETURN
        END
