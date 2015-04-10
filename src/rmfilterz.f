
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


        SUBROUTINE RMFILTERZ(SLICE,AMASK,NSAM,NROW,ISLDIM,FND,III)
       
C       PROGRAM TO APPLY MASK TO THE FOURIER TRANSFORM OF THE 2D RADON TRANSFORM

        DIMENSION SLICE(ISLDIM),AMASK(ISLDIM)
        COMMON /DEBUG1/DEBFLAG
        CHARACTER*(*) FND
        LOGICAL DEBFLAG
       
       
        IDIR=1
        CALL FMRS_2(SLICE,NSAM,NROW,IDIR)

        NSAMF=NSAM+2
        IF(DEBFLAG)  CALL WRITEPICT(SLICE,NSAMF,NROW,1,51,4)
        IF(III.EQ.1) CALL WRITEPICTN(SLICE,NSAMF,NROW,1,51,FND)

C       APPLY MASK:
        DO I=1,ISLDIM
           SLICE(I)=SLICE(I)*AMASK(I)
        ENDDO

        IF(DEBFLAG) CALL WRITEPICT(SLICE,NSAMF,NROW,1,51,5)

        IDIR=-1
        CALL FMRS_2(SLICE,NSAM,NROW,IDIR)

        IF(DEBFLAG) CALL WRITEPICT(SLICE,NSAMF,NROW,1,51,6)

        RETURN
        END
