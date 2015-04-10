

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

        SUBROUTINE RMFILTY(VOLIN,SLICE,TEMP,IDTEMP,NSAM,NROW,NSLICE,RAD,
     $                     PFRO,PTO,PINC,FND,FNM)

C       PROGRAM TO FILTER Z-SLICES IN A 3D RADON TRANSFORM

        COMMON /DEBUG1/DEBFLAG
       DIMENSION VOLIN(NSAM,NROW,NSLICE),SLICE(NSAM*NSLICE),TEMP(IDTEMP)


        LOGICAL DEBFLAG
        CHARACTER*(*) FND,FNM

        NSAM2=NSAM-4
        NSAMF=NSAM-2
        ISLDIM=NSAMF*NROW

        DO I=1,NROW
           ANG=(I-1)*PINC
           C=COS(ANG)

           CALL CUTYSLICE(VOLIN,SLICE,NSAM,NROW,NSLICE,NSAMF,I)

           NSL=1
           IF(DEBFLAG) CALL WRITEPICT(SLICE,NSAMF,NSLICE,NSL,51,1)
           CALL RMMAKEMASK(TEMP,NSAMF,NSLICE,C,RAD)

           IF(DEBFLAG) CALL WRITEPICT(TEMP,NSAMF,NROW,NSL,51,2)
           IF(I.EQ.1) CALL WRITEPICTN(TEMP,NSAMF,NROW,NSL,51,FNM)

           CALL RMFILTERZ(SLICE,TEMP,NSAM2,NSLICE,ISLDIM,FND,I)
           IF (DEBFLAG) CALL WRITEPICT(SLICE,NSAMF,NROW,NSL,51,3)

           CALL PUTYSLICE(VOLIN,SLICE,NSAM,NROW,NSLICE,NSAMF,I)
        ENDDO

        RETURN
        END
