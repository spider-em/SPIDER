
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


        SUBROUTINE PUTYSLICE(VOL,YSLICE,NSAM,NROW,NSLICE,NSAMF,IPOS)
      
        DIMENSION VOL(NSAM,NROW,NSLICE),YSLICE(NSAMF*NSLICE)
C       PROGRAM TO PUT Y-SLICE INTO A 3D RADON OR RADON FOURIER FILE.
C       THE LINE CAN BE TRUNCATED AT NSAMF.
        DO I=1,NSLICE
           INDS=(I-1)*NSAMF
           DO K=1,NSAM-4
              VOL(K,IPOS,I) = YSLICE(INDS+K)
           ENDDO
        ENDDO
        RETURN
        END
