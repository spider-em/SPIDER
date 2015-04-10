
C++*********************************************************************
C
C REFORM.F                          CREATED        JAN 91
C                                   Y BUG          NOV 2003 ArDean Leith
C                                   ALLOCATED      APR 2010 ArDean Leith
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
C    REFORM(LUNIN,LUNOUT,NSAM,NSAM1,NSAM2,NSAMS,
C           NROW,NROW1,NROW2,NROWS,NSLICE,NSLICE1,NSLICE2,NSLICES,
C           AXIS,IANG,IRTFLG)
C
C    PURPOSE:   REFORM AN IMAGE STAGE BY ROTATION OF 90, 180 OR 270 
C               DEGREES AROUND  THE X,Y, OR Z AXIS.  CENTER OF ROTATION
C               IS ALWAYS CENTER OF THE RELEVANT SLICE
C               ROTATIONS ARE CLOCKWISE WHEN LOOKING ALONG AXES.  THE
C               X AXIS + TO RIGHT,  Y AXIS + UP THE SCREEN, Z AXIS + OUT
C              OF THE SCREEN!
C
C    PARAMETERS:       
C                     
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************


        SUBROUTINE REFORM(LUNIN,LUNOUT,NSAM,NSAM1,NSAM2,NSAMS,
     &           NROW,NROW1,NROW2,NROWS,NSLICE,NSLICE1,NSLICE2,NSLICES,
     &           AXIS,IANG,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON /IOBUF/ BUF(NBUFSIZ)

        CHARACTER (LEN=1) :: AXIS
        REAL, ALLOCATABLE :: BUF2(:)

        IF (AXIS .EQ. 'Z') THEN
C          ROTATE AROUND Z AXIS ----------------------------------- Z

           ALLOCATE(BUF2(NSAM*NROW), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
               CALL ERRT(46,'BUF2',NSAM*NROW)
               GOTO 9999
           ENDIF

           IF (IANG .EQ. 90 .OR. IANG .EQ. 270) THEN
C             ROTATE 90 OR -90 AROUND Z AXIS

              IT5 = NSAMS * (NROWS + NSAM1 - 1)

              DO ISLICE1 = NSLICE1, NSLICE2
                ISLICE2  = ISLICE1 - NSLICE1 + 1

                DO IROW1 = NROW1,NROW2
                  IREC1  = (ISLICE1 - 1) * NROW + IROW1
                  CALL REDLIN(LUNIN,BUF,NSAM,IREC1)

                  IF (IANG .EQ. 270) THEN
                     ISAM2 = NSAMS - (IROW1 - NROW1)
                     IT    = ISAM2 - (NSAM1 * NSAMS)
                  ELSE
                     ISAM2 = IROW1 - NROW1 + 1
                     IT3   = IT5 + ISAM2
                  ENDIF                  

                  DO ISAM1 = NSAM1,NSAM2
                     IF (IANG .EQ. 270) THEN
C                      IROW2  = ISAM1 - NSAM1 + 1
C                      IPT    = (IROW2 - 1) * NSAMS + ISAM2
C                      IPT    = (ISAM1 - NSAM1) * NSAMS + ISAM2
C                      IPT    = ISAM1 * NSAMS - NSAM1 * NSAMS + ISAM2
                       IPT    = ISAM1 * NSAMS + IT
                    ELSE
C                      IROW2  = NROWS - (ISAM1 - NSAM1)
C                      IPT    = (IROW2 - 1) * NSAMS + ISAM2
C                      IPT    = (NROWS - ISAM1 + NSAM1 - 1) * NSAMS + ISAM2
C                      IPT    = (NROWS + NSAM1 - ISAM1 - 1) * NSAMS + ISAM2
C                      IPT    = NROWS*NSAMS + NSAM1*NSAMS - ISAM1*NSAMS - NSAMS +
C                               ISAM2
C                      IT3    = NROWS*NSAMS + NSAM1*NSAMS - NSAMS + ISAM2
C                      IT5    = NSAMS * (NROWS + NSAM1 - 1)
C                      IT3    = IT5 + ISAM2
                       IPT    = IT3 - ISAM1 * NSAMS
                    ENDIF
                    BUF2(IPT) = BUF(ISAM1)
                  ENDDO
                ENDDO

C               OUTPUT THE REFORMED SLICE
                DO IROW2 = 1,NROWS
                  IPT    = (IROW2 - 1)   * NSAMS + 1
                  IREC2  = (ISLICE2 - 1) * NROWS + IROW2
                  CALL WRTLIN(LUNOUT,BUF2(IPT),NSAMS,IREC2)
                ENDDO
              ENDDO

           ELSEIF (IANG .EQ. 180) THEN
C             ROTATE 180 DEGREES AROUND Z AXIS ! -------------------- 

              IT7        = NSAMS + NSAM1
              DO ISLICE1 = NSLICE1,NSLICE2
                ISLICE2  = ISLICE1 - NSLICE1 + 1
                DO IROW1 = NROW1,NROW2
                  IREC1  = (ISLICE1 -1) * NROW + IROW1
                  IROW2  = NROWS - (IROW1 - NROW1)
                  IREC2  = (ISLICE2 -1) * NROWS + IROW2
                  CALL REDLIN(LUNIN,BUF,NSAM,IREC1)

C                 INVERT THE ROW
                  DO ISAM1 = NSAM1,NSAM2
C*                   ISAM2 = NSAMS - (ISAM1 - NSAM1)
C*                   ISAM2 = NSAMS - ISAM1 + NSAM1
C*                   IT7   = NSAMS + NSAM1
C*                   ISAM2 = IT7   - ISAM1
C*                   BUF2(ISAM2) = BUF(ISAM1)
                     BUF2(IT7 - ISAM1) = BUF(ISAM1)
                  ENDDO
                  CALL WRTLIN(LUNOUT,BUF2,NSAMS,IREC2)
                ENDDO
              ENDDO                    
           ENDIF 

        ELSEIF (AXIS .EQ. 'X') THEN  ! ---------------------------- X

           DO ISLICE1 = NSLICE1,NSLICE2
             DO IROW1 = NROW1,NROW2
               IREC1 = (ISLICE1-1) * NROW + IROW1
               CALL REDLIN(LUNIN,BUF,NSAM,IREC1)

               IF (IANG .EQ. 90) THEN
                 ISLICE2 = IROW1 - NROW1 + 1
                 IROW2   = NROWS - (ISLICE1 - NSLICE1)

               ELSEIF (IANG .EQ. 180) THEN
                 ISLICE2 = NSLICES - (ISLICE1 - NSLICE1)
                 IROW2   = NROWS   - (IROW1   - NROW1)

               ELSEIF (IANG .EQ. 270) THEN
                 ISLICE2 = NSLICES - (IROW1   - NROW1)
                 IROW2   = ISLICE1 - NSLICE1  + 1
               ENDIF

               IREC2 = (ISLICE2 - 1) * NROWS + IROW2
               CALL WRTLIN(LUNOUT,BUF(NSAM1),NSAMS,IREC2)

             ENDDO
           ENDDO

        ELSE IF (AXIS .EQ. 'Y') THEN ! ---------------------------- Y

           ALLOCATE(BUF2(NSAM*NSLICE), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
               CALL ERRT(46,'BUF2',NSAM)
               GOTO 9999
           ENDIF

           IF (IANG .EQ. 180) THEN  ! --------------------------- Y 180

            IT = NSAMS + NSAM1
            DO ISLICE1 = NSLICE1,NSLICE2
C             INVERT THE SLICES
              ISLICE2  = NSLICES - (ISLICE1 - NSLICE1)

              DO IROW1 = NROW1,NROW2
                IROW2  = IROW1 - NROW1 + 1
                IREC1  = (ISLICE1 - 1) * NROW + IROW1
                CALL REDLIN(LUNIN,BUF,NSAM,IREC1)

C               MIRROR THE COLUMS ALONG THE ROW
                DO ISAM1 = NSAM1,NSAM2
                   BUF2(IT - ISAM1) = BUF(ISAM1)
                ENDDO

                IREC2  = (ISLICE2 - 1) * NROWS + IROW2
                CALL WRTLIN(LUNOUT,BUF2(NSAM1),NSAMS,IREC2)
              ENDDO
            ENDDO

          ELSEIF (IANG .EQ. 90 .OR. IANG .EQ. 270) THEN ! 
          
C          ROTATE AROUND  Y AXIS, RELATIVELY SLOW!!  ------------ Y 90

           DO IROW1  = NROW1,NROW2
C            FIND CORRESPONDING CURRENT OUTPUT SLICE 
             IROW2   = IROW1 - NROW1 + 1

             DO ISLICE1 = NSLICE1,NSLICE2

C              READ ONE RECORD FROM SOURCE FILE
               IREC1    = (ISLICE1 - 1) * NROW + IROW1
               CALL REDLIN(LUNIN,BUF,NSAM,IREC1)

               IF (IANG .EQ. 270) THEN
                  ISAM2   = ISLICE1 - NSLICE1 + 1
               ELSE
                  ISAM2   = NSAMS   - (ISLICE1 - NSLICE1)
               ENDIF

               ISAMT    = ISAM2
               IP1      = - NSAMS + ISAMT
               IP2      = - NSAM1 + 1

               DO ISAM1 = NSAM1,NSAM2
C                CREATE NEW X-Z IMAGE PLANE
C                ISLICE2 = ISAM1 - NSAM1 + 1
C                ISLICE2 = ISAM1 + IP2
C                IROWT   = ISLICE2
                 IROWT   = ISAM1 + IP2

C                BUF2((IROWT-1) * NSAMS + ISAMT)      = BUF(ISAM1)
C                BUF2((IROWT * NSAMS - NSAMS + ISAMT) = BUF(ISAM1)
                 BUF2(IROWT * NSAMS + IP1)            = BUF(ISAM1)
               ENDDO
             ENDDO

C            WRITE OUT THE REFORMED X-Y PLANE
             IPT    = 1 

             DO ISLICE2 = 1,NSLICES
                IREC2 = (ISLICE2 - 1) * NROWS + IROW2
                CALL WRTLIN(LUNOUT,BUF2(IPT),NSAMS,IREC2)
                IPT   = IPT + NSAMS
             ENDDO
           ENDDO
          ENDIF
        ENDIF

9999    IF (ALLOCATED(BUF2)) DEALLOCATE (BUF2) 
               
	RETURN
	END

