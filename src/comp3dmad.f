
C ++********************************************************************
C  COMP3DMAD                                                           *
C                                                                      *
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
C                                                                      *
C  PURPOSE: CALCULATES MEAN RELATIVE ERROR AND  DISCREPANCY            *
C           BETWEEN TWO VOLUMES WITHIN AN (OPTIONAL) MASK.             *
C                                                                      *
C  AUTHOR: M.RADERMACHER                                               *
C                                                                      *
C***********************************************************************

      SUBROUTINE COMP3DMAD(LUN,LUN1,LUN2,LUN3) 
                         
      IMPLICIT DOUBLE PRECISION (A-H,O-X)
                               
      INCLUDE 'CMBLOCK.INC' 
      INCLUDE 'CMLIMIT.INC' 

      CHARACTER (LEN=MAXNAM)          :: FIL2,FIL1,FIL3,FIL4
      REAL, ALLOCATABLE, DIMENSION(:) :: AIMG
                                                          
      CHARACTER *1  :: NULL,YN                                             
      INTEGER       :: NVOX                                                
      LOGICAL       :: FMASK,FOUT,FDIFMASK                                       
      REAL          :: RA,RCON,RDIF,RDIS 
                                                          
      NULL = CHAR(0)
                                                      
      WRITE(NOUT,*) ' COMPARISON OF TWO 3D ARRAYS'
 
      CALL OPFILEC(0,.TRUE.,FIL1,LUN,'O',IFORM1,NSAM,NROW,NSLICE,
     &           MAXIM,'FIRST',.FALSE.,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9100

      MAXIM = 0                         
      CALL OPFILEC(0,.TRUE.,FIL2,LUN1,'O',IFORM2,NSAM2,NROW2,NSLIC2,
     &           MAXIM,'SECOND',.FALSE.,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9100

      IF(NSAM2.NE.NSAM.OR.NROW2.NE.NROW.OR.NSLIC2.NE.NSLICE) GOTO 9100
      FMASK    = .TRUE.
      FDIFMASK = .TRUE.
      CALL FILERD(FIL3,NLET2,NULL,'MASK',IRTFLG)
      IF (IRTFLG .GT. 0) GOTO 9999
      IF (FIL3(1:1) .EQ. '*') FMASK = .FALSE.

      IF (FMASK) THEN 
         MAXIM = 0 
         CALL OPFILEC(0,.FALSE.,FIL3,LUN2,'O',IFORMM,NSAM3,NROW3,NSLIC3,
     &           MAXIM,'XXXX',.FALSE.,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9100

         IF (NSAM3.NE.NSAM.OR.NROW3.NE.NROW.OR.NSLIC3.NE.NSLICE)  
     &       GOTO 9100
  
         CALL RDPRMC(YN,NA,.TRUE., 
     &      'APPLY MASK TO BOTH, SCALING AND DIFFERENCE/ERROR? (Y/N)', 
     &      NULL,IRT) 
         IF (YN .EQ. 'N') FDIFMASK=.FALSE. 
      ENDIF  
                                                           
      FOUT = .TRUE. 
      CALL FILERD(FIL4,NLET2,NULL,'OUTPUT DIFFERENCE',IRTFLG) 
      IF (IRTFLG .GT. 0) GOTO 9999                                       
      IF (FIL4(1:1) .EQ. '*') FOUT = .FALSE.                                  
      IF (FOUT) THEN 
         MAXIM = 0 
         CALL OPFILEC(0,.FALSE.,FIL4,LUN3,'U',IFORM1,NSAM,NROW,NSLICE,
     &                MAXIM,'XXXX',.TRUE.,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9100
      ENDIF
      S1=0.
      S2=0.
      S3=0.
      S4=0.
      SUM1=0.
      SUM2=0.
      NVOX=0
      ALLOCATE (AIMG(4*NSAM), STAT=IRTFLG)
      IF (IRTFLG.NE.0) THEN
         CALL ERRT(46,'DR DIFF, AIMG',IER)
         GOTO 9999
      ENDIF

      DO ISL=1,NSLICE
         IOFF=(ISL-1)*NROW
         DO  I = 1,NROW
            CALL REDLIN(LUN,AIMG,NSAM,I+IOFF)
            CALL REDLIN(LUN1,AIMG(NSAM+1),NSAM,I+IOFF)
            IF(FMASK) CALL REDLIN(LUN2,AIMG(2*NSAM+1),NSAM,I+IOFF)
            DO  K = 1,NSAM
               IF (FMASK) THEN
                  IF (AIMG(2*NSAM+K).LT..5) CYCLE
               ENDIF                                                             
               NVOX = NVOX+1
               B    = AIMG(K)
               C    = AIMG(NSAM+K)
               P1   = B*C
               P2   = C*C
               S1   = S1+P1
               S2   = S2+P2
               S3   = S3+C
               S4   = S4+B
C              IF ((B.NE.0.OR.C.NE.0) .AND. NSAM.LE.16)
C     $        WRITE(NOUT,103)ISL,I,K,B,C,P1,P2,S1,S2,S3,S4
C103            FORMAT(1H ,'SLICE: ',I3,' , ROW: ',I3,' , COL: ',I3,/
C     $         ' VOXEL 1: ',F16.4,' , VOXEL 2: ',F16.4, /
C     $         '     AD1: ',F16.4,'       AD2: ',F16.4, /
C     $         '      S1: ',F16.4,'        S2: ',F16.4, /
C     $         '      S3: ',F16.4,'        S4: ',F16.4)
           ENDDO
         ENDDO
      ENDDO
      IF (NVOX .EQ. 0) THEN
        WRITE(NOUT,601)
601     FORMAT(//,' THERE ARE NO VOXELS PASSED WITH THE MASK PROVIDED'/
     &   ' A COMPARISON WOULD BE MEANINGLESS AND THEREFORE HAS NOT'/              
     &   ' BEEN DONE. IT WOULD BE ADVISABLE TO CHECK YOUR MASK.'//)
        DEALLOCATE (AIMG)
        GOTO 9999
      ENDIF

      A     = (S1-S3*S4/NVOX)/(S2-S3*S3/NVOX)
      CONST = (A*S3-S4)/NVOX
      BMEAN = S4/NVOX
      SUM1  = 0.
      SUM2  = 0.
      DNOM  = 0.
      DNUM  = 0.
      DO  ISL=1,NSLICE
         IOFF = (ISL-1)*NROW
         DO  I = 1,NROW
            CALL REDLIN(LUN,AIMG,NSAM,I+IOFF)
            CALL REDLIN(LUN1,AIMG(NSAM+1),NSAM,I+IOFF)
            IF (FMASK) CALL REDLIN(LUN2,AIMG(2*NSAM+1),NSAM,I+IOFF)

            DO  K = 1,NSAM
               IF (FOUT) AIMG(3*NSAM+K) = 0.
               IF (FMASK .AND. FDIFMASK) THEN
                  IF (AIMG(2*NSAM+K) .LT. .1) CYCLE
               ENDIF                                                         
               B     = AIMG(K)
               C     = AIMG(NSAM+K)
               BD    = A*C-B-CONST
               BD1   = DABS(BD)
               BD2   = DABS(B-BMEAN)
               SUM1  = SUM1+BD1
               SUM2  = SUM2+BD2
               BD1SQ = BD1*BD1
               DNOM  = DNOM+(B-BMEAN)*(B-BMEAN)
               DNUM  = DNUM+BD1SQ  
                                                 
               IF (FOUT) AIMG(3*NSAM+K) = BD
            ENDDO
            IF (FOUT) CALL WRTLIN(LUN3,AIMG(3*NSAM+1),NSAM,I+IOFF)
         ENDDO
      ENDDO

      DIF   = SUM1 / SUM2
      DISCR = DSQRT(DNUM / DNOM)

      WRITE(NOUT,*) ' '
      WRITE(NOUT,100) A, CONST
100   FORMAT('  SCALING OF SECOND VOLUME: VOXEL_NEW = VOXEL_OLD *',
     &          1pg14.6, ' - ',1pg14.6,/)
 
      WRITE(NOUT,101) DIF, DISCR
101   FORMAT('  MEAN RELATIVE ERROR: ',F10.5,/,                                   
     &       '  DISCREPANCY:    ',F15.5,/)

      RA = A
      RCON = CONST
      RDIF = DIF
      RDIS = DISCR

      CALL REG_SET_NSEL(1,4,RA,RCON, RDIF,RDIS,0.0,IRTFLG)
      DEALLOCATE (AIMG)
      GOTO 9999 
                                                        
9100  CALL ERRT(IER,'DR DIFF',NE)
                                        
9999  CLOSE (LUN)                                                       
      CLOSE (LUN1)                                                      
      CLOSE (LUN2)                                                      
      CLOSE (LUN3)                                                      
      END                                                               
       
C     IF ((B.NE.0.OR.C.NE.0).AND.NSAM.LE.16)
C     $ WRITE(NOUT,104)ISL,I,K,B,C,BD1,BD2,SUM1,SUM2                     
C104   FORMAT(' SLICE: ',I3,' , ROW: ',I3,' , COL: ',I3,/             
C     $' VOXEL 1: ',F16.4,' , VOXEL 2: ',F16.4, /                        
C     $'     BD1: ',F16.4,'       BD2: ',F16.4, /                        
C     $'    SUM1: ',F16.4,'      SUM2: ',F16.4)
