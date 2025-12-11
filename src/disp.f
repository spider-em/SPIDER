

C++*********************************************************************
C
C  DISP.F  -- EXTRACTED FROM COPYTOJPED           Jan 2020 ArDean Leith
C          -- IRTFLG 654321                       Oct 2025 ArDean Leith
C          -- REMOVE MRC EXTENSION                Oct 2025 ArDean Leith
C
C **********************************************************************
C=* AUTHOR: A. LEITH                                                   *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2025  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email:                                                             *
C=*                                                                    *
C=* SPIDER is free software; you can redistribute it and/or            *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* SPIDER is distributed in the hope that it will be useful,          *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* merchantability or fitness for a particular purpose.  See the      *
C=* GNU General Public License (www.gnu.org/licenses) for details.     *
C=*                                                                    *
C **********************************************************************
C
C  DISP()
C  DISP_RELION()
C
C  DISP:         COPIES SPIDER IMAGE TO JPEG USING IMAGEMAGICK 
C                THEN USES  SYSTEM COMMAND TO DISPLAY IMAGE. 
C
C  DISP_RELION:  INVOKES RELION_DIPLAY USING A SYSTEM COMMAND TO 
C                DISPLAY IMAGE(S). 
C
C  CALLED BY:    UTIL2
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE DISP()

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM) :: FILOLD,FILNEW
        CHARACTER(LEN=MAXNAM) :: OPTIONS
        CHARACTER(LEN=160)    :: COMLIN
 
        CHARACTER(LEN=1)      :: NULL = CHAR(0)
        LOGICAL               :: VERBOSET,WANTOUT,SAYIT
        INTEGER               :: NX,NY,NZ,MAXIM,ITYPE,IRTFLG,NLET,NLETC
        INTEGER               :: lnblnkn,NLETO,IGOEXT,NDUM

        INTEGER, PARAMETER    :: LUN1   = 14 
        INTEGER, PARAMETER    :: LUN2   = 15 
        INTEGER, PARAMETER    :: IDELAY = 3
          
        IF (FCHAR(4:4) == 'R') THEN
           CALL DISP_RELION()
           RETURN
        ENDIF


C       OPEN INPUT FILE, WHOLE STACK NOT ALLOWED
        MAXIM = 0
        ITYPE = 0   ! ?
        CALL OPFILEC(0,.TRUE.,FILOLD,LUN1,'O',ITYPE,
     &               NX,NY,NZ,MAXIM,'SPIDER OR MRC INPUT~6',
     &               .FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (NZ > 1) THEN
           CALL ERRT(101,'CAN NOT DISPLAY A VOLUME',NDUM)
           RETURN
        ENDIF

#if defined (SP_DBUGIO)
        nletc = lnblnkn(filold)
        write(3,*)' In disp; nletc,filold: ',nletc,filold(:nletc)
#endif

        IF (IMAMI .NE. 1) 
     &      CALL NORM3(LUN1,NX,NY,NZ,FMAX,FMIN,AV)

C       IF MRC AND ORIGIN IS 'LL' FLIP ORIGIN FOR DISPLAY
        CALL LUNFLIPORG_MRC(LUN1,NX,NY,NZ,.FALSE.,IRTFLG)

C       CREATE JPEG OUTPUT FILE NAME
        NLET   = lnblnkn(FILOLD)
        IGOEXT = SCAN(FILOLD,'.')
        IF(IGOEXT > 0 .AND. IGOEXT < NLET) NLET = IGOEXT - 1

#if defined (SP_DBUGIO)
        write(3,*)' In disp; nlet,filold(:nlet): ',nlet,filold(:nlet)
#endif

        FILNEW = FILOLD(:NLET) // '.jpg'
        NLET   = NLET + 4

       ! FILNEW = 'jnk_temp_for_display.jpg'
       ! NLET   = lnblnkn(FILNEW)

#if defined (SP_DBUGIO)
        write(3,*)' In disp; nlet,filnew(:nlet): ', nlet,filnew(:nlet)
#endif

C       CONVERT IMAGE FILE INTO JPG FILE
C       COPYTOJPG HAS A DELAY NEAR END SINCE AN INTERMEDIATE FILE .gray
C       IS CREATED TO HANDLE FLOATING POINTS AND NEGATIVE VALUES

        !!!!!!!!!!!!!VERBOSET = .FALSE.
        VERBOSET = .TRUE.
        CALL COPYTOJPG(LUN1,LUN2,FILNEW,NX,NY,NZ, VERBOSET,IDELAY)

        IRTFLG = 654321   ! KEEP LOWERCASE AND ACCEPT <CR>
        CALL RDPRMC(OPTIONS,NLETO,.TRUE.,
     &          'IMAGEMAGICK DISPLAY OPTIONS (or <CR>)',NULL,IRTFLG)

        WRITE(COMLIN,90) OPTIONS(1:NLETO),FILNEW(1:NLET)
90      FORMAT( ' display ', A,' ',A, ' &' )

#if defined (SP_DBUGIO)
        write(3,*)' In disp; comlin: ',comlin
        !sayit = .TRUE.
#endif

C       SAYIT: ECHO COMLIN TO OUTPUT
        SAYIT = .FALSE.
        CALL CSVMS(COMLIN,SAYIT,IRTFLG)

        END


C       ----------------- DISP_RELION ---------------------------------

        SUBROUTINE DISP_RELION()

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM) :: FILELINE
        CHARACTER(LEN=MAXNAM) :: OPTIONS
        CHARACTER(LEN=160)    :: COMLIN
 
        CHARACTER(LEN=1)      :: NULL = CHAR(0)
        INTEGER               :: IRTFLG
        INTEGER               :: NLET1,NLET2
        INTEGER               :: lnblnkn

        INTEGER, PARAMETER    :: IDELAY = 3
          
        IRTFLG = -999   ! KEEP LOWERCASE
        CALL RDPRMC(FILELINE,NLET1,.FALSE.,
     &              'SPIDER OR MRC IMAGE OR STACK',
     &              NULL,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IRTFLG = -999   ! KEEP LOWERCASE
        CALL RDPRMC(OPTIONS,NLET2,.FALSE.,
     &              'RELION DISPLAY OPTIONS (or <CR>)',
     &              NULL,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        WRITE(COMLIN,90) FILELINE(1:NLET1), OPTIONS(1:NLET2)
90      FORMAT( ' relion_display --i ', A,' ',A, ' &' )

C       DO NOT ECHO COMLIN
        CALL CSVMS(COMLIN,.FALSE.,IRTFLG)

        END





#ifdef NEVER
gyan 877% relion_display -h
+++ RELION: command line arguments (with defaults for optional ones between parantheses) +++
====== General options ===== 
                             --i () : Input STAR file, image or stack
                      --gui (false) : Use this to provide all other parameters through a GUI
           --display (rlnImageName) : Metadata label to display
                         --table () : Name of the table to read from in the input STAR file
                        --scale (1) : Relative scale
                        --black (0) : Pixel value for black (default is auto-contrast)
                        --white (0) : Pixel value for white (default is auto-contrast)
               --sigma_contrast (0) : Set white and black pixel values this many times the 
                                      image stddev from the mean
         --read_whole_stack (false) : Read entire stacks at once (to speed up when many 
                                      images of each stack are displayed)
  --show_fourier_amplitudes (false) : Show amplitudes of 2D Fourier transform?
  --show_fourier_phase_angles (false) : Show phase angles of 2D Fourier transforms?
====== Multiviewer options ===== 
                          --col (5) : Number of columns
             --apply_orient (false) : Apply the orientation as stored in the input STAR file angles 
                                       and offsets
                    --ori_scale (1) : Relative scale for viewing individual images in multiviewer
            --sort (EMDL_UNDEFINED) : Metadata label to sort images on
              --random_sort (false) : Use random order in the sorting
                  --reverse (false) : Use reverse order (from high to low) in the sorting
                    --class (false) : Use this to analyse classes in input model.star file
                     --regroup (-1) : Number of groups to regroup saved particles from 
                                      selected classes in (default is no regrouping)
               --allow_save (false) : Allow saving of selected particles or class averages
                       --fn_imgs () : Name of the STAR file in which to save selected images.
                      --fn_parts () : Name of the STAR file in which to save particles from 
                                      selected classes.
      --max_nr_parts_per_class (-1) : Select maximum this number of particles from each selected classes.
                 --recenter (false) : Recenter the selected images to the center-of-mass of 
                                      all positive pixel values. 
               --max_nr_images (-1) : Only show this many images (default is show all)
====== Picking options ===== 
                     --pick (false) : Pick coordinates in input image
           --pick_start_end (false) : Pick start-end coordinates in input image
                        --coords () : STAR file with picked particle coordinates
            --particle_radius (100) : Particle radius in pixels
                      --lowpass (0) : Lowpass filter (in A) to filter micrograph before displaying
                     --highpass (0) : Highpass filter (in A) to filter micrograph before displaying
                      --angpix (-1) : Pixel size (in A) to calculate lowpass filter
                    --color_star () : STAR file with a column for red-blue coloring (a subset of) 
                                      the particles
                   --color_label () : MetaDataLabel to color particles on (e.g. rlnParticleSelectZScore)
                        --blue (1.) : Value of the blue color
                         --red (0.) : Value of the red color
                         --verb (1) : Verbosity
                          --version : Print RELION version and exit
#endif

