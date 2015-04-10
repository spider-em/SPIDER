
C ++********************************************************************
C                                                                      *
C                                                                      *
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
C                                                                      *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C NORCDF
C  Package : DATAPAC
C A Fortran subroutine library for probability distribution, density, percent
C point, and sparsity function evaluation; random number generation;
C line-printer plotting - histograms, scatter diagrams, probability plots; data
C manipulation; general statistical analysis; time series analysis; polynomial
C regression; ANOVA. (Approximately 170 subroutines.)
C 
C Type         : Subprogram library
C 
C Portablity   : Portable with some effort
C 
C Availability : Public domain
C
C     PURPOSE--THIS SUBROUTINE COMPUTES THE CUMULATIVE DISTRIBUTION
C              FUNCTION VALUE FOR THE NORMAL (GAUSSIAN)
C              DISTRIBUTION WITH MEAN = 0 AND STANDARD DEVIATION = 1. 
C              THIS DISTRIBUTION IS DEFINED FOR ALL X AND HAS
C              THE PROBABILITY DENSITY FUNCTION
C              F(X) = (1/SQRT(2*PI))*EXP(-X*X/2). 
C MODIFIED by PP 11/7/96
C   rather  F(X) = (2/SQRT(PI))*EXP(-X*X)
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VALUE AT
C                                WHICH THE CUMULATIVE DISTRIBUTION
C                                FUNCTION IS TO BE EVALUATED.
C     OUTPUT ARGUMENTS--CDF    = THE SINGLE PRECISION CUMULATIVE
C                                DISTRIBUTION FUNCTION VALUE.
C     OUTPUT--THE SINGLE PRECISION CUMULATIVE DISTRIBUTION
C             FUNCTION VALUE CDF.
C     PRINTING--NONE.
C     RESTRICTIONS--NONE.
C     OTHER DATAPAC   SUBROUTINES NEEDED--NONE.
C     FORTRAN LIBRARY SUBROUTINES NEEDED--EXP.
C     MODE OF INTERNAL OPERATIONS--SINGLE PRECISION.
C     LANGUAGE--ANSI FORTRAN. 
C     REFERENCES--NATIONAL BUREAU OF STANDARDS APPLIED MATHEMATICS
C                 SERIES 55, 1964, PAGE 932, FORMULA 26.2.17.
C               --JOHNSON AND KOTZ, CONTINUOUS UNIVARIATE
C                 DISTRIBUTIONS--1, 1970, PAGES 40-111.
C     WRITTEN BY--JAMES J. FILLIBEN
C                 STATISTICAL ENGINEERING LABORATORY (205.03)
C                 NATIONAL BUREAU OF STANDARDS
C                 WASHINGTON, D. C. 20234
C                 PHONE:  301-921-2315
C     ORIGINAL VERSION--JUNE      1972. 
C     UPDATED         --SEPTEMBER 1975. 
C     UPDATED         --NOVEMBER  1975. 
C
C SUPPORT_ROUTINE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      FUNCTION ORCDF(X)

C     Z=X

      Z=X*1.4142136 
      IF(X.LT.0.0)Z=-Z
      T=1.0/(1.0+.2316419*Z)
C     CDF=1.0-((0.39894228040143  )*EXP(-0.5*Z*Z))*(.319381530*T+
      CDF=0.79788456*EXP(-0.5*Z*Z)*(.319381530*T+
     &          (-0.356563782)*T**2+1.781477937*T**3
     &          +(-1.821255978)*T**4+1.330274429*T**5)
      IF(X.LT.0.0)  CDF=1.0-CDF
      ORCDF=CDF

      END 
