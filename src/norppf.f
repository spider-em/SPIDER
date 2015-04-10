
C ++********************************************************************
C                                                                      *
C                                                                      *
C                                                                      *
C **********************************************************************
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C 
C SUPPORT_ROUTINE                                                                      *
C
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE NORPPF(P,PPF,IER)

C NORPPF
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
C     PURPOSE--THIS SUBROUTINE COMPUTES THE PERCENT POINT
C              FUNCTION VALUE FOR THE NORMAL (GAUSSIAN)
C              DISTRIBUTION WITH MEAN = 0 AND STANDARD DEVIATION = 1. 
C              THIS DISTRIBUTION IS DEFINED FOR ALL X AND HAS
C              THE PROBABILITY DENSITY FUNCTION
C              F(X) = (1/SQRT(2*PI))*EXP(-X*X/2). 
C              NOTE THAT THE PERCENT POINT FUNCTION OF A DISTRIBUTION 
C              IS IDENTICALLY THE SAME AS THE INVERSE CUMULATIVE
C              DISTRIBUTION FUNCTION OF THE DISTRIBUTION.
C     INPUT  ARGUMENTS--P      = THE SINGLE PRECISION VALUE 
C                                (BETWEEN 0.0 AND 1.0)
C                                AT WHICH THE PERCENT POINT 
C                                FUNCTION IS TO BE EVALUATED.
C     OUTPUT ARGUMENTS--PPF    = THE SINGLE PRECISION PERCENT
C                                POINT FUNCTION VALUE.
C     OUTPUT--THE SINGLE PRECISION PERCENT POINT
C             FUNCTION VALUE PPF.
C     PRINTING--NONE UNLESS AN INPUT ARGUMENT ERROR CONDITION EXISTS. 
C     RESTRICTIONS--P SHOULD BE BETWEEN 0.0 AND 1.0, EXCLUSIVELY.
C     OTHER DATAPAC   SUBROUTINES NEEDED--NONE.
C     FORTRAN LIBRARY SUBROUTINES NEEDED--SQRT, ALOG.
C     MODE OF INTERNAL OPERATIONS--SINGLE PRECISION.
C     LANGUAGE--ANSI FORTRAN. 
C     REFERENCES--ODEH AND EVANS, THE PERCENTAGE POINTS
C                 OF THE NORMAL DISTRIBUTION, ALGORTIHM 70, 
C                 APPLIED STATISTICS, 1974, PAGES 96-97.
C               --EVANS, ALGORITHMS FOR MINIMAL DEGREE
C                 POLYNOMIAL AND RATIONAL APPROXIMATION,
C                 M. SC. THESIS, 1972, UNIVERSITY 
C                 OF VICTORIA, B. C., CANADA.
C               --HASTINGS, APPROXIMATIONS FOR DIGITAL
C                 COMPUTERS, 1955, PAGES 113, 191, 192.
C               --NATIONAL BUREAU OF STANDARDS APPLIED MATHEMATICS
C                 SERIES 55, 1964, PAGE 933, FORMULA 26.2.23.
C               --FILLIBEN, SIMPLE AND ROBUST LINEAR ESTIMATION
C                 OF THE LOCATION PARAMETER OF A SYMMETRIC
C                 DISTRIBUTION (UNPUBLISHED PH.D. DISSERTATION,
C                 PRINCETON UNIVERSITY), 1969, PAGES 21-44, 229-231.
C               --FILLIBEN, 'THE PERCENT POINT FUNCTION',
C                 (UNPUBLISHED MANUSCRIPT), 1970, PAGES 28-31.
C               --JOHNSON AND KOTZ, CONTINUOUS UNIVARIATE
C                 DISTRIBUTIONS--1, 1970, PAGES 40-111.
C               --THE KELLEY STATISTICAL TABLES, 1948.
C               --OWEN, HANDBOOK OF STATISTICAL TABLES,
C                 1962, PAGES 3-16.
C               --PEARSON AND HARTLEY, BIOMETRIKA TABLES
C                 FOR STATISTICIANS, VOLUME 1, 1954,
C                 PAGES 104-113.
C     COMMENTS--THE CODING AS PRESENTED BELOW
C               IS ESSENTIALLY IDENTICAL TO THAT
C               PRESENTED BY ODEH AND EVANS
C               AS ALGORTIHM 70 OF APPLIED STATISTICS.
C               THE PRESENT AUTHOR HAS MODIFIED THE
C               ORIGINAL ODEH AND EVANS CODE WITH ONLY
C               MINOR STYLISTIC CHANGES.
C             --AS POINTED OUT BY ODEH AND EVANS
C               IN APPLIED STATISTICS,
C               THEIR ALGORITHM REPRESENTES A
C               SUBSTANTIAL IMPROVEMENT OVER THE
C               PREVIOUSLY EMPLOYED
C               HASTINGS APPROXIMATION FOR THE
C               NORMAL PERCENT POINT FUNCTION--
C               THE ACCURACY OF APPROXIMATION
C               BEING IMPROVED FROM 4.5*(10**-4)
C               TO 1.5*(10**-8).
C     WRITTEN BY--JAMES J. FILLIBEN
C                 STATISTICAL ENGINEERING LABORATORY (205.03)
C                 NATIONAL BUREAU OF STANDARDS
C                 WASHINGTON, D. C. 20234
C                 PHONE:  301-921-2315
C     ORIGINAL VERSION--JUNE      1972. 
C     UPDATED         --SEPTEMBER 1975. 
C     UPDATED         --NOVEMBER  1975. 
C     UPDATED         --OCTOBER   1976. 
C MODIFIED BY PP
C
C---------------------------------------------------------------------
C     CHECK THE INPUT ARGUMENTS FOR ERRORS
C
      IF(P.LE.0.0.OR.P.GE.1.0)  THEN
        IER=1
C    1 FORMAT(1H ,115H***** FATAL ERROR--THE FIRST  INPUT ARGUMENT TO THE
C     1 NORPPF SUBROUTINE IS OUTSIDE THE ALLOWABLE (0,1) INTERVAL *****)
      ELSE
        IER=0
C
C-----START POINT-----------------------------------------------------
C
      IF(P.EQ.0.5)  THEN
      PPF=0.0
      ELSE
      R=P 
      IF(P.GT.0.5)R=1.0-R
      T=SQRT(-2.0*ALOG(R))
      ANUM=((((T*(-.453642210148E-4)-.204231210245E-1)*T
     &  -.342242088547)*T-1.0)*T-.322232431088)
      ADEN=((((T*.38560700634E-2+.103537752850)*T
     &  +.531103462366)*T+.588581570495)*T+.993484626060E-1)
      PPF=T+(ANUM/ADEN)
      IF(P.LT.0.5)PPF=-PPF
      ENDIF
      ENDIF
      END 


