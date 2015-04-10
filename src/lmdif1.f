
C++*********************************************************************
C
C LMDIF1(FCN,M,N,X,FVEC,TOL,INFO, IWA,WA,LWA)
C
C CALL TREE:   'SP' -->  DIFF1O --> LATTICE --> SOLVE
C                               --> WFTCIRC --> LMDIF1  --> LMDIF
C
C PURPOSE: MINIMIZE THE SUM OF THE SQUARES OF M NONLINEAR
C          FUNCTIONS IN N VARIABLES BY A MODIFICATION OF THE
C          LEVENBERG-MARQUARDT ALGORITHM. THIS IS DONE BY USING THE MORE
C          GENERAL LEAST-SQUARES SOLVER LMDIF. USER MUST PROVIDE A
C          SUBROUTINE WHICH CALCULATES THE FUNCTIONS. THE JACOBIAN IS
C          THEN CALCULATED BY A FORWARD-DIFFERENCE APPROXIMATION.
c
c     The subroutine statement is:
c
c       subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
c
c     Where:
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program, and should be written as follows.
c
c         subroutine       fcn(m,n,x,fvec,iflag)
c         integer          m,n,iflag
c         double precision x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of lmdif1.
c         in this case set iflag to a negative integer.
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an array of length n. on input x must contain
c         an initial estimate of the solution vector. on output x
c         contains the final estimate of the solution vector.
c
c       fvec is an output array of length m which contains
c         the functions evaluated at the output x.
c
c       tol is a nonnegative input variable. termination occurs
c         when the algorithm estimates either that the relative
c         error in the sum of squares is at most tol or that
c         the relative error between x and the solution is at
c         most tol.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  algorithm estimates that the relative error
c                   in the sum of squares is at most tol.
c
c         info = 2  algorithm estimates that the relative error
c                   between x and the solution is at most tol.
c
c         info = 3  conditions for info = 1 and info = 2 both hold.
c
c         info = 4  fvec is orthogonal to the columns of the
c                   jacobian to machine precision.
c
c         info = 5  number of calls to fcn has reached or
c                   exceeded 200*(n+1).
c
c         info = 6  tol is too small. no further reduction in
c                   the sum of squares is possible.
c
c         info = 7  tol is too small. no further improvement in
c                   the approximate solution x is possible.
c
c       iwa is an integer work array of length n.
c
c       wa  is a work array of length lwa.
c
c       lwa is a positive integer input variable not less than
c         m*n+5*n+m.
c
c     Subprograms called:
c
c       User-supplied ...... fcn
c
c       Minpack-supplied ... lmdif
c
c     Argonne National Laboratory. Minpack Project. March 1980.
c     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
c
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE LMDIF1(FCN,M,N,X,FVEC,TOL,INFO,IWA,WA,LWA)

      IMPLICIT NONE

      INTEGER          :: M,N,INFO,LWA
      INTEGER          :: IWA(N)
      DOUBLE PRECISION :: TOL
      DOUBLE PRECISION :: X(N),FVEC(M),WA(LWA)

      EXTERNAL         :: FCN

      INTEGER          :: MAXFEV,MODE,MP5N,NFEV,NPRINT
      DOUBLE PRECISION :: EPSFCN,FACTOR,FTOL,GTOL,XTOL,ZERO
     
      DATA  FACTOR,ZERO /1.0D2, 0.0D0/

      INFO = 0

C     CHECK THE INPUT PARAMETERS FOR ERRORS.
      IF (N   <= 0 .OR. M < N .OR. TOL < ZERO .OR.
     $    LWA < (M*N + 5*N + M))  RETURN

      MAXFEV = 200 * (N + 1)
      FTOL   = TOL
      XTOL   = TOL
      GTOL   = ZERO
      EPSFCN = ZERO
      MODE   = 1
      NPRINT = 0
      MP5N   = M + 5 * N

      CALL LMDIF(FCN,M,N,X,FVEC,FTOL,XTOL,GTOL,MAXFEV,EPSFCN,WA(1),
     *           MODE,FACTOR,NPRINT,INFO,NFEV,WA(MP5N+1),M,IWA,
     *           WA(N+1),WA(2*N+1),WA(3*N+1),WA(4*N+1),WA(5*N+1))

      IF (INFO == 8) INFO = 4

      END


      
