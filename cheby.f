!
      SUBROUTINE CHEBY (M, ND, X, NP, Y)
!     Last change: 20 Jan 1990
!     Purpose    : to find the values Y(1),...Y(NP)
!                  of the NDthe derivative of the Cheby poly
!                  of degree M at the points
!                  X(1),...X(NP)
      INTEGER ONE,M,NP,IA,MM,J,ND,NJ 
      DOUBLE PRECISION A(500),AX(500),X(NP),Y(NP)
      IA = 500
      MM = M + 1
      ONE = 1
      IF (MM .GT. IA) GO TO 70
      IF (M .EQ. 0) GO TO 50
      DO 10 J = 1, M
        A(J) = 0D0
   10 CONTINUE 
      A(MM) = 1D0
      IF (ND .EQ. 0) GO TO 40
      DO 30 NJ = 1, ND
        CALL DERIV2(A, AX, MM)
        DO 20 J = 1, MM
          A(J) = AX(J)
   20   CONTINUE 
   30 CONTINUE 
   40 CALL EVAL2(A, MM, X, NP, Y)
      RETURN 
   50 DO 60 J = 1, NP
        IF (ND .EQ. 0) Y(J) = 1D0
        IF (ND .GT. 0) Y(J) = 0D0
   60 CONTINUE 
      RETURN 
   70 WRITE (6,80)
   80 FORMAT (' Sub CHEBY, M too big')
      RETURN 
      END 
! 
      SUBROUTINE EVAL2 (A, M, X, NX, Y)
!     Last change: 20 Jan 1990
!     Purpose    : to compute the value of the Cheby series
!                  A(1),...A(M) of length M, where A(1)....A(M)
!                  are respectively the coeffs of T(0)...T(M-1),
!                  at the points X(j) j=1,...NX.
!                  The output is Y(j) j=1,...NX, which is the
!                  value of the series at these NX points.
      DOUBLE PRECISION A(M),X(NX),Y(NX),B(0:502),XA,XB,XX
      INTEGER M,NX,IA,M1,M2,J,I
      IA = 502
      M1 = M + 1
      M2 = M + 2
      IF (M2 .GT. IA) GO TO 30
      DO 20 I = 1, NX
        XX = X(I)
        XA = 2D0 * XX - 1D0
        XB = 2D0 * XA
        B(M2) = 0D0
        B(M1) = 0D0
        DO 10 J = M - 1, 0, -1
          B(J + 1) = XB * B(J + 2) - B(J + 3) + A(J + 1)
   10   CONTINUE 
        Y(I) = B(1) - XA * B(2)
   20 CONTINUE 
      RETURN 
   30 WRITE (6,40)
   40 FORMAT (' M too big in Sub EVAL2')
      RETURN 
      END 
!
      SUBROUTINE DERIV2 (A, AX, M)
!     Last change: 20 Jan 1990
!     Purpose    :   given the Cheby series A(1),...A(M)
!                    which represents the function
!                    f=Sum(j=1..M) A(j)*T(j-1)
!                    the subroutines computes the coeffs
!                    AX(1).......AX(M) of the function df/dx
!                    Note that AX(M)=0
! 
      DOUBLE PRECISION A(M),AX(M)
      INTEGER          M,J
      AX(M) = 0D0
      IF (M .LE. 0) GO TO 70
      IF (M .GT. 4) GO TO 50
      IF (M .EQ. 1) GO TO 10
      IF (M .EQ. 2) GO TO 20
      IF (M .EQ. 3) GO TO 30
      IF (M .EQ. 4) GO TO 40
   10 AX(1) = 0D0
      RETURN 
   20 AX(1) = 2D0 * A(2)
      RETURN 
   30 AX(1) = 2D0 * A(2)
      AX(2) = 8D0 * A(3)
      RETURN 
   40 AX(1) = 2D0 * A(2) + 6D0 * A(4)
      AX(2) = 8D0 * A(3)
      AX(3) = 12D0 * A(4)
      RETURN 
   50 AX(M - 1) = - 4D0 * (1D0 - M) * A(M)
      AX(M - 2) = - 4D0 * (2D0 - M) * A(M - 1)
      DO 60 J = M - 3, 2, -1
        AX(J) = 4D0 * J * A(J + 1) + AX(J + 2)
   60 CONTINUE 
      AX(1) = 2D0 * A(2) + AX(3) / 2D0
      RETURN 
   70 WRITE (6,80)
   80 FORMAT (' M .LE. 0 in Sub DERIV2')
      RETURN 
      END 
! 
      SUBROUTINE CMESH (N, X, NKIND)
!     Last changed: 20 Jan 1990
!     Purpose :     Chebyshev mesh:
!                   it finds N points X(j) j=1,...N which are the
!                   zeroes of the Cheby polynomial of first kind Tn
!                   or of the Cheby polynomial of second kind Un,
!                   over the interval 0 < x < 1
      DOUBLE PRECISION X(N),PI,y
      INTEGER N,IV,I,NKIND
      PI = 4D0 * DATAN (1D0)
      IF (NKIND .EQ. 1) GO TO 10
      IF (NKIND .EQ. 2) GO TO 30
!     Zeroes of the Chebyshev polynomials of the first kind 
   10 DO 20 I = 1, N
        IV = N + 1 - I
        X(I) = 0.5D0 * (1D0+DCOS((2*IV - 1)*PI/(2D0*N)))
   20 CONTINUE 
      RETURN 
!     Zeroes of the Chebyshev polynomials of the second kind 
   30 DO 40 I = 1, N
        Y = DCOS(PI*I/(N + 1))
        X(I) = 0.5D0 * (Y + 1D0)
   40 CONTINUE 
      RETURN 
      END 
