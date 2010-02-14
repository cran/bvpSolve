c the ELASTICA problem:
c -------- elastica.f -> elastica.dll ------
c compile in R with: system("g77 -shared -o elastica.dll elastica.f")
c or with system("R CMD SHLIB elastica.f")

c  The differential system coded up in a suitable form is as follows:

      SUBROUTINE fsub(NCOMP,X,Z,F,RPAR,IPAR)
      IMPLICIT NONE
      INTEGER NCOMP, IPAR  , I
      DOUBLE PRECISION F, Z, RPAR, X
      DIMENSION Z(*),F(*)
      DIMENSION RPAR(*), IPAR(*)

      F(1)=cos(Z(3))
      F(2)=sin(Z(3))
      F(3)=Z(4)
      F(4)=Z(5)*cos(Z(3))
      F(5)=0

      RETURN
      END

c The analytic Jacobian for the F-function:

      SUBROUTINE dfsub(NCOMP,X,Z,DF,RPAR,IPAR)
      IMPLICIT NONE
      INTEGER NCOMP, IPAR, I, J
      DOUBLE PRECISION X, Z, DF, RPAR
      DIMENSION Z(*),DF(NCOMP,*)
      DIMENSION RPAR(*), IPAR(*)
      CHARACTER (len=50) str

      DO I=1,5
         DO J=1,5
            DF(I,J)=0.D0
         END DO
      END DO

C     dF1/dZ3
      DF(1,3)=-sin(Z(3))
C     dF2/dZ3
      DF(2,3)=cos(Z(3))
C     dF3/dZ4
      DF(3,4)=1.0D0
C     dF4/dZ3
      DF(4,3)=-Z(5)*sin(Z(3))
C     dF4/dZ4
      DF(4,4)=1.0D0
C     dF4/dZ5
      DF(4,5)=cos(Z(3))

      RETURN
      END

c The boundary conditions can be coded up in Fortran 77:

      SUBROUTINE gsub(I,NCOMP,Z,G,RPAR,IPAR)
      IMPLICIT NONE
      INTEGER I, NCOMP, IPAR
      DOUBLE PRECISION Z, RPAR, G
      DIMENSION Z(*)
      DIMENSION RPAR(*), IPAR(*)

c

C     I == the boundary condition "number".
C     The conditions at the left are enumerated first,
C     then the ones at the right.
C     The order of left or right conditions does not matter,
C     but we must be consistent when defining the jacobian
C     of the boundary conditions!

C     We have specified 3 left bcs, and 5 bcs total.
C     This means that:
C     BC(1) = x(0) = 0
C     BC(2) = y(0) = 0
C     BC(3) = kappa(0) = 0
C     BC(4) = y(0.5) = 0
C     BC(5) = phi(0.5) = -pi/2

      IF (I.EQ.1) G=Z(1)
      IF (I.EQ.2) G=Z(2)
      IF (I.EQ.3) G=Z(4)
      IF (I.EQ.4) G=Z(2)
      IF (I.EQ.5) G=Z(3)+1.5707963267948966192313216916397514D0

      RETURN
      END

c The analytic Jacobian for the G-function:

      SUBROUTINE dgsub(I,NCOMP,Z,DG,RPAR,IPAR)
      IMPLICIT NONE
      INTEGER I, NCOMP, IPAR
      DOUBLE PRECISION Z, DG, RPAR
      DIMENSION Z(*),DG(*)
      DIMENSION RPAR(*), IPAR(*)

C     I == the boundary condition "number".
C     The conditions at the left are enumerated first,
C     then the ones at the right.
C     The order of left or right conditions does not matter,
C     but we must be consistent when defining the jacobian
C     of the boundary conditions!

C     We have specified 3 left bcs, and 5 bcs total.
C     This means that:
C     BC(1) = x(0) = 0
C     BC(2) = y(0) = 0
C     BC(3) = kappa(0) = 0
C     BC(4) = y(0.5) = 0
C     BC(5) = phi(0.5) = -pi/2

      DG(1)=0.D0
      DG(2)=0.D0
      DG(3)=0.D0
      DG(4)=0.D0
      DG(5)=0.D0

C     dG1/dZ1
      IF (I.EQ.1) DG(1)=1.D0
C     dG2/dZ2
      IF (I.EQ.2) DG(2)=1.D0
C     dG3/dZ4
      IF (I.EQ.3) DG(4)=1.D0
C     dG4/dZ2
      IF (I.EQ.4) DG(2)=1.D0
C     dG5/dZ3
      IF (I.EQ.5) DG(3)=1.D0

      RETURN
      END
