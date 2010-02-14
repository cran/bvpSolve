c Initialiser for parameter common block
      SUBROUTINE initbnd(bvpparms)
      EXTERNAL bvpparms
      
      DOUBLE PRECISION parms(2)
      COMMON / pars / parms

       CALL bvpparms(2, parms)

      END
      
c derivative function
      SUBROUTINE funbnd(NCOMP,X,Y,F,RPAR,IPAR)
      IMPLICIT NONE
      INTEGER NCOMP, IPAR(*), I
      DOUBLE PRECISION F(2), Y(2), RPAR(*), X

      DOUBLE PRECISION a, p
      COMMON / pars / a, p

        F(1)= Y(2)
        F(2)= - a * p *Y(1)/(p+ x*x)**2

      END

c The analytic Jacobian for the derivative-function:
      SUBROUTINE dfbnd(NCOMP,X,Y,DF,RPAR,IPAR)
      IMPLICIT NONE
      INTEGER NCOMP, IPAR(*), I, J
      DOUBLE PRECISION X, Y(2), DF(2,2), RPAR(*)

      DOUBLE PRECISION a, p
      COMMON / pars / a, p

        DF(1,1)=0.D0
        DF(1,2)=1.D0
        DF(2,1)= - a *p /(p+x*x)**2
        DF(2,2)=0.D0

      END

c The boundary conditions:
      SUBROUTINE gbnd(I,NCOMP,Y,G,RPAR,IPAR)
      IMPLICIT NONE
      INTEGER I, NCOMP, IPAR(*)
      DOUBLE PRECISION Y(2), RPAR(*), G

      DOUBLE PRECISION a, p
      COMMON / pars / a, p

        IF (I.EQ.1) THEN
          G=Y(1) + 0.1 / sqrt(p+0.01)
        ELSE IF (I.EQ.2) THEN
          G=Y(1) - 0.1 / sqrt(p+0.01)
        ENDIF
      
      END

c The analytic Jacobian for the boundaries:
      SUBROUTINE dgbnd(I,NCOMP,Y,DG,RPAR,IPAR)
      IMPLICIT NONE
      INTEGER I, NCOMP, IPAR(*)
      DOUBLE PRECISION Y(2), DG(2), RPAR(*)

        DG(1)=1.D0
        DG(2)=0.D0

      END
