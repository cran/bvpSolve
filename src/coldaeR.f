c ===================================================================================
* karline: to make this code compatible with R:
* 1. change all write statements into rprint statements
* 2. changed all ( ,1) declarations into (,*)
* 3. changed interface to fsub, gsub to make it compatible with TWPBVPC
* 3. added rpar, ipar to calls of guess
* 4. added counters
* 5. renamed ipar -> iset and solutn -> guess
* 6. renamed comon blocks
c ===================================================================================

C.... *****************************************************************
C  this package solves boundary value problems for
C  ordinary differential equations with constraints,
C  as described below. for more see [1].
C
C  COLDAE is a modification of the package COLNEW by bader and
C  ascher [4], which in turn is a modification of COLSYS by ascher,
C  christiansen and russell [3]. it does what colnew does plus
C  optionally solving semi-explicit differential-algebraic equations
C  with index at most 2.
C**********************************************************************
C----------------------------------------------------------------------
C                            p a r t  1
C        main storage allocation and program control subroutines
C----------------------------------------------------------------------
C
      SUBROUTINE COLDAE (NCOMP, NY, M, ALEFT, ARIGHT, ZETA, ISET, LTOL,
     1                   TOL, FIXPNT, ISPACE, FSPACE, IFLAG,
     2                   FSUB, DFSUB, GSUB, DGSUB, GUESS, RPAR, IPAR,
     3                   ICOUNT)
C
C
C**********************************************************************
C
C     written by
C                  uri ascher and ray spiteri,
C                            department of computer science,
C                            university of british columbia,
C                            vancouver, b. c., canada   v6t 1z2
C
C**********************************************************************
C
C     purpose
C
C     this package solves a multi-point boundary value
C     problem for a mixed order system of ode-s with constraints
C     given by
C
C          (m(i))
C         u       =  f  ( x; z(u(x)), y(x) )   i = 1, ... ,ncomp
C          i          i
C
C             0   =  f  ( x; z(u(x)), y(x) )   i = ncomp+1,...,ncomp+ny
C                     i
C
C                                          aleft .lt. x .lt. aright,
C
C
C         g  ( zeta(j); z(u(zeta(j))) ) = 0   j = 1, ... ,mstar
C          j
C                                    mstar = m(1)+m(2)+...+m(ncomp),
C
C
C         where                          t                       t
C               u = (u , u , ... ,u     ) , y = (y ,y , ... ,y  )
C                     1   2        ncomp          1  2        ny
C
C         is the exact solution vector: u are the differential solution
C         components and y are the algebraic solution components.
C
C                (mi)
C               u     is the mi=m(i) th  derivative of u
C                i                                      i
C
C                                  (1)        (m1-1)       (mncomp-1)
C               z(u(x)) = ( u (x),u  (x),...,u    (x),...,u      (x) )
C                            1     1          1            ncomp
C
C                f (x,z(u),y)   is a (generally) nonlinear function of
C                 i
C                             z(u)=z(u(x)) and y=y(x).
C
C                g (zeta(j);z(u))  is a (generally) nonlinear function
C                 j
C                               used to represent a boundary condition.
C
C         the boundary points satisfy
C               aleft .le. zeta(1) .le. .. .le. zeta(mstar) .le. aright
C
C         the orders mi of the differential equations satisfy
C                            1 .le. m(i) .le. 4.
C
C         regarding the dae, note that:
C         i)  with ny=0, the code is essentially identical to the ode
C             code  colnew
C         ii) no explicit checking of the index of the problem
C             is provided. if the index is &gt; 2 then the code will
C             not work well.
C         iii) the constraints are treated like de-s of order 0 and
C             correspondingly the approximation to y is sought in
C             a piecewise discontinuous polynomial space.
C         iv) the number of boundary conditions required is independent
C             of the index. it is the user's responsibility to ensure
C             that these conditions are consistent with the constraints.
C             the conditions at the left end point aleft  must include
C             a subset equivalent to specifying the index-2
C             constraints there.
C         v)  for an index-2 problem in hessenberg form, the projected
C             collocation method of ascher and petzold [2] is used.
C         vi) if the constraints are of a mixed type (and possibly
C             a mixed index) then coldae can transform and project
C             appropriately -- see description of ISET(12).
C
C**********************************************************************
C
C     method
C
C        the method used to approximate the solution u is
C     collocation at gaussian points, requiring m(i)-1 continuous
C     derivatives in the i-th component, i = 1, ..., ncomp.
C     here, k is the number of collocation points (stages) per
C     subinterval and is chosen such that k .ge. max m(i).
C     a runge-kutta-monomial solution representation is utilized.
C     for hessenberg index-2 daes, a projection on the constraint
C     manifold at each interval's end is used [2].
C
C     references
C
C     [1] u. ascher and r. spiteri,
C         collocation software for boundary value differential-
C         algebraic equations,
C         siam j. scient. stat. comput., to appear.
C
C     [2] u. ascher and l. petzold,
C         projected implicit runge-kutta methods for differential-
C         algebraic equations,
C         siam j. num. anal. 28 (1991), 1097-1120.
C
C     [3] u. ascher, j. christiansen and r.d. russell,
C         collocation software for boundary-value odes,
C         acm trans. math software 7 (1981), 209-222.
C
C     [4] g. bader and u. ascher,
C         a new basis implementation for a mixed order
C         boundary value ode solver,
C         siam j. scient. stat. comput. 8 (1987), 483-500.
C
C     [5] u. ascher, j. christiansen and r.d. russell,
C         a collocation solver for mixed order
C         systems of boundary value problems,
C         math. comp. 33 (1979), 659-679.
C
C     [6] u. ascher, j. christiansen and r.d. russell,
C         colsys - a collocation code for boundary
C         value problems,
C         lecture notes comp.sc. 76, springer verlag,
C         b. childs et. al. (eds.) (1979), 164-185.
C
C     [7] c. deboor and r. weiss,
C         solveblok: a package for solving almost block diagonal
C         linear systems,
C         acm trans. math. software 6 (1980), 80-87.
C
C
C**********************************************************************
C
C     ***************     input to coldae     ***************
C
C     variables
C
C     ncomp - no. of differential equations   (ncomp .le. 20)
C
C     ny   - no. of constraints     (ny .le. 20)
C
C     m(j) - order of the j-th differential equation
C            ( mstar = m(1) + ... + m(ncomp) .le. 40 )
C
C     aleft - left end of interval
C
C     aright - right end of interval
C
C     zeta(j) - j-th side condition point (boundary point). must
C               have  zeta(j) .le. zeta(j+1). all side condition
C               points must be mesh points in all meshes used,
C               see description of ISET(11) and fixpnt below.
C
C     ISET - an integer array dimensioned at least 11.
C            a list of the parameters in ISET and their meaning follows
C            some parameters are renamed in coldae; these new names are
C            given in parentheses.
C
C     ISET(1)     ( = nonlin )
C             = 0 if the problem is linear
C             = 1 if the problem is nonlinear
C
C     ISET(2) = no. of collocation points per subinterval  (= k )
C               where max m(i) .le.  k .le. 7 . if ISET(2)=0 then
C               coldae sets  k = max ( max m(i)+1, 5-max m(i) )
C
C     ISET(3) = no. of subintervals in the initial mesh  ( = n ).
C               if ISET(3) = 0 then coldae arbitrarily sets n = 5.
C
C     ISET(4) = no. of solution and derivative tolerances.  ( = ntol )
C               we require  0 .lt. ntol .le. mstar.
C
C     ISET(5) = dimension of fspace.     ( = ndimf )
C
C     ISET(6) = dimension of ispace.     ( = ndimi )
C
C     ISET(7) -  output control ( = iprint )
C              = -1 for full diagnostic printout
C              = 0 for selected printout
C              = 1 for no printout
C
C     ISET(8)     ( = iread )
C             = 0 causes coldae to generate a uniform initial mesh.
C             = 1 if the initial mesh is provided by the user.  it
C                 is defined in fspace as follows:  the mesh
C                 aleft=x(1).lt.x(2).lt. ... .lt.x(n).lt.x(n+1)=aright
C                 will occupy  fspace(1), ..., fspace(n+1). the
C                 user needs to supply only the interior mesh
C                 points  fspace(j) = x(j), j = 2, ..., n.
C             = 2 if the initial mesh is supplied by the user
C                 as with ISET(8)=1, and in addition no adaptive
C                 mesh selection is to be done.
C
C     ISET(9)     ( = iguess )
C             = 0 if no initial guess for the solution is
C                 provided.
C             = 1 if an initial guess is provided by the user
C                 in subroutine  guess.
C             = 2 if an initial mesh and approximate solution
C                 coefficients are provided by the user in  fspace.
C                 (the former and new mesh are the same).
C             = 3 if a former mesh and approximate solution
C                 coefficients are provided by the user in fspace,
C                 and the new mesh is to be taken twice as coarse;
C                 i.e.,every second point from the former mesh.
C             = 4 if in addition to a former initial mesh and
C                 approximate solution coefficients, a new mesh
C                 is provided in fspace as well.
C                 (see description of output for further details
C                 on iguess = 2, 3, and 4.)
C
C     ISET(10)= -1 if the first relax factor is RSTART
C        (use for an extra sensitive nonlinear problem only)
C          =  0 if the problem is regular
C          =  1 if the newton iterations are not to be damped
C                 (use for initial value problems).
C          =  2 if we are to return immediately upon  (a) two
C                 successive nonconvergences, or  (b) after obtaining
C                 error estimate for the first time.
C
C     ISET(11)= no. of fixed points in the mesh other than aleft
C               and aright. ( = nfxpnt , the dimension of fixpnt)
C               the code requires that all side condition points
C               other than aleft and aright (see description of
C               zeta ) be included as fixed points in fixpnt.
C
C     ISET(12)    ( = index )
C                 this parameter is ignored if ny=0.
C             = 0 if the index of the dae is not as per one of the
C                 following cases
C             = 1 if the index of the dae is 1. in this case the
C                 ny x ny jacobian matrix of the constraints with
C                 respect to the algebraic unknowns, i.e.
C                 df(i,j), i=ncomp+1,...,ncomp+ny,
C                          j=mstar+1,...,mstar+ny
C                 (see description of dfsub below)
C                 is nonsingular wherever it is evaluated. this
C                 allows usual collocation to be safely used.
C             = 2 if the index of the dae is 2 and it is in Hessenberg
C                 form. in this case the
C                 ny x ny jacobian matrix of the constraints with
C                 respect to the algebraic unknowns is 0, and the
C                 ny x ny matrix  CB  is nonsingular, where
C                 C(i,j) = df(i+ncomp, m(j)), i=1,...,ny,
C                                             j=1,...,ncomp
C                 B(i,j) = df(i,j+mstar),     i=1,...,ncomp,
C                                             j=1,...,ny
C                 the projected collocation method described in [2]
C                 is then used.
C             in case of ISET(12)=0 and ny &gt; 0, coldae determines the
C             appropriate projection needed at the right end of each
C             mesh subinterval using SVD. this is the most expensive
C             and most general option.
C
C     ltol  -  an array of dimension ISET(4). ltol(j) = l  specifies
C              that the j-th tolerance in  tol  controls the error
C              in the l-th component of z(u).   also require that
C              1.le.ltol(1).lt.ltol(2).lt. ... .lt.ltol(ntol).le.mstar
C
C     tol    - an array of dimension ISET(4). tol(j) is the
C              error tolerance on the ltol(j) -th component
C              of z(u). thus, the code attempts to satisfy
C              for j=1,...,ntol  on each subinterval
C              abs(z(v)-z(u))       .le. tol(j)*abs(z(u))       +tol(j)
C                            ltol(j)                     ltol(j)
C
C              if v(x) is the approximate solution vector.
C
C     fixpnt - an array of dimension ISET(11).   it contains
C              the points, other than aleft and aright, which
C              are to be included in every mesh.
C
C     ispace - an integer work array of dimension ISET(6).
C              its size provides a constraint on nmax,
C              the maximum number of subintervals. choose
C              ISET(6) according to the formula
C                      ISET(6)  .ge.  nmax*nsizei
C                where
C                      nsizei = 3 + kdm
C                with
C                      kdm = kdy + mstar  ;  kdy = k * (ncomp+ny) ;
C                      nrec = no. of right end boundary conditions.
C
C
C     fspace - a real work array of dimension ISET(5).
C              its size provides a constraint on nmax.
C              choose ISET(5) according to the formula
C                      ISET(5)  .ge.  nmax*nsizef
C                where
C                      nsizef = 4 + 3 * mstar + (5+kdy) * kdm +
C                              (2*mstar-nrec) * 2*mstar +
C                              ncomp*(mstar+ny+2) + kdy.
C
C
C     iflag - the mode of return from coldae.
C           = 1 for normal return
C           = 0 if the collocation matrix is singular.
C           =-1 if the expected no. of subintervals exceeds storage
C               specifications.
C           =-2 if the nonlinear iteration has not converged.
C           =-3 if there is an input data error.
C
C
C**********************************************************************
C
C     *************    user supplied subroutines   *************
C
C
C     the following subroutines must be declared external in the
C     main program which calls coldae.
C
C
C     fsub  - name of subroutine for evaluating f(x,z(u(x)),y(x)) =
C                               t
C             (f ,...,f        )  at a point x in (aleft,aright).  it
C               1      ncomp+ny
C             should have the heading
C
C                       subroutine fsub (x , z , y, f)
C
C             where f is the vector containing the value of fi(x,z(u),y)
C             in the i-th component and y(x) = (y(1),...y(ny)) and
C                                       z(u(x))=(z(1),...,z(mstar))
C             are defined as above under  purpose .
C
C
C     dfsub - name of subroutine for evaluating the jacobian of
C             f(x,z(u),y) at a point x.  it should have the heading
C
C                       subroutine dfsub (x , z , y, df)
C
C             where z(u(x)) and y(x) are defined as for fsub and the
C             (ncomp+ny) by (mstar+ny)
C             array  df  should be filled by the partial deriv-
C             atives of f, viz, for a particular call one calculates
C                          df(i,j) = dfi / dzj, i=1,...,ncomp+ny
C                                               j=1,...,mstar,
C                          df(i,mstar+j) = dfi / dyj, i=1,...,ncomp+ny
C                                               j=1,...,ny.
C
C
C     gsub  - name of subroutine for evaluating the i-th component of
C             g(x,z(u(x))) = g (zeta(i),z(u(zeta(i)))) at a point x =
C                             i
C             zeta(i) where 1.le.i.le.mstar. it should have the heading
C
C                       subroutine gsub (i , z , g)
C
C             where z(u) is as for fsub, and i and g=g  are as above.
C                                                     i
C             note that in contrast to f in  fsub , here
C             only one value per call is returned in g.
C             also, g is independent of y and, in case of a higher
C             index dae (index = 2), should include the constraints
C             sampled at aleft  plus independent additional constraints,
C             or an equivalent set.
C
C
C     dgsub - name of subroutine for evaluating the i-th row of
C             the jacobian of g(x,z(u(x))).  it should have the heading
C
C                       subroutine dgsub (i , z , dg)
C
C             where z(u) is as for fsub, i as for gsub and the mstar-
C             vector dg should be filled with the partial derivatives
C             of g, viz, for a particular call one calculates
C                   dg(i,j) = dgi / dzj      j=1,...,mstar.
C
C
C     guess - name of subroutine to evaluate the initial
C             approximation for  z(u(x)), y(x) and for dmval(u(x))= vector
C             of the mj-th derivatives of u(x). it should have the
C             heading
C
C                       subroutine guess (x , z , y, dmval)
C
C             note that this subroutine is needed only if using
C             ISET(9) = 1, and then all  mstar  components of z,
C             ny components of  y
C             and  ncomp  components of  dmval  should be specified
C             for any x,  aleft .le. x .le. aright .
C
C
C**********************************************************************
C
C     ************   use of output from coldae   ************
C
C                 ***   solution evaluation   ***
C
C     on return from coldae, the arrays fspace and ispace
C     contain information specifying the approximate solution.
C     the user can produce the solution vector  z( u(x) ), y(x)  at
C     any point x, aleft .le. x .le. aright, by the statement,
C
C           call appsln (x, z, y, fspace, ispace)
C
C     when saving the coefficients for later reference, only
C     ispace(1),...,ispace(8+ncomp)    and
C     fspace(1),...,fspace(ispace(8))    need to be saved as
C     these are the quantities used by appsln.
C
C
C                 ***   simple continuation   ***
C
C
C     a formerly obtained solution can easily be used as the
C     first approximation for the nonlinear iteration for a
C     new problem by setting   (iguess =) ISET(9) = 2, 3 or 4.
C
C     if the former solution has just been obtained then the
C     values needed to define the first approximation are
C     already in ispace and fspace.
C     alternatively, if the former solution was obtained in a
C     previous run and its coefficients were saved then those
C     coefficients must be put back into
C     ispace(1),..., ispace(8+ncomp)    and
C     fspace(1),..., fspace(ispace(8)).
C
C     for ISET(9) = 2 or 3 set ISET(3) = ispace(1) ( = the
C     size of the previous mesh ).
C
C     for ISET(9) = 4 the user specifies a new mesh of n subintervals
C     as follows.
C     the values in  fspace(1),...,fspace(ispace(8))  have to be
C     shifted by n+1 locations to  fspace(n+2),..,fspace(ispace(8)+n+1)
C     and the new mesh is then specified in fspace(1),..., fspace(n+1).
C     also set ISET(3) = n.
C
C
C**********************************************************************
C
C     ***************      package subroutines      ***************
C
C     the following description gives a brief overview of how the
C     procedure is broken down into the subroutines which make up
C     the package called  coldae . for further details the
C     user should refer to documentation in the various subroutines
C     and to the references cited above.
C
C     the subroutines fall into four groups:
C
C part 1 - the main storage allocation and program control subr
C
C     coldae - tests input values, does initialization and breaks up
C              the work areas, fspace and ispace, into the arrays
C              used by the program.
C
C     contrl - is the actual driver of the package. this routine
C              contains the strategy for nonlinear equation solving.
C
C     skale  - provides scaling for the control
C              of convergence in the nonlinear iteration.
C
C
C part 2 - mesh selection and error estimation subroutines
C
C     consts - is called once by  coldae  to initialize constants
C              which are used for error estimation and mesh selection.
C
C     newmsh - generates meshes. it contains the test to decide
C              whether or not to redistribute a mesh.
C
C     errchk - produces error estimates and checks against the
C              tolerances at each subinterval
C
C
C part 3 - collocation system set-up subroutines
C
C     lsyslv - controls the set-up and solution of the linear
C              algebraic systems of collocation equations which
C              arise at each newton iteration.
C
C     gderiv - is used by lsyslv to set up the equation associated
C              with a side condition point.
C
C     vwblok - is used by lsyslv to set up the equation(s) associated
C              with a collocation point.
C
C     gblock - is used by lsyslv to construct a block of the global
C              collocation matrix or the corresponding right hand
C              side.
C
C
C part 4 - service subroutines
C
C     appsln - sets up a standard call to  approx .
C
C     approx - evaluates a piecewise polynomial solution.
C
C     rkbas  - evaluates the mesh independent runge-kutta basis
C
C     vmonde - solves a vandermonde system for given right hand
C              side
C
C     horder - evaluates the highest order derivatives of the
C              current collocation solution used for mesh refinement.
C
C
C part 5 - linear algebra  subroutines
C
C     to solve the global linear systems of collocation equations
C     constructed in part 3,  coldae  uses a column oriented version
C     of the package  solveblok originally due to de boor and weiss.
C
C     to solve the linear systems for static parameter condensation
C     in each block of the collocation equations, the linpack
C     routines  dgefa and  dgesl  are included. but these
C     may be replaced when solving problems on vector processors
C     or when solving large scale sparse jacobian problems.
C
C----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION M(*), ZETA(*), ISET(*), LTOL(*), TOL(*), DUMMY(1),
     1    FIXPNT(*), ISPACE(*), FSPACE(*), RPAR(*), IPAR(*), ICOUNT(*)
      DIMENSION DUMMY2(840)
C
      COMMON /DAEOUT/ PRECIS, IOUT, IPRINT
      COMMON /COLLOC/ RHO(7), COEF(49)
      COMMON /DAEORD/ K, NC, NNY, NCY, MSTAR, KD, KDY, MMAX, MT(20)
      COMMON /DAEAPR/ N, NOLD, NMAX, NZ, NDMZ
      COMMON /DAEMSH/ MSHFLG, MSHNUM, MSHLMT, MSHALT
      COMMON /DAESID/ TZETA(40), TLEFT, TRIGHT, IZETA, IDUM
      COMMON /DAENLN/ NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX
      COMMON /DAEEST/ TTL(40), WGTMSH(40), WGTERR(40), TOLIN(40),
     1                ROOT(40), JTOL(40), LTTOL(40), NTOL
C
      EXTERNAL FSUB, DFSUB, GSUB, DGSUB, GUESS

      integer nfunc, njac, nstep, nbound, njacbound
      common/CDAEdiag/nfunc, njac, nstep, nbound, njacbound

C
C*********************************************************************
C
C     the actual subroutine coldae serves as an interface with
C     the package of subroutines referred to collectively as
C     coldae. the subroutine serves to test some of the input
C     parameters, rename some of the parameters (to make under-
C     standing of the coding easier), to do some initialization,
C     and to break the work areas fspace and ispace up into the
C     arrays needed by the program.
C
C**********************************************************************
C     intialise counters
      nfunc = 0
      njac = 0
      nstep = 0
      nbound = 0
      njacbound = 0

C
C...  specify machine dependent output unit  iout  and compute machine
C...  dependent constant  precis = 100 * machine unit roundoff
C
      IOUT = 6
      PRECIS = 1.D0
   10 PRECIS = PRECIS / 2.D0
      PRECP1 = PRECIS + 1.D0
      IF ( PRECP1 .GT. 1.D0 )                       GO TO 10
      PRECIS = PRECIS * 100.D0
C
C...  in case incorrect input data is detected, the program returns
C...  immediately with iflag=-3.
C
      IFLAG = -3
      NCY = NCOMP + NY
      IF ( NCOMP .LT. 0 .OR. NCOMP .GT. 20 )        RETURN
      IF ( NY .LT. 0 .OR. NY .GT. 20 )              RETURN
      IF ( NCY .LT. 1 .OR. NCY .GT. 40 )            RETURN
      DO 20 I=1,NCOMP
       IF ( M(I) .LT. 1 .OR. M(I) .GT. 4 )        RETURN
   20 CONTINUE
C
C...  rename some of the parameters and set default values.
C
      NONLIN = ISET(1)
      K = ISET(2)
      N = ISET(3)
      IF ( N .EQ. 0 )  N = 5
      IREAD = ISET(8)
      IGUESS = ISET(9)
      IF ( NONLIN .EQ. 0 .AND. IGUESS .EQ. 1 )  IGUESS = 0
      IF ( IGUESS .GE. 2 .AND. IREAD .EQ. 0 )   IREAD = 1
      ICARE = ISET(10)
      NTOL = ISET(4)
      NDIMF = ISET(5)
      NDIMI = ISET(6)
      NFXPNT = ISET(11)
      IPRINT = ISET(7)
      INDEX = ISET(12)
      IF (NY .EQ. 0) INDEX = 0
      MSTAR = 0
      MMAX = 0
      DO  30 I = 1, NCOMP
       MMAX = MAX0 ( MMAX, M(I) )
       MSTAR = MSTAR + M(I)
       MT(I) = M(I)
   30 CONTINUE
      IF ( K .EQ. 0 )   K = MAX0( MMAX + 1 , 5 - MMAX )
      DO 40 I = 1, MSTAR
   40 TZETA(I) = ZETA(I)
      DO 50 I = 1, NTOL
       LTTOL(I) = LTOL(I) 
   50 TOLIN(I) = TOL(I)
      TLEFT = ALEFT
      TRIGHT = ARIGHT
      NC = NCOMP
      NNY = NY
      KD = K * NCOMP
      KDY = K * NCY
C
C...  print the input data for checking.
C
C Karline toggled a lot of printing off
C ...
C
C...  check for correctness of data
C
      IF ( K .LT. 0 .OR. K .GT. 7 )                 RETURN
      IF ( N .LT. 0 )                               RETURN
      IF ( IREAD .LT. 0 .OR. IREAD .GT. 2 )         RETURN
      IF ( IGUESS .LT. 0 .OR. IGUESS .GT. 4 )       RETURN
      IF ( ICARE .LT. -1 .OR. ICARE .GT. 2 )       RETURN
      IF ( INDEX .LT. 0 .OR. INDEX .GT. 2 )         RETURN
      IF ( NTOL .LT. 0 .OR. NTOL .GT. MSTAR )       RETURN
      IF ( NFXPNT .LT. 0 )                          RETURN
      IF ( IPRINT .LT. (-1) .OR. IPRINT .GT. 1 )    RETURN
      IF ( MSTAR .LT. 0 .OR. MSTAR .GT. 40 )        RETURN
      IP = 1
      DO 100 I = 1, MSTAR
      IF ( DABS(ZETA(I) - ALEFT) .LT. PRECIS .OR.
     1     DABS(ZETA(I) - ARIGHT) .LT. PRECIS )     GO TO 100
   90 IF ( IP .GT. NFXPNT )                         RETURN
      IF ( ZETA(I) - PRECIS .LT. FIXPNT(IP) )     GO TO 95
      IP = IP + 1
      GO TO 90
   95 IF ( ZETA(I) + PRECIS .LT. FIXPNT(IP) )       RETURN
  100 CONTINUE
C
C...  set limits on iterations and initialize counters.
C...  limit = maximum number of newton iterations per mesh.
C...  see subroutine  newmsh  for the roles of  mshlmt , mshflg ,
C...  mshnum , and  mshalt .
C
      MSHLMT = 3
      MSHFLG = 0
      MSHNUM = 1
      MSHALT = 1
      LIMIT = 40
C
C...  compute the maxium possible n for the given sizes of
C...  ispace  and  fspace.
C
      NREC = 0
      DO 110 I = 1, MSTAR
       IB = MSTAR + 1 - I
       IF ( ZETA(IB) .GE. ARIGHT )  NREC = I
  110 CONTINUE
      NFIXI = MSTAR
      NSIZEI = 3 + KDY + MSTAR
      NFIXF = NREC * (2*MSTAR) + 5 * MSTAR + 3
      NSIZEF = 4 + 3 * MSTAR + (KDY+5) * (KDY+MSTAR) +
     1(2*MSTAR-NREC) * 2*MSTAR + (MSTAR+NY+2)*NCOMP + KDY
      NMAXF = (NDIMF - NFIXF) / NSIZEF
      NMAXI = (NDIMI - NFIXI) / NSIZEI
      IF ( IPRINT .LT. 1 )  THEN
       CALL Rprinti1('The maximum number of subintervals is min',nmaxF)
       CALL Rprinti1('The maximum number allowed from ispace', NMAXI)
      ENDIF
      NMAX = MIN0( NMAXF, NMAXI )
      IF ( NMAX .LT. N )                            RETURN
      IF ( NMAX .LT. NFXPNT+1 )                     RETURN
      IF (NMAX .LT. 2*NFXPNT+2 .AND. IPRINT .LT. 1) THEN
      CALL Rprint('Insufficient space to double mesh for err. estimate')
      ENDIF
C
C...  generate pointers to break up  fspace  and  ispace .
C
      LXI = 1
      LG = LXI + NMAX + 1
      LXIOLD = LG + 2*MSTAR * (NMAX * (2*MSTAR-NREC) + NREC)
      LW     = LXIOLD + NMAX + 1
      LV     = LW + KDY**2 * NMAX
      LFC    = LV + MSTAR * KDY * NMAX
      LZ     = LFC + (MSTAR+NY+2) * NCOMP * NMAX
      LDMZ   = LZ + MSTAR * (NMAX + 1)
      LDMV   = LDMZ + KDY* NMAX
      LDELZ  = LDMV + KDY * NMAX
      LDELDZ = LDELZ + MSTAR * (NMAX + 1)
      LDQZ   = LDELDZ + KDY * NMAX
      LDQDMZ = LDQZ + MSTAR * (NMAX + 1)
      LRHS   = LDQDMZ + KDY * NMAX
      LVALST = LRHS   + KDY * NMAX + MSTAR
      LSLOPE = LVALST + 4 * MSTAR * NMAX
      LACCUM = LSLOPE + NMAX
      LSCL   = LACCUM + NMAX + 1
      LDSCL  = LSCL + MSTAR * (NMAX + 1)
      LPVTG = 1
      LPVTW = LPVTG + MSTAR * (NMAX + 1)
      LINTEG = LPVTW + KDY * NMAX
C
C...  if  iguess .ge. 2, move  xiold, z, and  dmz  to their proper
C...  locations in  fspace.
C
      IF ( IGUESS .LT. 2 )                          GO TO 160
      NOLD = N
      IF (IGUESS .EQ. 4)  NOLD = ISPACE(1)
      NZ = MSTAR * (NOLD + 1)
      NDMZ = KDY * NOLD
      NP1 = N + 1
      IF ( IGUESS .EQ. 4 )  NP1 = NP1 + NOLD + 1
      DO 120 I=1,NZ
  120 FSPACE( LZ+I-1 )  =  FSPACE( NP1+I )
      IDMZ = NP1 + NZ
      DO 125 I=1,NDMZ
  125 FSPACE( LDMZ+I-1 )  =  FSPACE( IDMZ+I )
      NP1 = NOLD + 1
      IF ( IGUESS .EQ. 4 )                          GO TO 140
      DO 130 I=1,NP1
  130 FSPACE( LXIOLD+I-1 )  =  FSPACE( LXI+I-1 )
      GO TO 160
  140 DO 150 I=1,NP1
  150 FSPACE( LXIOLD+I-1 )  =  FSPACE( N+1+I )
  160 CONTINUE
C
C...  initialize collocation points, constants, mesh.
C
      CALL CONSTS_DAE ( K, RHO, COEF )
      IF (NY.EQ.0) THEN
       NYCB = 1
      ELSE
       NYCB = NY
      ENDIF
      CALL NEWMSH_DAE (3+IREAD, FSPACE(LXI), FSPACE(LXIOLD), DUMMY,
     1             DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, NFXPNT, FIXPNT,
     2           DUMMY2, DFSUB, DUMMY2, DUMMY2, NCOMP, NYCB,RPAR,IPAR)
C
C...  determine first approximation, if the problem is nonlinear.
C
      IF (IGUESS .GE. 2)                            GO TO 230
      NP1 = N + 1
      DO 210 I = 1, NP1
  210 FSPACE( I + LXIOLD - 1 ) = FSPACE( I + LXI - 1 )
      NOLD = N
      IF ( NONLIN .EQ. 0  .OR. IGUESS .EQ. 1 )      GO TO 230
C
C...  system provides first approximation of the solution.
C...  choose z(j) = 0  for j=1,...,mstar.
C
      DO 220 I=1, NZ
  220 FSPACE( LZ-1+I ) = 0.D0
      DO 225 I=1, NDMZ
  225 FSPACE( LDMZ-1+I ) = 0.D0
  230 CONTINUE
      IF (IGUESS .GE. 2)  IGUESS = 0
      CALL CONTRL_DAE(FSPACE(LXI),FSPACE(LXIOLD),FSPACE(LZ),
     +     FSPACE(LDMZ),FSPACE(LDMV),
     1     FSPACE(LRHS),FSPACE(LDELZ),FSPACE(LDELDZ),FSPACE(LDQZ),
     2     FSPACE(LDQDMZ),FSPACE(LG),FSPACE(LW),FSPACE(LV),FSPACE(LFC),
     3     FSPACE(LVALST),FSPACE(LSLOPE),FSPACE(LSCL),FSPACE(LDSCL),
     4     FSPACE(LACCUM),ISPACE(LPVTG),ISPACE(LINTEG),ISPACE(LPVTW),
     5     NFXPNT,FIXPNT,IFLAG,FSUB,DFSUB,GSUB,DGSUB,GUESS,RPAR,IPAR)
C
C...  prepare output
C
      ISPACE(1) = N
      ISPACE(2) = K
      ISPACE(3) = NCOMP
      ISPACE(4) = NY
      ISPACE(5) = MSTAR
      ISPACE(6) = MMAX
      ISPACE(7) = NZ + NDMZ + N + 2
      K2 = K * K
      ISPACE(8) = ISPACE(7) + K2 - 1
      DO 240 I = 1, NCOMP
  240 ISPACE(8+I) = M(I)
      DO 250 I = 1, NZ
  250 FSPACE( N+1+I ) = FSPACE( LZ-1+I )
      IDMZ = N + 1 + NZ
      DO 255 I = 1, NDMZ
  255 FSPACE( IDMZ+I ) = FSPACE( LDMZ-1+I )
      IC = IDMZ + NDMZ
      DO 258 I = 1, K2
  258 FSPACE( IC+I ) = COEF(I)
  259 icount(1) = nfunc
      icount(2) = njac
      icount(3) = nbound
      icount(4) = njacbound
      icount(5) = nstep
  
      RETURN
C----------------------------------------------------------------------
      END
      SUBROUTINE CONTRL_DAE(XI, XIOLD, Z, DMZ, DMV, RHS, DELZ, DELDMZ,
     1           DQZ, DQDMZ, G, W, V, FC, VALSTR, SLOPE, SCALE, DSCALE,
     2           ACCUM, IPVTG, INTEGS, IPVTW, NFXPNT, FIXPNT, IFLAG,
     3           FSUB, DFSUB, GSUB, DGSUB, GUESS,RPAR,IPAR)
C
C**********************************************************************
C
C   purpose
C     this subroutine is the actual driver.  the nonlinear iteration
C     strategy is controlled here ( see [6] ). upon convergence, errchk
C     is called to test for satisfaction of the requested tolerances.
C
C   variables
C
C     check  - maximum tolerance value, used as part of criteria for
C              checking for nonlinear iteration convergence
C     relax  - the relaxation factor for damped newton iteration
C     relmin - minimum allowable value for relax  (otherwise the
C              jacobian is considered singular).
C     rlxold - previous relax
C     rstart - initial value for relax when problem is sensitive
C     ifrz   - number of fixed jacobian iterations
C     lmtfrz - maximum value for ifrz before performing a reinversion
C     iter   - number of iterations (counted only when jacobian
C              reinversions are performed).
C     xi     - current mesh
C     xiold  - previous mesh
C     ipred  = 0  if relax is determined by a correction
C            = 1  if relax is determined by a prediction
C     ifreez = 0  if the jacobian is to be updated
C            = 1  if the jacobian is currently fixed (frozen)
C     iconv  = 0  if no previous convergence has been obtained
C            = 1  if convergence on a previous mesh has been obtained
C     icare  =-1  no convergence occurred (used for regular problems)
C            = 0  a regular problem
C            = 1  no damped newton
C            = 2  used for continuation (see description of ISET(10)
C                 in coldae).
C     rnorm  - norm of rhs (right hand side) for current iteration
C     rnold  - norm of rhs for previous iteration
C     anscl  - scaled norm of newton correction
C     anfix  - scaled norm of newton correction at next step
C     anorm  - scaled norm of a correction obtained with jacobian fixed
C     nz     - number of components of  z  (see subroutine approx)
C     ndmz   - number of components of  dmz  (see subroutine approx)
C     imesh  - a control variable for subroutines newmsh and errchk
C            = 1  the current mesh resulted from mesh selection
C                 or is the initial mesh.
C            = 2  the current mesh resulted from doubling the
C                 previous mesh
C
C**********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XI(*), XIOLD(*), Z(*), DMZ(*), RHS(*), DMV(*)
      DIMENSION G(*), W(*), V(*), VALSTR(*), SLOPE(*), ACCUM(*)
      DIMENSION DELZ(*), DELDMZ(*), DQZ(*), DQDMZ(*) , FIXPNT(*)
      DIMENSION DUMMY(1), SCALE(*), DSCALE(*), FC(*), DF(800)
      DIMENSION FCSP(40,60), CBSP(20,20),RPAR(*),IPAR(*)
      DIMENSION INTEGS(*), IPVTG(*), IPVTW(*)
C
      COMMON /DAEOUT/ PRECIS, IOUT, IPRINT
      COMMON /DAEORD/ K, NCOMP, NY, NCY, MSTAR, KD, KDY, MMAX, M(20)
      COMMON /DAEAPR/ N, NOLD, NMAX, NZ, NDMZ
      COMMON /DAEMSH/ MSHFLG, MSHNUM, MSHLMT, MSHALT
      COMMON /DAESID/ ZETA(40), ALEFT, ARIGHT, IZETA, IDUM
      COMMON /DAENLN/ NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX
      COMMON /DAEEST/ TOL(40), WGTMSH(40), WGTERR(40), TOLIN(40),
     1                ROOT(40), JTOL(40), LTOL(40), NTOL
C
      EXTERNAL FSUB, DFSUB, GSUB, DGSUB, GUESS
C
C...  constants for control of nonlinear iteration
C
      RELMIN = 1.D-3
      RSTART = 1.D-2
      LMTFRZ = 4
C
C...  compute the maximum tolerance
C
      CHECK = 0.D0
      DO 10 I = 1, NTOL
   10   CHECK = DMAX1 ( TOLIN(I), CHECK )
      IMESH = 1
      ICONV = 0
      IF ( NONLIN .EQ. 0 ) ICONV = 1
      ICOR = 0
      NOCONV = 0
      MSING = 0
      ISING = 0
C
C...  the main iteration begins here .
C...  loop 20 is executed until error tolerances are satisfied or
C...  the code fails (due to a singular matrix or storage limitations)
C
   20      CONTINUE
C
C...       initialization for a new mesh
C
       ITER = 0
       IF ( NONLIN .GT. 0 )                     GO TO 50
C
C...       the linear case.
C...       set up and solve equations
C
       CALL LSYSLV_DAE (MSING, XI, XIOLD, DUMMY, DUMMY, Z, DMZ, G,
     1          W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 0,
     2          FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING,RPAR,IPAR)
C
C...       check for a singular matrix
C
       IF (ISING .NE. 0) THEN
         IF ( IPRINT .LT. 1 )  THEN
          CALL Rprint('Singular projection matrix due to index > 2')
         ENDIF
           IFLAG = 0
         RETURN
       END IF
       IF ( MSING .EQ. 0 )                      GO TO 400
   30      IF ( MSING .LT. 0 )                      GO TO 40
       IF ( IPRINT .LT. 1 )  THEN
        CALL Rprint('A local elimination matrix is singular')
       ENDIF 
          GO TO 460
   40      IF ( IPRINT .LT. 1 ) THEN
            CALL Rprint('The global BVP-matrix is singular')
           ENDIF 
       IFLAG = 0
       RETURN
C
C...       iteration loop for nonlinear case
C...       define the initial relaxation parameter (= relax)
C
   50      RELAX = 1.D0
C
C...       check for previous convergence and problem sensitivity
C
       IF ( ICARE .EQ. (-1) )  RELAX = RSTART
           IF ( ICARE .EQ. 1 )     RELAX = 1.D0
       IF ( ICONV .EQ. 0 )                      GO TO 160
C
C...       convergence on a previous mesh has been obtained.    thus
C...       we have a very good initial approximation for the newton
C...       process.    proceed with one full newton and then iterate
C...       with a fixed jacobian.
C
       IFREEZ = 0
C
C...       evaluate right hand side and its norm  and
C...       find the first newton correction
C
       CALL LSYSLV_DAE (MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
     1          W, V, FC, RHS, DQDMZ, INTEGS, IPVTG, IPVTW, RNOLD, 1,
     2          FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING,RPAR,IPAR)
C
       IF ( IPRINT .LT. 0 )  THEN
        CALL Rprinti1('Fixed Jacobian iteration ', ITER)
        CALL Rprintd1('Norm (RHS) = ', RNOLD)
       ENDIF 
          
       GO TO 70
C
C...       solve for the next iterate .
C...       the value of ifreez determines whether this is a full
C...       newton step (=0) or a fixed jacobian iteration (=1).
C
   60      IF ( IPRINT .LT. 0 )  THEN
        CALL Rprinti1('Fixed Jacobian iteration ', ITER)
        CALL Rprintd1('Norm (RHS) = ', RNOLD)
           ENDIF   
       RNOLD = RNORM
       CALL LSYSLV_DAE (MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
     1          W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM,
     2    3+IFREEZ, FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING,RPAR,IPAR)
C
C...       check for a singular matrix
C
   70      IF ( MSING .NE. 0 )                      GO TO 30
       IF (ISING .NE. 0) THEN
         IF ( IPRINT .LT. 1 ) THEN
        CALL Rprint('Singular projection matrix due to index > 2')
         ENDIF 
           
         IFLAG = 0
         RETURN
       END IF
       IF ( IFREEZ .EQ. 1 )                     GO TO 80
C
C...       a full newton step
C
       ITER = ITER + 1
       IFRZ = 0
   80      CONTINUE
C
C...       update   z and dmz , compute new  rhs  and its norm
C
       DO 90 I = 1, NZ
         Z(I) = Z(I) + DELZ(I)
   90      CONTINUE
       DO 100 I = 1, NDMZ
         DMZ(I) = DMZ(I) + DELDMZ(I)
  100      CONTINUE
       CALL LSYSLV_DAE (MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
     1          W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 2,
     2          FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING,RPAR,IPAR)
C
C...       check monotonicity. if the norm of  rhs  gets smaller,
C...       proceed with a fixed jacobian; else proceed cautiously,
C...       as if convergence has not been obtained before (iconv=0).
C
       IF ( RNORM .LT. PRECIS )                 GO TO 390
       IF ( RNORM .GT. RNOLD )                  GO TO 130
       IF ( IFREEZ .EQ. 1 )                     GO TO 110
       IFREEZ = 1
       GO TO 60
C
C...       verify that the linear convergence with fixed jacobian
C...       is fast enough.
C
  110      IFRZ = IFRZ + 1
       IF ( IFRZ .GE. LMTFRZ )       IFREEZ = 0
       IF ( RNOLD .LT. 4.D0*RNORM )  IFREEZ = 0
C
C...       check convergence (iconv = 1).
C
       DO 120 IT = 1, NTOL
         INZ = LTOL(IT)
         DO 120 IZ = INZ, NZ, MSTAR
           IF ( DABS(DELZ(IZ)) .GT.
     1           TOLIN(IT) * (DABS(Z(IZ)) + 1.D0))  GO TO 60
  120      CONTINUE
C
C...       convergence obtained
C
       IF ( IPRINT .LT. 1 ) THEN
        CALL Rprinti1('Convergence after iteration: ', ITER)
       ENDIF 
       GO TO 400
C
C...      convergence of fixed jacobian iteration failed.
C
  130      IF ( IPRINT .LT. 0 ) THEN
        CALL Rprinti1('Fixed Jacobian iteration ', ITER)
        CALL Rprintd1('Norm (RHS) = ', RNOLD)
        CALL Rprint('Switch to damped Newton iteration')
           ENDIF 
              
       ICONV = 0
       IF ( ICARE .NE. 1 )   RELAX = RSTART
       DO 140 I = 1, NZ
         Z(I) = Z(I) - DELZ(I)
  140      CONTINUE
       DO 150 I = 1, NDMZ
         DMZ(I) = DMZ(I) - DELDMZ(I)
  150      CONTINUE
C
C...       update old mesh
C
       NP1 = N + 1
       DO 155 I = 1, NP1
  155        XIOLD(I) = XI(I)
       NOLD = N
C
       ITER = 0
C
C...       no previous convergence has been obtained. proceed
C...       with the damped newton method.
C...       evaluate rhs and find the first newton correction.
C
  160      IF(IPRINT .LT. 0)  THEN
            CALL Rprint('Full damped Newton iteration')
           ENDIF 

             CALL LSYSLV_DAE (MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
     1          W, V, FC, RHS, DQDMZ, INTEGS, IPVTG, IPVTW, RNOLD, 1,
     2          FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING,RPAR,IPAR)
C
C...       check for a singular matrix
C
       IF ( MSING .NE. 0 )                       GO TO 30
       IF (ISING .NE. 0) THEN
         IF ( IPRINT .LT. 1 ) THEN
          CALL Rprint('Singular projection matrix due to index > 2')
         ENDIF 
         IFLAG = 0
         RETURN
       END IF
C
C...       bookkeeping for first mesh
C
       IF ( IGUESS .EQ. 1 )  IGUESS = 0
C
C...       find initial scaling
C
       CALL SKALE_DAE (N, MSTAR, KDY, Z, DMZ, XI, SCALE, DSCALE)
           RLXOLD = RELAX
       IPRED = 1
       GO TO 220
C
C...       main iteration loop
C
  170      RNOLD = RNORM
       IF ( ITER .GE. LIMIT )                   GO TO 430
C
C...       update scaling
C
       CALL SKALE_DAE (N, MSTAR, KDY, Z, DMZ, XI, SCALE, DSCALE)
C
C...       compute norm of newton correction with new scaling
C
       ANSCL = 0.D0
       DO 180 I = 1, NZ
         ANSCL = ANSCL + (DELZ(I) * SCALE(I))**2
 180   CONTINUE
       DO 190 I = 1, NDMZ
         ANSCL = ANSCL + (DELDMZ(I) * DSCALE(I))**2
  190      CONTINUE
       ANSCL = DSQRT(ANSCL / FLOAT(NZ+NDMZ))
C
C...       find a newton direction
C
       CALL LSYSLV_DAE (MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ, G,
     1          W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 3,
     2          FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING,RPAR,IPAR)
C
C...       check for a singular matrix
C
       IF ( MSING .NE. 0 )                      GO TO 30
       IF (ISING .NE. 0) THEN
         IF ( IPRINT .LT. 1 ) THEN
          CALL Rprint('Singular projection matrix due to index > 2' )
         ENDIF 
         IFLAG = 0
         RETURN
       END IF
C
C...       predict relaxation factor for newton step.
C
           IF (ICARE .NE. 1)  THEN
       ANDIF = 0.D0
       DO 200 I = 1, NZ
         ANDIF = ANDIF + ((DQZ(I) - DELZ(I)) * SCALE(I))**2
  200      CONTINUE
       DO 210 I = 1, NDMZ
         ANDIF = ANDIF + ((DQDMZ(I) - DELDMZ(I)) * DSCALE(I))**2
  210      CONTINUE
       ANDIF = DSQRT(ANDIF/FLOAT(NZ+NDMZ) + PRECIS)
       RELAX = RELAX * ANSCL / ANDIF
       IF ( RELAX .GT. 1.D0 )  RELAX = 1.D0
           RLXOLD = RELAX
       IPRED = 1
           END IF
  220      ITER = ITER + 1
C
C...       determine a new  z and dmz  and find new  rhs  and its norm
C
       DO 230 I = 1, NZ
         Z(I) = Z(I)  +  RELAX * DELZ(I)
  230      CONTINUE
       DO 240 I = 1, NDMZ
         DMZ(I) = DMZ(I)  +  RELAX * DELDMZ(I)
  240      CONTINUE
  250      CALL LSYSLV_DAE (MSING, XI, XIOLD, Z, DMZ, DQZ, DQDMZ, G,
     1          W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 2,
     2          FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING,RPAR,IPAR)
C
C...       compute a fixed jacobian iterate (used to control relax)
C
       CALL LSYSLV_DAE (MSING, XI, XIOLD, Z, DMZ, DQZ, DQDMZ, G,
     1          W, V, FC, RHS, DUMMY, INTEGS, IPVTG, IPVTW, RNORM, 4,
     2          FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING,RPAR,IPAR)
C
C...       find scaled norms of various terms used to correct relax
C
       ANORM = 0.D0
       ANFIX = 0.D0
       DO 260 I = 1, NZ
         ANORM = ANORM  +  (DELZ(I) * SCALE(I))**2
         ANFIX = ANFIX  +  (DQZ(I) * SCALE(I))**2
  260      CONTINUE
       DO 270 I = 1, NDMZ
         ANORM = ANORM  +  (DELDMZ(I) * DSCALE(I))**2
         ANFIX = ANFIX  +  (DQDMZ(I) * DSCALE(I))**2
  270      CONTINUE
       ANORM = DSQRT(ANORM / FLOAT(NZ+NDMZ))
       ANFIX = DSQRT(ANFIX / FLOAT(NZ+NDMZ))
       IF ( ICOR .EQ. 1 )                         GO TO 280
       IF (IPRINT .LT. 0) THEN
       CALL Rprintid('Iteration = , Relaxation factor = ',ITER, RELAX)
       CALL Rprintd2('Norm of scaled RHS changes from, to',ANORM,ANFIX)
       CALL Rprintd2('Norm OF RHS changes from, to', RNOLD, RNORM)
       ENDIF 
       GO TO 290
  280      IF (IPRINT .LT. 0) THEN
       CALL Rprintd1('Relaxation factor corrected to  ', RELAX)
       CALL Rprintd2('Norm of scaled RHS changes from, to',ANORM,ANFIX)
       CALL Rprintd2('Norm or RHS changes from, to', RNOLD, RNORM)
           ENDIF 

  290      ICOR = 0
C
C...       check for monotonic decrease in  delz and deldmz.
C
       IF (ANFIX.LT.PRECIS .OR. RNORM.LT.PRECIS)  GO TO 390
       IF ( ANFIX .GT. ANORM .AND. ICARE .NE. 1)  GO TO 300
C
C...       we have a decrease.
C...       if  dqz  and dqdmz  small, check for convergence
C
       IF ( ANFIX .LE. CHECK )                    GO TO 350
C
C...       correct the predicted  relax  unless the corrected
C...       value is within 10 percent of the predicted one.
C
       IF ( IPRED .NE. 1 )                        GO TO 170
  300      IF ( ITER .GE. LIMIT )                     GO TO 430
           IF ( ICARE .EQ. 1 )                        GO TO 170
C
C...       correct the relaxation factor.
C
       IPRED = 0
       ARG = (ANFIX/ANORM - 1.D0) / RELAX + 1.D0
       IF ( ARG .LT. 0.D0 )                       GO TO 170
       IF (ARG .LE. .25D0*RELAX+.125D0*RELAX**2)  GO TO 310
       FACTOR = -1.D0 + DSQRT (1.D0+8.D0 * ARG)
       IF ( DABS(FACTOR-1.D0) .LT. .1D0*FACTOR )  GO TO 170
       IF ( FACTOR .LT. 0.5D0 )  FACTOR = 0.5D0
       RELAX = RELAX / FACTOR
       GO TO 320
  310      IF ( RELAX .GE. .9D0 )                     GO TO 170
       RELAX = 1.D0
  320      ICOR = 1
       IF ( RELAX .LT. RELMIN )                   GO TO 440
       FACT = RELAX - RLXOLD
       DO 330 I = 1, NZ
        Z(I) = Z(I)  +  FACT * DELZ(I)
  330      CONTINUE
       DO 340 I = 1, NDMZ
         DMZ(I) = DMZ(I)  +  FACT * DELDMZ(I)
  340      CONTINUE
       RLXOLD = RELAX
       GO TO 250
C
C...       check convergence (iconv = 0).
C
  350      CONTINUE
       DO 360 IT = 1, NTOL
         INZ = LTOL(IT)
         DO 360 IZ = INZ, NZ, MSTAR
           IF ( DABS(DQZ(IZ)) .GT.
     1         TOLIN(IT) * (DABS(Z(IZ)) + 1.D0) )   GO TO 170
  360      CONTINUE
C
C...       convergence obtained
C
       IF ( IPRINT .LT. 1 ) THEN 
       CALL Rprinti1('Convergence after iteration ',ITER)
       ENDIF 
C
C...       since convergence obtained, update  z and dmz  with term
C...       from the fixed jacobian iteration.
C
       DO 370 I = 1, NZ
         Z(I) = Z(I)  +  DQZ(I)
  370      CONTINUE
       DO 380 I = 1, NDMZ
         DMZ(I) = DMZ(I)  +  DQDMZ(I)
  380      CONTINUE
  390      IF ( (ANFIX .LT. PRECIS .OR. RNORM .LT. PRECIS)
     1          .AND. IPRINT .LT. 1 ) THEN
       CALL Rprinti1('Convergence after iteration ',ITER)
           ENDIF 
       ICONV = 1
       IF ( ICARE .EQ. (-1) )  ICARE = 0
C
C...       if full output has been requested, print values of the
C...       solution components   z  at the meshpoints and  y  at
C...       collocation points.
C...       check for error tolerance satisfaction
C
  400      IFIN = 1
       IF (IMESH .EQ. 2) CALL ERRCHK_DAE (XI, Z, DMZ, VALSTR, IFIN)
       IF ( IMESH .EQ. 1 .OR.
     1          IFIN .EQ. 0 .AND. ICARE .NE. 2)     GO TO 460
       IFLAG = 1
       RETURN
C
C...       diagnostics for failure of nonlinear iteration.
C
  430      IF ( IPRINT .LT. 1 ) THEN
       CALL Rprinti1('NO convergence after iteration ',ITER)
           ENDIF 
       GO TO 450
  440       IF( IPRINT .LT. 1 )    THEN
      CALL Rprintd1('NO convergence, relaxation factor too small',RELAX)
      CALL Rprintd1('Should not be less than', RELMIN)
       ENDIF
  450      IFLAG = -2
       NOCONV = NOCONV + 1
       IF ( ICARE .EQ. 2 .AND. NOCONV .GT. 1 )  RETURN
       IF ( ICARE .EQ. 0 )  ICARE = -1
C
C...       update old mesh
C
  460      NP1 = N + 1
       DO 470 I = 1, NP1
  470        XIOLD(I) = XI(I)
       NOLD = N
C
C...       pick a new mesh
C...       check safeguards for mesh refinement
C
       IMESH = 1
       IF ( ICONV .EQ. 0 .OR. MSHNUM .GE. MSHLMT
     1                       .OR. MSHALT .GE. MSHLMT )  IMESH = 2
       IF ( MSHALT .GE. MSHLMT .AND.
     1        MSHNUM .LT. MSHLMT )  MSHALT = 1
       IF (NY.EQ.0) THEN
           NYCB = 1
       ELSE
           NYCB = NY
       ENDIF

       CALL NEWMSH_DAE (IMESH, XI, XIOLD, Z, DMZ, DMV, VALSTR,
     1                  SLOPE, ACCUM, NFXPNT, FIXPNT, DF, DFSUB,
     2            FCSP, CBSP, NCOMP, NYCB,RPAR,IPAR)
C
C...       exit if expected n is too large (but may try n=nmax once)
C
       IF ( N .LE. NMAX )                       GO TO 480
       N = N / 2
       IFLAG = -1
       IF ( ICONV .EQ. 0 .AND. IPRINT .LT. 1 ) THEN
        CALL Rprint('NO convergence')
       ENDIF 
       IF ( ICONV .EQ. 1 .AND. IPRINT .LT. 1 ) THEN 
       CALL Rprint('Probably tolerance too stringent or nmax too small')
       ENDIF 
       RETURN
  480      IF ( ICONV .EQ. 0 )  IMESH = 1
C          IF ( ICARE .EQ. 1 )  ICONV = 0
      GO TO 20
C     ---------------------------------------------------------------
      END

      SUBROUTINE SKALE_DAE(N, MSTAR, KDY, Z, DMZ, XI, SCALE, DSCALE)
C
C**********************************************************************
C
C   purpose
C            provide a proper scaling of the state variables, used
C            to control the damping factor for a newton iteration [4].
C
C   variables
C
C            n      = number of mesh subintervals
C            mstar  = number of unknomns in z(u(x))
C            kdy     = number of unknowns in dmz per mesh subinterval
C            z      = the global current solution vector
C            dmz    = the global current highest derivs vector
C            xi     = the current mesh
C            scale  = scaling vector for z
C            dscale = scaling vector for dmz
C
C**********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(MSTAR,*), SCALE(MSTAR,*), DMZ(KDY,N), DSCALE(KDY,N)
      DIMENSION XI(*), BASM(5)
C
      COMMON /DAEORD/ K, NCOMP, NY, NCY, ID1, KD, ID3, MMAX, M(20)
C
      BASM(1) = 1.D0
      DO 50 J=1,N
      IZ = 1
      H = XI(J+1) - XI(J)
      DO 10 L = 1, MMAX
      BASM(L+1) = BASM(L) * H / FLOAT(L)
  10    CONTINUE
      DO 40 ICOMP = 1, NCOMP
      SCAL = (DABS(Z(IZ,J)) + DABS(Z(IZ,J+1))) * .5D0 + 1.D0
      MJ = M(ICOMP)
      DO 20 L = 1, MJ
        SCALE(IZ,J) = BASM(L) / SCAL
        IZ = IZ + 1
  20      CONTINUE
      SCAL = BASM(MJ+1) / SCAL
      DO 30 IDMZ = ICOMP, KDY, NCY
        DSCALE(IDMZ,J) = SCAL
  30      CONTINUE
  40    CONTINUE
      DO 45 ICOMP = 1+NCOMP,NCY
      SCAL = 1.D0 / (DABS(DMZ(ICOMP,J)) + 1.D0)
      DO 45 IDMZ = ICOMP, KDY, NCY
        DSCALE(IDMZ,J) = SCAL
  45    CONTINUE
  50  CONTINUE
      NP1 = N + 1
      DO 60 IZ = 1, MSTAR
      SCALE(IZ,NP1) = SCALE(IZ,N)
  60  CONTINUE
      RETURN
      END
C----------------------------------------------------------------------
C                            p a r t  2
C          mesh selection, error estimation, (and related
C          constant assignment) routines -- see [5], [6]
C----------------------------------------------------------------------
C
      SUBROUTINE NEWMSH_DAE(MODE, XI, XIOLD, Z, DMZ, DMV, VALSTR,
     1                   SLOPE, ACCUM, NFXPNT, FIXPNT, DF, DFSUB,
     2             FC, CB, NCOMP, NYCB,RPAR,IPAR)
C
C**********************************************************************
C
C   purpose
C            select a mesh on which a collocation solution is to be
C            determined
C
C                           there are 5 possible modes of action:
C            mode = 5,4,3 - deal mainly with definition of an initial
C                           mesh for the current boundary value problem
C                 = 2,1   - deal with definition of a new mesh, either
C                           by simple mesh halving or by mesh selection
C            more specifically, for
C            mode = 5  an initial (generally nonuniform) mesh is
C                      defined by the user and no mesh selection is to
C                      be performed
C                 = 4  an initial (generally nonuniform) mesh is
C                      defined by the user
C                 = 3  a simple uniform mesh (except possibly for some
C                      fixed points) is defined; n= no. of subintervals
C                 = 1  the automatic mesh selection procedure is used
C                      (see [5] for details)
C                 = 2  a simple mesh halving is performed
C
C**********************************************************************
C
C   variables
C
C            n      = number of mesh subintervals
C            nold   = number of subintervals for former mesh
C            xi     - mesh point array
C            xiold  - former mesh point array
C            mshlmt - maximum no. of mesh selections which are permitted
C                     for a given n before mesh halving
C            mshnum - no. of mesh selections which have actually been
C                     performed for the given n
C            mshalt - no. of consecutive times ( plus 1 ) the mesh
C                     selection has alternately halved and doubled n.
C                     if mshalt .ge. mshlmt then  contrl  requires
C                     that the current mesh be halved.
C            mshflg = 1  the mesh is a halving of its former mesh
C                       (so an error estimate has been calculated)
C                   = 0  otherwise
C            iguess - ISET(9) in subroutine coldae.  it is used
C                     here only for mode=5 and 4, where
C                   = 2 the subroutine sets xi=xiold.  this is
C                       used e.g. if continuation is being per-
C                       formed, and a mesh for the old differen-
C                       tial equation is being used
C                   = 3 same as for =2, except xi uses every other
C                       point of xiold (so mesh xiold is mesh xi
C                       halved)
C                   = 4 xi has been defined by the user, and an old
C                       mesh xiold is also available
C                       otherwise, xi has been defined by the user
C                       and we set xiold=xi in this subroutine
C            slope  - an approximate quantity to be equidistributed for
C                     mesh selection (see [5]), viz,
C                             .                        (k+mj)
C                     slope(i)=     max   (weight(l) *u      (xi(i)))
C                               1.le.l.le.ntol         j
C
C                     where j=jtol(l)
C            slphmx - maximum of slope(i)*(xiold(i+1)-xiold(i)) for
C                     i = 1 ,..., nold.
C            accum  - accum(i) is the integral of  slope  from  aleft
C                     to  xiold(i).
C            valstr - is assigned values needed in  errchk  for the
C                     error estimate.
C            fc     - you know
C**********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D1(40), D2(40), SLOPE(*), ACCUM(*), VALSTR(*), DMV(*)
      DIMENSION XI(*), XIOLD(*), Z(*), DMZ(*), FIXPNT(*), DUMMY(1)
      DIMENSION FC(NCOMP,60), ZVAL(40), YVAL(40), A(28), DF(NCY,*)
      DIMENSION CB(NYCB,NYCB), IPVTCB(40), BCOL(40), U(400), V(400)
      DIMENSION RPAR(*),IPAR(*)
      EXTERNAL DFSUB
C
      COMMON /COLLOC/ RHO(7), COEF(49)
      COMMON /DAEOUT/ PRECIS, IOUT, IPRINT
      COMMON /DAEORD/ K, NCDUM, NY, NCY, MSTAR, KD, KDY, MMAX, M(20)
      COMMON /DAEAPR/ N, NOLD, NMAX, NZ, NDMZ
      COMMON /DAEMSH/ MSHFLG, MSHNUM, MSHLMT, MSHALT
      COMMON /DAENLN/ NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX
      COMMON /DAESID/  ZETA(40), ALEFT, ARIGHT, IZETA, IDUM
      COMMON /DAEBAS/ B(28), ACOL(28,7), ASAVE(28,4)
      COMMON /DAEEST/ TOL(40), WGTMSH(40), WGTERR(40), TOLIN(40),
     1                ROOT(40), JTOL(40), LTOL(40), NTOL

      integer nfunc, njac, nstep, nbound, njacbound
      common/CDAEdiag/nfunc, njac, nstep, nbound, njacbound
C
      NFXP1 = NFXPNT +1
      GO TO (180, 100, 50, 20, 10), MODE
C
C...  mode=5   set mshlmt=1 so that no mesh selection is performed
C
   10 MSHLMT = 1
C
C...  mode=4   the user-specified initial mesh is already in place.
C
   20 IF ( IGUESS .LT. 2 )                          GO TO 40
C
C...  iguess=2, 3 or 4.
C
      NOLDP1 = NOLD + 1
      IF ( IGUESS .NE. 3 )                          GO TO 40
C
C...  if iread ( ISET(8) ) .ge. 1 and iguess ( ISET(9) ) .eq. 3
C...  then the first mesh is every second point of the
C...  mesh in  xiold .
C
      N = NOLD /2
      I = 0
      DO 30 J = 1, NOLD, 2
       I = I + 1
   30 XI(I) = XIOLD(J)
   40 CONTINUE
      NP1 = N + 1
      XI(1) = ALEFT
      XI(NP1) = ARIGHT
      GO TO 320
C
C...  mode=3   generate a (piecewise) uniform mesh. if there are
C...  fixed points then ensure that the n being used is large enough.
C
   50 IF ( N .LT. NFXP1 )  N = NFXP1
      NP1 = N + 1
      XI(1) = ALEFT
      ILEFT = 1
      XLEFT = ALEFT
C
C...  loop over the subregions between fixed points.
C
      DO 90 J = 1, NFXP1
       XRIGHT = ARIGHT
       IRIGHT = NP1
       IF ( J .EQ. NFXP1 )                      GO TO 60
       XRIGHT = FIXPNT(J)
C
C...       determine where the j-th fixed point should fall in the
C...       new mesh - this is xi(iright) and the (j-1)st fixed
C...       point is in xi(ileft)
C
       NMIN = (XRIGHT-ALEFT) / (ARIGHT-ALEFT) * FLOAT(N) + 1.5D0
       IF (NMIN .GT. N-NFXPNT+J)  NMIN = N - NFXPNT + J
       IRIGHT = MAX0 (ILEFT+1, NMIN)
   60      XI(IRIGHT) = XRIGHT
C
C...       generate equally spaced points between the j-1st and the
C...       j-th fixed points.
C
       NREGN = IRIGHT - ILEFT - 1
       IF ( NREGN .EQ. 0 )                      GO TO 80
       DX = (XRIGHT - XLEFT) / FLOAT(NREGN+1)
       DO 70 I = 1, NREGN
   70      XI(ILEFT+I) = XLEFT  +  FLOAT(I) * DX
   80      ILEFT = IRIGHT
       XLEFT = XRIGHT
   90 CONTINUE
      GO TO 320
C
C...  mode=2  halve the current mesh (i.e. double its size)
C
  100 N2 = 2 * N
C
C...  check that n does not exceed storage limitations
C
      IF ( N2 .LE. NMAX )                           GO TO 120
C
C...  if possible, try with n=nmax. redistribute first.
C
      IF ( MODE .EQ. 2 )                            GO TO 110
      N = NMAX / 2
      GO TO 220
  110 IF ( IPRINT .LT. 1 ) THEN 
        CALL Rprint('Expected N too large')
       ENDIF 
      N = N2
      RETURN
C
C...  calculate the old approximate solution values at
C...  points to be used in  errchk  for error estimates.
C...  if  mshflg  =1 an error estimate was obtained for
C...  for the old approximation so half the needed values
C...  will already be in  valstr .
C
  120 IF ( MSHFLG .EQ. 0 )                          GO TO 140
C
C...  save in  valstr  the values of the old solution
C...  at the relative positions 1/6 and 5/6 in each subinterval.
C
      KSTORE = 1
      DO 130 I = 1, NOLD
      HD6 = (XIOLD(I+1) - XIOLD(I)) / 6.D0
      X = XIOLD(I) + HD6
      CALL APPROX_DAE (I, X, VALSTR(KSTORE), DUMMY, ASAVE(1,1),
     +         DUMMY, XIOLD, NOLD, Z, DMZ,
     1         K, NCOMP, NY, MMAX, M, MSTAR, 4, DUMMY, 0)
      X = X + 4.D0 * HD6
      KSTORE = KSTORE  +  3 * MSTAR
      CALL APPROX_DAE (I, X, VALSTR(KSTORE), DUMMY, ASAVE(1,4),
     +         DUMMY, XIOLD, NOLD, Z, DMZ,
     1         K, NCOMP, NY, MMAX, M, MSTAR, 4, DUMMY, 0)
      KSTORE = KSTORE  +  MSTAR
  130 CONTINUE
      GO TO 160
C
C...  save in  valstr  the values of the old solution
C...  at the relative positions 1/6, 2/6, 4/6 and 5/6 in
C...  each subinterval.
C
  140 KSTORE = 1
      DO 150 I = 1, N
       X = XI(I)
       HD6 = (XI(I+1) - XI(I)) / 6.D0 
       DO 150 J = 1, 4
       X = X + HD6
       IF ( J.EQ.3 )  X = X + HD6
       CALL APPROX_DAE (I, X, VALSTR(KSTORE), DUMMY, ASAVE(1,J),
     +          DUMMY, XIOLD, NOLD, Z, DMZ,
     1          K, NCOMP, NY, MMAX, M, MSTAR, 4, DUMMY, 0)
       KSTORE = KSTORE  +  MSTAR
  150 CONTINUE
  160 MSHFLG = 0
      MSHNUM = 1
      MODE = 2
C
C...  generate the halved mesh.
C
      J = 2
      DO 170 I = 1, N
       XI(J) = (XIOLD(I) + XIOLD(I+1)) / 2.D0
       XI(J+1) = XIOLD(I+1)
  170 J = J + 2
      N = N2
      GO TO 320
C
C...  mode=1  we do mesh selection if it is deemed worthwhile
C
  180 IF ( NOLD .EQ. 1 )                            GO TO 100
      IF ( NOLD .LE. 2*NFXPNT )                     GO TO 100
C
C...  we now project DMZ for mesh selection strategy, if required
C...  but set DMV = DMZ in case it is not

      IDMZ = 1
      DO 183 I = 1, NOLD
       DO 182 KK = 1, K
        DO 181 J = 1, NCY
           DMV(IDMZ) = DMZ(IDMZ)
           IDMZ = IDMZ + 1
181         CONTINUE
182      CONTINUE
183   CONTINUE

      IF (INDEX.NE.1 .AND. NY .GT. 0) THEN
       IDMZ = 1
       DO 500 I=1, NOLD
          XI1 = XIOLD(I+1)
          CALL APPROX_DAE (I, XI1, ZVAL, YVAL, A,
     1                     COEF, XIOLD, NOLD, Z, DMZ,
     1                     K, NCOMP, NY, MMAX, M, MSTAR, 3, DUMMY, 1)
         CALL DFSUB (NCY, XI1, ZVAL, YVAL, DF,RPAR,IPAR)
c francesca: added
         njac = njac + 1

C...         if index=2, form projection matrices directly
C...         otherwise use svd to define appropriate projection

         IF (INDEX.EQ.0) THEN
           CALL PRJSVD (FC,DF,CB,U,V,NCOMP,NCY,NY,IPVTCB,ISING,2)
         ELSE

C...        form cb

           DO 212 J = 1, NY
          DO 212 J1 = 1, NY
             FACT = 0.0D0
             ML = 0
             DO 211 L = 1, NCOMP
            ML = ML + M(L)
            FACT = FACT + DF(J+NCOMP,ML)*DF(L,MSTAR+J1)
211                  CONTINUE
             CB(J,J1) = FACT
212            CONTINUE

C...           decompose cb

           CALL DGEFA (CB, NY, NY, IPVTCB, ISING)
           IF (ISING.NE.0) RETURN

C...           form columns of fc

           ML = 0
           DO 215 L = 1, NCOMP
          ML = ML + M(L)
          DO 213 J1 = 1, NY
             BCOL(J1) = DF(J1+NCOMP,ML)
213               CONTINUE

          CALL DGESL(CB, NY, NY, IPVTCB, BCOL, 0)

          DO 215 J1 = 1, NCOMP
             FACT = 0.0D0
             DO 214 J = 1, NY
            FACT = FACT + DF(J1, J+MSTAR)*BCOL(J)
214                  CONTINUE
             FC(J1,L) = FACT
C                 CONTINUE
215            CONTINUE

         ENDIF

C..         finally, replace fc with the true projection SR = I - fc

           DO 217 J = 1, NCOMP
          DO 216 L = 1, NCOMP
             FC(J,L) = -FC(J,L)
             IF (J.EQ.L) FC(J,L) = FC(J,L) + 1.0D0
216               CONTINUE
217            CONTINUE

C...        project DMZ for the k collocation points, store in DMV

        DO 221 KK = 1, K

           DO 219 J = 1, NCOMP
          FACT = 0.0D0
          DO 218 L = 1, NCOMP
             FACT = FACT + FC(J,L)*DMZ(IDMZ+L-1)
218               CONTINUE
          DMV(IDMZ+J-1) = FACT
219            CONTINUE

           IDMZ = IDMZ + NCY
221         CONTINUE

500      CONTINUE

      ENDIF
C
C...  the first interval has to be treated separately from the
C...  other intervals (generally the solution on the (i-1)st and ith
C...  intervals will be used to approximate the needed derivative, but
C...  here the 1st and second intervals are used.)
C
      I = 1
      HIOLD = XIOLD(2) - XIOLD(1)
      CALL HORDER_DAE (1, D1, HIOLD, DMV, NCOMP, NCY, K)
      IDMZ = IDMZ + (NCOMP + NY) * K
      HIOLD = XIOLD(3) - XIOLD(2)
      CALL HORDER_DAE (2, D2, HIOLD, DMV, NCOMP, NCY, K)
      ACCUM(1) = 0.D0
      SLOPE(1) = 0.D0
      ONEOVH = 2.D0 / ( XIOLD(3) - XIOLD(1) )
      DO 190 J = 1, NTOL
       JJ = JTOL(J)
       JZ = LTOL(J)
  190 SLOPE(1) = DMAX1(SLOPE(1),(DABS(D2(JJ)-D1(JJ))*WGTMSH(J)*
     1           ONEOVH / (1.D0 + DABS(Z(JZ)))) **ROOT(J))
      SLPHMX = SLOPE(1) * (XIOLD(2) - XIOLD(1))
      ACCUM(2) = SLPHMX
      IFLIP = 1
C
C...  go through the remaining intervals generating  slope
C...  and  accum .
C
      DO 210 I = 2, NOLD
       HIOLD = XIOLD(I+1) - XIOLD(I)
       IF ( IFLIP .EQ. -1 )
     1        CALL HORDER_DAE ( I, D1, HIOLD, DMV, NCOMP, NCY, K)
       IF ( IFLIP .EQ. 1 )
     1        CALL HORDER_DAE ( I, D2, HIOLD, DMV, NCOMP, NCY, K)
       ONEOVH = 2.D0 / ( XIOLD(I+1) - XIOLD(I-1) )
       SLOPE(I) = 0.D0
C
C...       evaluate function to be equidistributed
C
       DO 200 J = 1, NTOL
         JJ = JTOL(J)
         JZ = LTOL(J)  +  (I-1) * MSTAR
         SLOPE(I) = DMAX1(SLOPE(I),(DABS(D2(JJ)-D1(JJ))*WGTMSH(J)*
     1                  ONEOVH / (1.D0 + DABS(Z(JZ)))) **ROOT(J))
  200      CONTINUE
C
C...       accumulate approximate integral of function to be
C...       equidistributed
C
       TEMP = SLOPE(I) * (XIOLD(I+1)-XIOLD(I))
       SLPHMX = DMAX1(SLPHMX,TEMP)
       ACCUM(I+1) = ACCUM(I) + TEMP
       IFLIP = - IFLIP
  210 CONTINUE
      AVRG = ACCUM(NOLD+1) / FLOAT(NOLD)
      DEGEQU = AVRG / DMAX1(SLPHMX,PRECIS)
C
C...  naccum=expected n to achieve .1x user requested tolerances
C
      NACCUM = ACCUM(NOLD+1) + 1.D0
      IF ( IPRINT .LT. 0 ) THEN 
      CALL Rprintd1('Mesh degree of equidistribution =',DEGEQU)
      CALL Rprinti1('Prediction for required N = ', NACCUM)
       ENDIF 

C
C...  decide if mesh selection is worthwhile (otherwise, halve)
C
      IF ( AVRG .LT. PRECIS )                       GO TO 100
      IF ( DEGEQU .GE. .5D0 )                       GO TO 100
C
C...  nmx assures mesh has at least half as many subintervals as the
C...  previous mesh
C
      NMX = MAX0 ( NOLD+1, NACCUM ) / 2
C
C...  this assures that halving will be possible later (for error est)
C
      NMAX2 = NMAX / 2
C
C...  the mesh is at most halved
C
      N = MIN0 ( NMAX2, NOLD, NMX )
  220 NOLDP1 = NOLD + 1
      IF ( N .LT. NFXP1 )  N = NFXP1
      MSHNUM = MSHNUM + 1
C
C...  if the new mesh is smaller than the old mesh set mshnum
C...  so that the next call to  newmsh  will produce a halved
C...  mesh. if n .eq. nold / 2 increment mshalt so there can not
C...  be an infinite loop alternating between n and n/2 points.
C
      IF ( N .LT. NOLD )  MSHNUM = MSHLMT
      IF ( N .GT. NOLD/2 )  MSHALT = 1
      IF ( N .EQ. NOLD/2 )  MSHALT = MSHALT + 1
      MSHFLG = 0
C
C...  having decided to generate a new mesh with n subintervals we now
C...  do so, taking into account that the nfxpnt points in the array
C...  fixpnt must be included in the new mesh.
C
      IN = 1
      ACCL = 0.D0
      LOLD = 2
      XI(1) = ALEFT
      XI(N+1) = ARIGHT
      DO 310 I = 1, NFXP1
       IF ( I .EQ. NFXP1 )                      GO TO 250
       DO 230 J = LOLD, NOLDP1
         LNEW = J
         IF ( FIXPNT(I) .LE. XIOLD(J) )         GO TO 240
  230      CONTINUE
  240      CONTINUE
       ACCR = ACCUM(LNEW) + (FIXPNT(I)-XIOLD(LNEW))*SLOPE(LNEW-1)
       NREGN = (ACCR-ACCL) / ACCUM(NOLDP1) * FLOAT(N) - .5D0
       NREGN = MIN0(NREGN, N - IN - NFXP1 + I)
       XI(IN + NREGN + 1) = FIXPNT(I)
       GO TO 260
  250      ACCR = ACCUM(NOLDP1)
       LNEW = NOLDP1
       NREGN = N - IN
  260      IF ( NREGN .EQ. 0 )                      GO TO 300
       TEMP = ACCL
       TSUM = (ACCR - ACCL) / FLOAT(NREGN+1)
       DO 290 J = 1, NREGN
         IN = IN + 1
         TEMP = TEMP + TSUM
         DO 270 L = LOLD, LNEW
           LCARRY = L
           IF ( TEMP .LE. ACCUM(L) )            GO TO 280
  270        CONTINUE
  280        CONTINUE
         LOLD = LCARRY
  290      XI(IN) = XIOLD(LOLD-1) + (TEMP - ACCUM(LOLD-1)) /
     1     SLOPE(LOLD-1)
  300      IN = IN + 1
       ACCL = ACCR
       LOLD = LNEW
  310 CONTINUE
      MODE = 1
  320 CONTINUE
      NP1 = N + 1

      NZ   = MSTAR * (N + 1)
      NDMZ = KDY * N
      RETURN
C----------------------------------------------------------------
      END

      SUBROUTINE CONSTS_DAE (K, RHO, COEF)
C
C**********************************************************************
C
C   purpose
C            assign (once) values to various array constants.
C
C   arrays assigned during compilation:
C     cnsts1 - weights for extrapolation error estimate
C     cnsts2 - weights for mesh selection
C              (the above weights come from the theoretical form for
C              the collocation error -- see [5])
C
C   arrays assigned during execution:
C     wgterr - the particular values of cnsts1 used for current run
C              (depending on k, m)
C     wgtmsh - gotten from the values of cnsts2 which in turn are
C              the constants in the theoretical expression for the
C              errors. the quantities in wgtmsh are 10x the values
C              in cnsts2 so that the mesh selection algorithm
C              is aiming for errors .1x as large as the user
C              requested tolerances.
C     jtol   - components of differential system to which tolerances
C              refer (viz, if ltol(i) refers to a derivative of u(j),
C              then jtol(i)=j)
C     root   - reciprocals of expected rates of convergence of compo-
C              nents of z(j) for which tolerances are specified
C     rho    - the k collocation points on (0,1)
C     coef   -
C     acol  -  the runge-kutta coefficients values at collocation
C              points
C
C**********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RHO(7), COEF(K,*), CNSTS1(28), CNSTS2(28), DUMMY(1)
C
      COMMON /DAEORD/ KDUM, NCOMP, NY, NCY, MSTAR, KD, KDY, MMAX, M(20)
      COMMON /DAEBAS/ B(28), ACOL(28,7), ASAVE(28,4)
      COMMON /DAEEST/ TOL(40), WGTMSH(40), WGTERR(40), TOLIN(40),
     1                ROOT(40), JTOL(40), LTOL(40), NTOL
C
      DATA CNSTS1 /    .25D0,     .625D-1,  7.2169D-2, 1.8342D-2,
     1     1.9065D-2, 5.8190D-2, 5.4658D-3, 5.3370D-3, 1.8890D-2,
     2     2.7792D-2, 1.6095D-3, 1.4964D-3, 7.5938D-3, 5.7573D-3,
     3     1.8342D-2, 4.673D-3,  4.150D-4,  1.919D-3,  1.468D-3,
     4     6.371D-3,  4.610D-3,  1.342D-4,  1.138D-4,  4.889D-4,
     5     4.177D-4,  1.374D-3,  1.654D-3,  2.863D-3  /
      DATA CNSTS2 /   1.25D-1,   2.604D-3,  8.019D-3,  2.170D-5,
     1     7.453D-5,  5.208D-4,  9.689D-8,  3.689D-7,  3.100D-6,
     2     2.451D-5,  2.691D-10, 1.120D-9,  1.076D-8,  9.405D-8,
     3     1.033D-6,  5.097D-13, 2.290D-12, 2.446D-11, 2.331D-10,
     4     2.936D-9,  3.593D-8,  7.001D-16, 3.363D-15, 3.921D-14,
     5     4.028D-13, 5.646D-12, 7.531D-11, 1.129D-9  /
C
C...  assign weights for error estimate
C
      KOFF = K * ( K + 1 ) / 2
      IZ = 1
      DO 10 J = 1, NCOMP
       MJ = M(J)
       DO 10 L = 1, MJ
         WGTERR(IZ) = CNSTS1(KOFF - MJ + L)
         IZ = IZ + 1
   10 CONTINUE
C
C...  assign array values for mesh selection: wgtmsh, jtol, and root
C
      JCOMP = 1
      MTOT = M(1)
      DO 40 I = 1, NTOL
       LTOLI = LTOL(I)
   20      CONTINUE
       IF ( LTOLI .LE. MTOT )                   GO TO 30
       JCOMP = JCOMP + 1
       MTOT = MTOT + M(JCOMP)
       GO TO 20
   30      CONTINUE
       JTOL(I) = JCOMP
       WGTMSH(I) = 1.D1 * CNSTS2(KOFF+LTOLI-MTOT) / TOLIN(I)
       ROOT(I) = 1.D0 / FLOAT(K+MTOT-LTOLI+1)
   40 CONTINUE
C
C...  specify collocation points
C
      GO TO (50,60,70,80,90,100,110), K
   50 RHO(1) = 0.D0
      GO TO 120
   60 RHO(2) = .57735026918962576451D0
      RHO(1) = - RHO(2)
      GO TO 120
   70 RHO(3) = .77459666924148337704D0
      RHO(2) = .0D0
      RHO(1) = - RHO(3)
      GO TO 120
   80 RHO(4) = .86113631159405257523D0
      RHO(3) = .33998104358485626480D0
      RHO(2) = - RHO(3)
      RHO(1) = - RHO(4)
      GO TO 120
   90 RHO(5) = .90617984593866399280D0
      RHO(4) = .53846931010568309104D0
      RHO(3) = .0D0
      RHO(2) = - RHO(4)
      RHO(1) = - RHO(5)
      GO TO 120
  100 RHO(6) = .93246951420315202781D0
      RHO(5) = .66120938646626451366D0
      RHO(4) = .23861918608319690863D0
      RHO(3) = -RHO(4)
      RHO(2) = -RHO(5)
      RHO(1) = -RHO(6)
      GO TO 120
  110 RHO(7) = .949107991234275852452D0
      RHO(6) = .74153118559939443986D0
      RHO(5) = .40584515137739716690D0
      RHO(4) = 0.D0
      RHO(3) = -RHO(5)
      RHO(2) = -RHO(6)
      RHO(1) = -RHO(7)
  120 CONTINUE
C
C...  map (-1,1) to (0,1) by  t = .5 * (1. + x)
C
      DO 130 J = 1, K
       RHO(J) = .5D0 * (1.D0 + RHO(J))
  130 CONTINUE
C
C...  now find runge-kutta coeffitients b, acol and asave
C...  the values of asave are to be used in  newmsh  and errchk .
C
      DO 140 J = 1, K
       DO 135 I = 1, K
  135      COEF(I,J) = 0.D0
       COEF(J,J) = 1.D0
       CALL VMONDE (RHO, COEF(1,J), K)
  140 CONTINUE
      CALL RKBAS ( 1.D0, COEF, K, MMAX, B, DUMMY, 0)
      DO 150 I = 1, K
       CALL RKBAS ( RHO(I), COEF, K, MMAX, ACOL(1,I), DUMMY, 0)
  150 CONTINUE
      CALL RKBAS ( 1.D0/6.D0, COEF, K, MMAX, ASAVE(1,1), DUMMY, 0)
      CALL RKBAS ( 1.D0/3.D0, COEF, K, MMAX, ASAVE(1,2), DUMMY, 0)
      CALL RKBAS ( 2.D0/3.D0, COEF, K, MMAX, ASAVE(1,3), DUMMY, 0)
      CALL RKBAS ( 5.D0/6.D0, COEF, K, MMAX, ASAVE(1,4), DUMMY, 0)
      RETURN
      END

      SUBROUTINE ERRCHK_DAE (XI, Z, DMZ, VALSTR, IFIN)
C
C**********************************************************************
C
C      purpose
C               determine the error estimates and test to see if the
C               error tolerances are satisfied.
C
C      variables
C        xi     - current mesh points
C        valstr - values of the previous solution which are needed
C                 for the extrapolation- like error estimate.
C        wgterr - weights used in the extrapolation-like error
C                 estimate. the array values are assigned in
C                 subroutine  consts.
C        errest - storage for error estimates
C        err    - temporary storage used for error estimates
C        z      - approximate solution on mesh xi
C        ifin   - a 0-1 variable. on return it indicates whether
C                 the error tolerances were satisfied
C        mshflg - is set by errchk to indicate to newmsh whether
C                 any values of the current solution are stored in
C                 the array valstr. (0 for no, 1 for yes)
C
C**********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ERR(40), ERREST(40), DUMMY(1)
      DIMENSION XI(*), Z(*), DMZ(*), VALSTR(*)
C
      COMMON /DAEOUT/ PRECIS, IOUT, IPRINT
      COMMON /DAEORD/ K, NCOMP, NY, NCY, MSTAR, KD, KDY, MMAX, M(20)
      COMMON /DAEAPR/ N, NOLD, NMAX, NZ, NDMZ
      COMMON /DAEMSH/ MSHFLG, MSHNUM, MSHLMT, MSHALT
      COMMON /DAEBAS/ B(28), ACOL(28,7), ASAVE(28,4)
      COMMON /DAEEST/ TOL(40), WGTMSH(40), WGTERR(40), TOLIN(40),
     1                ROOT(40), JTOL(40), LTOL(40), NTOL
C
C...  error estimates are to be generated and tested
C...  to see if the tolerance requirements are satisfied.
C
      IFIN = 1
      MSHFLG = 1
      DO 10 J = 1, MSTAR
   10   ERREST(J) = 0.D0
      DO 60 IBACK = 1, N
       I = N + 1 - IBACK
C
C...       the error estimates are obtained by combining values of
C...       the numerical solutions for two meshes.
C...       for each value of iback we will consider the two
C...       approximations at 2 points in each of
C...       the new subintervals.  we work backwards through
C...       the subinterval so that new values can be stored
C...       in valstr in case they prove to be needed later
C...       for an error estimate. the routine  newmsh
C...       filled in the needed values of the old solution
C...       in valstr.
C
       KNEW = ( 4 * (I-1) + 2 ) * MSTAR + 1
       KSTORE = ( 2 * (I-1) + 1 ) * MSTAR + 1
       X = XI(I) +  (XI(I+1)-XI(I)) * 2.D0 / 3.D0
       CALL APPROX_DAE (I, X, VALSTR(KNEW), DUMMY, ASAVE(1,3),
     +            DUMMY, XI, N, Z, DMZ,
     1            K, NCOMP, NY, MMAX, M, MSTAR, 4, DUMMY, 0)
       DO 20 L = 1,MSTAR
         ERR(L) = WGTERR(L) * DABS(VALSTR(KNEW) -
     1       VALSTR(KSTORE))
         KNEW = KNEW + 1
         KSTORE = KSTORE + 1
   20      CONTINUE
       KNEW = ( 4 * (I-1) + 1 ) * MSTAR + 1
       KSTORE = 2 * (I-1) * MSTAR + 1
       X = XI(I) +  (XI(I+1)-XI(I)) / 3.D0
       CALL APPROX_DAE (I, X, VALSTR(KNEW), DUMMY, ASAVE(1,2),
     +            DUMMY, XI, N, Z, DMZ,
     1            K, NCOMP, NY, MMAX, M, MSTAR, 4, DUMMY, 0)
       DO 30 L = 1,MSTAR
         ERR(L) = ERR(L) + WGTERR(L) * DABS(VALSTR(KNEW) -
     1       VALSTR(KSTORE))
         KNEW = KNEW + 1
         KSTORE = KSTORE + 1
   30      CONTINUE
C
C...       find component-wise maximum error estimate
C
       DO 40 L = 1,MSTAR
         ERREST(L) = DMAX1(ERREST(L),ERR(L))
   40      CONTINUE
C
C...       test whether the tolerance requirements are satisfied
C...       in the i-th interval.
C
       IF ( IFIN .EQ. 0 )                       GO TO 60
       DO 50 J = 1, NTOL
         LTOLJ = LTOL(J)
         LTJZ = LTOLJ  +  (I-1) * MSTAR
       IF ( ERR(LTOLJ) .GT.
     1          TOLIN(J) * (DABS(Z(LTJZ))+1.D0) )  IFIN = 0
   50      CONTINUE
   60 CONTINUE

      RETURN
C--------------------------------------------------------------
      END
C---------------------------------------------------------------------
C                            p a r t  3
C          collocation system setup routines
C---------------------------------------------------------------------
C
      SUBROUTINE LSYSLV_DAE(MSING, XI, XIOLD, Z, DMZ, DELZ, DELDMZ,
     1           G, W, V, FC, RHS, DMZO, INTEGS, IPVTG, IPVTW, RNORM,
     2        MODE, FSUB, DFSUB, GSUB, DGSUB, GUESS, ISING,RPAR,IPAR)
C*********************************************************************
C
C   purpose
C         this routine controls the set up and solution of a linear
C      system of collocation equations.
C         the matrix  g  is cast into an almost block diagonal
C      form by an appropriate ordering of the columns and solved
C      using the package of de boor-weiss [7] modified.
C      the matrix is composed of n blocks. the i-th block has the size
C                  integs(1,i) * integs(2,i).
C      it contains in its last rows the linearized collocation
C      equations, condensed as described in [4],
C      and the linearized side conditions corresponding to
C      the i-th subinterval.  integs(3,i)  steps of gaussian
C      elimination are applied to it to achieve a  partial plu
C      decomposition.  the right hand side vector is put into  rhs
C      and the solution vector is returned in  delz and deldmz.
C      note that the presence of algebraic solution components
C      does not affect the structure (size) of g -- only the contents
C      of the blocks (and the size of deldmz) changes.
C
C         lsyslv operates according to one of 5 modes:
C      mode = 0 - set up the collocation matrices  v , w , g
C                 and the right hand side  rhs ,  and solve.
C                 (for linear problems only.)
C      mode = 1 - set up the collocation matrices  v , w , g
C                 and the right hand sides  rhs  and  dmzo ,
C                 and solve. also set up  integs .
C                 (first iteration of nonlinear problems only).
C      mode = 2 - set up  rhs  only and compute its norm.
C      mode = 3 - set up  v, w, g  only and solve system.
C      mode = 4 - perform forward and backward substitution only
C                 (do not set up the matrices nor form the rhs).
C
C   variables
C
C      ig,izeta  - pointers to g,zeta respectively
C                       (necessary to keep track of blocks of g
C                       during matrix manipulations)
C      idmz,irhs,iv,iw - pointers to  rhs,v,w rspectively
C      df    - partial derivatives of f from dfsub
C      rnorm - euclidean norm of rhs
C      lside - number of side conditions in current and previous blocks
C      iguess = 1 when current soln is user specified via  guess
C             = 0 otherwise
C      dmzo  - an array used to project the initial solution into
C              the current pp-space, when mode=1.
C              in the case mode=1 the current solution iterate may not
C              be in the right space, being defined by an arbitrary
C              user's guess or as a pp on a different mesh.
C              when forming collocation equations we are using values
C              of z, y and dmval at collocation points and of z at
C              boundary points. at the end of lsyslv (with mode=1)
C              a similar projection used to obtain the corrections
C              delz and deldmz is used to obtain the projected initial
C              iterate z and dmz.
C      fc    - an array used to store projection matrices for
C              the case of projected collocation
C
C
C*********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  Z(*), DMZ(*), DELZ(*), DELDMZ(*), XI(*), XIOLD(*)
      DIMENSION  G(*), W(*), V(*), RHS(*), DMZO(*), DUMMY(1), Y(1)
      DIMENSION  INTEGS(3,*), IPVTG(*), IPVTW(*), YVAL(20)
      DIMENSION  ZVAL(40), F(40), DGZ(40), DMVAL(20), DF(800), AT(28)
      DIMENSION  FC(*), CB(400), IPVTCB(20),RPAR(*),IPAR(*)
C
      COMMON /DAEOUT/ PRECIS, IOUT, IPRINT
      COMMON /COLLOC/ RHO(7), COEF(49)
      COMMON /DAEORD/ K, NCOMP, NY, NCY, MSTAR, KD, KDY, MMAX, M(20)
      COMMON /DAESID/ ZETA(40), ALEFT, ARIGHT, IZETA, IZSAVE
      COMMON /DAEAPR/ N, NOLD, NMAX, NZ, NDMZ
      COMMON /DAENLN/ NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX
      COMMON /DAEBAS/ B(28), ACOL(28,7), ASAVE(28,4)
      integer nfunc, njac, nstep, nbound, njacbound
      common/CDAEdiag/nfunc, njac, nstep, nbound, njacbound
C
      EXTERNAL DFSUB, DGSUB
C
      IF (NY.EQ.0) THEN
       NYCB = 1
      ELSE
       NYCB = NY
      ENDIF
      INFC = (MSTAR+NY) * NCOMP
      M1 = MODE + 1

C..Francesca Mazzia Added initialization

      DO  I=1,MSTAR
        ZVAL(I) = 0.D0
      END DO
      DO  I=1,NY
         YVAL(I) = 0.D0
      ENDDO

      GO TO (10, 30, 30, 30, 310), M1
C
C...  linear problem initialization
C
   10 DO 20 I=1,MSTAR
   20   ZVAL(I) = 0.D0
      DO 25 I=1,NY
   25   YVAL(I) = 0.D0
C
C...  initialization
C
   30 IDMZ = 1
      IDMZO = 1
      IRHS = 1
      IG = 1
      IW = 1
      IV = 1
      IFC = 1
      IZETA = 1
      LSIDE = 0
      IOLD = 1
      NCOL = 2 * MSTAR
      RNORM = 0.D0
      IF ( MODE .GT. 1 )                            GO TO 80
C
C...  build integs (describing block structure of matrix)
C
      DO 70 I = 1,N
       INTEGS(2,I) = NCOL
       IF (I .LT. N)                            GO TO 40
       INTEGS(3,N) = NCOL
       LSIDE = MSTAR
       GO TO 60
   40      INTEGS(3,I) = MSTAR
   50      IF( LSIDE .EQ. MSTAR )                   GO TO 60
       IF ( ZETA(LSIDE+1) .GE. XI(I)+PRECIS )   GO TO 60
       LSIDE = LSIDE + 1
       GO TO 50
   60      NROW = MSTAR + LSIDE
   70      INTEGS(1,I) = NROW
   80 CONTINUE
      IF ( MODE .EQ. 2 )                            GO TO 90
C
C...  zero the matrices to be computed
C
      LW = KDY * KDY * N
      DO 84 L = 1, LW
   84   W(L) = 0.D0
C
C...  the do loop 290 sets up the linear system of equations.
C
  90  CONTINUE
C karline  
      nstep = nstep + 1  
      DO 290 I=1, N
C
C...       construct a block of  a  and a corresponding piece of  rhs.
C
       XII = XI(I)
       H = XI(I+1) - XI(I)
       NROW = INTEGS(1,I)
C
C...       go thru the ncomp collocation equations and side conditions
C...       in the i-th subinterval
C
  100      IF ( IZETA .GT. MSTAR )                   GO TO 140
       IF ( ZETA(IZETA) .GT. XII + PRECIS )      GO TO 140
C
C...       build equation for a side condition.
C
       IF ( MODE .EQ. 0 )                       GO TO 110
       IF ( IGUESS .NE. 1 )                     GO TO 102
C
C...       case where user provided current approximation
C
       CALL GUESS (XII, ZVAL, YVAL, DMVAL)
       GO TO 110
C
C...       other nonlinear case
C
  102      IF ( MODE .NE. 1 )                       GO TO 106
C.. Francesca Mazzia  y--> yval ?
       CALL APPROX_DAE (IOLD, XII, ZVAL, Y, AT, COEF, XIOLD, NOLD,
     1          Z, DMZ, K, NCOMP, NY, MMAX, M, MSTAR, 2, DUMMY, 0)
       GO TO 110
C.. Francesca Mazzia  y--> yval  ?
  106    CALL APPROX_DAE (I, XII, ZVAL, Y, AT, DUMMY, XI, N, Z, DMZ,
     1                  K, NCOMP, NY, MMAX, M, MSTAR, 1, DUMMY, 0)
  108      IF ( MODE .EQ. 3 )                       GO TO 120
C
C...       find  rhs  boundary value.
C
  110      CALL GSUB (IZETA, MSTAR, ZVAL, GVAL,RPAR,IPAR)
           nbound = nbound + 1  
       RHS(NDMZ+IZETA) = -GVAL
       RNORM = RNORM + GVAL**2
       IF ( MODE .EQ. 2 )                       GO TO 130
C
C...       build a row of  a  corresponding to a boundary point
C
  120      CALL GDERIV_DAE(G(IG), NROW, IZETA, ZVAL, DGZ, 1, DGSUB,
     & RPAR,IPAR)
  130      IZETA = IZETA + 1
       GO TO 100
C
C...       assemble collocation equations
C
  140      DO 220 J = 1, K
         HRHO = H * RHO(J)
         XCOL = XII + HRHO
C
C...         this value corresponds to a collocation (interior)
C...         point. build the corresponding  ncy  equations.
C
         IF ( MODE .EQ. 0 )                     GO TO 200
         IF ( IGUESS .NE. 1 )                   GO TO 160
C
C...         use initial approximation provided by the user.
C
         CALL GUESS (XCOL, ZVAL, YVAL, DMZO(IRHS) )
         GO TO 170
C
C...         find  rhs  values
C
  160        IF ( MODE .NE. 1 )                     GO TO 190
         CALL APPROX_DAE (IOLD, XCOL, ZVAL, YVAL, AT, COEF,
     +            XIOLD, NOLD, Z, DMZ,
     1            K, NCOMP, NY, MMAX, M, MSTAR, 2, DMZO(IRHS), 2)
C
  170        CALL FSUB (NCY, XCOL, ZVAL, YVAL, F,RPAR,IPAR)
c karline: added           
             nfunc = nfunc  + 1
  
         DO 175 JJ = NCOMP+1,NCY
  175          DMZO(IRHS+JJ-1) = 0.D0
         DO 180 JJ = 1, NCY
           VALUE = DMZO(IRHS) - F(JJ)
           RHS(IRHS) = - VALUE
           RNORM = RNORM + VALUE**2
           IRHS = IRHS + 1
  180        CONTINUE
         GO TO 210
C
C...         evaluate former collocation solution
C
C.. Francesca Mazzia  y--> yval ?
  190      CALL APPROX_DAE (I, XCOL, ZVAL, Y, ACOL(1,J), COEF, XI, N,
     1            Z, DMZ, K, NCOMP, NY, MMAX, M, MSTAR, 4, DUMMY, 0)
         IF ( MODE .EQ. 3 )                     GO TO 210
C
C...         fill in  rhs  values (and accumulate its norm).
C
         CALL FSUB (NCY,XCOL, ZVAL, DMZ(IRHS+NCOMP), F,RPAR,IPAR)
c karline: added           
             nfunc = nfunc  + 1
		 
         DO 195 JJ = 1, NCY
           VALUE = F(JJ)
           IF (JJ .LE. NCOMP) VALUE = VALUE - DMZ(IRHS)
           RHS(IRHS) =  VALUE
           RNORM = RNORM + VALUE**2
           IRHS = IRHS + 1
  195        CONTINUE
         GO TO 220
C
C...         the linear case
C
  200        CALL FSUB (NCY, XCOL, ZVAL, YVAL, RHS(IRHS),RPAR,IPAR)
c karline: added           
             nfunc = nfunc  + 1
  
         IRHS = IRHS + NCY
C
C...         fill in ncy rows of  w and v
C
  210        CALL VWBLOK_DAE(XCOL, HRHO, J, W(IW), V(IV), IPVTW(IDMZ),
     1            KDY, ZVAL, YVAL, DF, ACOL(1,J), DMZO(IDMZO),
     2            NCY, DFSUB, MSING,RPAR,IPAR)
         IF ( MSING .NE. 0 )                    RETURN
  220      CONTINUE
C
C...       build global bvp matrix  g
C
        IF (INDEX .NE. 1 .AND. NY .GT. 0) THEN
C
C...          projected collocation: find solution at xi(i+1)
C
          XI1 = XI(I+1)
          IF (MODE .NE. 0) THEN
           IF (IGUESS .EQ. 1) THEN
        CALL GUESS (XI1, ZVAL, YVAL, DMVAL)
           ELSE
        IF (MODE .EQ. 1) THEN
         CALL APPROX_DAE (IOLD, XI1, ZVAL, YVAL, AT,COEF,
     +            XIOLD, NOLD, Z, DMZ, K, NCOMP, NY, MMAX,
     +            M, MSTAR, 2, DUMMY, 1)
         IF (I .EQ. N)
     +            CALL APPROX_DAE (NOLD+1, XI1, ZVAL, YVAL, AT,COEF,
     +            XIOLD, NOLD, Z, DMZ, K, NCOMP, NY, MMAX,
     +            M, MSTAR, 1, DUMMY, 0)
        ELSE
         CALL APPROX_DAE (I, XI1, ZVAL, YVAL, AT,COEF,
     +            XI, N, Z, DMZ, K, NCOMP, NY, MMAX,
     +            M, MSTAR, 3, DUMMY, 1)
         CALL APPROX_DAE (I+1, XI1, ZVAL, YVAL, AT,COEF,
     +            XI, N, Z, DMZ, K, NCOMP, NY, MMAX,
     +            M, MSTAR, 1, DUMMY, 0)
        END IF
           END IF
          END IF
C
C...          find rhs at next mesh point (also for linear case)
C
          CALL FSUB (NCY,XI1, ZVAL, YVAL, F,RPAR,IPAR)
c karline: added           
             nfunc = nfunc  + 1
		  
        END IF
C
        CALL GBLOCK_DAE(H, G(IG), NROW, IZETA, W(IW), V(IV), KDY,
     2                  DUMMY, DELDMZ(IDMZ), IPVTW(IDMZ), 1, MODE,
     +                  XI1, ZVAL, YVAL, F, DF, CB, IPVTCB,
     +            FC(IFC), DFSUB, ISING, NCOMP, NYCB, NCY,RPAR,IPAR)
       IF (ISING .NE. 0)                       RETURN
       IF ( I .LT. N )                          GO TO 280
       IZSAVE = IZETA
  240      IF ( IZETA .GT. MSTAR )                  GO TO 290
C
C...       build equation for a side condition.
C
       IF ( MODE .EQ. 0 )                       GO TO 250
       IF ( IGUESS .NE. 1 )                     GO TO 245
C
C...       case where user provided current approximation
C
       CALL GUESS (ARIGHT, ZVAL, YVAL, DMVAL)
       GO TO 250
C
C...       other nonlinear case
C
  245      IF ( MODE .NE. 1 )                       GO TO 246
C.. Francesca Mazzia  y--> yval  ?
       CALL APPROX_DAE (NOLD+1, ARIGHT, ZVAL, Y, AT, COEF,
     +          XIOLD, NOLD, Z, DMZ,
     1          K, NCOMP, NY, MMAX, M, MSTAR, 1, DUMMY, 0)
       GO TO 250
C.. Francesca Mazzia  y--> yval  ?
  246      CALL APPROX_DAE (N+1, ARIGHT, ZVAL, Y, AT, COEF, XI, N,
     1       Z, DMZ, K, NCOMP, NY, MMAX, M, MSTAR, 1, DUMMY, 0)
  248      IF ( MODE .EQ. 3 )                       GO TO 260
C
C...       find  rhs  boundary value.
C
  250      CALL GSUB (IZETA, MSTAR, ZVAL, GVAL,RPAR,IPAR)
       nbound = nbound + 1  
       RHS(NDMZ+IZETA) = - GVAL
       RNORM = RNORM + GVAL**2
       IF ( MODE .EQ. 2 )                       GO TO 270
C
C...       build a row of  a  corresponding to a boundary point
C
  260   CALL GDERIV_DAE(G(IG), NROW, IZETA+MSTAR, ZVAL, DGZ, 2, DGSUB,
     &                   RPAR,IPAR )
  270      IZETA = IZETA + 1
       GO TO 240
C
C...       update counters -- i-th block completed
C
  280      IG = IG + NROW * NCOL
       IV = IV + KDY * MSTAR
       IW = IW + KDY * KDY
       IDMZ = IDMZ + KDY
       IF ( MODE .EQ. 1 )  IDMZO = IDMZO + KDY
       IFC = IFC + INFC+2*NCOMP
  290 CONTINUE
C
C...       assembly process completed
C
      IF ( MODE .EQ. 0 .OR. MODE .EQ. 3 )           GO TO 300
      RNORM = DSQRT(RNORM / FLOAT(NZ+NDMZ))
      IF ( MODE .EQ. 2 )                            RETURN
C
C...  solve the linear system.
C
C...  matrix decomposition
C
  300 CALL FCBLOK (G, INTEGS, N, IPVTG, DF, MSING)
C
C...  check for singular matrix
C
      MSING = - MSING
      IF( MSING .NE. 0 )                            RETURN
C
C...  perform forward and backward substitution .
C
  310 CONTINUE
      DO 311 L = 1, NDMZ
      DELDMZ(L) = RHS(L)
  311 CONTINUE
      IZ = 1
      IDMZ = 1
      IW = 1
      IFC = 1
      IZET = 1
      DO 320 I=1, N
       NROW = INTEGS(1,I)
       IZETA = NROW + 1 - MSTAR
       IF ( I .EQ. N ) IZETA = IZSAVE
  322    IF ( IZET .EQ. IZETA )                     GO TO 324
       DELZ(IZ-1+IZET) = RHS(NDMZ+IZET)
       IZET = IZET + 1
       GO TO 322
  324    H = XI(I+1) - XI(I)
       CALL GBLOCK_DAE(H, G(1), NROW, IZETA, W(IW), V(1), KDY,
     1                DELZ(IZ), DELDMZ(IDMZ), IPVTW(IDMZ), 2, MODE,
     +                XI1, ZVAL, YVAL, FC(IFC+INFC), DF, CB,
     +      IPVTCB, FC(IFC), DFSUB, ISING, NCOMP, NYCB, NCY,RPAR,IPAR)
       IZ = IZ + MSTAR
       IDMZ = IDMZ + KDY
       IW = IW + KDY * KDY
       IFC = IFC + INFC+2*NCOMP
       IF ( I .LT. N )                            GO TO 320
  326    IF ( IZET .GT. MSTAR )                     GO TO 320
       DELZ(IZ-1+IZET) = RHS(NDMZ+IZET)
       IZET = IZET + 1
       GO TO 326
  320 CONTINUE
C
C...  perform forward and backward substitution for mode=0,2, or 3.
C
      CALL SBBLOK (G, INTEGS, N, IPVTG, DELZ)
C
C...  finally find deldmz
C
      CALL DMZSOL (KDY, MSTAR, N, V, DELZ, DELDMZ)
C
      IF ( MODE .NE. 1 )                            RETURN
C
C...  project current iterate into current pp-space
C
      DO 321 L = 1, NDMZ
      DMZ(L) = DMZO(L)
  321 CONTINUE
      IZ = 1
      IDMZ = 1
      IW = 1
      IFC = 1
      IZET = 1
      DO 350 I=1, N
       NROW = INTEGS(1,I)
       IZETA = NROW + 1 - MSTAR
       IF ( I .EQ. N ) IZETA = IZSAVE
  330    IF ( IZET .EQ. IZETA )                     GO TO 340
       Z(IZ-1+IZET) = DGZ(IZET)
       IZET = IZET + 1
       GO TO 330
  340    H = XI(I+1) - XI(I)
       CALL GBLOCK_DAE(H, G(1), NROW, IZETA, W(IW), DF, KDY,
     1                Z(IZ), DMZ(IDMZ), IPVTW(IDMZ), 2, MODE,
     +                  XI1, ZVAL, YVAL, FC(IFC+INFC+NCOMP),
     +                  DF, CB, IPVTCB, FC(IFC), DFSUB, ISING,
     +            NCOMP, NYCB, NCY,RPAR,IPAR)
       IZ = IZ + MSTAR
       IDMZ = IDMZ + KDY
       IW = IW + KDY * KDY
       IFC = IFC + INFC+2*NCOMP
       IF ( I .LT. N )                            GO TO 350
  342    IF ( IZET .GT. MSTAR )                     GO TO 350
        Z(IZ-1+IZET) = DGZ(IZET)
        IZET = IZET + 1
       GO TO 342
  350 CONTINUE
      CALL SBBLOK (G, INTEGS, N, IPVTG, Z)
C
C...  finally find dmz
C
      CALL DMZSOL (KDY, MSTAR, N, V, Z, DMZ)
C
      RETURN
      END

      SUBROUTINE GDERIV_DAE(GI,NROW,IROW,ZVAL,DGZ,MODE,DGSUB,RPAR,IPAR)
C
C**********************************************************************
C
C   purpose:
C
C      construct a collocation matrix row according to mode:
C      mode = 1  -  a row corresponding to a initial condition
C                   (i.e. at the left end of the subinterval).
C      mode = 2  -  a row corresponding to a condition at  aright.
C
C   variables:
C
C      gi     - the sub-block of the global bvp matrix in
C               which the equations are to be formed.
C      nrow   - no. of rows in gi.
C      irow   - the row in gi to be used for equations.
C      zval   - z(xi)
C      dg     - the derivatives of the side condition.
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GI(NROW,*), ZVAL(*), DGZ(*), DG(40),RPAR(*),IPAR(*)
C
      COMMON /DAEORD/ KDUM, NDUM, NY, NCY, MSTAR, KD, KDY, MMAX, M(20)
      COMMON /DAESID/ ZETA(40), ALEFT, ARIGHT, IZETA, IDUM
      COMMON /DAENLN/ NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX
      integer nfunc, njac, nstep, nbound, njacbound
      common/CDAEdiag/nfunc, njac, nstep, nbound, njacbound
	  
C
C...  zero jacobian dg
C
      DO 10 J=1,MSTAR
   10   DG(J) = 0.D0
C
C...  evaluate jacobian dg
C
      CALL DGSUB (IZETA, mstar, ZVAL, DG,RPAR,IPAR)
      NJACBOUND = NJACBOUND + 1 
C
C...  evaluate  dgz = dg * zval  once for a new mesh
C
      IF (NONLIN .EQ. 0 .OR. ITER .GT. 0)           GO TO 30
      DOT = 0.D0
      DO 20 J = 1, MSTAR
   20   DOT = DOT  +  DG(J) * ZVAL(J)
      DGZ(IZETA) = DOT
C
C...  branch according to  m o d e
C
   30 IF ( MODE .EQ. 2 )                            GO TO 50
C
C...  provide coefficients of the j-th linearized side condition.
C...  specifically, at x=zeta(j) the j-th side condition reads
C...  dg(1)*z(1) + ... +dg(mstar)*z(mstar) + g = 0
C
C
C...  handle an initial condition
C
      DO 40 J = 1, MSTAR
      GI(IROW,J) =  DG(J)
   40 GI(IROW,MSTAR+J) = 0.D0
      RETURN
C
C...  handle a final condition
C
   50 DO 60 J= 1, MSTAR
      GI(IROW,J) = 0.D0
   60 GI(IROW,MSTAR+J) = DG(J)
      RETURN
      END
      SUBROUTINE VWBLOK_dae(XCOL, HRHO, JJ, WI, VI, IPVTW, KDY, ZVAL,
     1              YVAL, DF, ACOL, DMZO, NCY, DFSUB, MSING,RPAR,IPAR)
C
C**********************************************************************
C
C   purpose:
C
C      construct a group of  ncomp  rows of the matrices  wi  and  vi.
C      corresponding to an interior collocation point.
C
C
C   variables:
C
C      xcol   - the location of the collocation point.
C      jj     - xcol is the jj-th of k collocation points
C               in the i-th subinterval.
C      wi,vi  - the i-th block of the collocation matrix
C               before parameter condensation.
C      kdy    - no. of rows in vi and wi .
C      zval   - z(xcol)
C      yval   - y(xcol)
C      df     - the jacobian at xcol .
C      jcomp  - counter for the component being dealt with.
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WI(KDY,*), VI(KDY,*), ZVAL(*), DMZO(*), DF(NCY,*)
      DIMENSION IPVTW(*),  HA(7,4), ACOL(7,4), BASM(5), YVAL(*)
      DIMENSION RPAR(*),IPAR(*)
      COMMON /DAEORD/ K, NCOMP, NY, NDM, MSTAR, KD, KDYM, MMAX, M(20)
      COMMON /DAENLN/ NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX
      integer nfunc, njac, nstep, nbound, njacbound
      common/CDAEdiag/nfunc, njac, nstep, nbound, njacbound

	  
C
C...  initialize  wi
C
      I1 = (JJ-1) * NCY
      DO 10 ID = 1+I1, NCOMP+I1
      WI(ID,ID) = 1.D0
   10 CONTINUE
C
C...  calculate local basis
C
   30        FACT = 1.D0
         DO 35 L=1,MMAX
        FACT = FACT * HRHO / FLOAT(L)
        BASM(L) = FACT
        DO 35 J=1,K
           HA(J,L) = FACT * ACOL(J,L)
  35         CONTINUE
C
C... zero jacobian
C
      DO 40 JCOL = 1, MSTAR+NY
      DO 40 IR = 1, NCY
   40   DF(IR,JCOL) = 0.D0
C
C...  build ncy rows for interior collocation point x.
C...  the linear expressions to be constructed are:
C...   (m(id))
C...  u     -  df(id,1)*z(1) - ... - df(id,mstar)*z(mstar) -
C...   id
C...        -  df(id,mstar+1)*u(1) - ... - df(id,mstar+ny)*y(ny)
C...  for id = 1 to ncy  (m(id)=0 for id &gt; ncomp).
C
      CALL DFSUB (NCY, XCOL, ZVAL, YVAL, DF,RPAR,IPAR)
      njac = njac +1	  
      I0 = (JJ-1) * NCY
      I1 = I0 + 1
      I2 = I0 + NCY
C
C...  evaluate  dmzo = dmzo - df * (zval,yval)  once for a new mesh
C
      IF (NONLIN .EQ. 0 .OR. ITER .GT. 0)          GO TO 60
      DO 50 J = 1, MSTAR+NY
      IF (J .LE. MSTAR) THEN
      FACT = - ZVAL(J)
      ELSE
      FACT = - YVAL(J-MSTAR)
      END IF
      DO 50 ID = 1, NCY
      DMZO(I0+ID) = DMZO(I0+ID)  +  FACT * DF(ID,J)
  50  CONTINUE
C
C...  loop over the  ncomp  expressions to be set up for the
C...  current collocation point.
C
   60 DO 70 J = 1, MSTAR
      DO 70 ID = 1, NCY
      VI(I0+ID,J) = DF(ID,J)
   70 CONTINUE
      JN = 1
      DO 140 JCOMP = 1, NCOMP
       MJ = M(JCOMP)
       JN = JN + MJ
       DO 130 L = 1, MJ
        JV = JN - L
        JW = JCOMP
        DO 90 J = 1, K
          AJL = - HA(J,L)
          DO 80 IW = I1, I2
         WI(IW,JW) = WI(IW,JW)  +  AJL * VI(IW,JV)
   80         CONTINUE
   90       JW = JW + NCY
        LP1 = L + 1
        IF ( L .EQ. MJ )                        GO TO 130
        DO 110 LL = LP1, MJ
          JDF = JN - LL
          BL = BASM(LL-L)
          DO 100 IW = I1, I2
        VI(IW,JV) = VI(IW,JV)  +  BL * VI(IW,JDF)
  100         CONTINUE
  110       CONTINUE
  130    CONTINUE
  140 CONTINUE
C
C...  loop for the algebraic solution components
C
      DO 150 JCOMP = 1,NY
      JD = NCOMP+JCOMP
      DO 150 ID = 1,NCY
      WI(I0+ID,I0+JD) = -DF(ID,MSTAR+JCOMP)
  150 CONTINUE
      IF ( JJ .LT. K )                          RETURN
C
C   ...decompose the wi block and solve for the mstar columns of vi
C
C
C...  do parameter condensation
C
      MSING = 0
      CALL DGEFA  (WI, KDY, KDY, IPVTW, MSING)
C
C...   check for singularity
C
      IF ( MSING .NE. 0 )                         RETURN
      DO 250 J= 1,MSTAR
       CALL DGESL  (WI, KDY, KDY, IPVTW, VI(1,J), 0)
  250 CONTINUE
      RETURN
      END

      SUBROUTINE GBLOCK_DAE(H, GI, NROW, IROW, WI, VI, KDY,
     1                   RHSZ, RHSDMZ, IPVTW, MODE, MODL,
     +                   XI1, ZVAL, YVAL, F, DF, CB, IPVTCB,
     +             FC, DFSUB, ISING, NCOMP, NYCB, NCY,RPAR,IPAR)
C
C**********************************************************************
C
C   purpose:
C
C      construct collocation matrix rows according to mode:
C      mode = 1  -  a group of  mstar    rows corresponding
C                   to a mesh interval.
C           = 2  -  compute condensed form of rhs
C      modl = mode of lsyslv
C
C   variables:
C
C      h      - the  local stepsize.
C      gi     - the sub-block of the collocation matrix in
C               which the equations are to be formed.
C      wi     - the sub-block of noncondensed collocation equations,
C               left-hand side part.
C      vi     - the sub-block of noncondensed collocation equations,
C               right-hand side part.
C      rhsdmz - the inhomogenous term of the uncondensed collocation
C               equations.
C      rhsz   - the inhomogenous term of the condensed collocation
C               equations.
C      nrow   - no. of rows in gi.
C      irow   - the first row in gi to be used for equations.
C      xi1    - next mesh point (right end)
C      zval,yval - current solution at xi1
C      fc     - projection matrices
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION HB(7,4), BASM(5)
      DIMENSION GI(NROW,*), WI(*), VI(KDY,*)
      DIMENSION RHSZ(*), RHSDMZ(*), IPVTW(*),RPAR(*),IPAR(*)
      DIMENSION ZVAL(*), YVAL(*), F(*), DF(NCY,*), CB(NYCB,NYCB),
     +          IPVTCB(*), FC(NCOMP,*), BCOL(40), U(400), V(400)
C
      COMMON /DAEORD/  K, NCD, NY, NCYD, MSTAR, KD, KDUM, MMAX, M(20)
      COMMON /DAEBAS/ B(7,4), ACOL(28,7), ASAVE(28,4)
      COMMON /DAENLN/ NONLIN, ITER, LIMIT, ICARE, IGUESS, INDEX
      integer nfunc, njac, nstep, nbound, njacbound
      common/CDAEdiag/nfunc, njac, nstep, nbound, njacbound
C
C...  compute local basis
C
      FACT = 1.D0
      BASM(1) = 1.D0
      DO 30 L=1,MMAX
       FACT = FACT * H / FLOAT(L)
       BASM(L+1) = FACT
       DO 20 J=1,K
   20       HB(J,L) = FACT * B(J,L)
   30 CONTINUE
C
C...  branch according to  m o d e
C
      GO TO (40, 120 ), MODE
C
C...  set right gi-block to identity
C
   40 CONTINUE
      IF (MODL .EQ. 2)                          GO TO 110
      DO 60 J = 1, MSTAR
      DO 50 IR = 1, MSTAR
      GI(IROW-1+IR,J) = 0.D0
   50     GI(IROW-1+IR,MSTAR+J) = 0.D0
   60   GI(IROW-1+J,MSTAR+J) = 1.D0
C
C...  compute the block gi
C
      IR = IROW
      DO 100 ICOMP = 1, NCOMP
       MJ = M(ICOMP)
       IR = IR + MJ
       DO 90 L = 1, MJ
        ID = IR - L
        DO 80 JCOL = 1, MSTAR
           IND = ICOMP
           RSUM = 0.D0
           DO 70 J = 1, K
          RSUM = RSUM  -  HB(J,L) * VI(IND,JCOL)
   70             IND = IND + NCY
           GI(ID,JCOL) = RSUM
   80       CONTINUE
        JD = ID - IROW
        DO 85 LL = 1, L
           GI(ID,JD+LL) = GI(ID,JD+LL) - BASM(LL)
   85       CONTINUE
   90    CONTINUE
  100 CONTINUE

      IF (INDEX .EQ. 1 .OR. NY .EQ. 0)              RETURN
C
C...  projected collocation
C...  set up projection matrix and update gi-block
C
      CALL DFSUB (NCY, XI1, ZVAL, YVAL, DF,RPAR,IPAR)
      njac = njac + 1
C
C...  if index=2 then form projection matrices directly
C...  otherwise use svd to define appropriate projection
C
      IF (INDEX .EQ. 0) THEN
      CALL PRJSVD (FC,DF,CB,U,V,NCOMP,NCY,NY,IPVTCB,ISING,1)
      IF (ISING .NE. 0)                           RETURN
      ELSE
C
C...  form  cb
C
      DO 102 I=1,NY
      DO 102 J=1,NY
      FACT = 0
      ML = 0
      DO 101 L=1,NCOMP
        ML = ML + M(L)
        FACT = FACT + DF(I+NCOMP,ML) * DF(L,MSTAR+J)
  101     CONTINUE
      CB(I,J) = FACT
  102 CONTINUE
C
C...  decompose cb
C
      CALL DGEFA  (CB, NY, NY, IPVTCB, ISING)
      IF (ISING .NE. 0)                           RETURN
C
C...  form columns of fc
C
      DO 105 J=1,MSTAR+NY

      IF (J .LE. MSTAR) THEN
      DO 103 I=1,NY
  103      BCOL(I) = DF(I+NCOMP,J)
      ELSE
       DO 203 I=1,NY
  203      BCOL(I) = 0.0D0
       BCOL(J-MSTAR) = 1.D0
      END IF

      CALL DGESL  (CB, NY, NY, IPVTCB, BCOL, 0)

      DO 105 I=1,NCOMP
      FACT = 0.0D0
      DO 104 L=1,NY
        FACT = FACT + DF(I,L+MSTAR) * BCOL(L)
  104     CONTINUE
      FC(I,J) = FACT
  105 CONTINUE
C
      END IF
C
C...  update gi
C
      DO 108 J = 1,MSTAR
      DO 107 I=1,NCOMP
      FACT = 0
      DO 106 L=1,MSTAR
        FACT = FACT + FC(I,L) * GI(IROW-1+L,J)
  106     CONTINUE
      BCOL(I) = FACT
  107   CONTINUE
      ML = 0
      DO 108 I = 1,NCOMP
      ML = ML + M(I)
      GI(IROW-1+ML,J) = GI(IROW-1+ML,J) - BCOL(I)
  108 CONTINUE
C
C...  prepare extra rhs piece; two if new mesh
C
  110 IF (INDEX .EQ. 1 .OR. NY .EQ. 0 )           RETURN
      DO 115 JCOL=1,2
      DO 112 I=1,NCOMP
      FACT = 0
      DO 111  L=1,NY
        FACT = FACT + FC(I,L+MSTAR) * F(L+NCOMP)
  111     CONTINUE
      FC(I,JCOL+MSTAR+NY) = FACT
  112   CONTINUE
C
      IF (MODL .NE. 1 .OR. JCOL .EQ. 2)         RETURN
      DO 113 I = 1+NCOMP,NY+NCOMP
  113     F(I) = 0
      DO 114 J=1,MSTAR
      FACT = -ZVAL(J)
      DO 114 I = 1+NCOMP,NY+NCOMP
        F(I) = F(I) + DF(I,J) * FACT
  114   CONTINUE
  115 CONTINUE
      RETURN
C
C...  compute the appropriate piece of  rhsz
C
  120 CONTINUE
      CALL DGESL  (WI, KDY, KDY, IPVTW, RHSDMZ, 0)
      IR = IROW
      DO 140 JCOMP = 1, NCOMP
       MJ = M(JCOMP)
       IR = IR + MJ
       DO 130 L = 1, MJ
        IND = JCOMP
        RSUM = 0.D0
        DO 125 J = 1, K
           RSUM = RSUM  +  HB(J,L) * RHSDMZ(IND)
  125       IND = IND + NCY
        RHSZ(IR-L) = RSUM
  130    CONTINUE
  140 CONTINUE
      IF (INDEX .EQ. 1 .OR. NY .EQ. 0)              RETURN
C
C...  projected collocation
C...  calculate projected rhsz
C
      DO 160 I=1,NCOMP
      FACT = 0
      DO 150 L=1,MSTAR
      FACT = FACT + FC(I,L) * RHSZ(L+IROW-1)
  150   CONTINUE
      BCOL(I) = FACT
  160 CONTINUE
      ML = 0
      DO 170 I = 1,NCOMP
      ML = ML + M(I)
      RHSZ(IROW-1+ML) = RHSZ(IROW-1+ML) - BCOL(I) - F(I)
  170 CONTINUE
C
      RETURN
      END
      SUBROUTINE PRJSVD (FC,DF,D,U,V,NCOMP,NCY,NY,IPVTCB,ISING,MODE)
C
C**********************************************************************
C
C   purpose:
C
C      construct projection matrices in  fc  in the general case
C      where the problem may have a higher index but is not in
C      a Hessenberg index-2 form
C
C      calls the linpack routine  dsvdc  for a singular value
C      decomposition.
C
C      mode = 1 - called from gblock
C           = 2 - called from newmsh
C                 (then fc consists of only ncomp columns)
C
C**********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FC(NCOMP,*), DF(NCY,*), D(NY,NY), U(NY,NY), V(NY,NY)
      DIMENSION WORK(20), S(21), E(20), IPVTCB(*)
      COMMON /DAEORD/  K, NCD, NYD, NCYD, MSTAR, KD, KDUM, MMAX, M(20)
      COMMON /DAEOUT/ PRECIS, IOUT, IPRINT
      COMMON /DAEEST/ TOL(40), WGTMSH(40), WGTERR(40), TOLIN(40),
     1                ROOT(40), JTOL(40), LTOL(40), NTOL
c
C...  compute the maximum tolerance
C
      CHECK = 0.D0
      DO 5 I = 1, NTOL
    5   CHECK = DMAX1 ( TOLIN(I), CHECK )
C
C...  construct d and find its svd
C
      DO 10 I=1,NY
      DO 10 J=1,NY
      D(I,J) = DF(I+NCOMP,J+MSTAR)
   10 CONTINUE
      JOB = 11
      CALL DSVDC (D,NY,NY,NY,S,E,U,NY,V,NY,WORK,JOB,INFO)
C
C...  determine rank of d
C
      S(NY+1) = 0
      IRANK = 0
   20 IF (S(IRANK+1) .LT. CHECK)                    GO TO 30
      IRANK = IRANK + 1
      GO TO 20
C
C...  if d has full rank then no projection is needed
C
   30 IF (IRANK .EQ. NY) THEN
      DO 35 I=1,NCOMP
      DO 35 J=1,MSTAR+NY
   35       FC(I,J) = 0.D0
      RETURN
      ELSE
C
C...  form projected cb
C
      IR = NY-IRANK
      DO 50 I=1,NY
      DO 50 J=1,NY
      FACT = 0
      ML = 0
      DO 40 L=1,NCOMP
        ML = ML + M(L)
        FACT = FACT + DF(I+NCOMP,ML) * DF(L,MSTAR+J)
   40     CONTINUE
      D(I,J) = FACT
   50 CONTINUE
      DO 70 I=1,NY
      DO 60 J=1,IR
      WORK(J) = 0
      DO 60 L=1,NY
        WORK(J) = WORK(J) + D(I,L)*V(L,J+IRANK)
   60   CONTINUE
      DO 70 J=1,IR
      D(I,J) = WORK(J)
   70 CONTINUE
      DO 90 I=1,IR
      DO 80 J=1,IR
      WORK(J) = 0
      DO 80 L=1,NY
        WORK(J) = WORK(J) + U(L,I+IRANK)*D(L,J)
   80   CONTINUE
      DO 90 J=1,IR
       D(I,J) = WORK(J)
   90 CONTINUE
C
C...  decompose projected cb
C
      CALL DGEFA  (D, NY, IR, IPVTCB, ISING)
      IF (ISING .NE. 0)                           RETURN
C
C...  form columns of fc
C
      DO 130 J=MSTAR+1,MSTAR+NY
      DO 100 I=1,IR
  100     WORK(I) = U(J-MSTAR,I+IRANK)
      CALL DGESL  (D, NY, IR, IPVTCB, WORK, 0)
      DO 110 I=1,NY
      U(J-MSTAR,I) = 0
      DO 110 L=1,IR
        U(J-MSTAR,I) = U(J-MSTAR,I) + V(I,L+IRANK)*WORK(L)
  110   CONTINUE
      DO 130 I=1,NCOMP
      FACT = 0
      DO 120 L=1,NY
        FACT = FACT + DF(I,MSTAR+L)*U(J-MSTAR,L)
  120     CONTINUE
      FC(I,J) = FACT
  130 CONTINUE
C
      IF (MODE .EQ. 1) THEN
C
      DO 150 I=1,NCOMP
      DO 150 J=1,MSTAR
      FACT = 0
      DO 140 L=1,NY
        FACT = FACT +  FC(I,L+MSTAR) * DF(L+NCOMP,J)
  140     CONTINUE
      FC(I,J) = FACT
  150 CONTINUE
C
      ELSE
C
      DO 160 I=1,NCOMP
      MJ = 0
      DO 160 J=1,NCOMP
      MJ = MJ + M(J)
      FACT = 0
      DO 155  L=1,NY
        FACT = FACT +  FC(I,L+MSTAR) * DF(L+NCOMP,MJ)
  155     CONTINUE
      FC(I,J) = FACT
  160 CONTINUE
      END IF
C
      END IF
      RETURN
      END
C----------------------------------------------------------------------
C                             p a r t  4
C               polynomial and service routines
C----------------------------------------------------------------------
C
      SUBROUTINE APPSLN_DAE(X, Z, Y, FSPACE, ISPACE)
C
C*****************************************************************
C
C     purpose
C
C           set up a standard call to  approx  to evaluate the
C           approximate solution  z = z( u(x) ), y = y(x)  at a
C           point x (it has been computed by a call to  coldae ).
C           the parameters needed for  approx  are retrieved
C           from the work arrays  ispace  and  fspace .
C
C*****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Z(*), Y(*), FSPACE(*), ISPACE(*), A(28), DUMMY(1)
      IS6 = ISPACE(7)
      IS5 = ISPACE(1) + 2
      IS4 = IS5 + ISPACE(5) * (ISPACE(1) + 1)
      I = 1
      CALL APPROX_DAE (I,X,Z,Y,A, FSPACE(IS6), FSPACE(1), ISPACE(1),
     1             FSPACE(IS5), FSPACE(IS4), ISPACE(2), ISPACE(3),
     2             ISPACE(4), ISPACE(6), ISPACE(9), ISPACE(5), 2,
     3             DUMMY, 1)
      RETURN
      END
      SUBROUTINE APPROX_DAE(I, X, ZVAL, YVAL, A, COEF, XI, N, Z, DMZ,K,
     1                   NCOMP, NY, MMAX, M, MSTAR, MODE, DMVAL, MODM)
C
C**********************************************************************
C
C   purpose
C                                    (1)       (m1-1)     (mncomp-1)
C           evaluate z(u(x))=(u (x),u (x),...,u  (x),...,u  (x)      )
C                              1     1         1          mncomp
C           as well as optionally y(x) and dmval(x) at one point x.
C
C   variables
C     a      - array of mesh independent rk-basis coefficients
C     xi     - the current mesh (having n subintervals)
C     z      - the current solution vector (differential components).
C              it is convenient to imagine z as a two-dimensional
C              array with dimensions mstar x (n+1). then
C              z(j,i) = the jth component of z at the ith mesh point
C     dmz    - the array of mj-th derivatives of the current solution
C              plus algebraic solution components at collocation points
C              it is convenient to imagine dmz as a 3-dimensional
C              array with dimensions ncy x k x n. then
C              dmz(l,j,i) = a solution value at the jth collocation
C              point in the ith mesh subinterval: if l &lt;= ncomp then
C              dmz(l,j,i) is the ml-th derivative of ul, while if
C              l &gt; ncomp then dmz(l,j,i) is the value of the current
C              (l-ncomp)th component of y at this collocation point
C     mode   - determines the amount of initialization needed
C            = 4  forms z(u(x)) using z, dmz and ha
C            = 3  as in =4, but computes local rk-basis
C            = 2  as in =3, but determines i such that
C                       xi(i) .le. x .lt. xi(i+1) (unless x=xi(n+1))
C            = 1  retrieve  z=z(u(x(i)))  directly
C     modm   = 0  evaluate only zval
C            = 1  evaluate also yval
C            = 2  evaluate in addition dmval
C   output
C     zval   - the solution vector z(u(x)) (differential components)
C     yval   - the solution vector y(x)  (algebraic components)
C     dmval  - the mth derivatives of u(x)
C
C**********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ZVAL(*), DMVAL(*), XI(*), M(*), A(7,*), DM(7)
      DIMENSION Z(*), DMZ(*), BM(4), COEF(*), YVAL(*)
C
      COMMON /DAEOUT/ PRECIS, IOUT, IPRINT
C
      GO TO (10, 30, 80, 90), MODE
C
C...  mode = 1 , retrieve  z( u(x) )  directly for x = xi(i).
C
   10 X  = XI(I)
      IZ = (I-1) * MSTAR
      DO 20 J = 1, MSTAR
      IZ = IZ + 1
      ZVAL(J) = Z(IZ)
   20 CONTINUE
      RETURN
C
C...  mode = 2 ,  locate i so  xi(i) .le. x .lt. xi(i+1)
C
   30 CONTINUE
      IF ( X .GE. XI(1)-PRECIS .AND. X .LE. XI(N+1)+PRECIS )
     1                                              GO TO 40
      IF (IPRINT .LT. 1) THEN 
      CALL Rprintd3('Domain Error In Approx, X, Aleft, Aright ',
     +  X, Xi(1), Xi(N+1))
       ENDIF 

      IF ( X .LT. XI(1) )  X = XI(1)
      IF ( X .GT. XI(N+1) )  X = XI(N+1)
   40 IF ( I .GT. N .OR. I .LT. 1 )  I = (N+1) / 2
      ILEFT = I
      IF ( X .LT. XI(ILEFT) )                       GO TO 60
      DO 50 L = ILEFT, N
       I = L
       IF ( X .LT. XI(L+1) )                    GO TO 80
   50 CONTINUE
      GO TO 80
   60 IRIGHT = ILEFT - 1
      DO 70 L = 1, IRIGHT
       I = IRIGHT + 1 - L
       IF ( X .GE. XI(I) )                      GO TO 80
   70 CONTINUE
C
C...  mode = 2 or 3 , compute mesh independent rk-basis.
C
   80 CONTINUE
      S = (X - XI(I)) / (XI(I+1) - XI(I))
      CALL RKBAS ( S, COEF, K, MMAX, A, DM, MODM )
C
C...  mode = 2, 3, or 4 , compute mesh dependent rk-basis.
C
   90 CONTINUE
      BM(1) = X - XI(I)
      DO 95 L = 2, MMAX
       BM(L) = BM(1) / FLOAT(L)
   95 CONTINUE
C
C...  evaluate  z( u(x) ).
C
  100 IR = 1
      NCY = NCOMP + NY
      IZ = (I-1) * MSTAR + 1
      IDMZ = (I-1) * K * NCY
      DO 140 JCOMP = 1, NCOMP
      MJ = M(JCOMP)
      IR = IR + MJ
      IZ = IZ + MJ
      DO 130 L = 1, MJ
         IND = IDMZ + JCOMP
         ZSUM = 0.D0
         DO 110 J = 1, K
           ZSUM = ZSUM  +  A(J,L) * DMZ(IND)
  110          IND = IND + NCY
         DO 120 LL = 1, L
           LB = L + 1 - LL
  120          ZSUM = ZSUM * BM(LB)  +  Z(IZ-LL)
  130     ZVAL(IR-L) = ZSUM
  140 CONTINUE
      IF ( MODM .EQ. 0 )                            RETURN
C
C...  for modm = 1 evaluate  y(j) = j-th component of y.
C
      DO 150 JCOMP = 1, NY
  150    YVAL(JCOMP) = 0.D0
      DO 170 J = 1, K
       IND = IDMZ + (J-1)*NCY + NCOMP + 1
       FACT = DM(J)
       DO 160 JCOMP = 1, NY
        YVAL(JCOMP) = YVAL(JCOMP)  +  FACT * DMZ(IND)
        IND = IND + 1
  160    CONTINUE
  170 CONTINUE
      IF ( MODM .EQ. 1 )                            RETURN
C
C...  for modm = 2 evaluate  dmval(j) = mj-th derivative of uj.
C
      DO 180 JCOMP = 1, NCOMP
  180    DMVAL(JCOMP) = 0.D0
      DO 200 J = 1, K
       IND = IDMZ + (J-1)*NCY + 1
       FACT = DM(J)
       DO 190 JCOMP = 1, NCOMP
        DMVAL(JCOMP) = DMVAL(JCOMP)  +  FACT * DMZ(IND)
        IND = IND + 1
  190    CONTINUE
  200 CONTINUE
      RETURN
C--------------------------------------------------------------------
      END
C
      SUBROUTINE HORDER_DAE(I, UHIGH, HI, DMZ, NCOMP, NCY, K)
C
C**********************************************************************
C
C   purpose
C           determine highest order (piecewise constant) derivatives
C           of the current collocation solution
C
C   variables
C     hi     - the stepsize, hi = xi(i+1) - xi(i)
C     dmz    - vector of mj-th derivative of the solution
C     uhigh  - the array of highest order (piecewise constant)
C              derivatives of the approximate solution on
C              (xi(i),xi(i+1)), viz,
C                          (k+mj-1)
C              uhigh(j) = u   (x)    on (xi(i),xi(i+1))
C                          j
C
C**********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION UHIGH(*), DMZ(*)
C
      COMMON /COLLOC/ RHO(7), COEF(49)
C
      DN = 1.D0 / HI**(K-1)
C
C...  loop over the ncomp solution components
C
      DO 10 ID = 1, NCOMP
       UHIGH(ID) = 0.D0
   10 CONTINUE
      KIN = 1
      IDMZ = (I-1) * K * NCY + 1
      DO 30 J = 1, K
       FACT = DN * COEF(KIN)
       DO 20 ID = 1, NCOMP
        UHIGH(ID) = UHIGH(ID)  +  FACT * DMZ(IDMZ)
        IDMZ = IDMZ + 1
   20    CONTINUE
       KIN = KIN + K
   30 CONTINUE
      RETURN
      END

