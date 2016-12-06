c-----------------------------------------------------------------------
* karline: to make this code compatible with R:
* 1. change all write statements into rprint statements
* 2. changed all ( ,1) declarations into (,*)
* 3. added rpar, ipar to calls of guess
* 4. TODO added counters
* 5. renamed ipar -> iset and solutn -> guess
c all subroutines renamed, with 'sys' in front.
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
c                         p a r t  1
c        main storage allocation and program control subroutines
c-----------------------------------------------------------------------
c
      SUBROUTINE COLSYS (ncomp, m, aleft, aright, zeta, iset, ltol,
     1           tol, fixpnt, ispace, fspace, iflag, fsub,
     2           dfsub, gsub, dgsub, guess, RPAR, IPAR,
     3           icount)
c
c
C small changes by Karline Soetaert to make it compatible with R
c***********************************************************************
c
c     purpose
c
c     subroutine colsys solves a multi-point boundary value
c     problem for a mixed order system of ode-s given by
c
c           (m(i))
c          u       =  f  ( x; z(u(x)) )      i = 1, ... ,ncomp
c           i          i
c
c                                           aleft .lt. x .lt. aright,
c
c
c          g  ( zeta(j); z(u(zeta(j))) )= 0  j = 1, ... ,mstar
c           j
c                                    mstar=m(1)+m(2)+...+m(ncomp),
c
c
c          where                          t
c                u = (u , u , ... ,u     )  is the exact solution vector
c                      1   2        ncomp
c
c                 (mi)
c                u     is the mi=m(i) th  derivative of u
c                 i                                      i
c
c                                   (1)        (m1-1)       (mncomp-1) t
c                z(u(x)) = ( u (x),u  (x),...,u    (x),...,u      (x) )
c                             1     1          1            ncomp
c
c                f (x,z(u)) is a (generally) nonlinear function of
c                 i
c                           z(u)=z(u(x)).
c
c                g (zeta(j);z(u)) is a (generally) nonlinear function
c                 j
c                           used to represent a boundary condition.
c
c         the boundary points satisfy
c                aleft .le. zeta(1) .le. .. .le. zeta(mstar) .le. aright
c
c         the orders mi of the differential equations satisfy
c                m1 .le. m2 .le. ... .le. mncomp .le. 4.
c
c
c***********************************************************************
c
c     written by
c                  u. ascher,
c                            department of computer science,
c                            university of british columbia,
c                            vancouver, b. c., canada   v6t 1w5
c                  j. christiansen  and
c                  r. d. russell,
c                            mathematics department,
c                            simon fraser university,
c                            burnaby, b. c., canada  v5a 1s6
c
c***********************************************************************
c
c     method
c
c        the method used to approximate the solution u is
c     collocation at gaussian points, using b-splines of
c     order k+mi and continuity mi-1 in the i-th component,
c     i = 1, ..., ncomp. here, k is the number of collocation
c     points per subinterval and is chosen such that k .ge. m(ncomp).
c
c     main references
c
c     (1) u. ascher, j. christiansen and r.d. russell,
c
c        a collocation solver for mixed order
c        systems of boundary value problems
c
c        tech. rep. 77-13, dept. computer sc., univ. b.c.,
c        vancouver, canada.  also in math. comp. 33 (1979), 659-679.
c
c     (2) u. ascher, j. christiansen and r.d. russell,
c
c        colsys - a collocation code for boundary
c        value problems
c
c        proc. conf. for codes for bvp-s in ode-s,
c        houston, texas, 1978.
c
c     other references
c
c     (3) u. ascher and r. d. russell
c
c        evaluation of b-splines for solving systems
c        of boundary value problems
c
c        tech. rep. 77-14, dept. computer sc., univ. b.c.,
c        vancouver, canada.
c
c     (4) c. deboor and r. weiss
c
c        solveblok: a package for solving almost block diagonal
c        linear systems, with applications to spline approximation
c        and the numerical solution of ordinary differential equations
c
c        mrc tech report 1625, university of wisconsin - madison
c
c
c     (5) r. d. russell and j. christiansen
c
c        adaptive mesh selection strategies for
c        solving boundary value problems
c
c        siam j. numer. anal. 7(1978), 59-80.
c
c***********************************************************************
c
c     ***************     input to colsys     ***************
c
c     variables
c
c     ncomp - no. of differential equations   (ncomp .le. 20)
c
c     m(j) - order of the j-th differential equation ( m(j).le.m(j+1)
c            and mstar = m(1) + ... + m(ncomp) .le. 40 )
c
c     aleft - left end of interval
c
c     aright - right end of interval
c
c     zeta(j) - j-th side condition point (boundary point). must
c               have  zeta(j) .le. zeta(j+1)
c
c     iset - an integer array dimensioned at least 11.
c            a list of the parameters in iset and their meaning follows.
c            some parameters are renamed in colsys; these new names are
c            given in parentheses.
c
c     iset(1)     ( = nonlin )
c             = 0 if the problem is linear
c             = 1 if the problem is nonlinear
c
c     iset(2) = no. of collocation points per subinterval  (= k )
c               where m(ncomp) .lt.  k .le. 7 . if iset(2)=0 then
c               colsys sets  k = max ( m(ncomp)+1, 5-m(ncomp) )
c
c     iset(3) = no. of subintervals in the initial mesh  ( = n ).
c               if iset(3) = 0 then colsys arbitrarily sets n = 5.
c
c     iset(4) = no. of solution and derivative tolerances.  ( = ntol )
c               we require  0 .lt. ntol .le. mstar.
c
c     iset(5) = dimension of fspace.     ( = ndimf )
c
c     iset(6) = dimension of ispace.     ( = ndimi )
c
c     iset(7) -  output control ( = iprint )
c              = -1 for full diagnostic printout  - karline:should never be -1
c              = 0 for selected printout
c              = 1 for no printout
c
c     iset(8)     ( = iread )
c             = 0 causes colsys to generate a uniform initial mesh.
c             = 1 if the initial mesh is provided by the user.  it
c                 is defined in fspace as follows:  the mesh
c                 aleft=x(1).lt.x(2).lt. ... .lt.x(n).lt.x(n+1)=aright
c                 will occupy  fspace(1), ..., fspace(n+1). the
c                 user needs to supply only the interior mesh
c                 points  fspace(j) = x(j), j = 2, ..., n.
c             = 2 if the initial mesh is supplied by the user
c                 as with iset(8)=1, and in addition no adaptive
c                 mesh selection is to be done.
c
c     iset(9)     ( = iguess )
c             = 0 if no initial guess for the solution is
c                 provided.
c             = 1 if an initial guess is provided by the user
c                 in subroutine  guess.
c             = 2 if an initial mesh and approximate solution
c                 coefficients are provided by the user in  fspace.
c                 (the former and new mesh are the same).
c             = 3 if a former mesh and an approximate solution
c                 coefficients are provided by the user in fspace,
c                 and the new mesh is to be taken twice as coarse;
c                 i.e.,every second point from the former mesh.
c             = 4 if in addition to a former initial mesh and an
c                 approximate solution coefficients, a new mesh
c                 is provided in fspace as well.
c                 (see description of output for further details
c                 on iguess = 2, 3, and 4.)
c
c     iset(10)= 0 if the problem is regular
c             = 1 if the first relax factor is =rstart, and the
c                 nonlinear iteration does not rely on past covergence
c                 (use for an extra sensitive nonlinear problem only).
c             = 2 if we are to return immediately upon  (a) two
c                 successive nonconvergences, or  (b) after obtaining
c                 error estimate for the first time.
c
c     iset(11)= no. of fixed points in the mesh other than
c               aleft and aright. ( = nfxpnt , the dimension of fixpnt)
c
c     ltol  -  an array of dimension iset(4). ltol(j) = l  specifies
c              that the j-th tolerance in  tol  controls the error
c              in the l-th component of z(u).   also require that
c              1.le.ltol(1).lt.ltol(2).lt. ... .lt.ltol(ntol).le.mstar
c
c     tol    - an array of dimension iset(4). tol(j) is the
c              error tolerance on the ltol(j) -th component
c              of z(u). thus, the code attempts to satisfy
c              for j=1,...,ntol  on each subinterval
c              abs(z(v)-z(u))       .le. tol(j)*abs(z(u))       +tol(j)
c                            ltol(j)                     ltol(j)
c
c              if v(x) is the approximate solution vector.
c
c     fixpnt - an array of dimension iset(11).   it contains
c              the points, other than aleft and aright, which
c              are to be included in every mesh.
c
c     ispace - an integer work array of dimension iset(6).
c              its size provides a constraint on nmax,
c              the maximum number of subintervals. choose
c              iset(6) according to the formula
c                      iset(6)  .ge.  nmax*nsizei
c                where
c                      nsizei = 3 + kdm - nrec
c                with
c                      kdm = kd + mstar  ;  kd = k * ncomp ;
c                      nrec = no. of right end boundary conditions.
c
c
c     fspace - a real work array of dimension iset(5).
c              its size provides a constraint on nmax.
c              choose iset(5) according to the formula
c                      iset(5)  .ge.  nmax*nsizef
c                where
c                      nsizef = 4 + k + 2 * kd + (4+2*k) * mstar +
c                              (kdm-nrec) * (kdm+1).
c
c
c     iflag - the mode of return from colsys.
c           = 1 for normal return
c           = 0 if the collocation matrix is singular.
c           =-1 if the expected no. of subintervals exceeds storage
c               specifications.
c           =-2 if the nonlinear iteration has not converged.
c           =-3 if there is an input data error.
c
c
c***********************************************************************
c
c     *************    user supplied subroutines   *************
c
c
c     the following subroutines must be declared external in the
c     main program which calls colsys.
c
c
C    karline: changed from FSUB (x , z , f)
C     FSUB  - name of subroutine for evaluating f(x,z(u(x))) =
C                            t
C             (f ,...,f     )  at a point x in (aleft,aright).  it
C               1      ncomp
C             should have the heading
c                  subroutine FSUB (mstar, x , z , f, rpar, ipar)
C
C             where f is the vector containing the value of fi(x,z(u))
C             in the i-th component and                            t
C                                       z(u(x))=(z(1),...,z(mstar))
C             is defined as above under  purpose .
C
C
C     karline: changed from subroutine dFSUB (x , z , df)
C     DFSUB - name of subroutine for evaluating the jacobian of
C             f(x,z(u)) at a point x.  it should have the heading
C                dfsub (mstar, x , z , df,rpar, ipar)
C              
C             where z(u(x)) is defined as for FSUB and the (ncomp) by
C             (mstar) array df should be filled by the partial deriv-
C             atives of f, viz, for a particular call one calculates
C                                df(i,j) = dfi / dzj, i=1,...,ncomp
C                                                     j=1,...,mstar.
C
C
C     karline: changed from subroutine GSUB (i , z , g)
C     GSUB  - name of subroutine for evaluating the i-th component of
C             g(x,z(u(x))) = g (zeta(i),z(u(zeta(i)))) at a point x =
C                             i
C             zeta(i) where 1.le.i.le.mstar. it should have the heading
C
C                      subroutine GSUB (i , mstar,  z , g, rpar, ipar)
C
C             where z(u) is as for FSUB, and i and g=g  are as above.
C                                                     i
C             note that in contrast to f in  FSUB , here
C             only one value per call is returned in g.
C
C
C     karline: changed from subroutine dGSUB (i , z , dg)
C     dGSUB - name of subroutine for evaluating the i-th row of
C             the jacobian of g(x,u(x)).  it should have the heading
C
C                      subroutine dGSUB(i, mstar, z, dg, rpar, ipar)
C
C             where z(u) is as for FSUB, i as for GSUB and the mstar-
C             vector dg should be filled with the partial derivatives
C             of g, viz, for a particular call one calculates
C                   dg(i,j) = dgi / dzj      j=1,...,mstar.
C
C
c     guess- name of subroutine to evaluate the initial
c             approximation for  z(u(x)) and for dmval(u(x))= vector
c             of the mj-th derivatives of u(x). it should have the
c             heading
c
c                       subroutine guess (x , z , dmval, rpar, ipar)
c
c             note that this subroutine is needed only if using
c             iset(9) = 1, and then all  mstar  components of z
c             and  ncomp  components of  dmval  should be specified
c             for any x,  aleft .le. x .le. aright .
c
c
c***********************************************************************
c
c     ************   use of output from colsys   ************
c
c                 ***   solution evaluation   ***
c
c     on return from colsys, the arrays fspace and ispace
c     contain information specifying the approximate solution.
c     the user can produce the solution vector  z( u(x) )  at
c     any point x, aleft .le. x .le. aright, by the statement,
c
c           call sysappsln (x, z, fspace, ispace)
c
c     when saving the coefficients for later reference, only
c     ispace(1),...,ispace(7+ncomp)    and
c     fspace(1),...,fspace(ispace(7))    need to be saved as
c     these are the quantities used by sysappsln.
c
c
c                 ***   simple continuation   ***
c
c
c     a formerly obtained solution can easily be used as the
c     first approximation for the nonlinear iteration for a
c     new problem by setting   (iguess =) iset(9) = 2, 3 or 4.
c
c     if the former solution has just been obtained then the
c     values needed to define the first approximation are
c     already in ispace and fspace.
c     alternatively, if the former solution was obtained in a
c     previous run and its coefficients were saved then those
c     coefficients must be put back into
c     ispace(1),..., ispace(7+ncomp)    and
c     fspace(1),..., fspace(ispace(7)).
c
c     for iset(9) = 2 or 3 set iset(3) = ispace(1) ( = the
c     size of the previous mesh ).
c
c     for iset(9) = 4 the user specifies a new mesh of n subintervals
c     as follows.
c     the values in  fspace(1),...,fspace(ispace(7))  have to be
c     shifted by n+1 locations to  fspace(n+2),..,fspace(ispace(7)+n+1)
c     and the new mesh is then specified in fspace(1),..., fspace(n+1).
c     also set iset(3) = n.
c
c
c***********************************************************************
c
c     ***************      package subroutines      ***************
c
c     the following description gives a brief overview of how the
c     procedure is broken down into the subroutines which make up
c     the package called  colsys . for further details the
c     user should refer to documentation in the various subroutines
c     and to the references cited above.
c
c     the subroutines fall into four groups:
c
c part 1 - the main storage allocation and program control subroutines.
c
c     colsys - tests input values, does initialization and breaks up
c              the work areas, fspace and ispace, into the arrays
c              used by the program.
c
c     syscontrl - is the actual driver of the package. this routine
c              contains the strategy for nonlinear problems.
c
c
c part 2 - mesh selection and error estimation subroutines
c
c     sysconsts - is called once by  colsys  to initialize constants
c              which are used for error estimation and mesh selection.
c
c     sysnewmsh - generates meshes. it contains the test to decide
c              whether or not to redistribute a mesh.
c
c     syserrchk - produces error estimates and checks against the
c              tolerances at each subinterval
c
c
c part 3 - collocation system set-up subroutines
c
c     syslsyslv - controls the set-up and solution of the linear
c              algebraic systems of collocation equations which
c              arise at each newton iteration.
c
c     sysbldblk - is used by syslsyslv to set up the equation(s) associated
c              with a side condition point or a collocation point.
c
c
c part 4 - b-spline subroutines
c
c     sysappsln - sets up a standard call to  sysapprox .
c
c     sysapprox - evaluates a piecewise polynomial solution.
c
c     sysbspfix - evaluates the mesh independent b-splines
c              (i.e. the fixed b-splines)
c
c     sysbspvar - evaluates the mesh dependent b-splines (i.e. the
c              varying b-splines)
c
c     sysbspder - generates values for the derivatives needed to set
c              up the collocation equations.
c
c     sysappdif - generates a divided difference table from the b-spline
c              coefficients for a collocation solution. the table
c              is used in  sysapprox .
c
c     syshorder - evaluates the highest order derivatives of the
c              current collocation solution used for mesh refinement.
c
c
c     to solve the linear systems of collocation equations
c     constructed in part 3,  colsys  uses the package  solveblok
c     of de boor and weiss (to appear in toms).
c
c
c----------------------------------------------------------------------
      implicit real(kind=8) (a-h,o-z)
      common /order/ k,nc,mstar,kd,kdm,mnsum,mt(20)
      common /appr/ n,nold,nmax,nalpha,mshflg,mshnum,mshlmt,mshalt
      common /side/ tzeta(40),tleft,tright,izeta
      common /nonln/ precis,nonlin,iter,limit,icare,iprint,iguess,ifreez
      common /eqord/  ind(5), ineq(20), mnd(5), nd, neq
      common /errors/ ttl(40),wgtmsh(40),tolin(40),root(40),
     1       jtol(40),lttol(40),ntol
      external fsub, dfsub, gsub, dgsub, guess
      dimension m(*), zeta(*), iset(*), ltol(*), tol(*),
     1     fixpnt(*), ispace(*), fspace(*), ipar(*), rpar(*), icount(*)

c karline: added this 'dimenioning'
      dimension dum1(1),dum2(1),dum3(1),dum4(1)
c karline: added counters
      integer nfunc, njac, nstep, nbound, njacbound
      common/coldiag/nfunc, njac, nstep, nbound, njacbound

c***********************************************************************
c
c     the actual subroutine colsys serves as an interface with
c     the package of subroutines referred to collectively as
c     colsys. the subroutine serves to test some of the input
c     parameters, rename some of the parameters (to make under-
c     standing of the coding easier), to do some initialization,
c     and to break the work areas fspace and ispace up into the
c     arrays needed by the program.
c
c***********************************************************************

C     intialise counters
      nfunc = 0
      njac = 0
      nstep = 0
      nbound = 0
      njacbound = 0

c
c...  compute machine
c...  dependent constant  precis = 100 * machine unit roundoff
c
      precis = 1.d0
   10 precis = precis / 2.d0
      precp1 = precis + 1.d0
      if (precp1 .gt. 1.d0)                         go to 10
      precis = precis * 100.d0
c
c...  in case incorrect input data is detected, the program returns
c...  immediately with iflag=-3.
c
      iflag = -3
      if (ncomp .lt. 1 .or. ncomp .gt. 20)          return
      if (m(1) .lt. 1 .or. m(ncomp) .gt. 4)         return
      if (ncomp .eq. 1)                             go to 30
      do 20 i=2,ncomp
           if (m(i-1) .gt. m(i))                    return
   20 continue
   30 continue
c
c...  rename some of the parameters and set default values.
c
      nonlin = iset(1)
      k = iset(2)
      if (k .eq. 0)   k = max0( m(ncomp)+1, 5-m(ncomp) )
      n = iset(3)
      if (n .eq. 0)  n = 5
      iread = iset(8)
      iguess = iset(9)
      if (nonlin .eq. 0 .and. iguess .eq. 1) iguess = 0
      if (iguess .ge. 2 .and. iread .eq. 0)  iread = 1
      icare = iset(10)
      ntol = iset(4)
      ndimf = iset(5)
      ndimi = iset(6)
      nfxpnt = iset(11)
      iprint = iset(7)
      mstar = 0
      mnsum = 0
      do  40 i=1,ncomp
           mnsum = mnsum + m(i)**2
           mstar = mstar + m(i)
   40 Continue
      do 50 i=1,ncomp
        mt(i) = m(i)
   50 Continue
      do 60 i=1,mstar
        tzeta(i) = zeta(i)
   60 Continue
      do 70 i=1,ntol
           lttol(i) = ltol(i)
           tolin(i) = tol(i)
   70 Continue
      tleft = aleft
      tright = aright
      nc = ncomp
      kd = k * ncomp
      kdm = kd + mstar
c
c...  print the input data for checking.
c



c
c...  check for correctness of data
c
      if (k .lt. 0 .or. k .gt. 7)                   return
      if (n .lt. 0)                                 return
      if (iread .lt. 0 .or. iread .gt. 2)           return
      if (iguess .lt. 0 .or. iguess .gt. 4)         return
      if (icare .lt. 0 .or. icare .gt. 2)           return
      if (ntol .lt. 0 .or. ntol .gt. mstar)         return
      if (nfxpnt .lt. 0)                            return
      if (iprint .lt. (-1) .or. iprint .gt. 1)      return
      if (mstar .lt. 0 .or. mstar .gt. 40)          return
c
c...  set limits on iterations and initialize counters.
c...  limit = maximum number of newton iterations per mesh.
c...  see subroutine  sysnewmsh  for the roles of  mshlmt , mshflg ,
c...  mshnum , and  mshalt .
c
      mshlmt = 3
      mshflg = 0
      mshnum = 1
      mshalt = 1
      limit = 40
c
c...  compute the maxium possible n for the given sizes of
c...  ispace  and  fspace.
c
      nrec = 0
      do 110 ii=1,mstar
           i = mstar + 1 - ii
           if (zeta(i) .lt. aright)                 go to 110
           nrec = ii
  110 continue
      nfixi = nrec
      nsizei = 3 + kdm - nrec
      nfixf = nrec * (kdm+1) + 2 * mnsum + 2 * mstar + 3
      nsizef = 4 + k + 2 * kd + (4+2*k) * mstar +
     1(kdm-nrec) * (kdm+1)
      nmaxf = (ndimf - nfixf) / nsizef
      nmaxi = (ndimi - nfixi) / nsizei
      if (iprint .lt. 1) then
       CALL Rprinti1('The maximum number of subintervals is min',nmaxF)
       CALL Rprinti1('The maximum number allowed from ispace', NMAXI)
      endif 
      nmax = min0(nmaxf,nmaxi)
      if (nmax .lt. n)                              return
      if (nmax .lt. nfxpnt+1)                       return
      if (nmax .lt. 2*nfxpnt+2  .and.  iprint .lt. 1) then
      CALL Rprint('Insufficient space to double mesh for err estimate')
      endif 
c
c...  generate pointers to break up  fspace  and  ispace .
c
      lxi = 1
      la = lxi + nmax + 1
      lxiold = la + kdm * (nmax * (kdm-nrec) + nrec)
      lxij = lxiold + nmax + 1
      lalpha = lxij + k * nmax
      ldlpha = lalpha + nmax * kd + mstar
      lelpha = ldlpha + nmax * kd + mstar
      laldif = lelpha + nmax * k * mstar + mnsum
      lrhs = laldif + nmax * k * mstar + mnsum
      lvalst = lrhs + nmax * (kdm - nrec) + nrec
      lslope = lvalst + 4 * mstar * nmax
      laccum = lslope + nmax
      lipiv = 1
      linteg = lipiv + (lvalst - lrhs)
c
c...  if  iguess .ge. 2, move  xiold  and  aldif  to their proper
c...  locations in  fspace.
c
      if (iguess .lt. 2)                            go to 160
      nold = n
      if (iguess .eq. 4)  nold = ispace(1)
      naldif =  nold * k * mstar + mnsum
      np1 = n + 1
      if (iguess .eq. 4)  np1 = np1 + nold + 1
      do 120 i=1,naldif
        fspace( laldif+i-1 )  =  fspace( np1+i )
  120 Continue
      np1 = nold + 1
      if (iguess .eq. 4)                            go to 140
      do 130 i=1,np1
        fspace( lxiold+i-1 )  =  fspace( lxi+i-1 )
  130 Continue
      go to 160
  140 do 150 i=1,np1
        fspace( lxiold+i-1 )  =  fspace( n+1+i )
  150 Continue
  160 continue
c
c...  initialize collocation points, constants, mesh.
c
      call sysconsts
c karline: this is a strange call; sysnewmsh expects arrays of certain length; here dum..
      call sysnewmsh (3+iread, fspace(lxi), fspace(lxiold),
     1     fspace(lxij), dum1, dum2, dum3, dum4,
     2     nfxpnt, fixpnt)
c
c...  determine which are the different order equations and
c...  put these orders in  mnd , also generate the pointers
c...  ind  and  ineq  which will be used in  sysbspder .
c
      ind(1) = 1
      mnd(1) = m(1)
      nd = 1
      neq = 0
      ig = (m(1)+1) * (m(1)+k) + 1
      if (ncomp .le. 1)                             go to 200
      do 190 j=2,ncomp
           mj = m(j)
           if (mj .eq. m(j-1))                      go to 170
           nd = nd + 1
           ind(nd) = ig
           mnd(nd) = mj
           go to 180
  170      neq = neq + 1
           ineq(neq) = ig
  180      ig = ig + (mj+1) * (mj+k)
  190 continue
      ind(nd+1) =ind(nd) + ig
  200 continue
c
c...  determine first approximation, if the problem is nonlinear.
c
      if (iguess .ge. 2)                            go to 230
      np1 = n + 1
      do 210 i = 1,np1
        fspace(i + lxiold - 1) = fspace(i + lxi - 1)
  210 Continue
      nold = n
      if (nonlin .eq. 0 .or. iguess .eq. 1)         go to 230
c
c...  system provides first approximation of the solution.
c...  choose  z(j) = 0   for j=1,..,mstar.
c
      do 220 i = 1,nalpha
        fspace(i + lalpha - 1) = 0.d0
  220 Continue
      call sysappdif (fspace(laldif), fspace(lalpha), fspace(lxi),
     1     n, k, nc, mt, mstar)
  230 continue
      if (iguess .ge. 2)  iguess = 0
      call syscontrl (fspace(lxi),fspace(lxiold),fspace(lxij),
     1     fspace(lalpha),fspace(laldif),fspace(lrhs),
     2     fspace(ldlpha), fspace(lelpha),
     3     fspace(la),fspace(lvalst),fspace(lslope),
     4     fspace(laccum),ispace(lipiv),ispace(linteg),
     5     nfxpnt,fixpnt,iflag,fsub,dfsub,gsub,dgsub,
     6     guess,rpar,ipar)
c...  prepare output
      ispace(1) = n
      ispace(2) = k
      ispace(3) = ncomp
      ispace(4) = mstar
      naldif = n * k * mstar + mnsum
      ispace(5) = naldif
      ispace(6) = naldif + n + 2
      ispace(7) = ispace(6) + 65
      do 240 i=1,ncomp
        ispace(7+i) = m(i)
  240 Continue
      do 250 i=1,naldif
        fspace(n+1+i) = fspace(laldif-1+i)
  250 Continue
      icount(1) = nfunc
      icount(2) = njac
      icount(3) = nbound
      icount(4) = njacbound
      icount(5) = nstep
      return
c-----------------------------------------------------------------------
      end
      subroutine syscontrl(xi, xiold, xij, alpha, aldif, rhs,
     1           dalpha, ealpha, a, valstr, slope,
     2           accum, ipiv, integs, nfxpnt, fixpnt, iflag,
     3           fsub, dfsub, gsub, dgsub, guess,rpar,ipar)
c
c**********************************************************************
c
c   purpose
c     this subroutine is the actual driver.  the nonlinear iteration
c     strategy is controlled here ( see (2) ).  upon convergence, syserrchk
c     is called to test for satisfaction of the requested tolerances.
c
c   variables
c
c     check  - maximum tolerance value, used as part of criteria for
c              checking for nonlinear iteration convergence
c     relax  - the relaxation factor for damped newton iteration
c     relmin - minimum allowable value for relax  (otherwise the
c              jacobian is considered singular).
c     rlxold - previous relax
c     rstart - initial value for relax when problem is sensitive
c     ifrz   - number of fixed jacobian iterations
c     lmtfrz - maximum value for ifrz before performing a reinversion
c     iter   - number of iterations (counted only when jacobian
c              reinversions are performed).
c     xi     - current mesh
c     xiold  - previous mesh
c     ipred  = 0  if relax is determined by a correction
c            = 1  if relax is determined by a prediction
c     ifreez = 0  if the jacobian is to be inverted
c            = 1  if the jacobian is currently fixed (frozen)
c     icon   = 0  if no previous convergence has been obtained
c            = 1  if convergence on a previous mesh has been obtained
c     icare  =-1  no convergence occurred (used for regular problems)
c            = 0  a regular problem
c            = 1  a sensitive problem
c            = 2  used for continuation (see description of iset(10)
c                 in colsys).
c     rnorm  - norm of rhs (right hand side) for current iteration
c     rnold  - norm of rhs for previous iteration
c     anscl  - scaled norm of newton correction
c     anfix  - scaled norm of newton correction at next step
c     anorm  - scaled norm of a correction obtained with jacobian fixed
c     naldif - number of components of aldif (see subroutine approx)
c     imesh  - a control variable for subroutines sysnewmsh and syserrchk
c
c***********************************************************************
c
      implicit real(kind=8) (a-h,o-z)
      external fsub, dfsub, gsub, dgsub, guess
      dimension xi(*), xiold(*), xij(*), alpha(*), aldif(*), rhs(*)
      dimension a(*), valstr(*), slope(*), accum(*), ipiv(*), integs(*)
      dimension dalpha(*), ealpha(*) , fixpnt(*),rpar(*),ipar(*)
      common /order/ k,ncomp,mstar,kd,kdm,mnsum,m(20)
      common /appr/ n,nold,nmax,nalpha,mshflg,mshnum,mshlmt,mshalt
      common /side/ zeta(40),aleft,aright,izeta
      common /nonln/ precis,nonlin,iter,limit,icare,iprint,iguess,ifreez
      common /eqord/  ind(5), ineq(20), mnd(5), nd, neq
      common /errors/ tol(40),wgtmsh(40),tolin(40),root(40),
     1       jtol(40),ltol(40),ntol
c
c...  constants for control of nonlinear iteration
c
      relmin = 1.d-3
      rstart = 1.d-2
      lmtfrz = 4
      anscl = 1.d0
      ifrz = 0
c
c...  compute the maximum tolerance
c
      check = 0.d0
      do 10 i=1,ntol
        check = dmax1 (tolin(i), check )
   10 Continue
      falpha = dfloat(nalpha)
      imesh = 1
      icon = 0
      if (nonlin .eq. 0) icon=1
      icor = 0
      lconv = 0
c
c...  the main iteration begins here
c...  loop 20 is executed until error tolerances are satisfied or
c...  the code fails (due to a singular matrix or storage limitations)
c
   20      continue
c
c...       initialization for a new mesh
c
           iter = 0
           naldif = n * k * mstar + mnsum
           if (nonlin .gt. 0)                       go to 60
c
c...       the linear case.
c...       set up and solve equations
c
           call syslsyslv (iflag, xi, xiold, xij, alpha, aldif, rhs,
     1          ealpha, a, ipiv, integs, rnorm, 0, fsub,
     2          dfsub, gsub, dgsub, guess,rpar,ipar)
c
c...       check for a singular matrix
c
           if (iflag .ne. 0)                        go to 40
   30      if (iprint .lt. 1) then
             CALL Rprint('The matrix is singular')
           endif 

           return
c
c...       update the old mesh
c
   40      np1 = n + 1
           do 50 i=1,np1
             xiold(i) = xi(i)
   50      Continue
           nold = n
c
c...       prepare table of divided differences and call  syserrchk
c
           call sysappdif (aldif, alpha, xi, n, k, ncomp, m, mstar)
           go to 450
c
c...       iteration loop for nonlinear case
c...       define the initial relaxation parameter (= relax)
c
   60      relax = 1.d0
c
c...       check for previous convergence and problem sensitivity
c
           if (icare .eq. 1 .or. icare .eq. (-1))  relax = rstart
           if (icon .eq. 0)                         go to 140
c
c...       convergence on a previous mesh has been obtained.    thus
c...       we have a very good initial approximation for the newton
c...       process.    proceed with one full newton and then iterate
c...       with a fixed jacobian.
c
           ifreez = 0
c
c...       evaluate right hand side and its norm
c
           call syslsyslv (iflag, xi, xiold, xij, dalpha, aldif, rhs,
     1          alpha, a, ipiv, integs, rnorm, 1, fsub,
     2          dfsub, gsub, dgsub, guess,rpar,ipar)
c
c...       solve for the next iterate .
c...       the value of ifreez determines whether this is a full
c...       newton step (=0) or a fixed jacobian iteration (=1).
c
           if (iprint .lt. 0  .and.  iter .eq. 0)  then
             CALL Rprint('Fixed Jacobian iterations')
           endif 

   70      if (iprint .lt. 0)  then
      CALL Rprintid('Iteration = , Norm (RHS) = ', Iter, Rnorm)
           endif 

           rnold = rnorm
           call syslsyslv (iflag, xi, xiold, xij, dalpha, aldif, rhs,
     1          alpha, a, ipiv, integs, rnorm, 2+ifreez ,
     2          fsub, dfsub, gsub, dgsub, guess,rpar,ipar)
c
c...       check for a singular matrix
c
           if (iflag .eq. 0)                        go to 30
           if (ifreez .eq. 1)                       go to 90
c
c...       a full newton step
c
           iter = iter + 1
           ifrz = 0
c
c...       update the old mesh.
c
           np1 = n + 1
           do 80 i=1,np1
             xiold(i) = xi(i)
   80      Continue
           nold = n
   90      continue
c
c...       update   alpha , compute new  rhs  and its norm
c
           do 100 i=1,nalpha
             alpha(i) = alpha(i) + dalpha(i)
  100      Continue
           call sysappdif (aldif, alpha, xi, n, k, ncomp, m, mstar)
           call syslsyslv (iflag, xi, xiold, xij, dalpha, aldif, rhs,
     1          alpha, a, ipiv, integs, rnorm, 1, fsub,
     2          dfsub, gsub, dgsub, guess,rpar,ipar)
c
c...       check monotonicity. if the norm of  rhs  gets smaller,
c...       proceed with a fixed jacobian; else proceed cautiously,
c...       as if convergence has not been obtained before (icon=0).
c
           if (rnorm .lt. precis)                   go to 405
           if (rnorm .le. rnold)                    go to 120
           if (iprint .lt. 0)  then
      CALL Rprintid('Iteration = , Norm (RHS) = ', Iter, Rnorm)
             CALL Rprint('Switch to damped Newton iteration')
           endif 

           icon = 0
           relax = rstart
           do 110 i=1,nalpha
             alpha(i) = alpha(i) - dalpha(i)
  110      Continue
           call sysappdif (aldif, alpha, xi, n, k, ncomp, m, mstar)
           iter = 0
           go to 140
  120      if (ifreez .eq. 1)                       go to 130
           ifreez = 1
           go to 70
c
c...       verify that the linear convergence with fixed jacobian
c...       is fast enough.
c
  130      ifrz = ifrz + 1
           if (ifrz .ge. lmtfrz) ifreez = 0
           if (rnold .lt. 4.d0*rnorm) ifreez = 0
           go to 300
c
c...       no previous convergence has been obtained. proceed
c...       with the modified newton method.
c...       evaluate rhs.
c
  140      if(iprint .lt. 0)  then
             CALL Rprint('Full damped Newton iteration')
           endif 
       call syslsyslv (iflag, xi, xiold, xij, dalpha, aldif, rhs,
     1          alpha, a, ipiv, integs, rnorm, 1, fsub,
     2          dfsub, gsub, dgsub, guess,rpar,ipar)
c
c...       find a newton direction
c
  150      rnold = rnorm
           if (iter .ge. limit)                     go to 420
           call syslsyslv (iflag, xi, xiold, xij, dalpha, aldif, rhs,
     1          alpha, a, ipiv, integs, rnorm, 2, fsub,
     2          dfsub, gsub, dgsub, guess,rpar,ipar)
c
c...       check for a singular matrix
c
           if (iflag .eq. 0)                        go to 30
           if (iter .gt. 0)                         go to 170
c
c...       bookkeeping for first mesh
c
           if ( iguess .eq. 1) iguess = 0
c
c...       update the old mesh
c
           np1 = n + 1
           do 160 i=1,np1
             xiold(i) = xi(i)
  160      Continue
           nold = n
           go to 190
  170      continue
c
c...       predict relaxation factor for newton step.
c
           andif = 0.d0
           do 180 i=1,nalpha
             andif = andif + (ealpha(i) - dalpha(i))**2
     1     / (alpha(i)*alpha(i) + precis)
  180      Continue
           relax = relax * anscl / dmax1( dsqrt(andif/falpha),
     1     precis)
           if (relax .gt. 1.d0)  relax = 1.d0
  190      rlxold = relax
           ipred = 1
           iter = iter + 1
c
c...       determine a new  alpha  and find new  rhs  and its norm           
c
           do 200 i=1,nalpha
             alpha(i) = alpha(i) + relax * dalpha(i)
  200      Continue
  210      call sysappdif (aldif, alpha, xi, n, k, ncomp, m, mstar)
           call syslsyslv (iflag, xi, xiold, xij, dalpha, aldif, rhs,
     1          alpha, a, ipiv, integs, rnorm, 1, fsub,
     2          dfsub, gsub, dgsub, guess,rpar,ipar)
c
c...       compute a fixed jacobian iterate (used to control relax)
c
           call syslsyslv (iflag, xi, xiold, xij, ealpha, aldif, rhs,
     1          alpha, a, ipiv, integs, rnorm, 3, fsub,
     2          dfsub, gsub, dgsub, guess,rpar,ipar)
c
c...       find scaled norms of various terms used to correct relax
c
           anorm = 0.d0
           anfix = 0.d0
           anscl = 0.d0
           do 220 i=1,nalpha
             anscl = anscl + dalpha(i) * dalpha(i) /
     1       (alpha(i)*alpha(i) + precis)
             scale = alpha(i) - relax*dalpha(i)
             scale = 1.d0 / (scale*scale + precis)
             anorm = anorm + dalpha(i) * dalpha(i) * scale
             anfix = anfix + ealpha(i) * ealpha(i) * scale
  220      Continue
           anorm = dsqrt(anorm / falpha)
           anfix = dsqrt(anfix / falpha)
           anscl = dsqrt(anscl / falpha)
           if (icor .eq. 1)                         go to 230
           if (iprint .lt. 0)  then
      CALL Rprintid('Iteration = , Relaxation factor = ', Iter, Relax)
      CALL Rprintd2('Norm of scaled RHS changes from, to ',Anorm,Anfix)
      CALL Rprintd2('Norm of RHS changes from, to ', Rnold, Rnorm)
           endif 

           go to 240
  230      if (iprint .lt. 0) then
      CALL Rprintd1('Relaxation factor corrected to ', Relax)
      CALL Rprintd2('Norm of scaled RHS changes from, to ',Anorm,Anfix)
      CALL Rprintd2('Norm of RHS changes from, to ', Rnold, Rnorm)
           endif 

  240      icor = 0
c
c...       check for monotonic decrease in  dalpha.
c
           if (anfix.lt.precis .or. rnorm.lt.precis)go to 405
           if (anfix .gt. anorm)                    go to 250
c
c...       we have a decrease. if  ealpha  small, check convergence
c
           if (anfix .le. check)                    go to 290
c
c...       correct the predicted  relax  unless the corrected
c...       value is within 10 percent of the predicted one.
c
           if (ipred .ne. 1)                        go to 150
  250      if (iter .ge. limit)                     go to 420
c
c...       correct the relaxation factor.
c
           ipred = 0
           arg = (anfix/anorm - 1.d0) / relax + 1.d0
           if (arg .lt. 0.d0)                         go to 150
           if (arg .le. .25d0*relax+.125d0*relax**2 ) go to 260
           factor = -1.d0 + dsqrt (1.d0+8.d0 * arg)
           if ( dabs(factor-1.d0) .lt. .1d0*factor )  go to 150
           relax = relax / factor
           go to 270
  260      if (relax .ge. .9d0)                       go to 150
           relax = 1.d0
  270      icor = 1
           if (relax .lt. relmin)                     go to 430
           do 280 i=1,nalpha
            alpha(i) = alpha(i) + (relax-rlxold) * dalpha(i)
  280      Continue
           rlxold = relax
           go to 210
c
c...       check convergence.
c...       compute divided difference tables for correction terms.
c
  290      call sysappdif (a, ealpha, xi, n, k, ncomp, m, mstar)
           go to 310
c
c...       if  icon = 1 then also save  a.
c
  300      call sysappdif (ealpha, dalpha, xi, n, k, ncomp, m, mstar)
  310      continue
           inn = 0
           jcol = 0
           jinit = 1
           do 380 i = 1, ntol
             jend = jtol(i) - 1
             if (jend .lt. jinit)                   go to 330
             do 320 j = jinit, jend
               mj = m(j)
               nalphj = n * k + mj
               jcol = jcol + mj
               inn = inn + mj * nalphj
  320        continue
  330        jinit = jend + 1
             nalphj = n * k + m(jinit)
             inn1 = inn
             jcol1 = jcol + 1
  340        if (jcol1 .eq. ltol(i))                go to 350
             inn1 = inn1 + nalphj
             jcol1 = jcol1 + 1
             go to 340
  350        iinit = jcol1 - jcol
c
c...         check that tolerances are satisfied for b-spline coeffs.
c
             do 370 ii = iinit, nalphj
               in = inn1 + ii
               if (icon .eq. 1)                     go to 360
               if (dabs(a(in)) .gt. tolin(i) *
     1              (dabs(aldif(in)) +1.d0))        go to 410
               go to 370
  360          if (dabs(ealpha(in)) .gt. tolin(i) *
     1              (dabs(aldif(in)) + 1.d0))       go to 410
  370        continue
  380      continue
c
c...       convergence obtained
c
           if (iprint .lt. 1) then 
      CALL Rprinti1('Convergence after iteration ',Iter)
           endif 

           if (icon .eq. 1)                         go to 450
c
c...       since convergence obtained, update coeffs with term from
c...       the fixed jacobian iteration.
c
           do 390 i=1,naldif
             aldif(i) = aldif(i) + a(i)
  390      Continue
           do 400 i=1,nalpha
             alpha(i) = alpha(i) + ealpha(i)
  400      Continue
  405    if ((anfix.lt.precis.or.rnorm.lt.precis).and.iprint.lt.1) then
      CALL Rprinti1('Convergence after iteration ',Iter)
           endif 

           icon = 1
           if (icare .eq. (-1))  icare = 0
           go to 450
c
c...       no convergence. repeat
c
  410      if ( icon .eq. 0)                        go to 150
           go to 70
c
c...       diagnostics for failure of nonlinear iteration.
c
  420      if(iprint .lt. 1) then
      CALL Rprinti1('NO convergence after iteration ',Iter)
           endif 
           go to 440
  430      if(iprint .lt. 1) then
      CALL Rprintd1('NO convergence, Relaxation factor too small',Relax)
      CALL Rprintd1('Should not be less than ', Relmin)
           endif 
  440      iflag = -2
           lconv = lconv + 1
           if (icare .eq. 2 .and. lconv .gt. 1)     return
           if (icare .eq. 0)  icare = -1
           go to 460
c
c...       check for error tolerance satisfaction
c
  450      call syserrchk(imesh,xiold,aldif,valstr,mstar,ifin)
           if (imesh .eq. 1 .or. ifin .eq. 0 .and.
     1          icare .ne. 2)                       go to 460
           iflag = 1
           return
c
c...       pick a new mesh
c...       check safeguards for mesh refinement
c
  460      imesh = 1
           if (icon .eq. 0 .or. mshnum .ge. mshlmt
     1     .or. mshalt .ge. mshlmt)  imesh = 2
           if (mshalt .ge. mshlmt .and. mshnum .lt. mshlmt)
     1     mshalt = 1
           call sysnewmsh(imesh, xi, xiold, xij, aldif, valstr,
     1          slope, accum, nfxpnt, fixpnt)
c
c...       exit if expected n is too large (but may try n=nmax once)
c
           if (n .le. nmax)                         go to 470
           n = n / 2
           iflag = -1
           if (icon .eq. 0 .and. iprint .lt. 1) then
             CALL Rprint('NO convergence')
           endif 
           if (icon .eq. 1 .and. iprint .lt. 1) then
      CALL Rprint('Probably tolerances too stringent or nmax too small')
           endif 
           return
  470      if (icon .eq. 0)  imesh = 1
           if (icare .eq. 1)  icon = 0
      go to 20
c     ---------------------------------------------------------------
      end
c-----------------------------------------------------------------------
c                    p a r t  2
c          mesh selection, error estimation, (and related
c          constant assignment) routines -- see (1), (2), (5)
c-----------------------------------------------------------------------
c
      subroutine sysnewmsh (mode, xi, xiold, xij, aldif, valstr,
     1           slope, accum, nfxpnt, fixpnt)
c
c***********************************************************************
c
c   purpose
c            select a mesh on which a collocation solution is to be
c            determined
c
c                           there are 5 possible modes of action:
c            mode = 5,4,3 - deal mainly with definition of an initial
c                           mesh for the current boundary value problem
c                 = 2,1   - deal with definition of a new mesh, either
c                           by simple mesh halving or by mesh selection
c            more specifically, for
c            mode = 5  an initial (generally nonuniform) mesh is
c                      defined by the user and no mesh selection is to
c                      be performed
c                 = 4  an initial (generally nonuniform) mesh is
c                      defined by the user
c                 = 3  a simple uniform mesh (except possibly for some
c                      fixed points) is defined; n= no. of subintervals
c                 = 1  the automatic mesh selection procedure is used
c                      (see (1) and (5) for details)
c                 = 2  a simple mesh halving is performed
c
c***********************************************************************
c
c   variables
c
c            n      = number of mesh subintervals
c            nold   = number of subintervals for former mesh
c            xi     - mesh point array
c            xiold  - former mesh point array
c            mshlmt - maximum no. of mesh selections which are permitted
c                     for a given n before mesh halving
c            mshnum - no. of mesh selections which have actually been
c                     performed for the given n
c            mshalt - no. of consecutive times ( plus 1 ) the mesh
c                     selection has alternately halved and doubled n.
c                     if mshalt .ge. mshlmt then  syscontrl  requires
c                     that the current mesh be halved.
c            mshflg =  1  the mesh is a halving of its former mesh
c                       (so an error estimate has been calculated)
c                   = 0  otherwise
c            iguess - iset(9) in subroutine colsys.  it is used
c                     here only for mode=5 and 4, where
c                   = 2 the subroutine sets xi=xiold.  this is
c                       used e.g. if continuation is being per-
c                       formed, and a mesh for the old differen-
c                       tial equation is being used
c                   = 3 same as for =2, except xi uses every other
c                       point of xiold (so mesh xiold is mesh xi
c                       halved)
c                   = 4 xi has been defined by the user, and an old
c                       mesh xiold is also available
c                       otherwise, xi has been defined by the user
c                       and we set xiold=xi in this subroutine
c            slope  - an approximate quantity to be equidistributed for
c                     mesh selection (see (1)), viz,
c                             .                        (k+mj)
c                     slope(i)=     max   (weight(l) *u      (xi(i)))
c                               1.le.l.le.ntol         j
c
c                     where j=jtol(l)
c            slphmx - maximum of slope(i)*(xiold(i+1)-xiold(i)) for
c                     i = 1 ,..., nold.
c            accum  - accum(i) is the integral of  slope  from  aleft
c                     to  xiold(i).
c            valstr - is assigned values needed in  syserrchk  for the
c                     error estimate.
c***********************************************************************
c
      implicit real(kind=8) (a-h,o-z)
      common /order/k,ncomp,mstar,kd,kdm,mnsum,m(20)
      common /appr/ n,nold,nmax,nalpha,mshflg,mshnum,mshlmt,mshalt
      common /errors/ tol(40),wgtmsh(40),tolin(40),root(40),
     1       jtol(40),ltol(40),ntol
      common /colloc/ rho(7),wgterr(40)
      common /side/  zeta(40), aleft, aright, izeta
      common /nonln/ precis,nonlin,iter,limit,icare,iprint,iguess,ifreez
      common /bsplin/ vncol(66,7), vnsave(66,5), vn(66)
      dimension d1(40), d2(40), zv(40), slope(*), accum(*), valstr(*)
c karline: added dimension of dummy double precision 'dumm'
      dimension xi(*), xiold(*), xij(*), aldif(*), fixpnt(*), dumm(1)
c
      nfxp1 = nfxpnt +1

      if (mode .eq. 1) then
        goto  180
      else if (mode .eq. 2) then
        goto  100
      else if (mode .eq. 3) then
        goto  50
      else if (mode .eq. 4) then
        goto  20
      else if (mode .eq. 5) then
        goto  10
      endif 
C      go to (180, 100, 50, 20, 10), mode
c
c...  mode=5   set mshlmt=1 so that no mesh selection is performed
c
   10 mshlmt = 1
c
c...  mode=4   the user-specified initial mesh is already in place.
c
   20 if (iguess .lt. 2)                            go to 40
c
c...  iguess=2, 3 or 4.
c
      noldp1 = nold + 1

      if (iguess .ne. 3)                            go to 40
c
c...  if iread ( iset(8) ) .ge. 1 and iguess ( iset(9) )
c...  .eq. 3 then the first mesh is every second point of the
c...  mesh in  xiold .
c
      n = nold /2
      i = 0
      do 30 j = 1, nold, 2
           i = i + 1
           xi(i) = xiold(j)
   30  Continue
   40 continue
      np1 = n + 1
      xi(1) = aleft
      xi(np1) = aright
      go to 320
c
c...  mode=3   generate a (piecewise) uniform mesh. if there are
c...  fixed points then ensure that the n being used is large enough.
c
   50 if ( n .lt. nfxp1 ) n = nfxp1
      np1 = n + 1
      xi(1) = aleft
      ileft = 1
      xleft = aleft
c
c...  loop over the subregions between fixed points.
c
      do 90 j = 1,nfxp1
           xright = aright
           iright = np1
           if ( j .eq. nfxp1 )                      go to 60
           xright = fixpnt(j)
c
c...       determine where the j-th fixed point should fall in the
c...       new mesh - this is xi(iright) and the (j-1)st fixed
c...       point is in xi(ileft)
c
           nmin = int((xright-aleft)/(aright-aleft)*dfloat(n) + 1.5d0)
           if (nmin .gt. n-nfxpnt+j)  nmin = n - nfxpnt + j
           iright = max0 (ileft+1, nmin)
   60      xi(iright) = xright
c
c...       generate equally spaced points between the j-1st and the
c...       j-th fixed points.
c
           nregn = iright - ileft - 1
           if ( nregn .eq. 0 )                      go to 80
           dx = (xright - xleft) / dfloat(nregn+1)
           do 70 i = 1,nregn
             xi(ileft+i) = xleft + dfloat(i) * dx
   70      Continue
   80      ileft = iright
           xleft = xright
   90 continue
      go to 320
c
c...  mode=2  halve the current mesh (i.e. double its size)
c
  100 n2 = 2 * n
c
c...  check that n does not exceed storage limitations
c
      if (n2 .le. nmax)                             go to 120
c
c...  if possible, try with n=nmax. redistribute first.
c
      if (mode .eq. 2)                              go to 110
      n = nmax / 2
      go to 220
  110 if (iprint .lt. 1)  then 
             CALL Rprint('Expected N too large')
      endif 
      n = n2
      return
c
c...  calculate the old approximate solution values at
c...  points to be used in  syserrchk  for error estimates.
c...  if  mshflg  =1 an error estimate was obtained for
c...  for the old approximation so half the needed values
c...  will already be in  valstr .
c
  120 if (mshflg .eq. 0)                            go to 140
c
c...  save in  valstr  the values of the old solution
c...  at the relative positions 1/6 and 5/6 in each subinterval.
c
      kstore = 1
      do 130 i = 1,nold
           hd6 = (xiold(i+1) - xiold(i)) / 6.d0
           x = xiold(i) + hd6
           call sysapprox (i, x, valstr(kstore), vnsave(1,2), xiold,
     1          nold, aldif, k, ncomp, m, mstar, 3,dumm,0)
           x = x + 4.d0 * hd6
           kstore = kstore + 3 * mstar
           call sysapprox (i, x, valstr(kstore), vnsave(1,5),
     1          xiold, nold, aldif, k, ncomp, m, mstar, 3,dumm,0)
           kstore = kstore + mstar
  130 continue
      go to 160
c
c...  save in  valstr  the values of the old solution
c...  at the relative positions 1/6, 2/6, 4/6 and 5/6 in
c...  each subinterval.
c
  140 kstore = 1
      do 155 i = 1,n
           x = xi(i)
           hd6 = (xi(i+1) - xi(i)) / 6.d0
           do 150 j = 1,4
             x = x + hd6
             if ( j.eq.3 ) x = x + hd6
             call sysapprox (i, x, valstr(kstore), vnsave(1,j+1),
     1            xiold, nold, aldif, k, ncomp, m, mstar, 3,
     2            dumm,0)
             kstore = kstore + mstar
  150 continue
  155 continue
  160 mshflg = 0
      mshnum = 1
      mode = 2
c
c...  generate the halved mesh.
c
      j = 2
      do 170 i = 1,n
           xi(j) = (xiold(i) + xiold(i+1)) / 2.d0
           xi(j+1) = xiold(i+1)
           j = j + 2
  170  Continue
      n = n2
      go to 320
c
c...  mode=1  we do mesh selection if it is deemed worthwhile
c
  180 if ( nold .eq. 1 )                            go to 100
      if (nold .le. 2*nfxpnt)                       go to 100
c
c...  the first interval has to be treated separately from the
c...  other intervals (generally the solution on the (i-1)st and ith
c...  intervals will be used to approximate the needed derivative, but
c...  here the 1st and second intervals are used.)
c
      i = 1
      call syshorder (1, d1, xiold, aldif)
      call syshorder (2, d2, xiold, aldif)
      call sysapprox (i, xiold(i), zv, vnsave(1,1), xiold, nold,
     1     aldif, k, ncomp, m, mstar, 3, dumm, 0)
      accum(1) = 0.d0
      slope(1) = 0.d0
      oneovh = 2.d0 / ( xiold(3) - xiold(1) )
      do 190 j = 1,ntol
           jj = jtol(j)
           jv = ltol(j)
       slope(1) = dmax1(slope(1),(dabs(d2(jj)-d1(jj))*wgtmsh(j)*
     1  oneovh / (1.d0 + dabs(zv(jv)))) **root(j))
  190 Continue
      slphmx = slope(1) * (xiold(2) - xiold(1))
      accum(2) = slphmx
      iflip = 1
c
c...  go through the remaining intervals generating  slope
c...  and  accum .
c
      do 210 i = 2,nold
           if ( iflip .eq. (-1) ) call syshorder ( i, d1, xiold, aldif)
           if ( iflip .eq. 1 ) call syshorder ( i, d2, xiold, aldif)
           call sysapprox (i, xiold(i), zv, vnsave(1,1), xiold, nold,
     1          aldif, k, ncomp, m, mstar, 3, dumm, 0)
           oneovh = 2.d0 / ( xiold(i+1) - xiold(i-1) )
           slope(i) = 0.d0
c
c...       evaluate function to be equidistributed
c
           do 200 j = 1,ntol
             jj = jtol(j)
             jv = ltol(j)
           slope(i) = dmax1(slope(i),(dabs(d2(jj)-d1(jj))*wgtmsh(j)*
     1     oneovh / (1.d0 + dabs(zv(jv)))) **root(j))
  200      Continue
c
c...       accumulate approximate integral of function to be
c...       equidistributed
c
           temp = slope(i) * (xiold(i+1)-xiold(i))
           slphmx = dmax1(slphmx,temp)
           accum(i+1) = accum(i) + temp
           iflip = - iflip
  210  Continue
      avrg = accum(nold+1) / dfloat(nold)
      degequ = avrg / dmax1(slphmx,precis)
c
c...  naccum=expected n to achieve .1x user requested tolerances
c
      naccum = int(accum(nold+1) + 1.d0)
      if (iprint .lt. 0)  then
      CALL Rprintd1('Mesh info, degree of equidistribution = ',Degequ)
      CALL Rprinti1('Prediction for required N =  ', Naccum)
      endif 
c
c...  decide if mesh selection is worthwhile (otherwise, halve)
c
      if (avrg .lt. precis)                         go to 100
      if (degequ .ge. .5d0)                         go to 100
c
c...  nmx assures mesh has at least half as many subintervals as the
c...  previous mesh
c
      nmx = max0 (nold+1, naccum) / 2
c
c...  this assures that halving will be possible later (for error est)
c
      nmax2 = nmax / 2
c
c...  the mesh is at most halved
c
      n = min0 (nmax2, nold, nmx)
  220 noldp1 = nold + 1
      if (n .lt. nfxp1) n=nfxp1
      mshnum = mshnum + 1
c
c...  if the new mesh is smaller than the old mesh set mshnum
c...  so that the next call to  sysnewmsh  will produce a halved
c...  mesh. if n .eq. nold / 2 increment mshalt so there can not
c...  be an infinite loop alternating between n and n/2 points.
c
      if (n .lt. nold) mshnum = mshlmt
      if (n .gt. nold/2)  mshalt = 1
      if (n .eq. nold/2)  mshalt = mshalt + 1
      mshflg = 0
c
c...  having decided to generate a new mesh with n subintervals we now
c...  do so, taking into account that the nfxpnt points in the array
c...  fixpnt must be included in the new mesh.
c
      in = 1
      accl = 0.d0
      lold =2
      xi(1) = aleft
      xi(n+1) = aright
      do 310 i = 1, nfxp1
           if (i .eq. nfxp1)                        go to 250
           lnew=lold
           do 230 j = lold, noldp1
             lnew = j
             if (fixpnt(i) .le. xiold(j))           go to 240
  230      continue
  240      continue
           accr = accum(lnew) + (fixpnt(i)-xiold(lnew))*slope(lnew-1)
           nregn = int((accr-accl) / accum(noldp1) * dfloat(n) - .5d0)
           nregn = min0(nregn, n - in - nfxp1 + i)
           xi(in + nregn + 1) = fixpnt(i)
           go to 260
  250      accr = accum(noldp1)
           lnew = noldp1
           nregn = n - in
  260      if (nregn .eq. 0)                        go to 300
           temp = accl
           tsum = (accr - accl) / dfloat(nregn+1)
           do 290 j = 1, nregn
             in = in + 1
             temp = temp + tsum
             lcarry = lold
             do 270 l = lold, lnew
               lcarry = l
               if (temp .le. accum(l))              go to 280
  270        continue
  280        continue
             lold = lcarry
             xi(in) = xiold(lold-1) + (temp - accum(lold-1)) /
     1     slope(lold-1)
  290     Continue
  300      in = in + 1
           accl = accr
           lold = lnew
  310 continue
      mode = 1
  320 continue
c
c...  regardless of how the new mesh is chosen, the new collocation
c...  points xij in (aleft,aright) are generated here
c
      k2 = 1
      do 335 i = 1,n
           h = (xi(i+1) - xi(i)) / 2.d0
           xm = (xi(i+1) + xi(i)) / 2.d0
           do 330 j = 1,k
             xij(k2) = rho(j) * h + xm
             k2 = k2 + 1
  330 continue
  335 continue
      np1 = n + 1

      nalpha = n * k * ncomp + mstar
      return
c----------------------------------------------------------------
      end
      subroutine sysconsts
c
c***********************************************************************
c
c   purpose
c            assign (once) values to various array constants.
c
c   arrays assigned during compilation:
c     cnsts1 - weights for extrapolation error estimate
c     cnsts2 - weights for mesh selection
c              (the above weights come from the theoretical form for
c              the collocation error -- see (5))
c
c   arrays assigned during execution:
c     wgterr - the particular values of cnsts1 used for current run
c              (depending on k, m)
c     wgtmsh - gotten from the values of cnsts2 which in turn are
c              the constants in the theoretical expression for the
c              errors. the quantities in wgtmsh are 10x the values
c              in cnsts2 so that the mesh selection algorithm
c              is aiming for errors .1x as large as the user
c              requested tolerances.
c     jtol   - components of differential system to which tolerances
c              refer (viz, if ltol(i) refers to a derivative of u(j),
c              then jtol(i)=j)
c     root   - reciprocals of expected rates of convergence of compo-
c              nents of z(j) for which tolerances are specified
c     rho    - the k collocation points on (-1,1)
c     vncol  - the mesh independent b-splines values at collocation
c              points
c
c***********************************************************************
c
      implicit real(kind=8) (a-h,o-z)
      common /colloc/ rho(7),wgterr(40)
      common /order/ k,ncomp,mstar,kd,kdm,mnsum,m(20)
      common /bsplin/ vncol(66,7), vnsave(66,5), vn(66)
      common /errors/ tol(40),wgtmsh(40),tolin(40),root(40),
     1       jtol(40),ltol(40),ntol
      common /nonln/ precis,nonlin,iter,limit,icare,iprint,iguess,ifreez
      dimension cnsts1(28), cnsts2(28)
      data cnsts1 / .25d0,   .625d-1, .72169d-1,    1.8342d-2,
     1     1.9065d-2, 5.8190d-2,   5.4658d-3, 5.3370d-3, 1.8890d-2,
     2     2.7792d-2,   1.6095d-3, 1.4964d-3, 7.5938d-3, 5.7573d-3,
     3     1.8342d-2,  4.673d-3, 4.150d-4, 1.919d-3, 1.468d-3,
     4     6.371d-3, 4.610d-3,  1.342d-4, 1.138d-4, 4.889d-4,
     5     4.177d-4, 1.374d-3, 1.654d-3, 2.863d-3  /
      data cnsts2 / 1.25d-1,    2.604d-3,  8.019d-3,  2.170d-5,
     1     7.453d-5, 5.208d-4, 9.689d-8, 3.689d-7, 3.100d-6, 2.451d-5,
     2     2.691d-10, 1.120d-9,   1.076d-8,  9.405d-8,  1.033d-6,
     3     5.097d-13, 2.290d-12, 2.446d-11, 2.331d-10, 2.936d-9,
     4     3.593d-8,  7.001d-16, 3.363d-15, 3.921d-14, 4.028d-13,
     5     5.646d-12, 7.531d-11, 1.129d-9  /
c
c...  assign weights for error estimate
c
      koff = k * ( k + 1 ) / 2
      iextra = 1
      do 15 j = 1,ncomp
           lim = m(j)
           do 10 l = 1,lim
             wgterr(iextra) = cnsts1(koff - lim + l)
             iextra = iextra + 1
   10 continue
   15 continue
c
c...  assign array values for mesh selection: wgtmsh, jtol, and root
c
      jcomp = 1
      mtot = m(1)
      do 40 i=1,ntol
           ltoli = ltol(i)
   20      continue
           if (ltoli .le. mtot)                     go to 30
           jcomp = jcomp + 1
           mtot = mtot + m(jcomp)
           go to 20
   30      continue
           jtol(i) = jcomp
           wgtmsh(i) = 1.d1 * cnsts2(koff+ltoli-mtot) / tolin(i)
           root(i) = 1.d0 / dfloat(k+mtot-ltoli+1)
   40 continue
c
c...  specify collocation points
c
      if (k .eq. 1) then
        goto  50
      else if (k .eq. 2) then
        goto  60
      else if (k .eq. 3) then
        goto  70
      else if (k .eq. 4) then
        goto  80
      else if (k .eq. 5) then
        goto  90
      else if (k .eq. 6) then
        goto  100
      else if (k .eq. 7) then
        goto  110
      endif

C      go to (50,60,70,80,90,100,110),k
   50 rho(1) = 0.d0
      go to 120
   60 rho(2) = .57735026918962576451d0
      rho(1) = - rho(2)
      go to 120
   70 rho(3) = .77459666924148337704d0
      rho(2) = .0d0
      rho(1) = - rho(3)
      go to 120
   80 rho(1) = -.86113631159405257523d0
      rho(2) = -.33998104358485626480d0
      rho(3) = - rho(2)
      rho(4) = - rho(1)
      go to 120
   90 rho(5) = .90617984593866399280d0
      rho(4) = .53846931010568309104d0
      rho(3) = .0d0
      rho(2) = - rho(4)
      rho(1) = - rho(5)
      go to 120
  100 rho(6) = .93246951420315202781d0
      rho(5) = .66120938646626451366d0
      rho(4) = .23861918608319690863d0
      rho(3) = -rho(4)
      rho(2) = -rho(5)
      rho(1) = -rho(6)
      go to 120
  110 rho(7) = .949107991234275852452d0
      rho(6) = .74153118559939443986d0
      rho(5) = .40584515137739716690d0
      rho(4) = 0.d0
      rho(3) = -rho(5)
      rho(2) = -rho(6)
      rho(1) = -rho(7)
  120 continue
c
c...  put mesh independent b-splines values at collocation point
c...  rho(j) into vncol(*,j), j=1,...,k.
c
      do 130 j=1,k
           arg = .5d0 * (1.d0 - rho(j))
           call sysbspfix (arg, vncol(1,j), k, ncomp, m)
  130 continue
c
c...  put mesh independent b-splines values at the points in unit in-
c...  terval 0, 1/6, 1/3, 2/3, 5/6  into vnsave.  these values are to
c...  be used in   sysnewmsh   and   syserrchk .
c
      call sysbspfix (1.d0, vnsave(1,1), k, ncomp, m)
      call sysbspfix (5.d0/6.d0, vnsave(1,2), k, ncomp, m)
      call sysbspfix (2.d0/3.d0, vnsave(1,3), k, ncomp, m)
      call sysbspfix (1.d0/3.d0, vnsave(1,4), k, ncomp, m)
      call sysbspfix (1.d0/6.d0, vnsave(1,5), k, ncomp, m)
      return
      end
      subroutine syserrchk(imesh,xiold,aldif,valstr,mstar,ifin)
c
c***********************************************************************
c
c      purpose
c               determine the error estimates and test to see if the
c               error tolerances are satisfied.
c
c      variables
c        xiold  - current mesh points
c        valstr - values of the previous solution which are needed
c                 for the extrapolation- like error estimate.
c        wgterr - weights used in the extrapolation-like error
c                 estimate. the array values are assigned in
c                 subroutine  sysconsts.
c        errest - storage for error estimates
c        err    - temporary storage used for error estimates
c Francesca eliminated unused argment
c        work   - space to be used to store values of z at the
c                 mesh points for printout. its dimension is
c                 mstar * nmax.
c        z      - approximate solution on mesh xi
c        ifin   - a 0-1 variable. if imesh = 2 then on return it
c                 indicates whether the error tolerances were satisfied.
c        imesh  = 1  the current mesh resulted from mesh selection
c                    or is the initial mesh.
c               = 2  the current mesh resulted from doubling the
c                    previous mesh
c        mshflg - is set by syserrchk to indicate to sysnewmsh whether
c                 any values of the current solution are stored in
c                 the array valstr. (0 for no, 1 for yes)
c
c**********************************************************************
c
      implicit real(kind=8) (a-h,o-z)
      dimension err(40),z(40),errest(40)
      common /order/k,ncomp,mstr,kd,kdm,mnsum,m(20)
      common /appr/ n,nold,nmax,nalpha,mshflg,mshnum,mshlmt,mshalt
      common /side/  zeta(40), aleft, aright, izeta
      common /errors/ tol(40),wgtmsh(40),tolin(40),root(40),
     1       jtol(40),ltol(40),ntol
      common /colloc/ rho(7),wgterr(40)
      common /nonln/ precis,nonlin,iter,limit,icare,iprint,iguess,ifreez
      common /bsplin/ vncol(66,7), vnsave(66,5), vn(66)
C Karline: added dumm(1)      
      dimension xiold(*), aldif(*), valstr(*),dumm(1)
c
      ifin = 1
      noldp1 = nold + 1
c
c...  if full output has been requested, print values of the
c...  solution components   z  at the meshpoints.
c
c karline: toggled off printing

      if (imesh.eq.1)                               return
c
c...  imesh = 2 so error estimates are to be generated and tested
c...  to see if the tolerance requirements are satisfied.
c
      do 40 j = 1,mstar
        errest(j) = 0.d0
   40 Continue
      do 100 iback = 1,nold
           i = nold +1 -iback
c
c...       the error estimates are obtained by combining values of
c...       the numerical solutions for two meshes.
c...       for each value of iback we will consider the two
c...       approximations at 2 points in each of
c...       the new subintervals.  we work backwards through
c...       the subinterval so that new values can be stored
c...       in valstr in case they prove to be needed later
c...       for an error estimate. the routine  sysnewmsh
c...       filled in the needed values of the old solution
c...       in valstr.
c
           mshflg = 1
           do 50 j = 1,mstar
             z(j) = 0.d0
             err(j) = 0.d0
   50      Continue
           do 65 j = 1,2
             jj = 5 - j
             knew = ( 4 * (i-1) + 3 - j ) * mstar + 1
             kstore = ( 2 * (i-1) + 2 - j ) * mstar + 1
             x = xiold(i) + dfloat(3-j)/3.d0*(xiold(i+1)-xiold(i))
             call sysapprox (i, x, valstr(knew), vnsave(1,jj), xiold,
     1            nold, aldif, k, ncomp, m, mstar, 3,dumm,0)
             do 60 l = 1,mstar
               err(l) = err(l) + wgterr(l) * dabs(valstr(knew) -
     1         valstr(kstore))
               z(l) = z(l) + .5d0 * dabs(valstr(knew))
               knew = knew + 1
               kstore = kstore + 1
   60      continue
   65      continue
c
c...       test whether the tolerance requirements are satisfied
c...       in the i-th interval.
c
           if (ifin .eq. 0)                         go to 80
           do 70 j = 1, ntol
             ltolj = ltol(j)
            if ( err(ltolj) .gt. tolin(j) * (z(ltolj)+1.d0) ) ifin = 0
   70      Continue
   80      do 90 l = 1,mstar
             errest(l) = dmax1(errest(l),err(l))
   90      Continue
  100 continue

C karline: toggled off lot of printing
      return
c--------------------------------------------------------------
      end
c
c-----------------------------------------------------------------------
c                    p a r t  3
c          collocation system setup routines -- see (1)
c-----------------------------------------------------------------------
c
      subroutine syslsyslv (iflag, xi, xiold, xij, alpha, aldif,
     1           rhs, alpho, a, ipiv, integs, rnorm,
     2           mode, fsub, dfsub, gsub, dgsub, guess,rpar,ipar)
c*********************************************************************
c
c   purpose
c         this routine controls the set up and solution of a linear
c      system of collocation equations.
c         the matrix  a  is cast into an almost block diagonal
c      form by an appropriate ordering of the columns and solved
c      using the package of de boor-weiss (4). the matrix is composed
c      of n blocks. the i-th block has the size
c                  integs(1,i) * integs(2,i).
c      it contains in its last rows the linearized collocation equa-
c      tions (both bundary conditions and differential equations )
c      corresponding to the i-th subinterval.  integs(3,i)  steps of
c      gaussian elimination are applied to it to achieve a  plu
c      decomposition.  the right hand side vector is put into  rhs
c      and the solution vector is returned in  alpha.
c
c         syslsyslv operates according to one of 5 modes:
c      mode = 0 - set up both  a  and  rhs , and solve system
c                 (for linear problems).
c      mode = 1 - set up  rhs  only and compute its norm.
c      mode = 2 - set up  a  only and solve system.
c      mode = 3 - perform forward and backward substitution (do not set
c                 up a nor form the rhs).
c
c         for the first iteration on a particular mesh,
c      integs  is computed.  also, the initial  alpha  on
c      the new mesh is computed.
c
c   variables
c
c      irhs,ia,izeta,ialpho - pointers to rhs,a,zeta,alpho respectively
c                             (necessary to keep track of blocks of a
c                             during matrix manipulations)
c      alpho - b-spline coeffs for previous solution
c      dg    - partial derivatives of g from dgsub
c      df    - partial derivatives of f from dfsub
c      rnorm - euclidean norm of rhs
c      lside - number of side conditions in current and previous blocks
c      icolc - pointer to current collocation point array xij
c      id    - (another) pointer for rhs
c      iguess = 1 when current soln is user specified
c             = 0 otherwise
c
c*********************************************************************
      implicit real(kind=8) (a-h,o-z)
      common /order/ k, ncomp, mstar, kd, kdm, mnsum, m(20)
      common /side/  zeta(40), aleft, aright, izeta
      common /bsplin/  vncol(66,7), vnsave(66,5), vn(66)
      common /appr/ n,nold,nmax,nalpha,mshflg,mshnum,mshlmt,mshalt
      common /nonln/ precis,nonlin,iter,limit,icare,iprint,iguess,ifreez
      common /hi/   dn1, dn2, dn3
      external dfsub, dgsub
      dimension  alpho(*), xi(*), xiold(*), xij(*), alpha(*)
C Karline: added dimension dummy(1)   
      dimension aldif(*), rhs(*), a(*), ipiv(*), integs(3,*),dummy(1)
      dimension z(40), f(40), df(800), dmval(20)
      dimension rpar(*), ipar(*)

      integer nfunc, njac, nstep, nbound, njacbound
      common/coldiag/nfunc, njac, nstep, nbound, njacbound
c
      m1 = mode + 1
      if (m1 .eq. 1) then
        goto  10
      else if (m1 .eq. 2 .or. m1 .eq. 3) then
        goto  30
      else if (m1 .eq. 4) then
        goto  310
      endif
      
C      go to (10, 30, 30, 310), m1
c
c...  linear problem initialization
c
   10 do 20 i=1,mstar
        z(i) = 0.d0
   20 Continue
c
c...  initialization
c
   30 irhs = 0
      ia = 1
      izeta = 1
      lside = 0
      rnorm = 0.d0
      ialpho = 0
      if (iter .ge. 1 .or. mode .eq. 2)             go to 80
c
c...  build integs (describing block structure of matrix)
c
      do 70 i = 1,n
           integs(2,i) = kdm
           if (i .lt. n)                            go to 40
           integs(3,i) = kdm
           lside = mstar
           go to 60
   40      integs(3,i) = kd
   50      if( lside .eq. mstar )                   go to 60
           if ( zeta(lside+1) .ge. xi(i+1) )        go to 60
           lside = lside + 1
           go to 50
   60      nrow = kd + lside
           integs(1,i) = nrow
   70  Continue
c
c...  the do loop 290 sets up the linear system of equations.
c
   80 CONTINUE
C karline
      nstep = nstep + 1

      do 290 i=1,n
           xil = xi(1)
           if (i .gt. 1) xil = xi(i-1)
           xir = xi(n+1)
           if (i .lt. n) xir = xi(i+2)
           dn1 = 1.d0 / (xi(i+1) - xil)
           dn2 = 1.d0 / (xi(i+1) - xi(i))
           dn3 = 1.d0 / (xir - xi(i))
c
c...       construct a block of  a  and a corresponding piece of  rhs.
c
           nrow = integs(1,i)
           ii = i
           icolc = (i-1) * k
           id = irhs + izeta - 1
c
c...       go thru the ncomp collocation equations and side conditions
c...       in the i-th subinterval
c
           do 270  ll=1,k
             xx = xij (icolc+ll)
  100        if ( izeta .gt. mstar )                go to 160
             if ( zeta(izeta) .ge. xx )             go to 160
c
c...         build equation for a side condition.
c
  110        id = id + 1
             ialpho = ialpho + 1
             if (mode .eq. 0)                       go to 130
             if (iguess .ne. 1)                     go to 120
c
c...         case where user provided current approximation
c
             call guess (zeta(izeta), z, dmval, rpar, ipar)
      if (mode .eq. 1) then
        goto  130
      else if (mode .eq. 2) then
        goto  140
      endif
C             go to (130, 140), mode
c
c...         other nonlinear case
c
  120        call sysapprox (ii, zeta(izeta), z, vn, xiold, nold,
     1            aldif, k, ncomp, m, mstar, 1, dummy, 0)
             if (mode .eq. 2)                       go to 140
c
c...         find  rhs  boundary value.
c
  130        call gsub (izeta, mstar, z, g,rpar,ipar)
c karline
             nbound = nbound + 1
             rhs(id) = -g
             rnorm = rnorm + g**2
             if (mode .eq. 1)                       go to 150
c
c...         build a row of  a  corresponding to a boundary point
c
  140        call sysbldblk (i, zeta(izeta), ll, a(ia), nrow, id-irhs, 
     1       z, df, ncomp, xi, alpho, ialpho, 1, dfsub, dgsub,rpar,ipar)
  150        izeta = izeta + 1
c
c...         check for other side conditions.
c
             if (izeta .gt. mstar .and.
     1           zeta(mstar) .ge. dmin1(xx,aright)) go to 280
             if (xx .gt. xi(n+1))                   go to 260
             go to 100
c
c...         this value corresponds to a collocation (interior)
c...         point. build the corresponding  ncomp  equations.
c
  160        if (iguess .ne. 1) THEN
             if (m1 .eq. 1) then
               goto  210
             else if (m1 .eq. 2) then
               goto  170
             else if (m1 .eq. 3) then
              goto  230
             endif

C  160        if (iguess .ne. 1)  go to (210, 170, 230), m1              
             ENDIF
c
c...         use initial approximation provided by the user.
c
             call guess(xx, z, dmval, rpar, ipar)
      if (mode .eq. 1) then
        goto  190
      else if (mode .eq. 2) then
        goto  250
      endif
C             go to (190, 250), mode
c
c...         find  rhs  values
c
  170        if (iter .ge. 1 )                      go to 180
             call sysapprox (ii, xx, z, vn, xiold, nold, aldif, k,
     1            ncomp, m, mstar, 1, dmval, 1)
             go to 190
  180        call sysapprox (i, xx, z, vncol(1,ll), xiold, nold,
     1            aldif, k, ncomp, m, mstar, 3, dmval, 1)
  190        continue
             call fsub (mstar, xx, z, f,rpar,ipar)
c karline
             nfunc = nfunc +1
c
c...         fill in  rhs  values (and accumulate its norm).
c
             do 200 j=1,ncomp
               id = id + 1
               value = dmval(j) - f(j)
               rhs(id) = -value
               rnorm = rnorm + value**2
               if (iter .ge. 1)                     go to 200
               ialpho = ialpho + 1
               alpho(ialpho) = dmval(j)
  200        continue
             go to 260
c
c...         the linear case
c
  210        continue
             call fsub (mstar, xx, z, f,rpar,ipar)
             nfunc = nfunc +1

             do 220 j=1,ncomp
               id = id + 1
               rhs(id) = f(j)
  220        Continue
             id = id - ncomp
             go to 250
c
c...         evaluate former collocation soln for mode=2
c
  230        if (iter .ge. 1 )                      go to 240
             call sysapprox (ii, xx, z, vn, xiold, nold, aldif, k,
     1            ncomp, m, mstar, 1, dummy, 0)
             go to 250
  240        call sysapprox (i, xx, z, vncol(1,ll), xiold, nold,
     1            aldif, k, ncomp, m, mstar, 3, dummy, 0)
c
c...         fill in ncomp rows of  a
c
  250        call sysbldblk (i, xx, ll, a(ia), nrow, id-irhs+1, z,
     1       df, ncomp, xi, alpho, ialpho, 2, dfsub, dgsub,rpar,ipar)
             id = id + ncomp
c
c...         prepare to set up side conditions for last subinterval
c
  260        if (ll .lt. k)                         go to 270
             if (i .lt. n .or. izeta .gt. mstar)    go to 280
             xx = xi(n+1) + 1.d0
             go to 110
  270      continue
c
c...       update counters -- i-th block completed
c
  280      irhs = irhs + nrow
           ia = ia + nrow * kdm
  290 continue
      if (mode .ne. 1)                              go to 300
      rnorm = dsqrt(rnorm / dfloat(nalpha))
      return
c
c...  solve the linear system.
c
c...  matrix decomposition
c
  300 call sysfcblok (a, integs, n, ipiv, alpha, iflag)
c
c...  check for singular matrix
c
      if(iflag .eq. 0)                              return
c
c...  perform forward and backward substitution for mode=0,2, or 3.
c
  310 call syssbblok (a, integs, n, ipiv, rhs, alpha)
c
c...  find the coefficients alpha of the initial sysapprox if necessary.
c
      if (iter .ge. 1 .or. mode .ne. 2)             return
      ialpho = 0
      irhs = 0
      isto = 0
      do 325 i=1,n
           nrow = integs(1,i)
           irhs = irhs + isto
           istart = isto + 1
           isto = nrow - kd
           do 320 j=istart,nrow
             irhs = irhs + 1
             ialpho = ialpho + 1
             rhs(irhs) = rhs(irhs) + alpho(ialpho)
  320 continue
  325 Continue
      call syssbblok (a, integs, n, ipiv, rhs, alpho)
      do 330 i=1,nalpha
        alpho(i) = alpho(i) - alpha(i)
  330 Continue
      return
      end
      subroutine sysbldblk (i, x, ll, q, nrow, nc, z, df, ncomp,
     1           xi, alpho, ialpho, mode, dfsub, dgsub,rpar,ipar)
c
c**********************************************************************
c
c   purpose:
c
c      construct collocation matrix rows according to mode:
c      mode = 1  -  a row corresponding to a side condition.
c      mode = 2  -  a group of   ncomp   rows corresponding
c                   an interior collocation point.
c
c   variables:
c
c      alpho  - used only on the first iteration for nonlinear
c               problems when the first approximation is other
c               than a b-spline representation on the current mesh.
c               a right hand side is being built up in  alpho which,
c               when the inverted collocation matrix is applied to it,
c               will produce a first approximation on the current mesh
c               in terms of b- splines so the step-length algorithm
c               in  syscontrl  can operate.
c      x      - the collocation or side condition point.
c      i      - the subinterval containing x
c      ll     - if x is a collocation point then it is the ll-th
c               of k collocation points on the i-th subinterval.
c      q      - the sub-block of the collocation matrix in
c               which the equations are to be formed.
c      nrow   - no. of rows in q.
c      nc     - the first row in q to be used for equations.
c      z      - z(x)
c      dg     - the derivatives of the side condition.
c      df     - the jacobian at x.
c      id     - the row of q being constructed.
c      basef  - values and derivatives of the b-spline basis
c               for each of the components.
c      jcomp  - counter for the component being dealt with.
c      l      - counter for the b-splines representing u(jcomp).
c      j      - counter for the lowest m(jcomp) derivatives of
c               bsplines representing u     .
c                                      jcomp
c
c**********************************************************************
      implicit real(kind=8) (a-h,o-z)
      common /appr/ n,nold,nmax,nalpha,mshflg,mshnum,mshlmt,mshalt
      common /order/  k, nd, mstar, kd, kdm, mnsum, m(20)
      common /side/   zeta(40), aleft, aright, izeta
      common /nonln/ precis,nonlin,iter,limit,icare,iprint,iguess,ifreez
      common /bsplin/  vncol(66,7), vnsave(66,5), vn(66)
      dimension q(nrow,*), z(*), df(ncomp,*),rpar(*),ipar(*)
      dimension xi(*), basef(620), alpho(*), dg(40)
      integer nfunc, njac, nstep, nbound, njacbound
      common/coldiag/nfunc, njac, nstep, nbound, njacbound

c
      nk = nc
      if (mode .eq. 2)  nk = nc + ncomp - 1
      do 15 j=nc,nk
           do 10 l=1,kdm
             q(j,l)=0.d0
   10      Continue
   15 Continue
c
c...  branch according to  m o d e
c
      IF (MODE .EQ. 1) THEN
        GOTO  20
      ELSE IF (MODE .EQ. 2) THEN
        GOTO  130
      ENDIF
c      go to (20, 130), mode
c
c...  x is a boundary point
c
   20 call sysbspder (vn, xi, n, x, i, basef, 2)
c
c...  provide coefficients of the j-th linearized side condition.
c...  specifically, at x=zeta(j) the j-th side condition reads
c...  dg(1)*z(1) + ... +dg(mstar)*z(mstar) + g = 0
c
      call dgsub (izeta, mstar, z, dg, rpar,ipar)
      njacbound = njacbound + 1

      if (iter .ge. 1 .or. nonlin .eq. 0)           go to 40
      value = 0.d0
      do 30 j=1,mstar
         value = value + dg(j) * z(j)
   30 Continue
      alpho(ialpho) =  value
   40 iq = 0
      iqm = mstar
      idg = 0
      ibasef = 0
      id = nc
c
      do 120 jcomp=1,ncomp
           mj = m(jcomp)
           mj1 = mj + 1
           kmj = k - mj
c
c...       incorporate the values and derivatives for
c...       the b-splines which are nonzero on the preceeding
c...       subinterval.
c
           do 60 l=1,mj
             do 50 j=1,mj
               q(id, iq+l) = q(id, iq+l) + dg(idg+j) *
     1       basef(ibasef + j)
   50        Continue
            ibasef = ibasef + mj1
   60      Continue
c
c...       the b-splines which are nonzero on the current
c...       subinterval only.
c
           if (kmj .le. 0)                          go to 90
           do 80 l=1,kmj
             do 70 j=1,mj
                q(id, iqm+l) = q(id, iqm+l) + dg(idg+j) *
     1       basef(ibasef+j)
   70        Continue
             ibasef = ibasef + mj1
   80      Continue
c
c...       the b-splines which are nonzero on the succeeding
c...       subinterval as well.
c
   90      do 110 l=1,mj
             do 100 j=1,mj
               q(id, iq+kd+l) = q(id, iq+kd+l) + dg(idg+j) *
     1       basef(ibasef+j)
  100        Continue
             ibasef = ibasef + mj1
  110      Continue
c
           idg =idg + mj
           iq = iq + mj
           iqm = iqm + kmj
  120 continue
      return
c
c...  build ncomp rows for interior collocation point x.
c...  the linear expressions to be constructed are:
c...   (m(jj))
c...  u     -  df(jj,1)*z(1) - ... - df(jj,mstar)*z(mstar)
c...   jj
c...  for jj = 1 to ncomp.
c
  130 call sysbspder (vncol(1,ll), xi, n, x , i, basef, 3)
      call dfsub (mstar,x, z, df,rpar,ipar)
      njac = njac + 1

c
c...  loop over the  ncomp  expressions to be set up for the
c...  current collocation point.
c
      do 240 jj=1,ncomp
           if (iter .ge. 1 .or. nonlin .eq. 0)      go to 150
           ialpho = ialpho + 1
           value = 0.d0
           do 140 j=1,mstar
             value = value + df(jj,j) * z(j)
  140      Continue
           alpho(ialpho) = alpho(ialpho) - value
  150      id = jj + nc - 1
           iq=0
           iqm=mstar
           idf=0
           ibasef=0
c
c...       note that if jj .eq. jcomp an entry has to be made for the
c...       m(jcomp)-th derivative of the jcomp-th component.
c
           do 230 jcomp=1,ncomp
             mj = m(jcomp)
             mj1 = mj + 1
             kmj = k - mj
c
c...         use the b-splines which are nonzero on the preceeding
c...         subinterval.
c
             do 170 l=1,mj
               if (jcomp . eq. jj) q(id, iq+l) = basef(ibasef+mj1)
               do 160 j=1,mj
                 q(id, iq+l) = q(id, iq+l) - df(jj, idf+j) *
     1         basef(ibasef+j)
  160          Continue
               ibasef = ibasef + mj1
  170        Continue
c
c...         the b-splines which are nonzero on the current
c...         subinterval only.
c
             if (kmj .le. 0)                        go to 200
             do 190 l=1,kmj
               if (jcomp .eq. jj) q(id, iqm+l) = basef(ibasef+mj1)
               do 180 j=1,mj
                 q(id, iqm+l) = q(id, iqm+l) -
     1         df(jj, idf+j) * basef(ibasef+j)
  180          Continue
               ibasef = ibasef + mj1
  190        Continue
c
c...         the b-splines which are nonzero on the succeeding
c...         subinterval as well.
c
  200        do 220 l=1,mj
               if (jcomp .eq. jj) q(id, iq+kd+l) = basef(ibasef+mj1)
               do 210 j=1,mj
                 q(id, iq+kd+l) = q(id, iq+kd+l) - df(jj, idf+j) *
     1         basef(ibasef+j)
  210        Continue
             ibasef = ibasef + mj1
  220    Continue
c
             idf = idf + mj
             iq = iq + mj
             iqm = iqm + kmj
  230      continue
  240 continue
      return
      end
c
c-----------------------------------------------------------------------
c                    p a r t  4
c          b-spline routines -- see (3)
c-----------------------------------------------------------------------
c
      SUBROUTINE sysappsln (x, z, fspace, ispace)
c
c*****************************************************************
c
c     purpose
c
c           set up a standard call to  sysapprox  to evaluate the
c           approximate solution  z = z( u(x) )  at a point x
c           (it has been computed by a call to  colsys ).
c           the parameters needed for  sysapprox  are retrieved
c           from the work arrays  ispace  and  fspace .
c
c*****************************************************************
c
      implicit real(kind=8) (a-h,o-z)
C Karline: added dumm(1)      
      dimension z(*), fspace(*), ispace(*), dumm(1)
      is6 = ispace(6) + 1
      is5 = ispace(1) + 2
      i = 1
      call sysapprox (i, x, z, fspace(is6), fspace, ispace(1),
     1     fspace(is5), ispace(2), ispace(3), ispace(8),
     2     ispace(4), 1, dumm, 0)
      return
      end
      subroutine sysapprox (i, x, z, vn, xi, n, aldif, k, ncomp,
     1           m, mstar, mode, dmval, modhi)
c
c***********************************************************************
c
c   purpose
c                                    (1)       (m1-1)     (mncomp-1)
c           evaluate z(u(x))=(u (x),u (x),...,u  (x),...,u  (x)      )
c                              1     1         1          mncomp
c           at one point x.
c           if modhi=1, evaluate  mj-th  derivatives too.
c
c   variables
c     vn     - triangular array of b-spline values filled in by
c              routines sysbspfix and sysbspvar
c     xi     - the current mesh (having n subintervals)
c     aldif  - the array of divided differences of the current
c              solution vectors coefficients alpha (and previously
c              determined in the routine sysappdif)
c     mode   - determines the amount of initialization needed
c            = 5  forms z(u(x)) using aldif and vn
c            = 3  as in =5, but finishes filling in vn using sysbspvar
c            = 2  as in =3, but starts filling in vn using sysbspfix
c            = 1  as in =2, but determines i such that
c                           xi(i) .le. x .lt. xi(i+1) (unless x=xi(n+1))
c            = 4  a special case which only determines i as above
c     dmval  - array of mj-th derivatives of the solution components
c              uj (evaluated if modhi=1)
c
c***********************************************************************
c
      implicit real(kind=8) (a-h,o-z)
      common /nonln/ precis,nonlin,iter,limit,icare,iprint,iguess,ifreez
      common /side/  zeta(40), aleft, aright, izeta
      dimension z(*), vn(*), xi(*), aldif(*), m(*), dmval(*)
c
      dnk2  = 0.0d0
      ivnhi = 1
      IF (MODE .EQ. 1) THEN
        GOTO  10
      ELSE IF (MODE .EQ. 2) THEN
        GOTO  60
      ELSE IF (MODE .EQ. 3) THEN
        GOTO  70
      ELSE IF (MODE .EQ. 4) THEN
        GOTO  10
      ELSE IF (MODE .EQ. 5) THEN
        GOTO  80
      ENDIF
C      go to (10, 60, 70, 10, 80), mode
c
c...  mode = 1 or 4,  locate i so  xi(i) .le. x .lt. xi(i+1)
c
   10 continue
      if (x .ge. xi(1)-precis .and. x .le. xi(n+1)+precis)
     1                                              go to 20
      if (iprint .lt. 1) then
      CALL Rprintd3('Domain error in Approx, X, Aleft, Aright ',
     +  X, Xi(1), Xi(N+1))
      endif 

      if (x .lt. xi(1)) x = xi(1)
      if (x .gt. xi(n+1)) x = xi(n+1)
   20 if (i .gt. n .or. i .lt. 1) i = (n+1) / 2
      ileft = i
      if (x .lt. xi(ileft))                         go to 40
      do 30 l=ileft,n
           i = l
           if (x .lt. xi(l+1))                      go to 60
   30 continue
      go to 60
   40 iright = ileft - 1
      do 50 l=1,iright
           i = iright + 1 - l
           if (x .ge. xi(i))                        go to 60
   50 continue
   60 if (mode .eq. 4)                              return
c
c...  mode = 1 or 2   begin filling in vn using sysbspfix .
c...  compute mesh independent splines.
c
      rhox = (xi(i+1) - x) / (xi(i+1) - xi(i))
      call sysbspfix (rhox, vn, k, ncomp, m)
c
c...  mode = 1, 2, or 3   finish filling in vn using sysbspvar
c
   70 call sysbspvar (i, x, vn, xi, n, k, ncomp, m)
c
c...  mode .ne. 4  determine z(u(x))
c
   80 do 90 l=1,mstar
        z(l) = 0.d0
   90 Continue
      indif = 0
      k5 = 1
      if (modhi .eq. 0)                             go to 110
c
c...  initialize for subsequent evaluation of  mj-th  derivatives.
c
      ivnhi = k * (k-1) / 2
      dnk2 = dfloat(k) / (xi(i+1) - xi(i))
      incomp = 0
      do 100 j=1,ncomp
        dmval(j) = 0.d0
  100 Continue
c
c...  evaluate  z( u(x) ).
c
  110 do 150 j = 1, ncomp
           mj = m(j)
           nalphj = n * k + mj
           kmr = k + mj
           ivn = kmr * (kmr - 1) / 2
           do 130 nr = 1, mj
             left = i * k + mj - kmr
             do 120 l = 1, kmr
               leftpl = left + l
               z(k5) = z(k5) + aldif(indif+leftpl) * vn(ivn+l)
  120        Continue
             kmr = kmr - 1
             ivn = ivn - kmr
             k5 = k5 + 1
             indif = indif + nalphj
  130      Continue
           if (modhi .eq. 0)                        go to 150
c
c...       evaluate  dmval(j) = mj-th derivative of uj.
c
           incomp = incomp + (mj-1) * nalphj
           left = (i-1) * k + mj - 1
           do 140 l = 1, k
              dmval(j) = dmval(j) + dnk2 * (aldif(incomp+left+l+1) -
     1     aldif(incomp+left+l)) * vn(ivnhi+l)
  140      Continue
           incomp = incomp + nalphj
  150 continue
      return
c--------------------------------------------------------
      end
      subroutine sysbspfix (rhox, vn, k, ncomp, m)
c
c**********************************************************************
c
c   purpose
c             evaluate the mesh independent bsplines at one point
c
c
c   variables
c     vn     - triangular array of b-spline values at x for orders
c              1 to k+m(ncomp) where xi(i) .le. x .le. xi(i+1) , column
c              j has length j and contains the j-th order b-spline
c              values and begins in location i + j*(j-1)/2.  values
c              not computed here are computed in sysbspvar.
c     rhox   = (xi(i+1)-x)/(xi(i+1)-xi(i))
c
c***********************************************************************
c
      implicit real(kind=8) (a-h,o-z)
      dimension vn(*), m(*)
      xrho = 1.d0 - rhox
      ivn = 0
c
c...  compute first group of mesh independent b-spline values
c
      vn(1) = 1.d0
      do 20 l=1,k
           ivn = ivn + l
           vnp = 0.d0
           do 10 j=1,l
             rep =  vn(ivn-l+j)
             vn(ivn+j) = vnp + rep * rhox
             vnp = rep * xrho
   10      Continue
          vn(ivn+l+1) = vnp
   20  Continue
c
c...  compute second group of mesh independent b-spline values
c
      md1 = m(ncomp) - 1
      if (md1 .le. 0)                               return
      do 40 l=1,md1
           ivn = ivn + k + l
           inc = l + 2
           vnp = vn(ivn+1-k) * xrho
           if (k .lt. inc)                          return
           do 30 j=inc,k
             rep = vn(ivn-l-k+j)
             vn(ivn+j) = vnp + rep  * rhox
             vnp = rep * xrho
   30      Continue
         vn(ivn+k+1) = vnp
   40 Continue
      return
      end
      subroutine sysbspvar (i, x, vn, xi, n, k, ncomp, m)
c
c***********************************************************************
c
c   purpose
c            evaluate the mesh dependent b-splines at one point x
c
c   variables
c     vn    - triangular array of values of b-splines of orders 1
c             to k+m(ncomp)  (described in sysbspfix)
c     x     - satisfies xi(i) .le. x .le. xi(i+1)
c
c**********************************************************************
c
      implicit real(kind=8) (a-h,o-z)
      dimension vn(*), xi(*), m(*)
      md1 = m(ncomp) -1
      if(md1 .le. 0)                                return
      xil = xi(1)
      if (i .gt. 1) xil = xi(i-1)
      xir = xi(n+1)
      if (i .lt. n) xir = xi(i+2)
      rho1 = (xi(i+1) - x) / (xi(i+1) - xi(i))
      rho2 = (xi(i+1) - x) / (xi(i+1) - xil)
      rho3 = (xir - x) / (xir - xi(i))
      xrho1 = 1.d0 - rho1
      xrho2 = 1.d0 - rho2
      xrho3 = 1.d0 - rho3
      ivn = k * (k+1) / 2
c
c...  recursively compute b-spline values.
c
      do 30 l=1,md1
           ivn = ivn + k + l
           vnp = 0.d0
           do 10 j=1,l
             rep = vn(ivn-l-k+j)
             vn(ivn+j) = vnp + rep  * rho2
             vnp = rep * xrho2
   10      Continue
           vn(ivn+l+1) = vnp + rho1 * vn(ivn-k+1)
           vnp = vn(ivn-l) * xrho1
           do 20 j=1,l
             rep = vn(ivn+j-l)
             vn(ivn+k+j) = vnp + rep * rho3
             vnp = rep * xrho3
   20      Continue
           vn(ivn+k+l+1) = vnp
   30   Continue
      return
      end
      subroutine sysbspder (vn, xmesh, n, x, i, basef, mode)
c
c***********************************************************************
c
c   purpose
c           evaluate the derivatives of the b-splines of appropriate
c           orders at one point x (used to set up the
c           collocation equations.)
c
c   variables
c
c     vn     - the triangular array of b-spline values calculated in
c              sysbspfix and sysbspvar
c     basef - b-spline derivatives needed to set up collocation
c             equations, viz, derivatives of orders 0,1,...,mj of
c             b-splines of order k+mj  (j=1,...,ncomp).  these
c             values are found using vn, alphd, and alphn (see below).
c     alphd  - array of divided differences corresponding to deriva-
c              tives of b-splines of order k+mncomp
c     alphn  - same as alphd, but for other order b-splines
c     alphdo - divided differences of one lower order, used to deter-
c              mine alphd
c     alphno - divided differences of one lower order, used to deter-
c              mine alphn
c      nd    - the no. of differential equations of distinct orders
c              (so no. of other differential equations =neq =ncomp-nd)
c     mnd    - the distinct orders of these nd differential equations
c     xmesh  - current mesh, with xmesh(i) .le. x .lt. xmesh(i+1)  (unle
c              x=xmesh(n+1)
c     mode   - determines the amount of initialization needed
c            = 4  compute the array basef
c            = 3  as in =4, but fill in subinterval dependent values
c                 of vn using sysbspvar
c            = 2  as in =3, but fill in subinterval independent values
c                 of vn using sysbspfix
c            = 1  as in =2, but calculate certain subinterval depen-
c                 dent constants
c
c***********************************************************************
c
      implicit real(kind=8) (a-h,o-z)
      common /order/ k, ncomp, mstar, kd, kdm, mnsum, m(20)
      common /hi/    dn1, dn2, dn3
      common /eqord/ ind(5), ineq(20), mnd(5), nd, neq
      dimension basef(*), vn(*), xmesh(*)
      dimension alphd(80), alphdo(80), alphn(280) , alphno(280)
c
      IF (MODE .EQ. 1) THEN
        GOTO  10
      ELSE IF (MODE .EQ. 2) THEN
        GOTO  20
      ELSE IF (MODE .EQ. 3) THEN
        GOTO  30
      ELSE IF (MODE .EQ. 4) THEN
        GOTO  40
      ENDIF
C      go to (10, 20, 30, 40), mode
c
c...  mode = 1   compute subinterval dependent constants
c
   10 xil = xmesh(1)                                                                  
      if (i .gt. 1)  xil = xmesh(i-1)
      xir = xmesh(n+1)
      if (i .lt. n)  xir = xmesh(i+2)
      dn1 = 1.d0 / (xmesh(i+1) - xil)
      dn2 = 1.d0 / (xmesh(i+1) - xmesh(i))
      dn3 = 1.d0 / (xir - xmesh(i))
c
c...  mode = 2   compute subinterval independent b-splines
c
   20 rhox = (xmesh(i+1) - x) * dn2
      call sysbspfix (rhox, vn, k, ncomp, m)
c
c...  mode = 3   compute subinterval dependent b-splines
c
   30 call sysbspvar (i, x, vn, xmesh, n, k, ncomp, m)
c
c...  mode = 4
c
   40 md = mnd(nd)
      kmd = k + md
      kmd1 = kmd + 1
      md1 = md + 1
      md2m2 = md * 2 - 2
      md2m1 = md2m2 + 1
      inl = kmd * 2
c
c...  initialize arrays alphdo and alphno
c
      do 50 j=1,kmd
           alphdo (j)  = 0.d0
           alphdo(j+kmd) = 1.d0
   50 Continue
      kup = kmd * md
      do 60 j=1,kup
        alphdo(j+inl) = 0.d0
   60 Continue
      ndm1 = nd - 1
      nrest = md2m2  - k
      inn = 0
      if (nrest .le. 0)                             go to 100
      if (nd .eq. 1)                                go to 100
      inl = 2 * md2m2
      do 90 nn = 1,ndm1
           mn2 = mnd(nn) + 2
           do 70 j = 1,md2m2
             alphno(j+inn) = 0.d0
             alphno(j+inn+md2m2) = 1.d0
   70      Continue
           kup = md2m2 * mnd(nn)
           do 80 j=1,kup
             alphno(j+inn+inl) = 0.d0
   80      Continue
          inn = inn + mn2 * md2m2
   90 Continue
  100 inns = inn
c
c...  initialize b-spline derivative values basef
c
      do 125 j=1,nd
           k1 = ind(j)
           mj = mnd(j)
           kmj = k + mj
           mj1 = mj + 1
           ivn = kmj * (kmj-1) / 2
           do 120 l=1,kmj
             basef(k1) = vn(ivn+l)
             do 110 jj=1,mj
                basef(k1+jj) = 0.d0
  110        Continue
          k1 = k1 + mj1
  120    Continue
  125  Continue
c
c...  for each derivative nr do loop 310
c
      do 310 nr=1,md
           nr1 = nr + 1
           mdr = md - nr
           k1 = ind(nd) + nr
           kmdr = k + mdr
           ivn = kmdr * (kmdr-1) /2
           if (mdr .eq. 0)                          go to 150
c
c...       first, determine nr(th) derivative of b-splines
c...       corresponding to the highest order solution component
c...       (i.e. of order mncomp=md).
c
           do 140 j=1,mdr
             jr = j + nr
             jin = jr + nr1 * kmd
             jink = jin + k
             do 130 l=j,jr
               jin1 = jin - kmd1
               jink1 = jink - kmd1
               alphd(jin) = dn1 * (alphdo(jin) - alphdo(jin1))
               alphd(jink) = dn3 * (alphdo(jink) - alphdo(jink1))
               in = k1 + (l-1) * md1
               basef(in) = basef(in) + alphd(jin) * vn(ivn+j)
               in = in + k * md1
               basef(in) = basef(in) + alphd(jink) * vn(ivn+j+k)
               jin = jin - kmd
               jink = jink - kmd
  130        Continue
  140      continue
  150      mdr1 = mdr + 1
           if ( mdr1 .gt. k)                        go to 180
           do 170 j = mdr1,k
             jr = j + nr
             jin = jr + nr1 * kmd
             do 160 l = j,jr
               jin1 = jin - kmd1
               alphd(jin) = dn2 * (alphdo(jin) - alphdo(jin1))
               in = k1 + (l-1) * md1
               basef(in) = basef(in) + alphd(jin) * vn(ivn+j)
               jin = jin - kmd
  160        Continue
  170      continue
  180      continue
           if (nd .eq. 1)                           go to 230
           inn = inns
c
c...       now determine nr(th) derivative basef for b-splines
c...       corresponding to all other solution components (nn)
c
           do 220 nn=1,ndm1
             nj = nd - nn
             mj = mnd(nj)
             inn = inn - (mj+2) * md2m2
             if (nr .gt. mj)                        go to 230
             kmjr = k + mj - nr
             k1 = ind(nj)+ nr
             ivn = kmjr * (kmjr-1) / 2
             mj1 = mj + 1
             jr1 = kmjr - md + 1
             jr1 = min0 (jr1, md-1)
c
c...         compute portion of b-spline derivative values (basef)
c...         using divided differences previously calculated for the
c...         highest order solution component in alphd.
c
             do 195 j=1,jr1
               jr = j + nr
               jin = jr + nr1 * kmd + md - mj
               do 190 l=j,jr
                 in = k1 + (l-1) * mj1
                 basef(in) = basef(in) + alphd(jin) * vn(ivn+j)
                 jin = jin - kmd
  190          continue
  195        Continue
             do 205 j=md,kmjr
               jr = j + nr
               jin = jr + nr1 * kmd
               do 200 l=j,jr
                 in = k1 + (l-1) * mj1
                 basef(in) = basef(in) + alphd(jin) * vn(ivn+j)
                 jin = jin - kmd
  200          continue
  205        Continue
c
c...         finish computing b-spline derivative values using the
c...         new nr(th) divided differences alphn
c
             jr2 = md2m2 - kmjr
             if (jr2 .le. 0)                        go to 220
             do 215 jj=1,jr2
               j = jj + jr1
               jr = j + nr
               jin = jr + nr1 * md2m2 + inn
               do 210 l=j,jr
                 jin1 = jin - md2m1
                 alphn(jin) = dn2 * (alphno(jin) - alphno(jin1))
                 in = k1 + (l-1) * mj1
                 basef(in) = basef(in) + alphn(jin) * vn(ivn+j)
                 jin = jin - md2m2
  210        continue
  215       Continue
  220      continue
  230      continue
c
c...       save nr(th) divided difference values, alphd and alphn,
c...       to be used to determine the next higher order divided
c...       differences, by storing them in alphdo and alphno
c
           if (nr .eq. md)                          go to 300
           nr2 = nr + 2
           inj = nr
           do 245 l=2,nr2
             inj = inj + kmd
             do 240 j=1,kmdr
                alphdo(j+inj) = alphd(j+inj)
  240        Continue
  245      Continue
           if (nd .eq. 1)                           go to 300
           if (nrest .le. 0)                        go to 300
           inn = 0
           do 290 nn = 1,ndm1
             mn = mnd(nn)
             if (mn .le. nr)                        go to 280
             kmnr = k + mn - nr
             jr1 = min0 (kmnr-md+1, md-1)
             inj = nr + inn
             inl = nr + md - mn
             do 255 l=2,nr2
               inj = inj + md2m2
               inl = inl + kmd
               do 250 j=1,jr1
                 alphno(inj+j) = alphd(inl+j)
  250          Continue
  255        continue
             mup = min0 (kmnr, md2m2)
             inj = nr + inn
             inl = nr
             do 265 l=2,nr2
               inj = inj + md2m2
               inl = inl + kmd
               do 260 j =md,mup
                 alphno(inj+j) = alphd(inl+j)
  260          Continue
  265        Continue
             jr2 = md2m2 - kmnr
             if (jr2 .le. 0)                        go to 280
             inj = nr + inn
             do 275 l=2,nr2
               inj = inj + md2m2
               do 270 jj = 1,jr2
                 jin = inj + jj + jr1
                 alphno(jin) = alphn(jin)
  270          Continue
  275        Continue
  280        inn = inn + (mn+2) * md2m2
  290      continue
  300      continue
  310 continue
c
c...  properly normalize basef values
c
      do 325 j=1,nd
           in = ind(j)
           icons = 1
           mj = mnd(j)
           kmj = k + mj
           mj1 = mj + 1
           do 322 nr = 1,mj
             icons = icons * (kmj-nr)
             in = in + 1
             do 320 l=1,kmj
               lbasef = in + (l-1) * mj1
               basef(lbasef) = basef(lbasef) * dfloat(icons)
  320 continue
  322 Continue
  325 Continue
c
c...  copy basef values corresponding to equal order solution components
c
      if (neq .eq. 0)                               return
      jd = 1
      do 360 j=1,neq
           in1 = ineq(j)
  330      if (in1 .lt. ind(jd+1))                  go to 340
           jd = jd + 1
           go to 330
  340      mj = mnd(jd)
           ntot = (k+mj)*(1+mj)
           in2 = ind(jd)
           do 350 l=1,ntot
             basef(in1-1+l) = basef(in2-1+l)
  350      Continue
  360 continue
      return
      end
      subroutine sysappdif (aldif, alpha, xi, n, k, ncomp, m, mstar)
c
c***********************************************************************
c
c   purpose
c           compute a divided difference table based upon the vector
c           of solution components
c
c   variables
c     alpha  - vector of solution coefficients (for all components)
c              corresponding to the mesh xi(1),...,xi(n+1)
c     aldif  - the divided difference array based upon alpha, viz,
c              aldif(i,r,j) = (r-1)st divided difference of alpha
c                             corresponding to u (x), for
c                                               j
c                             i=r,...,k+n+mj; r=1,...,mj; j=1,...,ncomp
c
c***********************************************************************
c
      implicit real(kind=8) (a-h,o-z)
      dimension aldif(*), alpha(*), xi(*), m(*)
      kd = k * ncomp
      incomp = 0
      k3 = 0
      k4 = 0
c
c...  construct the difference table for each component.
c
      do 130 j=1,ncomp
           mj = m(j)
           kmj = k - mj
           mjm1 = mj - 1
           kmr = k + mj
           nalphj = n * k + mj
           inn = incomp
           k1 = mstar
           k2 = kd
           k5 = inn + 1
c
c...       copy alpha into the first rows (nr=0) of aldif
c
           do 10 l=1,mj
             aldif(k5) = alpha(k3+l)
             k5 = k5 + 1
   10      Continue
           do 50 i = 1, n
             if (kmj .eq. 0)                        go to 30
             do 20 l = 1, kmj
               aldif(k5) = alpha(k1+k4+l)
               k5 = k5 + 1
   20        Continue
   30        do 40 l = 1, mj
               aldif(k5) = alpha(k2+k3+l)
               k5 = k5 + 1
   40        Continue
             k1 = k1 + kd
             k2 = k2 + kd
   50      continue
c
c...       for each derivative nr compute divided differences
c
           if (mjm1 .eq. 0)                         go to 120
           do 110 nr = 1, mjm1
             inn1 = inn + nalphj
             kmr = kmr - 1
             mjr = mj - nr
             kmjr = k - mjr
             xip1 = xi(1)
             dnk2 = dfloat(kmr) / (xi(2) - xip1)
c
c...         for xi(1),xi(2), the divided difference is a special case
c
             do 60 l=1,nr
               aldif(inn1+l) = 0.d0
   60        Continue
             do 70 l = nr, mjm1
               l1 = l + 1
             aldif(inn1+l1) = (aldif(inn+l1) - aldif(inn+l)) *
     1       dnk2
   70        Continue
             ibeg1 = mj
             ibeg2 = k + nr
c
c...         now the divided difference calculations for xi(i),xi(i+1),
c...                                                     i=1,...,n
c
             do 100 i = 1, n
               xii = xip1
               xip1 = xi(i+1)
               dnk1 = dfloat(kmr) / (xip1 - xii)
               if (i .lt. n) dnk2 = dfloat(kmr) / (xi(i+2) - xii)
               if (i .eq. n) dnk2 = dnk1
c
c...           the actual calculations involve two loops
c
               do 80 l = 1, kmjr
                 l1 = ibeg1 + l
               aldif(inn1+l1) = (aldif(inn+l1) - aldif(inn+l1-1)) *
     1         dnk1
   80           Continue     
               do 90 l = 1, mjr
                 l1 = ibeg2 + l
              aldif(inn1+l1) = (aldif(inn+l1) - aldif(inn+l1-1)) *
     1         dnk2
   90     Continue
               ibeg1 = ibeg1 + k
               ibeg2 = ibeg2 + k
  100        continue
             inn = inn1
  110      continue
  120      continue
           k3 = k3 + mj
           k4 = k4 + kmj
           incomp = incomp + nalphj * mj
  130 continue
      return
      end
      subroutine syshorder (i, uhigh, xiold, aldif)
c
c***********************************************************************
c
c   purpose
c           determine highest order (piecewise constant) derivatives
c           of the current collocation solution
c
c   variables
c     aldif  - divided differences of the solution coefficients alpha
c     uhigh  - the array of highest order (piecewise constant)
c              derivatives of the sysapproximate solution on
c              (xiold(i),xiold(i+1)), viz,
c                          (k+mj-1)
c              uhigh(j) = u   (x)    on (xiold(i),xiold(i+1))
c                          j
c
c***********************************************************************
c
      implicit real(kind=8) (a-h,o-z)
      common /appr/ n,nold,nmax,nalpha,mshflg,mshnum,mshlmt,mshalt
      common /order/ k,ncomp,mstar,kd,kdm,mnsum,m(20)
      dimension uhigh(*) , ar(20), arm1(20)
      dimension aldif(*), xiold(*)
c
      dn2 = 1.d0 / (xiold(i+1) - xiold(i))
      incomp = 0
c
c...  loop through the ncomp solution components
c
      do 50 j = 1,ncomp
           mj = m(j)
           nalphj = k * nold + mj
           kpmj = k + mj
           kmr = k + 1
           mjm1 = mj - 1
           incomp = incomp + mjm1 * nalphj
           left = i * k + mj - kmr
c
c...       further divided differences of the appropriate aldif
c...       (viz. of the (mj-1)st divided differences of the alpha) are
c...       calculated to obtain the (k+mj-1)st divided difference
c
           do 10 l=1,kmr
             leftpl = left + l
             arm1(l+mj-1) = aldif(incomp+leftpl)
   10      Continue
           incomp = incomp + nalphj
           kpmj1 = kpmj - 1
           do 40 nr = mj,kpmj1
             kmr = kmr - 1
             dnk2 = dn2 * dfloat(kmr)
             do 20 l = 1,kmr
               ar(l+nr) = dnk2 * (arm1(l+nr) - arm1(l+nr-1))
   20        Continue
             do 30 l=nr,kpmj1
               arm1(l+1) = ar(l+1)
   30        Continue
   40      continue
           uhigh(j) = ar(kpmj)
   50 continue
      return
      end
c-----------------------------------------------------------------------
c          for convenience of the user we list here the package
c          solveblok  of de boor - weiss (4), used in  colsys.
c-----------------------------------------------------------------------
c
      subroutine sysfcblok (bloks, integs, nbloks, ipivot, scrtch,iflag)
c
c******************************************************************
c
c    calls subroutines  sysfactrb  and  sysshiftb .
c
c     sysfcblok  supervises the plu factorization with pivoting of
c     scaled rows of the almost block diagonal matrix stored in the
c     arrays  bloks  and  integs .
c
c     sysfactrb = subprogram which carries out steps 1,...,last of gauss
c            elimination (with pivoting) for an individual block.
c     sysshiftb = subprogram which shifts the remaining rows to the top of
c            the next block
c
c     parameters
c      bloks   an array that initially contains the almost block diagona
c            matrix  a  to be factored, and on return contains the com-
c            puted factorization of  a .
c      integs  an integer array describing the block structure of  a .
c      nbloks  the number of blocks in  a .
c      ipivot  an integer array of dimension   sum (integs(1,n) ; n=1,
c            ...,nbloks) which, on return, contains the pivoting stra-
c            tegy used.
c      scrtch  work area required, of length  max (integs(1,n) ; n=1,
c            ...,nbloks).
c      iflag   output parameter;
c            = 0  in case matrix was found to be singular.
c            otherwise,
c            = (-1)**(number of row interchanges during factorization)
c
c***********************************************************************
c
      integer integs(3,nbloks),ipivot(*),iflag, i,index,indexb,indexn,
     1        last,ncol,nrow
      double precision bloks(*),scrtch(*)
      iflag = 1
      indexb = 1
      indexn = 1
      i = 1
c
c...  loop over the blocks.  i  is loop index
c
   10      index = indexn
           nrow = integs(1,i)
           ncol = integs(2,i)
           last = integs(3,i)
c
c...       carry out elimination on the i-th block until next block
c...       enters, i.e., for columns 1,...,last  of i-th block.
c
           call sysfactrb ( bloks(index), ipivot(indexb), scrtch, nrow,
     1          ncol, last, iflag)
c
c...       check for having reached a singular block or the last block
c
           if (iflag .eq. 0 .or. i .eq. nbloks)     return
           i = i+1
           indexn = nrow*ncol + index
c
c...       put the rest of the i-th block onto the next block
c
           call sysshiftb ( bloks(index), ipivot(indexb), nrow, ncol,
     1          last, bloks(indexn), integs(1,i), integs(2,i) )
           indexb = indexb + nrow
      go to 10
      end
      subroutine sysfactrb ( w, ipivot, d, nrow, ncol, last, iflag )
c
c***********************************************************************
c
c     adapted from p.132 of  element.numer.analysis  by conte-de boor
c
c     constructs a partial plu factorization, corresponding to steps
c      1,..., last   in gauss elimination, for the matrix  w  of
c      order ( nrow ,  ncol ), using pivoting of scaled rows.
c
c     parameters
c       w       contains the (nrow,ncol) matrix to be partially factored
c               on input, and the partial factorization on output.
c       ipivot  an integer array of length nrow containing a record of
c               the pivoting strategy used; row ipivot(i) is used
c               during the i-th elimination step, i=1,...,last.
c       d       a work array of length nrow used to store row sizes
c               temporarily.
c       nrow    number of rows of w.
c       ncol    number of columns of w.
c       last    number of elimination steps to be carried out.
c       iflag   on output, equals iflag on input times (-1)**(number of
c               row interchanges during the factorization process), in
c               case no zero pivot was encountered.
c               otherwise, iflag = 0 on output.
c
c***********************************************************************
c
      integer ipivot(nrow),ncol,last,iflag, i,ipivi,ipivk,j,k,kp1
      double precision w(nrow,ncol),d(nrow), awikdi,colmax,ratio,rowmax
      double precision dabs,dmax1
c
c...  initialize ipivot, d
c
      do 20 i=1,nrow
           ipivot(i) = i
           rowmax = 0.d0
           do 10 j=1,ncol
             rowmax = dmax1(rowmax, dabs(w(i,j)))
   10      Continue
           if (rowmax .eq. 0.d0)                    go to 90
           d(i) = rowmax
   20   Continue
c
c...  gauss elimination with pivoting of scaled rows, loop over
c...  k=1,.,last
c
      k = 1
c
c...  as pivot row for k-th step, pick among the rows not yet used,
c...  i.e., from rows ipivot(k),...,ipivot(nrow), the one whose k-th
c...  entry (compared to the row size) is largest. then, if this row
c...  does not turn out to be row ipivot(k), redefine ipivot(k) ap-
c...  propriately and record this interchange by changing the sign
c...  of  iflag .
c
   30      ipivk = ipivot(k)
           if (k .eq. nrow)                         go to 80
           j = k
           kp1 = k+1
           colmax = dabs(w(ipivk,k))/d(ipivk)
c
c...       find the (relatively) largest pivot
c
           do 40 i=kp1,nrow
             ipivi = ipivot(i)
             awikdi = dabs(w(ipivi,k))/d(ipivi)
             if (awikdi .le. colmax)                go to 40
             colmax = awikdi
             j = i
   40      continue
           if (j .eq. k)                            go to 50
           ipivk = ipivot(j)
           ipivot(j) = ipivot(k)
           ipivot(k) = ipivk
           iflag = -iflag
   50      continue
c
c...       if pivot element is too small in absolute value, declare
c...       matrix to be noninvertible and quit.
c
           if (dabs(w(ipivk,k))+d(ipivk) .le. d(ipivk))
     1                                              go to 90
c
c...       otherwise, subtract the appropriate multiple of the pivot
c...       row from remaining rows, i.e., the rows ipivot(k+1),...,
c...       ipivot(nrow), to make k-th entry zero. save the multiplier
c...       in its place.
c
           do 65 i=kp1,nrow
             ipivi = ipivot(i)
             w(ipivi,k) = w(ipivi,k)/w(ipivk,k)
             ratio = -w(ipivi,k)
             do 60 j=kp1,ncol
               w(ipivi,j) = ratio*w(ipivk,j) + w(ipivi,j)
   60        Continue
   65   Continue
           k = kp1
c
c...       check for having reached the next block.
c
           if (k .le. last)                         go to 30
      return
c
c...  if  last  .eq. nrow , check now that pivot element in last row
c...  is nonzero.
c
   80 if(dabs(w(ipivk,nrow))+d(ipivk) .gt. d(ipivk)) return
c
c...  singularity flag set
c
   90 iflag = 0
      return
      end
      subroutine sysshiftb ( ai, ipivot, nrowi, ncoli, last,
     1           ai1, nrowi1, ncoli1 )
c
c**********************************************************************
c
c     shifts the rows in current block, ai, not used as pivot rows, if
c     any, i.e., rows ipivot(last+1),...,ipivot(nrowi), onto the first
c     mmax = nrow-last rows of the next block, ai1, with column last+j
c     of ai  going to column j , j=1,...,jmax=ncoli-last. the remaining
c     columns of these rows of ai1 are zeroed out.
c
c                                picture
c
c          original situation after         results in a new block i+1
c          last = 2 columns have been       created and ready to be
c          done in sysfactrb (assuming no      factored by next sysfactrb
c          interchanges of rows)            call.
c                      1
c                 x  x 1x  x  x           x  x  x  x  x
c                      1
c                 0  x 1x  x  x           0  x  x  x  x
c     block i          1                       ---------------
c     nrowi = 4   0  0 1x  x  x           0  0 1x  x  x  0  01
c     ncoli = 5        1                       1             1
c     last = 2    0  0 1x  x  x           0  0 1x  x  x  0  01
c     -------------------------------          1             1   new
c                      1x  x  x  x  x          1x  x  x  x  x1  block
c                      1                       1             1   i+1
c     block i+1        1x  x  x  x  x          1x  x  x  x  x1
c     nrowi1= 5        1                       1             1
c     ncoli1= 5        1x  x  x  x  x          1x  x  x  x  x1
c     -------------------------------          1-------------1
c                      1
c
c***********************************************************************
c
      integer ipivot(nrowi),last, ip,j,jmax,jmaxp1,m,mmax
      double precision ai(nrowi,ncoli),ai1(nrowi1,ncoli1)
      mmax = nrowi - last
      jmax = ncoli - last
      if (mmax .lt. 1 .or. jmax .lt. 1)             return
c
c...  put the remainder of block i into ai1
c
      do 15 m=1,mmax
           ip = ipivot(last+m)
           do 10 j=1,jmax
             ai1(m,j) = ai(ip,last+j)
   10      Continue
   15 Continue 
      if (jmax .eq. ncoli1)                         return
c
c...  zero out the upper right corner of ai1
c
      jmaxp1 = jmax + 1
      do 30 j=jmaxp1,ncoli1
           do 20 m=1,mmax
             ai1(m,j) = 0.d0
   20      Continue
   30 Continue
      return
      end
      subroutine syssbblok ( bloks, integs, nbloks, ipivot, b, x )
c
c**********************************************************************
c
c     calls subroutines  syssubfor  and  syssubbak .
c
c     supervises the solution (by forward and backward substitution) of
c     the linear system  a*x = b  for x, with the plu factorization of
c     a  already generated in  sysfcblok .  individual blocks of
c     equations are solved via  syssubfor  and  syssubbak .
c
c    parameters
c       bloks, integs, nbloks, ipivot    are as on return from sysfcblok.
c       b       the right side, stored corresponding to the storage of
c               the equations. see comments in s l v b l k for details.
c       x       solution vector
c
c***********************************************************************
c
      integer integs(3,nbloks),ipivot(*), i,index,indexb,indexx,j,last,
     1        nbp1,ncol,nrow
      double precision bloks(*),b(*),x(*)
c
c...  forward substitution pass
c
      index = 1
      indexb = 1
      indexx = 1
      do 10 i=1,nbloks
           nrow = integs(1,i)
           last = integs(3,i)
           call syssubfor ( bloks(index), ipivot(indexb), nrow, last,
     1          b(indexb), x(indexx) )
           index = nrow*integs(2,i) + index
           indexb = indexb + nrow
           indexx = indexx + last
   10  Continue
c
c...  back substitution pass
c
      nbp1 = nbloks + 1
      do 20 j=1,nbloks
           i = nbp1 - j
           nrow = integs(1,i)
           ncol = integs(2,i)
           last = integs(3,i)
           index = index - nrow*ncol
           indexb = indexb - nrow
           indexx = indexx - last
       call syssubbak ( bloks(index), ipivot(indexb), nrow, ncol,
     1     last, x(indexx) )
   20 Continue
      return
      end
      subroutine syssubfor ( w, ipivot, nrow, last, b, x )
c
c***********************************************************************
c
c     carries out the forward pass of substitution for the current
c     block, i.e., the action on the right side corresponding to the
c     elimination carried out in  sysfactrb  for this block.
c        at the end, x(j) contains the right side of the transformed
c     ipivot(j)-th equation in this block, j=1,...,nrow. then, since
c     for i=1,...,nrow-last, b(nrow+i) is going to be used as the right
c     side of equation  i  in the next block (shifted over there from
c     this block during factorization), it is set equal to x(last+i)
c     here.
c
c    parameters
c       w, ipivot, nrow, last  are as on return from sysfactrb.
c       b(j)   is expected to contain, on input, the right side of j-th
c              equation for this block, j=1,...,nrow.
c       b(nrow+j) contains, on output, the appropriately modified right
c              side for equation j in next block, j=1,...,nrow-last.
c       x(j)   contains, on output, the appropriately modified right
c              side of equation ipivot(j) in this block, j=1,...,last
c              (and even for j=last+1,...,nrow).
c
c***********************************************************************
c
      integer ipivot(nrow), ip,jmax,k
      double precision w(nrow,last),b(*),x(nrow),sum
      ip = ipivot(1)
      x(1) = b(ip)
      if (nrow .eq. 1)                              go to 40
      do 20 k=2,nrow
           ip = ipivot(k)
           jmax = min0(k-1,last)
           sum = 0.d0
           do 10 j=1,jmax
             sum = w(ip,j)*x(j) + sum
   10      Continue
         x(k) = b(ip) - sum
   20 Continue
c
c...  transfer modified right sides of equations ipivot(last+1),...,
c...  ipivot(nrow) to next block.
c
      nrowml = nrow - last
      if (nrowml .eq. 0)                            go to 40
      lastp1 = last+1
      do 30 k=lastp1,nrow
        b(nrowml+k) = x(k)
   30 Continue
   40 return
      end
      subroutine syssubbak ( w, ipivot, nrow, ncol, last, x )
c
c***********************************************************************
c
c     carries out backsubstitution for current block.
c
c    parameters
c       w, ipivot, nrow, ncol, last  are as on return from sysfactrb.
c       x(1),...,x(ncol)  contains, on input, the right side for the
c               equations in this block after backsubstitution has been
c               carried up to but not including equation ipivot(last).
c               means that x(j) contains the right side of equation ipi-
c               vot(j) as modified during elimination, j=1,...,last,
c               while for j .gt. last, x(j) is already a component of
c               the solution vector.
c       x(1),...,x(ncol) contains, on output, the components of the
c               solution corresponding to the present block.
c
c**********************************************************************
c
      integer ipivot(nrow),last,  ip,j,k,kp1
      double precision w(nrow,ncol),x(ncol), sum
      k = last
      ip = ipivot(k)
      sum = 0.d0
      if (k .eq. ncol)                              go to 30
      kp1 = k+1
   10 do 20 j=kp1,ncol
        sum = w(ip,j)*x(j) + sum
   20 continue
   30 x(k) = (x(k) - sum)/w(ip,k)
      if (k .eq. 1)                                 return
      kp1 = k
      k = k-1
      ip = ipivot(k)
      sum = 0.d0
      go to 10
      end
