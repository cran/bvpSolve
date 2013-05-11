
c ===================================================================================
c main driver for twpbvplc, written by Jeff Cash and Francesca Mazzia
c ===================================================================================

       subroutine twpbvplc(ncomp, nlbc, aleft, aright,
     *       nfxpnt, fixpnt, ntol, ltol, tol,
     *       linear, givmsh, giveu, nmsh,
     *       nxxdim, xx, nudim, u, nmax,
     *       lwrkfl, wrk, lwrkin, iwrk, precis,
     *       fsub, dfsub, gsub, dgsub,
     *       ckappa1,gamma1,sigma,ckappa,
     *       ckappa2,rpar,ipar,iflbvp,liseries,iseries,indnms,
     *       full, useC,  nmguess, xguess, nygdim,yguess, iset)

*     OUTPUT
*
*     IFLBVP = 0   SUCCESFULL TERMINATION
*
*     IFLBVP = -1 (SUCCESSFULL TERMINATION BUT CONDITIONING PARAMETER
*                    NOT STABILIZED (ONLY IF USE_C = .true.)
*
*     IFLBVP = 1   terminated too many mesh points
*     IFLBVP = 2   terminated too many meshes (maximum number of
*                       meshes is liseries )
*     IFLBVP = 3   terminated ill conditioned problem
*                            (ONLY IF USE_C = .true.)
*     IFLBVP = 4   terminated invalid input
*

* WORKSPACE IN INPUT
*
* WRK workspace WRK(LWRKFL)
*   LWRKFL >= 24*NCOMP+16*NCOMP*NCOMP+14*NCOMP*NXXDIM
*                 +5*NCOMP*NCOMP*NXXDIM+4*NXXDIM)
* IWRK workspace IWRK(LWRKIN)
*    LWRKIN >= 2*NCOMP*NXXDIM+2*NXXDIM+3*NCOMP)
*
*  The subroutine twpbvpl is intended to solve two-point boundary
*  value problems.
* References:
*
*  CASH, J. R. AND MAZZIA, F. 2006. Hybrid mesh selection algorithms
*  based on conditioning for two-point boundary value problems.
*  JNAIAM J. Numer. Anal. Ind. Appl. Math. 1, 1, 81-90.
*
*  CASH, J. R. AND MAZZIA, F. 2009. Conditioning and Hybrid Mesh
*   Selection Algorithms For Two Point Boundary Value Problems.
*   Scalable Computing: Practice and Experience 10, 4, 347-361.
*
*   Revision History
*
* revision  July 10,  2006
*   added rpar and ipar in the functions
*   DFSUB, DGSUB, FSUB, GSUB
*   changed the name of the variable double in ddouble
*
*
* revision 31 August 2004
* This is a modified version of twpbvp that uses the conditioning
* in the mesh selection and lobatto methods.
*
*     New subroutines not included in the old version:
*           condestim
*           moncond
*           selcondmsh
*           selconderrmsh
*           smpselcondmsh
*
*     Updates subroutines:
*           bvpsol
*           fail4
*           fail6
*           newteq_l
*           mshref
*
*     Auxiliary function not used by the old version
*           donest
*           abdnrm
*
*     The common block algprs contains two more variable
*     logical  use_c, comp_c
*     common/algprs/ nminit,  iprint, idum,  use_c, comp_c
*
*     Revision February 2011 (Francesca Mazzia)
*          changed the input parameters in order to have the same
*          input of twpbvpc that works with R
*          All the subroutine in common with twpbvpc are in a companion
*          file called  twpbvpa.f
*


      implicit double precision (a-h,o-z)
      integer nudim, nygdim, lwrkfl, lwrkin
      dimension rpar(*),ipar(*), iset(*), precis(3)
      dimension fixpnt(*), ltol(*), tol(*)
      dimension xx(*), u(nudim,*), xguess(*), yguess(nygdim,*)
      dimension wrk(lwrkfl), iwrk(lwrkin)
      dimension iseries(*)
      logical linear, givmsh, giveu, full, useC
      external fsub
      external dfsub
      external  gsub
      external  dgsub

      integer nminit, iprint, idum
      logical  use_c, comp_c, giv_u
      common/algprs/ nminit, iprint, idum, use_c, comp_c
C Karline: added
      integer ureset
      common/gu/ giv_u, ureset
      integer nfunc, njac, nstep, nbound, njacbound
      common/diagnost/nfunc, njac, nstep, nbound, njacbound
      intrinsic abs
      intrinsic  min
      parameter ( zero = 0.0d+0 )


      use_c = useC

c Karline: set this = 1 => always computed...
      COMP_C =  .TRUE.
      giv_u = giveu
cF     the block data is now in bvpsol
cF      block data
c
*  This block data routine initializes nminit (the initial number
*  of mesh points), uval0 (the initial value for the trial
*  solution) to their default values, use_c (if the conditioning is used
*  in the mesh selection), and comp_c (if the conditining parameters
*  are computed)

      nminit = 11
      iprint = -1

c initialise counters
      ureset = 0
      nfunc = 0
      njac = 0
      nstep = 0
      nbound = 0
      njacbound = 0

      if (full) then
        iprint = 1
      else
        iprint = -1
      endif

*  Check for invalid input parameters.  If any parameters are
*  invalid, exit with the flag iflbvp set to 4.


      iflbvp = 4
      if (ncomp .le. 0)  return
      if (nlbc .lt. 0 .or. nlbc .gt. ncomp) return
      if (aleft .ge. aright) return

      if (nfxpnt .lt. 0)  return
      if (givmsh .and. nmsh .lt. nfxpnt+2) return
      if (givmsh .and. xx(1) .ne. aleft) return
C     SCMODIFIED add an extra condition to avoid accessing xx(0)
      if (nmsh .gt. 0) then
        if (givmsh .and. xx(nmsh) .ne. aright) return
      end if
      if (nfxpnt .gt. 0) then
         if (fixpnt(1) .le. aleft) return
         if (fixpnt(nfxpnt) .ge. aright) return
         do 50 i = 1, nfxpnt-1
            if (fixpnt(i+1) .le. fixpnt(i)) return
 50             continue
      endif

      if (ntol .lt. 1) return
      do 60 i = 1, ntol
         if (ltol(i) .lt. 0 .or. ltol(i) .gt. ncomp) return
         if (tol(i) .le. zero) return
 60       continue

      if (giveu .and. .not. givmsh) return

      if (use_c .and. .not. comp_c) return

      if (nudim .le. 0) return
      if (lwrkfl .le. 0 .or. lwrkin .le. 0) return


*  Calculate maximum number of mesh points possible with the
*  given floating-point and integer workspace.

      isp = lwrkfl  - 2*ntol - 24*ncomp - 14*ncomp*ncomp
      iden = 5*ncomp*ncomp + 14*ncomp + 4
      nmax1 = isp/iden

      isp = lwrkin - 3*ncomp
      nmax2 = isp/(2*ncomp+2)

      nmax = min(nmax1,nmax2)
* nmax from workspace
      nmax = min(nmax, nxxdim)
* nmax from size of u and xx

      if (iprint .ge. 0) then
       CALL Rprinti1('Nmax from workspace =',nmax)
      end if

      if (nmax .le. 1) return


*  Partition floating point workspace.

      irhs = 1
      lrhs = ncomp*nmax
* 1 ncomp*nmax
      itpblk = irhs + lrhs
      ltpblk = ncomp*nlbc
* 1 ncomp*ncomp
      ibtblk = itpblk + ltpblk
      lbtblk = ncomp*(ncomp - Nlbc)
* 2 ncomp*ncomp
      iajac = ibtblk + lbtblk
      lajac = 2*ncomp*ncomp*nmax
* 2 ncomp*ncomp*nmax
      ibhold = iajac + lajac
      lbhold = ncomp*ncomp*nmax
* 3 ncomp*ncomp*nmax
      ichold = ibhold + lbhold
      lchold = ncomp*ncomp*nmax
* 4 ncomp*ncomp*nmax
      ifval = ichold + lchold
      lfval = ncomp*nmax
* 2 ncomp*nmax
      idef = ifval + lfval
      ldef = ncomp*(nmax-1)
* 3 ncomp*nmax
      idefex = idef + ldef
*      ldefex = ncomp*(nmax-1)

*  def6 uses the same space as defexp

      idef6 = idefex
      ldef6 = ncomp*(nmax-1)
* 4 ncomp*nmax
      idefim = idef6 + ldef6
*      ldefim = ncomp*(nmax-1)

*  def8 uses the same space as defimp

      idef8 = idefim
      ldef8 = ncomp*(nmax-1)
* 5 ncomp*nmax
      iusve = idef8 + ldef8
      lusve = ncomp*nmax
* 6 ncomp*nmax
      iuold = iusve + lusve
      luold = ncomp*nmax
* 7 ncomp*nmax
      itmrhs = iuold + luold
      ltmrhs = ncomp*nmax
* 8 ncomp*nmax
      irhtri = itmrhs + ltmrhs
      lrhtri = ncomp*nmax
* 9 ncomp*nmax
      idelu = irhtri + lrhtri
      ldelu = ncomp*nmax
* 10 ncomp*nmax
      ixmer = idelu + ldelu
*      lxmer = ncomp*nmax

*  rerr occupies the same space as xmerit
      irerr = ixmer
      lrerr = ncomp*nmax
* 11 ncomp*nmax
      iutri = irerr + lrerr
      lutri = ncomp*nmax
* 12 ncomp*nmax
      iermx = iutri + lutri
      lermx = nmax
* 1  nmax
      irtdc = iermx + lermx
      lrtdc = nmax
* 2  nmax
      ixxold = irtdc + lrtdc
      lxxold = nmax
* 3  nmax
      iuint = ixxold + lxxold
      luint = ncomp
* 1 ncomp
      iftmp = iuint + luint
      lftmp = ncomp
* 2 ncomp
      idgtm = iftmp + lftmp
      ldgtm = ncomp
* 3 ncomp
      idftm1 = idgtm + ldgtm
      ldftm1 = ncomp*ncomp
* 3 ncomp*ncomp
      idftm2 = idftm1 + ldftm1
      ldftm2 = ncomp*ncomp
* 4 ncomp*ncomp
      itmp = idftm2 + ldftm2
      ltmp = ncomp*20
* 23 ncomp
      idsq = itmp + ltmp
      ldsq = ncomp*ncomp
* 5 ncomp*ncomp
      idexr = idsq + ldsq
      ldexr = ncomp
* 24 ncomp
      ietst6 = idexr + ldexr
      letst6 = ntol
* 1 ntol
      ietst8 = ietst6 + letst6
      letst8 = ntol
* 2 ntol
      iamg = ietst8 + letst8
      lamg = ncomp*nmax
* 13 ncomp*nmax
      ic1 = iamg + lamg
      lc1 = ncomp*ncomp*nmax
* 5 ncomp*ncomp*nmax

      iwrkrhs = ic1+lc1
      lwrkrhs = ncomp*nmax
* 14 ncomp*nmax
      ir4 = iwrkrhs + lwrkrhs
      lr4 = nmax
* 4  nmax
      idhold = ir4 + lr4
      ldhold = 9*ncomp*ncomp
* 14 ncomp*ncomp
      ilast = idhold +  ldhold


*  Partition integer workspace.

      iiref = 1
      liref = nmax

      iihcom = iiref + liref
      lihcom = nmax

      iipvbk = iihcom + lihcom
      lipvbk = ncomp*nmax

      iipvlu = iipvbk + lipvbk
      lipvlu = 3*ncomp

      iisign = iipvlu + lipvlu
      lisign = ncomp*nmax

      if (iprint .eq. 1) then
        CALL Rprinti1('Integer workspace', iisign+lisign)
      end if

      call bvpsol_l(ncomp, nmsh, nlbc, aleft, aright,
     *   nfxpnt, fixpnt, ntol, ltol, tol, nmax, linear,
     *   giveu, givmsh, xx, nudim, u,
     *   wrk(idefex), wrk(idefim), wrk(idef), wrk(idelu),
     *   wrk(irhs), wrk(ifval),
     *   wrk(itpblk), wrk(ibtblk), wrk(iajac), wrk(ibhold),
     *   wrk(ichold), wrk(idhold), iwrk(iipvbk), iwrk(iipvlu),
     *   iwrk(iisign), wrk(iuint), wrk(iftmp), wrk(itmrhs),
     *   wrk(idftm1), wrk(idftm2), wrk(idgtm),
     *   wrk(iutri), wrk(irhtri), wrk(ixmer),
     *   wrk(ixxold), wrk(iuold), wrk(iusve),
     *   wrk(itmp), wrk(idsq), wrk(idexr), wrk(irtdc),
     *   wrk(irerr), wrk(ietst6), wrk(ietst8), wrk(iermx),
     *   iwrk(iihcom), iwrk(iiref), wrk(idef6), wrk(idef8),
     *   fsub,dfsub,gsub,dgsub,iflbvp,
     *   wrk(iamg),wrk(ic1),wrk(iwrkrhs),
     *   ckappa1,gamma1,sigma,ckappa,ckappa2,wrk(ir4),rpar,ipar,
     *    liseries,iseries,indnms,precis,
     *   nmguess,xguess, nygdim, yguess )


      iset(1) = nfunc
      iset(2) = njac
      iset(3) = nbound
      iset(4) = njacbound
      iset(5) = nstep
      iset(6) = ureset
C karline moved return statement to here...
      return

      end

c ===================================================================================

      subroutine bvpsol_l(ncomp, nmsh, nlbc, aleft, aright,
     *   nfxpnt, fixpnt,
     *   ntol, ltol, tol, nmax, linear, giveu, givmsh,
     *   xx, nudim, u, defexp, defimp, def, delu, rhs, fval,
     *   topblk, botblk, ajac, bhold, chold, dhold,
     *   ipvblk, ipivlu,isign,
     *   uint, ftmp, tmprhs, dftmp1, dftmp2, dgtm,
     *   utrial, rhstri, xmerit, xxold, uold, usave,
     *   tmp, dsq, dexr, ratdc, rerr,
     *   etest6, etest8, ermx, ihcomp, irefin,
     *   def6, def8, fsub, dfsub, gsub, dgsub, iflbvp,
     *   amg, c1, wrkrhs,ckappa1,gamma1,sigma,ckappa,ckappa2,r4,
     *   rpar,ipar,liseries,iseries,indnms,precis,
     *   nmguess,xguess, nygdim, yguess )

      implicit double precision (a-h,o-z)

      dimension rpar(*), ipar(*), precis(3)
      dimension  fixpnt(*), ltol(ntol), tol(ntol)
      dimension  xx(*), u(nudim, *), xguess(*), yguess(nygdim,*)
      dimension  defexp(ncomp,*), defimp(ncomp,*), def(ncomp,*)
      dimension  delu(ncomp, *), rhs(*), fval(ncomp,*)
      dimension  topblk(nlbc,*), botblk(ncomp-nlbc,*)
      dimension  ajac(ncomp, 2*ncomp, *)
      dimension  bhold(ncomp, ncomp, *), chold(ncomp, ncomp, *)
      dimension  ipivlu(*), ipvblk(*), isign(*)
      dimension  uint(ncomp), ftmp(ncomp)
      dimension  dgtm(ncomp), tmprhs(*)
      dimension  dftmp1(ncomp, ncomp), dftmp2(ncomp, ncomp)
      dimension  utrial(ncomp,*), rhstri(*)
      dimension  xmerit(ncomp, *)
      dimension  xxold(*), uold(ncomp,*), usave(ncomp,*)
      dimension  tmp(ncomp,*)
      dimension  dsq(ncomp,ncomp), dexr(ncomp)
      dimension  ratdc(*), rerr(ncomp,*)
      dimension  etest6(*), etest8(*), ermx(*)
      dimension  ihcomp(*), irefin(*)
      dimension  def6(ncomp,*), def8(ncomp,*)
      dimension  amg(*), c1(ncomp,*), wrkrhs(*)
      Dimension  Dhold(3*ncomp,*),df(ncomp,ncomp)
      Dimension  r4(*), iseries(*)

      logical linear, giveu, givmsh, ddouble

      external fsub
      external dfsub
      external gsub
      external dgsub

      common/mchprs/flmin, flmax, epsmch

      logical  use_c, comp_c, Chstif
      integer nminit, iprint, idum
      common/algprs/ nminit, iprint, idum,  use_c, comp_c
      common/monpar/ sfatt_alpha, sfatt_r3, sfatt_r1r3
      common/newt/greps
      intrinsic max

      logical smooth, succes, strctr, trst6, reaft6
      logical onto6, onto8, ludone, rhsgiv, maxmsh
      logical first4, first8
      logical nodouble, forcedouble
      logical reposs

      logical mchset
      save mchset
      logical frscal
c      save frscal


      logical stab_kappa, stab_gamma, stab_cond, stiff_cond, ill_cond
      logical stab_kappa1, ill_cond_newt, stab_sigma, comparekappa
      logical errok

      parameter (zero = 0.0d+0, one = 1.0d+0)
      parameter (third = 0.33d+0, fourth = 0.25d+0)
      parameter (quan6 = 0.5d+0 )
      parameter (itcondmax = 5)
      parameter (power = 1.0d+0/6.0d+0)
*  blas: dload
*  double precision d1mach


      data mchset/.true./
      data fxfct/10.0d+0/
      data maxmsh/.false./
      data itcond/0/
      data frscal/.true./
c      frscal = .true.
      if (mchset) then
c Karline: use precis instead of d1mach

         flmin = precis(1)
         flmax = precis(2)
         epsmch = precis(3)
         mchset = .false.
      endif






*  The routine stcons calculates integration constants stored in
*  labeled common consts.

      Call Stcon1
      Call Stcon2

*  Set up arrays for the error tests.


      if (.not. linear) then
         call dload(ntol, one, etest6, 1)
      else
         do 10 i = 1, ntol
            etest6(i) = one/max(quan6, tol(i)**third)
 10             continue
      endif

      call dload(ntol, 10d0, etest6, 1)


      do i=1,ncomp
         do j=1,ncomp
            df(i,j)=0.0d+0
         enddo
      enddo
      nmold = 1
      smooth = .false.
      strctr = .false.
      trst6 = .true.
      reaft6 = .false.
      numbig = 0
      nummed = 0
      first4 = .true.
      first8 = .true.
      onto6  = .false.
      maxmsh = .false.
      Chstif = .true.
      ddouble = .false.


      if (comp_c) then
*     initialize parameter for the conditioning estimation
      gamma1old  = flmax
      gamma1     = flmax
      ckappa1old = flmax
      ckappa1    = flmax
      ckappa     = flmax
      ckappaold  = flmax
      ckappa2    = flmax
      sigma      = flmax
      sigmaold   = flmax
      stiff_cond = .false.
      stab_cond  = .false.
      ill_cond   = .false.
      endif
*     initialize parameter for the conditioning estimation
      if (linear) then
         sfatt_alpha = 1e-5
         sfatt_r3 = 1d-5
      else
         sfatt_alpha = 1e-5
         sfatt_r3 = 1d-5
      end if

      sfatt_r1r3 = 0.65d0
      greps=10d0


      tolmin = flmax
      do i=1,ntol
         tolmin = min(tol(i),tolmin)
      end do
*  If givmsh is .true., the initial number of mesh points must be
*  provided by the user in nmsh, and the mesh points must be
*  contained in the array xx (of dimension nmsh).
*  Otherwise, nmsh is set to its default value, and a
*  uniform initial mesh is created.



      if (.not. giveu .and. .not. givmsh) then
         nmsh = nminit
           if (nmsh .ge. nmax ) then
            maxmsh = .true.
            if (iprint .eq. 1) THEN
             CALL Rprint('Initial mesh greater than nmax')
            ENDIF

            goto 900
         end if
         if (nmsh .lt. nfxpnt+2) nmsh = nfxpnt + 2
         call unimsh(nmsh, aleft, aright, nfxpnt, fixpnt, xx)
      endif
        if (nmsh .ge. nmax ) then
            maxmsh = .true.

            if (iprint .eq. 1) THEN
             CALL Rprint('Initial mesh greater than nmax')
            ENDIF

            goto 900
         end if

       if (.not. giveu) call initu(ncomp, nmsh, xx, nudim, u,
     *    nygdim, nmguess,xguess, yguess)
      indnms = 0
      indnmsold = 0
      nfail4=0
***** top of logic for 4th order solution ****

 400   continue
       if (indnmsold.ne.nmsh) then
       indnms = indnms + 1
       iseries(indnms) = nmsh
       indnmsold = nmsh
       endif

       if (indnms .ge. liseries) then
          goto 1900
      end if
      If(Maxmsh) goto 900

      if (iprint .eq. 1) THEN
       CALL Rprinti1('Start 4th order, nmsh', nmsh)
      ENDIF

*  Set the def (deferred correction) array to zero.

      call mtload(ncomp, nmsh-1, zero, ncomp, def)
      iorder = 4

*  The routine fneval calls fsub at the mesh xx and the
*  solution u, and saves the values in the array fval.

      call fneval(ncomp, nmsh, xx, nudim, u, fval, fsub,rpar,ipar)


*  Try to compute a 4th order solution by solving a system of nonlinear
*  equations.

      if (linear) then
         ludone = .false.

          call lineq( ncomp, nmsh, nlbc,
     *    ludone, xx, nudim, u, def,
     *    delu, rhs, fval,
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,
     *    ajac, topblk, botblk, bhold, chold, ipvblk,
     *    fsub, dfsub, gsub, dgsub, iflnwt,rpar,ipar)

*  Call fneval to evaluate the fval array at the new solution u.
*  (Such a call is not necessary for the nonlinear case because
*  fval is called within newteq for the new u.)

         call fneval(ncomp, nmsh, xx, nudim, u, fval, fsub,rpar,ipar)

      else

         rhsgiv = .false.
         call newteq(ncomp, nmsh, nlbc,
     *    rhsgiv, ntol, ltol, tol,
     *    xx, nudim, u, def,
     *    delu, rhs, fval,
     *    utrial, rhstri,
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,
     *    ajac, topblk, botblk, bhold, chold, ipvblk,
     *    fsub, dfsub, gsub, dgsub, itnwt, iflnwt,rpar,ipar,
     *    frscal)

      endif
*
*  these flags are used in the mesh selection strategy
*



      if (iflnwt .eq. 0) then

c       COMPUTE ESTIMATIONS OF CONDITIONING NUMBERS: norms of inverse
c       jacobian matrix
c       BY BRUGNANO & TRIGIANTE, AND HIGHAM

        N =nmsh*ncomp
        ninter=nmsh-1
        if (comp_c) then
          gamma1old = gamma1
          ckappa1old = ckappa1
          ckappaold = ckappa
          sigmaold=sigma
            call CONDESTIM(aleft,aright,nmsh,ncomp,N,xx,topblk,nlbc,
     *       ncomp, ajac, ncomp,2*ncomp,ninter,botblk,ncomp-nlbc,
     *   ipvblk,isign,amg,c1,wrkrhs,ckappa1,gamma1,sigma,ckappa,ckappa2)


           if (iprint .ge. 0) then
            CALL Rprintd1('stiffness = ', sigma)
            CALL Rprintd1('gamma1    = ', gamma1)
            CALL Rprintd1('kappa1    = ', ckappa1)
            CALL Rprintd1('kappa     = ', ckappa)           
            CALL Rprintd1('kappa2    = ', ckappa2)
           end if

          stab_sigma = abs(sigmaold-sigma)/(sigma).lt.5d-2
     *      .and. sigma .lt. flmax

          stab_kappa = abs(ckappaold-ckappa)/(ckappa).lt.5d-2
     *      .and. ckappa .lt. flmax

          stab_kappa1 = abs(ckappa1old-ckappa1)/(ckappa1).lt.5d-2
     *      .and. ckappa1 .lt. flmax .and. gamma1.lt.flmax

          stab_gamma = abs(gamma1old-gamma1)/(gamma1).lt.5d-2
     *     .and. gamma1.lt.flmax .and. ckappa1 .lt. flmax

          stab_cond = stab_kappa  .and. stab_kappa1 .and. stab_gamma

c          stab_cond = stab_kappa1 .and. stab_gamma




        stiff_cond = (( (sigma .ge. 1d1 )))
        ill_cond   = ckappa2 .ge.  1d16 .and. ckappa2 .lt. flmax
        ill_cond_newt = ckappa2 .ge.  1d10 .and. ckappa2 .lt. flmax
          if (ill_cond .and. use_c) goto 2000

           if (iprint .eq. 1) then
             CALL Rprintd1('stab_sigma = ',stab_sigma)
             CALL Rprintd1('stab_kappa = ', stab_kappa)
             CALL Rprintd1('stab_kappa1 = ', stab_kappa1)
             CALL Rprintd1('stab_gamma = ', stab_gamma)
             CALL Rprintd1('stiff_cond = ', stiff_cond)
             CALL Rprinti1('ill_cond   = ', ill_cond)
           end if
        end if
c endif if (comp_c)
c
c  The subroutine dfexcl_l substitute conv4
c
         Call dfexcl_l(Ncomp, Nmsh, Xx, Nudim, U,def8,Def, Linear, Fval,
     *        Tmp, Fsub, Dfsub, Df, Ipivlu, Dhold,
     *        Ntol, Ltol, Tol, Jflag,rpar,ipar)

         if (reaft6) then
            onto6=.true.
            goto 408
         endif
      else

       if (comp_c) then
         if (iflnwt .ne. -1) then
           gamma1old = gamma1
           ckappa1old = ckappa1
           ckappaold = ckappa
           sigmaold=sigma
           N =nmsh*ncomp
           ninter=nmsh-1
          call CONDESTIM(aleft,aright,nmsh,ncomp,N,xx,topblk,nlbc,
     *       ncomp, ajac, ncomp,2*ncomp,ninter,botblk,ncomp-nlbc,
     *   ipvblk,isign,amg,c1,wrkrhs,ckappa1,gamma1,sigma,ckappa,ckappa2)




           if (iprint .ge. 0) then
            CALL Rprintd1('stiffness = ', sigma)
            CALL Rprintd1('gamma1    = ', gamma1)
            CALL Rprintd1('kappa1    = ', ckappa1)
            CALL Rprintd1('kappa     = ', ckappa)           
            CALL Rprintd1('kappa2    = ', ckappa2)
           end if

          stab_sigma = abs(sigmaold-sigma)/(sigma).lt.5d-2
     *      .and. sigma .lt. flmax

         stab_kappa = abs(ckappaold-ckappa)/(ckappa).lt.5d-2
     *      .and. ckappa .lt. flmax

         stab_kappa1 = abs(ckappa1old-ckappa1)/(ckappa1).lt.5d-2
     *      .and. ckappa1 .lt. flmax .and. gamma1.lt.flmax

         stab_gamma = abs(gamma1old-gamma1)/(gamma1).lt.5d-2
     *     .and. gamma1.lt.flmax .and. ckappa1 .lt. flmax

         stab_cond = stab_kappa  .and. stab_kappa1 .and. stab_gamma
c         stab_cond = stab_kappa1 .and. stab_gamma

         stiff_cond = (( (sigma .ge. 1.0d1  )))
         ill_cond   = ckappa2 .ge.  1d16 .and. ckappa2 .lt. flmax
         ill_cond_newt = ckappa2 .ge.  1d10 .and. ckappa2 .lt. flmax


         if (ill_cond .and. use_c) goto 2000


           if (iprint .eq. 1) then
             CALL Rprintd1('stab_sigma = ',stab_sigma)
             CALL Rprintd1('stab_kappa = ', stab_kappa)
             CALL Rprintd1('stab_kappa1 = ', stab_kappa1)
             CALL Rprintd1('stab_gamma = ', stab_gamma)
             CALL Rprintd1('stiff_cond = ', stiff_cond)
             CALL Rprinti1('ill_cond   = ', ill_cond)
           end if
        end if
       end if
c end if if(comp_c)


         succes = .false.
         onto6 = .false.
         reaft6 = .false.
         nfail4 = nfail4+1
C karline: added , after ipar
         call fail4( ncomp, nmsh, nlbc, ntol, ltol,
     *           xx, nudim, u, rhs, linear, nmax,
     *           nmold, xxold, uold, ratdc,
     *           iorder, iflnwt, itnwt, ddouble, maxmsh,
     *           numbig, nummed,wrkrhs,amg,stab_cond,stiff_cond,
     *           ill_cond_newt,nfail4,
     *           nfxpnt, fixpnt, irefin,itcond,itcondmax,rpar,ipar,
     *           nmguess,xguess,nygdim,yguess)


c note: ratdc in subroutines fail4, fail6 does not related to
c       ratdc=dfexmx/defimp , they only use the storage
         goto 400

      endif


      If(Jflag.eq.1) then


c      if (nodouble .and. .not. forcedouble) then
c         call selcondmsh(ncomp, nmsh,
c    *     nfxpnt, fixpnt,  nmax, xx,  irefin,
c    *     nmold, xxold, ddouble, maxmsh,r4,amg)
c          ddouble = .false.
c          itcond=itcond+1
c       else
         Call Dblmsh(Nmsh, Nmax, Xx, Nmold, Xxold, Maxmsh)
          itcond = 0
c       endif
         Call Matcop(Nudim, Ncomp, Ncomp, Nmold, U, Uold)
         call interp(ncomp, nmsh, xx, nudim, u,ncomp,
     *                        nmold, xxold, uold)

         If(Maxmsh) goto 900
         goto 400
      Else

       nodouble = ( (stiff_cond .and. .not. stab_cond)
     * .and. (use_c))

      forcedouble = .false.

       if (use_c .and. itcond .eq. itcondmax) then
             itcond = 0

       endif
c    find where biggest deferred correction is
c        If(Chstif) Then
           smaldef=1.0D+40
           Bigdef=0.0D+0
           Icmph = 1
           Ix = 1
           Do 405 Iv=1,Nmsh-1
           Do 405 Iu = 1,Ntol
             Ipoint = Ltol(Iu)
             Holdef=Abs(Def8(Ipoint,Iv))
             if(smaldef.gt.holdef) smaldef=holdef
             If(Holdef.gt.Bigdef) Then
              Bigdef = Holdef
              Icmph=Ipoint
              Ixx = Iv
              intol = Iu
             Endif
 405       Continue
c   Biggest deferred correction is in component Icmph and
c   at the mesh interval Ix.
c   Now compute an explicit deferred correction for this.
         Call expl(Ncomp,Nmsh,xx,Nudim,U,dgtm,fval,fsub,Ixx,rpar,ipar)

           Ix=Ixx
           siz=abs(dgtm(Icmph))
           rat=siz/bigdef

          If(Rat.gt.50.0D+0.and.bigdef.gt.dsqrt(tol(Icmph)).and.
     +      siz.gt.1.0D+0) Then

        if (use_c) then
         if ((stiff_cond)) then
             drat = bigdef/
     *               (max(one, abs(u(icmph,Ix)))*tol(intol))

c            numadd = drat**power
             numadd = 15
            call smpselcondmsh(ncomp, nmsh,
     *        nfxpnt, fixpnt,  nmax, xx,  irefin,Ix,numadd,
     *        nmold, xxold, ddouble, maxmsh,r4,amg)
              itcond=itcond+1
         else
             drat = bigdef/
     *               (max(one, abs(u(icmph,Ix)))*tol(intol))

c            numadd = drat**power
             numadd = 15
                call smpmsh (nmsh, nmax, xx, ix, numadd,
     *             nmold, xxold, maxmsh)
         endif
       else
             drat = bigdef/
     *               (max(one, abs(u(icmph,Ix)))*tol(intol))
            numadd = 15
c            numadd = drat**power
            call smpmsh (nmsh, nmax, xx, Ix, numadd,
     *             nmold, xxold, maxmsh)
       end if

cf the solution is  interpolated, using the new  solution
         Call Matcop(Nudim, Ncomp, Ncomp, Nmold, U, Uold)
         Call Interp(Ncomp, Nmsh, Xx, Nudim,U,Ncomp,Nmold, Xxold, Uold)
         If(Maxmsh) goto 900

             goto 400
           else
c             Chstif = .false.
             onto6=.true.
             if (linear.and.ddouble) reposs=.true.
           endif
c        Endif
c endif of Chstif=true

c       Iorder = 6
c       if (iflnwt.ne.0) then
c         call mshref(ncomp, nmsh, nlbc, ntol, ltol,
c     *             iorder, rhs, ratdc,
c     *             nmax, xx, nmold, xxold, ddouble, maxmsh,
c     *             numbig, nummed,
c     *             amg,stab_cond,stiff_cond,
c     *             r4, nfxpnt,fixpnt, irefin,itcond,itcondmax)
c       elseif((.not.onto6)) then
c         Call Matcop(Nudim, Ncomp, Ncomp, Nmold, U, Uold)
c         Call Interp(Ncomp, Nmsh, Xx, Nudim, U, Nmold, Xxold, Uold)
c         goto 400
c       endif
       Endif
 408     continue
      Call Matcop(Ncomp, Ncomp, Ncomp, Nmsh-1, Def, Def6)

c      if (succes) then
c          iflbvp = 0
c          return
      if (maxmsh) then
          go to 900
      elseif (.not. onto6)  then
          go to 400
       endif

*  To reach here, onto6 must be .true.

**** logic for 6th order ****

c karline removed

*  Save the 4th order solution on this mesh in uold.

      call matcop(nudim, ncomp, ncomp, nmsh, u, uold)
*  Copy the current mesh into the xxold array.
c      nmold = nmsh
c      call dcopy(nmold, xx, 1, xxold, 1)

      iorder = 6

      if (linear) then
         call lineq( ncomp, nmsh, nlbc,
     *    ludone, xx, nudim, u, def,
     *    delu, rhs, fval,
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,
     *    ajac, topblk, botblk, bhold, chold, ipvblk,
     *    fsub, dfsub, gsub, dgsub, iflnwt,rpar,ipar)

      else

         call fixjac(ncomp, nmsh, nlbc,
     *     iorder, ntol, ltol, tol,
     *     xx, nudim, u, def, def, delu,
     *     rhs, fval, utrial, rhstri,
     *     rnsq, uint, ftmp, tmprhs,
     *     ajac, topblk, botblk, ipvblk,
     *     fsub, gsub, iflnwt,rpar,ipar)

*  If the fixed Jacobian iterations fail but rnsq is small,
*  try a Newton procedure.  Set rhsgiv to indicate that
*  the right-hand side and fval have already been evaluated
*  at the given u.

         if (iflnwt .eq. -3 .and. rnsq .lt. fxfct*epsmch    ) then
            rhsgiv = .true.

            call newteq(ncomp, nmsh, nlbc,
     *           rhsgiv, ntol, ltol, tol,
     *           xx, nudim, u, def,
     *           delu, rhs, fval,
     *           utrial, rhstri,
     *           uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,
     *           ajac, topblk, botblk, bhold, chold, ipvblk,
     *           fsub, dfsub, gsub, dgsub, iter, iflnwt,rpar,ipar,
     *              frscal)

         endif
      endif

      if (iflnwt .eq. 0) then

C karline: what about 'linear' -> rmoved
c         call conv6(ncomp, nmsh, ntol, ltol, tol,
c     *             nudim, u, uold, etest6, err6,
c     *             trst6, onto8, reaft6, linear, succes)
         call conv6(ncomp, nmsh, ntol, ltol, tol,
     *             nudim, u, uold, etest6, err6,
     *             trst6, onto8, reaft6, succes)

      else

         onto8 = .false.

         call fail6( ncomp, nmsh, nlbc, ntol, ltol, tol,
     *              nfxpnt, fixpnt, iorder, nmax,
     *              xx, nudim, u, rhs, usave, xxold, uold, nmold,
     *              ihcomp, irefin,
     *              rerr, ermx, ratdc,
     *              reaft6, ddouble, succes, maxmsh,
     *              numbig, nummed,
     *              wrkrhs,amg, stab_cond,ckappa1,gamma1,ckappa,
     *              stiff_cond,itcond, itcondmax)



      endif
      If(Maxmsh) goto 900
      if (succes) then

         if (iprint .ne. -1 .and. comp_c .and. use_c) then
           if (ill_cond) then
             CALL Rprint('The problem is ill-conditioned, ')
             CALL Rprint('The solution could be inaccurate')
           end if
         end if
         iflbvp = 0
         if (comp_c .and. use_c) then
            if (.not. stab_kappa .and. indnms .gt. 1) then
                 iflbvp = -1
                if  (iprint .ne. -1) then
        CALL Rprint('The conditioning parameters did not stabilise, ')
        CALL Rprint('The solution could be inaccurate')
                end if
            end if
            if (ill_cond) iflbvp = -2
         end if
         return
       elseif (.not. onto8) then
         go to 400
       endif

***** logic for trying to calculate 8th order solution *****

      if (iprint .eq. 1)  then
             CALL Rprint('Start 8th order')
      end if

      call matcop(nudim, ncomp, ncomp, nmsh, u, uold)
*  Copy the current mesh into the xxold array.
c      nmold = nmsh
c      call dcopy(nmold, xx, 1, xxold, 1)


*  Save the old deferred correction vector def in def6.

      call matcop(ncomp, ncomp, ncomp, nmsh-1, def, def6)

*  For linear problems, calculate the fval array for the
*  new solution u.

      if (linear) call fneval(ncomp, nmsh, xx, nudim, u, fval, fsub,
     *       rpar,ipar)

*  Calculate 8th order deferred corrections (the array def8).

      Call df8cal_l (Ncomp, Nmsh, Xx, Nudim, U, Fval, Def8, Linear,
     *             Tmp, Fsub, Dfsub, Df, Ipivlu, Dhold,
     *             Ntol, Ltol, Tol,JC,rpar,ipar)

      If(Jc.eq.1) then

      if (use_c) then
       nodouble = ((stiff_cond) .and. .not.stab_cond )
      else
        nodouble =.false.
      end if
      forcedouble = .false.

       if (iprint .eq. 1) then
         CALL Rprinti2('Forcedouble, itcond', nodouble, itcond)
       end if

       if (use_c .and.  itcond .eq. itcondmax) then
            itcond = 0
            forcedouble = .true.
       endif

       if (nodouble .and. .not. forcedouble) then
          call selcondmsh(ncomp, nmsh,
     *     nfxpnt, fixpnt,  nmax, xx,  irefin,
     *     nmold, xxold, ddouble, maxmsh,r4,amg)
           ddouble = .false.
           itcond=itcond+1
        else
         Call Dblmsh(Nmsh, Nmax, Xx, Nmold, Xxold, Maxmsh)
          itcond = 0
        endif
         If(Maxmsh) goto 900
         Call Matcop(Nudim, Ncomp, Ncomp, Nmold, Uold, U)
         Call Interp(Ncomp, Nmsh, Xx, Nudim, U,Ncomp,Nmold, Xxold, Uold)
        goto 400
      endif


*  For linear problems, the def array is the def8 array.
*  For nonlinear problems, add the def8 array to the
*  already-calculated def array.

      if (linear) then
         call matcop(ncomp, ncomp, ncomp, nmsh-1, def8, def)
      else
         call maxpy(ncomp, nmsh-1, one, def8, ncomp, def)

      endif

      iorder = 8

      if (linear) then
         call lineq( ncomp, nmsh, nlbc,
     *    ludone, xx, nudim, u, def,
     *    delu, rhs, fval,
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,
     *    ajac, topblk, botblk, bhold, chold, ipvblk,
     *    fsub, dfsub, gsub, dgsub, iflnwt,rpar,ipar)

      else

         call fixjac(ncomp, nmsh, nlbc,
     *     iorder, ntol, ltol, tol,
     *     xx, nudim, u, def, def8, delu,
     *     rhs, fval, utrial, rhstri,
     *     rnsq, uint, ftmp, tmprhs,
     *     ajac, topblk, botblk, ipvblk,
     *     fsub, gsub, iflnwt,rpar,ipar)

*  If the fixed Jacobian iterations fail but rnsq is small,
*  try a Newton procedure.  Set rhsgiv to indicate that
*  the right-hand side and fval have already been evaluated
*  at the given u.

         if (iflnwt .eq. -3 .and. rnsq .lt. fxfct*epsmch ) then
            rhsgiv = .true.


            call newteq(ncomp, nmsh, nlbc,
     *           rhsgiv, ntol, ltol, tol,
     *           xx, nudim, u, def,
     *           delu, rhs, fval,
     *           utrial, rhstri,
     *           uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,
     *           ajac, topblk, botblk, bhold, chold, ipvblk,
     *           fsub, dfsub, gsub, dgsub, iter, iflnwt,rpar,ipar,
     *            frscal)

         endif
      endif

      if (iflnwt .eq. 0) then

         errok=.false.

         call conv8( ncomp, nmsh, ntol, ltol, tol,
     *              nfxpnt, fixpnt, linear, nmax,
     *              xx, nudim, u, def, def6, def8, uold,
     *              ihcomp, irefin, ermx, err6,
     *              etest8, strctr,
     *              ddouble, nmold, xxold, maxmsh, succes, first8,
     *              wrkrhs,amg, stab_cond,ckappa1,gamma1,ckappa,
     *              stiff_cond,rpar,ipar,nmguess, xguess,nygdim, yguess)



      else

         succes = .false.
         call  fail8(ncomp, nmsh, nfxpnt, fixpnt, nmax,
     *      ntol, ltol, tol, nmold,
     *      xx, nudim, u, def6, xxold, uold,
     *      ihcomp, irefin, ermx, ddouble, maxmsh,
     *       wrkrhs,amg, stiff_cond, stab_cond)

      endif

      if (maxmsh) then
         go to 900
      elseif (.not. succes) then
         go to 400
      endif

*  Successful termination.


      if (iprint .ne. -1 .and. comp_c .and. use_c ) then
        if ( ill_cond ) then
             CALL Rprint('The problem is ill-conditioned, ')
             CALL Rprint('The solution could be inaccurate')
        end if
      end if
       iflbvp = 0
       if (comp_c .and. use_c) then
            if (.not. stab_kappa .and. indnms.gt.1) then
                iflbvp = -1
                if  (iprint .ne. -1) then
       CALL Rprint('The conditioning parameters did not stabilise, ')
       CALL Rprint('The solution could be inaccurate')
                end if
             end if
             if (ill_cond) iflbvp=-2
        end if



      return

 900   continue



      nmsh = nmold
      call dcopy(nmsh, xxold, 1, xx, 1)
c      call matcop(nudim, ncomp, ncomp, nmsh, uold, u)
*  Copy the current mesh into the xxold array.
c      nmold = nmsh
c      call dcopy(nmold, xx, 1, xxold, 1)

* Error exit---too many mesh points.

      iflbvp = 1
      if (iprint .ge. 0   ) then
             CALL Rprint('Terminated, too many mesh points')
      end if
      if (nmsh .gt. nmold .and. use_c .and.
     *    abs(ckappaold-ckappa)/(ckappa) .lt. nmsh*1d-16 ) then
      	 CALL Rprint('Try with a less stringent tolerance ')
      end if
        if (linear .and. ill_cond .and. use_c) then
             CALL Rprint('The problem is ill-conditioned, ')
             CALL Rprint('Try with a less stringent tolerance')
        end if
        if (.not.linear .and. ill_cond .and. use_c ) then
             CALL Rprint('The problem is ill-conditioned, ')
             CALL Rprint('Try with a less stringent tolerance ')
             CALL Rprint(' or with a different initial guess' )
         end if

      return
 1900   continue

* Error exit---too many meshes  .return

      iflbvp = 2
      if (iprint .ne.-1) then
       CALL Rprinti1('Terminated too many meshes, nmsh ',indnms)
      end if

      if (iprint .ne. -1 .and. use_c ) then
        if (linear .and. ill_cond ) then
             CALL Rprint('The problem is ill-conditioned, ')
             CALL Rprint('Try with a less stringent tolerance')
        end if
        if (.not.linear .and. ill_cond )then
             CALL Rprint('The problem is ill-conditioned, ')
             CALL Rprint('Try with a less stringent tolerance ')
             CALL Rprint(' or with a different initial guess' )
         end if
      end if
        return
 2000  continue

* Error exit -- ill_cond
         iflbvp = 3
         if (iprint .ne. -1) then
          if (linear) then
             CALL Rprint('The problem is ill-conditioned, ')
             CALL Rprint('Try with a less stringent tolerance')
          else
             CALL Rprint('The problem is ill-conditioned, ')
             CALL Rprint('Try with a less stringent tolerance ')
             CALL Rprint(' or with a different initial guess' )
          end if
         end if
         return
      end


c ===================================================================================

      Subroutine dfexcl_l (Ncomp, Nmsh, Xx, Nudim, U, Def8,Def6, Linear,
     *                 Fval, Tmp, Fsub, Dfsub, Df, Ip, Dhold,
     *                 Ntol, Ltol, Tol,JC,rpar,ipar)

      Implicit Double Precision (A-H, O-Z)
      Integer Ncomp, Nmsh
      dimension rpar(*),ipar(*)
      Dimension Xx(Nmsh), U(Nudim,Nmsh), Fval(Ncomp, Nmsh)
      Dimension Def6(Ncomp,Nmsh-1), Def8(Ncomp,Nmsh-1),Tmp(Ncomp,*)
      Dimension Df(Ncomp,Ncomp), Ltol(Ntol), Tol(Ntol)
      Dimension Ip(2*Ncomp), Dhold(2*Ncomp,2*Ncomp)
      External Fsub
      external Dfsub

      Parameter ( One = 1.0d+0, Two = 2.0d+0 )
C Karline: moved this statement upward
      integer nminit, iprint, idum
      Logical Linear,use_c,comp_c

      Common /Cons1/ A21,A22,A23,A24,A31,A32,A33,A34,C1,C2,
     +       C16,C26,C123,C223,C14,C24
      Common /Algprs/nminit,Iprint,idum,use_c,comp_c
      Common /Flags/ Ifinal,Iback,Iprec
      integer nfunc, njac, nstep, nbound, njacbound
      common/diagnost/nfunc, njac, nstep, nbound, njacbound


c       St1 --> Tmp(ncomp,10)
c       St2 --> Tmp(ncomp,11)
c       St3 --> Tmp(ncomp,12)
*  Given The Nmsh Mesh Points Xx, The Estimated Solution
*  U And The Array Fval Of Function Values At (Xx(Im), U(*,Im)),
*  Im = 1,...,Nmsh, dfexcl_l Calculates Sixth-Order Implicit
*  Deferred Correction, Stored In The Array Def6, Indexed
*  Over The Components And Mesh Intervals.

*  The Array Tmp Is Workspace For 4 Intermediate Vectors Of
*  Dimension Ncomp.
      nstep = nstep + 1
      Do 90 Im = 1, Nmsh-1
         Hmsh = Xx(Im+1) - Xx(Im)
         C16h=C16/Hmsh
         C26h=C26/Hmsh
         Do 10 Ic = 1, Ncomp
            Fvim = Fval(Ic,Im)
            Fvim1 = Fval(Ic,Im+1)
            Uim = U(Ic,Im)
            Uim1 = U(Ic,Im+1)
            Tmp(Ic,3) = C16h*(Uim1-Uim)+C123*Fvim1+C14*Fvim
            Tmp(Ic,4) = C26h*(Uim1-Uim)+C223*Fvim1+C24*Fvim
c
c     Put cubic Hermite approximations to the unknowns in
c     Tmp(Ic,3) and Tmp(Ic,4).
c
            Tmp(Ic,10) = (Uim+Uim1)/Two
            Tmp(Ic,11) = A21*Fvim + A24*Fvim1
            Tmp(Ic,12) = A31*Fvim + A34*Fvim1
 10                        Continue

         Xxc1 = Xx(Im)+C1*Hmsh
         Xxc2 = Xx(Im)+C2*Hmsh

         Do 60 Nit = 1, 10

         Do 20 Ic = 1,Ncomp
            Tmp3 = Tmp(Ic,3)
            Tmp4 = Tmp(Ic,4)
            Tmp(Ic,1)  = Tmp(Ic,10) +
     +                 Hmsh*(Tmp(Ic,11) + A22*Tmp3 + A23*Tmp4)
            Tmp(Ic,2)  = Tmp(Ic,10) +
     +                 Hmsh*(Tmp(Ic,12) + A32*Tmp3 + A33*Tmp4)
 20                        Continue

         Call Fsub (Ncomp,Xxc1,Tmp(1,1),Tmp(1,5),rpar,ipar)
         Call Fsub (Ncomp,Xxc2,Tmp(1,2),Tmp(1,6),rpar,ipar)
         nfunc = nfunc + 2

         Call Dfsub (Ncomp,Xxc1,Tmp(1,1),Df(1,1),rpar,ipar)
         Do 30 I = 1, Ncomp
            Tmp(I,5) = Tmp(I,5)-Tmp(I,3)
            Tmp(I,6) = Tmp(I,6)-Tmp(I,4)
            Do 25 J = 1,Ncomp
               Dfij = Hmsh*Df(I,J)
               Dhold(I,J) = -A22*Dfij
               Dhold(I,J+Ncomp) = -A23*Dfij
 25         Continue
 30      Continue

         Call Dfsub(Ncomp,Xxc2,Tmp(1,2),Df,rpar,ipar)
         Do 35 I = 1, Ncomp
            Do 32 J = 1, Ncomp
               Dfij = Hmsh*Df(I,J)
               Dhold(I+Ncomp,J) = -A32*Dfij
               Dhold(I+Ncomp,J+Ncomp) = -A33*Dfij
 32         Continue
 35      Continue
         njac  = njac + 2

         Do 40 I = 1,Ncomp
           Dhold(I,I) = Dhold(I,I) + One
           Dhold(I+Ncomp,I+Ncomp) = Dhold(I+Ncomp,I+Ncomp) + One
 40      Continue

         Call Lufac(2*Ncomp,2*Ncomp,Dhold,Ip,Ier)
         Call Lusol(2*Ncomp,2*Ncomp,Dhold,Ip,Tmp(1,5),Tmp(1,7))

         Do 45 I = 1,Ncomp
            Tmp(I,3) = Tmp(I,3) + Tmp(I,7)
            Tmp(I,4) = Tmp(I,4) + Tmp(I,8)
 45      Continue
         Jc = 0
         If (Linear) Goto 70


         Do 50 I = 1, Ntol
           Ii = Ltol(I)
           Er = Tol(I)/Hmsh
           If (Abs(Tmp(Ii,7)) .gt. Er*Max(One,Abs(Tmp(Ii,3)))  .or.
     *         Abs(Tmp(Ii,8)) .gt. Er*Max(One,Abs(Tmp(Ii,4)))) Jc = 1
 50      Continue

         If (Jc .eq. 0) Goto 70

 60      Continue

         if (iprint.eq.1) then
             CALL Rprint('NO convergence of corrections')
         end if

         Return


 70      Continue

         Do 80 Ic = 1, Ncomp
            Def6(Ic,Im) = (Hmsh/12.d+0)*(Fval(Ic,Im)+
     *              5.d+0*(Tmp(Ic,3)+Tmp(Ic,4))+Fval(Ic,Im+1))-
     *              U(Ic,Im+1)+U(Ic,Im)

 80      Continue
      do 85 ic=1,ncomp
      tmp(ic,5)=def6(ic,im)
      tmp(ic,6)=def6(ic,im)
 85         continue
      call lusol(2*ncomp,2*ncomp,dhold,Ip,tmp(1,5),tmp(1,7))
      do ic=1,ncomp
         def8(ic,Im)=tmp(ic,7)
      end do


 90   Continue

      Return

      End

c ===================================================================================

      subroutine expl (ncomp, nmsh, xx, nudim, u, defexp, fval,
     *                fsub,Im,rpar,ipar)

      implicit double precision (a-h, o-z)
      integer ncomp, nmsh
      dimension rpar(*),ipar(*)
      dimension xx(nmsh), u(nudim,nmsh), fval(ncomp, nmsh)
C     SCMODIFIED: increased the number of dimensions
C      dimension t1(2),t2(2),t3(2),t4(2)
      dimension t1(ncomp),t2(ncomp),t3(ncomp),t4(ncomp)
      dimension defexp(Ncomp)
      external fsub
      integer nfunc, njac, nstep, nbound, njacbound
      common/diagnost/nfunc, njac, nstep, nbound, njacbound

      logical use_c,comp_c
      integer nminit, iprint, idum

      common/algprs/ nminit,  iprint, idum, use_c,comp_c
      parameter ( half = 0.5d+0, fourth = 0.25d+0, thfrth= 0.75d+0 )
      a5=5.0D+0/32.0D+0
      b5=27.0D+0/32.0D+0
      c5=9.0D+0/64.0D+0
      d5=3.0D+0/64.0D+0
      e5=5.0D+0/24.0D+0
      f5=2.0D+0/3.0D+0
      a6=7.0D+0/90.0D+0
      b6=16.0D+0/45.0D+0
      c6=2.0D+0/15.0D+0
      hmsh = xx(im+1) - xx(im)
      Do 10 Ic=1,Ncomp
            t1(ic) = (a5*u(ic, im+1) + b5*u(ic, im))
     *         + hmsh*(c5*fval(ic,im)- d5*fval(ic,im+1))
            t2(ic) = (b5*u(ic,im+1) + a5*u(ic,im))
     *         + hmsh*(-c5*fval(ic,im+1) + d5*fval(ic,im))

 10   continue
      call fsub (ncomp, xx(im) + fourth*hmsh, t1,
     *          t3,rpar,ipar)
      call fsub (ncomp, xx(im) + thfrth*hmsh, t2,
     *          t4,rpar,ipar)

      Do 20 Ic=1,Ncomp
            t1(ic) = half*(u(ic,im+1) + u(ic,im))
     *          + e5*hmsh*(fval(ic,im+1) - fval(ic,im))
     *          - f5*hmsh*(t4(ic) - t3(ic))
 20   continue

      call fsub(ncomp, half*(xx(im) + xx(im+1)), t1,
     *          t2,rpar,ipar)
      nfunc = nfunc+3

      do 30 Ic=1,Ncomp
            defexp(ic) = hmsh*(a6*(fval(ic,im+1) + fval(ic,im))
     *           + b6*(t3(ic) + t4(ic)) + c6*t2(ic))
     *           - u(ic,im+1) + u(ic,im)
      Au=Max(abs(u(Ic,Im)),abs(u(ic,Im+1)))
      au=0.0D+0
      Defexp(ic)=defexp(ic)/max(1.0D+0,Au)

 30   continue

      return
      end

c ===================================================================================

      Subroutine df8cal_l(Ncomp, Nmsh, Xx, Nudim, U, Fval, Def8, Linear,
     *                   Tmp, Fsub, Dfsub, Df, Ip, Dhold, Ntol,
     *                   Ltol, Tol,JC,rpar,ipar)

*   Given The Mesh Points Xx, The Solution U, And The Function
*   Values Fval, df8cal_l Computes Eighth-Order Deferred Corrections,
*   Which Are Stored In Def8.
*   The Array Tmp Is Workspace For 9 Intermediate Vectors.

      Implicit Double Precision (A-H, O-Z)
      Integer Ncomp, Nmsh
      dimension rpar(*),ipar(*)
      Dimension Xx(Nmsh), U(Nudim,Nmsh), Fval(Ncomp,Nmsh)
      Dimension Def8(Ncomp, Nmsh-1), Tmp(Ncomp,*)
      Dimension Df(Ncomp,Ncomp), Ltol(Ntol), Tol(Ntol)
      Dimension Ip(3*Ncomp), Dhold(3*Ncomp,3*Ncomp)

      External Fsub
      External Dfsub


      Parameter ( One = 1.0d+0, Two = 2.0d+0 )

      Common/Cons2/ A21,A22,A23,A24,A25,A31,A32,A34,A35,A41,A42,
     +      A43,A44,A45,B1,B2,B3,C1,C2,C3,C16,C26,C36,C123,C223,
     +      C323,C14,C24,C34
      logical use_c,comp_c
      integer nminit, iprint, idum

      Common /Algprs/nminit,Iprint,idum,use_c,comp_c


      Logical Linear
      integer nfunc, njac, nstep, nbound, njacbound
      common/diagnost/nfunc, njac, nstep, nbound, njacbound

      JC=0
c      St1 --> Tmp(ncomp,13)
c      St2 --> Tmp(ncomp,14)
c      St3 --> Tmp(ncomp,15)
c      St4 --> Tmp(ncomp,16)
      Do 110 Im = 1, Nmsh-1
         Hmsh = Xx(Im+1) - Xx(Im)

         C16h=C16/Hmsh
         C26h=C26/Hmsh
         C36h=C36/Hmsh

         Do 10 Ic = 1, Ncomp
             Fvim = Fval(Ic,Im)
             Fvim1 = Fval(Ic,Im+1)
             Uim = U(Ic,Im)
             Uim1 = U(Ic,Im+1)
             Tmp(Ic,4) = C16h*(Uim1-Uim)+C123*Fvim1+C14*Fvim
             Tmp(Ic,5) = C26h*(Uim1-Uim)+C223*Fvim1+C24*Fvim
             Tmp(Ic,6) = C36h*(Uim1-Uim)+C323*Fvim1+C34*Fvim
             Tmp(Ic,13) = (Uim+Uim1)/Two
             Tmp(Ic,14) = A21*Fvim + A25*Fvim1
             Tmp(Ic,15) = A31*Fvim + A35*Fvim1
             Tmp(Ic,16) = A41*Fvim + A45*Fvim1
 10                          Continue

         Xxc1 = Xx(Im)+C1*Hmsh
         Xxc2 = Xx(Im)+C2*Hmsh
         Xxc3 = Xx(Im)+C3*Hmsh

         Do 80 Nit = 1, 10

            Do 20 Ic = 1, Ncomp
              Tmp4 = Tmp(Ic,4)
              Tmp5 = Tmp(Ic,5)
              Tmp6 = Tmp(Ic,6)
              Tmp(Ic,1) = Tmp(Ic,13) + Hmsh*(Tmp(Ic,14) + A22*Tmp4
     +            + A23*Tmp5  + A24*Tmp6)
              Tmp(Ic,2) = Tmp(Ic,13) + Hmsh*(Tmp(Ic,15)
     +            + A32*Tmp4 + A34*Tmp6)
              Tmp(Ic,3) = Tmp(Ic,13) + Hmsh*(Tmp(Ic,16) + A42*Tmp4
     +            + A43*Tmp5  + A44*Tmp6)
 20                               Continue

            Call Fsub (Ncomp,Xxc1,Tmp(1,1),Tmp(1,7),rpar,ipar)
            Call Fsub (Ncomp,Xxc2,Tmp(1,2),Tmp(1,8),rpar,ipar)
            Call Fsub (Ncomp,Xxc3,Tmp(1,3),Tmp(1,9),rpar,ipar)
            nfunc = nfunc+3

            Call Dfsub(Ncomp,Xxc1,Tmp(1,1),Df,rpar,ipar)
            Do 30 I = 1, Ncomp
               Tmp(I,7) = Tmp(I,7)-Tmp(I,4)
               Tmp(I,8) = Tmp(I,8)-Tmp(I,5)
               Tmp(I,9) = Tmp(I,9)-Tmp(I,6)
               Do 25 J = 1,Ncomp
                  Dfij = Hmsh*Df(I,J)
                  Dhold(I,J) = -A22*Dfij
                  Dhold(I,J+Ncomp) = -A23*Dfij
                  Dhold(I,J+2*Ncomp) = -A24*Dfij
 25            Continue
 30         Continue

            Call Dfsub(Ncomp,Xxc2,Tmp(1,2),Df,rpar,ipar)
            Do 40 I = 1, Ncomp
                Do 35 J = 1, Ncomp
                   Dfij = Hmsh*Df(I,J)
                   Dhold(I+Ncomp,J) = -A32*Dfij
                   Dhold(I+Ncomp,J+Ncomp) = 0.d+0
                   Dhold(I+Ncomp,J+2*Ncomp) = -A34*Dfij
 35             Continue
 40         Continue

            Call Dfsub(Ncomp,Xxc3,Tmp(1,3),Df,rpar,ipar)
            Do 50 I = 1, Ncomp
                 Do 45 J =  1, Ncomp
                    Dfij = Hmsh*Df(I,J)
                    Dhold(I+2*Ncomp,J) = -A42*Dfij
                    Dhold(I+2*Ncomp,J+Ncomp) = -A43*Dfij
                    Dhold(I+2*Ncomp,J+2*Ncomp) = -A44*Dfij
 45              Continue
 50         Continue

            Do 60 I = 1, Ncomp
               Dhold(I,I) = Dhold(I,I) + One
               Dhold(I+Ncomp,I+Ncomp) = Dhold(I+Ncomp,I+Ncomp) + One
               Dhold(I+2*Ncomp,I+2*Ncomp) =
     *              Dhold(I+2*Ncomp,I+2*Ncomp) + One
 60         Continue

            Call Lufac(3*Ncomp,3*Ncomp,Dhold,Ip,Ier)
            Call Lusol(3*Ncomp,3*Ncomp,Dhold,Ip,Tmp(1,7),Tmp(1,10))

            Do 65 I = 1,Ncomp
               Tmp(I,4) = Tmp(I,4) + Tmp(I,10)
               Tmp(I,5) = Tmp(I,5) + Tmp(I,11)
               Tmp(I,6) = Tmp(I,6) + Tmp(I,12)
 65         Continue
            njac  = njac+ 3
            Jc = 0
            If (Linear) Goto 90

            Do 70 I = 1, Ntol
               Ii = Ltol(I)
               Er = Tol(I)/Hmsh
               If (Abs(Tmp(Ii,10)) .gt. Er*Max(One,Abs(Tmp(Ii,4))) .or.
     *           Abs(Tmp(Ii,11)) .gt. Er*Max(One,Abs(Tmp(Ii,5))) .or.
     *           Abs(Tmp(Ii,12)) .gt. Er*Max(One,Abs(Tmp(Ii,6)))) Jc = 1
 70                                 Continue

            If (Jc .eq. 0) Goto 90

 80         Continue

         if (iprint.eq.1) then
             CALL Rprint('NO convergence of 8th order defcors')
         end if

      Return
 90   Continue


         Do 100 Ic = 1, Ncomp
           Def8(Ic,Im) = Hmsh*(B1*(Fval(Ic,Im)+Fval(Ic,Im+1))+
     *                   B2*(Tmp(Ic,4)+Tmp(Ic,6))+B3*Tmp(Ic,5))-
     *                   U(Ic,Im+1)+U(Ic,Im)
 100     Continue

 110     Continue

      Return

      End

c ===================================================================================

      Subroutine Stcon1
      Implicit Double Precision(A-H,O-Z)
      Common/Cons1/A21,A22,A23,A24,A31,A32,A33,A34,C1,C2,
     +      C16,C26,C123,C223,C14,C24

      Parameter ( One = 1.0d+0, Two = 2.0d+0,  Three = 3.0d+0 )
      Parameter ( Four = 4.0d+0, Five = 5.0d+0,  Six = 6.0d+0 )

      Rt5 = Sqrt(5.0d+0)

      A21 = (Six + Rt5)/120.0d0
      A22 = -Rt5/120.0d0
      A23 = (-13.d0*Rt5)/120.0d0
      A24 = (-Six + Rt5)/120.d0

      A31 = (Six-Rt5)/120.0d0
      A32 = (13.0d0*Rt5)/120.0d0
      A33 = Rt5 / 120.0d0
      A34 = (-Six - Rt5)/120.d0

      C1 = (Five - Rt5)/10.0d0
      C2 = (Five + Rt5)/10.0d0

      C12 = C1*C1
      C22 = C2*C2

      C16 = Six*(C1 - C12)
      C26 = Six*(C2 - C22)

      C123 = Three*C12 - Two*C1
      C223 = Three*C22 - Two*C2

      C14 = One - Four*C1 + Three*C12
      C24 = One - Four*C2 + Three*C22

      Return
      End

c ===================================================================================

      Subroutine Stcon2
      Implicit Double Precision(A-H,O-Z)
      Common/Cons2/ A21,A22,A23,A24,A25,A31,A32,A34,A35,A41,A42,
     +      A43,A44,A45,B1,B2,B3,C1,C2,C3,C16,C26,C36,C123,C223,
     +      C323,C14,C24,C34

      Parameter ( One = 1.0d+0, Two = 2.0d+0, Three = 3.0d+0 )
      Parameter ( Four = 4.0d+0, Six = 6.0d+0 )

      Rt21 = Sqrt(21.0d+0)

      A21 = One/28.d0 + Three*Rt21/1960.d0
      A22 = -Rt21/280.d0
      A23 = -32.d0*Rt21/735.d0
      A24 = -23.d0*Rt21/840.d0
      A25 = -One/28.d0 + Three*Rt21/1960.d0

      A31 = One/64.d0
      A32 = 7.d0*Rt21/192.d0
      A34 = -7.d0*Rt21/192.d0
      A35 = -One/64.d0

      A41 = One/28.d0 - Three*Rt21/1960.d0
      A42 = 23.d0*Rt21/840.d0
      A43 = 32.d0*Rt21/735.d0
      A44 = Rt21/280.d0
      A45 = -(One/28.d0) - Three*Rt21/1960.d0

      B1 = One/20.0d0
      B2 = 49.0d0/180.d0
      B3 = 16.0d0/45.d0

      C1 = One/Two - Rt21/14.d0
      C2 = One/Two
      C3 = One/Two + Rt21/14.d0

      C12 = C1*C1
      C22 = C2*C2
      C32 = C3*C3

      C16 = Six*(C1 - C12)
      C26 = Six*(C2- C22)
      C36 = Six*(C3 - C32)

      C123 = Three*C12 - Two*C1
      C223 = Three*C22 - Two*C2
      C323 = Three*C32 - Two*C3

      C14 = One - Four*C1 + Three*C12
      C24 = One - Four*C2 + Three*C22
      C34 = One - Four*C3 + Three*C32

      Return
      End

c ===================================================================================

