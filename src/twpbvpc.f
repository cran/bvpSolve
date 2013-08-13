
c ===================================================================================
c main driver for twpbvpc, written by Jeff Cash and Francesca Mazzia
c with small adaptations to make it work with R by Karline Soetaert
c ===================================================================================

      subroutine twpbvpc(ncomp, nlbc, aleft, aright,
     *       nfxpnt, fixpnt, ntol, ltol, tol,
     *       linear, givmsh, giveu, nmsh,
     *       nxxdim, xx, nudim, u, nmax,
     *       lwrkfl, wrk, lwrkin, iwrk, precis,
     *       fsub, dfsub, gsub, dgsub,
     *       ckappa1,gamma1,sigma,ckappa,
     *       ckappa2,rpar,ipar,iflbvp,liseries,iseries,indnms,
     *       full, useC, nmguess, xguess,  nygdim, yguess, iset)
C
*     OUTPUT
*
*     IFLBVP = 0   SUCCESSFULL TERMINATION
*
*     IFLBVP = -1 (SUCCESSFULL TERMINATION BUT CONDITIONING PARAMETERS
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
* WRK workspace WRK(LWRKFL)
*      LWRKFL>=16*NCOMP+3*NCOMP*NCOMP+16*NCOMP*NXXDIM
*             +5*NCOMP*NCOMP*NXXDIM+3*NXXDIM
* IWRK workspace IWRK(LWRKIN)
*      LWRKIN>=2*NCOMP*NXXDIM+2*NXXDIM+NCOMP
*
*  The subroutine twpbvpc is intended to solve two-point boundary
*  value problems.
*
*
* References:
*  Cash, J. R.; Mazzia, F.
*   A new mesh selection algorithm, based on conditioning,
*   for two-point boundary value codes.
*   J. Comput. Appl. Math.  184  (2005), no. 2, 362--381.
*
*  CASH, J. R. AND MAZZIA, F. 2009. Conditioning and Hybrid Mesh
*   Selection Algorithms For Two Point Boundary Value Problems.
*   Scalable Computing: Practice and Experience 10, 4, 347-361.
*   Revision History
*
* revision  July 3,  2006
*   added rpar and ipar in the functions
*   DFSUB, DGSUB, FSUB, GSUB
*   changed the name of the variable double in ddouble
*
* revision   September 20, 2004
*   This is a modified version of twpbvp that uses the conditioning
*   in the mesh selection.
*
*     New subroutines not included in the old version:
*           condestim
*           estimkappa
*           moncond
*           selcondmsh
*           selconderrmsh
*           smpselcondmsh
*
*     Updates subroutines:
*           bvpsol
*           conv4
*           fail4
*           fail6
*           conv8
*           fail8
*           newteq
*           mshref
*           wtchdg
*
*     Auxiliary function not used by the old version
*           donest
*           abdnrm
*
*     The common block algprs contains two more variable
*     logical use_c, comp_c
*     common/algprs/ nminit, iprint, idum,uval0, use_c, comp_c
*
* The starting subroutine is twpbvp by J.R. Cash and M.H. Wright
*
*
c ===================================================================================
* karline: to make this code compatible with R:
* 1. change all write(6,...) -> rprint
* 2. add initu (in file twpbvpa.f)
* 3. pass precisions, in 3-valued vector 'precis', as passed form C calling routine
*         - do not use d1mach of original FORTRAN code
* 4. add argument useC for conditioning or not
* 5. add arguments xguess and yguess (used if givu=TRUE)
* 6. got rid of 'pdebug'
* 7. added iset, to contain several 'counters'
c ===================================================================================

*

      implicit double precision (a-h,o-z)
      dimension rpar(*),ipar(*), iset(*), precis(3)
      dimension fixpnt(*), ltol(*), tol(*)
      dimension xx(*), u(nudim,*), xguess(*), yguess(nygdim,*)
      dimension wrk(lwrkfl), iwrk(lwrkin)
      dimension iseries(*)
      logical linear, givmsh, giveu, full, useC
      external fsub
      external dfsub
      external gsub
      external dgsub

      integer nminit, iprint, idum
      logical  use_c, comp_c, giv_u
      common/algprs/ nminit, iprint, idum, use_c, comp_c
      integer nmguess, ureset

      common/gu/ giv_u, ureset

c Karline: diagnostic properties
      integer nfunc, njac, nstep, nbound, njacbound
      common/diagnost/nfunc, njac, nstep, nbound, njacbound
      intrinsic abs
      intrinsic min

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
   50    continue
      endif

      if (ntol .lt. 1) return
      do 60 i = 1, ntol
         if (ltol(i) .lt. 0 .or. ltol(i) .gt. ncomp) return
         if (tol(i) .le. zero) return
   60 continue

      if (giveu .and. .not. givmsh) return

      if (use_c .and. .not. comp_c) return

      if (nudim .le. 0) return
      if (lwrkfl .le. 0 .or. lwrkin .le. 0) return


*  Calculate maximum number of mesh points possible with the
*  given floating-point and integer workspace.

      isp = lwrkfl  - 2*ntol - 14*ncomp - 3*ncomp*ncomp
      iden = 5*ncomp*ncomp + 16*ncomp + 3
      nmax1 = isp/iden

      isp = lwrkin - ncomp
      nmax2 = isp/(2*ncomp+2)

      nmax = min(nmax1,nmax2)
* nmax from workspace
      nmax = min(nmax, nxxdim)
* nmax from size of u and xx

      if (iprint .ge. 0) THEN
       CALL Rprinti1('Nmax from workspace =',  nmax)
      ENDIF

      if (nmax .le. 1) return


*  Partition floating point workspace.

      irhs = 1
      lrhs = ncomp*nmax
* 1 ncomp*nmax
      itpblk = irhs + lrhs
      ltpblk = ncomp*nlbc
* 2 ncomp*nmax
      ibtblk = itpblk + ltpblk
      lbtblk = ncomp*(ncomp - Nlbc)
* 3 ncomp*nmax
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
* 4 ncomp*nmax
      idef = ifval + lfval
      ldef = ncomp*(nmax-1)
* 5 ncomp*nmax
      idefex = idef + ldef
      ldefex = ncomp*(nmax-1)
* 6 ncomp*nmax
*  def6 uses the same space as defexp

      idef6 = idefex
      ldef6 = ncomp*(nmax-1)

      idefim = idef6 + ldef6
      ldefim = ncomp*(nmax-1)
* 7 ncomp*nmax
*  def8 uses the same space as defimp

      idef8 = idefim
      ldef8 = ncomp*(nmax-1)

      iusve = idef8 + ldef8
      lusve = ncomp*nmax
* 8 ncomp*nmax
      iuold = iusve + lusve
      luold = ncomp*nmax
* 9 ncomp*nmax
      itmrhs = iuold + luold
      ltmrhs = ncomp*nmax
* 10 ncomp*nmax
      irhtri = itmrhs + ltmrhs
      lrhtri = ncomp*nmax
* 11 ncomp*nmax
      idelu = irhtri + lrhtri
      ldelu = ncomp*nmax
* 12 ncomp*nmax
      ixmer = idelu + ldelu
      lxmer = ncomp*nmax

*  rerr occupies the same space as xmerit
      irerr = ixmer
      lrerr = ncomp*nmax
* 13 ncomp*nmax
      iutri = irerr + lrerr
      lutri = ncomp*nmax
* 14 ncomp*nmax
      iermx = iutri + lutri
      lermx = nmax
* 1 nmax
      irtdc = iermx + lermx
      lrtdc = nmax
* 2 nmax
      ixxold = irtdc + lrtdc
      lxxold = nmax
* 3 nmax
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
* 1 ncomp*ncomp
      idftm2 = idftm1 + ldftm1
      ldftm2 = ncomp*ncomp
* 2 ncomp*ncomp
      itmp = idftm2 + ldftm2
      ltmp = ncomp*8
* 11 ncomp
      idsq = itmp + ltmp
      ldsq = ncomp*ncomp
* 3 ncomp*ncomp
      idexr = idsq + ldsq
      ldexr = ncomp
* 12 ncomp
      ietst6 = idexr + ldexr
      letst6 = ntol
* 1 ntol
      ietst8 = ietst6 + letst6
      letst8 = ntol
* 2 ntol
      iamg = ietst8 + letst8
      lamg = ncomp*nmax
* 15 ncomp*nmax
      ic1 = iamg + lamg
      lc1 = ncomp*ncomp*nmax
* 5 ncomp*ncomp*nmax
      idelta0 = ic1 + lc1
      ldelta0 = ncomp*2
* 14 ncomp
      iwrkrhs = idelta0+ldelta0
      lwrkrhs = ncomp*nmax
* 16 ncomp*nmax
      ilast = iwrkrhs +  lwrkrhs


*  Partition integer workspace.

      iiref = 1
      liref = nmax

      iihcom = iiref + liref
      lihcom = nmax

      iipvbk = iihcom + lihcom
      lipvbk = ncomp*nmax

      iipvlu = iipvbk + lipvbk
      lipvlu = ncomp

      iisign = iipvlu + lipvlu
      lisign = ncomp*nmax

      if (iprint .eq. 1) then
        Call Rprinti1 ('Integer workspace', iisign+lisign)
      end if


c ksks: add precis as argument: machine precision...
      call bvpsol(ncomp, nmsh, nlbc, aleft, aright,
     *   nfxpnt, fixpnt, ntol, ltol, tol, nmax, linear,
     *   giveu, givmsh, xx, nudim, u,
     *   wrk(idefex), wrk(idefim), wrk(idef), wrk(idelu),
     *   wrk(irhs), wrk(ifval),
     *   wrk(itpblk), wrk(ibtblk), wrk(iajac), wrk(ibhold),
     *   wrk(ichold), iwrk(iipvbk), iwrk(iipvlu),iwrk(iisign),
     *   wrk(iuint), wrk(iftmp), wrk(itmrhs),
     *   wrk(idftm1), wrk(idftm2), wrk(idgtm),
     *   wrk(iutri), wrk(irhtri), wrk(ixmer),
     *   wrk(ixxold), wrk(iuold), wrk(iusve),
     *   wrk(itmp), wrk(idsq), wrk(idexr), wrk(irtdc),
     *   wrk(irerr), wrk(ietst6), wrk(ietst8), wrk(iermx),
     *   iwrk(iihcom), iwrk(iiref), wrk(idef6), wrk(idef8),
     *   fsub,dfsub,gsub,dgsub,iflbvp,
     *   wrk(iamg),wrk(ic1),wrk(idelta0),wrk(iwrkrhs),
     *   ckappa1,gamma1,ckappa,
     *   ckappa2, sigma,rpar, ipar,liseries,iseries,indnms,precis,
     *   nmguess,xguess, nygdim, yguess )

      iset(1) = nfunc
      iset(2) = njac
      iset(3) = nbound
      iset(4) = njacbound
      iset(5) = nstep
      iset(6) = ureset

      return
      end


c ksks: add precis as argument: machine precision...

      subroutine bvpsol(ncomp, nmsh, nlbc, aleft, aright,
     *   nfxpnt, fixpnt,
     *   ntol, ltol, tol, nmax, linear, giveu, givmsh,
     *   xx, nudim, u, defexp, defimp, def, delu, rhs, fval,
     *   topblk, botblk, ajac, bhold, chold, ipvblk, ipivlu,isign,
     *   uint, ftmp, tmprhs, dftmp1, dftmp2, dgtm,
     *   utrial, rhstri, xmerit, xxold, uold, usave,
     *   tmp, dsq, dexr, ratdc, rerr,
     *   etest6, etest8, ermx, ihcomp, irefin,
     *   def6, def8, fsub, dfsub, gsub, dgsub, iflbvp,
     *   amg, c1,delta0, wrkrhs,ckappa1,gamma1,ckappa,
     *   ckappa2,sigma,rpar,ipar,liseries,iseries,indnms,precis,
     *   nmguess,xguess, nygdim, yguess)

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
      dimension  tmp(ncomp,8)
      dimension  dsq(ncomp,ncomp), dexr(ncomp)
      dimension  ratdc(*), rerr(ncomp,*)
      dimension  etest6(*), etest8(*), ermx(*)
      dimension  ihcomp(*), irefin(*)
      dimension  def6(ncomp,*), def8(ncomp,*)
      dimension  amg(*), c1(ncomp,*), wrkrhs(*), delta0(ncomp,2)
      dimension  iseries(*)

      logical linear, giveu, givmsh, ddouble

      external fsub
      external dfsub
      external gsub
      external dgsub

      common/mchprs/flmin, flmax, epsmch

      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c
      common/monpar/ sfatt_alpha, sfatt_r3, sfatt_r1r3
      common/newt/greps
      intrinsic max

      logical smooth, succes, strctr, trst6, reaft6
      logical onto6, onto8, ludone, rhsgiv, maxmsh
      logical first4, first8

      logical mchset
      save mchset

      logical frscal
c      save frscal

c karline: added ill_cond_newt
      logical stab_kappa, stab_gamma, stab_cond, stiff_cond, ill_cond
      logical stab_kappa1, ill_cond_newt, comparekappa

      parameter (zero = 0.0d+0, one = 1.0d+0)
      parameter (third = 0.33d+0, fourth = 0.25d+0)
      parameter (quan6 = 0.1d+0 )
      parameter (itcondmax = 10)

      data mchset/.true./
      data fxfct/10.0d+0/
      data maxmsh/.false./

      frscal = .true.
      if (mchset) then
c Karline: use precis instead of d1mach

         flmin = precis(1)
         flmax = precis(2)
         epsmch = precis(3)
         mchset = .false.
      endif

*  The routine stcons calculates integration constants stored in
*  labeled common consts.

      call stcons

*  Set up arrays for the error tests.


      if (.not. linear) then
         call dload(ntol, one, etest6, 1)
      else
         do 10 i = 1, ntol
            etest6(i) = one/max(quan6, tol(i)**third)
   10    continue
      endif
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


*     initialize parameter for the conditioning estimation
      if (linear) then
         sfatt_alpha = 0.08
         sfatt_r3 = 1d-5
      else
         sfatt_alpha = 0.08
         sfatt_r3 = 1d-5
      end if

      sfatt_r1r3 = 0.5d0
      greps = 1.0d0

      itcond  = 0
      gamma1old  = flmax
      gamma1     = flmax
      ckappa1old  = flmax
      ckappa1     = flmax
      ckappa      = flmax
      ckappaold   = flmax
      ckappa2     = flmax
      stiff_cond = .false.
      stab_cond  = .false.
      ill_cond   = .false.

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
         if (nmsh .lt. nfxpnt+2) nmsh = nfxpnt + 2
         if (nmsh .ge. nmax ) then
            maxmsh = .true.
            if (iprint .eq. 1) THEN
             CALL Rprint('Initial mesh greater than nmax')
            ENDIF

            goto 900
         end if

         call unimsh(nmsh, aleft, aright, nfxpnt, fixpnt, xx)
      end if

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

cf       number of failure of Newton iteration
      nfail4=0
***** top of logic for 4th order solution ****

  400 continue
      if (iprint .eq. 1) THEN
       CALL Rprinti1('Start 4th order, nmsh', nmsh)
      ENDIF


       if (indnmsold.ne.nmsh) then
       indnms = indnms + 1
       iseries(indnms) = nmsh
       indnmsold = nmsh
       endif

       if (indnms .ge. liseries) then
          if (iprint .eq. 1) then
           CALL Rprinti1('Terminated too many meshes, nmsh', nmsh)
          end if
          goto 1900
       end if

*  Set the def (deferred correction) array to zero.

      call mtload(ncomp, nmsh-1, zero, ncomp, def)
      iorder = 4

*  The routine fneval calls fsub at the mesh xx and the
*  solution u, and saves the values in the array fval.

      call fneval(ncomp, nmsh, xx, nudim, u, fval, fsub, rpar, ipar)


*  Try to compute a 4th order solution by solving a system of nonlinear
*  equations.

      if (linear) then
         ludone = .false.

          call lineq( ncomp, nmsh, nlbc,
     *    ludone, xx, nudim, u, def,
     *    delu, rhs, fval,
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,
     *    ajac, topblk, botblk, bhold, chold, ipvblk,
     *    fsub, dfsub, gsub, dgsub, iflnwt, rpar, ipar)

*  Call fneval to evaluate the fval array at the new solution u.
*  (Such a call is not necessary for the nonlinear case because
*  fval is called within newteq for the new u.)

         call fneval(ncomp, nmsh, xx, nudim, u, fval, fsub, rpar, ipar)

      else

         rhsgiv = .false.
         call newteq(ncomp, nmsh, nlbc,
     *    rhsgiv, ntol, ltol, tol,
     *    xx, nudim, u, def,
     *    delu, rhs, fval,
     *    utrial, rhstri,
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,
     *    ajac, topblk, botblk, bhold, chold, ipvblk,
     *    fsub, dfsub, gsub, dgsub, itnwt, iflnwt, rpar, ipar, frscal)

      endif
*
*  these flags are used in the mesh selection strategy
*

      if (iflnwt .eq. 0) then


c       COMPUTE ESTIMATIONS OF CONDITION NUMBERS KAPPA1 and GAMMA1
c       related to perturbation to the boundary conditions
c

        N =nmsh*ncomp
        ninter=nmsh-1
        if (comp_c) then
          gamma1old = gamma1
          ckappa1old = ckappa1
          ckappaold = ckappa

           call CONDESTIM(aleft,aright,nmsh,ncomp,N,xx,topblk,nlbc,
     *      ncomp, ajac, ncomp,2*ncomp,ninter,botblk,ncomp-nlbc,
     *   ipvblk,isign,amg,c1,wrkrhs,ckappa1,gamma1,sigma,ckappa,ckappa2)


          if (iprint .ge. 0) then
           CALL Rprintd1('stiffness = ', sigma)
           CALL Rprintd1('gamma1    = ', gamma1)
           CALL Rprintd1('kappa1    = ', ckappa1)
           CALL Rprintd1('kappa     = ', ckappa)
           CALL Rprintd1('kappa2    = ', ckappa2)
          end if

          stab_kappa = abs(ckappaold-ckappa)/(ckappa).lt.5d-2
     *      .and. ckappa .lt. flmax

          stab_kappa1 = abs(ckappa1old-ckappa1)/(ckappa1).lt.5d-2
     *      .and. ckappa1 .lt. flmax .and. gamma1.lt.flmax

          stab_gamma = abs(gamma1old-gamma1)/(gamma1).lt.5d-2
     *     .and. gamma1.lt.flmax .and. ckappa1 .lt. flmax

          stab_cond = stab_kappa .and. stab_kappa1 .and. stab_gamma

          stiff_cond = (( (sigma .ge.  10d0  )))
          ill_cond   = (( (ckappa2 .ge.  1d016  )))

c karline: added ill_con_newt
        ill_cond_newt = ckappa2 .ge.  1d10 .and. ckappa2 .lt. flmax

           if (iprint .eq. 1) then
             CALL Rprintl1('stab_kappa = ', stab_kappa)
             CALL Rprintl1('stab_kappa1 = ', stab_kappa1)
             CALL Rprintl1('stab_gamma = ', stab_gamma)
             CALL Rprintl1('stiff_cond = ', stiff_cond)
             CALL Rprintl1('ill_cond   = ', ill_cond)
           end if
           if (ill_cond .and. use_c) goto 2000

         end if
c endif if (comp_c)


         call conv4( ncomp, nmsh, ntol, ltol, tol, linear, nmax,
     *           xx, nudim, u, defexp, defimp, def, fval,
     *           tmp, bhold, chold, dsq, dexr, usave,
     *           ratdc, rerr, ipivlu, nmold, xxold,
     *           smooth, reaft6, onto6, strctr, trst6, ddouble ,
     *           fsub, maxmsh, succes, first4,
     *           amg,stab_cond,ckappa,stiff_cond,wrkrhs,
     *           nfxpnt, fixpnt, irefin, rpar,ipar,
     *           nmguess,xguess,nygdim,yguess)


      else

       if (comp_c) then
         if (iflnwt .ne. -1) then
           gamma1old = gamma1
           ckappa1old = ckappa1
           ckappaold  =ckappa
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

          stab_kappa = abs(ckappaold-ckappa)/(ckappa).lt.5d-2
     *      .and. ckappa .lt. flmax

          stab_kappa1 = abs(ckappa1old-ckappa1)/(ckappa1).lt.5d-2
     *      .and. ckappa1 .lt. flmax .and. gamma1.lt.flmax

         stab_gamma = abs(gamma1old-gamma1)/(gamma1).lt.5d-2
     *     .and. gamma1.lt.flmax .and. ckappa1 .lt. flmax

         stab_cond = stab_kappa .and. stab_kappa1 .and. stab_gamma
         stiff_cond =  (( (sigma .ge. 5d0  )))
         ill_cond   = (( (ckappa2 .ge.  1d016  )))


         if (ill_cond .and. use_c) goto 2000
         if (iprint .eq. 1) then
             CALL Rprintl1('stab_sigma = ',stab_sigma)
             CALL Rprintl1('stab_kappa = ', stab_kappa)
             CALL Rprintl1('stab_kappa1 = ', stab_kappa1)
             CALL Rprintl1('stab_gamma = ', stab_gamma)
             CALL Rprintl1('stiff_cond = ', stiff_cond)
             CALL Rprintl1('ill_cond   = ', ill_cond)
         end if
         end if
       end if
c end if if(comp_c)

         succes = .false.
         onto6 = .false.
         reaft6 = .false.
         nfail4 = nfail4+1
c karline: added ill_cond_newt
         call fail4( ncomp, nmsh, nlbc, ntol, ltol,
     *       xx, nudim, u, rhs, linear, nmax,
     *       nmold, xxold, uold, ratdc,
     *       iorder, iflnwt, itnwt, ddouble , maxmsh,
     *       numbig,nummed,wrkrhs,amg,stab_cond,stiff_cond,
     *       ill_cond_newt,nfail4,nfxpnt,fixpnt, irefin,itcond,
     *       itcondmax,rpar,ipar,nmguess,xguess,nygdim, yguess)


      endif

      if (succes) then
          if (iprint .ne. -1 .and. comp_c .and. use_c) then
            if ( ill_cond) THEN
             CALL Rprint('The problem is ill-conditioned, ')
             CALL Rprint('The solution could be inaccurate')
            ENDIF

          end if

        iflbvp=0
        if (comp_c .and. use_c) then
            if (.not. stab_kappa .and. indnms .gt. 1) then
                iflbvp = -1
                if  (iprint .ne. -1) THEN
           CALL Rprint('The conditioning parameters did not stabilise')
           CALL Rprint('The solution could be inaccurate')
                ENDIF

             end if
        end if


          return
      elseif (maxmsh) then
          go to 900
      elseif (.not. onto6)  then
          go to 400
      endif

*  To reach here, onto6 must be .true.

**** logic for 6th order ****

      if (iprint .eq. 1) then
             CALL Rprint('Start 6th order')
      ENDIF


*  Save the 4th order solution on this mesh in uold.
      call matcop(nudim, ncomp, ncomp, nmsh, u, uold)
*  Copy the current mesh into the xxold array.
c     nmold = nmsh
c     call dcopy(nmold, xx, 1, xxold, 1)


      iorder = 6

      if (linear) then
         call lineq( ncomp, nmsh, nlbc,
     *    ludone, xx, nudim, u, def,
     *    delu, rhs, fval,
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,
     *    ajac, topblk, botblk, bhold, chold, ipvblk,
     *    fsub, dfsub, gsub, dgsub, iflnwt, rpar, ipar)

      else

         call fixjac(ncomp, nmsh, nlbc,
     *     iorder, ntol, ltol, tol,
     *     xx, nudim, u, def, def, delu,
     *     rhs, fval, utrial, rhstri,
     *     rnsq, uint, ftmp, tmprhs,
     *     ajac, topblk, botblk, ipvblk,
     *     fsub, gsub, iflnwt, rpar, ipar)

*  If the fixed Jacobian iterations fail but rnsq is small,
*  try a Newton procedure.  Set rhsgiv to indicate that
*  the right-hand side and fval have already been evaluated
*  at the given u.

         if (iflnwt .eq. -3 .and. rnsq .lt. fxfct*epsmch) then
            rhsgiv = .true.

            call newteq(ncomp, nmsh, nlbc,
     *           rhsgiv, ntol, ltol, tol,
     *           xx, nudim, u, def,
     *           delu, rhs, fval,
     *           utrial, rhstri,
     *           uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,
     *           ajac, topblk, botblk, bhold, chold, ipvblk,
     *           fsub, dfsub, gsub, dgsub, iter, iflnwt,
     *           rpar, ipar, frscal)

         endif
      endif

      if (iflnwt .eq. 0) then

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
     *              reaft6, ddouble , succes, maxmsh,
     *              numbig, nummed,
     *              wrkrhs,amg, stab_cond,ckappa1,gamma1,ckappa,
     *              stiff_cond,itcond, itcondmax)


      endif

      if (succes) then

         if (iprint .ne. -1 .and. comp_c .and. use_c) then
           if (  ill_cond ) THEN
             CALL Rprint('The problem is ill-conditioned, ')
             CALL Rprint('The solution could be inaccurate')
           ENDIF

         end if

        iflbvp=0
        if (comp_c .and. use_c) then
            if (.not. stab_kappa .and. indnms .gt. 1) then
                iflbvp = -1
                if  (iprint .ne. -1) THEN
           CALL Rprint('The conditioning parameters did not stabilise')
           CALL Rprint('The solution could be inaccurate')
                ENDIF

             end if
        end if


         return
       elseif (maxmsh) then
         go to 900
       elseif (.not. onto8) then
         go to 400
       endif

***** logic for trying to calculate 8th order solution *****

      if (iprint .eq. 1) THEN
             CALL Rprint('Start 8th order')
      ENDIF


      call matcop(nudim, ncomp, ncomp, nmsh, u, uold)
*  Copy the current mesh into the xxold array. KSKS THIS WAS TOGGLED OFF - RESET IT....
      nmold = nmsh
      call dcopy(nmold, xx, 1, xxold, 1)


*  Save the old deferred correction vector def in def6.

      call matcop(ncomp, ncomp, ncomp, nmsh-1, def, def6)

*  For linear problems, calculate the fval array for the
*  new solution u.

      if (linear) call fneval(ncomp, nmsh, xx, nudim, u, fval, fsub,
     *       rpar, ipar)

*  Calculate 8th order deferred corrections (the array def8).

      call df8cal (ncomp, nmsh, xx, nudim, u, fval, def8,
     *      tmp, fsub, rpar, ipar)

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
     *    fsub, dfsub, gsub, dgsub, iflnwt, rpar, ipar)

      else

         call fixjac(ncomp, nmsh, nlbc,
     *     iorder, ntol, ltol, tol,
     *     xx, nudim, u, def, def8, delu,
     *     rhs, fval, utrial, rhstri,
     *     rnsq, uint, ftmp, tmprhs,
     *     ajac, topblk, botblk, ipvblk,
     *     fsub, gsub, iflnwt, rpar, ipar)

*  If the fixed Jacobian iterations fail but rnsq is small,
*  try a Newton procedure.  Set rhsgiv to indicate that
*  the right-hand side and fval have already been evaluated
*  at the given u.

         if (iflnwt .eq. -3 .and. rnsq .lt. fxfct*epsmch) then
            rhsgiv = .true.

            call newteq(ncomp, nmsh, nlbc,
     *           rhsgiv, ntol, ltol, tol,
     *           xx, nudim, u, def,
     *           delu, rhs, fval,
     *           utrial, rhstri,
     *           uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,
     *           ajac, topblk, botblk, bhold, chold, ipvblk,
     *           fsub, dfsub, gsub, dgsub, iter, iflnwt,
     *           rpar, ipar, frscal)


         endif
      endif

      if (iflnwt .eq. 0) then

         call conv8( ncomp, nmsh, ntol, ltol, tol,
     *              nfxpnt, fixpnt, linear, nmax,
     *              xx, nudim, u, def, def6, def8, uold,
     *              ihcomp, irefin, ermx, err6,
     *              etest8, strctr,
     *              ddouble , nmold, xxold, maxmsh, succes, first8,
     *              wrkrhs,amg, stab_cond,ckappa1,gamma1,ckappa,
     *     stiff_cond,rpar,ipar,nmguess, xguess,nygdim, yguess)

      else

         succes = .false.
         call  fail8(ncomp, nmsh, nfxpnt, fixpnt, nmax,
     *      ntol, ltol, tol, nmold,
     *      xx, nudim, u, def6, xxold, uold,
     *      ihcomp, irefin, ermx, ddouble , maxmsh,
     *       wrkrhs,amg, stiff_cond, stab_cond)

      endif

      if (maxmsh) then
         go to 900
      elseif (.not. succes) then
         go to 400
      endif


*  Successful termination.

      if (iprint .ne. -1 .and. comp_c .and. use_c) then
        if ( ill_cond) THEN
             CALL Rprint('The problem is ill-conditioned, ')
             CALL Rprint('The solution could be inaccurate')
      ENDIF

      end if
       iflbvp=0
       if (comp_c .and. use_c) then
            if (.not. stab_kappa .and. indnms .gt. 1) then
                iflbvp = -1
                if  (iprint .ne. -1) THEN
           CALL Rprint('The conditioning parameters did not stabilise')
           CALL Rprint('The solution could be inaccurate')
                ENDIF
             end if
       end if


      return


 1900   continue

* Error exit---too many meshes  .return

      iflbvp = 2
      if (iprint .ne. -1) then
        CALL Rprinti1('Terminated too many meshes, nmsh ', indnms)
        if (linear .and. ill_cond) THEN
             CALL Rprint('The problem is ill-conditioned')
             CALL Rprint('Try with a less stringent tolerance')
       ENDIF
        if (.not.linear .and. ill_cond) THEN
      CALL Rprint('The problem is ill-conditioned, try with a less')
      CALL Rprint('stringent tolerance or with different initial guess')
        ENDIF
      end if
      return
 2000 continue

         iflbvp = 3
         if (iprint .ne. -1) then
           if (linear) then
             CALL Rprint('The problem is ill-conditioned')
             CALL Rprint('Try with a less stringent tolerance')
           else
      CALL Rprint('The problem is ill-conditioned, try with a less')
      CALL Rprint('stringent tolerance or with different initial guess')
           end if
         end if
         return
  900 continue

* Error exit---too many mesh points.

      iflbvp = 1
      if (iprint .ne. -1) then
             CALL Rprint('Terminated, too many mesh points')
      end if
      if (iprint .ne. -1 .and. use_c ) then
        if (linear .and. ill_cond) THEN
             CALL Rprint('The problem is ill-conditioned')
             CALL Rprint('Try with a less stringent tolerance')
        ENDIF

        if (.not.linear .and.ill_cond) THEN
      CALL Rprint('The problem is ill-conditioned, try with a less')
      CALL Rprint('stringent tolerance or with different initial guess')
        ENDIF

      end if
      return


      end



