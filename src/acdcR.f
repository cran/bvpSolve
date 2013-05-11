C Francesca: September 2011, changed again the exit strategy
c francesca: add precis as argument: machine precision...
c changed acinitu, added in input more information about xguess,uguess
c
c ===================================================================================
c karline: some subroutines renamed by adding ac in front
C karline: changed the exit strategy if epsmin was changed
C karline: rewrote all write statements, removed all formats
c ===================================================================================


c ===================================================================================
c acinitu resets u after re-meshing for linear problems or for nonlinear problems
c when interpolation of the old solution is not used.
c it interpolates between (Xguess,Yguess), if these are inputted
c otherwise sets to constant value...
c ===================================================================================


      subroutine acinitu(ncomp,nmsh,xx,nudim,u,
     +                       nugdim,nmguess,xguess,uguess)
      implicit double precision (a-h,o-z)
      dimension xx(*), u(nudim, *), xguess(*), uguess(nugdim,*)

      logical use_c, comp_c, giv_u
      integer nmguess, ureset
      common/algprs/ nminit, iprint, idum, use_c, comp_c
      Common/acAlgprs/Maxcon,Itsaim,Uval0
      common/acgu/ giv_u, ureset



      ureset = ureset + 1

      IF (giv_u) THEN
       if (iprint .ne. -1) then
        CALL Rprint('acinitu = 0.0')
       endif

       call acinterp(ncomp, nmsh, xx, nudim, u,
     *                  nugdim,nmguess, xguess, uguess)

      ELSE
       if (iprint .ne. -1) then
        CALL Rprintd1('acinitu', Uval0)
       endif
       call mtload(ncomp, nmsh, Uval0, nudim, u)
      ENDIF

      return
      end


c November 2010
c Francesca: subroutine acdc, corrected a lot of bugs
c        modified all the linear algebra routines
c        and added the computation of the conditioning parameters
c
c   changed the mening of IFLBVP
c
c     IFLBVP = 0   SUCCESSFULL TERMINATION
c     IFLBVP = -1  Epsmin changed, successfull termination
*     IFLBVP = 1   terminated, final problem not solved
*     IFLBVP = 2   terminated too many continuation steps
*     IFLBVP = 3   terminated ill conditioned problem
*                            (ONLY IF USE_C = .true.)
*     IFLBVP = 4   terminated invalid input
*
c added statistics in icount
c      icount(1) = total number of function evaluation
c      icount(2) = total number of jacobian evaluation
c      icount(3) = total number of boundary condition evaluation
c      icount(4) = total number of jacobian of boundary condition evaluation
c      icount(5) = total number of deferred corrections steps
c      icound(6) = total number of continuation steps
c      icount(7) = total number of successfull continuation steps
c      icount(8) = maximum mesh used
c
c     ckappa1, gamma1,sigma,ckappa,ckappa2 = conditioning parameters
c
c karline: added 'Full' to set output level
c          added 'useC' for specification of conditioning
      Subroutine acdc(Ncomp, Nlbc, Nucol, Aleft, Aright, Nfxpnt, Fixpnt,
     +            Ntol, Ltol, Tol, Linear, Givmsh, Giveu,
     +            Full,nmshguess, xguess, nugdim, uguess,Nmsh, Xx,
     +            Nudim, U, Nmax, Lwrkfl, Wrk, Lwrkin, Iwrk, Giveps,
     +            Eps, Epsmin, acfsub, acdfsub, acgsub, acdgsub,
     +            ckappa1,gamma1,sigma,ckappa,ckappa2,rpar,ipar,icount,
     +            precis, useC, Iflbvp)

      Implicit Double Precision (A-H,O-Z)
      DIMENSION RPAR(*), IPAR(*)
      Dimension Fixpnt(*), Ltol(*), Tol(*), icount(*)
      Dimension Xx(*), U(Nudim,*), Xguess(*), Uguess(Nugdim,*)
      Dimension Wrk(Lwrkfl), Iwrk(Lwrkin), precis(3)
      Dimension Phi(3), E(3), Npr(2), Pmax(2), Hord(2), StiffSig(3)

      Logical Linear, Givmsh, Giveu, Giveps, eps_changed, Full, useC
      External acfsub
      External acdfsub
      External acgsub
      External acdgsub
c Francesca: added use_c and comp_c
      LOGICAL use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c
      Common/acAlgprs/Maxcon,Itsaim,Uval0
      Common /acFlags/ Ifinal,Iatt,Iback,Iprec,Iucond
      Common /acConvg/ Nits
      Common /acMshvar/ Hsml,Npr,Pmax,Hord
c Francesca: added counters
      integer nfunc, njac, nstep, nbound, njacbound
      common /Mcoldiag/nfunc, njac, nstep, nbound, njacbound, maxmesh

      Parameter ( Zero = 0.0d+0, One = 1.0d+0, Two = 2.0d+0 )
      Parameter ( Three = 3.0d+0, Third = 1.0d+0/3.0d+0, Huge = 1.d+30 )

C.... The First Part Of This Subroutine Is Used For Checking Input
C.... Parameters And Dividing Up The Floating Point And Integer
C.... Workspaces. The Rest Of This Subroutine Is Devoted To The
C.... Continuation Algorithm.
c KS:   initialise block data
       Nminit = 7
       Iprint = 0
       Itsaim = 7
       Maxcon = 100
       Uval0 = 0.0d+0
       use_C = useC
       comp_C = .TRUE.

C     intialise counters
      nfunc = 0
      njac = 0
      nstep = 0
      nbound = 0
      njacbound = 0

C.... Output Details Of The Problem
      IF (Full) THEN
        iprint = 1
      ELSE
        iprint = -1
      ENDIF
C karline: Removed a lot of unnecessary printing

C.... Check For Invalid Input Parameters.  If Any Parameters Are
C.... Invalid, Exit With The Flag Iflbvp Set To 4

      Iflbvp = 4
      If (Ncomp .le. 0)  Return
      If (Nlbc .lt. 0 .or. Nlbc .gt. Ncomp) Return
      If (Nucol .le. 0) Return
      If (Aleft .ge. Aright) Return

      If (Nfxpnt .lt. 0)  Return
      If (Maxcon .lt. 0)  Return
      If (.not. Givmsh) Nmsh = 1

      If (Givmsh .and. Nmshguess .lt. Nfxpnt+2) Return
      If (Givmsh .and. nmshguess .gt. nucol ) Return
      If (Givmsh .and. xguess(1) .ne. Aleft) Return
      If (Givmsh .and. xguess(Nmsh) .ne. Aright) Return
      If (Nfxpnt .gt. 0) Then
         If (Fixpnt(1) .le. Aleft) Return
         If (Fixpnt(Nfxpnt) .ge. Aright) Return
         Do 10 I = 1, Nfxpnt-1
            If (Fixpnt(I+1) .le. Fixpnt(I)) Return
   10    Continue
      Endif

      If (Ntol .lt. 1) Return
      Do 20 I = 1, Ntol
         If (Ltol(I) .lt. 0 .or. Ltol(I) .gt. Ncomp) Return
         If (Tol(I) .le. Zero) Return
   20 Continue

      If (Giveu .and. .not. Givmsh) Return

      If (Nudim .le. 0) Return
      If (Lwrkfl .le. 0 .or. Lwrkin .le. 0) Return
      If (Itsaim .lt. 1) Return
      If (.not. Giveps) Eps = 0.5d0
      If (Eps .le. Zero .or. Eps .ge. One) Return
      If (Epsmin .gt. Eps .or. Epsmin .le. Zero) Return
      If (Iprint .ge. 0) then
       CALL Rprintd1('The initial value of Epsilon is ',Eps)
       CALL Rprintd1('The desired final value of Epsilon is ',Epsmin)
      end if
C.... Calculate Maximum Number Of Mesh Points Possible With The
C.... Given Floating-Point And Integer Workspace.


      isp = lwrkfl  - 2*ntol - 19*ncomp - 14*ncomp*ncomp
      iden = 5*ncomp*ncomp + 14*ncomp + 6
      nmax1 = isp/iden

      isp = lwrkin - 3*ncomp
      nmax2 = isp/(2*ncomp+2)

      nmax = min(nmax1,nmax2)
* nmax from workspace
      nmax = min(nmax, nucol)
* nmax from size of u and xx
      If (Iprint .eq. 1) then
      CALL Rprinti1('The largest mesh size from workspace, nmax =',Nmax)
      end if
      Nmax = Min(Nmax, Nucol)
      If (Iprint .eq. 1 .and. Nmax .eq. Nucol) then
        CALL Rprinti1('Nmax = nucol =', Nmax)
      end if
      If (nmsh .gt. nmax)    Return

      If (Nmax .le. 1) Return

C.... Partition Floating Point Workspace.

      Irhs = 1
      Lrhs = Ncomp*Nmax
* 1 ncomp*nmax
      Itpblk = Irhs + Lrhs
      Ltpblk = Ncomp*Nlbc
* 1 ncomp*ncomp
      Ibtblk = Itpblk + Ltpblk
      Lbtblk = Ncomp*(ncomp - Nlbc)
* 2 ncomp*ncomp
      Iajac = Ibtblk + Lbtblk
      Lajac = 2*Ncomp*Ncomp*Nmax
* 2 ncomp*ncomp*nmax
      Ibhold = Iajac + Lajac
      Lbhold = Ncomp*Ncomp*Nmax
* 3 ncomp*ncomp*nmax
      Ichold = Ibhold + Lbhold
      Lchold = Ncomp*Ncomp*Nmax
* 4 ncomp*ncomp*nmax
      Ifval = Ichold + Lchold
      Lfval = Ncomp*Nmax
* 2 ncomp*nmax
      Idef = Ifval + Lfval
      Ldef = Ncomp*(Nmax-1)
* 3 ncomp*nmax
      idef6 = idef
      ldef6 = ncomp*(nmax-1)
* 4 ncomp*nmax
      Idef8 = Idef6 + Ldef6
      Ldef8 = Ncomp*(Nmax-1)
* 5 ncomp*nmax
      Iuold = Idef8 + Ldef8
      Luold = Ncomp*Nmax
* 6 ncomp*nmax
      Itmrhs = Iuold + Luold
      Ltmrhs = Ncomp*Nmax
* 7 ncomp*nmax
      Irhtri = Itmrhs + Ltmrhs
      Lrhtri = Ncomp*Nmax
* 8 ncomp*nmax
      Idelu = Irhtri + Lrhtri
      Ldelu = Ncomp*Nmax
* 9 ncomp*nmax
      Ixmer = Idelu + Ldelu
      Lxmer = Ncomp*Nmax
* 10 ncomp*nmax
      Iutri = Ixmer + Lxmer
      Lutri = Ncomp*Nmax
* 11 ncomp*nmax
      Iermx = Iutri + Lutri
      Lermx = Nmax
* 1 nmax
      Irtdc = Iermx + Lermx
      Lrtdc = Nmax
* 2 nmax
      Ixxold = Irtdc + Lrtdc
      Lxxold = Nmax
* 3 nmax
      Iuint = Ixxold + Lxxold
      Luint = Ncomp
* 1 ncomp
      Iftmp = Iuint + Luint
      Lftmp = Ncomp
* 2 ncomp
      Idgtm = Iftmp + Lftmp
      Ldgtm = Ncomp
* 3 ncomp
      Idftm1 = Idgtm + Ldgtm
      Ldftm1 = Ncomp*Ncomp
* 3 ncomp*ncomp
      Idftm2 = Idftm1 + Ldftm1
      Ldftm2 = Ncomp*Ncomp
* 4 ncomp*ncomp
      Itmp = Idftm2 + Ldftm2
      Ltmp = Ncomp*16
* 19 ncomp
      Idsq = Itmp + Ltmp
      Ldsq = Ncomp*Ncomp
* 5 ncomp*ncomp
      Ietst6 = Idsq + Ldsq
      Letst6 = Ntol
* 1 ntol
      Ietst8 = Ietst6 + Letst6
      Letst8 = Ntol
* 2 ntol
      Iextra = Ietst8 + Letst8
      Lextra = 9*Ncomp*Ncomp
* 14 ncomp*ncomp
      Iextra1 = Iextra + Lextra
      Lextra1 = Nmax
* 4 nmax
      iamg = Iextra1 + Lextra1
      lamg = ncomp*nmax
* 12 ncomp*nmax
      ic1 = iamg + lamg
      lc1 = ncomp*ncomp*nmax
* 5 ncomp*ncomp*nmax

      iwrkrhs = ic1+lc1
      lwrkrhs = ncomp*nmax
* 13 ncomp*nmax
      ir4 = iwrkrhs + lwrkrhs
      lr4 = nmax
* 5 nmax
      Imsh =  ir4 + lr4
      Lmsh = Nmax
* 6 nmax

      If (Linear) Then
         Ilast = Imsh + Lmsh
      Else
         Isol = Imsh + Lmsh
         Lsol = Ncomp*Nmax
* 14 ncomp*nmax
         Ilast = Isol + Lsol
      Endif

C.... Partition Integer Workspace.

      Iiref = 1
      Liref = Nmax

      Iihcom = Iiref + Liref
      Lihcom = Nmax

      Iipvbk = Iihcom + Lihcom
      Lipvbk = Ncomp*Nmax

      Iipvlu = Iipvbk + Lipvbk
      Lipvlu = 3*Ncomp


      iisign = iipvlu + lipvlu
      lisign = ncomp*nmax

      if (iprint .eq. 1) then
        CALL Rprinti1('Integer workspace ', iisign+lisign)
      end if

C.... *****************************************************************
C.... Initialisation And Explanation Of Variables Used In The
C.... Continuation Strategy Begins Here.
C.... *****************************************************************

C.... Initialise Extrapolation Variables

      Do 30 J = 1,3
        Phi(J) = Zero
        E(J) = Zero
 30   Continue

      Emin = One/Epsmin
      Ep = One/Eps


C.... The Maximum Monitor Function Values On The First Two Meshes
C.... For A Given Step Are Stored In Pmax(1) And Pmax(2).

      Pmax(1) = Zero
      Pmax(2) = Zero

C.... Nits Counts The Number Of Newton Iterations Required For
C.... Convergence On The First Mesh For A Given Continuation Problem

      Nits = 0

C.... Ifinal Is A Flag Indicating Whether The Final Problem Is Being
C.... Attempted. If Iback = 1 Then We Do Not Backtrack. Iprec Is A
C.... Flag That Is Set To One When We Believe That The A Calculation
C.... Is Being Attempted Which May Be Beyond The Bounds Imposed By
C.... The Precision Of The Machine.
C.... Epsp Is The Value Of Eps Beyond Which The Machine Precision Is
C.... (Possibly) Not Sufficient To Solve The Given Problem.

      Ifinal = 0
      Iback = 0
      Iprec = 0

      Epsp = Zero
C francesca
C  maxmesh is the max mesh used
C  nstep is the total  number of deferred correction steps
      maxmesh = 0
      nstep   = 0
c  francesca eps_changed is true if epsmin is changed
      eps_changed = .false.

C.... Istep Is A Flag Which Limits The Size Of The Continuation Steps
C.... After Experiencing A Failure For Some Step.

      Istep = 0

C.... Nc Counts The Total Number Of Continuation Steps Taken.
C.... Nss Counts The Number Of Successful Continuation Steps Taken.

      Nc = 0
      Nss = 0

C.... Idc Counts The Number Of Consecutive Steps For Which The Maximum
C.... Value Of The Monitor Function Is Decreasing.
C.... Iextrap Is A Flag That Indicates Whether Monitor Function
C.... Extrapolation Is Possible.
C.... Istuk Is A Flag Which Is Set To One When The Continuation
C.... Algorithm Is Backtracking.

      Idc = 0
      Iextrap = 0
      Istuk = 0

C.... Ilin Is A Flag Indicating Whether Linear Or Quadratic
C.... Extrapolation Should Be Used.

      Ilin = -1

C.... Hrat Is A Variable Utilised When We Have A Moving Layer.

      Hrat = One

C.... The Following Four Variables Are Used In The Parameter
C.... Selection Process.

      Amax = Huge
      Bmax = Zero
      Cmax = Huge
      Estep = Huge

      Iusep = 0

C .... Initialize Xx and U

      If (.not. Giveu .and. .not. Givmsh) Then
         Nmsh = Nminit
         If (Nmsh .lt. Nfxpnt+2) Nmsh = Nfxpnt + 2
         Call Unimsh(Nmsh, Aleft, Aright, Nfxpnt, Fixpnt, Xx)
      Endif

      If (Givmsh) Then
        Nmsh = Nmshguess
        Call Dcopy(Nmsh, xguess, 1, Xx, 1)
      Endif

c      If (.not. Giveu) Then
c         Call Initu(Ncomp, Nmsh, Nudim, U)
c
c      else
c         Call Matcop(Nudim, Nudim, Ncomp, Nmsh, Uguess, U)
c      endif

      call acinitu(ncomp, nmsh, xx, nudim, u,
     +            nugdim,nmshguess,xguess,Uguess)

C.... Calculate Phimax, Which Is The Largest Monitor Function Value
C.... That We Believe Is Permissible In Double Precision For The Given
C.... Tolerance.

      Tolmax = One
      Phimax = One/D1mach(3)
      Do 40 J = 1,Ntol
        Tolmax = Min(Tol(J),Tolmax)
 40   Continue
      Phimax = Phimax*Tolmax
      Phiaim = Zero
      Phialt = Zero
      Epold = Zero
      Hsmlold = Zero
      Htpv = Zero

C.... We Predict That For Epsilon = 1 The Maximum Value Of The
C.... Monitor Function Is Minimised By Phi(3). We Will Use This Value
C.... When Performing Extrapolation.

      E(3) = One
      Phi(3) = One/(Aright-Aleft)



C.... *****************************************************************
C.... End Of Initialisation Of Continuation Variables.
C.... *****************************************************************

C.... *****************************************************************
C.... The Algorithm Loops Back To Line 50 For Every Continuation Step.
C.... *****************************************************************

 50   Eps = Max(One/Ep,Epsmin)
      Ep = One/Eps


C.... Solve For Epsmin Exactly

      If (Eps .lt. (1.00001d0*Epsmin)) Ifinal = 1

C.... Increase The Continuation Step Counter Nc.

      Nc = Nc+1
      If (Iprint .ge. 0) then
        CALL Rprinti1('Continuation step ',Nc)
        CALL Rprintd1('Epsilon = ',Eps)
      end if

C.... If We Have Reached The Maximum Number Of Continuation Steps
C.... Then This Will Be The Last Problem We Attempt.

      If (Nc .eq. Maxcon) Then
         If (Iprint .ge. 0) then
           CALL Rprint('This is the final continuation step')
         end if
         Iback = 1
         Ifinal = 1
         eps_changed=.true.
      Endif

C.... Iatt Is A Flag That Governs Whether Mesh Selection Is Performed
C.... Or Control Is Returned To The Driver To Select A New Parameter.


      Iucond = 0
      Iatt = -1
      Iflbvp = 0

C.... Attempt To Solve The Latest Continuation Problem.

      CALL acBvpsol(Ncomp, Nmsh, Nlbc, Aleft, Aright,
     *   Nfxpnt, Fixpnt, Ntol, Ltol, Tol, Nmax, Linear,
     *   Giveu, Givmsh, nmshguess, xguess,nugdim,
     *   uguess, Xx, Nudim, U,
     *   Wrk(Idef), Wrk(Idelu),
     *   Wrk(Irhs), Wrk(Ifval),
     *   Wrk(Itpblk), Wrk(Ibtblk), Wrk(Iajac), Wrk(Ibhold),
     *   Wrk(Ichold), Iwrk(Iipvbk), Iwrk(Iipvlu),
     *   Wrk(Iuint), Wrk(Iftmp), Wrk(Itmrhs),
     *   Wrk(Idftm1), Wrk(Idftm2), Wrk(Idgtm),
     *   Wrk(Iutri), Wrk(Irhtri), Wrk(Ixmer),
     *   Wrk(Ixxold), Wrk(Iuold),
     *   Wrk(Itmp), Wrk(Idsq), Wrk(Irtdc),
     *   Wrk(Ietst6), Wrk(Ietst8), Wrk(Iermx),
     *   Iwrk(Iihcom), Iwrk(Iiref), wrk(idef6), Wrk(Idef8),
     *   Wrk(Iextra), Wrk(Iextra1),
     *   acfsub, acdfsub, acgsub, acdgsub, Iflbvp,
     *   wrk(iamg),wrk(ic1),wrk(iwrkrhs),
     *   ckappa1,gamma1,sigma,ckappa,ckappa2,wrk(ir4),
     *   precis,Eps, rpar,ipar)

C.... *****************************************************************
C.... The Logic For A Successful Continuation Step Starts Here
C.... *****************************************************************

      If (Iflbvp .eq. 0) Then


         Nss = Nss + 1

C.... If The Problem Eps = Epsmin Has Been Solved Succesfully Then
C.... We May Finish. Ifinal = 1 When Eps = Epsmin.

         If (Ifinal .eq. 1) Then
            If (Iprint .ge. 0) Then
        CALL Rprinti1('Continuation step ',Nc)
        CALL Rprintd1('Epsilon = ', Eps)
        CALL Rprint('Tolerances satisfied for final problem')
        CALL Rprintd1('Epsilon = ', Eps)
C Karline: removed lot of print statements
            Endif

            if (eps_changed) iflbvp = -1
C KARLINE: CHANGED
            GOTO 100
         Endif

C.... When Bactracking The Program Should Find A Problem That It Can
C.... Eventually Solve. The Relationship Between The Last Two
C.... Successful Epsilon Values Is Used To Restrict Future Epsilon
C.... Values. This Is The Purpose Of Estep. If Istep = 1, Then We
C.... Must Restrict Future Values Of Ep.

         If (Istuk .eq. 1) Then
            Istep = 1
            Estep = Min(Ep/E(3),Estep)
         Endif
         Istuk = 0

C.... Calculate The Best Approximation To The Maximum Value Of The
C.... Monitor Function By Extrapolating The Maximum Value On The
C.... First Two Meshes. The Best Approximation Is Stored In Phit.

         If (Iextrap .eq. 0) Then
            H1 = Hord(1)/Pmax(1)
            H2 = Hord(2)/Pmax(2)
            If (H1/H2 .gt. 1.1d0) Then
               C1h = (Pmax(2) - Pmax(1))/(H2-H1)
               Phit = Pmax(1) - C1h*H1
            Else
               Phit = Pmax(2)
            Endif
         Endif

C.... It Is Possible That We Are Attempting To Solve Problems
C.... That Are Beyond The Bounds Imposed By The Machine Precision.
C.... The User Is Warned Of This Possibility.

         If (Iprec .eq. 0 .and. Phit .gt. Phimax) Then
            If (Iprint .ge. 0) then
      CALL Rprint(' ** Machine precision (possibly) not sufficient for')
        CALL Rprintd1(' Epsilon less than ',Eps)
            end if
            Epsp = Eps
            Iprec = 1
         Endif

C.... Save Details Of Last Problem Solved Successfully In Case We
C.... Need To Backtrack

         Nnewold = Nmsh
         CALL acMeshdet(Wrk(Imsh), Xx, Nmsh)
         If (.not. Linear) CALL acSoldet(Wrk(Isol), U, Ncomp, Nmsh,
     +        Ncomp, Nudim)

Cf   save details for the next step
c          Nmshguess = Nmsh
c          Call Dcopy(Nmsh, Xx, 1, Xguess, 1)
c          Call Matcop(Nudim, Nudim, Ncomp, Nmsh, U, Uguess)

C.... If Iextrap Equals 1 Then We Have Abandoned Monitor Function
C.... Extrapolation As Our Parameter Selection Basis.

         If (Iextrap .eq. 1) Then
            E(1) = E(2)
            E(2) = E(3)
            E(3) = Ep
            If (Istep .eq. 1) Then
               Ep = Estep*Ep
            Else
               Ep = Emin
            Endif
            Goto 70
         Endif

C.... If We Have A Decreasing Monitor Function Twice Consecutively
C.... Then Extrapolation Will Not Work. The Flag Iextrap Equals 1
C.... When Extrapolation Does Not Work.

         If (Phit .le. Phi(3)) Then
            If (Idc .eq. 1 .or. Phit .lt. Phi(2)) Then
               Iextrap = 1
               E(1) = E(2)
               E(2) = E(3)
               E(3) = Ep

               If (Istep .eq. 1) Then
                  Ep = Estep*Ep
               Else
                  Ep = Emin
               Endif
               Goto 70
            Else
               Idc = 1
               E(3) = Ep
               Phi(3) = Phit
            Endif

C.... Otherwise Update Extrapolation Data

         Else
            Idc = 0
            Do 60 J = 1,2
               E(J) = E(J+1)
               Phi(J) = Phi(J+1)
 60         Continue
            E(3) = Ep
            Phi(3) = Phit
         Endif

C.... The Variable Irest Is A Flag Which Informs Us Whether Our
C.... Extrapolation Procedure Is Accurate Or Not. Irest = 1 Means
C.... Extrapolation Is Not Reliable And So We Should Restrict The
C.... Next Parameter Step.

         Irest = 0

C.... Check Errors In Extrapolation Procedure And Calculate Hmult
C.... Which Restricts The Size Of Parameter Steps If Extrapolation
C.... Is Working Poorly.

         If (Nss .eq. 2 .and. Phit .gt. Two*Phiaim) Then
            Irest = 1
            Hh = Ep - Epold
            Errp = Abs(Phit - Phiaim)
            Fx = Errp/Hh**2
            Hmax = Sqrt(0.2d0*Phiaim/Fx)
            Hmult = (Epold+Hmax)/Epold
         Endif


         If (Nss .gt. 2 .and. Phit .gt. Two*Phiaim) Then
            Irest = 1

            Hh = Ep - Epold
            Errp = Abs(Phit - Phiaim)
            Errp1 = Abs(Phit - Phialt)
            If (Ilin .eq. -1) Then
               Fx = Errp/Hh**2
               Hmax = Sqrt(0.2d0*Phiaim/Fx)
               Hmult = (Epold+Hmax)/Epold
               Fx = Errp1/Hh**3
               Hmax = (0.2d0*Phiaim/Fx)**(Third)
               Hmult1 = (Epold+Hmax)/Epold
            Else
               Fx = Errp/Hh**3
               Hmax = (0.2d0*Phiaim/Fx)**(Third)
               Hmult = (Epold+Hmax)/Epold
               Fx = Errp1/Hh**2
               Hmax = Sqrt(0.2d0*Phiaim/Fx)
               Hmult1 = (Epold+Hmax)/Epold
            Endif
         Endif
         Epold = Ep

C.... If The Extrapolation Procedure Is Particularly Inaccurate Then
C.... We Need To Restrict It Permanently By Setting Istep = 1.

         Pinacc = Three*Max(Pmax(1),Phiaim,Phialt)
         If (Nss .gt. 2 .and. Phit .gt. Pinacc) Then
            Istep = 1
            Pestep = Max(Hmult,Hmult1)
            Estep = Min(Estep,Pestep)
         Endif

C.... Decide Whether Linear Or Quadratic Extrapolation Is Best.

         If (Nss .gt. 2) Then
            Per1 = Abs(Phiaim-Phit)
            Per2 = Abs(Phialt-Phit)
            If (Per2 .lt. Per1) Then
               Ilin = -Ilin
               If (Phit .gt. Two*Phialt) Then
                  Hmult = Hmult1
               Else
                  Irest = 0
               Endif
            Endif
         Endif

C.... If The Continuation Steps Are Becoming Very Small Then We Have
C.... Reached Our Final Problem

         Dele = E(3)-E(2)
         If (Dele .lt. 0.01d0*E(3) .or. (Nmsh .gt. (3*Nmax/4))) Then
            Ep = E(3)
            Epsmin = One/Ep
            If (Iprint .ge. 0) Then
               If (Dele .lt. 0.01d0*E(3)) Then
        CALL Rprint('Continuation steps too small')
        CALL Rprintd1('Change Eps to ', Epsmin)
               Else
        CALL Rprint('Storage limit being approached')
        CALL Rprintd1('Change Eps to ', Epsmin)
               Endif
            Endif
            Emin = Ep
            Iback = 1
            eps_changed=.true.
            iflbvp=-1
C KARLINE: ADDED THAT - NOW GOTO 100, NOT TO 70 which crashes R!
C            Goto 70
            GOTO 100
         Endif

C.... The Following Section Of Code Calculates The Desired Value
C.... (Phiaim) That We Would Like The Maximum Value Of The Monitor
C.... Function To Take.


         Itru = 0
         If (Phit .gt. Pmax(1) .and. Pmax(2) .ne. Phit) Itru = 1
         If (Itru .eq. 1) Amax = -Phit*Phit/(Two*C1h)
         Bmax = Phit*H1
         If (Npr(1) .gt. 3*Npr(2)/2) Then
            Ccm = Max(Bmax-Three,Bmax/1.5d0)
            Cmax = Min(Cmax,Ccm)
         Endif
         If (Itru .eq. 0) Amax = Max(1.5d0*Bmax,Bmax+Three)
         If (Nss .ge. 2) Hrat = Max(H1/Hsmlold,One)
         Hsmlold = Hsml
         Bbm = Max(1.5d0*Bmax,Bmax+Three)
         If (.not. Linear) Then
            Fmax = Dble(Itsaim)*Bmax/Dble(Nits)
            If (Fmax .gt. Bmax) Fmax = Max(Fmax,Bbm)
            If (Fmax .lt. Bmax) Fmax = Max(Bmax-Three,Bmax/1.5d0,Fmax)
         Endif
         Htot = Min(Amax,Bbm,Cmax)
         If (Linear) Fmax = Htot
         If (Nss .gt. 1 .and. Htot .ne. Cmax) Then
            If (Itru .eq. 0 .or. Amax .gt. Htpv) Htot = Max(Htot,Htpv)
            Htot = Min(Htot,Fmax,Cmax)
         Endif
         Htpv = Htot
         Phiaim = Htot/(Hsml*Hrat)
         Phiaim = Max(1.05d0*Phit,Phiaim)


C.... Quadratic And Linear Extrapolation To Find The Value Of Ep That
c.... Corresponds To Phiaim

         If (Nss .ne. 1) Then
            F1 = (Phi(2)-Phi(1))/(E(2)-E(1))
            F2 = (Phi(3)-Phi(2))/(E(3)-E(2))
            Ff2 = (F2-F1)/(E(3)-E(1))

            If (Ff2 .gt. 0) Then
               Qa = Ff2
               Qb = F1-(E(1)+E(2))*Ff2
               Qc = Phi(1)-E(1)*(F1-E(2)*Ff2)-Phiaim
               Qd = Qb**2-4.d0*Qa*Qc
               If (Ilin .eq. 1) Then
                  Ep = Min(Emin,(-Qb+Sqrt(Qd))/(Two*Qa))
                  If (Ep .eq. Emin) Phiaim = Phi(1)+
     +                 (Ep-E(1))*(F1+(Ep-E(2))*Ff2)
                  Phialt = Phi(2)+F2*(Ep-E(2))
               Else
                  Ep = Min(Emin,E(2)+(Phiaim-Phi(2))/F2)
                  If (Ep .eq. Emin) Phiaim = Phi(2)+F2*(Ep-E(2))
                  Phialt = Phi(1)+(Ep-E(1))*(F1+(Ep-E(2))*Ff2)
               Endif
            Else
               G1 = (E(2)-E(1))/(Phi(2)-Phi(1))
               G2 = (E(3)-E(2))/(Phi(3)-Phi(2))
               Gg2 = (G2-G1)/(Phi(3)-Phi(1))
               If (Ilin .eq. 1) Then
                  Ep = Min(Emin,
     +                 E(1)+(Phiaim-Phi(1))*(G1+(Phiaim-Phi(2))*Gg2))
                  If (Ep .eq. Emin) Then
                     Qa = Gg2
                     Qb = G1-(Phi(1)+Phi(2))*Gg2
                     Qc = E(1)-Phi(1)*(G1-Phi(2)*Gg2)-Ep
                     Qd = Qb**2-4.d0*Qa*Qc
                     Phiaim = (-Qb+Sqrt(Qd))/(Two*Qa)
                  Endif
                  Phialt = Phi(2)+F2*(Ep-E(2))
               Else
                  Ep = Min(Emin,E(2)+G2*(Phiaim-Phi(2)))
                  If (Ep .eq. Emin) Phiaim = Phi(2)+F2*(Ep-E(2))
                  Qa = Gg2
                  Qb = G1-(Phi(1)+Phi(2))*Gg2
                  Qc = E(1)-Phi(1)*(G1-Phi(2)*Gg2)-Ep
                  Qd = Qb**2-4.d0*Qa*Qc
                  Phialt = (-Qb+Sqrt(Qd))/(Two*Qa)
               Endif
            Endif

C.... Linear Extrapolation

         Else

            F2 = (Phi(3)-Phi(2))/(E(3)-E(2))
            Ep = Min(Emin,E(2)+(Phiaim-Phi(2))/F2)
         Endif



C.... Extrapolation May Not Be Working Very Well - If Not, Adjust
C.... Continuation Stepsize And Recalculate Phiaim.

         If ((Irest .eq. 1 .or. Istep .eq. 1) .and. Nss .ne. 1) Then
            If (Irest .eq. 1 .and. Istep .eq. 1) Then
               Ealt = Min(Epold*Hmult,Epold*Estep)
            Elseif (Irest .eq. 1) Then
               Ealt = Epold*Hmult
            Else
               Ealt = Epold*Estep
            Endif
            If (Ep .gt. Ealt) Then
               Ep = Min(Emin,Ealt)

               Philin = Phi(2)+F2*(Ep-E(2))
               If (Ff2 .gt. 0) Then
                  Phiqua = Phi(1)+(Ep-E(1))*(F1+(Ep-E(2))*Ff2)
               Else
                  Qa = Gg2
                  Qb = G1-(Phi(1)+Phi(2))*Gg2
                  Qc = E(1)-Phi(1)*(G1-Phi(2)*Gg2)-Ep
                  Qd = Qb**2-4.d0*Qa*Qc
                  Phiqua = (-Qb+Sqrt(Qd))/(Two*Qa)
               Endif
               If (Ilin .eq. 1) Then
                  Phiaim = Phiqua
                  Phialt = Philin
               Else
                  Phiaim = Philin
                  Phialt = Phiqua
               Endif
            Endif
         Endif

 70      Continue

C.... *****************************************************************
C.... End Of Logic For A Successful Continuation Step And
C.... Beginning Of Logic For An Unsuccessful Continuation Step.
C.... *****************************************************************

      Else


C.... If A Problem Is Not Solved Successfully Then Print Details And
C.... Select An Alternative Value For The Continuation Parameter.
C.... This Process Is Known As Backtracking. Istuk = 1 Indicates That
C.... We Are Backtracking.



         Istuk = 1

C.... If Iback = 1 Then The Final Problem Has Not Been Solved.
C.... In This Case We Stop.

         If (Iback .eq. 1) Then
C            If (Iprint .ge. 0) Then          Karline: toggled this off
               If (Epsp .ne. Zero) Then
        CALL Rprintd1('Final problem not solved, Epsilon = ',Eps)
        CALL Rprintd1('Try running problem again with Eps >', Eps)
        CALL Rprintd1('Machine prec. not Sufficient for Eps < ', Epsp)
               Else
        CALL Rprintd1('Final problem not solved, Epsilon = ',Eps)
        CALL Rprintd1('Try running problem again with Eps >', Eps)
               Endif
C            Endif
            if (nc .eq. Maxcon) Iflbvp=2
C KARLINE CHANGED THAT
            GOTO 100
         Endif

C.... If Iprec = 2, Then We Know That We Cannot Define A Mesh On
C.... Which The Current Eps Value Can Be Solved To The Requested
C.... Tolerances. We Alter Epsmin Accordingly.

         If (Iprec .eq. 2) Then
            Epsmin = One/(Max((Ep+E(3))/Two,0.9d0*Ep))
            eps_changed = .true.
         Endif

C.... Insert Details For Backtracking

         If (Iprint .ge. 0) then
         CALL Rprint('Failed step-Backtracking for larger value of Eps')
         end if
         Ifinal = 0
         Ep = (Ep+E(3))/Two

         If ((Ep-E(3)) .lt. 0.01d0*E(3)) Then
            Ep = E(3)
            Epsmin = One/Ep
            If (Iprint .ge. 0) then
        CALL Rprint('Continuation steps too small')
        CALL Rprintd1('Change Eps to ',Epsmin)
            end if
            Emin = Ep
            Iback = 1
C FRANCESCA deleted iflbvp=0
C KARLINE: ADDED THE GOTO... and iflag changed and epschanged...
            eps_changed = .true.
C FRANCESCA deleted iflbvp=-1
            GOTO 100
C KARLINE: CHANGES TILL HERE
         Endif

         If (Nss .eq. 1) Then
            Phiaim = Phi(2)+F2*(Ep-E(2))
         Else If (Nss .gt. 1) Then
            If (Iextrap .ne. 1) Then
               Philin = Phi(2)+F2*(Ep-E(2))
               If (Ff2 .gt. 0) Then
                  Phiqua = Phi(1)+(Ep-E(1))*(F1+(Ep-E(2))*Ff2)
               Else
                  Qa = Gg2
                  Qb = G1-(Phi(1)+Phi(2))*Gg2
                  Qc = E(1)-Phi(1)*(G1-Phi(2)*Gg2)-Ep
                  Qd = Qb**2-4.d0*Qa*Qc
                  Phiqua = (-Qb+Sqrt(Qd))/(Two*Qa)
               Endif

               If (Ilin .eq. 1) Then
                  Phiaim = Phiqua
                  Phialt = Philin
               Else
                  Phiaim = Philin
                  Phialt = Phiqua
               Endif
            Endif
         Endif

C.... Re-Insert Details From Last Problem Successfully Solved.

         If (Nss .gt. 0) Then
            Nmsh = Nnewold
            CALL acMeshdet(Xx, Wrk(Imsh), Nmsh)
            If (.not. Linear) CALL acSoldet(U, Wrk(Isol), Ncomp, Nmsh,
     +           Nudim, Ncomp)

Cf   save details for the next step
c            Nmshguess = Nmsh
c            Call Dcopy(Nmsh, Xx, 1, Xguess, 1)
c            Call Matcop(Nudim, Nudim, Ncomp, Nmsh, U, Uguess)
         Endif

      Endif

C.... *****************************************************************
C.... End Of Logic For An Unsuccessful Continuation Step.
C.... *****************************************************************

C.... Pass On Mesh And Solution Details To The Next Continuation
C.... Problem.

      If (Nss .gt. 0) Then
         Givmsh = .true.
         If (.not. Linear) Giveu = .true.
      Endif

C.... *****************************************************************
C.... The Program Loops Back To Line 50 For Every Continuation Step.
C.... *****************************************************************

      Goto 50

C KARLINE: ADDED THAT...
100   CONTINUE
      icount(1) = nfunc
      icount(2) = njac
      icount(3) = nbound
      icount(4) = njacbound
      icount(5) = nstep
      icount(6) = nc
      icount(7) = nss
      icount(8) = maxmesh
      Return

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      End

C------------------


      SUBROUTINE acMeshdet(X, X1, Nmsh)

      Implicit Double Precision (A-H,O-Z)

C.... This Subroutine Is Used For Saving And Re-Inserting Meshes

      Dimension X(*), X1(*)

      Do 10, J = 1,Nmsh
         X(J) = X1(J)
 10   Continue

      Return
      End

      SUBROUTINE acSoldet(U, U1, Ncomp, Nmsh, Idim1, Idim2)

      Implicit Double Precision (A-H,O-Z)

C.... This Subroutine Is Used For Saving And Re-Inserting Solutions.

      Dimension U(Idim1,*), U1(Idim2,*)

      Do 14, J = 1,Ncomp
         Do 12 K = 1,Nmsh
            U(J,K) = U1(J,K)
 12      Continue
 14   Continue

      Return
      End



* File Bvps.f

      SUBROUTINE acBvpsol(Ncomp, Nmsh, Nlbc, Aleft, Aright,
     *   Nfxpnt, Fixpnt, Ntol, Ltol, Tol, Nmax, Linear, Giveu,
     *   Givmsh, nmshguess, xguess,nugdim,uguess,
     *   Xx, Nudim, U, Def, Delu, Rhs, Fval,
     *   Topblk, Botblk, Ajac, Bhold, Chold, Ipvblk, Ipivlu,
     *   Uint, Ftmp, Tmprhs, Dftmp1, Dftmp2, Dgtm,
     *   Utrial, Rhstri, Xmerit, Xxold, Uold, Tmp, Dsq, Ratdc,
     *   Etest6, Etest8, Ermx, Ihcomp, Irefin, Def6, Def8, Dhold,
     *   Voldmsh, acfsub, acdfsub, acgsub, acdgsub, Iflbvp,
     *   amg, c1, wrkrhs,ckappa1,gamma1,sigma,ckappa,ckappa2,r4,
     *    precis,Eps, rpar,ipar)

      Implicit Double Precision (A-H,O-Z)
      DIMENSION RPAR(*), IPAR(*)
      Dimension  Fixpnt(*), Ltol(Ntol), Tol(Ntol)
      Dimension  Xx(*), U(Nudim, *),xguess(*),uguess(Nugdim,*)
      Dimension  Def(Ncomp,*), precis(3)
      Dimension  Delu(Ncomp, *), Rhs(*), Fval(Ncomp,*)
      Dimension  Topblk(Nlbc,*), Botblk(Ncomp-nlbc,*)
      Dimension  Ajac(Ncomp, 2*Ncomp, *)
      Dimension  Bhold(Ncomp, Ncomp, *), Chold(Ncomp, Ncomp, *)
      Dimension  Ipivlu(*), Ipvblk(*)
      Dimension  Uint(Ncomp), Ftmp(Ncomp), Dgtm(Ncomp), Tmprhs(*)
      Dimension  Dftmp1(Ncomp, Ncomp), Dftmp2(Ncomp, Ncomp)
      Dimension  Utrial(Ncomp,*), Rhstri(*), Xmerit(Ncomp, *)
      Dimension  Xxold(*), Uold(Ncomp,*), Tmp(Ncomp,16)
      Dimension  Dsq(Ncomp,Ncomp), Ratdc(*)
      Dimension  Etest6(*), Etest8(*), Ermx(*)
      Dimension  Ihcomp(*), Irefin(*)
      Dimension  Def6(Ncomp,*), Def8(Ncomp,*)
      Dimension  Dhold(3*Ncomp,3*Ncomp), Voldmsh(*)
      dimension  amg(*), c1(ncomp,*), wrkrhs(*), r4(*)
      Logical Linear, Giveu, Givmsh, dDouble, Errok

      External acfsub
      External acdfsub
      External acgsub
      External acdgsub
      common/Mcoldiag/nfunc, njac, nstep, nbound, njacbound, maxmesh
      Common/Mchprs/ flmin, Flmax, Epsmch
      Common /acFlags/ Ifinal,Iatt,Iback,Iprec,Iucond
      Common /acConvg/ Nits
      LOGICAL use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c
      Common/acAlgprs/Maxcon,Itsaim,Uval0

      logical stab_kappa, stab_gamma, stab_cond, stiff_cond, ill_cond
      logical stab_kappa1, ill_cond_newt, stab_sigma, comparekappa

      Intrinsic Max

      Logical Succes, Strctr
      Logical Ludone, Rhsgiv, Maxmsh

      Logical Mchset
      Save Mchset

      Parameter (Zero = 0.0d+0, One = 1.0d+0)
      Parameter (Third = 0.33d+0)
      Parameter (Quan6 = 0.1d+0 )

*  Blas: Dload
*  Double Precision D1mach

      Data Mchset/.true./
      Data Fxfct/10.0d+0/
      Data Maxmsh/.false./

*  The Parameter Inumb Is A Counter Used To Limit To Three The Number
*  Of Mesh Selections Performed For The Final Continuation Problem.

      Inumb = 0

      if (mchset) then
c Karline: use precis instead of d1mach

         flmin = precis(1)
         flmax = precis(2)
         epsmch = precis(3)
         mchset = .false.
      endif


*  The Routines Stcon1 And Stcon2 Calculate Integration Constants Stored
*  In Labeled Common Consts.

      CALL acStcon1
      CALL acStcon2



*  Set Up Arrays For The Error Tests.


      If (.not. Linear) Then
         CALL Dload(Ntol, One, Etest6, 1)
      Else
         Do 10 I = 1, Ntol
            Etest6(I) = One/Max(Quan6, Tol(I)**Third)
   10    Continue
      Endif

      Nmold = 1
      Strctr = .false.

*
*  If Givmsh Is .true., The Initial Number Of Mesh Points Must Be
*  Provided By The User In Nmsh, And The Mesh Points Must Be
*  Contained In The Array Xx (Of Dimension Nmsh).
*  Otherwise, Nmsh Is Set To Its Default Value, And A
*  Uniform Initial Mesh Is Created.

cf :jan 2011 moved outside bvpsol
c      If (.not. Giveu .and. .not. Givmsh) Then
c         Nmsh = Nminit
c         If (Nmsh .lt. Nfxpnt+2) Nmsh = Nfxpnt + 2
c         CALL acUnimsh(Nmsh, Aleft, Aright, Nfxpnt, Fixpnt, Xx)
c      Endif

c      If (Givmsh) Then
c        Nmsh = Nmshguess
c        Call Dcopy(Nmsh, xguess, 1, Xx, 1)
c      Endif

c      If (.not. Giveu) Then
c         Call Initu(Ncomp, Nmsh, Nudim, U)
c      else
c         Call Matcop(Nudim, Nudim, Ncomp, Nmsh, Uguess, U)
c      endif





      if (comp_c) then
*     initialize parameter for the conditioning estimation
      gamma1old  = flmax
      gamma1     = flmax
      ckappa1old = flmax
      ckappa1    = flmax
      ckappa     = flmax
      ckappaold  = flmax
      ckappa2    = flmax
      sigma      = 1

      endif

      stiff_cond = .false.
      ill_cond_newt =  .false.



***** Top Of Logic For 4th Order Solution ****

  400 Continue
      If (Iprint .ge. 0) then
        CALL Rprinti1('Number of points of the new mesh ',  Nmsh)
      end if
      maxmesh = max(maxmesh,nmsh)
      nstep   = nstep +1
*  Set The Def (Deferred Correction) Array To Zero.

      CALL Mtload(Ncomp, Nmsh-1, Zero, Ncomp, Def)
      Iorder = 4

*  The Routine Fneval Calls acfsub At The Mesh Xx And The
*  Solution U, And Saves The Values In The Array Fval.

      CALL acFneval(Ncomp,Nmsh,Xx,Nudim,U, Fval,acfsub, Eps,rpar,ipar)

*  Try To Compute A 4th Order Solution By Solving A System Of Nonlinear
*  Equations.


      If (Linear) Then
         Ludone = .false.

         CALL acLineq( Ncomp, Nmsh, Nlbc,
     *    Ludone, Xx, Nudim, U, Def,
     *    Delu, Rhs, Fval,
     *    Uint, Ftmp, Dftmp1, Dftmp2, Dgtm, Tmprhs,
     *    Ajac, Topblk, Botblk, Bhold, Chold, Ipvblk,
     *    acfsub, acdfsub, acgsub, acdgsub, Iflnwt, Eps,rpar,ipar)

*  CALL acFneval To Evaluate The Fval Array At The New Solution U.
*  (Such A Call Is Not Necessary For The Nonlinear Case Because
*  Fval Is Called Within Newteq For The New U.)

         CALL acFneval(Ncomp,Nmsh,Xx,Nudim,U,Fval,acfsub,Eps,rpar,ipar)

      Else
         Rhsgiv = .false.
         CALL acNewteq(Ncomp, Nmsh, Nlbc,
     *        Rhsgiv, Ntol, Ltol, Tol,
     *        Xx, Nudim, U, Def,
     *        Delu, Rhs, rnsq, Fval,
     *        Utrial, Rhstri,
     *        Uint, Ftmp, Dftmp1, Dftmp2, Dgtm, Tmprhs, Xmerit,
     *        Ajac, Topblk, Botblk, Bhold, Chold, Ipvblk,
     *        acfsub,acdfsub,acgsub,acdgsub,Itnwt,Iflnwt,Eps,rpar,ipar)
         If (Iatt .eq. -1) Nits = Max(1,Itnwt)
         If (Iprint .eq. 0) then
        CALL Rprinti1('Newton converged after iteration ', Itnwt)
         end if
      Endif



      If (Iflnwt .eq. 0) Then

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
            CALL acCONDESTIM(aleft,aright,nmsh,ncomp,N,xx,topblk,nlbc,
     *       ncomp, ajac, ncomp,2*ncomp,ninter,botblk,ncomp-nlbc,
     *       ipvblk,amg,c1,wrkrhs,ckappa1,gamma1,sigma,ckappa,ckappa2)


          if (iprint .ge. 0) then
        CALL Rprintd1('stiffness = ', sigma)
        CALL Rprintd1('gamma1    = ', gamma1)
        CALL Rprintd1('kappa1    = ', ckappa1)
        CALL Rprintd1('kappa     = ', ckappa)
        CALL Rprintd1('kappa2    = ',ckappa2)
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




        stiff_cond = (( (sigma .ge. 1d2 )))
        ill_cond   = ckappa2 .ge.  1d16 .and. ckappa2 .lt. flmax
        ill_cond_newt = ckappa2 .ge.  1d7 .and. ckappa2 .lt. flmax
     *        .and. gamma1 .ge.  1d5 .and. gamma1 .lt. flmax

          if (ill_cond .and. use_c) goto 2000

          if (iprint .ge. 1) then
        CALL Rprintd1('stab_sigma = ', stab_sigma)
        CALL Rprintd1('stab_kappa = ', stab_kappa)
        CALL Rprintd1('stab_kappa1 = ',stab_kappa1)
        CALL Rprintd1('stab_gamma = ', stab_gamma)
        CALL Rprintd1('stiff_cond = ', stiff_cond)
        CALL Rprintd1('ill_cond   = ', ill_cond)
          end if
        end if
c endif if (comp_c)



         CALL acDfexcl( Ncomp, Nmsh, Xx, Nudim, U, Def, Linear, Fval,
     *        Tmp, acfsub, acdfsub, Dsq, Ipivlu, Dhold,
     *        Ntol, Ltol, Tol, JC,Eps,rpar,ipar )

      Else

         Iflbvp = 1
         Return

      Endif

**** Logic For 6th Order ****

      If (Iprint .eq. 1) then
        CALL Rprint('Start 6th order')
      end if

*  Save The 4th Order Solution On This Mesh In Uold.

      Call Matcop(Nudim, Ncomp, Ncomp, Nmsh, U, Uold)
      Iorder = 6

      If (Linear) Then
        CALL acLineq( Ncomp, Nmsh, Nlbc,
     *    Ludone, Xx, Nudim, U, Def,
     *    Delu, Rhs, Fval,
     *    Uint, Ftmp, Dftmp1, Dftmp2, Dgtm, Tmprhs,
     *    Ajac, Topblk, Botblk, Chold, Bhold, Ipvblk,
     *    acfsub, acdfsub, acgsub, acdgsub, Iflnwt, Eps,rpar,ipar)

      Else
            Rhsgiv = .false.
        CALL acFixjac(Ncomp, Nmsh, Nlbc,
     *    Iorder, Ntol, Ltol, Tol,
     *    Xx, Nudim, U, Def, Def, Delu,
     *    Rhs, Fval, Utrial, Rhstri,
     *    Rnsq, Uint, Ftmp, Tmprhs,
     *    Ajac, Topblk, Botblk, Ipvblk,
     *    acfsub, acgsub, Iflnwt, Eps,rpar,ipar)

*  If The Fixed Jacobian Iterations Fail But Rnsq Is Small,
*  Try A Newton Procedure.  Set Rhsgiv To Indicate That
*  The Right-Hand Side And Fval Have Already Been Evaluated
*  At The Given U.

        If (Iflnwt .eq. -3 .and. Rnsq .lt. Fxfct*Epsmch) Rhsgiv = .true.
        If (Iflnwt .ne. 0) Then
           CALL acNewteq(Ncomp, Nmsh, Nlbc,
     *          Rhsgiv, Ntol, Ltol, Tol,
     *          Xx, Nudim, U, Def,
     *          Delu, Rhs, rnsq, Fval,
     *          Utrial, Rhstri,
     *          Uint, Ftmp, Dftmp1, Dftmp2, Dgtm, Tmprhs, Xmerit,
     *          Ajac, Topblk, Botblk, Bhold, Chold, Ipvblk,
     *          acfsub,acdfsub,acgsub,acdgsub,Iter,Iflnwt,Eps,rpar,ipar)
        Endif

      Endif


      If (Iflnwt .ne. 0) Then

         Iflbvp = 1
         Return

      Elseif (Ifinal .eq. 1) Then

         CALL acErrest (Ncomp, Nmsh, Ntol, Ltol, Tol,
     *        Nudim, U, Uold, Etest6, Errok)
         If (Errok) Then
            Iflbvp = 0
            Return
         Endif

      Endif

***** Logic For Trying To Calculate 8th Order Solution *****

      If (Iprint .eq. 1) then
         CALL Rprint('Start 8th order')
      end if

      Call Matcop(Nudim, Ncomp, Ncomp, Nmsh, U, Uold)


      call matcop(ncomp, ncomp, ncomp, nmsh-1, def, def6)
*  For Linear Problems, Calculate The Fval Array For The
*  New Solution U.

      If (Linear)
     +  CALL acFneval(Ncomp,Nmsh,Xx,Nudim,U,Fval,acfsub, Eps,rpar,ipar)

*  Calculate 8th Order Deferred Corrections (The Array Def8).

      CALL acDf8cal (Ncomp, Nmsh, Xx, Nudim, U, Fval, Def8, Linear,
     *             Tmp, acfsub, acdfsub, Dsq, Ipivlu, Dhold,
     *             Ntol, Ltol, Tol, JC, Eps,rpar,ipar)


*  For Linear Problems, The Def Array Is The Def8 Array.
*  For Nonlinear Problems, Add The Def8 Array To The
*  Already-Calculated Def Array.

      If (Linear) Then
         CALL Matcop(Ncomp, Ncomp, Ncomp, Nmsh-1, Def8, Def)
      Else
         CALL Maxpy(Ncomp, Nmsh-1, One, Def8, Ncomp, Def)
      Endif

      Iorder = 8

      If (Linear) Then
        CALL acLineq( Ncomp, Nmsh, Nlbc,
     *   Ludone, Xx, Nudim, U, Def,
     *   Delu, Rhs, Fval,
     *   Uint, Ftmp, Dftmp1, Dftmp2, Dgtm, Tmprhs,
     *   Ajac, Topblk, Botblk, Chold, Bhold, Ipvblk,
     *   acfsub, acdfsub, acgsub, acdgsub, Iflnwt, Eps,rpar,ipar)
      Else
              Rhsgiv = .false.
        CALL acFixjac(Ncomp, Nmsh, Nlbc,
     *    Iorder, Ntol, Ltol, Tol,
     *    Xx, Nudim, U, Def, Def8, Delu,
     *    Rhs, Fval, Utrial, Rhstri,
     *    Rnsq, Uint, Ftmp, Tmprhs,
     *    Ajac, Topblk, Botblk, Ipvblk,
     *    acfsub, acgsub, Iflnwt, Eps,rpar,ipar)

*  If The Fixed Jacobian Iterations Fail But Rnsq Is Small,
*  Try A Newton Procedure.  Set Rhsgiv To Indicate That
*  The Right-Hand Side And Fval Have Already Been Evaluated
*  At The Given U.

        If (Iflnwt .eq. -3 .and. Rnsq .lt. Fxfct*Epsmch) Rhsgiv = .true.
        if_fixjac = Iflnwt
        If (Iflnwt .ne. 0) Then
           CALL acNewteq(Ncomp, Nmsh, Nlbc,
     *          Rhsgiv, Ntol, Ltol, Tol,
     *          Xx, Nudim, U, Def,
     *          Delu, Rhs,rnsq, Fval,
     *          Utrial, Rhstri,
     *          Uint, Ftmp, Dftmp1, Dftmp2, Dgtm, Tmprhs, Xmerit,
     *          Ajac, Topblk, Botblk, Bhold, Chold, Ipvblk,
     *          acfsub,acdfsub,acgsub,acdgsub,Iter,Iflnwt,Eps,rpar,ipar)
        Endif
      Endif



      If (Iflnwt .eq. 0) Then

c         CALL acConv8( Ncomp, Nmsh, Ntol, Ltol, Tol,
c     *              Nfxpnt, Fixpnt, Linear, Nmax,
c     *              Xx, Nudim, U, Def8, Uold,
c     *              Ihcomp, Irefin, Ermx,
c     *              Etest8, Strctr, Ratdc, Voldmsh,
c     *              Double, Nmold, Xxold, Maxmsh, Succes)

           CALL acConv8( Ncomp, Nmsh, Ntol, Ltol, Tol,
     *              Nfxpnt, Fixpnt, Linear, Nmax,
     *              Xx, Nudim, U,nmshguess,xguess,nugdim,uguess,
     *              Def6, Def8, Uold,
     *              Ihcomp, Irefin, Ermx,
     *              Etest8, Strctr, Ratdc, Voldmsh,
     *              dDouble, Nmold, Xxold, Maxmsh, Succes,
     *       r4, amg,stab_cond,ckappa1,gamma1,ckappa,stiff_cond,
     *            ill_cond_newt,if_fixjac,rpar,ipar)




         If (Iprec .eq. 2) Then
            If (Iprint .ge. 0) then
        CALL Rprint('Mesh cannot be defined within the bounds imposed')
        CALL Rprint('by the machine precision')
            end if
            Iflbvp = 1
            Return
         Endif
      Else

         Iflbvp = 1
         Return

      Endif

      If (Maxmsh) Go To 900

      Iflbvp = 0

      If (Ifinal .eq. 1) Then
        If (Succes) Then
           Return
        Elseif (Iatt .ge. 1) Then
           If (Nmsh .lt. Nmax/2) Then
              Inumb = Inumb+1
           Else
              Iflbvp = 1
              Return
           Endif
        Endif
        If (Iback .ne. 1 .and. Inumb .eq. 3) Then
           Iflbvp = 1
           Return
        Endif
      Elseif (Iatt .eq. 0) Then
        Return
      Endif

      Iatt = Iatt+1

      Goto 400

  900 Continue

* Error Exit---Too Many Mesh Points.

      Iflbvp = 1

      Return

 2000  continue

* Error exit -- ill_cond
         iflbvp = 3
         if (linear .and. iprint .ge. 0) then
          CALL Rprint('The problem is ill-conditioned, ')
          CALL Rprint('Try with a less stringent tolerance')
         else
          CALL Rprint('The problem is ill-conditioned, ')
          CALL Rprint('Try with a less stringent tolerance ')
          CALL Rprint(' or with a different initial guess' )
         end if
         return
      End



      SUBROUTINE acConv8( Ncomp, Nmsh, Ntol, Ltol, Tol,
     *              Nfxpnt, Fixpnt, Linear, Nmax,
     *              Xx, Nudim, U,nmshguess,xguess,nugdim,uguess,
     *              Def6, Def8, Uold,
     *              Ihcomp, Irefin, Ermx,
     *              Etest8, Strctr, Ratdc, Voldmsh,
     *              dDouble, Nmold, Xxold, Maxmsh, Succes,
     *       r4, amg,stab_cond,ckappa1,gamma1,ckappa,stiff_cond,
     *               ill_cond_newt,if_fixjac, rpar,ipar)

      Implicit Double Precision (A-H,O-Z)
      Dimension Ltol(Ntol), Tol(Ntol)
      dimension rpar(*), ipar(*)
      Dimension Fixpnt(*), Etest8(Ntol)
      Dimension Xx(*), U(Nudim,*), xguess(*), Uguess(Nugdim,*)
      Dimension def6(ncomp,*),  Def8(Ncomp,*), Uold(Ncomp,*)
      Dimension Ihcomp(*), Irefin(*), amg(*), r4(*)
      Dimension Ermx(*), Xxold(*), Ratdc(*), Voldmsh(*)
      Logical Linear, Strctr, dDouble, Maxmsh, Succes
***bugfix 2jul01
      Save Nvold
      LOGICAL use_c, comp_c
      LOGICAL stiff_cond  ,stab_cond,  ill_cond_newt

      common/algprs/ nminit, iprint, idum, use_c, comp_c
      Common/acAlgprs/Maxcon,Itsaim,Uval0
      Common /acFlags/ Ifinal,Iatt,Iback,Iprec,Iucond
      Intrinsic Max

      Logical Errok

      Parameter (One = 1.0d+0, Fourth = 0.25d+0, Quan8 = 0.025d+0)

*  Blas: Dload

*  The Newton Iteration Converged For The 8th Order Solution.

      If (Iprint .eq. 1) then
        CALL Rprint('Conv8')
      end if

      Succes = .false.
      Maxmsh = .false.
      If (Ifinal .eq. 1) Then
         If (.not. Linear) Then
            Call Dload(Ntol, One, Etest8, 1)
         Else
            Do 10 I = 1, Ntol
               Etest8(I) = One/Max(Quan8, Tol(I)**Fourth)
 10         Continue
         Endif

*  Check Estimated Error.  For A Nonlinear Problem, All Components
*  Of Etest8 (The Ratios Used In Testing The Error) Are Set To One.
*  For A Linear Problem, The Components Of Etest8 Are In General
*  Larger Than One.  But If Strctr Is .true. And The Number Of Mesh
*  Points Decreased, We Set The Elements Of Etest8 To One (Which
*  Makes A Stricter Test).

         If (Linear .and. Strctr .and. Nmsh .lt. Nmold)
     *        Call Dload(Ntol,One, Etest8, 1)


         CALL acErrest (Ncomp, Nmsh, Ntol, Ltol, Tol,
     *        Nudim, U, Uold, Etest8, Errok)
         If (Errok) Then
            Succes = .true.
            Return
         Endif
      Endif

      If (Ifinal .eq. 1 .and. Iatt .ge. 1 ) Then

         CALL acDblmsh (Nmsh, Nmax, Xx, Nmold, Xxold, Maxmsh)
         If (.not. Maxmsh .and. Iprec .ne. 2) Then
***bugfix 2jul01
*           Call Matcop(Nudim, Ncomp, Ncomp, Nmsh, U, Uold)
            Call Matcop(Nudim, Ncomp, Ncomp, Nmold, U, Uold)
            Call acInterp(Ncomp, Nmsh, Xx, Nudim, U,ncomp,
     *           Nmold, Xxold, Uold)
         Endif
         Return
      Endif

*  Perform Selective Mesh Refinement Based On The 8th Order Deferred
*  Corrections.  The Value Of Ipow Indicates That The Error Estimate
*  Is Of Order 6.  Then, For A Nonlinear Problem, Interpolate The
*  Latest Solution Onto The New Mesh.

      Ipow = 6


*  NB: The array def8 will be overwritten by selmsh.


      if (use_c .and. stiff_cond .and. .not. ill_cond_newt) then
*  NB: The array def8 will be overwritten by selconderrmsh.
            CALL acselconderrmsh(Ncomp, Nmsh, Ntol, Ltol, Tol,
     *     Nfxpnt, Fixpnt, Ipow, Nmax, Nvold,
     *     Xx, Nudim, U, Def8, Irefin, Ihcomp,
     *     Nmold, Xxold, Ermx, dDouble, Maxmsh, Ratdc, Voldmsh,
     *     r4,amg,stab_cond,stiff_cond,linear)
      else

          CALL acSelmsh(Ncomp, Nmsh, Ntol, Ltol, Tol,
     *     Nfxpnt, Fixpnt, Ipow, Nmax, Nvold,
     *     Xx, Nudim, U, Def8, Irefin, Ihcomp,
     *     Nmold, Xxold, Ermx, dDouble, Maxmsh, Ratdc, Voldmsh)


      end if

*  The Array Def8 Will Be Overwritten By Selmsh.


      If (.not. Maxmsh .and. Iprec .ne. 2) Then
         If (Linear .and. (Ifinal .eq. 1 .or. Iatt .ne. 0)) Then
c              Call Initu(Ncomp, Nmsh, Nudim, U)
              call acinitu(ncomp, nmsh, xx, nudim, u,
     +            nugdim,nmshguess,xguess,uguess)
         Else
***bugfix 2jul01
*           Call Matcop(Nudim, Ncomp, Ncomp, Nmsh, U, Uold)
            Call Matcop(Nudim, Ncomp, Ncomp, Nmold, U, Uold)
            Call acInterp(Ncomp, Nmsh, Xx, Nudim, U,Ncomp,
     *           Nmold, Xxold, Uold)
         Endif
      Endif

      Return

      End

* File Defs.f


      SUBROUTINE acdfexcl (Ncomp, Nmsh, Xx, Nudim, U, Def6, Linear,
     *                 Fval, Tmp, acfsub, acdfsub, Df, Ip, Dhold,
     *                 Ntol, Ltol, Tol,JC,eps,rpar,ipar)

      Implicit Double Precision (A-H, O-Z)
      Integer Ncomp, Nmsh
      dimension rpar(*),ipar(*)
      Dimension Xx(Nmsh), U(Nudim,Nmsh), Fval(Ncomp, Nmsh)
      Dimension Def6(Ncomp,Nmsh-1), Tmp(Ncomp,*)
      Dimension Df(Ncomp,Ncomp), Ltol(Ntol), Tol(Ntol)
      Dimension Ip(2*Ncomp), Dhold(2*Ncomp,2*Ncomp)
      Dimension St1(200), St2(200), St3(200)
      External acfsub
      External acdfsub

      Parameter ( One = 1.0d+0, Two = 2.0d+0 )


      Common /acCons1/ A21,A22,A23,A24,A31,A32,A33,A34,C1,C2,
     +       C16,C26,C123,C223,C14,C24


      Logical Linear,use_c,comp_c

      common/algprs/ nminit, iprint, idum, use_c, comp_c
      Common/acAlgprs/Maxcon,Itsaim,Uval0
      Common /acFlags/ Ifinal,Iatt,Iback,Iprec,Iucond
      common/Mcoldiag/nfunc, njac, nstep, nbound, njacbound, maxmesh
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

         CALL acfsub (Ncomp,Xxc1,Tmp(1,1),Tmp(1,5),eps,rpar,ipar)
         nfunc=nfunc+1
         Call acfsub (Ncomp,Xxc2,Tmp(1,2),Tmp(1,6),eps,rpar,ipar)
         nfunc=nfunc+1
         Call acdfsub (Ncomp,Xxc1,Tmp(1,1),Df(1,1),eps,rpar,ipar)
         njac=njac+1
         Do 30 I = 1, Ncomp
            Tmp(I,5) = Tmp(I,5)-Tmp(I,3)
            Tmp(I,6) = Tmp(I,6)-Tmp(I,4)
            Do 25 J = 1,Ncomp
               Dfij = Hmsh*Df(I,J)
               Dhold(I,J) = -A22*Dfij
               Dhold(I,J+Ncomp) = -A23*Dfij
 25         Continue
 30      Continue

         Call acdfsub(Ncomp,Xxc2,Tmp(1,2),Df,eps,rpar,ipar)
         njac=njac+1
         Do 35 I = 1, Ncomp
            Do 32 J = 1, Ncomp
               Dfij = Hmsh*Df(I,J)
               Dhold(I+Ncomp,J) = -A32*Dfij
               Dhold(I+Ncomp,J+Ncomp) = -A33*Dfij
 32         Continue
 35      Continue

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
           Er = Tol(I)
           If (Abs(Tmp(Ii,7)) .gt. Er*Max(One,Abs(Tmp(Ii,3)))  .or.
     *         Abs(Tmp(Ii,8)) .gt. Er*Max(One,Abs(Tmp(Ii,4)))) Jc = 1
 50      Continue

         If (Jc .eq. 0) Goto 70

 60      Continue


         if (iprint.eq.1) then
          CALL Rprint('NO convergence of corrections')
         end if


         If (Iprec .eq. 0) Then
            If (Iprint .ge. 0) then
      CALL Rprint('** Warning - possibly approaching machine precision')
        CALL Rprintd1('beyond Epsilon  = ', Eps)
            end if
            Iprec = 1
         Endif

 70      Continue

         Do 80 Ic = 1, Ncomp
            Def6(Ic,Im) = (Hmsh/12.d+0)*(Fval(Ic,Im)+
     *              5.d+0*(Tmp(Ic,3)+Tmp(Ic,4))+Fval(Ic,Im+1))-
     *              U(Ic,Im+1)+U(Ic,Im)
 80      Continue
c      do 85 ic=1,ncomp
c      tmp(ic,5)=def6(ic,im)
c      tmp(ic,6)=def6(ic,im)
c 85         continue
c      call lusol(2*ncomp,2*ncomp,dhold,Ip,tmp(1,5),tmp(1,7))
c      do ic=1,ncomp
c         def8(ic,Im)=tmp(ic,7)
c      end do


 90   Continue

      Return


      End



      SUBROUTINE acdf8cal(Ncomp, Nmsh, Xx, Nudim, U, Fval, Def8,Linear,
     *                   Tmp, acfsub, acdfsub, Df, Ip, Dhold, Ntol,
     *                   Ltol, Tol,JC,eps,rpar,ipar)

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

      External acfsub
      External acdfsub


      Parameter ( One = 1.0d+0, Two = 2.0d+0 )

      Common/acCons2/ A21,A22,A23,A24,A25,A31,A32,A34,A35,A41,A42,
     +      A43,A44,A45,B1,B2,B3,C1,C2,C3,C16,C26,C36,C123,C223,
     +      C323,C14,C24,C34
      logical use_c,comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c
      Common/acAlgprs/Maxcon,Itsaim,Uval0
      Common /acFlags/ Ifinal,Iatt,Iback,Iprec,Iucond
      common/Mcoldiag/nfunc, njac, nstep, nbound, njacbound, maxmesh
      Logical Linear

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

            Call acfsub (Ncomp,Xxc1,Tmp(1,1),Tmp(1,7),eps,rpar,ipar)
            nfunc=nfunc+1
            Call acfsub (Ncomp,Xxc2,Tmp(1,2),Tmp(1,8),eps,rpar,ipar)
            nfunc=nfunc+1
            Call acfsub (Ncomp,Xxc3,Tmp(1,3),Tmp(1,9),eps,rpar,ipar)
            nfunc=nfunc+1
            Call acdfsub(Ncomp,Xxc1,Tmp(1,1),Df,eps,rpar,ipar)
            njac=njac+1
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

            Call acdfsub(Ncomp,Xxc2,Tmp(1,2),Df,eps,rpar,ipar)
            njac=njac+1
            Do 40 I = 1, Ncomp
                Do 35 J = 1, Ncomp
                   Dfij = Hmsh*Df(I,J)
                   Dhold(I+Ncomp,J) = -A32*Dfij
                   Dhold(I+Ncomp,J+Ncomp) = 0.d+0
                   Dhold(I+Ncomp,J+2*Ncomp) = -A34*Dfij
 35             Continue
 40         Continue

            Call acdfsub(Ncomp,Xxc3,Tmp(1,3),Df,eps,rpar,ipar)
            njac=njac+1
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

            Jc = 0
            If (Linear) Goto 90

            Do 70 I = 1, Ntol
               Ii = Ltol(I)
               Er = Tol(I)/Hmsh
               Er = Tol(I)
               If (Abs(Tmp(Ii,10)) .gt. Er*Max(One,Abs(Tmp(Ii,4))) .or.
     *           Abs(Tmp(Ii,11)) .gt. Er*Max(One,Abs(Tmp(Ii,5))) .or.
     *           Abs(Tmp(Ii,12)) .gt. Er*Max(One,Abs(Tmp(Ii,6)))) Jc = 1
 70                                 Continue

            If (Jc .eq. 0) Goto 90

 80         Continue

         if (iprint.eq.1) then
          CALL Rprint('NO convergence of 8th order defcors')
         end if



         If (Iprec .eq. 0) Then
            If (Iprint .ge. 0)  then
      CALL Rprint('** Warning -possibly approaching machine precision')
      CALL Rprintd1('beyond Epsilon  = ',Eps)
             end if
            Iprec = 1
         Endif


 90   Continue

         Do 100 Ic = 1, Ncomp
           Def8(Ic,Im) = Hmsh*(B1*(Fval(Ic,Im)+Fval(Ic,Im+1))+
     *                   B2*(Tmp(Ic,4)+Tmp(Ic,6))+B3*Tmp(Ic,5))-
     *                   U(Ic,Im+1)+U(Ic,Im)
 100     Continue

 110     Continue

      Return

      End



      SUBROUTINE acStcon1
      Implicit Double Precision(A-H,O-Z)
      Common/acCons1/A21,A22,A23,A24,A31,A32,A33,A34,C1,C2,
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

      SUBROUTINE acStcon2
      Implicit Double Precision(A-H,O-Z)
      Common/acCons2/ A21,A22,A23,A24,A25,A31,A32,A34,A35,A41,A42,
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

* File Eqs.f




      SUBROUTINE acfixjac(ncomp, nmsh, nlbc,
     *    iorder, ntol, ltol, tol,
     *    xx, nudim, u, defcor, defnew, delu,
     *    rhs, fval, utrial, rhstri,
     *    rnsq, uint, ftmp, tmprhs,
     *    ajac, topblk, botblk, ipivot,
     *    acfsub, acgsub, iflag,eps,rpar,ipar)

* Fixed Jacobian iterations.

      implicit double precision (a-h,o-z)
      dimension rpar(*),ipar(*)
      dimension  ltol(ntol), tol(ntol)
      dimension  xx(nmsh), u(nudim,nmsh), defcor(ncomp,nmsh-1)
      dimension  defnew(ncomp,nmsh-1), delu(ncomp,nmsh)
      dimension  rhs(ncomp*nmsh), fval(ncomp,nmsh)
      dimension  utrial(ncomp,nmsh), rhstri(ncomp*nmsh)
      dimension  uint(ncomp), ftmp(ncomp), tmprhs(ncomp*nmsh)
      dimension  ajac(ncomp, 2*ncomp, nmsh-1)
      dimension  topblk(nlbc,ncomp), botblk(ncomp-nlbc,ncomp)
      dimension  ipivot(ncomp*nmsh)
      logical    better

      external   acfsub
      external   acgsub

      logical  use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c
      Common/acAlgprs/Maxcon,Itsaim,Uval0
      Common/Mchprs/ Flmin,Flmax, Epsmch
      Common /acFlags/ Ifinal,Iatt,Iback,Iprec,Iucond

      intrinsic  abs
      intrinsic  max

*  blas: dcopy, dssq

      parameter ( one    = 1.0d+0 )
      parameter ( xlarge = 1.0d+6, huge = 1.0d+20, lmtfrz = 8)
      parameter ( rngrow = 16.0d+0, rfact = 100.0d+0 )
      parameter ( tolfct = 0.1d+0 )


*  The iteration scheme uses a fixed Jacobian matrix to solve for
*  correction vectors, once there has been convergence of the Newton
*  iterations on this mesh.   It is assumed that the LU
*  factors of the Jacobian have been left unaltered since
*  their calculation.

      if (iprint .eq. 1) then
        CALL Rprint('Fixed Jacobian iterations')
      end if
      ninter = nmsh - 1
      rnold = flmax
      isize=nmsh*ncomp

*  Evaluate the right-hand side rhstri at the initial solution u by
*  adding the new deferred corrections to the already-calculated
*  rhs vector.

      call dcopy(nlbc, rhs, 1, rhstri, 1)
      ind = nlbc
      do 10 im = 1, ninter
      do 10 ic = 1, ncomp
         ind = ind + 1
         rhstri(ind) = rhs(ind) + defnew(ic, im)
 10       continue
      ind = ninter*nmsh + nlbc + 1
      call dcopy(ncomp-nlbc, rhs, 1, rhstri, 1)

      call dssq  ( nmsh*ncomp, rhstri, 1, scale, sumsq )
      rnsq = (scale**2)*sumsq

      iter = 0

*  If the initial right-hand side is too large, do not even attempt to
*  solve the nonlinear equations.

      if (rnsq.gt.huge .or.
     *      (iorder.eq. 8 .and. rnsq.gt.xlarge)) then
         if (iprint .eq. 1) then
        CALL Rprintd1('Large residual, rnsq =', rnsq)
         end if
         iflag = -2
         return
      end if
      call dcopy(ncomp*nmsh, rhstri, 1, rhs, 1)

*  Statement 100 is the top of the iteration loop.

 100   continue

*  If rnsq is sufficiently small, terminate immediately.

      if (iprint .eq. 1) then
        CALL Rprintid('iter, rnsq', iter, rnsq)
      end if
      if (rnsq .le. epsmch) then
         iflag = 0
         return
      endif

      iter = iter + 1

*  Solve for the step delu by solving a system involving the fixed
*  Jacobian (whose LU factors are saved).  Copy the rhs array into
*  tmprhs, which is then overwritten by blkslv.

      call dcopy(ncomp*nmsh, rhs, 1, tmprhs, 1)
      call dcopy(ncomp*nmsh,tmprhs,1,delu,1)
      job = 0
      call crslve(Topblk, Nlbc, Ncomp, Ajac, Ncomp, 2*Ncomp,
     + Ninter, Botblk, Ncomp-Nlbc, Ipivot, Delu, job)

*  Compute the trial point utrial by adding delu to u.

      call matcop( nudim, ncomp, ncomp, nmsh, u, utrial )
      call maxpy( ncomp, nmsh, one, delu, ncomp, utrial )

*  compute the right-hand side vector rhstri and its squared
*  two-norm at the trial point.

      rnold = rnsq
      CALL acfneval(ncomp,nmsh,xx,ncomp,utrial,fval,acfsub,
     *   eps,rpar,ipar)
      CALL acrhscal (ncomp, nmsh, nlbc, xx, ncomp, utrial, defcor,
     *   acfsub, acgsub, rhstri, rnsq, fval, ftmp, uint,eps,rpar,ipar)

*  If rnsq strictly decreased, update the solution vector u
*  and the right-hand side rhs.

      better = .false.
      if (rnsq .lt. rnold) then
         better = .true.
         call matcop( ncomp, nudim, ncomp, nmsh, utrial, u )
         call dcopy( ncomp*nmsh, rhstri, 1, rhs, 1 )
      endif

*  Stop the fixed Jacobian iterations if there have been too
*  many iterations, or if rnsq has not decreased by a factor
*  of at least rngrow.

      if (iter .ge. lmtfrz .or. rnsq .gt. (rnold/rngrow) ) then


         if (better) then

*  Setting iflag to -3 signals that, although the fixed Jacobian
*  iterations did not succeed, the current point was an improvement
*  on the previous one.  Hence, if we switch to a Newton procedure,
*  the right-hand side does not need to be recalculated.

            iflag = -3
         else
            iflag = -2
         endif
         if (iprint .eq. 1) then
        CALL Rprinti1('Failure of fixed Jacobian, iflag =',iflag)
         end if
         return
      endif


*  Test For Convergence Using The Ratio Abs((Change In U)/Max(U,1)).


      do 150 im = 1, nmsh
      do 150 it = 1, ntol
         itol = ltol(it)
c Francesca using rhs instead of delu
c         er = abs(rhs( itol+(im-1)*ncomp))/max(abs(u(itol,im)), one)
         Er = Abs(Delu(Itol,Im))/Max(Abs(U(Itol,Im)), One)
         if (er .gt. tolfct*tol(it).and. er .gt. epsmch
     +              ) go to 100
 150      continue

*  To exit from the loop here, the convergence tests have
*  been passed.

      if (iprint .ge. 0) then
        CALL Rprintid('Fixed Jacobian convergence',iter, rnsq)
      end if
      iflag = 0
      return
      end

************************************************************************





      SUBROUTINE aclineq( ncomp, nmsh, nlbc,
     *    ludone, xx, nudim, u, defcor,
     *    delu, rhs, fval,
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,
     *    ajac, topblk, botblk, bhold, chold, ipivot,
     *    acfsub, acdfsub, acgsub, acdgsub, iflag,eps,rpar,ipar)
***********************************************************************
*     call by: bvpsol
*     calls to: lnrhs, jaccal, dcopy, colrow, dload, dload, clrslve
*           maxpy
**********************************************************************
      implicit double precision (a-h,o-z)
      DIMENSION RPAR(*), IPAR(*)
      dimension  xx(nmsh), u(nudim, nmsh), defcor(ncomp, nmsh-1)
      dimension  delu(ncomp, nmsh), rhs(ncomp*nmsh)
      dimension  fval(ncomp,nmsh), uint(ncomp), ftmp(ncomp)
      dimension  topblk(nlbc,*), botblk(ncomp-nlbc,*)
      dimension  dftmp1(ncomp, ncomp), dftmp2(ncomp, ncomp)
      dimension  dgtm(ncomp), tmprhs(ncomp*nmsh)
      dimension  ajac(ncomp, 2*ncomp, nmsh-1),
     *               bhold(ncomp, ncomp, nmsh-1),
     *               chold(ncomp, ncomp, nmsh-1)
      dimension  ipivot(ncomp*nmsh)
      integer    job
      logical    ludone
      external   acfsub
      external   acdfsub
      external   acgsub
      external   acdgsub
      LOGICAL use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c
      Common/acAlgprs/Maxcon,Itsaim,Uval0

*  blas: dcopy, dload

      parameter  ( one = 1.0d+0, zero = 0.0d+0 )

*  The routine lineq calculates the Newton step for a linear
*  problem.  The Newton step is exact unless the Jacobian
*  matrix is singular.

      isize=nmsh*ncomp
      ninter = nmsh - 1

      if (.not. ludone) then

*  Compute the right-hand side vector rhs.
         iflag = 0
         CALL aclnrhs (ncomp, nmsh, nlbc, xx, nudim, u,
     *       acfsub,acgsub,rhs,rnsq,fval,ftmp, uint,eps,rpar,ipar)

*  If the Jacobian for this mesh has not previously been
*  calulated and factorized successfully, CALL acjaccal.
*  The block-structured Jacobian matrix is stored in three
*  matrices (topblk, ajac, and botblk).
*  The matrices bhold and chold are also calculated in jaccal,
*  and are saved for later use in outer routines.

         CALL acjaccal (ncomp, nmsh, nlbc,
     *      xx, nudim, u, fval, dgtm, dftmp1, dftmp2, uint,
     *      ajac, topblk, botblk, bhold, chold,
     *      acdfsub, acdgsub,eps,rpar,ipar)

*  Call blkdcm to calculate the LU factors of the Jacobian.
*  The factors are overwritten on the matrices topblk, ajac and botblk.
*  Interchanges are represented in the integer array ipivot.

      call dcopy(ncomp*nmsh,rhs,1,tmprhs,1)
      call dcopy(ncomp*nmsh,tmprhs,1,delu,1)


      job = 0
      call colrow(isize, topblk,nlbc, ncomp, ajac, ncomp, 2*ncomp,
     +   ninter, botblk, ncomp-nlbc, ipivot, delu, iflag, job)

         ludone = .true.

*  Copy the rhs into the temporary vector tmprhs, which will be
*  overwritten by blkslv.


      else

*  The right-hand side is the deferred correction array,
*  padded with zeros at the boundary conditions.
        iflag = 0
         call dload(nlbc, zero, tmprhs(1), 1)
         do 100 im = 1, ninter
            loc = (im-1)*ncomp + nlbc + 1
            call dcopy(ncomp, defcor(1,im), 1, tmprhs(loc), 1)
 100            continue
         nrhs = ninter*ncomp + nlbc + 1
         call dload(ncomp-nlbc, zero, tmprhs(nrhs), 1)
         call dcopy(ncomp*nmsh,tmprhs,1,delu,1)
      job = 0

      CALL crslve(topblk,nlbc,ncomp,ajac,ncomp,2*ncomp,ninter,
     +   botblk,ncomp-nlbc,ipivot,delu,job)

      endif

*  Since the problem is linear, the Newton step  is exact.  The
*  new u array is obtained by adding delu to u.

      CALL maxpy ( ncomp, nmsh, one, delu, nudim, u )

c      iflag = 0
      return

      end




******************************************************************
      SUBROUTINE acnewteq(ncomp, nmsh, nlbc,
     *    rhsgiv, ntol, ltol, tol,
     *    xx, nudim, u, defcor,
     *    delu, rhs, rnsq, fval,
     *    utrial, rhstri,
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,
     *    ajac, topblk, botblk, bhold, chold, ipivot,
     *    acfsub, acdfsub, acgsub, acdgsub, iter, iflag,
     *    eps,rpar,ipar)
******************************************************************
*     call by: bvpsol
*     calls to:
******************************************************************
      implicit double precision (a-h,o-z)
      DIMENSION RPAR(*), IPAR(*)
      dimension  ltol(*), tol(*), xx(*)
      dimension  fval(ncomp,*)
      dimension  u(nudim, *), delu(ncomp, *), utrial(ncomp,nmsh)
      dimension  rhs(ncomp*nmsh),  defcor(ncomp,*)
      dimension  ajac(ncomp, 2*ncomp, *)
      dimension  topblk(nlbc,*), botblk(ncomp-nlbc,*)
      dimension  ftmp(*), uint(*), dgtm(ncomp)
      dimension  dftmp1(ncomp, ncomp), dftmp2(ncomp, ncomp)
      dimension  bhold(ncomp, ncomp, *), chold(ncomp, ncomp, *)
      dimension  ipivot(*)
      dimension  rhstri(ncomp*nmsh)
      dimension  tmprhs(ncomp*nmsh), xmerit(ncomp, nmsh)
      logical rhsgiv

      LOGICAL use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c
      Common/acAlgprs/Maxcon,Itsaim,Uval0
      external   acfsub
      external   acdfsub
      external   acgsub
      external   acdgsub

      parameter  ( zero   = 0.0d+0, one    = 1.0d+0 )
      parameter ( two = 2.0d+0, half = 0.5d+0, fourth = 0.25d+0 )
      parameter ( tenth = 0.1d+0, ten = 10.0d+0, hund = 100.0d+0 )


      common/mchprs/ flmin, flmax, epsmch

      logical gtpdeb, imprvd, braktd, crampd, extrap, vset, wset
      save  gtpdeb, mfsrch, epsaf, epsag, eta, rmu, tolabs, alfmax
      save  tolrel, toltny

      intrinsic  abs
      intrinsic  max

*  blas: dcopy

      logical frscal
      save frscal
      parameter (cnvfct = 0.1d+0 )
      data  frscal/.true./
c frscal has been updated in bvpsol
c
      data  alfsml/1.0d-4/,  alfmax/1.1d+0/
      data  imerit/1/, lmtnwt/39/
      data  shrfct/100.0d+0/, stpfct/2.0d+0/
      data  gtpdeb/.false./, mfsrch/5/
      data  eta/.999999d+0/, rmu/1.0d-6/


*  The routine newteq_l performs Newton iterations with a line
*  search, to solve the nonlinear equations.


*  Set up constants if this is the first call to newteq_l.

      if (frscal) then
         frscal = .false.
         epsaf = epsmch
         epsag = epsmch
         tolabs = epsmch
         tolrel = epsmch
         toltny = epsmch
      endif
      ninter = nmsh - 1

      if (iprint .eq. 1) then
        CALL Rprint('Start Newton iterations')
      end if
*  A Newton method with line search and watchdog safeguarding
*  is used to try to solve the nonlinear equations.

*  Initialize iter (the counter of Newton iterations) and alfold
*  (the step taken at the previous iteration).

      iter = -1
      alfold = one
      alfa = zero


      if (.not. rhsgiv) then

*  If necessary, evaluate the right-hand side at the initial u.
          CALL acfneval(ncomp, nmsh, xx, nudim, u, fval,acfsub,
     *        eps,rpar,ipar)

          CALL acrhscal (ncomp, nmsh, nlbc, xx, nudim, u, defcor,
     *      acfsub, acgsub, rhs, rnsq, fval, ftmp, uint,eps,rpar,ipar)
      endif



      rnsqtr = rnsq

*  At any given Newton iteration, rnprev is the value of rnsq at
*  the immediately preceding Newton iteration.

      rnprev = flmax
      rnbest = flmax
      if ( iprint .ge. 0) then
        CALL Rprint(' iter  , alfa , merit,   rnsq')
      end if
*  Initialize counter of watchdog iterations.

      itwtch = 0

*  Statement 100 is the top of the Newton iteration loop.

 100   continue

      iter = iter + 1

      if (iprint .eq. 1) then
        CALL Rprinti1('Newton iteration ', iter)
      end if
*  If there have been too many Newton iterations, terminate.

      if (iter .ge. lmtnwt) then
         if (iprint .ge. 0) then
          CALL rprint('Too many Newton iterations')
         end if
         iflag = -2
         return
      endif

*  The vector rhs is the right-hand side at the current iterate,
*  and rnsq is its squared two-norm.
*  Perform watchdog tests, using the unscaled merit function (rnsq)
*  as the watchdog function.  The routine wtchdg updates rnbest
*  and itwtch.  If iflwat is not zero, this sequence of Newton
*  iterations is terminated.

      iflwat = 0

      CALL acwtchdg ( iter, rnsq, rnbest, rnprev, itwtch,
     *                alfold, iflwat )

      if (iflwat .ne. 0) then
         if (iprint .ge. 0) then
        CALL Rprinti1('Watchdog tests fail, iter = ', iter)
         end if
         iflag = -3
         return
      endif

*  Watchdog tests are passed.  Proceed with the Newton iteration.
*  Call jaccal to evaluate the block-structured Jacobian matrix,
*  which is stored in three matrices (topblk, ajac, and botblk).
*  The matrices bhold and chold are saved for use in later
*  calculations in the outer routine.


*  If rnsq is sufficiently small, terminate immediately.
*  Note that the stored Jacobian does not correspond exactly
*  to the final point.

      if (rnsq .le. epsmch) then
         if (iprint .ge. 0)  then
        CALL Rprintid('Convergence, iter , and  rnsq =', iter, rnsq)
         end if
         iflag = 0
         return
      endif

      CALL acjaccal (ncomp, nmsh, nlbc,
     *    xx, nudim, u, fval, dgtm, dftmp1, dftmp2, uint,
     *    ajac, topblk, botblk, bhold, chold,
     *    acdfsub, acdgsub,eps,rpar,ipar)

*  blkdcm is called to calculate the LU factors of the Jacobian,
*  which are overwritten on topblk, ajac and botblk.
*  Interchanges are represented in the integer array ipivot.

*   Solve for the Newton step delu.  Copy the rhs array into tmprhs,
*   which is then overwritten by blkslv.

      call dcopy(ncomp*nmsh, rhs, 1, tmprhs, 1)
      call dcopy(ncomp*nmsh,tmprhs,1,delu,1)
      job = 0
      call colrow(Nmsh*ncomp,topblk,nlbc,ncomp,ajac,ncomp,2*ncomp,
     +   ninter,botblk,ncomp-nlbc,ipivot,delu,iflag,job)

      if (iprint .ge. 0 .and. iflag.ne.0) then
        CALL Rprinti1('Singular Jacobian, iter= ',iter)
      end if
      if (iflag .ne. 0) then
         iflag = -1
         return
      end if
*     the jacobian is singular


*  If imerit = 1, the line search is based on the scaled merit function,
*  the squared two-norm of the solution xmerit of the linear system
*      (Jacobian)*xmerit = rhs,
*  where (Jacobian) is the Jacobian at the current Newton iterate.
*  Thus the initial value of the scaled merit function is simply
*  the squared two-norm of the Newton step delu itself.

      if (imerit.eq.1) then
         call mssq( ncomp, nmsh, delu, xmscal, xmsq )
         fmtry = (xmscal**2)*xmsq
      else
*  The unscaled merit function is simply the squared two-norm of rhs.
         fmtry = rnsq
      end if

c  fa and oldg represent the merit function and its gradient
c  at the initial point of the line search.

      fa = fmtry
      oldg = -two*fa
      alfa = zero
      if (iprint .eq. 1) then
        CALL Rprintd3('alfa, merit, rnsq', alfa, fmtry, rnsq)
      end if
*  On the first Newton iteration, the initial trial step is unity.
*  On subsequent iterations, the initial step is not allowed to
*  be more than the factor stpfct larger than the final step at
*  the immediately preceding iteration.


      alfa = one
      if(stpfct*alfold .lt. one) alfa = stpfct*alfold

      if (alfa .lt. alfsml) alfa = alfsml

      fmold = fa
      inform = -1

*  Statement 150 is the top of the inner line search iteration.
*  The line search routine getptq has been altered so that it
*  terminates with an indication of success as soon as a
*  strictly lower value of the merit function is found.  Note that
*  this is a much less strict requirement than the usual sufficient
*  decrease conditions.


 150   continue
      iwr = 6

      CALL acgetptq (gtpdeb, mfsrch, iwr, alfmax, alfsml, alfuzz,
     *      epsaf, epsag,
     *      eta, fmtry, fmold, oldg, rmu, tolabs, tolrel, toltny,
     *      imprvd, inform, nfsrch, alfa, alfbst, fbest,
     *      braktd, crampd, extrap, vset, wset, nsamea, nsameb,
     *      alin, blin, fa, factor, fv, fw, xtry, xv, xw)



*  inform = 1, 2 or 3 indicates success in finding an acceptable point.
*  inform = 4 means alfmax is too small (this should never happen here,
*  since alfmax is set always to 1.1).
*  inform = 5 means that a decrease was not achieved for any step
*  greater than alfsml.
*  inform = 6 means a better point could not be found (the minimum
*  probably lies too close to alfa=0).
*  inform = 7 means that the gradient at alfa=0 (oldg) is positive
*  (this cannot happen here, since oldg=-two*fa, and fa is a non-negative
*  number)

      if (iprint.eq.1) then
        CALL Rprintid('Inform, alfa after getptq',inform, alfa)
      end if
      if (inform .eq. 5) then
          iflag = -5
          return
      elseif (inform .eq. 4 .or. inform .eq. 7) then
         iflag = -4
         return
      elseif (inform .eq. 0) then

*  inform = 0 means that a new function value should be obtained
*  with the step alfa.
*  We may override alfa from getptq by requiring that the step is not
*  allowed to decrease by more than a factor of shrfct during
*  a line search iteration.
*
         if (alfa .lt. alfold/shrfct) alfa = alfold/shrfct
         alfold = alfa

*  Define the next iterate utrial = u + alfa*delu.
*  Call fneval and rhscal to evaluate the right-hand side
*  rhstri at utrial.
*  The vector rhstri is stored separately, and rhs is overwritten
*  only when an improved point is found.

         call matcop ( nudim, ncomp, ncomp, nmsh, u, utrial)
         call maxpy ( ncomp, nmsh, alfa, delu, ncomp, utrial )
         CALL acfneval(ncomp, nmsh, xx, ncomp, utrial, fval,acfsub,
     *        eps,rpar,ipar)

         CALL acrhscal (ncomp, nmsh, nlbc, xx, ncomp, utrial, defcor,
     *      acfsub,acgsub,rhstri,rnsqtr,fval,ftmp, uint,eps,rpar,ipar)

         fmold = fmtry
         if (imerit .eq. 1) then

*  Solve a linear system to obtain the 2-d array xmerit whose squared
*  norm is the scaled merit function.   The LU factors of the Jacobian
*  have already been calculated by blkdcm.
*  Copy rhstri into tmprhs, which is overwritten by blkslv.

      call dcopy(ncomp*nmsh, rhstri, 1, tmprhs, 1)
      call dcopy(ncomp*nmsh,tmprhs,1,xmerit,1)
      job = 0
      call crslve(topblk, nlbc, ncomp, ajac, ncomp, 2*ncomp, ninter,
     +   botblk, ncomp-nlbc, ipivot, xmerit,job)
      call mssq( ncomp, nmsh, xmerit, xscale, xsolsq )
      fmtry = (xscale**2)*xsolsq
         else

*  The unscaled merit function is the squared two-norm of the right-hand
*  side.
            fmtry = rnsqtr
         end if
         if (iprint .eq. 1) then
        CALL Rprintd3('alfa, merit, rnsq', alfa, fmtry, rnsqtr)
         end if
         go to 150
      endif

*  To reach here, inform must be 1, 2, 3, or 6, and the line search
*  has found a strictly lower value of the merit function.
*  Store the new Newton iterate in u, and the corresponding rhs
*  vector in rhs.

      rnprev = rnsq
      rnsq = rnsqtr


      call matcop (ncomp, nudim, ncomp, nmsh, utrial, u)
      call dcopy(ncomp*nmsh, rhstri, 1, rhs, 1)
      if (iprint .ge. 0) then
        CALL Rprintid3(' ',iter, alfa, fmtry, rnsq)
      end if
*  Now test for convergence using the ratio of the Newton step
*  for each component with max(1, abs(current solution estimate)).
*  If the test fails for any element of u, branch back to the
*  top of the Newton iteration.

      do 160 im = 1, nmsh
      do 160 it = 1, ntol
         icmp = ltol(it)
         er = abs(delu(icmp,im))/max(abs(u(icmp,im)), one)
         if (er .gt. cnvfct*tol(it)) go to 100
 160      continue


      if (iprint .ge. 0) then
      CALL Rprintid('Convergence, iter = ,  rnsq =', iter+1, rnsq)
      end if
      iflag = 0

*  To fall through the above loop, the termination test for a
*  sufficiently small delu is satisfied.
*  Note that the stored Jacobian and its factorization do not
*  correspond to the final solution.


      return

      end



      SUBROUTINE acwtchdg ( iter, wmerit, wmbest, wmprev,
     *      itwtch, alfold, iflag )

*  Logic for watchdog tests.

      implicit double precision (a-h,o-z)
      parameter ( itonew = 5, itwtmx = 8, grfct = 100.0d+0 )
      parameter ( half = 0.5d+0 )

*  Perform watchdog tests in two forms:
*  (1) to determine whether a sufficient decrease in the
*  watchdog merit function has occurred within the most recent
*  sequence of itwtmx iterations;
*  (2) to determine whether the watchdog merit function has increased
*  too much in a single iteration after itonew Newton iterations
*  have been performed.  This allows the merit function to increase
*  wildly only during the first itonew iterations.

*  wmbest is the smallest watchdog merit function achieved in this
*  sequence of Newton iterations.
*  wmprev is the watchdog merit function from the immediately
*  preceding Newton iteration.

*  itwtch counts the number of iterations without an improvement
*  in the unscaled merit function.

      iflag = 0
      if (wmerit .le. wmbest) then

*  The current watchdog merit function is the best.

         wmbest = wmerit
         itwtch = 0
         return
      endif

*  The current merit function is not the best.

      itwtch = itwtch + 1

*  Do not apply watchdog tests if (1) the previous step alfold
*  exceeds 1/2, or (2) the watchdog merit function decreased in
*  the immediately preceding iteration and itwtch does not
*  exceed twice its maximum.

      if (alfold .ge. half) return
      if (wmerit .le. wmprev .and. itwtch .le. 2*itwtmx) return


*  If more than itwtmx iterations have occurred without
*  an overall improvement in the watchdog merit function,
*  signal for termination.

      if (itwtch .ge. itwtmx) then
         iflag = -1

*  If a too-large increase in the watchdog merit function
*  compared to the best value occurred, and iter .ge. itonew,
*  signal for termination.

      elseif (iter .ge. itonew .and.
     *          wmerit .gt. grfct*wmbest) then
          iflag = -1
      endif
      return
      end


* File Evals.f

      SUBROUTINE acFneval(Ncomp, Nmsh, Xx, Nudim, U, Fval, acfsub,
     +                                          Eps,rpar,ipar)
      Implicit Double Precision (A-H,O-Z)
      Dimension rpar(*),ipar(*)
      Dimension Xx(Nmsh), U(Nudim,Nmsh), Fval(Ncomp,Nmsh)
      common/Mcoldiag/nfunc, njac, nstep, nbound, njacbound, maxmesh
      External acfsub

*  Fneval Evaluates The Function Values (From acfsub) For
*  A Given Mesh Xx And Array U, And Stores The Values
*  In The Array Fval.

      Call acfsub (Ncomp, Xx(1), U(1,1), Fval(1,1), Eps,rpar,ipar)
      nfunc=nfunc+1
      Do 50 Im = 1, Nmsh-1
       Call acfsub(Ncomp,Xx(Im+1),U(1,Im+1),Fval(1,Im+1),Eps,rpar,ipar)
         nfunc=nfunc+1
   50 Continue
      Return
      End


      SUBROUTINE acjaccal (ncomp, nmsh, nlbc,
     *   xx, nudim, u, fval,
     *   dgtm, dftm1, dftm2, uint,
     *   ajac, topblk, botblk, bhold, chold,
     *   acdfsub, acdgsub,eps,rpar,ipar)

      implicit double precision (a-h,o-z)
      DIMENSION RPAR(*), IPAR(*)
      dimension xx(nmsh), u(nudim,nmsh), fval(ncomp,nmsh)
      dimension dgtm(ncomp)
      dimension dftm1(ncomp, ncomp), dftm2(ncomp, ncomp),
     *             uint(ncomp)
      dimension ajac(ncomp, 2*ncomp, nmsh-1)
      dimension topblk(nlbc, ncomp), botblk(ncomp-nlbc,ncomp)
      dimension bhold(ncomp, ncomp, nmsh-1),
     *             chold(ncomp, ncomp, nmsh-1)
      common/Mcoldiag/nfunc, njac, nstep, nbound, njacbound, maxmesh
      external  acdfsub
      external  acdgsub

      logical pdebug, use_c, comp_c


*  blas: dcopy, ddot

      parameter ( zero = 0.0d+0, half = 0.5d+0, eighth = 0.125d+0 )
      parameter ( four = 4.0d+0, six = 6.0d+0 )
      parameter ( one = 1.0d+0, three = 3.0d+0, twelve = 12.0d+0 )


      ninter = nmsh - 1

      do 110 i = 1, nlbc
         call acdgsub (i, ncomp, u(1,1), dgtm,eps,rpar,ipar)
         njacbound=njacbound+1
         call dcopy(ncomp, dgtm(1), 1, topblk(i,1), nlbc)
 110      continue

      call acdfsub (ncomp, xx(1), u(1,1), dftm1(1,1),eps,rpar,ipar)
      njac=njac+1
*  on entry to jaccal, the array fval contains the function values
*  at (xx(im), u(ic,im)), ic=1,...,ncomp and im = 1,...,nmsh,
*  calculated by a preceding call of rhscal with the same xx and u
*  arrays.

      do 200 im = 1, ninter

         hmsh = xx(im+1) - xx(im)

         do 120 ic = 1, ncomp
            uint(ic) = half*(u(ic,im) + u(ic,im+1))
     *         - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))
 120            continue
         xhalf = half*(xx(im+1) + xx(im))
         call acdfsub (ncomp, xhalf, uint, dftm2(1,1),eps,rpar,ipar)
         njac=njac+1
         do 140 ic = 1, ncomp
            do 130 jc = 1, ncomp
               dsq = ddot(ncomp, dftm2(ic,1), ncomp,
     *                       dftm1(1,jc), 1)
               ajac(ic,jc,im) = -hmsh*(dftm1(ic,jc)/six
     *             + dftm2(ic,jc)/three + hmsh*dsq/twelve)
 130                  continue
            ajac(ic,ic,im) = ajac(ic,ic,im) - one
 140            continue

         call acdfsub(ncomp,xx(im+1),u(1,im+1),dftm1(1,1),eps,rpar,ipar)
         njac=njac+1
         do 170 ic = 1, ncomp
            do 160 jc = 1, ncomp
               dsq = ddot(ncomp, dftm2(ic,1), ncomp,
     *                         dftm1(1,jc), 1)
               ajac(ic,jc+ncomp,im) = -hmsh*(dftm1(ic,jc)/six
     *               + dftm2(ic,jc)/three - hmsh*dsq/twelve)
 160                  continue
            call dcopy(ncomp, ajac(ic,ncomp+1,im), ncomp,
     *                   chold(ic,1,im), ncomp)
            call dcopy(ncomp, dftm1(ic,1), ncomp,
     *                   bhold(ic,1,im), ncomp)
            ajac(ic,ic+ncomp,im) = ajac(ic,ic+ncomp,im) + one
            chold(ic,ic,im) = ajac(ic,ic+ncomp,im)
 170            continue


 200             continue
      do 220 i = nlbc+1, ncomp
         call acdgsub (i, ncomp, u(1, nmsh), dgtm,eps,rpar,ipar)
         njacbound=njacbound+1
         call dcopy(ncomp, dgtm(1), 1, botblk(i-nlbc,1), ncomp-nlbc)
 220      continue

      return

      end




      SUBROUTINE aclnrhs (ncomp, nmsh, nlbc,
     *   xx, nudim, u, acfsub, acgsub,
     *   rhs, rnsq, fval, ftmp, uint,eps,rpar,ipar)

       implicit double precision(a-h,o-z)

*  This SUBROUTINE acis designed to calculate the right-hand
*  side for linear problems.
      DIMENSION RPAR(*), IPAR(*)
      dimension xx(*), u(nudim,*)
      dimension rhs(*), fval(ncomp,*), ftmp(*), uint(*)
      external acfsub
      external acgsub

      parameter ( zero = 0.0d+0, half = 0.5d+0, eighth = 0.125d+0 )
      parameter ( one = 1.0d+0, four = 4.0d+0, six = 6.0d+0 )
      common/Mcoldiag/nfunc, njac, nstep, nbound, njacbound, maxmesh
      common/mchprs/ flmin, flmax, epsmch
      intrinsic abs

*  blas: dssq

      ninter = nmsh - 1
      rnsq = zero

*  first, process the left-hand boundary conditions.

      do 20 i = 1, nlbc
         call acgsub (i, ncomp, u(1,1), wg,eps,rpar,ipar)
         nbound=nbound+1
         rhs(i) = -wg
   20 continue

*  Next, process the interior mesh points.  The fval array
*  contains the function values from acfsub at xx and u.

      do 50 im = 1, ninter
         hmsh = xx(im+1) - xx(im)
         do 30 ic = 1, ncomp
            uint(ic) = half*(u(ic,im) + u(ic,im+1))
     *         - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))
   30    continue
         xhalf = half*(xx(im) + xx(im+1))
         call acfsub (ncomp, xhalf, uint, ftmp,eps,rpar,ipar)
         nfunc=nfunc+1
         loc = (im-1)*ncomp + nlbc
         do 40 ic = 1, ncomp
            rhs(loc+ic) = -u(ic,im+1) + u(ic,im)
     *        + hmsh*
     *            (fval(ic,im) + fval(ic,im+1) + four*ftmp(ic))/six
   40    continue
   50 continue

      nrhs = ninter*ncomp
      do 60 ii = nlbc+1, ncomp
         call acgsub (ii, ncomp, u(1,nmsh), wg,eps,rpar,ipar)
         nbound=nbound+1
         rhs(nrhs+ii) = -wg
   60 continue

      call dssq  ( nmsh*ncomp, rhs, 1, scale, sumsq )
      rnsq = (scale**2)*sumsq

      return
      end



      SUBROUTINE acrhscal (ncomp, nmsh, nlbc,
     *   xx, nudim, u, defcor,
     *   acfsub, acgsub,
     *   rhs, rnsq, fval, ftmp, uint,eps,rpar,ipar)

       implicit double precision(a-h,o-z)

*  This SUBROUTINE acconstructs the (ncomp*nmsh)-dimensional
*  vector rhs, which is the right-hand side of the Newton equations.
*  The ncomp by nmsh array fval is assumed to have been calculated
*  elsewhere by routine fneval.
      DIMENSION RPAR(*), IPAR(*)
      dimension  xx(nmsh), u(nudim,nmsh), defcor(ncomp,nmsh-1)
      dimension  rhs(ncomp*nmsh), fval(ncomp,nmsh)
      dimension  ftmp(ncomp), uint(ncomp)
      common/Mcoldiag/nfunc, njac, nstep, nbound, njacbound, maxmesh
      external   acfsub
      external   acgsub

      logical pdebug, use_c, comp_c

      common/mchprs/flmin, flmax, epsmch
      intrinsic abs

*  blas: dssq

      parameter ( zero = 0.0d+0, half = 0.5d+0, eighth = 0.125d+0 )
      parameter ( one = 1.0d+0, four = 4.0d+0, six = 6.0d+0 )

*  ninter is the number of intervals in the mesh (one less than the
*  number of mesh points)

      ninter = nmsh - 1
      rnsq = zero

*  First, process the left-hand boundary conditions.

      do 20 i = 1, nlbc
         call acgsub (i, ncomp, u(1,1), wg,eps,rpar,ipar)
         nbound=nbound+1
         rhs(i) = -wg
 20       continue

*  Next, process the interior mesh points.  The fval array
*  contains the function values from acfsub at xx and u.

      do 50 im = 1, ninter
         hmsh = xx(im+1) - xx(im)
         do 30 ic = 1, ncomp
            uint(ic) = half*(u(ic,im) + u(ic,im+1))
     *         - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))
 30             continue
         xhalf = half*(xx(im) + xx(im+1))
         call acfsub (ncomp, xhalf, uint, ftmp,eps,rpar,ipar)
         nfunc=nfunc+1
         loc = (im-1)*ncomp + nlbc
         do 40 ic = 1, ncomp
            rhs(loc+ic) = -u(ic,im+1) + u(ic,im) + defcor(ic,im)
     *        + hmsh*
     *            (fval(ic,im) + fval(ic,im+1) + four*ftmp(ic))/six

 40             continue
 50              continue

      nrhs = ninter*ncomp
      do 60 ii = nlbc+1, ncomp
         call acgsub (ii, ncomp, u(1,nmsh), wg,eps,rpar,ipar)
         nbound=nbound+1
         rhs(nrhs+ii) = -wg
 60       continue


      call dssq  ( nmsh*ncomp, rhs, 1, sscale, sumsq )
      rnsq = (sscale**2)*sumsq
*f

      return

      end




      SUBROUTINE acDblmsh (Nmsh, Nmax, Xx, Nmold, Xxold, Maxmsh)
      Implicit Double Precision (A-H,O-Z)
      Dimension Xx(*), Xxold(*)
      Logical Maxmsh
      LOGICAL use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c
      Common/acAlgprs/Maxcon,Itsaim,Uval0
      Common /acFlags/ Ifinal,Iatt,Iback,Iprec,Iucond
*  Blas: Dcopy

      Parameter (Half = 0.5d+0)

*  This Routine Is Used To Double The Mesh, I.e., Produce A Mesh
*  With Twice As Many Intervals In Which Each New Interval Is
*  Half The Corresponding Old Interval.

*  On Entry To Dblmsh, The Integer Nmsh And The Array Xx
*  Specify A Set Of Mesh Points Xx(1),..., Xx(Nmsh) (Assumed
*  To Be In Ascending Order).

*  If The Number Of Mesh Points In The Doubled Mesh Would
*  Exceed The Maximum Allowed Number Nmax, The Flag Maxmsh Is
*  Set To True, And We Exit Without Changing Any Other Parameters.

*  Otherwise, Nmold Is Set To The Old Number Of Mesh Points,
*  Xxold Is Set To The Old Mesh, Nmsh Is The New Number Of Mesh
*  Points, And Xx Contains The New Mesh Points.

      Nmold = Nmsh
      Iprec = Min(Iprec,1)
      Call Dcopy(Nmold, Xx, 1, Xxold, 1)

      Ninnew = 2*(Nmsh-1)
      Nmnew = Ninnew + 1
      If(Nmnew .ge. Nmax) Then
         If (Iprint .ge. 0)  then
      CALL Rprinti1(' Dblmsh.  Maximum Mesh Exceeded, Nmnew  = ',Nmnew)
         end if
         Maxmsh = .true.
         Return
      Endif
      Maxmsh = .false.

*  Loop Backwards Through The Old Mesh Points To Create The New Ones.

      Xx(Nmnew) = Xx(Nmsh)
      Do 100 I = Ninnew, 4, -2
         Id2 = I/2
         Xx(I) = Half*(Xx(I+1) + Xx(Id2))
         Xx(I-1) = Xx(Id2)
         If (Xx(I) .eq. Xx(I+1) .or. Xx(I) .eq. Xx(I-1)) Then
            Iprec = 2
            Return
         Endif
  100 Continue

*  Calculate The New Xx(2). Xx(1) Remains Unchanged.

      Xx(2) = Half*(Xx(3) + Xx(1))
      If (Xx(2) .eq. Xx(3) .or. Xx(2) .eq. Xx(1)) Then
         Iprec = 2
         Return
      Endif
      Nmsh = Nmnew
      Return
      End


      SUBROUTINE acSelmsh(Ncomp, Nmsh, Ntol, Ltol, Tol,
     *     Nfxpnt, Fixpnt, Ipow, Nmax, Nvold,
     *     Xx, Nudim, U, Ermeas, Irefin, Ihcomp,
     *     Nmold, Xxold, Ermx, dDouble, Maxmsh, Phiold, Voldmsh)

      Implicit Double Precision (A-H,O-Z)

      Dimension  Ltol(Ntol), Tol(Ntol), Fixpnt(*)
      Dimension  Xx(*), U(Nudim, *), Ermeas(Ncomp,*)
      Dimension  Irefin(Nmsh-1), Ihcomp(Nmsh-1)
      Dimension  Xxold(*), Ermx(*), Phiold(*), Voldmsh(*)
      Dimension  Npr(2), Pmax(2), Hord(2)
      Logical    dDouble, Maxmsh

      Intrinsic Abs
      Intrinsic Max
      Intrinsic Int

*  Blas: Dcopy

      Common /acFlags/ Ifinal,Iatt,Iback,Iprec,Iucond
      Common /acMshvar/ Hsml,Npr,Pmax,Hord

      Parameter  ( Zero = 0.0d+0, One = 1.0d+0, Onep1 = 1.1d+0 )
      Parameter  ( Erdcid = 5.0d+0 )
      Parameter  ( Phitst = 0.1d+0 )

*  The Routine Selmsh Performs Selective Mesh Refinement, Depending
*  On The Error Measure Ermeas.

      Maxmsh = .false.

      Frcpow = One/Ipow
      dDouble = .false.
      Nmold = Nmsh
      Ninter = Nmsh - 1
      Iprec = Min(Iprec,1)

*  Copy The Current Mesh Into The Xxold Array.

      If (Iatt .eq. -1) Nvold = Nmsh
      If (Iatt .eq. 0) Call Dcopy(Nvold, Xxold, 1, Voldmsh, 1)

      Call Dcopy(Nmold, Xx, 1, Xxold, 1)
      Ithres = 0
      Thres = One

*  On Input, The Array Ermeas Represents Some Error Measure Defined
*  Over The Components And Mesh Intervals (Not Mesh Points).
*  It Is Normalized In The Following Loop With Respect To The
*  Tolerance Array And The Current Solution.
*  The Value Errmax Gives The Maximum Normalized Error.

      Errmax = Zero
      Do 120 Im = 1, Ninter
         Ermx(Im) = Zero
         Do 110 It = 1, Ntol
            Jcomp = Ltol(It)
            Umin = Min(Abs(U(Jcomp,Im)),Abs(U(Jcomp,Im+1)))
            Denom = Tol(It)*Max(One, Umin)
            Ems = Ermeas(Jcomp,Im)
            Ermeas(Jcomp,Im) = Abs(Ems)/Denom
            Err = Ermeas(Jcomp, Im)
            If (Err .ge. Ermx(Im)) Then
                Ermx(Im) = Err
                Ihcomp(Im) = Jcomp
            Endif
  110    Continue
         Errmax = Max(Ermx(Im), Errmax)
  120 Continue

  200 Continue

*  For Each Interval Im,  The Integer Irefin(Im) Is Calculated
*  Based On An Equidistribution Procedure Involving The
*  Threshold Value Thres.  If Irefin(Im) > 1, We Add
*  Points To Interval Im.  If Irefin(Im) = 1, We Check Whether
*  Point Im Can Be Removed From The Mesh.

*  Nmest Is A Lower Bound On The Number Of Points In The New Mesh.
*  We Do Not Know In Advance How Many Points Will Be Removed,
*  So Nmest Is Computed By Assuming All Eligible Points Are Removed.

      Philrg = 0.d0
      Esum = 0.d0
      Nmest = Nmsh
      Mshchng = 0
      Do 220 Im = 1, Ninter
         Errim = Ermx(Im)**Frcpow
         If (Ermx(Im) .ge. Thres) Then
            Irefin(Im) = Int(Errim) + 1
            Nmest = Nmest + Irefin(Im) - 1
            Mshchng = 1
         Else
            Irefin(Im) = 1
            Nmest = Nmest - 1
         Endif
         Him = Xxold(Im+1)-Xxold(Im)
         Phiim = Errim/Him
         Esum = Esum+Errim
         If (Iatt .eq. -1) Then
            Phiold(Im) = Phiim
         Elseif (Iatt .eq. 0 .and. Phiim .gt. Philrg) Then
           Philrg = Phiim
           Hordlrg = Errim
           Imreg = Im
         Endif
  220 Continue
C


      If (Iatt .eq. 0) Then
        Hord(2) = Hordlrg
        Pmax(2) = Philrg
        Xloc1 = Xxold(Imreg)
        Xloc2 = Xxold(Imreg+1)
        Ichkpt = Nvold/2
        If (Xloc1 .lt. Voldmsh(Ichkpt)) Ichkpt = 1
        Xb = Voldmsh(Ichkpt)
        Do 225 J = Ichkpt, Nvold-1
          Xa = Xb
          Xb = Voldmsh(J+1)
          If (Xloc1 .ge. Xa .and. Xloc1 .lt. Xb) Then
            If (Xloc2-Xb .lt. (Xloc2-Xloc1)/2.d0) Then
              Pmax(1) = Phiold(J)
              Hord(1) = (Xb-Xa)*Pmax(1)
            Else
              Pmax(1) = Phiold(J+1)
              Hord(1) = (Voldmsh(J+2)-Xb)*Pmax(1)
            Endif
            Goto 228
          Endif
 225    Continue
 228    Hsml=(Xloc2-Xloc1)/Dble(Irefin(Imreg))
        Npr(2)=Int(Esum)
      Endif

      If (Iatt .eq. -1) Npr(1) = Int(Esum)
      If (Ifinal .eq. 1 .and. Mshchng .eq. 0) Iatt = 0
      If (Nmest .gt. Nmax) Go To 360

*  It Appears That We Can Perform The Desired Selective Mesh
*  Refinement.

*  Now Begin Running Through The Mesh, Adding And Possibly Deleting
*  Points As Indicated By The Irefin Array.

*  The Integer New Is A Count Of The Number Of Intervals In
*  The Tentative Mesh Being Generated By The Refinement Strategy.

      New = 1

*  The First Interval Is Treated As A Special Case, Since Xx(1)
*  Always Remains In The Mesh, And Cannot Be A Fixed Point.

      Rlen = Xxold(2) - Xx(1)
      Slen = Rlen
      If (Irefin(1) .gt. 1) Then
         Dx = Rlen/Irefin(1)
         Do 230 J = 2, Irefin(1)
            New = New + 1
            Xx(New) = Xx(1) + (J-1)*Dx

            If (Xx(New) .eq. Xx(New-1)) Then
               Iprec = 2
               Return
            Endif
  230    Continue
      Endif

*  The Fixed Points Specified By The Fixpnt Array Cannot Be
*  Removed From The Mesh.  The Value Fxnext Indicates The 'Next'
*  Fixed Point.  When No Further Fixed Points Remain To Be Processed
*  (Or If Nfxpnt = 0), Fxnext Is Set To A Value Strictly Larger Than
*  The Last Mesh Point, So That No Mesh Point Can Equal Fxnext.
*  This Way We Need To Compare Only The Values Of Xxold(I)
*  And Fxnext.

      Ifxcnt = 1
      If (Nfxpnt .eq. 0) Then
         Fxnext = Onep1*Abs(Xxold(Nmsh))
      Else
         Fxnext = Fixpnt(Ifxcnt)
      Endif

*  Jtkout Is A Counter Of The Number Of Consecutive Points That
*  Have Been Removed From The Mesh.

      Jtkout = 0
      Do 330 Im = 2, Ninter
         Rlold = Rlen
         Rlen = Xxold(Im+1) - Xxold(Im)

*  If Xxold(Im) Is The Next Fixed Point, It Cannot Be Removed
*  And So We Don'T Test Its Error Estimates.

         If(Xxold(Im) .eq. Fxnext)  Then

            Ifxcnt = Ifxcnt + 1
            If(Ifxcnt .gt. Nfxpnt) Then
               Fxnext = Onep1*Abs(Xxold(Nmsh))
            Else
               Fxnext = Fixpnt(Ifxcnt)
            Endif

         Elseif (Irefin(Im) .eq. 1) Then

*  If Xxold(Im) Is Not A Fixed Point And Irefin(Im) = 1, We May Wish
*  To Remove Point Im From The Mesh.

*  If We Are Considering Removing Points And Jtkout = 0, This
*  Is The First Point In A Possible Consecutive Set To Be Removed,
*  And We Initialize Phihat, Which Represents A Maximum Of
*  Certain Estimates.
*  If Jtkout Is Not Zero, Previous Points Contiguous To This
*  Point Have Been Removed, And Phihat Retains Its Previous Value.

            Slen = Slen + Rlen

            If (Jtkout .eq. 0) Then
                Ind1 = Ihcomp(Im-1)
                Phihat = Ermeas(Ind1,Im-1)/(Rlold**Ipow)
            Endif
            Phihat = Max(Phihat,
     *                 Ermeas(Ihcomp(Im),Im)/(Rlen**Ipow))
            Val1 = Phihat*(Slen**Ipow)
            If (Val1 .le. Phitst
     *             .and. Jtkout .lt. 4) Then

*  Increment The Counter Of Removed Points.
*  'Remove' The Mesh Point Xxold(Im) By Not Including It.

               Jtkout = Jtkout+1
               Go To 330
            Endif
*        End Of Logic For Irefin(Im) = 1.
         Endif

         Jtkout = 0
         New = New + 1
         Xx(New) = Xxold(Im)

         If (Xx(New) .eq. Xx(New-1)) Then
            Iprec = 2
            Return
         Endif
         If (Irefin(Im) .gt. 1) Then
            Dx = Rlen/Irefin(Im)
            Do 300 J = 2, Irefin(Im)
               New = New + 1

               if (New .ge. Nmax) goto 360
               Xx(New) = Xxold(Im) + (J-1)*Dx

               If (Xx(New) .eq. Xx(New-1)) Then
                  Iprec = 2
                  Return
               Endif
 300        Continue
         Endif
         Slen = Rlen

*  If The New Mesh Contains Too Many Points, Branch Out Of The
*  Loop To Try Alternative Strategies.

         If (New .gt. Nmax-1) Goto 360

  330 Continue

*  To End Up Here, We Have Processed The Entire Interval,
*  And Have Neither Exceeded The Specified Maximum Nor
*  Exceeded Three Times The Number Of Intervals In The Old
*  Mesh.  The Last Mesh Point Remains Unchanged.

      New = New + 1
      Xx(New) = Xxold(Nmsh)
      If (Xx(New) .eq. Xx(New-1)) Then
         Iprec = 2
         Return
      Endif
      Nmsh = New
      Maxmsh = .false.
      Return

  360 Continue

*  To Reach Here, The Number Of Mesh Points Created At Some Stage
*  Of The Refinement Process Was Larger Than The Maximum Permitted
*  Value Nmax.

*  Check Whether The Mesh Can Safely Be Doubled.

      If ((2*Nmsh-1) .lt. Nmax) Then

*  Double The Mesh.
         CALL acDblmsh (Nmsh, Nmax, Xx, Nmold, Xxold, Maxmsh)
         dDouble = .true.

*  If The Number Of Intervals Is Too Large And The Mesh Cannot Be
*  Doubled, Increase The Threshold Thres By A Factor Of Erdcid And
*  Try The Selective Refinement Again.
*  If This Happens Three Times Without Success Or If Thres Exceeds
*  Or Is Equal To Errmax, Stop.  (In This Case, We Know Already
*  That Doubling The Mesh Produces Too Many Points.)

      Elseif (Thres .lt. Errmax .and. Ithres .lt. 3) Then
         Ithres = Ithres + 1
         Thres = Erdcid*Thres
         If(Thres .gt. Errmax) Thres = Errmax
         Call Dcopy(Nmsh, Xxold, 1, Xx, 1)
         Go To 200
      Else
         Nmsh = 2*Nmsh - 1
         Maxmsh = .true.
      Endif
      Return

      End


      SUBROUTINE acUnimsh(Nmsh, Aleft, Aright, Nfxpnt, Fixpnt, Xx)
      Implicit Double Precision (A-H,O-Z)
      Integer  Nmsh, Nfxpnt
      Dimension Fixpnt(*), Xx(Nmsh)

      Intrinsic Max

*  Given A Left Endpoint Aleft, A Right Endpoint Aright,
*  A Set Of Nfxpnt Fixed Points Fixpnt(I), I = 1,...,Nfxpnt,
*  (Where Fixpnt(I) Is Different From Aleft And Aright For All I),
*  And An Initial Target Number Nmsh Of Mesh Points,
*  The SUBROUTINE acUnimsh Generates A Piecewise Uniform Mesh
*  Beginning At Aleft, Ending At Aright, And With Equally
*  Spaced Points Between Aleft And Fixpnt(1), Then Between
*  Fixpnt(1) And Fixpnt(2), ..., And Finally Between
*  Fixpnt(Nfxpnt) And Aright.  The Final Number Of Intervals
*  Is The Maximum Of Nfxpnt+2 And The Initial Value Of Nmsh.

*  In The Simplest Case When Nfxpnt = 0, Unimsh Generates A
*  Uniform Mesh With Nmsh Intervals In The Closed Interval
*  (Aleft, Aright).

*  On Exit, The Integer Nmsh Contains The Number Of Mesh Points
*  (Which Is The Maximum Of The Initial Nmsh And Nfxpnt).
*  The Array Xx (Of Dimension Nmsh) Contains The Mesh Points.

      If (Nfxpnt .eq. 0) Then

*  If There Are No Interior Fixed Points, The Spacing Is Uniform
*  Throughout The Interval.  Calculate The Spacing Dx
*  And Set Up The Xx Array.

        Ninter = Nmsh - 1

         Dx = (Aright - Aleft)/Ninter
         Do 10 I = 1, Ninter
            Xx(I) = Aleft + (I-1)*Dx
   10    Continue
         Xx(Nmsh) = Aright
         Return
      Endif

*  We Know That There Is At Least One Fixed Point Strictly Between
*  The Endpoints.

      If (Nmsh .lt. Nfxpnt+2)  Nmsh = Nfxpnt + 2
      Ninter = Nmsh - 1
      Xx(1) = Aleft
      Ileft = 1
      Xleft = Aleft
      Totint = Aright - Aleft
      Ndif = Ninter - Nfxpnt
      Do 50 J = 1, Nfxpnt + 1

*  Deal In Turn With The Subintervals Defined By The Interval
*  Boundaries And The Fixed  Points.

         If (J .lt. Nfxpnt+1) Then

*  The J-Th Fixed Point Is Xright.  Calculate Where It Should
*  Fall In The Mesh.

            Xright = Fixpnt(J)
            Nmin = Int(Ninter*(Xright-Aleft)/Totint + 1.5d+0)
            If (Nmin .gt. Ndif+J) Nmin = Ndif + J
            Iright = Max(Ileft+1, Nmin)
         Else
            Xright = Aright
            Iright = Nmsh
         Endif

*  Npt Is The Number Of Equally Spaced Points That Should
*  Lie Strictly Between The (J-1)-Th And J-Th Fixed Points.

         Xx(Iright) = Xright
         Npt = Iright - Ileft - 1
         Dx = (Xright - Xleft)/(Npt + 1)
         Do 30 I = 1, Npt
            Xx(Ileft+I) = Xleft + I*Dx
   30    Continue
         Ileft = Iright
         Xleft = Xright
   50 Continue

      Return
C
      End


* File Rest.f

      SUBROUTINE acErrest (Ncomp, Nmsh, Ntol, Ltol, Tol,
     *   Nudim, U, Uold, Etest, Errok)

      Implicit Double Precision (A-H,O-Z)
      Dimension Ltol(Ntol), Tol(Ntol), U(Nudim,Nmsh),
     *              Uold(Ncomp,Nmsh), Etest(Ntol)
      Logical Errok

      Intrinsic Abs
      Intrinsic Max

      Parameter ( One = 1.0d+0 )


*  Given Current And Previous Solutions U And Uold On The Same
*  Mesh, Errest Calculates An Error Measure For Each
*  Component For Which A Tolerance Is Specified.
*  The Error Measure Is The Usual Relative Error Normalized
*  By Dividing By The Tolerance.

*  The Array Etest Specifies The Error Test To Be Applied To Each
*  Error Measure.

*  On Exit, The Logical Flag Errok
*   -- Is .false. If Any Of The Error Measures Exceeds The
*      Corresponding Value In The Array Etest
*   -- Is .true. If All Error Measures Are Less Than The
*      Corresponding Values Of Etest.

      Errok = .true.

      Do 10 Im = 1, Nmsh
      Do 10 It = 1, Ntol
         Icmp = Ltol(It)
         Er = U(Icmp,Im) - Uold(Icmp,Im)
         Denom = Max(One, Abs(Uold(Icmp,Im)))
         Errel = Abs(Er/(Tol(It)*Denom))
         If (Errel .gt. Etest(It)) Errok = .false.
   10 Continue

      Return
      End



      SUBROUTINE acgetptq( debug, mfsrch, nout, alfmax, alfsml, alfuzz,
     *                   epsaf, epsag, eta, ftry, oldf, oldg,
     *                   rmu, tolabs, tolrel, toltny, imprvd,
     *                   inform, nfsrch, alfa, alfbst, fbest,
     *                   braktd, crampd, extrap,vset,wset,nsamea,nsameb,
     *                   a, b, fa, factor, fv, fw, xtry, xv, xw )

      implicit double precision (a-h,o-z)
      logical            debug, imprvd
      logical            braktd, crampd, extrap, vset, wset
      integer            mfsrch, nout, inform, nfsrch, nsamea, nsameb
c
c  *********************************************************************
c  getptq  is a step-length algorithm for minimizing a function of one
c  variable.  it will be called repeatedly by a search routine whose
c  purpose is to estimate a point  alfa = alfbst  that minimizes some
c  function  f(alfa)  over the closed interval (0, alfmax).
c
c  getptq  requires the function  f(alfa)  (but not its gradient)
c  to be evaluated at various points within the interval.  new
c  step-length estimates are computed using quadratic interpolation with
c  safeguards.
c
c  reverse communication is used to allow the calling program to
c  evaluate  f.  some of the parameters must be set or tested
c  by the calling program.  the remainder would ordinarily be local
c  variables.
c
c
c  input parameters (relevant to the calling program)
c  --------------------------------------------------
c
c  debug         specifies whether detailed output is wanted.
c
c  inform        must be nonzero on the first entry (e.g., -1).
c                it will be altered by  getptq  for later entries.
c
c  mfsrch        is an upper limit on the number of times  getptq  is
c                to be entered consecutively with  inform = 0
c                (following an initial entry with  inform lt 0).
c
c  nout          is the file number to be used for printed output
c                if debug is true.
c
c  alfa          is the first estimate of the step length.  alfa  is
c                subsequently altered by  getptq  (see below).
c
c  alfmax        is the upper limit of the interval to be searched.
c
c  alfsml        is intended to prevent inefficiency when the optimum
c                step is very small, for cases where the calling
c                program would prefer to re-define  f(alfa).  alfsml is
c                allowed to be zero. early termination will occur if
c                getptq  determines that the optimum step lies
c                somewhere in the interval  (0, alfsml)  (but not if
c                alfmax .le. alfsml).
c
c  epsaf         is an estimate of the absolute precision in the
c                computed values of  f.
c
c  eta           controls the accuracy of the search.  it must lie
c                in the range   0.0  le  eta  lt  1.0.  decreasing
c                eta  tends to increase the accuracy of the search.
c
c  oldf          is the value of  f(0).
c
c  oldg          is an estimate of the gradient of  f  at  alfa = 0.
c                it should be non-positive.
c
c  rmu           controls what is meant by a significant decrease in  f.
c                the final  f(alfbst)  should lie on or below the line
c                      l(alfa)  =  oldf + alfa*rmu*oldg.
c                rmu  should be in the open interval (0, 0.5).
c                the value  rmu = 1.0d-4  is good for most purposes.
c
c  tolabs,tolrel define a function  tol(alfa) = tolrel*alfa + tolabs
c                such that if  f  has already been evaluated at step
c                alfa,  then it will not be evaluated at any point
c                closer than  tol(alfa).
c                these values may be reduced by  getptq  if they seem
c                to be too large.
c
c  toltny        is the smallest value that  tolabs  is allowed to be
c                reduced to.
c
c
c  output parameters (relevant to the calling program)
c  ---------------------------------------------------
c
c  imprvd        is true if the previous step  alfa  was the best
c                point so far.  any related quantities (e.g., arrays)
c                should be saved by the calling program before paying
c                attention to  inform.
c
c  inform = 0    means the calling program should evaluate
c                           ftry = f(alfa)
c                for the new trial step  alfa,  and then re-enter
c                getptq.
c
c  inform = 1    means the search has terminated successfully
c                with a step  alfbst  that is less than the
c                upper bound  alfmax.
c
c  inform = 2    means the search has terminated successfully
c                with a step  alfbst  that is equal to the
c                upper bound  alfmax.
c
c  inform = 3    means that the search failed to find a point of
c                sufficient decrease in  mfsrch  functions, but an
c                improved point was found.
c
c  inform = 4    means  alfmax  is so small that a search should
c                not have been done.
c
c  inform = 5    means that the search was terminated prematurely
c                because of the value of  alfsml  (see above).
c
c  inform = 6    means the search has failed to find a useful step.  if
c                the subroutine for the function and gradient has been
c                programmed correctly, this will usually occur if the
c                minimum lies very close to  alfa = 0  or the gradient
c                is not sufficiently accurate.
c
c  inform = 7    means that the value of  g(0) was positive on entry.
c
c  alfa          is the step at which the next function value must be
c                computed.
c
c  alfbst        should be accepted by the calling program as the
c                required step-length estimate, whenever  getptq
c                returns  inform = 1,  2  or  3.
c
c  fbest         will be the corresponding value of  f.
c
c
c  the following parameters retain information between entries
c  -----------------------------------------------------------
c
c  alfuzz        is such that, if the final  alfa  lies in the interval
c                (0,alfuzz)  and  abs( f(alfa)-oldf ) le epsaf,  alfa
c                cannot be guaranteed to be a point of sufficient
c                decrease.
c
c  braktd        is false if  f  has not been evaluated at the far end
c                of the interval of uncertainty.  in this case, the
c                point  b  will be at  alfmax + tol(alfmax).
c
c  crampd        is true if  alfmax  is very small (le tolabs).
c                if the search fails, this indicates that a zero
c                step should be taken.
c
c  extrap        is true if alfbst has moved at least once and  xv
c                lies outside the interval of uncertainty.  in this
c                case, extra safeguards are applied to allow for
c                instability in the polynomial fit.
c
c  vset          records whether a third-best point has been
c                determined.
c
c  wset          records whether a second-best point has been
c                determined.  it will always be true by the
c                time the convergence test is applied (label 300).
c
c  nsamea        is the number of consecutive times that the left-hand
c                end of the interval of uncertainty has remained the
c                same.
c
c  nsameb        similarly for the right-hand end.
c
c  a, b, alfbst  define the current interval of uncertainty.
c                the required minimum lies somewhere within the
c                closed interval  (alfbst + a, alfbst + b).
c
c  alfbst        is the best point so far.  it is strictly within the
c                the interval of uncertainty except when it lies at the
c                left-hand end when  alfbst  has not been moved.
c                hence we have    a le 0,   b gt 0.
c
c  fbest         is the value of  f  at the point  alfbst.
c
c  fa            is the value of  f  at the point  alfbst + a.
c
c  factor        controls the rate at which extrapolated estimates of
c                alfa  may expand into the interval of uncertainty.
c                factor is not used if the minimum has been bracketed
c                (i.e., when the variable  braktd  is true).
c
c  fv, fw        are the values of  f  at the points  alfbst + xv,
c                alfbst + xw.  they are not defined until  vset
c                or  wset  (respectively) is true.
c
c  ftry          is the value of  f  at the new point  alfbst + xtry.
c
c  xtry          is the trial point within the shifted interval (a, b).
c                the new trial function value must be computed at the
c                point  alfa  =  alfbst + xtry.
c
c  xv            is such that  alfbst + xv  is the third-best point.
c                it is not defined until  vset  is true.
c
c  xw            is such that  alfbst + xw  is the second-best point.
c                it is not defined until  wset  is true.
c                in some cases,  xw  will replace a previous  xw  that
c                has a lower function but has just been excluded from
c                the interval of uncertainty.
c
c
c  systems optimization laboratory, stanford university, california.
c  original version february 1982.  rev. may 1983.
c  *********************************************************************
c
      logical            closef, conv1, conv2, conv3, convrg
      logical            moved, sigdec, xinxw
      data               zero, point1, half/ 0.0d+0,  0.1d+0, 0.5d+0/
      data                one,  two,  five   / 1.0d+0, 2.0d+0,  5.0d+0/
      data                ten, eleven      /10.0d+0, 11.0d+0        /
c
c
c  local variables
c  ---------------
c
c  closef        is true if the worst function  fv  is within  epsaf
c                of  fbest  (up or down).
c
c  convrg        will be set to true if at least one of the convergence
c                conditions holds at  alfbst.
c
c  moved         is true if a better point has been found (alfbst gt 0).
c
c  sigdec        says whether  fbest  represents a significant decrease
c                in the function, compared to the initial value  oldf.
c
c  xinxw         is true if  xtry  is in  (xw,0)  or  (0,xw).
c  ---------------------------------------------------------------------
c

      imprvd = .false.
      if (inform .ne. -1) go to 100
c
c  ---------------------------------------------------------------------
c  first entry.  initialize various quantities, check input data and
c  prepare to evaluate the function at the initial step  alfa.
c  ---------------------------------------------------------------------
      nfsrch = 0
      alfbst = zero
      fbest  = oldf
      if (oldg   .gt.      zero) go to 970
      if (oldg   .ge. (- epsag)) go to 960
      if (alfmax .le.    toltny) go to 940
c
      braktd = .false.
      crampd = alfmax .le. tolabs
      extrap = .false.
      vset   = .false.
      wset   = .false.
      nsamea = 0
      nsameb = 0
      alfuzz = two*epsaf/(rmu*abs( oldg ))
      a      = zero
      b      = alfmax + (tolrel*alfmax + tolabs)
      fa     = oldf
      factor = five
      tol    = tolabs
      xtry   = alfa
      go to 800
c
c  ---------------------------------------------------------------------
c  subsequent entries.
c  the function has just been evaluated at  alfa = alfbst + xtry,
c  giving  ftry.
c  ---------------------------------------------------------------------
 100   nsamea = nsamea + 1
      nsameb = nsameb + 1
      xtry   = alfa - alfbst
      moved  = alfbst .gt. zero
c
c  check if  xtry  is in the interval  (xw,0)  or  (0,xw).
c
      xinxw  = .false.
      if (wset) xinxw =       zero .lt. xtry  .and.  xtry .le. xw
     *                  .or.    xw .le. xtry  .and.  xtry .lt. zero
c
c  see if the new step is better.
c
      deltaf = ftry   - oldf
      ctry   = deltaf - alfa*rmu*oldg
      if (alfa .le. alfuzz) sigdec = deltaf .le. (- epsaf)
      if (alfa .gt. alfuzz) sigdec = ctry   .le.    epsaf
      imprvd = sigdec  .and.  ( ftry - fbest ) .le. (- epsaf)
c
      if (.not. imprvd) go to 130
c
c  we seem to have an improvement.  the new point becomes the
c  origin and other points are shifted accordingly.
c
      if (.not. wset) go to 110
      xv     = xw - xtry
      fv     = fw
      vset   = .true.
 110   xw     = zero - xtry
      fw     = fbest
      wset   = .true.
      fbest  = ftry
      alfbst = alfa
      a      =    a - xtry
      b      =    b - xtry
      moved  = .true.
      extrap = .not. xinxw
c
c  decrease the length of the interval of uncertainty.
c
      if (xtry .lt. zero) go to 120
      a      = xw
      fa     = fw
      nsamea = 0
      go to 300
 120   b      = xw
      nsameb = 0
      braktd = .true.
      go to 300
c
c  the new function value is no better than the best point found so far.
c  the point  xtry  must be a new end point of the interval of
c  uncertainty.
c
 130   if (xtry .ge. zero) go to 140
      a      = xtry
      fa     = ftry
      nsamea = 0
      go to 150
 140   b      = xtry
      nsameb = 0
      braktd = .true.
c
c  the origin remains unchanged but  xtry  may qualify as  xw.
c
 150   if (.not. wset)   go to 160
      if ((ftry - fw) .gt. epsaf) go to 170
      xv     = xw
      fv     = fw
      vset   = .true.
 160   xw     = xtry
      fw     = ftry
      wset   = .true.
      if (moved) extrap = xinxw
      go to 300
c
c  ftry  is no better than  fbest  or  fw.  if the best point has not
c  been moved, there must be more than one minimum.
c
 170   if (moved) go to 175
      xw     = xtry
      fw     = ftry
      go to 300
c
c  ftry  is no better than  fbest  or  fw,  but  xtry  may become  xv.
c  extrap  has the value set in the previous entry.
c
 175   if (.not. vset) go to 180
      if ((ftry - fv) .gt. epsaf  .and.  extrap) go to 300
 180   if (xinxw) go to 190
      xv     = xtry
      fv     = ftry
      vset   = .true.
      go to 300
 190   if (vset) xw = xv
      if (vset) fw = fv
      xv     = xtry
      fv     = ftry
      vset   = .true.
c
c  ---------------------------------------------------------------------
c  check the termination criteria.
c  ---------------------------------------------------------------------
 300   tol    = tolrel*alfbst + tolabs
      deltaf = fbest  - oldf

      cbest  = deltaf - alfbst*rmu*oldg
      if (alfbst .le. alfuzz) sigdec = deltaf .le. (- epsaf)
      if (alfbst .gt. alfuzz) sigdec = cbest  .le.    epsaf
      closef = .false.
      if (vset) closef = abs( fbest - fv ) .le. epsaf
c
      conv1  = max( abs( a ), b )  .le.  (tol + tol)

*  conv2 changed by mhw, 20 sept 1992, to allow it to be
*  satified for any significant decrease in f
*      conv2  =  moved  .and.  sigdec
*     *                 .and.  abs( fa - fbest )  .le.  a*eta*oldg
      conv2 = moved .and. sigdec
      conv3  = closef  .and.  (sigdec  .or.
     *                        (.not. moved)  .and.  (b .le. alfuzz))
      convrg = conv1  .or.  conv2  .or.  conv3
c
      atrue  = alfbst + a
      btrue  = alfbst + b
      alfaw  = alfbst + xw
      gap    = b - a
      if (vset) alfav  = alfbst + xv
      if (convrg  .and.  moved) go to 910
c
c  exit if the step is too small.
c
      if (btrue   .lt.  alfsml) go to 950

      if (nfsrch  .ge.  mfsrch) go to 930
      if (.not. convrg) go to 400
c
c  a better point has not yet been found (the step  xw  is no better
c  than step  zero).  check that the change in  f  is consistent with a
c  perturbation in  x  of  tol, the estimate of the minimum spacing
c  constant.  if the change in  f  is larger than  epsaf,  the value
c  of  tol  is reduced.
c
      tol    = xw/ten
      tolabs = tol
      if (abs(fw - oldf) .gt. epsaf  .and.  tol .gt. toltny) go to 400
      if (crampd) go to 940
      go to 960
c
c  ---------------------------------------------------------------------
c  proceed with the computation of a trial step length.
c  the choices are...
c  1. parabolic fit using function values only.
c  2. damped parabolic fit if the regular fit appears to be
c     consistently over-estimating the distance to the minimum.
c  3. bisection, geometric bisection, or a step of  tol  if the
c     parabolic fit is unsatisfactory.
c  ---------------------------------------------------------------------
 400   xmidpt = half*(a + b)
      q      = zero
      s      = zero
c
c  ---------------------------------------------------------------------
c  fit a parabola.
c  ---------------------------------------------------------------------
c
c  check if there are two or three points for the parabolic fit.
c
      gw = (fw - fbest)/xw
      if (vset  .and.  moved) go to 450
c
c  only two points available.  use  fbest,  fw  and the derivative
c  oldg.
c
      if (.not. moved) s = oldg
      if (      moved) s = oldg - two*gw
      q = two*(oldg - gw)
      go to 600
c
c  three points available.  use  fbest,  fw  and  fv.
c
 450   gv = (fv - fbest)/xv
      s  = gv - (xv/xw)*gw
      q  = two*(gv - gw)
c
c  ---------------------------------------------------------------------
c  construct an artificial interval  (artifa, artifb)  in which the
c  new estimate of the step length must lie.  set a default value of
c  xtry  that will be used if the polynomial fit is rejected.  in the
c  following, the interval  (a,b)  is considered the sum of two
c  intervals of lengths  dtry  and  daux, with common end point at the
c  best point (zero).  dtry  is the length of the interval into which
c  the default  xtry  will be placed and  endpnt  denotes its non-zero
c  end point.  the magnitude of  xtry  is computed so that the exponents
c  of  dtry  and  daux  are approximately bisected.
c  ---------------------------------------------------------------------
 600   artifa = a
      artifb = b
      if (braktd) go to 610
c
c  the minimum has not been bracketed.  set an artificial upper bound
c  by expanding the interval  xw  by a suitable factor.
c
      xtry   = - factor*xw
      artifb =   xtry
      if (alfbst + xtry .lt. alfmax) factor = five*factor
      go to 700
c
c  the minimum has been bracketed.
c  if the gradient at the origin is being used for the
c  polynomial fit, the default  xtry  is one tenth of  xw.
c
 610   if (vset  .and.  moved) go to 620
      xtry   = xw/ten
      go to 700
c
c  three points exist in the interval of uncertainty.  check whether
c  the points are configured for an extrapolation or interpolation.
c
 620   if (extrap) go to 660
c
c  if the interpolation appears to be consistently over-estimating the
c  distance to the minimum,  damp the interpolation step.
c
      if (nsamea .lt. 3  .and.  nsameb .lt. 3) go to 630
      factor = factor / five
      s      = factor * s
      go to 640
 630   factor = one
c
c  the points are configured for an interpolation.  the artificial
c  interval will be just  (a,b).   set  endpnt  so that  xtry
c  lies in the larger of the intervals  (a,0)  and  (0,b).
c
 640    if (xmidpt .lt. zero) endpnt = a
      if (xmidpt .gt. zero) endpnt = b
c
c  if a bound has remained the same for three iterations, set  endpnt
c  so that  xtry  is likely to replace the offending bound.
c
      if (nsamea .ge. 3) endpnt = a
      if (nsameb .ge. 3) endpnt = b
      go to 680
c
c  the points are configured for an extrapolation.
c
 660   if (xw .lt. zero) endpnt = b
      if (xw .gt. zero) endpnt = a
c
c  compute the default value of  xtry.
c
 680   dtry = abs( endpnt )
      daux = gap - dtry
      if (daux .ge. dtry)   xtry = five*dtry*(point1 + dtry/daux)/eleven
      if (daux .lt. dtry)   xtry = half*sqrt( daux )*sqrt( dtry )
      if (endpnt .lt. zero) xtry = - xtry
c
c  if the points are configured for an extrapolation set the artificial
c  bounds so that the artificial interval lies strictly within  (a,b).
c  if the polynomial fit is rejected,  xtry  will remain at the relevant
c  artificial bound.
c
      if (extrap  .and.  xtry .le. zero) artifa = xtry
      if (extrap  .and.  xtry .gt. zero) artifb = xtry
c
c  ---------------------------------------------------------------------
c  the polynomial fits give  (s/q)*xw  as the new step.
c  reject this step if it lies outside  (artifa, artifb).
c  ---------------------------------------------------------------------
 700   if (q .eq. zero) go to 800
      if (q .lt. zero) s = - s
      if (q .lt. zero) q = - q
      if (s*xw .lt. q*artifa   .or.   s*xw .gt. q*artifb) go to 800
c
c  accept the polynomial fit.
c
      xtry = zero
      if (abs( s*xw ) .ge. q*tol) xtry = (s/q)*xw
c
c  ---------------------------------------------------------------------
c  test for  xtry  being larger than  alfmax  or too close to  a  or  b.
c  ---------------------------------------------------------------------
 800   if (braktd) go to 810
c
c  if the step is close to or larger than  alfmax,  replace it by
c  alfmax  (to force evaluation of the function at the boundary).
c
      alfa   = alfbst + xtry
      if (alfmax - alfa .gt. (tolrel*alfmax + tolabs)) go to 810
      braktd = .true.
      xtry   = alfmax - alfbst
      alfa   = alfmax
      go to 900
c
c  otherwise, the function must not be evaluated at a point too close
c  to  a  or  b.  (it has already been evaluated at both those points.)
c
 810   xmidpt = half*(a + b)
      if (xtry .gt. a + tol  .and.  xtry .lt. b - tol) go to 820
      if (xmidpt .gt. zero) xtry =   tol
      if (xmidpt .le. zero) xtry = - tol
c
c
c  f  must not be calculated too close to  alfbst.
c
 820   if (abs( xtry ) .lt. tol  .and.  xmidpt .lt. zero) xtry = - tol
      if (abs( xtry ) .lt. tol  .and.  xmidpt .ge. zero) xtry =   tol
      alfa   = alfbst + xtry
c
c  ---------------------------------------------------------------------
c  exit.
c  ---------------------------------------------------------------------
c
c  new function value required.
c
 900   inform = 0
      go to 990
c
c  convergence test satisfied.
c
 910   inform = 1
      if (alfa .eq. alfmax) inform = 2
      go to 990
c
c  mfsrch  function evaluations without sufficient decrease, but an
c  improved point was found.
c
 930   if (.not. moved) go to 960
      inform = 3
      go to 990
c
c  zero step (alfmax too small).
c
 940   inform = 4
      go to 990
c
c  premature termination.  the step is smaller than  alfsml.
c
 950   inform = 5
      go to 990
c
c  zero step (a sufficiently better point could not be found).
c
 960   inform = 6
      go to 990
c
c  zero step (positive gradient at the starting point).
c
 970   inform = 7
c
c  exit.
c
 990  return
c
c
c  end of getptq
      end



      SUBROUTINE acInterp(Ncomp, Nmsh, Xx, Nudim, U,
     *                    Nuold_dim,Nmold, Xxold, Uold)

      Implicit Double Precision (A-H, O-Z)
      Dimension Xx(*), U(Nudim,*), Xxold(*), Uold(Nuold_dim,*)

* Blas: Dcopy

      Parameter (Zero = 0.0d+0)

*  Interp Performs Piecewise Linear Interpolation Of The Old
*  Solution Uold At The Nmold Old Mesh Points Xxold Onto The Nmsh
*  New Mesh Points Xx, Producing An Interpolated Solution U.
*  Note That No Assumption Is Made That The New Mesh Has
*  More Points Than The Old, Nor That The New And Old Mesh
*  Points Are Related In A Specific Way (Except That Their First
*  And Last Points Are Identical).

*  By Construction, Xx(1) = Xxold(1).  Copy The First Ncomp
*  Components Of Uold Into Those Of U.


      Call Dcopy(Ncomp, Uold(1,1), 1, U(1,1), 1)

      I = 2
      Do 100 Im = 2, Nmsh-1

   50    Continue
         If (I .gt. Nmold) Return

*  Check Whether The Im-Th Point In The New Mesh Lies Strictly
*  To The Right Of, Or To The Left Of (Or Exactly On) The
*  I-Th Point In The Old Mesh.


         If (Xx(Im) .gt. Xxold(I)) Then
            I = I + 1
            Go To 50
         Else
            Xdif = Xxold(I) - Xx(Im)
            If (Xdif .eq. Zero) Then

*  Xx(Im) And Xxold(I) Are Identical.

               Call Dcopy(Ncomp, Uold(1,I), 1, U(1,Im), 1)
               I = I + 1
            Else
               Xint = Xxold(I) - Xxold(I-1)
               Xrat = Xdif/Xint
               Do 70 K = 1, Ncomp
                  U(K,Im) = Uold(K,I) + Xrat*(Uold(K,I-1)-Uold(K,I))
   70          Continue
            Endif
         Endif

  100 Continue
      Call Dcopy(Ncomp, Uold(1,Nmold), 1, U(1,Nmsh), 1)
      Return
      End







       SUBROUTINE acCONDESTIM(ALEFT,ARIGHT,NMSH,NCOMP,N,XX,TOPBLK,
     *            NRWTOP,NOVRLP,ARRAY,
     *          NRWBLK,NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,OMG,
     *          C1,WORK,KPPA,GAMMA,SIGMA,CKAPPA,CKAPPA2)

C     **************************************************************
C
C     COMPUTES THE FIRST AND LAST BLOCK COLUMN OF INVERSE MATRIX AND
C     ESTIMATE THE CONDITION NUMBERS KPPA, GAMMA
C
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        IMPLICIT NONE
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,WORK,C1,
     *          OMG,GAMMA,gamma1,KPPA,MINMG, GAMMAK,KPPAK,CKAPPA1,
     *          KPPAI, GAMMAI, CKAPPA,CKAPPA2,kappa1_n, kappa2_n,ckmax
        DOUBLE PRECISION BOMEGA1,BOMEGA2, SIGMA, SIGMAK,KPPAJ
        DOUBLE PRECISION ALEFT,ARIGHT,XX,  CSUM, ZNORM
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,IPVCD
        INTEGER NCOMP,NMSH,idmx,idmn,idamax,idomg,job
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),OMG(*),
     *          BOTBLK(NRWBOT,*),WORK(*),C1(ncomp,*),XX(*),
     *          IPVCD(*)
        INTEGER k,i,j,l, kn,indnm,itmax,its, indsk
        LOGICAL FIRSTCOL

        INTRINSIC SIGN
        DOUBLE PRECISION SIGN

        FIRSTCOL=.true.

        MINMG=0.0d+0
        KPPA=0.0d+0

        IF (FIRSTCOL) THEN


        do 300 k=1,ncomp
           do 330 l=1,N
              WORK(l)=0.0d0
              if (k.le.NRWTOP) then
                   WORK(k)=1.0d0
              else
                   WORK(N-ncomp+k)=1.0d0
              endif
 330       continue


           CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *                 NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,
     *                 WORK,0)


           do 340 l=1,N
               C1(k,l)=WORK(l)
 340       continue

 300    continue

        kn = ncomp
      ELSE

        kn = 2
        do  k=1,2
           do  l=1,N
              WORK(l)=0.0d0
           end do
           if (k==1) then
              do l=1, NRWTOP
                   WORK(l)=1.0d0/ncomp
              end do
              do l=N-ncomp+nrwtop+1,N
                   WORK(l)=1.0d0/ncomp
              enddo
            else
              j=1
              ckmax = 0.0d0
              do l=1, NRWTOP
                 WORK(l)=(-1.0d0)**(ncomp)*(1.0d0+(j-1.0d0)/(ncomp-1))
                 ckmax = max(abs(WORK(l)),ckmax)
                 j=j+1
              end do
              do l=N-ncomp+nrwtop+1,N
                 WORK(l)=(-1.0d0)**(ncomp)*(1.0d0+(j-1.0d0)/(ncomp-1))
                  j=j+1
                  ckmax = max(abs(WORK(l)),ckmax)
              enddo
              do l=1, NRWTOP
                 WORK(l)=WORK(l)/ckmax
              end do
              do l=N-ncomp+nrwtop+1,N
                 WORK(l)=WORK(l)/ckmax
              enddo



            end if

           CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *                 NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,
     *                 WORK,0)


           do  l=1,N
               C1(k,l)=WORK(l)
          end do

      end do

      end if

C
C     ESTIMATION OF CONDITION NUMBER BY BRUGNANO,TRIGIANTE & MAZZIA
C
c     infinity-norm
c
        l=1
        OMG(l) = 0d0
        j=1
        do  i = j,j+ncomp-1
            CSUM = 0d0
             do k=1,kn
               CSUM = CSUM + DABS(C1(k,i))
             end do
c             OMG(l) = DMAX1( OMG(l),CSUM)
             if (OMG(l).LE.CSUM) then
                  OMG(l)=CSUM
                  idmn=i
             endif
        end do
        KPPA = OMG(l)
        idmx=idmn
        GAMMA = 0.0d0
        l=l+1

        do 400 j=ncomp+1,N-ncomp+1,ncomp
           OMG(l) = 0d0
           do  i = j,j+ncomp-1
               CSUM = 0d0
               do k=1,kn
                 CSUM = CSUM + DABS(C1(k,i))
               end do
c               OMG(l) = DMAX1( OMG(l),CSUM)
               if (OMG(l).LT.CSUM) then
                  OMG(l)=CSUM
                  idmn=i
                endif
           end do
           if (OMG(l).GT.KPPA) then
               KPPA=OMG(l)
               idmx=idmn
           endif
           if (OMG(l).GT.OMG(l-1)) then
               gamma1=OMG(l)* (XX(l)-XX(l-1))
           else
               gamma1=OMG(l-1)* (XX(l)-XX(l-1))
           end if
           GAMMA=GAMMA + gamma1
           l=l+1
 400    continue
      GAMMA=GAMMA/(ARIGHT-ALEFT)



      GAMMAI =0d0
      KPPAI  =0d0
      SIGMA  =0d0

c      DO l=l,NBLOKS+1
c        OMG(l)=0.0d0
c      END DO
      DO k=1,kn
         I=2
         GAMMAK = 0d0
         DO l=1,ncomp
            WORK(l)=C1(k,l)
         END DO

         BOMEGA1 = 0d0
         DO l=1,ncomp
            BOMEGA1 = DMAX1(BOMEGA1,ABS(C1(k,l)))
         END DO
         KPPAK   = BOMEGA1
c        OMG(I-1)  = OMG(I-1)+BOMEGA1
         DO j=ncomp+1,N-ncomp+1,ncomp
               BOMEGA2 = 0d0
               DO l=1,ncomp
                  BOMEGA2 = DMAX1(BOMEGA2,ABS(C1(k,j+l-1)))
               END DO
c               OMG(I)  = OMG(I)+BOMEGA2
               IF (BOMEGA1 .GT. BOMEGA2) THEN
                  GAMMAK = GAMMAK + BOMEGA1*(XX(I)-XX(I-1))
               ELSE
                  GAMMAK = GAMMAK + BOMEGA2*(XX(I)-XX(I-1))
               END IF
               IF (BOMEGA2 .GT. KPPAK) THEN
                   KPPAK=BOMEGA2
               END IF
               BOMEGA1=BOMEGA2
               I=I+1
         END DO
           GAMMAK = GAMMAK/(ARIGHT-ALEFT)
           IF (GAMMAK .GT. GAMMAI) THEN
              GAMMAI = GAMMAK
           END IF
           SIGMAK = KPPAK/GAMMAK

           IF (SIGMAK .GT.  SIGMA) THEN
               SIGMA = SIGMAK
               indsk = k
           END IF
           IF (KPPAK .GT. KPPAI) THEN
              KPPAI = KPPAK
           END IF
       END DO

c nuovo OMG
c      DO l=l,NBLOKS+1
c        OMG(l)=0.0d0
c      END DO

c         k=indsk
c         I=2
c         BOMEGA1 = 0d0
c         DO l=1,ncomp
c            BOMEGA1 = DMAX1(BOMEGA1,ABS(C1(k,l)))
c         END DO
c         OMG(I-1)  = BOMEGA1
c         DO j=ncomp+1,N-ncomp+1,ncomp
c               BOMEGA2 = 0d0
c               DO l=1,ncomp
c                  BOMEGA2 = DMAX1(BOMEGA2,ABS(C1(k,j+l-1)))
c               END DO
c               OMG(I) = BOMEGA2
c               I=I+1
c         END DO

c fine nuovo OMG

       CKAPPA=KPPA
       CKAPPA1=KPPA
       CKAPPA2=0.0d0
       itmax = 0

 450    if (sigma .lt. 10d0 .and. itmax .eq. 0 ) then
         do l=1,N
            work(l)=1.0d0/dble(N+10)
          end do
          work(idmx)=(1.0d0+10.0d0)/dble(N+10)
          its = 1
       else
           do 380 l=1,N
              WORK(l)=0.0d0
 380       continue
           WORK(idmx)=1.0d0
           its = 0
       endif

       CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *                 NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,
     *                 WORK,1)

             DO i=1,NBLOKS
               DO j=1,NRWBLK
                  WORK( (i-1)*NRWBLK+NRWTOP+j ) =
     *             (XX(i+1)-XX(i))*WORK( (i-1)*NRWBLK+NRWTOP+j)
               ENDDO
             END DO
           kappa1_n = 0.0d0
           DO l=1,nrwtop
              kappa1_n = kappa1_n+ abs(WORK(l))
           END DO

           DO l=N-ncomp+nrwtop+1,N
              kappa1_n = kappa1_n+ abs(WORK(l))
           END DO
           kappa2_n = 0.0d0
           DO l=nrwtop+1,N-ncomp+nrwtop
              kappa2_n = kappa2_n+ abs(WORK(l))
           END DO

         IF ( kappa1_n + kappa2_n .ge. CKAPPA) THEN
           CKAPPA2 = kappa2_n
           CKAPPA = kappa1_n+kappa2_n
           CKAPPA1 = kappa1_n
         ELSEIF (itmax .ne. 0) THEN
            goto 500
         ELSE


         END IF

       if ( (SIGMA .lt. 10 .or. ckappa2 .GT. 0.1d0*KPPA)
     *                                   .and. itmax.le.2 ) then


          itmax = itmax+1
          DO I = 1,N
            WORK(I) = SIGN(1.0d0,WORK(I))
          END DO

            DO i=1,NBLOKS
               DO j=1,NRWBLK
                  WORK( (i-1)*NRWBLK+NRWTOP+j ) =
     *             (XX(i+1)-XX(i))*WORK( (i-1)*NRWBLK+NRWTOP+j)
               ENDDO
            END DO
            CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,
     *             WORK,0)


      ZNORM=0.0d0
      DO I=1,N
         IF ( ABS(WORK(I)) .GT. ZNORM) THEN
            ZNORM = ABS(WORK(I))
            indnm = I
         END IF
      END DO


      IF (ZNORM .GT. ABS(work(idmx)) .or. its .eq. 1) THEN
        idmx=indnm
        goto 450
      END IF

      END IF




 500     IF (.NOT. FIRSTCOL) THEN
             KPPA = ckappa1
             GAMMA = GAMMAI
         END IF

         RETURN
       END


      SUBROUTINE acSelconderrmsh(Ncomp, Nmsh, Ntol, Ltol, Tol,
     *     Nfxpnt, Fixpnt, Ipow, Nmax, Nvold,
     *     Xx, Nudim, U, Ermeas, Irefin, Ihcomp,
     *     Nmold, Xxold, Ermx, dDouble, Maxmsh, Phiold, Voldmsh,
     *     r4,amg,stab_cond, stiff_cond,linear)

      Implicit Double Precision (A-H,O-Z)

      Dimension  Ltol(Ntol), Tol(Ntol), Fixpnt(*)
      Dimension  Xx(*), U(Nudim, *), Ermeas(Ncomp,*)
      Dimension  Irefin(Nmsh-1), Ihcomp(Nmsh-1)
      Dimension  Xxold(*), Ermx(*), Phiold(*), Voldmsh(*)
      Dimension  Npr(2), Pmax(2), Hord(2)
      Dimension  r4(*), amg(*)
      Logical    dDouble, Maxmsh, stab_cond, stiff_cond, add,linear

      Intrinsic Abs
      Intrinsic Max
      Intrinsic Int

*  Blas: Dcopy

      Common /acFlags/ Ifinal,Iatt,Iback,Iprec,Iucond
      Common /acMshvar/ Hsml,Npr,Pmax,Hord

      Parameter  ( Zero = 0.0d+0, One = 1.0d+0, Onep1 = 1.1d+0 )
      Parameter  ( Erdcid = 5.0d+0 )
      Parameter  ( Phitst = 0.1d+0 )

*  The Routine Selmsh Performs Selective Mesh Refinement, Depending
*  On The Error Measure Ermeas.


      Maxmsh = .false.

      Frcpow = One/Ipow
      dDouble = .false.
      Nmold = Nmsh
      Ninter = Nmsh - 1
      Iprec = Min(Iprec,1)

*  Copy The Current Mesh Into The Xxold Array.

      If (Iatt .eq. -1) Nvold = Nmsh
      If (Iatt .eq. 0) Call Dcopy(Nvold, Xxold, 1, Voldmsh, 1)

      Call Dcopy(Nmold, Xx, 1, Xxold, 1)
      Ithres = 0


      Thres = One

*  On Input, The Array Ermeas Represents Some Error Measure Defined
*  Over The Components And Mesh Intervals (Not Mesh Points).
*  It Is Normalized In The Following Loop With Respect To The
*  Tolerance Array And The Current Solution.
*  The Value Errmax Gives The Maximum Normalized Error.

      Errmax = Zero
      Do 120 Im = 1, Ninter
         Ermx(Im) = Zero
         Do 110 It = 1, Ntol
            Jcomp = Ltol(It)
            Umin = Min(Abs(U(Jcomp,Im)),Abs(U(Jcomp,Im+1)))
            Denom = Tol(It)*Max(One, Umin)
            Ems = Ermeas(Jcomp,Im)
            Ermeas(Jcomp,Im) = Abs(Ems)/Denom
            Err = Ermeas(Jcomp, Im)
            If (Err .ge. Ermx(Im)) Then
                Ermx(Im) = Err
                Ihcomp(Im) = Jcomp
            Endif
  110    Continue
         Errmax = Max(Ermx(Im), Errmax)
  120 Continue




      call acmoncondmsh_l(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,nptcond,
     *                 r4,amg,linear)




cf   the conditioning and the error

 200  Continue

*  For Each Interval Im,  The Integer Irefin(Im) Is Calculated
*  Based On An Equidistribution Procedure Involving The
*  Threshold Value Thres.  If Irefin(Im) > 1, We Add
*  Points To Interval Im.  If Irefin(Im) = 1, We Check Whether
*  Point Im Can Be Removed From The Mesh.

*  Nmest Is A Lower Bound On The Number Of Points In The New Mesh.
*  We Do Not Know In Advance How Many Points Will Be Removed,
*  So Nmest Is Computed By Assuming All Eligible Points Are Removed.

      Philrg = 0.d0
      Esum = 0.d0
      Nmest = Nmsh
      Mshchng = 0
      Do 220 Im = 1, Ninter
         Errim = Ermx(Im)**Frcpow
         If (Ermx(Im) .ge. Thres) Then
            Irefin(Im) = Int(Errim) + 1
            Nmest = Nmest + Irefin(Im) - 1
            Mshchng = 1
         Else
            Irefin(Im) = 1
            Nmest = Nmest - 1
         Endif
         Him = Xxold(Im+1)-Xxold(Im)
         Phiim = Errim/Him
         Esum = Esum+Errim
         If (Iatt .eq. -1) Then
            Phiold(Im) = Phiim
         Elseif (Iatt .eq. 0 .and. Phiim .gt. Philrg) Then
           Philrg = Phiim
           Hordlrg = Errim
           Imreg = Im
         Endif

  220 Continue



c       if ( .not. stab_cond) then
c      if (nptcond .ge. 4 .and. .not. stab_cond)then
c      if (nptcond .ge. 4  .and.(.not. stab_cond  .or. .not. linear))then

      if (nptcond .ge. 4  )then

       do 221 im = 1, ninter
         if  (r4(im) .ge. fatt_r1r3) then
               irefin(im) = max(nptcond, irefin(im))
               if (nptcond .gt. irefin(im)) then
                    nmest = nmest + (nptcond-irefin(im)) - 1
               endif
         elseif (r4(im) .lt. fatt_r3
     *           .and. ifinal .eq. 0  ) then
            nmest = nmest -  max( 1 , floor(irefin(im)/20.0) ) +1
             irefin(im)  = max( 1 , 19*floor(irefin(im)/20.0) )

         endif
 221     continue

      end if

C
      If (Iatt .eq. 0) Then
        Hord(2) = Hordlrg
        Pmax(2) = Philrg
        Xloc1 = Xxold(Imreg)
        Xloc2 = Xxold(Imreg+1)
        Ichkpt = Nvold/2
        If (Xloc1 .lt. Voldmsh(Ichkpt)) Ichkpt = 1
        Xb = Voldmsh(Ichkpt)
        Do 225 J = Ichkpt, Nvold-1
          Xa = Xb
          Xb = Voldmsh(J+1)
          If (Xloc1 .ge. Xa .and. Xloc1 .lt. Xb) Then
            If (Xloc2-Xb .lt. (Xloc2-Xloc1)/2.d0) Then
              Pmax(1) = Phiold(J)
              Hord(1) = (Xb-Xa)*Pmax(1)
            Else
              Pmax(1) = Phiold(J+1)
              Hord(1) = (Voldmsh(J+2)-Xb)*Pmax(1)
            Endif
            Goto 228
          Endif
 225    Continue
 228    Hsml=(Xloc2-Xloc1)/Dble(Irefin(Imreg))
        Npr(2)=Int(Esum)
      Endif

      If (Iatt .eq. -1) Npr(1) = Int(Esum)
      If (Ifinal .eq. 1 .and. Mshchng .eq. 0) Iatt = 0
      If (Nmest .gt. Nmax) Go To 360

*  It Appears That We Can Perform The Desired Selective Mesh
*  Refinement.

*  Now Begin Running Through The Mesh, Adding And Possibly Deleting
*  Points As Indicated By The Irefin Array.

*  The Integer New Is A Count Of The Number Of Intervals In
*  The Tentative Mesh Being Generated By The Refinement Strategy.

      New = 1

*  The First Interval Is Treated As A Special Case, Since Xx(1)
*  Always Remains In The Mesh, And Cannot Be A Fixed Point.

      Rlen = Xxold(2) - Xx(1)
      Slen = Rlen
      If (Irefin(1) .gt. 1) Then
         Dx = Rlen/Irefin(1)
         Do 230 J = 2, Irefin(1)
            New = New + 1
            Xx(New) = Xx(1) + (J-1)*Dx
            If (Xx(New) .eq. Xx(New-1)) Then
               Iprec = 2
               Return
            Endif
  230    Continue
      Endif

*  The Fixed Points Specified By The Fixpnt Array Cannot Be
*  Removed From The Mesh.  The Value Fxnext Indicates The 'Next'
*  Fixed Point.  When No Further Fixed Points Remain To Be Processed
*  (Or If Nfxpnt = 0), Fxnext Is Set To A Value Strictly Larger Than
*  The Last Mesh Point, So That No Mesh Point Can Equal Fxnext.
*  This Way We Need To Compare Only The Values Of Xxold(I)
*  And Fxnext.

      Ifxcnt = 1
      If (Nfxpnt .eq. 0) Then
         Fxnext = Onep1*Abs(Xxold(Nmsh))
      Else
         Fxnext = Fixpnt(Ifxcnt)
      Endif

*  Jtkout Is A Counter Of The Number Of Consecutive Points That
*  Have Been Removed From The Mesh.

      Jtkout = 0
      Do 330 Im = 2, Ninter
         Rlold = Rlen
         Rlen = Xxold(Im+1) - Xxold(Im)

*  If Xxold(Im) Is The Next Fixed Point, It Cannot Be Removed
*  And So We Don'T Test Its Error Estimates.

         If(Xxold(Im) .eq. Fxnext)  Then

            Ifxcnt = Ifxcnt + 1
            If(Ifxcnt .gt. Nfxpnt) Then
               Fxnext = Onep1*Abs(Xxold(Nmsh))
            Else
               Fxnext = Fixpnt(Ifxcnt)
            Endif

         Elseif (Irefin(Im) .eq. 1) Then

*  If Xxold(Im) Is Not A Fixed Point And Irefin(Im) = 1, We May Wish
*  To Remove Point Im From The Mesh.

*  If We Are Considering Removing Points And Jtkout = 0, This
*  Is The First Point In A Possible Consecutive Set To Be Removed,
*  And We Initialize Phihat, Which Represents A Maximum Of
*  Certain Estimates.
*  If Jtkout Is Not Zero, Previous Points Contiguous To This
*  Point Have Been Removed, And Phihat Retains Its Previous Value.

            Slen = Slen + Rlen

            If (Jtkout .eq. 0) Then
                Ind1 = Ihcomp(Im-1)
                Phihat = Ermeas(Ind1,Im-1)/(Rlold**Ipow)
            Endif
            Phihat = Max(Phihat,
     *                 Ermeas(Ihcomp(Im),Im)/(Rlen**Ipow))
            Val1 = Phihat*(Slen**Ipow)
            If (Val1 .le. Phitst
     *         .and. Jtkout .lt. 4
     *         .and. r4(im) .lt. fatt_r3 ) Then

*  Increment The Counter Of Removed Points.
*  'Remove' The Mesh Point Xxold(Im) By Not Including It.

               Jtkout = Jtkout+1
               Go To 330
            Endif
*        End Of Logic For Irefin(Im) = 1.
         Endif

         Jtkout = 0
         New = New + 1
         Xx(New) = Xxold(Im)
         If (Xx(New) .eq. Xx(New-1)) Then
            Iprec = 2
            Return
         Endif
         If (Irefin(Im) .gt. 1) Then
            Dx = Rlen/Irefin(Im)
            Do 300 J = 2, Irefin(Im)
               New = New + 1
               if (New .ge. Nmax) goto 360
               Xx(New) = Xxold(Im) + (J-1)*Dx
               If (Xx(New) .eq. Xx(New-1)) Then
                  Iprec = 2
                  Return
               Endif
 300        Continue
         Endif
         Slen = Rlen

*  If The New Mesh Contains Too Many Points, Branch Out Of The
*  Loop To Try Alternative Strategies.


         If (New .gt. Nmax-1) Goto 360

  330 Continue

*  To End Up Here, We Have Processed The Entire Interval,
*  And Have Neither Exceeded The Specified Maximum Nor
*  Exceeded Three Times The Number Of Intervals In The Old
*  Mesh.  The Last Mesh Point Remains Unchanged.

      New = New + 1
      Xx(New) = Xxold(Nmsh)
      If (Xx(New) .eq. Xx(New-1)) Then
         Iprec = 2
         Return
      Endif
      Nmsh = New
      Maxmsh = .false.
      Return

  360 Continue

*  To Reach Here, The Number Of Mesh Points Created At Some Stage
*  Of The Refinement Process Was Larger Than The Maximum Permitted
*  Value Nmax.

*  Check Whether The Mesh Can Safely Be Doubled.

      If ((2*Nmsh-1) .lt. Nmax) Then

*  Double The Mesh.

         Call acDblmsh (Nmsh, Nmax, Xx, Nmold, Xxold, Maxmsh)
         dDouble = .true.

*  If The Number Of Intervals Is Too Large And The Mesh Cannot Be
*  Doubled, Increase The Threshold Thres By A Factor Of Erdcid And
*  Try The Selective Refinement Again.
*  If This Happens Three Times Without Success Or If Thres Exceeds
*  Or Is Equal To Errmax, Stop.  (In This Case, We Know Already
*  That Doubling The Mesh Produces Too Many Points.)

      Elseif (Thres .lt. Errmax .and. Ithres .lt. 3) Then
         Ithres = Ithres + 1
         Thres = Erdcid*Thres
         If(Thres .gt. Errmax) Thres = Errmax
         nptcond = nptcond/2
         Call Dcopy(Nmsh, Xxold, 1, Xx, 1)
         Go To 200
      Else
         Nmsh = 2*Nmsh - 1
         Maxmsh = .true.
      Endif
      Return
      End




      SUBROUTINE acmoncondmsh_l(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,
     *                        nptcond,r4,amg,linear)

      implicit double precision (a-h,o-z)

      dimension  xx(*)
      dimension  amg(*), r4(*)

      intrinsic abs
      intrinsic max
      intrinsic int

*  blas: dcopy
*  double precision dlog

      logical pdebug, use_c, comp_c, linear
      common/algprs/ nminit, iprint, idum, use_c, comp_c
      Common/acAlgprs/Maxcon,Itsaim,Uval0
      parameter  ( zero = 0.0d+0, one = 1.0d+0 )

* the function moncond compute the monitor function based on the
* conditioning parameters and the factor used to perform the mesh selection

       do i=1,nmsh-1
         r4(i) = (xx(i+1)-xx(i))*dabs(amg(i+1)- amg(i))
       end do



       r2 = r4(1)
       do i=2,nmsh-1
          r2 = r2 + r4(i)
       end do

      if (linear) then
         cfac=1e-10
      else
         cfac=1d-5
      endif


       do i=1,nmsh-1
c         r4(i) = r4(i)+(xx(i+1)-xx(i))*(r2/(xx(nmsh)-xx(1)))*cfac
          r4(i) = r4(i)+(r2/(xx(nmsh)-xx(1)))*cfac
       end do

       r1 = r4(1)
       do i=2,nmsh-1
          r1 = max(r1,r4(i))
       end do


       do i=1,nmsh-1
         r4(i) = r4(i)/r1
       end do


        r2m = r4(1)
        r1m = r4(1)
       do i=2,nmsh-1
          r1m = min(r1m,r4(i))
          r2m = r2m + r4(i)
       end do

       r3 = r2m/(nmsh-1)
       fatt_r3  = r3*1.0d-3
       fatt_r1r3= max(r3,0.25d0)
        nptm = 0
          nptr = 0

         r1m = one

        do i=1,nmsh-1
          if (r4(i) .ge. fatt_r1r3)  nptm = nptm + 1
          if (r4(i) .le. fatt_r3)  nptr = nptr + 1
          enddo


      if (nptm .LE. 1) then
         nptcond=14
      elseif (nptm .le. 2) then
         nptcond=10
      elseif (nptm .le. 4) then
         nptcond=8
      elseif (nptm .le. 8) then
         nptcond=6
      elseif (nptm .le. int(nmsh/20.0) ) then
         nptcond=4
      else
         nptcond=2
      endif
      end


