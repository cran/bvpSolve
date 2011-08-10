c ===================================================================================
c Routines in common to twpbvpC and twpbvpLC
c ===================================================================================

c Francesca 13/01/2011: modified newteq to work with rshgiv
c  modified external and intrinsic instructions

c karline: got rid of pdebug (common algprs)

c ===================================================================================
c print R-messages
c ===================================================================================

      subroutine rprint(msg)
      character (len=150) msg

            call dblepr(msg, 150, 0, 0)
      end subroutine

c ===================================================================================
c initu resets u after re-meshing for linear problems or for nonlinear problems
c when interpolation of the old solution is not used.
c it interpolates between (Xguess,Yguess), if these are inputted
c otherwise sets to constant value Uval0
c ===================================================================================


      subroutine initu(ncomp,nmsh,xx,nudim,u,
     +                       nugdim,nmguess,xguess,uguess)
      implicit double precision (a-h,o-z)
      dimension xx(*), u(nudim, *), xguess(*), uguess(nugdim,*)

      character(len=150) msg

      logical use_c, comp_c, giv_u
      integer nmguess, ureset
      common/algprs/nminit, iprint, idum, use_c, comp_c
      common/gu/ giv_u, ureset

      uval0 = 0.0d0

      ureset = ureset + 1

   99 format('initu',1pd15.5)

      IF (giv_u) THEN
       if (iprint .ne. -1) then
        write(msg,99) 0.0d0
        call Rprint(msg)
       endif

       call interp(ncomp, nmsh, xx, nudim, u,
     *                  nugdim,nmguess, xguess, uguess)

      ELSE
       if (iprint .ne. -1) then
        write(msg,99) Uval0
        call Rprint(msg)
       endif
        call mtload(ncomp, nmsh, Uval0, nudim, u)
      ENDIF

      return
      end

c ===================================================================================

      SUBROUTINE CONDESTIM(ALEFT,ARIGHT,NMSH,NCOMP,N,XX,TOPBLK,
     *            NRWTOP,NOVRLP,ARRAY,
     *          NRWBLK,NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,ISIGN,OMG,
     *          C1,WORK,KPPA,GAMMA,SIGMA,CKAPPA,CKAPPA2)

C     **************************************************************
C
C     COMPUTES THE FIRST AND LAST BLOCK COLUMN OF INVERSE MATRIX AND
C     ESTIMATE THE CONDITION NUMBERS KPPA, GAMMA
C
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        IMPLICIT NONE
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,WORK,C1,
     *          OMG,GAMMA,gamma1,KPPA,MINMG, GAMMAK,KPPAK,ckappa1,
     *          KPPAI, GAMMAI, CKAPPA,CKAPPA2,kappa1_n, kappa2_n,ckmax
        DOUBLE PRECISION BOMEGA1,BOMEGA2, SIGMA, SIGMAK
        DOUBLE PRECISION ALEFT,ARIGHT,XX,  CSUM,ZNORM,ALTSGN,TEMP
        INTEGER ISIGN(*)
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,IPVCD
        INTEGER NCOMP,NMSH,idmx,idmn,idamax,idomg,job
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),OMG(*),
     *          BOTBLK(NRWBOT,*),WORK(*),C1(ncomp,*),XX(*),
     *          IPVCD(*)
        INTEGER k,i,j,l, kn,INDNM,itmax,its
        LOGICAL FIRSTCOL
        DOUBLE PRECISION dasum
        EXTERNAL DASUM

        INTRINSIC SIGN
        DOUBLE PRECISION SIGN
      character(len=150) msg

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
             if (OMG(l).LT.CSUM) then
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
c estimation of the new value of sigma
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
           END IF
           IF (KPPAK .GT. KPPAI) THEN
              KPPAI = KPPAK
           END IF
       END DO

       CKAPPA = KPPA
       itmax=0

c estimation of the infinity norm using the Hager algorithm
c the initial vector depends on the value of sigma
c algorithm called kappastiffbvp

 450    if (sigma .lt. 10 .and. itmax .eq. 0 ) then
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

c solution of the traspose system

       CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *                 NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,
     *                 WORK,1)

c scaling

             DO i=1,NBLOKS
               DO j=1,NRWBLK
                  WORK( (i-1)*NRWBLK+NRWTOP+j ) =
     *             (XX(i+1)-XX(i))*WORK( (i-1)*NRWBLK+NRWTOP+j)
               ENDDO
             END DO
c computation of kappa
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

         IF ( kappa1_n + kappa2_n .gt. CKAPPA ) THEN
           CKAPPA2 = kappa2_n
           CKAPPA = kappa1_n+kappa2_n
           CKAPPA1 = kappa1_n
         ELSEIF (itmax .ne. 0) THEN
            goto 500
         END IF

c the new step if made only if sigma < 10 and kappa2 > kappa1/10
       if ((SIGMA .lt. 10 .or. ckappa2 .GT. 0.1d0*KPPA)
     *                                   .and. itmax.le.2 ) then
          itmax =itmax+1
c compute the right hand side
          DO ,I = 1,N
            WORK(I) = SIGN(1.0d0,WORK(I))
          END DO

c scaling
            DO i=1,NBLOKS
               DO j=1,NRWBLK
                  WORK( (i-1)*NRWBLK+NRWTOP+j ) =
     *             (XX(i+1)-XX(i))*WORK( (i-1)*NRWBLK+NRWTOP+j)
               ENDDO
            END DO
c solve the linear system
            CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,
     *             WORK,0)


c test for convergence

      ZNORM=0
      DO I=1,N
         IF ( ABS(WORK(I)) .GT. ZNORM) THEN
            ZNORM = ABS(WORK(I))
            indnm = I
         END IF
      END DO

      IF (ZNORM .GT. ABS(work(idmx))  .or. its .eq. 1) THEN
        idmx=indnm
        goto 450
      END IF

      end IF

 500     IF (.NOT. FIRSTCOL) THEN
             KPPA = ckappa1
             GAMMA = GAMMAI
         END IF

         RETURN
       END

c ===================================================================================

        SUBROUTINE ESTIMKAPPA(NMSH,NCOMP,N,XX,TOPBLK,
     *            NRWTOP,NOVRLP,ARRAY, NRWBLK,
     *          NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,ISIGN,C1,WORK,ckappa)

C     **************************************************************
C
C
C     ESTIMATE THE CONDITION NUMBER  kappa
C
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        IMPLICIT NONE
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,WORK,C1
        DOUBLE PRECISION XX
        INTEGER ISIGN(*)
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,IPVCD
        INTEGER NCOMP,NMSH,job
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),WORK(*),XX(*),C1(NCOMP,*),
     *          IPVCD(*)


        double precision ckappa
        INTEGER KASE,ISOLVE
        INTEGER k,i,j,l,jj
      character(len=150) msg


c
c       DO THE CONDITION NUMBER ESTIMATION BY HIGHAM:
C       (DONE IN COLROW) infinity-norm
c

        ISOLVE = 0
        KASE = 0
 55               CALL DONEST(N,C1,WORK,ISIGN,ckappa,KASE)

        IF (KASE .NE. 0) THEN
           ISOLVE = ISOLVE+1

           IF (KASE .eq. 1) THEN
               JOB = 1
           ELSE
               JOB = 0
           END IF

           IF (JOB .eq. 0) THEN
             DO i=1,NBLOKS
               DO j=1,NRWBLK
                  WORK( (i-1)*NRWBLK+NRWTOP+j ) =
     *             (XX(i+1)-XX(i))*WORK( (i-1)*NRWBLK+NRWTOP+j)
               ENDDO
             END DO
           END IF
           CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *                 NCLBLK,NBLOKS,BOTBLK,NRWBOT,IPVCD,
     *                 WORK,JOB)
           IF (JOB .eq. 1) THEN
             DO i=1,NBLOKS
               DO j=1,NRWBLK
                  WORK( (i-1)*NRWBLK+NRWTOP+j ) =
     *             (XX(i+1)-XX(i))*WORK( (i-1)*NRWBLK+NRWTOP+j)
               ENDDO
             END DO
           END IF
           GOTO 55
        END IF


      RETURN
      END

c ===================================================================================

      subroutine conv4( ncomp, nmsh, ntol, ltol, tol, linear, nmax,
     *             xx, nudim, u, defexp, defimp, def, fval,
     *             tmp, bhold, chold, dsq, dexr, usave,
     *             ratdc, rerr, ipivot, nmold, xxold,
     *             smooth, reaft6, onto6, strctr, trst6, ddouble ,
     *             fsub, maxmsh, succes, first4,
     *             amg,stab_cond,ckappa,stiff_cond,r4,
     *             nfxpnt, fixpnt,irefin,rpar,ipar,
     *             nmguess,xguess,nygdim,yguess)

      implicit double precision (a-h,o-z)
      dimension rpar(*), ipar(*)
      dimension ltol(ntol), ipivot(*)
      dimension xx(*), u(nudim,*), tol(ntol), xguess(*),yguess(nygdim,*)
      dimension defexp(ncomp,*), defimp(ncomp,*), def(ncomp,*)
      dimension fval(ncomp,*), tmp(ncomp,4)
      dimension bhold(ncomp,ncomp,*), chold(ncomp, ncomp,*)
      dimension dsq(ncomp,ncomp), dexr(ncomp), usave(ncomp,*)
      dimension xxold(*)
      dimension ratdc(*), rerr(ncomp,*), amg(*)
      dimension fixpnt(*), irefin(*), r4(*)

      logical linear, smooth, reaft6, onto6, strctr, trst6
      logical ddouble , succes, maxmsh, first4
      logical stab_cond, stiff_cond
      external fsub

      logical callrt, oscchk, reposs, savedu

      logical use_c, comp_c
      common/algprs/ nminit, iprint,idum,  use_c, comp_c
      common/mchprs/flmin, flmax, epsmch

      parameter (hundth = 1.0d+5, rerfct = 1.5d+0)
      parameter (power = 1.0d+0/6.0d+0, one = 1.0d+0)
      parameter (zero = 0.0d+0, huge = 1.0d+30)

      intrinsic abs
      intrinsic max

      save  dfold, oldrt1, savedu, reposs

      logical adjrer
      character(len=150) msg

*  The Newton iteration converged for a 4th order solution.

      if (iprint .eq. 1) THEN
        write(msg,901)
        CALL rprint(msg)
      ENDIF


      if (first4) then
         dfold = zero
         oldrt1 = huge
         savedu = .false.
         reposs = .false.
         first4 = .false.
      endif

      succes = .false.
      maxmsh = .false.


*  Compute the explicit deferred correction array defexp.
*  The parameter fval is an ncomp by nmsh array, calculated
*  by a previous call to fneval for this mesh and solution.

      call dfexcl( ncomp, nmsh, xx, nudim, u, defexp, fval,
     *                tmp, fsub, rpar, ipar )

      call matcop(ncomp, ncomp, ncomp, nmsh-1, defexp, def)

*  If the problem has been declared as smooth, or we might wish
*  to perform Richardson extrapolation after trying to calculate a
*  6th order solution, just go on to the 6th order logic (don't
*  bother to calculate implicit deferred corrections).

      if ( smooth .or. reaft6) then
         onto6 = .true.
         return
      endif

*  Compute the cheap implicit deferred correction array defimp.
*  The array chold must be unchanged from a call to jaccal
*  for this mesh and solution.
*  The temporary arrays dsq and dexr are calculated inside dfimcl.

      call dfimcl( ncomp, nmsh, defexp, chold, dsq, dexr,
     *                 ipivot, defimp)


*  Call dccal to calculate: dfexmx, the maximum-magnitude element
*  of defexp in components for which tolerances are specified;
*  incmp, inmsh, and intol, the indices of the component,
*  mesh interval, and tolerance of dfexmx; derivm, the
*  maximum-magnitude element of fval(incmp,*) for all mesh
*  points; dfimmx, the maximum-magnitude element of defimp in
*  component incmp; the ratios rat1 and rat2; and the array
*  ratdc (used in osc).
      dfctol = hundth*epsmch
      call dccal( ncomp, nmsh, ntol, ltol,
     *                  defexp, defimp, dfctol, fval,
     *                  ratdc, dfexmx, incmp, inmsh, intol,
     *                  derivm, dfimmx, rat1, rat2)

*  decid4 sets logical flags that determine the next step of the
*  algorithm.

      call decid4(linear, rat1, rat2, dfexmx, dfimmx,
     *     derivm, dfold, tol(intol), oldrt1,
     *     onto6, smooth, callrt, strctr, oscchk, ddouble , reposs)


      oldrt1 = rat1
      dfold = dfexmx
      if (callrt) then

*  ratcor calculates a more expensive implicit deferred correction.
*  The matrix bhold must be unchanged from the last call to jaccal.
*  If callrt is true, onto6 is always true also.

         call ratcor( ncomp, nmsh, xx, defimp, bhold, def)

      elseif (linear) then
         if (oscchk) then
            call osc (ncomp, nmsh, dfexmx, incmp,
     *            defexp, ratdc, ddouble , inmsh, onto6, trst6,
     *            smooth)

         elseif (reposs)  then

*  If reposs (`Richardson extrapolation possible') is true
*  for two successive 4th order solutions, a special termination
*  test may be used based on Richardson extrapolation.

*  If reposs is true for the first time, savedu will
*  be false; savedu can be true only when reposs is true for a
*  second consecutive mesh (a doubled version of the first).

            if (savedu) then

*  The flag savedu is .true. when the immediately preceding
*  converged 4th order solution was saved in the array usave,
*  and the mesh was then doubled.
*  In this case, the routine rerrvl is called to compute a
*  Richardson extrapolation (RE) error estimate remax.
*  The rerr array does not need to be adjusted, so adjrer is false.

               adjrer = .false.
               call rerrvl( ncomp, nmsh, nudim, u, usave, ntol,
     *              ltol, rerr, remax, itlmx, adjrer )

               if (remax .lt. rerfct*tol(itlmx)) then
                  succes = .true.
                  return
               endif
*           end of logic for savedu = .true.
            endif

*  Here, reposs is .true., but either we hadn't saved the previous
*  solution, or else the RE error estimate is not small.
*  Save u in usave, and set savedu to .true.
*  Set ddouble  to .true. to indicate that the mesh should be doubled.
*  Set .onto6. to .false. to return to the top of the 4th order
*  iteration loop.

      call matcop(nudim, ncomp, ncomp, nmsh, u, usave)
            ddouble= .true.
            savedu = .true.
            onto6 = .false.
*        end of logic for reposs = .true.
         endif
*     end of logic for linear
      endif

c      if (stiff_cond .and. stab_cond) onto6 = .true.

      if (.not. onto6) then

*  NB: onto6 can be false only for linear problems

         if (ddouble ) then
           if ((use_c) ) then
             if (stiff_cond .and. .not. stab_cond) then
               call selcondmsh(ncomp, nmsh,
     *            nfxpnt, fixpnt,  nmax, xx,  irefin,
     *            nmold, xxold, ddouble , maxmsh,r4,amg)
                  ddouble= .false.
              else
                call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
              end if
            else
               call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
           end if
         else



*  Refine the mesh near interval inmsh; the number of points
*  to be added depends on the relative size of dfexmx compared
*  to u and tol.

*  If the conditining parameter are not stabilized and the problem
*  is considered stiff  the points are added and removed using
*  the monitor function given by the conditioning parameters
*  otherwise  the old technique and the conditioning parameters
*  are used together.
*  new functions: selcondmsh, smpselcondmsh
*
             drat = dfexmx/
     *               (max(one, abs(u(incmp,inmsh)))*tol(intol))

            numadd = drat**power

         if ((use_c) .AND. (stiff_cond)) then

            call smpselcondmsh(ncomp, nmsh,
     *        nfxpnt, fixpnt,  nmax, xx,  irefin,inmsh,numadd,
     *        nmold, xxold, ddouble , maxmsh,r4,amg)
         else

             call smpmsh (nmsh, nmax, xx, inmsh, numadd,
     *             nmold, xxold, maxmsh)

         endif

      end if
      if (.not. maxmsh) call initu(ncomp, nmsh, xx, nudim,
     &   u, nygdim,nmguess,xguess, yguess)

*     end of logic for .not. onto6
      endif



      return

  901 format(1h ,'conv4')

      end

c ===================================================================================

      subroutine fail4( ncomp, nmsh, nlbc, ntol, ltol,
     *             xx, nudim, u, rhs, linear, nmax,
     *             nmold, xxold, uold, tmwork,
     *             iorder, iflnwt, itnwt, ddouble, maxmsh,
     *             numbig, nummed,r4,amg,stab_cond,stiff_cond,
     *             ill_cond_newt,nfail4,
     *             nfxpnt, fixpnt, irefin,itcond,itcondmax,rpar,ipar,
     *             nmguess,xguess,nygdim,yguess)


      implicit double precision (a-h,o-z)
      dimension rpar(*),ipar(*)
      dimension ltol(ntol)
      dimension xx(*), u(nudim, *), rhs(*),xguess(*), yguess(nygdim,*)
      dimension xxold(*), uold(ncomp, *), tmwork(*)
      logical linear, ddouble, maxmsh
      logical stab_cond,stiff_cond, ill_cond_newt
      dimension amg(*),r4(*), fixpnt(*), irefin(*)
      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum,  use_c, comp_c
      character(len=150) msg

*  The Newton procedure failed to obtain a 4th order solution.

      if (iprint .eq. 1) THEN
        write(msg,901)
        CALL Rprint(msg)
      ENDIF



      maxmsh = .false.


      if (iflnwt .eq. -1 ) then

*  iflnwt = -1 means that the Jacobian was considered singular.
*  (This is the only possible failure for a linear problem.)

         call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
         call initu(ncomp, nmsh, xx, nudim, u, nygdim,
     +                    nmguess,xguess, yguess)
      else
*  The routine mshref decides how to refine the mesh and then
*  performs the refinement, either by doubling or based on
*  the rhs vector at the best point.


         call mshref(ncomp, nmsh, nlbc, ntol, ltol,
     *             iorder, rhs, tmwork,
     *             nmax, xx, nmold, xxold, ddouble, maxmsh,
     *             numbig, nummed,
     *             amg,stab_cond,stiff_cond,
     *             r4, nfxpnt,fixpnt, irefin,itcond,itcondmax)





         if (.not. maxmsh) then
              if (linear  .or. itnwt .eq.0 .or.nfail4 .ge. 3 )  then
                call initu(ncomp, nmsh,xx,nudim,u,nygdim,
     +                       nmguess,xguess, yguess)
                nfail4 = 0
              else
*  Interpolate the partially converged solution.
               call matcop(nudim, ncomp, ncomp, nmold, u, uold)
               call interp(ncomp, nmsh, xx, nudim, u, ncomp,nmold,
     *                 xxold, uold)
              endif
         endif

*     End of logic for failure because of some reason other than
*     singularity.
      endif

      return
 901   format(1h ,'fail4')
      end

      subroutine conv6(ncomp, nmsh, ntol, ltol, tol,
     *             nudim, u, uold, etest6, err6,
     *             trst6, onto8, reaft6, succes)

      implicit none
      integer ncomp, nmsh, ntol, nudim
      integer ltol(ntol)
      double precision u(nudim,*), tol(ntol)
      double precision uold(ncomp,*), etest6(*),err6
      logical trst6, onto8, reaft6, succes

      integer nminit, iprint, idum
      logical use_c, comp_c
      common/algprs/ nminit, iprint,idum,  use_c, comp_c

      logical errok
      character(len=150) msg

*  The Newton iteration converged for a 6th-order solution.

      if (iprint .eq. 1) THEN
        write(msg,901)
        CALL rprint(msg)
      ENDIF


      succes = .false.

*  The logical flag reaft6 is true only when the 6th order solution
*  failed.  Since the 6th order solution converged, reaft6 is false.

      reaft6 = .false.
      onto8 = .true.

* Calculate the error estimates for the 4th order solution.

      call errest (ncomp, nmsh, ntol, ltol, tol,
     *          nudim, u, uold, etest6, err6, errok)

      if (trst6 .and. errok) then
         succes = .true.
         return
      endif

      return
  901 format(1h ,'conv6')
      end

c ===================================================================================

      subroutine fail6( ncomp, nmsh, nlbc, ntol, ltol, tol,
     *              nfxpnt, fixpnt, iorder, nmax,
     *              xx, nudim, u, rhs, usave, xxold, uold, nmold,
     *              ihcomp, irefin, rerr, ermx, tmwork,
     *              reaft6, ddouble , succes, maxmsh,
     *              numbig, nummed,
     *             r4,amg, stab_cond,ckappa1,gamma1,ckappa,stiff_cond,
     *            itcond, itcondmax)
      implicit double precision (a-h,o-z)
      dimension ltol(ntol)
      dimension fixpnt(*), tol(*), rhs(*)
      dimension xx(*), u(nudim,*), xxold(*)
      dimension uold(ncomp,*), usave(ncomp,*)
      dimension ihcomp(*), irefin(*)
      dimension rerr(ncomp,*), ermx(*), tmwork(*)
      logical reaft6, ddouble , succes, maxmsh
      logical adjrer
      logical stab_cond, stiff_cond
      dimension amg(*), r4(*)

      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c

*  blas: dcopy

      parameter (eight = 8.0d+0)

*  nmpt is the standard number of mesh points to be added to selected
*  intervals when the mesh is altered based on the distribution
*  in size of the rhs vector.
      parameter ( nmpt = 15)
      character(len=150) msg

*  Non-convergence of 6th order.

      if (iprint .eq. 1) THEN
        write(msg,901)
        CALL rprint(msg)
      ENDIF

      succes = .false.
      maxmsh = .false.

*  NB: the problem must be nonlinear.  Linear problems will either
*  fail to converge for 4th order, or else, once they've converged
*  for 4th order, must converge for 6th order.

*  Restore the u array to the previous 4th order solution.

      call matcop(ncomp, nudim, ncomp, nmold, uold, u)

      if (reaft6 .and. iprint .ge. 0) THEN
        write(msg,9999)
        CALL rprint(msg)
      ENDIF
 9999 format(1h ,'in fail6, reaft6is true')
      if (.not.reaft6 .and. iprint .ge. 0) THEN
        write(msg,9998)
        CALL rprint(msg)
      ENDIF
 9998 format(1h ,'in fail6, not reaft6')
      if (ddouble.and. iprint .ge. 0) THEN
        write(msg,9997)
        CALL rprint(msg)
      ENDIF
 9997 format(1h ,'in fail6, ddouble  is true')
      if (.not.ddouble.and. iprint.ge.0 ) THEN
        write(msg,9996)
        CALL rprint(msg)
      ENDIF

 9996 format(1h ,'in fail6, not double')

      if (.not. reaft6 .or. .not. ddouble ) then

*  Here, either
*  (1) the mesh for which this 6th order solution failed
*  is not a doubled version of the immediately preceding mesh, or
*  (2) for the immediately preceding mesh, it is not true
*  that the 4th order converged and the 6th order failed.

*  Setting reaft6 to .true. signals that Richardson extrapolation
*  may be possible if the next 6th order solution fails.  When
*  reaft6 is true, the routine conv4 immediately sets onto6 to true.

         reaft6 = .true.

*  Save the current 4th order solution in usave.

         call matcop(nudim, ncomp, ncomp, nmsh, u, usave)

*  Use the distribution of sizes in the rhs vector to refine the mesh.



         call mshref(ncomp, nmsh, nlbc, ntol, ltol,
     *              iorder, rhs, tmwork,
     *              nmax, xx, nmold, xxold, ddouble , maxmsh,
     *              numbig, nummed,amg,stab_cond,stiff_cond,
     *              r4,nfxpnt, fixpnt,irefin,itcond,itcondmax)
         if (.not. maxmsh) call interp(ncomp, nmsh, xx, nudim, u,ncomp,
     *                        nmold, xxold, uold)

      else

*  Here, reaft6 and ddoubleare both true.  So for two consecutive
*  meshes, the 4th order converged and the 6th order failed,
*  and the second mesh is the ddoubleof the first.

*  Calculate an error estimate from Richardson extrapolation
*  with the current and previous 4th order solutions.
*  (usave is the 4th order solution saved from the previous (halved)
*  mesh.)
*  Set addrer to .true. to signal that the rerr array should
*  be adjusted.

         adjrer = .true.

         call rerrvl( ncomp, nmsh, nudim, u, usave, ntol, ltol,
     *          rerr, remax, itlmx, adjrer )
         if (iprint.eq.1) THEN
           write(msg,9994)
           CALL rprint(msg)
         ENDIF

 9994    format(1h ,'***in fail6')
         if (iprint.eq.1) THEN
           write(msg,9993) remax, eight*tol(itlmx)
           CALL rprint(msg)
         ENDIF

 9993    format(1h ,'remax',1pe14.4,5x,'8*tol',1pe14.4)
         if (remax .lt. eight*tol(itlmx)) then
            succes = .true.
         else

*  Richardson extrapolation did not give sufficient accuracy.
*  Perform selective mesh refinement on the OLD (NB: old!) mesh
*  and the old (saved) solution, using the error estimate from
*  Richardson extrapolation to guide where the mesh points are placed.
cf of the old mesh we take the values using incx = 2
            nmsh = 1 + (nmsh-1)/2
            call dcopy(nmsh, xxold, 2, xx, 1)
            call dcopy(nmsh, amg, 2, tmwork, 1)
            ipow = 4

*  The rerr array is overwritten by selmsh.
        if (use_c) then
          if (.not. stiff_cond ) then
             call selmsh(ncomp, nmsh, ntol, ltol, tol,
     *             nfxpnt, fixpnt, ipow, nmax,
     *             xx, ncomp, usave, rerr, irefin, ihcomp,
     *             nmold, xxold, ermx, ddouble , maxmsh)

          else
*  The rerr array is overwritten by selconderrmsh.
           call selconderrmsh(ncomp, nmsh, ntol, ltol, tol,
     *             nfxpnt, fixpnt, ipow, nmax,
     *             xx, ncomp, usave, rerr, irefin, ihcomp,
     *             nmold, xxold, ermx, ddouble , maxmsh,
     *             r4,tmwork,stab_cond)

          end if
        else
              call selmsh(ncomp, nmsh, ntol, ltol, tol,
     *             nfxpnt, fixpnt, ipow, nmax,
     *             xx, ncomp, usave, rerr, irefin, ihcomp,
     *             nmold, xxold, ermx, ddouble , maxmsh)
        end if

*  If ddoubleis false on exit from selmsh, the call to selmsh has
*  produced a different (non-doubled) mesh.   Interpolate the
*  saved solution (from the old mesh) onto the mesh newly created
*  by selmsh.
*  NB: Because ddouble is false, we won't try Richardson extrapolation
*  if the next 4th order converges and 6th order fails.

            if (.not. maxmsh) then
               if (.not. ddouble ) then
                          call interp(ncomp, nmsh, xx, nudim, u, ncomp,
     *                   nmold, xxold, usave)
                else

*  Selective mesh refinement based on the old mesh simply
*  produced the same mesh that we started with.  So now refine
*  starting with the doubled mesh (from before) and the solution.

                  reaft6 = .true.
*  Save the solution in usave in case we can carry out Richardson
*  extrapolation in the same circumstances next time.
                  call matcop(nudim, ncomp, ncomp, nmsh, u, usave)

*  Use the distribution of sizes in the rhs vector to refine the mesh.

                  call mshref(ncomp, nmsh, nlbc, ntol, ltol,
     *                   iorder, rhs, tmwork,
     *                   nmax, xx, nmold, xxold, ddouble , maxmsh,
     *                   numbig, nummed,amg,stab_cond,stiff_cond,
     *                   r4,nfxpnt, fixpnt,irefin,itcond,itcondmax)

                  if (.not. maxmsh)
     *                call interp(ncomp, nmsh, xx, nudim, u,ncomp,
     *                        nmold, xxold, usave)


*              end of logic for needing to refine (again) based on the
*              current mesh
               endif

*           end of logic for not too many mesh points
            endif

*        end of logic for failure of Richardson extrapolation
*        to produce a converged solution
         endif

*     end of logic for both reaft6 and ddouble being true
      endif

      return
  901 format(1h ,'fail6')
      end

c ===================================================================================

      subroutine conv8( ncomp, nmsh, ntol, ltol, tol,
     *              nfxpnt, fixpnt, linear, nmax,
     *              xx, nudim, u, def, def6, def8, uold,
     *              ihcomp, irefin, ermx, err6,
     *              etest8, strctr,
     *              ddouble , nmold, xxold, maxmsh, succes, first8,
     *              r4, amg,stab_cond,ckappa1,gamma1,ckappa,stiff_cond,
     *    rpar,ipar,nmguess, xguess,nygdim, yguess)

      implicit double precision (a-h,o-z)
      dimension rpar(*), ipar(*)
      dimension ltol(ntol), tol(ntol)
      dimension fixpnt(*)
      dimension etest8(ntol)
      dimension xx(*),u(nudim,*),def(ncomp,*),xguess(*),yguess(nygdim,*)
      dimension def6(ncomp,*), def8(ncomp,*), uold(ncomp,*)
      dimension ihcomp(*), irefin(*)
      dimension ermx(*), xxold(*), amg(*), r4(*)
      logical linear, strctr, ddouble , maxmsh, succes, first8
      logical stab_cond, stiff_cond

      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c

      intrinsic max

      logical errok

      parameter (one = 1.0d+0, fourth = 0.25d+0, quan8 = 0.025d+0)

      parameter ( efact  = 100.0d+0, huge = 1.0d+30 )

*  blas: dload

      save er6old, er8old
      character(len=150) msg

*  The Newton iteration converged for the 8th order solution.


      if (iprint .eq. 1) THEN
        write(msg,901)
        CALL rprint(msg)
      ENDIF


      if (first8) then
         er6old = huge
         er8old = huge
         first8 = .false.
      endif

      if (.not. linear) then
         call dload(ntol, one, etest8, 1)
      else
         do 10 i = 1, ntol
            etest8(i) = one/max(quan8, tol(i)**fourth)
   10    continue
      endif
      succes = .false.
      maxmsh = .false.

*  Check estimated error.  For a nonlinear problem, all components
*  of etest8 (the ratios used in testing the error) are set to one.
*  For a linear problem, the components of etest8 are in general
*  larger than one.  But if strctr is .true. and the number of mesh
*  points decreased, we set the elements of etest8 to one (which
*  makes a stricter test).

      if (linear .and. strctr .and. nmsh .lt. nmold)
     *   call dload(ntol, one, etest8, 1)



      call errest (ncomp, nmsh, ntol, ltol, tol,
     *          nudim, u, uold, etest8, err8, errok)

c      write(*,*) 'etest8', etest8(1),etest8(2), errok
      if (errok)  then
         succes = .true.
         return
      endif

*  At this point, the 8th order solution converged, but did not
*  satisfy the test for termination.

      if (nmsh .lt. nmold. and.
     *         err6 .gt. efact*er6old .and.
     *         err8 .gt. efact*er8old .or.
     *    nmsh .lt. 3*nmold .and. err8 .gt. er8old ) then

*  If the number of mesh points decreased and the errors in the
*  6th and 8th order solutions did not decrease sufficiently compared
*  to the previous mesh, double the mesh and go back to try to
*  calculate a 4th order solution.


         call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
         if (.not. maxmsh) then
            er6old = err6
            er8old = err8
c If the problem is not linear we use the old solution
c instead the the first value
c old code
c            call initu(ncomp, nmsh, xx, nudim, u,rpar,ipar)
c new code
            if (linear) then
               call initu(ncomp,nmsh,xx,nudim,u, nygdim,
     +                  nmguess,xguess, yguess)
            else
c we do not use u but uold
               call matcop(nudim, ncomp, ncomp, nmold, u, uold)
               call interp(ncomp, nmsh, xx, nudim, u,ncomp,
     *                       nmold, xxold, uold)
            endif
         endif
         return
      endif

*   Here, we know that
*   (1) the number of mesh points exceeds that for the previous mesh; or
*   (2) the number of mesh points decreased, the 6th and 8th order
*       errors did not satisfy the termination tests, but they did not
*       increase so much that the mesh needed to be doubled.

      er6old = err6
      er8old = err8
      if (err8 .le. err6) then


*  Perform selective mesh refinement based on the 8th order deferred
*  corrections.  The value of ipow indicates that the error estimate
*  is of order 6.  Then, for a nonlinear problem, interpolate the
*  latest solution onto the new mesh.

         ipow = 6

*  NB: The array def8 will be overwritten by selmsh.


        if (use_c) then
          if (.not. stiff_cond ) then
            call selmsh(ncomp, nmsh, ntol, ltol, tol,
     *            nfxpnt, fixpnt, ipow, nmax,
     *            xx, nudim, u, def8, irefin, ihcomp,
     *            nmold, xxold, ermx, ddouble , maxmsh)
          else


*  NB: The array def8 will be overwritten by selconderrmsh.
            call selconderrmsh(ncomp, nmsh, ntol, ltol, tol,
     *            nfxpnt, fixpnt, ipow, nmax,
     *            xx, nudim, u, def8, irefin, ihcomp,
     *            nmold, xxold, ermx, ddouble , maxmsh,r4,amg,stab_cond)
         end if
        else
            call selmsh(ncomp, nmsh, ntol, ltol, tol,
     *            nfxpnt, fixpnt, ipow, nmax,
     *            xx, nudim, u, def8, irefin, ihcomp,
     *            nmold, xxold, ermx, ddouble , maxmsh)
        end if



         if (.not. maxmsh) then
            if (linear) then
               call initu(ncomp, nmsh, xx,nudim,u,nygdim,
     *                  nmguess,xguess,yguess)
            else
               call matcop(nudim, ncomp, ncomp, nmold, u, uold)
               call interp(ncomp, nmsh, xx, nudim, u,ncomp,
     *                       nmold, xxold, uold)
            endif
         endif

      else

*  err8 is greater than err6

*  For a linear problem, set all elements of etest8 to one,
*  which makes the error test stricter.  (The elements of etest8
*  may have already been set to one earlier in this routine.)

         if (linear) call dload(ntol, one, etest8, 1)

*  Selectively refine the mesh using the old solution and the
*  6th order deferred correction.  Then, for a nonlinear prpblem,
*  interpolate the old solution onto the new mesh.

         ipow = 4


*  The array def6 will be overwritten by selmsh.
      if (use_c) then
         if (.not. stiff_cond ) then
          call selmsh(ncomp, nmsh, ntol, ltol, tol,
     *            nfxpnt, fixpnt, ipow, nmax,
     *            xx, ncomp, uold, def6, irefin, ihcomp,
     *            nmold, xxold, ermx, ddouble , maxmsh)
           else

            call selconderrmsh(ncomp, nmsh, ntol, ltol, tol,
     *            nfxpnt, fixpnt, ipow, nmax,
     *            xx, ncomp, uold, def6, irefin, ihcomp,
     *        nmold, xxold, ermx, ddouble , maxmsh,r4,amg,stab_cond)

         end if
       else
          call selmsh(ncomp, nmsh, ntol, ltol, tol,
     *            nfxpnt, fixpnt, ipow, nmax,
     *            xx, ncomp, uold, def6, irefin, ihcomp,
     *            nmold, xxold, ermx, ddouble , maxmsh)
        end if
         if (.not. maxmsh) then
            if (linear) then
               call initu(ncomp,nmsh,xx,nudim,u,nygdim,
     +               nmguess,xguess,yguess)
            else
c               call matcop(nudim, ncomp, ncomp, nmold, u, uold)
               call interp(ncomp, nmsh, xx, nudim, u,ncomp,
     *                       nmold, xxold, uold)
            endif
         endif
*     end of logic for err8 greater than err6
      endif

      return

  901 format(1h ,'conv8')
      end

c ===================================================================================

      subroutine fail8(ncomp, nmsh, nfxpnt, fixpnt, nmax,
     *      ntol, ltol, tol, nmold,
     *      xx, nudim, u, def6, xxold, uold,
     *      ihcomp, irefin, ermx, ddouble , maxmsh,
     *       r4,amg, stiff_cond, stab_cond)

      implicit double precision(a-h,o-z)
      dimension fixpnt(*), ltol(ntol), tol(ntol)
      dimension xx(*), u(nudim,*), def6(ncomp,*)
      dimension xxold(*), uold(ncomp,*)
      dimension ihcomp(*), irefin(*)
      dimension ermx(*)
      dimension amg(*), r4(*)
      logical ddouble , maxmsh
      logical stiff_cond, stab_cond

      logical use_c, comp_c
      common/algprs/ nminit, iprint,idum,use_c, comp_c
      character(len=150) msg

*  8th order solution did not converge (the problem must be nonlinear)

      ipow = 4

*  Selectively refine the mesh based on the 6th order deferred
*  correction and the old solution.

*  The def6 array is overwritten by selmsh.

        if ( use_c .and. stiff_cond) then


         call selconderrmsh(ncomp, nmsh, ntol, ltol, tol,
     *            nfxpnt, fixpnt, ipow, nmax,
     *            xx, ncomp, uold, def6, irefin, ihcomp,
     *      nmold, xxold, ermx, ddouble , maxmsh,r4,amg,stab_cond)
        else
         call selmsh(ncomp, nmsh, ntol, ltol, tol,
     *        nfxpnt, fixpnt, ipow, nmax,
     *        xx, ncomp, uold, def6, irefin, ihcomp,
     *        nmold, xxold, ermx, ddouble , maxmsh)
        end if

*  Interpolate to obtain the new initial solution.

      if (.not. maxmsh) then
         call interp(ncomp,nmsh,xx,nudim,u,ncomp, nmold, xxold, uold)
      endif

      return
      end

c ===================================================================================

      subroutine dccal( ncomp, nmsh, ntol, ltol,
     *                     defexp, defimp, dfctol, fval,
     *                     ratdc, dfexmx, incmp, inmsh, intol,
     *                     derivm, dfimmx, rat1, rat2)

      implicit double precision  (a-h,o-z)

      dimension  ltol(ntol)
      dimension  defexp(ncomp,nmsh-1), defimp(ncomp,nmsh-1)
      dimension  fval(ncomp,nmsh), ratdc(nmsh-1)

      logical use_c, comp_c
      common/algprs/ nminit, iprint,idum,  use_c, comp_c

      parameter ( zero = 0.0d+0, one = 1.0d+0 )
      parameter ( rtst = 50.0d+0, tstrat = 0.1d+0 )

      intrinsic abs
      intrinsic max
      character(len=150) msg

*  blas: idamax

*  Find dfexmx, the maximum-magnitude element of defexp
*  in components for which a tolerance is specified.
*  The component index of dfexmx (incmp), its mesh
*  interval index (inmsh), and its tolerance index (intol),
*  are output parameters.

      dfexmx = zero
      do 10 it = 1, ntol
         icmp = ltol(it)
         idmx = idamax(nmsh-1, defexp(icmp, 1), ncomp)
         dval = abs(defexp(icmp, idmx))
         if (dval .ge. dfexmx) then
            dfexmx = dval
            incmp = icmp
            inmsh = idmx
            intol = it
         endif
   10 continue

*  Find derivm (maximum-magnitude element of fval(incmp,*))
*  for all mesh points.

      idmx = idamax(nmsh, fval(incmp, 1), ncomp)
      derivm = abs(fval(incmp, idmx))

*  For component incmp, go through the mesh intervals to calculate
*  (1) dfimmx, the maximum implicit deferred correction;
*  (2) two crucial ratios, rat1 and rat2, used in deciding whether
*      to refine the mesh;
*  (3) the array ratdc of deferred-correction ratios (explicit to
*      implicit).

*  In defining rat1 and rat2, we consider only intervals for
*  which the explicit deferred correction (defexp) exceeds the
*  tolerance dfctol in magnitude.  If it does not, the associated
*  interval does not affect rat1 or rat2, and the value of ratdc
*  is taken as 1.

*  rat2 is the maximum-magnitude ratio of sufficiently large
*  defexp to the larger of defimp and dfctol.

*  rat1 is the maximum-magnitude ratio of sufficiently large
*  defexp to the larger in magnitude of defimp and dfctol, but only
*  for those values of defexp greater than tstrat*dfexmx.
*  Thus by construction rat1 is less than or equal to rat2.

      rat1 = zero
      rat2 = zero
      dfimmx = zero
      smtest = tstrat*dfexmx

      do 100 im = 1, nmsh-1
         texp = defexp(incmp, im)
         timp = defimp(incmp, im)
         dfimmx = max(dfimmx, abs(timp))
         abtexp = abs(texp)
         if (abtexp .le. dfctol) then
            ratdc(im) = one
         else
            if (abs(timp) .lt. dfctol) timp = dfctol
            ratdc(im) = texp/timp
            abrat = abs(ratdc(im))
            rat2 = max(rat2, abrat)
            if (abs(texp) .ge. smtest
     *                .and. abrat .ge. rat1)  rat1 = abrat
         endif
  100 continue

      return

      end

c ===================================================================================

      subroutine decid4(linear, rat1, rat2, dfexmx, dfimmx,
     *     derivm, dfold, tolval, oldrt1,
     *     onto6, smooth, callrt, strctr, oscchk, ddouble , reposs)


      implicit double precision (a-h,o-z)

      logical linear, onto6, smooth, callrt, strctr,
     *     oscchk, ddouble , reposs

      logical use_c, comp_c
      common/algprs/ nminit, iprint,idum,  use_c, comp_c

      logical stest

      parameter ( tenth = 0.1d+0, one = 1.0d+0, two = 2.0d+0 )
      parameter ( thrtwo = 32.0d+0 )
      parameter ( rtst = 50.0d+0, derval = 50.0d+0 )
      character(len=150) msg

*  decid4 evaluates information about the deferred corrections
*  and the nature of the problem, and sets various logical
*  flags that control the subsequent logic of the algorithm.

*  This code has been written for clarity, and NOT for efficiency.

      onto6 = .true.
      callrt = .false.
      smooth = .false.
      oscchk = .false.
      strctr = .false.
      reposs = .false.
      ddouble = .false.

*  rat2 is always greater than or equal to rat1.


      stest = .true.
      if (linear) stest = dfexmx .lt. tenth*dfold



      if (rat2 .lt. rtst) then
         if (stest) then
            smooth = .true.
         else
            oscchk = .true.
         endif

         return
      endif

*  We know now that rat2 .ge. rtst.

      thttol = thrtwo*tolval

      if (rat1 .lt. rtst .and. dfexmx .lt. thttol) then
          if (stest) then
             smooth = .true.
          else
             oscchk = .true.
          endif

          return
      endif

      if (rat1 .lt. rtst .and. dfexmx .ge. thttol) then
         callrt = .true.

         return
      endif

*  We know now that rat1 .ge. rtst (and that rat2 .ge. rtst).

      if (derivm .gt. derval .and. dfexmx .lt. thttol) then
         if (stest) then
            smooth = .true.
         else
            oscchk = .true.
         endif

         return
      endif

      if (derivm .gt. derval .and. dfexmx .gt. thttol) then
         if (dfimmx .lt. one) then
            callrt = .true.

         else
            strctr = .true.
            if (linear) then
               onto6 = .false.
               if (two*rat1 .ge. oldrt1) ddouble = .true.

*           end of logic for linear
            endif
*        end of logic for dfimmx .ge. one
         endif

         return
*     end of logic for derivm .gt. derval .and dfexmx .gt. thttol
      endif

*  To reach this point in the code, both of the following must hold:
*    rat1 .ge. rtst (which means that rat2 .ge. rtst)
*    derivm .le. derval

*  On linear problems, a special termination rule is tried if two
*  conditions hold:
*    (1) the 4th order solution has been computed on two consecutive
*        meshes, one of which is the double of the other
*        (this holds when ddouble  is .true.), and
*    (2) on both meshes, rat1 .ge. rtst, and derivm is small.  When
*        the conditions in (2) hold for a particular mesh, decid4
*        sets reposs to .true. to indicate that Richardson
*        extrapolation may be possible.
*  This set of tests is to take care of transients kept out by
*  initial conditions.

      if (linear) reposs = .true.
      return

      end

c ===================================================================================

      subroutine dfexcl (ncomp, nmsh, xx, nudim, u, defexp, fval,
     *               tmp, fsub, rpar, ipar)

      implicit double precision (a-h, o-z)
      integer ncomp, nmsh
      dimension rpar(*), ipar(*)
      dimension xx(nmsh), u(nudim,nmsh), fval(ncomp, nmsh)
      dimension  defexp(ncomp,nmsh-1), tmp(ncomp,4)
      external fsub

      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c

      parameter ( half = 0.5d+0, fourth = 0.25d+0, thfrth= 0.75d+0 )

      common /consts/ alp1, alp2, alp3, bet0, bet2, bet3, bet4,
     *    a1, b1, c1, d1, e2, f2, c2, d2, b2, a2,
     *    p3, q3, e3, f3, c3, d3, a3, b3, a4, p4, x4, e4, c4,
     *    a5, b5, c5, d5, e5, f5, a6, b6, c6
      integer nfunc, njac, nstep, nbound, njacbound
      common/diagnost/nfunc, njac, nstep, nbound, njacbound

      character(len=150) msg

*  Given the nmsh mesh points xx, the estimated solution
*  u and the array fval of function values at (xx(im), u(*,im)),
*  im = 1,...,nmsh, dfexcl calculates sixth-order explicit
*  deferred correction, stored in the array defexp, indexed
*  over the components and mesh intervals.

*  The array tmp is workspace for 4 intermediate vectors of
*  dimension ncomp.

      do 50 im = 1, nmsh-1
         hmsh = xx(im+1) - xx(im)
         do 10 ic = 1, ncomp
            tmp(ic,1) = (a5*u(ic, im+1) + b5*u(ic, im))
     *         + hmsh*(c5*fval(ic,im)- d5*fval(ic,im+1))
            tmp(ic,2) = (b5*u(ic,im+1) + a5*u(ic,im))
     *         + hmsh*(-c5*fval(ic,im+1) + d5*fval(ic,im))
   10    continue

         call fsub (ncomp, xx(im) + fourth*hmsh, tmp(1,1),
     *          tmp(1,3), rpar, ipar)
         call fsub (ncomp, xx(im) + thfrth*hmsh, tmp(1,2),
     *          tmp(1,4), rpar, ipar)

         do 20 ic = 1, ncomp
            tmp(ic,1) = half*(u(ic,im+1) + u(ic,im))
     *          + e5*hmsh*(fval(ic,im+1) - fval(ic,im))
     *          - f5*hmsh*(tmp(ic,4) - tmp(ic,3))
   20    continue

         call fsub(ncomp, half*(xx(im) + xx(im+1)), tmp(1,1),
     *          tmp(1,2), rpar, ipar)
         do 30 ic = 1, ncomp
            defexp(ic,im) = hmsh*(a6*(fval(ic,im+1) + fval(ic,im))
     *           + b6*(tmp(ic,3) + tmp(ic,4)) + c6*tmp(ic,2))
     *           - u(ic,im+1) + u(ic,im)
   30    continue

   50 continue
      nfunc = nfunc + (nmsh-1)*3
      nstep = nstep + 1
      return
      end

c ===================================================================================

      subroutine df8cal (ncomp, nmsh, xx, nudim, u, fval, def8,
     *      tmp, fsub, rpar, ipar)

*   Given the mesh points xx, the solution u, and the function
*   values fval, df8cal computes eighth-order deferred corrections,
*   which are stored in def8.
*   The array tmp is workspace for 8 intermediate vectors.

      implicit double precision (a-h, o-z)
      integer ncomp, nmsh
      dimension rpar(*),ipar(*)
      dimension xx(nmsh), u(nudim,nmsh), fval(ncomp,nmsh)
      dimension def8(ncomp, nmsh-1),  tmp(ncomp,8)
      external fsub

      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c

      parameter ( half = 0.5d+0, two = 2.0d+0 )
      parameter ( fc1 = 0.625d+0, fc2= 0.375d+0 )

      common /consts/ alp1, alp2, alp3, bet0, bet2, bet3, bet4,
     *    a1, b1, c1, d1, e2, f2, c2, d2, b2, a2,
     *    p3, q3, e3, f3, c3, d3, a3, b3, a4, p4, x4, e4, c4,
     *    a5, b5, c5, d5, e5, f5, a6, b6, c6
      integer nfunc, njac, nstep, nbound, njacbound
      common/diagnost/nfunc, njac, nstep, nbound, njacbound

      character(len=150) msg

      do 100 im = 1, nmsh-1

         hmsh = xx(im+1) - xx(im)

         do 10 ic = 1, ncomp
            tmp(ic,1) = a1*u(ic,im+1) + b1*u(ic,im)
     *         + hmsh*(c1*fval(ic,im+1) + d1*fval(ic,im))
            tmp(ic,2) = b1*u(ic,im+1) + a1*u(ic,im)
     *         - hmsh*(c1*fval(ic,im) + d1*fval(ic,im+1))
   10    continue

         call fsub(ncomp, xx(im) + fc1*hmsh, tmp(1,1), tmp(1,3),
     *    rpar, ipar)
         call fsub(ncomp, xx(im) + fc2*hmsh, tmp(1,2), tmp(1,4),
     *    rpar, ipar)
         do 20 ic = 1, ncomp
            tmp(ic,1) = a2*u(ic,im+1) + b2*u(ic,im)
     *         + hmsh*(c2*fval(ic,im+1) + d2*fval(ic,im)
     *         + e2*tmp(ic,3) + f2*tmp(ic,4))
            tmp(ic,2) = b2*u(ic,im+1) + a2*u(ic,im)
     *         - hmsh*(d2*fval(ic,im+1) + c2*fval(ic,im)
     *         + f2*tmp(ic,3) + e2*tmp(ic,4))
   20    continue

         call fsub(ncomp, xx(im) + (half + alp2)*hmsh, tmp(1,1),
     *               tmp(1,5), rpar, ipar)
         call fsub(ncomp, xx(im) + (half - alp2)*hmsh, tmp(1,2),
     *               tmp(1,6), rpar, ipar)
         do 30 ic = 1, ncomp
            tmp(ic,1) = a3*u(ic,im+1) + b3*u(ic,im)
     *         + hmsh*(c3*fval(ic,im+1) + d3*fval(ic,im)
     *         + e3*tmp(ic,3) + f3*tmp(ic,4)
     *         + p3*tmp(ic,5) + q3*tmp(ic,6))
            tmp(ic,2) = b3*u(ic,im+1) + a3*u(ic,im)
     *         - hmsh*(d3*fval(ic,im+1) + c3*fval(ic,im)
     *         + f3*tmp(ic,3) + e3*tmp(ic,4)
     *         + q3*tmp(ic,5) + p3*tmp(ic,6))
   30    continue

         call fsub (ncomp, xx(im) + (half + alp3)*hmsh, tmp(1,1),
     *                 tmp(1,7), rpar, ipar)
         call fsub (ncomp, xx(im) + (half - alp3)*hmsh, tmp(1,2),
     *                 tmp(1,8), rpar, ipar)
         do 40 ic = 1, ncomp
            tmp(ic,1) = a4*(u(ic,im+1) + u(ic,im))
     *         + hmsh*(c4*(fval(ic,im+1) - fval(ic,im))
     *         + e4*(tmp(ic,3) - tmp(ic,4))
     *         + x4*(tmp(ic,7) - tmp(ic,8)))
   40    continue

         call fsub (ncomp, xx(im) + half*hmsh, tmp(1,1), tmp(1,2),
     *       rpar, ipar)
         do 50 ic = 1, ncomp
            def8(ic,im) =
     *          hmsh*(bet0*(fval(ic,im) + fval(ic,im+1))
     *          + bet2*(tmp(ic,5) + tmp(ic,6))
     *          + bet3*(tmp(ic,7) + tmp(ic,8))
     *          + two*bet4*tmp(ic,2))
     *          - u(ic,im+1) + u(ic,im)
   50    continue

  100 continue
      nstep = nstep + 1
      nfunc = nfunc + (nmsh-1)*7
      return
      end

c ===================================================================================

      subroutine dfimcl( ncomp, nmsh, defexp, chold, dsq, dexr,
     *                 ipivot, defimp)
      implicit double precision (a-h,o-z)
      dimension defexp(ncomp, nmsh-1), chold(ncomp, ncomp, nmsh-1)
      dimension dsq(ncomp, ncomp), dexr(ncomp)
      dimension ipivot(ncomp), defimp(ncomp, nmsh-1)

      logical use_c, comp_c
      common/algprs/ nminit, iprint,idum, use_c, comp_c

*  blas: dcopy

      parameter (zero = 0.0d+0)
      character(len=150) msg

*  dfimcl calculates the rational deferred correction array,
*  which is indexed over the components and mesh intervals.

      call mtload(ncomp, nmsh-1, zero, ncomp, defimp)

      do 100 im = 1, nmsh-1

         call dcopy(ncomp, defexp(1,im), 1, dexr(1), 1 )
         do 50 ic = 1, ncomp
            call dcopy(ncomp, chold(1,ic,im), 1, dsq(1,ic), 1)
   50    continue
         call lufac (ncomp, ncomp, dsq, ipivot, ierlu)

         if (ierlu .eq. 0) then

            call lusol (ncomp, ncomp, dsq, ipivot, dexr,
     *          defimp(1, im))
         endif
  100 continue
      return
      end

c ===================================================================================

      subroutine osc (ncomp, nmsh, dfexmx, incmp,
     *     defcor, ratdc, ddouble , inmsh, onto6, trst6, smooth)

      implicit double precision (a-h, o-z)
      dimension defcor(ncomp, *), ratdc(*)

      logical ddouble , onto6, trst6, smooth

      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c

      intrinsic abs

      parameter ( zero = 0.0d+0, half = 0.5d+0 )
      parameter ( frac1 = 0.1d+0, frac2 = 1.0d-2 )
      character(len=150) msg


*  For linear problems, subroutine osc performs heuristic tests
*  to detect an oscillating solution.  the tests check whether
*  (significant) explicit and implicit deferred corrections have
*  different (componentwise) signs.

*  dfexmx is the maximum-magnitude explicit deferred correction,
*  and is known to occur in component incmp.

*  The array defcor contains the explicit deferred corrections.
*  The array ratdc contains the ratios of explicit to implicit
*  deferred corrections.

*  jsndif counts the number of differences in sign.

      jsndif = 0
      rmax = zero

*  allsum is the sum of the magnitudes of all deferred corrections,
*  smlsum is the sum of the magnitudes of the small deferred
*  corrections, and bigsum is the sum of the magnitudes of the
*  large deferred corrections.  Here, small is defined as less
*  than half of the maximum.

      ninter = nmsh - 1

      allsum = zero
      smlsum = zero
      bigsum = zero
      ibig = 0
      ism = 0

      do 30 im = 1, ninter
         abdef = abs(defcor(incmp,im))
         allsum = allsum + abdef
         if (abdef .lt. half*dfexmx) then
              ism = ism + 1
              smlsum = smlsum + abdef
         else
              ibig = ibig + 1
              bigsum = bigsum + abdef
         endif

*  The counter of sign differences is incremented if (1) ratdc is negative
*  (which means that the two deferred corrections have opposite
*  sign) and (2) the explicit deferred correction is not too small
*  relative to the maximum.

         if (ratdc(im).lt.zero .and. abdef.ge.frac2*dfexmx) then
            jsndif = jsndif + 1

*  If more than 4 sign differences have occurred, exit after setting
*  ddouble to .true., which signals that the mesh
*  should be doubled (i.e., twice as many intervals).

            if (jsndif.gt.4) then
               onto6 = .false.
               ddouble = .true.
               return
            endif
            if (abs(ratdc(im)).ge.rmax) then
               rmax = abs(ratdc(im))
               inmsh = im
            endif
         endif
   30 continue

      avsm = zero
      if (ism.gt.0) avsm = smlsum/ism
      avbg = zero
      if (ibig.gt.0) avbg = bigsum/ibig
      ave = allsum/ninter

      if (avsm.gt.frac1*avbg .or. ave.gt.half*avbg) then

*  The error appears to be uniformly large.
*  Signal that the 6th order solution should be calculated.
         onto6 = .true.

      elseif (jsndif.eq.0) then
*  If there were no sign changes, the problem appears to be smooth.
         smooth = .true.
         onto6 = .true.
      else

*  If the sign changed at between 1 and 4 points, don't go on to
*  6th order, and don't ever accept a 6th order solution even if the
*  error estimate at a later stage indicates that it is OK to do so.
*  Set ddouble to .false., to signal that the mesh will not necessarily
*  be doubled.

          ddouble = .false.
          onto6 = .false.
          trst6 = .false.
      endif
      return
      end

c ===================================================================================

      subroutine ratcor ( ncomp, nmsh, xx, defimp, bhold, dfrat)
      implicit double precision (a-h, o-z)
      dimension xx(nmsh), defimp(ncomp,nmsh-1)
      dimension dfrat(ncomp,nmsh-1), bhold(ncomp,ncomp,nmsh-1)

      character(len=150) msg
*  blas: ddot, dscal

      parameter (zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0)

      ninter = nmsh - 1
      do 10 im = 1, ninter
         hmsh = xx(im+1) - xx(im)
         call dscal(ncomp*ncomp, (-half*hmsh), bhold(1,1,im), 1)
   10 continue

      do 20 im = 1, ninter
      do 20 ic = 1, ncomp
         bhold(ic,ic,im) = bhold(ic,ic,im) + one
   20 continue

      do 30 im = 1, ninter
      do 30 ic = 1, ncomp
         dfrat(ic,im) = ddot(ncomp, bhold(ic,1,im), ncomp,
     *                        defimp(1,im), 1)
   30 continue
      return
      end

c ===================================================================================

      subroutine stcons

      implicit double precision  (a-h,o-z)

*  stcons computes constants needed in integration formulae
*  and stores them in a labeled common area.

      parameter ( one = 1.0d+0,  four = 4.0d+0, two = 2.0d+0 )
      parameter ( five = 5.0d+0, three = 3.0d+0 )
      parameter ( half = 0.5d+0, fourth = 0.25d+0 )

      common /consts/ alp1, alp2, alp3, bet0, bet2, bet3, bet4,
     *    a1, b1, c1, d1, e2, f2, c2, d2, b2, a2,
     *    p3, q3, e3, f3, c3, d3, a3, b3, a4, p4, x4, e4, c4,
     *    a5, b5, c5, d5, e5, f5, a6, b6, c6



      alp1 = one/8.0d+0
      alp1sq = one/64.0d+0
      alp2sq = five/28.0d+0
      alp3sq = five/84.0d+0
      alp2 = sqrt(alp2sq)
      alp3 = sqrt(alp3sq)
      bet0 = one/6.0d+0 - (10.0d+0 - 28.0d+0*(alp2sq + alp3sq))/
     *         (105.0d+0*(one - four*alp2sq)*
     *         (one - four*alp3sq))
      bet2 = -(28.0d+0*alp3sq - three)/
     *         (1680.0d+0*alp2sq*
     *         (one - four*alp2sq)*(alp2sq - alp3sq))
      bet3 = (28.0d+0*alp2sq - three)/
     *          (1680.0d+0*alp3sq*
     *          (one - four*alp3sq)*(alp2sq - alp3sq))
      bet4 = half - bet0 - bet2-bet3
      a1 = half*(one - alp1)*(two*alp1 + one)*
     *          (two*alp1 + one)
      b1 = half*(one + alp1)*
     *          (two*alp1 - one)*(two*alp1 - one)
      c1 = half*(alp1sq - fourth)*(two*alp1 + one)
      d1 = half*(alp1sq - fourth)*(two*alp1 - one)
      uu = alp2*((four*alp2sq - one)**2)/
     *       ((four*alp1sq - one)*(20.0d+0*alp1*alp1 - one))
      vv = ((four*alp2sq - one)**2)/
     *      (16.0d+0*alp1*(four*alp1sq - one))
      e2 = half*(uu + vv)
      f2 = half*(uu - vv)
      rr = half*(alp2*(four*alp2sq - one) +
     *      (one - 12.0d+0*alp1sq)*(e2 + f2))
      ss = fourth*(four*alp2sq - one) - two*alp1*(e2 - f2)
      c2 = half*(rr + ss)
      d2 = half*(rr - ss)
      ww = two*(alp2 - (c2 + d2 + e2 + f2))
      b2 = half*(one - ww)
      a2 = half*(one + ww)
      z1 = (three - 28.0d+0*alp3sq)/
     *         (1680.0d+0*alp2*(four*alp2sq - one)*
     *         (alp2*alp2 - alp3sq)*bet3)
      z2 = one/(105.0d+0*alp3*bet3*
     *    (20.0d+0*alp2sq - one)*(four*alp2sq - one))
      p3 = half*(z1 + z2)
      q3 = half*(z2 - z1)
      u1 = (alp3*((four*alp3sq - one)**2)-
     *       (p3 + q3)*(20.0d+0*alp2sq - one)
     *       *(four*alp2sq - one))/
     *        ((four*alp1sq - one)*(20.0d+0*alp1sq - one))
      v1 = (alp3sq*(one - two*alp3sq)-
     *         two*alp2*(one - four*alp2sq)*(p3 - q3)
     *        -one/8.0d+0)/(two*alp1*(one - four*alp1sq))
      e3 = half*(u1 + v1)
      f3 = half*(u1 - v1)
      r1 = half*(alp3*(four*alp3sq - one) +
     *         (e3 + f3)*(one - 12.0d+0*alp1sq) +
     *         (p3 + q3)*(one - 12.0d+0*alp2sq))
      s1 = alp3sq - fourth - two*alp1*(e3 - f3)
     *         - two*alp2*(p3 - q3)
      c3 = half*(r1 + s1)
      d3 = half*(r1 - s1)
      w1 = two*(alp3 - (c3 + d3 + e3 + f3 + p3 + q3))
      a3 = half*(one + w1)
      b3 = half*(one - w1)
      a4 = half
      p4 = 0.0d+0
      x4 = (three - 28.0d+0*alp2sq)/
     *        (3360.0d+0*alp3*bet4*(four*alp3sq - one)
     *        *(alp3sq - alp2sq))
      e4 = (0.125d+0 + four*alp2*p4*
     *         (one - four*alp2sq) +
     *         four*alp3*x4*(one - four*alp3sq))/
     *         (four*alp1*(four*alp1sq - one))
      c4 = -(0.125d+0 + two*alp1*e4 + two*alp2*p4 +
     *         two*alp3*x4)
      a5 = five/32.0d+0
      b5 = 27.0d+0/32.0d+0
      c5 = 9.0d+0/64.0d+0
      d5 = three/64.0d+0
      e5 = five/24.0d+0
      f5 = two/three
      a6 = 7.0d+0/90.0d+0
      b6 = 16.0d+0/45.0d+0
      c6 = two/15.0d+0
      return
      end

c ===================================================================================

      subroutine fixjac(ncomp, nmsh, nlbc,
     *    iorder, ntol, ltol, tol,
     *    xx, nudim, u, defcor, defnew, delu,
     *    rhs, fval, utrial, rhstri,
     *    rnsq, uint, ftmp, tmprhs,
     *    ajac, topblk, botblk, ipivot,
     *    fsub, gsub, iflag, rpar, ipar)

* Fixed Jacobian iterations.

      implicit double precision (a-h,o-z)
      dimension rpar(*), ipar(*)
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

      external fsub
      external gsub

      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c
      common/mchprs/flmin, flmax, epsmch

      intrinsic  abs
      intrinsic max

*  blas: dcopy, dssq

      parameter ( one    = 1.0d+0 )
      parameter ( xlarge = 1.0d+6, huge = 1.0d+30, lmtfrz = 8)
      parameter ( rngrow = 16.0d+0, rfact = 100.0d+0 )
      parameter ( tolfct = 0.1d+0 )
      character(len=150) msg

*  The iteration scheme uses a fixed Jacobian matrix to solve for
*  correction vectors, once there has been convergence of the Newton
*  iterations on this mesh.   It is assumed that the LU
*  factors of the Jacobian have been left unaltered since
*  their calculation.

      if (iprint .eq. 1) THEN
        write(msg,901)
        CALL rprint(msg)
      ENDIF

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
   10 continue
      ind = ninter*nmsh + nlbc + 1
      call dcopy(ncomp-nlbc, rhs, 1, rhstri, 1)

      call dssq  ( nmsh*ncomp, rhstri, 1, scale, sumsq )
      rnsq = (scale**2)*sumsq

      iter = 0

*  If the initial right-hand side is too large, do not even attempt to
*  solve the nonlinear equations.

      if (rnsq.gt.huge .or.
     *      (iorder.eq. 8 .and. rnsq.gt.xlarge)) then
         if (iprint .eq. 1) THEN
           write (msg,902) rnsq
           CALL rprint(msg)
         ENDIF

         iflag = -2
         return
      end if
      call dcopy(ncomp*nmsh, rhstri, 1, rhs, 1)

*  Statement 100 is the top of the iteration loop.

  100 continue

*  If rnsq is sufficiently small, terminate immediately.

      if (iprint .eq. 1) THEN
        write(msg,903) iter, rnsq
        CALL rprint(msg)
      ENDIF

      if (rnsq .le. 1.0d2*epsmch) then
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
      call fneval(ncomp, nmsh, xx, ncomp, utrial, fval, fsub,
     *   rpar, ipar)
      call rhscal (ncomp, nmsh, nlbc, xx, ncomp, utrial, defcor,
     *   fsub, gsub, rhstri, rnsq, fval, ftmp, uint, rpar, ipar)

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

      if (iter .ge. lmtfrz .or. rnsq .gt. (rnold/rngrow)
     *     .and. rnsq .gt. 1.0d2*epsmch ) then
         if (better) then

*  Setting iflag to -3 signals that, although the fixed Jacobian
*  iterations did not succeed, the current point was an improvement
*  on the previous one.  Hence, if we switch to a Newton procedure,
*  the right-hand side does not need to be recalculated.

            iflag = -3
         else
            iflag = -2
         endif
         if (iprint .eq. 1) THEN
           write(msg,904) iflag
           CALL rprint(msg)
         ENDIF

         return
      endif

*  Test for convergence using the ratio abs((change in u)/max(u,1)).

      do 150 im = 1, nmsh
      do 150 it = 1, ntol
         itol = ltol(it)
         er = abs(delu(itol,im))/max(abs(u(itol,im)), one)
         if (er .gt. tolfct*tol(it) .and. er .gt. 1.0d2*epsmch)go to 100
  150 continue

*  To exit from the loop here, the convergence tests have
*  been passed.

      if (iprint .ge. 0) THEN
        write(msg,905) iter, rnsq
        CALL rprint(msg)
      ENDIF


      iflag = 0
      return
  901 format(1h ,'fixed Jacobian iterations')
  902 format(1h ,'Large residual, rnsq =',1pe12.4)
  903 format(1h ,'iter, rnsq',i5,1pe11.3)
  904 format(1h ,'failure of fixed Jacobian, iflag =',i5)
  905 format(1h ,'fixed Jacobian convergence',i5,1pe11.3)
      end

c ===================================================================================

      subroutine lineq( ncomp, nmsh, nlbc,
     *    ludone, xx, nudim, u, defcor,
     *    delu, rhs, fval,
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs,
     *    ajac, topblk, botblk, bhold, chold, ipivot,
     *    fsub, dfsub, gsub, dgsub, iflag, rpar, ipar)
***********************************************************************
*     call by: bvpsol
*     calls to: lnrhs, jaccal, dcopy, colrow, dload, dload, clrslve
*           maxpy
**********************************************************************
      implicit double precision (a-h,o-z)
      dimension rpar(*), ipar(*)
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
      external   fsub
      external   dfsub
      external   gsub
      external   dgsub

      logical use_c, comp_c
      common/algprs/ nminit, iprint,idum, use_c, comp_c

*  blas: dcopy, dload

      parameter  ( one = 1.0d+0, zero = 0.0d+0 )
      character(len=150) msg

*  The routine lineq calculates the Newton step for a linear
*  problem.  The Newton step is exact unless the Jacobian
*  matrix is singular.

      isize=nmsh*ncomp
      ninter = nmsh - 1
      iflag = 0
      if (.not. ludone) then

*  Compute the right-hand side vector rhs.

         call lnrhs (ncomp, nmsh, nlbc, xx, nudim, u,
     *          fsub, gsub, rhs, rnsq, fval, ftmp, uint, rpar, ipar)


*  If the Jacobian for this mesh has not previously been
*  calulated and factorized successfully, call jaccal.
*  The block-structured Jacobian matrix is stored in three
*  matrices (topblk, ajac, and botblk).
*  The matrices bhold and chold are also calculated in jaccal,
*  and are saved for later use in outer routines.

         call jaccal (ncomp, nmsh, nlbc,
     *      xx, nudim, u, fval, dgtm, dftmp1, dftmp2, uint,
     *      ajac, topblk, botblk, bhold, chold,
     *      dfsub, dgsub,rpar,ipar)

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

         call dload(nlbc, zero, tmprhs(1), 1)
         do 100 im = 1, ninter
            loc = (im-1)*ncomp + nlbc + 1
            call dcopy(ncomp, defcor(1,im), 1, tmprhs(loc), 1)
  100    continue
         nrhs = ninter*ncomp + nlbc + 1
         call dload(ncomp-nlbc, zero, tmprhs(nrhs), 1)
         call dcopy(ncomp*nmsh,tmprhs,1,delu,1)
         job = 0

      call crslve(topblk,nlbc,ncomp,ajac,ncomp,2*ncomp,ninter,
     +   botblk,ncomp-nlbc,ipivot,delu,job)

      endif

*  Since the problem is linear, the Newton step  is exact.  The
*  new u array is obtained by adding delu to u.

      call maxpy ( ncomp, nmsh, one, delu, nudim, u )

c      iflag = 0
      return

  901 format(1h ,'Singular matrix')
      end

c ===================================================================================

      subroutine newteq(ncomp, nmsh, nlbc,
     *    rhsgiv, ntol, ltol, tol,
     *    xx, nudim, u, defcor,
     *    delu, rhs, fval,
     *    utrial, rhstri,
     *    uint, ftmp, dftmp1, dftmp2, dgtm, tmprhs, xmerit,
     *    ajac, topblk, botblk, bhold, chold, ipivot,
     *    fsub, dfsub, gsub, dgsub, iter, iflag, rpar, ipar, frscal)
******************************************************************
*     call by: bvpsol
*     calls to:
******************************************************************
      implicit double precision (a-h,o-z)
         dimension rpar(*), ipar(*)
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


      external   fsub
      external   dfsub
      external   gsub
      external   dgsub

      parameter  ( zero   = 0.0d+0, one    = 1.0d+0 )
      parameter ( two = 2.0d+0, half = 0.5d+0, fourth = 0.25d+0 )
      parameter ( tenth = 0.1d+0, ten = 10.0d+0, hund = 100.0d+0 )

      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c
      common/mchprs/flmin, flmax, epsmch

      logical gtpdeb, imprvd, braktd, crampd, extrap, vset, wset
      save  gtpdeb, mfsrch, epsaf, epsag, eta, rmu, tolabs, alfmax
      save  tolrel, toltny

      intrinsic  abs
      intrinsic  max

*  blas: dcopy

      logical frscal
c      save frscal
      parameter (cnvfct = 0.1d+0 )
c      data  frscal/.true./
c frscal has been updated in bvpsol
c
      data  alfsml/1.0d-4/,  alfmax/1.1d+0/
      data  imerit/1/, lmtnwt/39/
      data  shrfct/100.0d+0/, stpfct/2.0d+0/
      data  gtpdeb/.false./, mfsrch/5/
      data  eta/.999999d+0/, rmu/1.0d-6/

      character(len=150) msg

*  The routine newteq performs Newton iterations with a line
*  search, to solve the nonlinear equations.


*  Set up constants if this is the first call to newteq.

      if (frscal) then
         frscal = .false.
         epsaf = epsmch
         epsag = epsmch
         tolabs = epsmch
         tolrel = epsmch
         toltny = epsmch
      endif
      ninter = nmsh - 1

      if (iprint .eq. 1) THEN
        write(msg,901)
        CALL rprint(msg)
      ENDIF


*  A Newton method with line search and watchdog safeguarding
*  is used to try to solve the nonlinear equations.

*  Initialize iter (the counter of Newton iterations) and alfold
*  (the step taken at the previous iteration).

      iter = -1
      alfold = one
      alfa = zero



      if (.not. rhsgiv) then
*  If necessary, evaluate the right-hand side at the initial u.
         call fneval(ncomp, nmsh, xx, nudim, u, fval, fsub,
     *    rpar, ipar)
         call rhscal (ncomp, nmsh, nlbc, xx, nudim, u, defcor,
     *      fsub, gsub, rhs, rnsq, fval, ftmp, uint, rpar, ipar)
      endif


*  At any given Newton iteration, rnprev is the value of rnsq at
*  the immediately preceding Newton iteration.

      rnprev = flmax
      rnbest = flmax
      if (iprint .ge. 0) THEN
        write (msg,902)
        CALL rprint(msg)
      ENDIF


*  Initialize counter of watchdog iterations.

      itwtch = 0

*  Statement 100 is the top of the Newton iteration loop.

  100 continue

      iter = iter + 1

      if (iprint .eq. 1) THEN
        write(msg,910) iter
        CALL rprint(msg)
      ENDIF


*  If there have been too many Newton iterations, terminate.

      if (iter .ge. lmtnwt) then
         if (iprint .ge. 0) THEN
           write(msg,903)
           CALL rprint(msg)
         ENDIF

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

      call wtchdg ( iter, rnsq, rnbest, rnprev, itwtch,
     *                alfold, iflwat )

      if (iflwat .ne. 0) then
         if (iprint .ge. 0) THEN
           write(msg,904) iter
           CALL rprint(msg)
         ENDIF
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

      if (rnsq .le. epsmch .and. .not. comp_c ) then
         if (iprint .ge. 0) THEN
           write(msg,906) iter, rnsq
           CALL rprint(msg)
         ENDIF

         iflag = 0
         return
      endif

      call jaccal (ncomp, nmsh, nlbc,
     *    xx, nudim, u, fval, dgtm, dftmp1, dftmp2, uint,
     *    ajac, topblk, botblk, bhold, chold,
     *    dfsub, dgsub, rpar, ipar)

*  blkdcm is called to calculate the LU factors of the Jacobian,
*  which are overwritten on topblk, ajac and botblk.
*  Interchanges are represented in the integer array ipivot.

*   Solve for the Newton step delu.  Copy the rhs array into tmprhs,
*   which is then overwritten by blkslv.

      call dcopy(ncomp*nmsh, rhs, 1, tmprhs, 1)
      call dcopy(ncomp*nmsh,tmprhs,1,delu,1)
      job = 0
      call colrow(nmsh*ncomp,topblk,nlbc,ncomp,ajac,ncomp,2*ncomp,
     +   ninter,botblk,ncomp-nlbc,ipivot,delu,iflag,job)

      if (iprint .ge. 0 .and. iflag.ne.0) THEN
        write(msg,905) iter
        CALL rprint(msg)
      ENDIF

C KSKS to account for the NOT COMP_C this was added...
      if (rnsq .le. epsmch )    return


      if (iflag .ne. 0) return
*        the jacobian is singular


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
      if (iprint .eq. 1) THEN
        write (msg,908) alfa, fmtry, rnsq
        CALL rprint(msg)
      ENDIF


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


  150 continue
      iwr = 6
      call getptq (gtpdeb, mfsrch, iwr, alfmax, alfsml, alfuzz,
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
         call fneval(ncomp, nmsh, xx, ncomp, utrial, fval, fsub,
     *    rpar, ipar)
         call rhscal (ncomp, nmsh, nlbc, xx, ncomp, utrial, defcor,
     *      fsub, gsub, rhstri, rnsqtr, fval, ftmp, uint, rpar, ipar)

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
         if (iprint .eq. 1) THEN
           write (msg,908) alfa, fmtry, rnsqtr
           CALL rprint(msg)
         ENDIF

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
      if (iprint .ge. 0) THEN
        write(msg,909) iter, alfa, fmtry, rnsq
        CALL rprint(msg)
      ENDIF


*  Now test for convergence using the ratio of the Newton step
*  for each component with max(1, abs(current solution estimate)).
*  If the test fails for any element of u, branch back to the
*  top of the Newton iteration.

      do 160 im = 1, nmsh
      do 160 it = 1, ntol
         icmp = ltol(it)
         er = abs(delu(icmp,im))/max(abs(u(icmp,im)), one)
         if (er .gt. cnvfct*tol(it)) go to 100
  160 continue


      if (iprint .ge. 0) THEN
        write(msg, 906) iter+1, rnsq
        CALL rprint(msg)
      ENDIF

      iflag = 0

*  To fall through the above loop, the termination test for a
*  sufficiently small delu is satisfied.
*  Note that the stored Jacobian and its factorization do not
*  correspond to the final solution.


      return

  901 format(1h ,'start Newton iterations')
  902 format(1h ,' iter',
     *              7x,'alfa',6x,'merit',7x,'rnsq')
  903 format(1h ,'Too many Newton iterations')
  904 format(1h ,'Watchdog tests fail, iter =', i5)
  905 format(1h ,'Singular Jacobian, iter=',i5)
  906 format(1h ,'Convergence, iter =',i5,4x,'rnsq =',1pe12.3)
  908 format(1h ,'alfa, merit, rnsq',3(1pe11.3))
  909 format(1h ,i5,3(1pe11.3))
  910 format(1h ,'Newton iteration',i5)
      end

c ===================================================================================

      subroutine wtchdg ( iter, wmerit, wmbest, wmprev,
     *      itwtch, alfold, iflag )

*  Logic for watchdog tests.

      implicit double precision (a-h,o-z)
      parameter ( itonew = 5, itwtmx = 8, grfct = 100.0d+0 )
      parameter ( half = 0.5d+0 )
      character(len=150) msg

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

*      write(6,99) iter, wmerit, wmbest, wmprev
*      write(6,98) itwtch, alfold
*   99 format(1h ,'iter,wmer,wbest,wprev',i5,3(1pe15.5))
*   98 format(1h ,'itwtch,alfold',i5,1pe15.5)
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

c ===================================================================================

      subroutine fneval(ncomp, nmsh, xx, nudim, u, fval, fsub,
     *   rpar, ipar)
      implicit double precision (a-h,o-z)
      dimension rpar(*),ipar(*)
      dimension xx(nmsh), u(nudim,nmsh), fval(ncomp,nmsh)
      integer nfunc, njac, nstep, nbound, njacbound
      common/diagnost/nfunc, njac, nstep, nbound, njacbound

      external fsub

*  fneval evaluates the function values (from fsub) for
*  a given mesh xx and array u, and stores the values
*  in the array fval.

      call fsub (ncomp, xx(1), u(1,1), fval(1,1), rpar, ipar)

      do 50 im = 1, nmsh-1
         hmsh = xx(im+1) - xx(im)
         call fsub (ncomp, xx(im+1), u(1,im+1), fval(1,im+1),
     *     rpar, ipar)
   50 continue
      nstep = nstep + 1
      nfunc = nfunc + nmsh

      return
      end

c ===================================================================================

      subroutine jaccal (ncomp, nmsh, nlbc,
     *   xx, nudim, u, fval,
     *   dgtm, dftm1, dftm2, uint,
     *   ajac, topblk, botblk, bhold, chold,
     *   dfsub, dgsub, rpar, ipar)

      implicit double precision (a-h,o-z)
       dimension rpar(*), ipar(*)
      dimension xx(nmsh), u(nudim,nmsh), fval(ncomp,nmsh)
      dimension dgtm(ncomp)
      dimension dftm1(ncomp, ncomp), dftm2(ncomp, ncomp),
     *             uint(ncomp)
      dimension ajac(ncomp, 2*ncomp, nmsh-1)
      dimension topblk(nlbc, ncomp), botblk(ncomp-nlbc,ncomp)
      dimension bhold(ncomp, ncomp, nmsh-1),
     *             chold(ncomp, ncomp, nmsh-1)

      external  dfsub
      external  dgsub

      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c
      integer nfunc, njac, nstep, nbound, njacbound
      common/diagnost/nfunc, njac, nstep, nbound, njacbound

*  blas: dcopy, ddot

      parameter ( zero = 0.0d+0, half = 0.5d+0, eighth = 0.125d+0 )
      parameter ( four = 4.0d+0, six = 6.0d+0 )
      parameter ( one = 1.0d+0, three = 3.0d+0, twelve = 12.0d+0 )
      character(len=150) msg


      ninter = nmsh - 1

      do 110 i = 1, nlbc
         call dgsub (i, ncomp, u(1,1), dgtm,rpar,ipar)
         call dcopy(ncomp, dgtm(1), 1, topblk(i,1), nlbc)
  110 continue

      call dfsub (ncomp, xx(1), u(1,1), dftm1(1,1), rpar, ipar)

*  on entry to jaccal, the array fval contains the function values
*  at (xx(im), u(ic,im)), ic=1,...,ncomp and im = 1,...,nmsh,
*  calculated by a preceding call of rhscal with the same xx and u
*  arrays.

      do 200 im = 1, ninter

         hmsh = xx(im+1) - xx(im)

         do 120 ic = 1, ncomp
            uint(ic) = half*(u(ic,im) + u(ic,im+1))
     *         - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))
  120    continue
         xhalf = half*(xx(im+1) + xx(im))
         call dfsub (ncomp, xhalf, uint, dftm2(1,1), rpar, ipar)
         do 140 ic = 1, ncomp
            do 130 jc = 1, ncomp
               dsq = ddot(ncomp, dftm2(ic,1), ncomp,
     *                       dftm1(1,jc), 1)
               ajac(ic,jc,im) = -hmsh*(dftm1(ic,jc)/six
     *             + dftm2(ic,jc)/three + hmsh*dsq/twelve)
  130       continue
            ajac(ic,ic,im) = ajac(ic,ic,im) - one
  140    continue

         call dfsub (ncomp, xx(im+1), u(1,im+1), dftm1(1,1), rpar, ipar)
         do 170 ic = 1, ncomp
            do 160 jc = 1, ncomp
               dsq = ddot(ncomp, dftm2(ic,1), ncomp,
     *                         dftm1(1,jc), 1)
               ajac(ic,jc+ncomp,im) = -hmsh*(dftm1(ic,jc)/six
     *               + dftm2(ic,jc)/three - hmsh*dsq/twelve)
  160       continue
            call dcopy(ncomp, ajac(ic,ncomp+1,im), ncomp,
     *                   chold(ic,1,im), ncomp)
            call dcopy(ncomp, dftm1(ic,1), ncomp,
     *                   bhold(ic,1,im), ncomp)
            ajac(ic,ic+ncomp,im) = ajac(ic,ic+ncomp,im) + one
            chold(ic,ic,im) = ajac(ic,ic+ncomp,im)
  170    continue


  200 continue

      njac = njac + 1 + ninter*2

      do 220 i = nlbc+1, ncomp
         call dgsub (i, ncomp, u(1, nmsh), dgtm,rpar,ipar)
         call dcopy(ncomp, dgtm(1), 1, botblk(i-nlbc,1), ncomp-nlbc)
  220 continue

      njacbound = njacbound + ncomp

      return

      end

c ===================================================================================

      subroutine lnrhs (ncomp, nmsh, nlbc,
     *   xx, nudim, u, fsub, gsub,
     *   rhs, rnsq, fval, ftmp, uint, rpar, ipar)

       implicit double precision(a-h,o-z)

*  This subroutine is designed to calculate the right-hand
*  side for linear problems.
      dimension rpar(*),ipar(*)
      dimension xx(*), u(nudim,*)
      dimension rhs(*), fval(ncomp,*), ftmp(*), uint(*)
      external fsub
      external gsub

      parameter ( zero = 0.0d+0, half = 0.5d+0, eighth = 0.125d+0 )
      parameter ( one = 1.0d+0, four = 4.0d+0, six = 6.0d+0 )

      common/mchprs/flmin, flmax, epsmch
      integer nfunc, njac, nstep, nbound, njacbound
      common/diagnost/nfunc, njac, nstep, nbound, njacbound
      intrinsic abs
      character(len=150) msg

*  blas: dssq

      ninter = nmsh - 1
      rnsq = zero

*  first, process the left-hand boundary conditions.

      do 20 i = 1, nlbc
         call gsub (i, ncomp, u(1,1), wg,rpar,ipar)
         rhs(i) = -wg
   20 continue

*  Next, process the interior mesh points.  The fval array
*  contains the function values from fsub at xx and u.

      do 50 im = 1, ninter
         hmsh = xx(im+1) - xx(im)
         do 30 ic = 1, ncomp
            uint(ic) = half*(u(ic,im) + u(ic,im+1))
     *         - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))
   30    continue
         xhalf = half*(xx(im) + xx(im+1))
         call fsub (ncomp, xhalf, uint, ftmp, rpar, ipar)
         loc = (im-1)*ncomp + nlbc
         do 40 ic = 1, ncomp
            rhs(loc+ic) = -u(ic,im+1) + u(ic,im)
     *        + hmsh*
     *            (fval(ic,im) + fval(ic,im+1) + four*ftmp(ic))/six
   40    continue
   50 continue

      nfunc = nfunc + ninter

      nrhs = ninter*ncomp
      do 60 ii = nlbc+1, ncomp
         call gsub (ii, ncomp, u(1,nmsh), wg,rpar,ipar)
         rhs(nrhs+ii) = -wg
   60 continue

      call dssq  ( nmsh*ncomp, rhs, 1, scale, sumsq )
      rnsq = (scale**2)*sumsq

      nbound = nbound + ncomp
      return
      end

c ===================================================================================

      subroutine rhscal (ncomp, nmsh, nlbc,
     *   xx, nudim, u, defcor,
     *   fsub, gsub,
     *   rhs, rnsq, fval, ftmp, uint, rpar, ipar)

       implicit double precision(a-h,o-z)

*  This subroutine constructs the (ncomp*nmsh)-dimensional
*  vector rhs, which is the right-hand side of the Newton equations.
*  The ncomp by nmsh array fval is assumed to have been calculated
*  elsewhere by routine fneval.
              dimension rpar(*),ipar(*)
      dimension  xx(nmsh), u(nudim,nmsh), defcor(ncomp,nmsh-1)
      dimension  rhs(ncomp*nmsh), fval(ncomp,nmsh)
      dimension  ftmp(ncomp), uint(ncomp)
      external   fsub
      external   gsub

      logical use_c, comp_c
      common/algprs/ nminit, iprint,idum,  use_c, comp_c
      common/mchprs/flmin, flmax, epsmch
      intrinsic abs
      integer nfunc, njac, nstep, nbound, njacbound
      common/diagnost/nfunc, njac, nstep, nbound, njacbound

*  blas: dssq

      parameter ( zero = 0.0d+0, half = 0.5d+0, eighth = 0.125d+0 )
      parameter ( one = 1.0d+0, four = 4.0d+0, six = 6.0d+0 )
      character(len=150) msg

*  ninter is the number of intervals in the mesh (one less than the
*  number of mesh points)

      ninter = nmsh - 1
      rnsq = zero

*  First, process the left-hand boundary conditions.

      do 20 i = 1, nlbc
         call gsub (i, ncomp, u(1,1), wg,rpar,ipar)
         rhs(i) = -wg
   20 continue

*  Next, process the interior mesh points.  The fval array
*  contains the function values from fsub at xx and u.

      do 50 im = 1, ninter
         hmsh = xx(im+1) - xx(im)
         do 30 ic = 1, ncomp
            uint(ic) = half*(u(ic,im) + u(ic,im+1))
     *         - eighth*hmsh*(fval(ic,im+1) - fval(ic,im))
   30    continue
         xhalf = half*(xx(im) + xx(im+1))
         call fsub (ncomp, xhalf, uint, ftmp, rpar, ipar)
         loc = (im-1)*ncomp + nlbc
         do 40 ic = 1, ncomp
            rhs(loc+ic) = -u(ic,im+1) + u(ic,im) + defcor(ic,im)
     *        + hmsh*
     *            (fval(ic,im) + fval(ic,im+1) + four*ftmp(ic))/six

   40    continue
   50 continue
      nfunc = nfunc + ninter

      nrhs = ninter*ncomp
      do 60 ii = nlbc+1, ncomp
         call gsub (ii, ncomp, u(1,nmsh), wg,rpar,ipar)
         rhs(nrhs+ii) = -wg
   60 continue

      nbound = nbound + ncomp

      call dssq  ( nmsh*ncomp, rhs, 1, scale, sumsq )
      rnsq = (scale**2)*sumsq

      return
      end

c ===================================================================================

      subroutine dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
      implicit double precision (a-h,o-z)
      dimension xx(*), xxold(*)
      logical maxmsh

      logical use_c, comp_c
      common/algprs/ nminit, iprint,idum, use_c, comp_c

*  blas: dcopy

      parameter (half = 0.5d+0)
      character(len=150) msg


*  This routine is used to double the mesh, i.e., produce a mesh
*  with twice as many intervals in which each new interval is
*  half the corresponding old interval.

*  On entry to dblmsh, the integer nmsh and the array xx
*  specify a set of mesh points xx(1),..., xx(nmsh) (assumed
*  to be in ascending order).

*  If the number of mesh points in the doubled mesh would
*  exceed the maximum allowed number nmax, the flag maxmsh is
*  set to true, and we exit without changing any other parameters.

*  Otherwise, nmold is set to the old number of mesh points,
*  xxold is set to the old mesh, nmsh is the new number of mesh
*  points, and xx contains the new mesh points.

      nmold = nmsh
      call dcopy(nmold, xx, 1, xxold, 1)

      ninnew = 2*(nmsh-1)
      nmnew = ninnew + 1
      if(nmnew .ge. nmax) then
         if (iprint .ge. 0) THEN
           write(msg,901) nmnew
           CALL rprint(msg)
         ENDIF

         nmsh = nmold
         call dcopy(nmold, xxold, 1, xx, 1)
         maxmsh = .true.
         return
      endif
      maxmsh = .false.

*  Loop backwards through the old mesh points to create the new ones.

      xx(nmnew) = xx(nmsh)
      do 100 i = ninnew, 4, -2
         id2 = i/2
         xx(i) = half*(xx(i+1) + xx(id2))
         xx(i-1) = xx(id2)
  100 continue

*  Calculate the new xx(2). xx(1) remains unchanged.
      xx(2) = half*(xx(3) + xx(1))
      nmsh = nmnew
      if(iprint .ge. 0) THEN
        write(msg,902) nmsh
        CALL rprint(msg)
      ENDIF

      return
  901 format (1h , ' dblmsh.  maximum mesh exceeded, nmnew =', i8)
  902 format (1h , ' dblmsh.  the doubled mesh has ', i8,' points.')
      end

c ===================================================================================

      subroutine selmsh(ncomp, nmsh, ntol, ltol, tol,
     *     nfxpnt, fixpnt, ipow, nmax,
     *     xx, nudim, u, ermeas, irefin, ihcomp,
     *     nmold, xxold, ermx, ddouble , maxmsh)

      implicit double precision (a-h,o-z)

      dimension  ltol(ntol), tol(ntol), fixpnt(*)
      dimension  xx(*), u(nudim, *), ermeas(ncomp,*)
      dimension  irefin(nmsh-1), ihcomp(nmsh-1)
      dimension  xxold(*), ermx(*)
      logical    ddouble , maxmsh

      intrinsic abs
      intrinsic  max
      intrinsic  int
      character(len=150) msg

*  blas: dcopy
*  double precision dlog

      logical use_c, comp_c
      common/algprs/ nminit, iprint,idum, use_c, comp_c

      parameter  ( zero = 0.0d+0, one = 1.0d+0, onep1 = 1.1d+0 )
      parameter  ( erdcid = 5.0d+0 )
      parameter  ( phitst = 0.1d+0 )

      logical first
      save    first, rlndec
      data    first / .true. /

*  The routine selmsh performs selective mesh refinement, depending
*  on the error measure ermeas.

      if (first) then
         first = .false.
         rlndec = dlog(erdcid)
      endif

      maxmsh = .false.

      frcpow = one/ipow
      ddouble = .false.
      nmold = nmsh
      ninter = nmsh - 1

*  Copy the current mesh into the xxold array.


      call dcopy(nmold, xx, 1, xxold, 1)
      ithres = 0
      thres = one

*  On input, the array ermeas represents some error measure defined
*  over the components and mesh intervals (not mesh points).
*  It is normalized in the following loop with respect to the
*  tolerance array and the current solution.
*  The value errmax gives the maximum normalized error.

      errmax = zero
      do 120 im = 1, ninter
         ermx(im) = zero
         do 110 it = 1, ntol
            jcomp = ltol(it)
            denom = tol(it)*max(one, abs(u(jcomp,im)))
            ems = ermeas(jcomp,im)
            ermeas(jcomp,im) = abs(ems)/denom
            err = ermeas(jcomp, im)
            if (err .ge. ermx(im)) then
                ermx(im) = err
                ihcomp(im) = jcomp
            endif
  110    continue
         errmax = max(ermx(im), errmax)
  120 continue


      if (errmax .gt. zero .and. errmax .le. erdcid) then

*  If errmax > 0 and .le. erdcid, find the smallest integer exponent ii
*  such that (erdcid**ii)*errmax > erdcid.

         if(errmax .gt. one) then
            ii = 1
            decii = erdcid
         else
            ilg = -dlog(errmax)/rlndec
            ii = 2 + ilg
            decii = erdcid**ii
         endif

*  Multiply error measures by erdcid**ii.

         errmax = decii*errmax
         do 140 im = 1, ninter
            ermx(im) = decii*ermx(im)
            do 140 it = 1, ntol
               jcomp = ltol(it)
               ermeas(jcomp,im) = decii*ermeas(jcomp,im)
  140    continue
      endif

  200 continue

*  For each interval im,  the integer irefin(im) is calculated
*  based on an equidistrbution procedure involving the
*  threshold value thres.  If irefin(im) > 1, we add
*  points to interval im.  If irefin(im) = 1, we check whether
*  point im can be removed from the mesh.

*  nmest is a lower bound on the number of points in the new mesh.
*  We do not know in advance how many points will be removed,
*  so nmest is computed by assuming all eligible points are removed.

      nmest = nmsh
      do 220 im = 1, ninter
         if (ermx(im).ge.thres) then
            irefin(im) = int(ermx(im)**frcpow) + 1
            nmest = nmest + irefin(im) - 1
         else
            irefin(im) = 1
            nmest = nmest - 1
         endif
  220 continue

      if (nmest .gt. nmax) then

         go to 360

      elseif (nmest-1 .gt. 3*ninter) then

         call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
         ddouble = .true.
         return
      endif

*  It appears that we can perform the desired selective mesh
*  refinement.

*  Now begin running through the mesh, adding and possibly deleting
*  points as indicated by the irefin array.

*  The integer new is a count of the number of intervals in
*  the tentative mesh being generated by the refinement strategy.

      new = 1

*  The first interval is treated as a special case, since xx(1)
*  always remains in the mesh, and cannot be a fixed point.

      rlen = xxold(2) - xx(1)
      slen = rlen
      if (new + irefin(1)  .gt. nmax) goto 360
      if (irefin(1).gt.1) then
         dx = rlen/irefin(1)
         do 230 j = 2, irefin(1)
            new = new + 1
            xx(new) = xx(1) + (j-1)*dx
  230    continue
      endif

*  The fixed points specified by the fixpnt array cannot be
*  removed from the mesh.  The value fxnext indicates the 'next'
*  fixed point.  When no further fixed points remain to be processed
*  (or if nfxpnt = 0), fxnext is set to a value strictly larger than
*  the last mesh point, so that no mesh point can equal fxnext.
*  This way we need to compare only the values of xxold(i)
*  and fxnext.

      ifxcnt = 1
      if (nfxpnt .eq. 0) then
         fxnext = onep1*abs(xxold(nmsh))
      else
         fxnext = fixpnt(ifxcnt)
      endif

*  jtkout is a counter of the number of consecutive points that
*  have been removed from the mesh.

      jtkout = 0
      do 330 im = 2, ninter
         rlold = rlen
         rlen = xxold(im+1) - xxold(im)

*  If xxold(im) is the next fixed point, it cannot be removed
*  and so we don't test its error estimates.

         if(xxold(im) .eq. fxnext)  then

            ifxcnt = ifxcnt + 1
            if(ifxcnt .gt. nfxpnt) then
               fxnext = onep1*abs(xxold(nmsh))
            else
               fxnext = fixpnt(ifxcnt)
            endif

         elseif (irefin(im).eq.1) then

*  If xxold(im) is not a fixed point and irefin(im) = 1, we may wish
*  to remove point im from the mesh.

*  If we are considering removing points and jtkout = 0, this
*  is the first point in a possible consecutive set to be removed,
*  and we initialize phihat, which represents a maximum of
*  certain estimates.
*  If jtkout is not zero, previous points contiguous to this
*  point have been removed, and phihat retains its previous value.

            slen = slen + rlen

            if (jtkout .eq. 0) then
                ind1 = ihcomp(im-1)
                phihat = ermeas(ind1,im-1)/(rlold**ipow)
            endif
            phihat = max(phihat,
     *                 ermeas(ihcomp(im),im)/(rlen**ipow))
            val1 = phihat*(slen**ipow)
            if (val1 .le. phitst
     *             .and. jtkout .lt. 4) then

*  Increment the counter of removed points.
*  'Remove' the mesh point xxold(im) by not including it.

               jtkout = jtkout+1
               go to 330
            endif
*        end of logic for irefin(im) = 1.
         endif

         jtkout = 0
         new = new + 1
         xx(new) = xxold(im)

         if (new + irefin(im)  .gt. nmax) goto 360

         if (irefin(im) .gt. 1) then
            dx = rlen/irefin(im)
            do 300 j = 2, irefin(im)
              new = new + 1
              xx(new) = xxold(im) + (j-1)*dx
  300       continue
         endif
         slen = rlen

         if (new  .gt. nmax) then

*  If the new mesh contains too many points, branch out of the
*  loop to try alternative strategies.

            go to 360

         elseif (new .gt. 3*ninter)  then

*  Here, the new mesh does not exceed the specified maximum,
*  but has more than 3 times as many intervals as the old mesh.
*  Try doubling the mesh if possible.

            call dcopy(nmsh, xxold, 1, xx, 1)
            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
            ddouble = .true.
            return

          end if
  330 continue

*  To end up here, we have processed the entire interval,
*  and have neither exceeded the specified maximum nor
*  exceeded three times the number of intervals in the old
*  mesh.  The last mesh point remains unchanged.

      new = new + 1
      xx(new) = xxold(nmsh)
      nmsh = new
      maxmsh = .false.
      if (iprint .ge. 0) THEN
        write(msg,905) nmsh
        CALL rprint(msg)
      ENDIF

      return

  360 continue
*  To reach here, the number of mesh points created at some stage
*  of the refinement process was larger than the maximum permitted
*  value nmax.

*  Check whether the mesh can safely be doubled.

      if ((2*nmsh-1) .lt. nmax) then

*  Double the mesh.
         call  dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
         ddouble = .true.

*  If the number of intervals is too large and the mesh cannot be
*  doubled, increase the threshold thres by a factor of erdcid and
*  try the selective refinement again.
*  If this happens three times without success or if thres exceeds
*  or is equal to errmax, stop.  (In this case, we know already
*  that doubling the mesh produces too many points.)

      elseif (thres .lt. errmax .and. ithres .lt. 3) then
         ithres = ithres + 1
         thres = erdcid*thres
         if(thres .gt. errmax) thres = errmax
         call dcopy(nmsh, xxold, 1, xx, 1)
         go to 200
      else
c         nmsh = 2*nmsh - 1
         nmsh = nmold
         call dcopy(nmold, xxold, 1, xx, 1)
         maxmsh = .true.
      endif
      return

  902 format(1h ,'im, jcomp, ermeas, normalized er',2i5,2(1pe11.3))
  903 format(1h ,'errmax',1pe11.3)
  905 format(1h ,'selmsh.  new nmsh =',i8)
  910 format(1h ,'ihcomp',(10i5))
      end
c       selmsh

c ===================================================================================

       subroutine smpmsh (nmsh, nmax, xx, intref, numadd,
     *      nmold, xxold, maxmsh)

      implicit double precision (a-h,o-z)
      logical maxmsh
      dimension xx(*), xxold(*)

      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum,use_c, comp_c
      character(len=150) msg

*  blas: dcopy

*  The routine smpmsh performs simple mesh refinement by adding
*  points to one or three interval(s) in the region indicated
*  by the integer intref.
*  numadd gives the trial number of points to be added in each
*  interval.

      nmold = nmsh
      call dcopy(nmold, xx, 1, xxold, 1)

*  numadd is altered if necessary so that it lies between 4 and 49

      if(numadd .gt. 49) then
         numadd = 49
      elseif (numadd .lt. 4) then
         numadd = 4
      endif


      maxmsh = .false.
      if (intref .eq. 1) then

*  Add numadd points to the first interval if intref = 1.

         nmnew = nmsh + numadd
         if (nmnew .gt. nmax) then
            if (iprint .ge. 0) THEN
              write(msg,903) nmnew
              CALL rprint(msg)
            ENDIF
            maxmsh = .true.
            return
         endif

*  Renumber the later points in reverse order.

         nint = numadd + 1
         do 100 i = nmnew, numadd+2, -1
            xx(i) = xx(i-numadd)
  100    continue
         dx = (xx(2) - xx(1))/nint
         do 110 i = 2, nint
            xx(i) = xx(1) + (i-1)*dx
  110    continue

      elseif (intref .eq. nmsh-1) then

*  Add numadd points to the last interval if intref = nmsh-1.

         nmnew = nmsh + numadd
         if (nmnew .gt. nmax) then
            if (iprint .ge. 0) THEN
               write(msg,903) nmnew
               CALL rprint(msg)
            ENDIF

            maxmsh = .true.
            return
         endif
         nint = numadd + 1
         dx = (xx(nmsh) - xx(nmsh-1))/nint
         xx(nmnew) = xx(nmsh)
         do 200 i = nmsh, nmnew-1
            xx(i) = xx(nmsh-1) + (i-nmsh+1)*dx
  200    continue

      else

         if (numadd .gt. 9) numadd = 9

*  Here, intref lies between 2 and nmsh-2.  Add numadd points to
*  each of the three intervals intref-1, intref and intref+1.

         nmnew = nmsh + 3*numadd
         if (nmnew .gt. nmax) then
            if (iprint .ge. 0) THEN
              write(msg,903) nmnew
              CALL rprint(msg)
            ENDIF
            maxmsh = .true.
            return
         endif

*  noalt is the number of points at the right end of the interval
*  whose numerical values remain the same, but whose indices change.
*  nochsm is the smallest index in the new ordering of one of these
*  points.

         noalt = nmsh - intref - 1
         nochsm = nmnew - noalt + 1

*  Renumber the noalt unchanged points at the right end of the interval
*  (in reverse order).

         j = 0
         do 300 i = nmnew, nochsm, -1
            xx(i) = xx(nmsh-j)
            j = j + 1
  300    continue

*  Add numadd points to the three contiguous intervals.
*  The remaining points at the left end of the interval retain
*  their original indices, and are left unchanged.

         nint = numadd + 1
         innew = nochsm - nint
         do 320 i = intref+1, intref-1, -1
            xx(innew) = xx(i)
            dx = (xx(innew + nint) - xx(innew))/nint
            do 310 j = 1, numadd
               xx(innew + j) = xx(innew) + j*dx
  310       continue
            innew = innew - nint
  320    continue

      endif
      nmsh = nmnew

      if(iprint .ge. 0) THEN
        write(msg,904) nmsh
        CALL rprint(msg)
      ENDIF

      return
  903 format(1h , ' smpmsh.  maximum points exceeded, nmnew =',i6)
  904 format(1h ,'smpmsh, new nmsh =',i7)
      end

c ===================================================================================

      subroutine unimsh(nmsh, aleft, aright, nfxpnt, fixpnt, xx)
      implicit double precision (a-h,o-z)
      integer  nmsh, nfxpnt
      dimension fixpnt(*), xx(nmsh)

      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c

      intrinsic max
      character(len=150) msg


*  Given a left endpoint aleft, a right endpoint aright,
*  a set of nfxpnt fixed points fixpnt(i), i = 1,...,nfxpnt,
*  (where fixpnt(i) is different from aleft and aright for all i),
*  and an initial target number nmsh of mesh points,
*  the subroutine unimsh generates a piecewise uniform mesh
*  beginning at aleft, ending at aright, and with equally
*  spaced points between aleft and fixpnt(1), then between
*  fixpnt(1) and fixpnt(2), ..., and finally between
*  fixpnt(nfxpnt) and aright.  The final number of intervals
*  is the maximum of nfxpnt+2 and the initial value of nmsh.

*  In the simplest case when nfxpnt = 0, unimsh generates a
*  uniform mesh with nmsh intervals in the closed interval
*  (aleft, aright).

*  On exit, the integer nmsh contains the number of mesh points
*  (which is the maximum of the initial nmsh and nfxpnt).
*  The array xx (of dimension nmsh) contains the mesh points.

      if (iprint .ge. 0) THEN
        write(msg,901) nmsh
        CALL rprint(msg)
      ENDIF


      if (nfxpnt .eq. 0) then

*  If there are no interior fixed points, the spacing is uniform
*  throughout the interval.  Calculate the spacing dx
*  and set up the xx array.

        ninter = nmsh - 1

         dx = (aright - aleft)/ninter
         do 10 i = 1, ninter
            xx(i) = aleft + (i-1)*dx
   10    continue
         xx(nmsh) = aright
         return
      endif

*  We know that there is at least one fixed point strictly between
*  the endpoints.

      if (nmsh .lt. nfxpnt+2)  nmsh = nfxpnt + 2
      ninter = nmsh - 1
      xx(1) = aleft
      ileft = 1
      xleft = aleft
      totint = aright - aleft
      ndif = ninter - nfxpnt
      do 50 j = 1, nfxpnt + 1

*  Deal in turn with the subintervals defined by the interval
*  boundaries and the fixed  points.

         if (j .lt. nfxpnt+1) then

*  The j-th fixed point is xright.  Calculate where it should
*  fall in the mesh.

            xright = fixpnt(j)
            nmin = ninter*(xright-aleft)/totint + 1.5d+0
            if (nmin .gt. ndif+j) nmin = ndif + j
            iright = max(ileft+1, nmin)
         else
            xright = aright
            iright = nmsh
         endif

*  npt is the number of equally spaced points that should
*  lie strictly between the (j-1)-th and j-th fixed points.

         xx(iright) = xright
         npt = iright - ileft - 1
         dx = (xright - xleft)/(npt + 1)
         do 30 i = 1, npt
            xx(ileft+i) = xleft + i*dx
   30    continue
         ileft = iright
         xleft = xright
   50 continue

      return
c
  901 format (1h ,'unimsh.  nmsh =',i5)
      end

c ===================================================================================

      subroutine stats(len, elem, ebigst, esecnd, summod, index)
      implicit double precision (a-h, o-z)
      integer index
      dimension elem(len)

      intrinsic abs

      parameter  ( zero = 0.0d+0 )
      character(len=150) msg

*  Given the real array elem of length len, stats calculates
*  the following:
*      - summod, the sum of the magnitudes of the elements of elem;
*      - ebigst (the largest element in magnitude);
*      - index (the index in elem of ebigst); and
*      - esecnd (the second largest element in magnitude, strictly
*          less than ebigst unless both are zero).

      index = 1
      ebigst = zero
      esecnd = zero
      summod = zero

      do 100 i = 1, len
         elmod = abs(elem(i))
         summod = summod + elmod
         if (elmod .gt. ebigst) then
            esecnd = ebigst
            ebigst = elmod
            index = i
         elseif (elmod .gt. esecnd) then
            esecnd = elmod
         endif
  100 continue
      return
      end

c ===================================================================================

      subroutine mshref(ncomp, nmsh, nlbc, ntol, ltol,
     *              iorder, rhs, tmwork,
     *              nmax, xx, nmold, xxold, ddouble , maxmsh,
     *              numbig, nummed,amg,stab_cond,stiff_cond,r4,
     *              nfxpnt, fixpnt,irefin,itcond,itcondmax)

      implicit double precision (a-h,o-z)

      dimension  ltol(ntol)
      dimension  rhs(ncomp*nmsh), tmwork(nmsh-1)
      dimension  xx(nmsh), xxold(nmold)
      dimension  fixpnt(*), irefin(*),r4(*)
      logical    ddouble , maxmsh
      dimension  amg(*)
      logical stab_cond, stiff_cond, forcedouble
      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c

*  blas: dcopy
      logical nodouble
      parameter (two = 2.0d+0)
      parameter (bigfac = 10.0d+0, small = 1.0d-2, numpt = 14)
      character(len=150) msg


*  This routine performs calculations leading to a decision
*  about how the mesh will be refined, and then refines the mesh.

*  The choices for mesh refinement in this routine are either to
*  double the mesh (i.e., divide each existing interval in half),
*  or to add points to just a few intervals.

*  The decision is made based on two criteria:  (1) the distribution
*  of magnitudes of components of rhs (broadly speaking, if
*  the maximum component of rhs is much larger than the
*  average, points are added near the corresponding interval),
*  and (2) the history of previous mesh refinements (if points
*  have been added to only a few intervals already, this strategy
*  is abandoned for the moment and the mesh is doubled).

*  The decision is indicated by setting the logical flag double
*  and (if ddouble is .false.) the integer intref.
*  Setting ddouble to .true. means that the mesh should be doubled
*  (i.e., the new mesh should contain twice as many intervals).
*  The integer intref indicates the region of the mesh where the
*  points should be added (see the routine smpmsh), and numadd
*  indicates how many points are to be added.
*  The integers nummed and numbig represent running totals,
*  used in deciding on the mesh refinement strategy.

*  If iorder = 4, meaning that we were just performing a Newton
*  iteration for a 4th order solution, check all elements of rhs.
*  If iorder .gt. 4, signalling that we were trying Newton
*  iterations for order 6 or 8, check only elements of rhs
*  corresponding to components for which a tolerance is specified.

      nodouble = ((iorder.eq.4) .and.
     *       (stiff_cond .and. .not. stab_cond) .and. (use_c))
c      nodouble = nodouble
c     *  .or.((iorder.gt.4) .and. ( stiff_cond .and. .not. stab_cond)
c     *   .and. (use_c))



      forcedouble = .false.
      if (use_c) then
         if ( itcond .eq. itcondmax) then
            itcond = 0
            forcedouble = .true.
         else
            forcedouble = .false.
         endif
      endif


      ninter = nmsh-1
      nup = ncomp
      if (iorder .gt. 4) nup = ntol

*  Check the vector rhs for a non-negligible component
*  whose magnitude is significantly larger than the average.
*  (small defines negligible, and bigfac defines significantly larger.)

      do 50 ic = 1, nup
         icmp = ic
         if (iorder .gt. 4) icmp = ltol(ic)

*  For component icmp, examine the ninter elements of rhs not
*  corresponding to boundary conditions.

*  subroutine stats calculates rbigst and rsecnd (the first- and
*  second-largest elements of rhs in magnitude), and the index
*  intref of the interval in which the largest value occurs.
*  The value sumrhs is the sum of the magnitudes of the components
*  of rhs.

         indrhs = nlbc + icmp
*FM: changed ic with icmp in the previous line
*  Copy the elements of rhs corresponding to interior mesh
*  points for component icmp into a single vector tmwork.

         call dcopy(ninter, rhs(indrhs), ncomp, tmwork, 1)

         call stats(ninter, tmwork, rbigst, rsecnd,
     *                 sumrhs, intref)
         tstval = bigfac*(sumrhs-rbigst)/ninter
         if (rbigst .ge. small .and. rbigst .ge. tstval) go to 100
   50 continue

*  If we reach this point, no interval has a significantly larger
*  element than average.  Set counters and double the mesh.

      numbig = 0
      nummed = 0
      ddouble = .true.

cf the mesh if not doubled if the problem is stiff and the order is 4
       if (nodouble .and. .not. forcedouble) then
          call selcondmsh(ncomp, nmsh,
     *     nfxpnt, fixpnt,  nmax, xx,  irefin,
     *     nmold, xxold, ddouble , maxmsh,r4,amg)
           ddouble = .false.
           itcond = itcond + 1
        else
          call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
            itcond = 0
        end if
c      call dblmsh(nmsh, nmax, xx, nmold, xxold, maxmsh)
      return

*  To reach statement 100, it must be true that for some component
*  icmp,  rbigst is non-negligible and large relative to the
*  average.  intref indicates the region to which points may be added.

*  If too many specialized refinements (adding a few points to
*  a small number of intervals) have been made, signal that
*  the mesh should be doubled.
*  Otherwise, increment counters and add numadd points as indicated.

  100 continue
      if (rbigst .lt. two*rsecnd) nummed = nummed + 1
      numadd = numpt
      numbig = numbig + 1
      ddouble = .false.

      if (rbigst .le. bigfac*rsecnd .or. numbig .gt. 8) then
         numbig = 0
         nummed = nummed + 1
         if (nummed .ge. 4 .and. iorder .eq. 4) then
            ddouble = .true.
            nummed = 0
         elseif (nummed .ge. 8 .and. iorder .gt. 4) then
            ddouble = .true.
            nummed = 0
         endif
      endif

*  Refine the mesh.
      if (ddouble) then
cf the mesh if not doubled if the problem is stiff
         if  (nodouble .and. .not. forcedouble)  then
             numadd = numpt
c            call smpselcondmsh(ncomp, nmsh,
c     *       nfxpnt, fixpnt,  nmax, xx,  irefin,intref,numadd,
c     *       nmold, xxold, ddouble , maxmsh,r4,amg)
c               itcond = itcond + 1
           call selcondmsh(ncomp, nmsh,
     *        nfxpnt, fixpnt,  nmax, xx,  irefin,
     *        nmold, xxold, ddouble , maxmsh,r4,amg)
           ddouble = .false.
           itcond = itcond + 1
         else
           call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
             itcond  = 0
         end if
      else
cf if the problem is stiff  we use both the old technique
cf and the conditioning
cf
         if (nodouble .and. .not. forcedouble)  then
           numadd = numpt
           call smpselcondmsh(ncomp, nmsh,
     *       nfxpnt, fixpnt,  nmax, xx,  irefin,intref,numadd,
     *       nmold, xxold, ddouble , maxmsh,r4,amg)
               itcond = itcond + 1
        elseif (forcedouble .and. use_c) then
            ddouble = .true.
              itcond = 0
            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
        else
                numadd = numpt
                call smpmsh (nmsh, nmax, xx, intref, numadd,
     *          nmold, xxold, maxmsh)
        endif

      endif

      return
      end

c ===================================================================================

      subroutine errest (ncomp, nmsh, ntol, ltol, tol,
     *   nudim, u, uold, etest, errmax, errok)

      implicit double precision (a-h,o-z)
      dimension ltol(ntol), tol(ntol), u(nudim,nmsh),
     *              uold(ncomp,nmsh), etest(ntol)
      logical errok

      logical use_c, comp_c
      common/algprs/ nminit, iprint,idum, use_c, comp_c
      intrinsic abs
      intrinsic max

      parameter ( one = 1.0d+0, zero = 0.0d+0 )
      character(len=150) msg


*  Given current and previous solutions u and uold on the same
*  mesh, errest calculates an error measure for each
*  component for which a tolerance is specified.
*  The error measure is the usual relative error normalized
*  by dividing by the tolerance.  On exit, errsum is the
*  sum of these error measures.

*  The array etest specifies the error test to be applied to each
*  error measure.

*  On exit, the logical flag errok
*   -- is .false. if any of the error measures exceeds the
*      corresponding value in the array etest
*   -- is .true. if all error measures are less than the
*      corresponding values of etest.

cf      errsum = zero
      errmax = zero
      errok = .true.

      do 10 im = 1, nmsh
      do 10 it = 1, ntol
         icmp = ltol(it)
         er = u(icmp,im) - uold(icmp,im)
         denom = max(one, abs(uold(icmp,im)))
         errel = abs(er/(tol(it)*denom))
         errmax = max(errmax,  errel)
cf         errsum = errsum + errel
         if (errel .gt. etest(it)) errok = .false.
   10 continue

      return
      end

c ===================================================================================

      subroutine getptq( debug, mfsrch, nout, alfmax, alfsml, alfuzz,
     *                   epsaf, epsag, eta, ftry, oldf, oldg,
     *                   rmu, tolabs, tolrel, toltny, imprvd,
     *                   inform, nfsrch, alfa, alfbst, fbest,
     *                   braktd, crampd, extrap,vset,wset,nsamea,nsameb,
     *                   a, b, fa, factor, fv, fw, xtry, xv, xw )

      implicit double precision (a-h,o-z)
      logical            debug, imprvd
      logical            braktd, crampd, extrap, vset, wset
      integer            mfsrch, nout, inform, nfsrch, nsamea, nsameb
      character(len=150) msg
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
c     if (debug) THEN
c  write (msg, 1000) alfmax, oldf, oldg, tolabs,
c    *   alfuzz, epsaf, epsag, tolrel, crampd
c       CALL rprint(msg)
c     ENDIF

      go to 800
c
c  ---------------------------------------------------------------------
c  subsequent entries.
c  the function has just been evaluated at  alfa = alfbst + xtry,
c  giving  ftry.
c  ---------------------------------------------------------------------
  100 nsamea = nsamea + 1
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
c     if (debug) THEN
c  write (msg, 1100) alfa, ftry, ctry
c       CALL rprint(msg)
c     ENDIF

      if (.not. imprvd) go to 130
c
c  we seem to have an improvement.  the new point becomes the
c  origin and other points are shifted accordingly.
c
      if (.not. wset) go to 110
      xv     = xw - xtry
      fv     = fw
      vset   = .true.
  110 xw     = zero - xtry
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
  120 b      = xw
      nsameb = 0
      braktd = .true.
      go to 300
c
c  the new function value is no better than the best point found so far.
c  the point  xtry  must be a new end point of the interval of
c  uncertainty.
c
  130 if (xtry .ge. zero) go to 140
      a      = xtry
      fa     = ftry
      nsamea = 0
      go to 150
  140 b      = xtry
      nsameb = 0
      braktd = .true.
c
c  the origin remains unchanged but  xtry  may qualify as  xw.
c
  150 if (.not. wset)   go to 160
      if ((ftry - fw) .gt. epsaf) go to 170
      xv     = xw
      fv     = fw
      vset   = .true.
  160 xw     = xtry
      fw     = ftry
      wset   = .true.
      if (moved) extrap = xinxw
      go to 300
c
c  ftry  is no better than  fbest  or  fw.  if the best point has not
c  been moved, there must be more than one minimum.
c
  170 if (moved) go to 175
      xw     = xtry
      fw     = ftry
      go to 300
c
c  ftry  is no better than  fbest  or  fw,  but  xtry  may become  xv.
c  extrap  has the value set in the previous entry.
c
  175 if (.not. vset) go to 180
      if ((ftry - fv) .gt. epsaf  .and.  extrap) go to 300
  180 if (xinxw) go to 190
      xv     = xtry
      fv     = ftry
      vset   = .true.
      go to 300
  190 if (vset) xw = xv
      if (vset) fw = fv
      xv     = xtry
      fv     = ftry
      vset   = .true.
c
c  ---------------------------------------------------------------------
c  check the termination criteria.
c  ---------------------------------------------------------------------
  300 tol    = tolrel*alfbst + tolabs
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
c     if (debug) write (nout, 1200) atrue, btrue, gap, tol,
c    *   nsamea, nsameb, braktd, closef, imprvd, conv1, conv2, conv3,
c    *   extrap, alfbst, fbest, cbest, alfaw, fw
      if (vset) alfav  = alfbst + xv
c     if (debug  .and.  vset) THEN
c     write (msg, 1300) alfav, fv
c       CALL rprint(msg)
c     ENDIF

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
  400 xmidpt = half*(a + b)
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
c     if (debug) THEN
c  write (msg, 2100)
c       CALL rprint(msg)
c     ENDIF

      go to 600
c
c  three points available.  use  fbest,  fw  and  fv.
c
  450 gv = (fv - fbest)/xv
      s  = gv - (xv/xw)*gw
      q  = two*(gv - gw)
c     if (debug) THEN
c  write (msg, 2200)
c       CALL rprint(msg)
c     ENDIF

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
  600 artifa = a
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
  610 if (vset  .and.  moved) go to 620
      xtry   = xw/ten
c     if (debug) THEN
c  write (msg, 2400) xtry
c       CALL rprint(msg)
c     ENDIF

      go to 700
c
c  three points exist in the interval of uncertainty.  check whether
c  the points are configured for an extrapolation or interpolation.
c
  620 if (extrap) go to 660
c
c  if the interpolation appears to be consistently over-estimating the
c  distance to the minimum,  damp the interpolation step.
c
      if (nsamea .lt. 3  .and.  nsameb .lt. 3) go to 630
      factor = factor / five
      s      = factor * s
      go to 640
  630 factor = one
c
c  the points are configured for an interpolation.  the artificial
c  interval will be just  (a,b).   set  endpnt  so that  xtry
c  lies in the larger of the intervals  (a,0)  and  (0,b).
c
  640 if (xmidpt .lt. zero) endpnt = a
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
  660 if (xw .lt. zero) endpnt = b
      if (xw .gt. zero) endpnt = a
c
c  compute the default value of  xtry.
c
  680 dtry = abs( endpnt )
      daux = gap - dtry
      if (daux .ge. dtry)   xtry = five*dtry*(point1 + dtry/daux)/eleven
      if (daux .lt. dtry)   xtry = half*sqrt( daux )*sqrt( dtry )
      if (endpnt .lt. zero) xtry = - xtry
c     if (debug) THEN
c  write (msg, 2500) xtry, daux, dtry
c       CALL rprint(msg)
c     ENDIF

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
  700 if (q .eq. zero) go to 800
      if (q .lt. zero) s = - s
      if (q .lt. zero) q = - q
      if (s*xw .lt. q*artifa   .or.   s*xw .gt. q*artifb) go to 800
c
c  accept the polynomial fit.
c
      xtry = zero
      if (abs( s*xw ) .ge. q*tol) xtry = (s/q)*xw
c     if (debug) THEN
c  write (msg, 2600) xtry
c       CALL rprint(msg)
c     ENDIF

c
c  ---------------------------------------------------------------------
c  test for  xtry  being larger than  alfmax  or too close to  a  or  b.
c  ---------------------------------------------------------------------
  800 if (braktd) go to 810
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
  810 xmidpt = half*(a + b)
      if (xtry .gt. a + tol  .and.  xtry .lt. b - tol) go to 820
      if (xmidpt .gt. zero) xtry =   tol
      if (xmidpt .le. zero) xtry = - tol
c
c
c  f  must not be calculated too close to  alfbst.
c
  820 if (abs( xtry ) .lt. tol  .and.  xmidpt .lt. zero) xtry = - tol
      if (abs( xtry ) .lt. tol  .and.  xmidpt .ge. zero) xtry =   tol
      alfa   = alfbst + xtry
c
c  ---------------------------------------------------------------------
c  exit.
c  ---------------------------------------------------------------------
c
c  new function value required.
c
  900 inform = 0
      go to 990
c
c  convergence test satisfied.
c
  910 inform = 1
      if (alfa .eq. alfmax) inform = 2
      go to 990
c
c  mfsrch  function evaluations without sufficient decrease, but an
c  improved point was found.
c
  930 if (.not. moved) go to 960
      inform = 3
      go to 990
c
c  zero step (alfmax too small).
c
  940 inform = 4
      go to 990
c
c  premature termination.  the step is smaller than  alfsml.
c
  950 inform = 5
      go to 990
c
c  zero step (a sufficiently better point could not be found).
c
  960 inform = 6
      go to 990
c
c  zero step (positive gradient at the starting point).
c
  970 inform = 7
c
c  exit.
c
  990 continue
c 990 if (debug) THEN
c       write (msg, 3000)
c       CALL rprint(msg)
c     ENDIF

      return
c KSKSKS:
c1000 format(/ 31h alfmax  oldf    oldg    tolabs, 1p2e22.14, 1p2e16.8
c    *       / 31h alfuzz  epsaf   epsag   tolrel, 1p2e22.14, 1p2e16.8
c    *       / 31h crampd                        ,  l6)
c1100 format(/ 31h alfa    ftry    ctry          , 1p2e22.14, 1pe16.8)
c1200 format(/ 31h a       b       b - a   tol   , 1p2e22.14, 1p2e16.8
c    *       / 31h nsamea  nsameb  braktd  closef, 2i3, 2l6
c    *       / 31h imprvd  convrg  extrap        ,  l6, 3x, 3l1, l6
c    *       / 31h alfbst  fbest   cbest         , 1p2e22.14, 1pe16.8
c    *       / 31h alfaw   fw                    , 1p2e22.14)
c1300 format(  31h alfav   fv                    , 1p2e22.14 /)

c2100 format(30h parabolic fit,    two points.)
c2200 format(30h parabolic fit,  three points.)
c2400 format(31h exponent reduced.  trial point, 1p1e22.14)
c2500 format(31h geo. bisection. xtry,daux,dtry, 1p3e22.14)
c2600 format(31h polynomial fit accepted.  xtry, 1p1e22.14)
c3000 format(53h ---------------------------------------------------- /)
c
c  end of getptq
      end

c ===================================================================================

      subroutine interp(ncomp, nmsh, xx, nudim, u,nuold,
     *                    nmold, xxold, uold)

      implicit double precision (a-h, o-z)
      dimension xx(*), u(nudim,*), xxold(*), uold(nuold,*)
      logical  use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c

* blas: dcopy

      parameter (zero = 0.0d+0)
      character(len=150) msg

*  interp performs piecewise linear interpolation of the old
*  solution uold at the nmold old mesh points xxold onto the nmsh
*  new mesh points xx, producing an interpolated solution u.
*  Note that no assumption is made that the new mesh has
*  more points than the old, nor that the new and old mesh
*  points are related in a specific way (except that their first
*  and last points are identical).

*  By construction, xx(1) = xxold(1).  Copy the first ncomp
*  components of uold into those of u.

      call dcopy(ncomp, uold(1,1), 1, u(1,1), 1)

      i = 2
      do 100 im = 2, nmsh-1

   50    continue
         if (i .gt. nmold) return

*  Check whether the im-th point in the new mesh lies strictly
*  to the right of, or to the left of (or exactly on) the
*  i-th point in the old mesh.


         if (xx(im) .gt. xxold(i)) then
            i = i + 1
            go to 50
         else
            xdif = xxold(i) - xx(im)
            if (xdif .eq. zero) then

*  xx(im) and xxold(i) are identical.

               call dcopy(ncomp, uold(1,i), 1, u(1,im), 1)
               i = i + 1
            else
               xint = xxold(i) - xxold(i-1)
               xrat = xdif/xint
               do 70 k = 1, ncomp
                  u(k,im) = uold(k,i) + xrat*(uold(k,i-1)-uold(k,i))
   70          continue
            endif
         endif

  100 continue
      call dcopy(ncomp, uold(1,nmold), 1, u(1,nmsh), 1)
      return
      end

c ===================================================================================

      subroutine rerrvl( ncomp, nmsh, nudim, u, usvrex, ntol, ltol,
     *          rerr, remax, itlmx, adjrer )
      implicit double precision (a-h,o-z)
      dimension ltol(*)
      dimension u(nudim, *), usvrex(ncomp, *)
      dimension rerr(ncomp, *)
      logical adjrer

      intrinsic abs
      intrinsic max

      parameter ( one = 1.0d+0, zero = 0.0d+0 )
      character(len=150) msg

*  rerrvl is used in considering Richardson extrapolation.
*  The two solutions u and usvrex have a special relationship:
*  u corresponds to a doubled mesh, with twice as many
*  intervals, where each interval is half the size of that for
*  usvrex's mesh.   nmsh is the number of mesh points in the
*  mesh for u.

*  remax is the maximum relative error, and itlmx is the
*  index of the tolerance for which the maximum occurs.

*  The array rerr contains the absolute value of the difference
*  between u and usvrex at the common mesh points, but is defined
*  only for components for which an error tolerance is specified.

*  The logical variable adjrer is true on entry if the values in
*  rerr are to be adjusted for later use in selective mesh refinement.


      itlmx = 1
      remax = zero
      nmold = 1 + (nmsh-1)/2
      do 100 it = 1, ntol
         icmp = ltol(it)
         imnew = 1
         do 50 im = 1, nmold
            rerr(icmp, im) = abs(usvrex(icmp,im) - u(icmp,imnew))
            denom = max(one, abs(usvrex(icmp,im)))
            rerel = rerr(icmp,im)/denom
            if (rerel .gt. remax) then
               remax = rerel
               itlmx = it
            endif
            imnew = imnew + 2
   50    continue
  100 continue

      if (adjrer) then

*  Adjust the rerr array if it may be used later in selective
*  mesh refinement.

         do 150 it = 1, ntol
            icmp = ltol(it)
            do 150 im = 1, nmold - 1
               rerr(icmp,im) = max(rerr(icmp,im),
     *                              rerr(icmp, im+1))
  150    continue
      endif

      return
      end

c ===================================================================================

      double precision function dasum(n,dx,incx)
c#
c#     takes the sum of the absolute values.
c#     jack dongarra, linpack, 3/11/78.
c#     modified 3/93 to return if incx .le. 0.
c#     modified 12/3/93, array(1) declarations changed to array(*)
c#
      double precision dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
      character(len=150) msg
c#
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n .le. 0 .or. incx .le. 0 )return
      if(incx.eq.1)go to 20
c#
c#        code for increment not equal to 1
c#
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
c#
c#        code for increment equal to 1
c#
c#
c#        clean-up loop
c#
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end

c ===================================================================================

      subroutine selconderrmsh(ncomp, nmsh, ntol, ltol, tol,
     *     nfxpnt, fixpnt, ipow, nmax,
     *     xx, nudim, u, ermeas, irefin, ihcomp, nmold, xxold,
     *     ermx, ddouble , maxmsh,r4,amg,stab_cond)

      implicit double precision (a-h,o-z)

      dimension  ltol(ntol), tol(ntol), fixpnt(*)
      dimension  xx(*), u(nudim, *), ermeas(ncomp,*)
      dimension  irefin(nmsh-1), ihcomp(nmsh-1)
      dimension  xxold(*), ermx(*), amg(*), r4(nmsh)
      logical    ddouble , maxmsh, stab_cond

      intrinsic abs
      intrinsic max
      intrinsic int

*  blas: dcopy
*  double precision dlog

      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c

      parameter  ( zero = 0.0d+0, one = 1.0d+0, onep1 = 1.1d+0 )
      parameter  ( erdcid = 5.0d+0 )
      parameter  ( phitst = 0.1d+0 )

      logical first, add
      save    first, rlndec
      data    first / .true. /
      character(len=150) msg

*  The routine selconderrmsh performs selective mesh refinement, depending
*  on the error measure ermeas and the monitor function based on the
*  conditioning.

      if (first) then
         first = .false.
         rlndec = dlog(erdcid)
      endif

      maxmsh = .false.

      frcpow = one/ipow
      ddouble = .false.
      nmold = nmsh
      ninter = nmsh - 1

*  Copy the current mesh into the xxold array.


      call dcopy(nmold, xx, 1, xxold, 1)
      ithres = 0
      thres = one

*  On input, the array ermeas represents some error measure defined
*  over the components and mesh intervals (not mesh points).
*  It is normalized in the following loop with respect to the
*  tolerance array and the current solution.
*  The value errmax gives the maximum normalized error.


      errmax = zero
      do 120 im = 1, ninter
         ermx(im) = zero
         do 110 it = 1, ntol
            jcomp = ltol(it)
            denom = tol(it)*max(one, abs(u(jcomp,im)))
            ems = ermeas(jcomp,im)
            ermeas(jcomp,im) = abs(ems)/denom
            err = ermeas(jcomp, im)
            if (err .ge. ermx(im)) then
                ermx(im) = err
                ihcomp(im) = jcomp
            endif
  110    continue
         errmax = max(ermx(im), errmax)
  120 continue

c      write(*,*) 'tol', tol(1),tol(2)

      call moncondmsh(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,nptcond,r4,amg)

c      write(*,*) 'errmax', errmax, 'nptcond', nptcond
c
       if (.not. stab_cond .and. errmax .ge. 1e20 ) then
cf     *           .and. r1 .gt. 1.0d0 ) then
cf  only the conditioning
cf
            call  selcondmsh(ncomp, nmsh,
     *         nfxpnt, fixpnt,  nmax, xx,  irefin,
     *         nmold, xxold, ddouble , maxmsh,r4,amg)
      else
cf   the conditioning and the error

      if (errmax .gt. zero .and. errmax .le. erdcid) then

*  If errmax > 0 and .le. erdcid, find the smallest integer exponent ii
*  such that (erdcid**ii)*errmax > erdcid.

         if(errmax .gt. one) then
            ii = 1
            decii = erdcid
         else
            ilg = -dlog(errmax)/rlndec
            ii = 2 + ilg
            decii = erdcid**ii
         endif

*  Multiply error measures by erdcid**ii.

         errmax = decii*errmax
         do 140 im = 1, ninter
            ermx(im) = decii*ermx(im)
            do 140 it = 1, ntol
               jcomp = ltol(it)
               ermeas(jcomp,im) = decii*ermeas(jcomp,im)
  140    continue
      endif

  200 continue

*  For each interval im,  the integer irefin(im) is calculated
*  based on an equidistrbution procedure involving the
*  threshold value thres.  If irefin(im) > 1, we add
*  points to interval im.  If irefin(im) = 1, we check whether
*  point im can be removed from the mesh.

*  nmest is a lower bound on the number of points in the new mesh.
*  We do not know in advance how many points will be removed,
*  so nmest is computed by assuming all eligible points are removed.


      nmest = nmsh
      do 220 im = 1, ninter
         if (ermx(im).ge.thres) then
            irefin(im) = int(ermx(im)**frcpow) + 1
            nmest = nmest + irefin(im) - 1
         else
            irefin(im) = 1
            nmest = nmest - 1
         endif
  220 continue

c       if(nptcond.ge.4) then
c         add = .false.
c         nptcond = nptcond/2
c       do 221 im = 1, ninter-1
c         if  (max(r4(im), r4(im+1)) .gt. fatt_r1r3) then
c            if (.not. add) then
c              irefin(im) = max(nptcond, irefin(im))
c              nmest = nmest + nptcond - 1
c            endif
c              irefin(im+1) = max(nptcond, irefin(im+1))
c                nmest = nmest + nptcond - 1
c              add = .true.
c         else
c            irefin(im) = max(1, irefin(im))
c            irefin(im+1) = max(1, irefin(im+1))
c            nmest = nmest - 1
c            add = .false.
c         endif
c 221     continue

c        do 221 im = 1, ninter
c         if ( r4(im) .gt. fatt_r1r3) then
c              irefin(im) = max(nptcond, irefin(im))
cc              nmest = nmest + nptcond - 1
c         else
c              irefin(im) = max(1, irefin(im))
cc              nmest = nmest - 1
c         endif
c 221     continue
c        end if

      if (nmest .gt. nmax) then

         go to 360

      elseif (nmest-1 .gt. 3*ninter) then

         call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
         ddouble = .true.
         return
      endif

*  It appears that we can perform the desired selective mesh
*  refinement.

*  Now begin running through the mesh, adding and possibly deleting
*  points as indicated by the irefin array.

*  The integer new is a count of the number of intervals in
*  the tentative mesh being generated by the refinement strategy.

      new = 1

*  The first interval is treated as a special case, since xx(1)
*  always remains in the mesh, and cannot be a fixed point.

      rlen = xxold(2) - xx(1)
      slen = rlen
      if (new + irefin(1)  .gt. nmax) goto 360
      if (irefin(1).gt.1) then
         dx = rlen/irefin(1)
         do 230 j = 2, irefin(1)
            new = new + 1
            xx(new) = xx(1) + (j-1)*dx
  230    continue
      endif

*  The fixed points specified by the fixpnt array cannot be
*  removed from the mesh.  The value fxnext indicates the 'next'
*  fixed point.  When no further fixed points remain to be processed
*  (or if nfxpnt = 0), fxnext is set to a value strictly larger than
*  the last mesh point, so that no mesh point can equal fxnext.
*  This way we need to compare only the values of xxold(i)
*  and fxnext.

      ifxcnt = 1
      if (nfxpnt .eq. 0) then
         fxnext = onep1*abs(xxold(nmsh))
      else
         fxnext = fixpnt(ifxcnt)
      endif

*  jtkout is a counter of the number of consecutive points that
*  have been removed from the mesh.

      jtkout = 0
      do 330 im = 2, ninter
         rlold = rlen
         rlen = xxold(im+1) - xxold(im)

*  If xxold(im) is the next fixed point, it cannot be removed
*  and so we don't test its error estimates.

         if(xxold(im) .eq. fxnext)  then

            ifxcnt = ifxcnt + 1
            if(ifxcnt .gt. nfxpnt) then
               fxnext = onep1*abs(xxold(nmsh))
            else
               fxnext = fixpnt(ifxcnt)
            endif

         elseif (irefin(im).eq.1) then

*  If xxold(im) is not a fixed point and irefin(im) = 1, we may wish
*  to remove point im from the mesh.

*  If we are considering removing points and jtkout = 0, this
*  is the first point in a possible consecutive set to be removed,
*  and we initialize phihat, which represents a maximum of
*  certain estimates.
*  If jtkout is not zero, previous points contiguous to this
*  point have been removed, and phihat retains its previous value.

            slen = slen + rlen

            if (jtkout .eq. 0) then
                ind1 = ihcomp(im-1)
                phihat = ermeas(ind1,im-1)/(rlold**ipow)
            endif
            phihat = max(phihat,
     *                 ermeas(ihcomp(im),im)/(rlen**ipow))
            val1 = phihat*(slen**ipow)
            if (val1 .le. phitst
     *             .and. jtkout .lt. 4
     *             .and. r4(im) .lt. fatt_r1r3) then

*  Increment the counter of removed points.
*  'Remove' the mesh point xxold(im) by not including it.

               jtkout = jtkout+1
               go to 330
            endif
*        end of logic for irefin(im) = 1.
         endif

         jtkout = 0
         new = new + 1
         if (new  .gt. nmax) goto 360
         xx(new) = xxold(im)
         if (new + irefin(im)  .gt. nmax) goto 360
         if (irefin(im) .gt. 1) then
            dx = rlen/irefin(im)
            do 300 j = 2, irefin(im)
              new = new + 1
              xx(new) = xxold(im) + (j-1)*dx
  300       continue
         endif
         slen = rlen

         if (new .gt. nmax) then

*  If the new mesh contains too many points, branch out of the
*  loop to try alternative strategies.

            go to 360

         elseif (new .gt. 3*ninter)  then

*  Here, the new mesh does not exceed the specified maximum,
*  but has more than 3 times as many intervals as the old mesh.
*  Try doubling the mesh if possible.
            nmsh = nmold
            call dcopy(nmsh, xxold, 1, xx, 1)
            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
            ddouble = .true.
            return

         endif

  330 continue

*  To end up here, we have processed the entire interval,
*  and have neither exceeded the specified maximum nor
*  exceeded three times the number of intervals in the old
*  mesh.  The last mesh point remains unchanged.

      new = new + 1
      if (new   .gt. nmax) goto 360
      xx(new) = xxold(nmsh)
      nmsh = new
      maxmsh = .false.
      if (iprint .ge. 0) THEN
        write(msg,905) nmsh
        CALL rprint(msg)
      ENDIF

      return

  360 continue
*  To reach here, the number of mesh points created at some stage
*  of the refinement process was larger than the maximum permitted
*  value nmax.

*  Check whether the mesh can safely be doubled.

      if ((2*nmsh-1) .lt. nmax) then

*  Double the mesh.
         call  dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
         ddouble = .true.

*  If the number of intervals is too large and the mesh cannot be
*  doubled, increase the threshold thres by a factor of erdcid and
*  try the selective refinement again.
*  If this happens three times without success or if thres exceeds
*  or is equal to errmax, stop.  (In this case, we know already
*  that doubling the mesh produces too many points.)

      elseif (thres .lt. errmax .and. ithres .lt. 3) then
         ithres = ithres + 1
         thres = erdcid*thres
         if(thres .gt. errmax) thres = errmax
         call dcopy(nmsh, xxold, 1, xx, 1)
         go to 200
      else
c         nmsh = 2*nmsh - 1
         nmsh = nmold
         call dcopy(nmsh, xxold, 1, xx, 1)
         maxmsh = .true.
      endif

      end if
c   endif use the conditioning and the error
      return

  905 format(1h ,'selconderrmsh.  new nmsh =',i8)
  910 format(1h ,'ihcomp',(10i5))
      end

c ===================================================================================

      subroutine selcondmsh(ncomp, nmsh,
     *     nfxpnt, fixpnt,  nmax, xx,  irefin,
     *     nmold, xxold, ddouble , maxmsh, r4, amg)

      implicit double precision (a-h,o-z)

      dimension  fixpnt(*)
      dimension  xx(*)
      dimension  irefin(*)
      dimension  xxold(*),  amg(*), r4(nmsh)
      logical    ddouble , maxmsh, add

      intrinsic abs
      intrinsic max
      intrinsic int

*  blas: dcopy
*  double precision dlog

      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c

      parameter  ( zero = 0.0d+0, one = 1.0d+0, onep1 = 1.1d+0 )
      parameter  ( erdcid = 5.0d+0 )
      character(len=150) msg




*  The routine selcondmsh performs selective mesh refinement, depending
*  on the monitor function based on the conditioning parameters.


      maxmsh = .false.
      ddouble = .false.
      nmold = nmsh
      ninter = nmsh - 1

*  Copy the current mesh into the xxold array.


      call dcopy(nmold, xx, 1, xxold, 1)

      ithres = 0
      thres = one

*  On input, the array amg represents the conditioning vector defined
*  over the components and mesh intervals (not mesh points).
*  We compute the monitor function and the related parameters in the
*  followinf function

      call moncondmsh(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,nptcond,r4,amg)




*  For each interval im,  the integer irefin(im) is calculated
*  based on an equidistrbution procedure involving the
*  threshold value thres.  If irefin(im) > 1, we add
*  points to interval im.  If irefin(im) = 1, we check whether
*  point im can be removed from the mesh.

*  nmest is a lower bound on the number of points in the new mesh.
*  We do not know in advance how many points will be removed,
*  so nmest is computed by assuming all eligible points are removed.

      add = .false.
      nmest = nmsh
c      do 220 im = 1, ninter-1
c         if  (max(r4(im), r4(im+1)) .gt. fatt_r1r3) then
c            if (.not. add) then
c              irefin(im) = nptcond
c              nmest = nmest + nptcond - 1
c            endif
c              irefin(im+1) = nptcond
c                nmest = nmest + nptcond - 1
c              add = .true.
c         else
c            irefin(im) = 1
c            irefin(im+1) = 1
c            nmest = nmest - 1
c            add = .false.
c         endif

c  220 continue

      do 220 im = 1, ninter
         if ( r4(im) .gt. fatt_r1r3) then
              irefin(im) = nptcond
              nmest = nmest + nptcond - 1
         else
              irefin(im) = 1
              nmest = nmest - 1
         endif

  220 continue


      if (nmest .gt. nmax) then

         go to 360

      endif

*  It appears that we can perform the desired selective mesh
*  refinement.

*  Now begin running through the mesh, adding and possibly deleting
*  points as indicated by the irefin array.

*  The integer new is a count of the number of intervals in
*  the tentative mesh being generated by the refinement strategy.

      new = 1

*  The first interval is treated as a special case, since xx(1)
*  always remains in the mesh, and cannot be a fixed point.

      rlen = xxold(2) - xx(1)
      slen = rlen
      if (new + irefin(1)  .gt. nmax) goto 360
      if (irefin(1).gt.1) then
         dx = rlen/irefin(1)
         do 230 j = 2, irefin(1)
            new = new + 1
            xx(new) = xx(1) + (j-1)*dx
  230    continue
      endif

*  The fixed points specified by the fixpnt array cannot be
*  removed from the mesh.  The value fxnext indicates the 'next'
*  fixed point.  When no further fixed points remain to be processed
*  (or if nfxpnt = 0), fxnext is set to a value strictly larger than
*  the last mesh point, so that no mesh point can equal fxnext.
*  This way we need to compare only the values of xxold(i)
*  and fxnext.

      ifxcnt = 1
      if (nfxpnt .eq. 0) then
         fxnext = onep1*abs(xxold(nmsh))
      else
         fxnext = fixpnt(ifxcnt)
      endif

*  jtkout is a counter of the number of consecutive points that
*  have been removed from the mesh.

      jtkout = 0
      do 330 im = 2, ninter

         rlold = rlen
         rlen = xxold(im+1) - xxold(im)

*  If xxold(im) is the next fixed point, it cannot be removed
*  and so we don't test its error estimates.

         if(xxold(im) .eq. fxnext)  then

            ifxcnt = ifxcnt + 1
            if(ifxcnt .gt. nfxpnt) then
               fxnext = onep1*abs(xxold(nmsh))
            else
               fxnext = fixpnt(ifxcnt)
            endif

         elseif (irefin(im).eq.1) then

*  If xxold(im) is not a fixed point and irefin(im) = 1, we may wish
*  to remove point im from the mesh.

*  If we are considering removing points and jtkout = 0, this
*  is the first point in a possible consecutive set to be removed,
*  and we initialize phihat, which represents a maximum of
*  certain estimates.
*  If jtkout is not zero, previous points contiguous to this
*  point have been removed, and phihat retains its previous value.

            slen = slen + rlen

            if ( jtkout .lt. 2
     *             .and. r4(im) .le. fatt_r3 ) then

*  Increment the counter of removed points.
*  'Remove' the mesh point xxold(im) by not including it.

               jtkout = jtkout+1
               go to 330
            endif
*        end of logic for irefin(im) = 1.
         endif
         if (new+irefin(im) .gt. nmax) then

*  If the new mesh contains too many points, branch out of the
*  loop to try alternative strategies.

            go to 360
         end if

         jtkout = 0
         new = new + 1
         xx(new) = xxold(im)
         if (irefin(im) .gt. 1) then
            dx = rlen/irefin(im)
            do 300 j = 2, irefin(im)
              new = new + 1
              xx(new) = xxold(im) + (j-1)*dx
  300       continue
         endif
         slen = rlen

         if (new .gt. nmax) then

*  If the new mesh contains too many points, branch out of the
*  loop to try alternative strategies.

            go to 360

         elseif (new .gt. 3*ninter)  then

*  Here, the new mesh does not exceed the specified maximum,
*  but has more than 3 times as many intervals as the old mesh.
*  Try doubling the mesh if possible.
            nmsh = nmold
            call dcopy(nmsh, xxold, 1, xx, 1)
            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
            ddouble = .true.
            return


         endif



  330 continue

*  To end up here, we have processed the entire interval,
*  and have neither exceeded the specified maximum nor
*  exceeded three times the number of intervals in the old
*  mesh.  The last mesh point remains unchanged.

      new = new + 1
      if (new   .gt. nmax) goto 360
      xx(new) = xxold(nmsh)
      nmsh = new
      maxmsh = .false.
      if (iprint .ge. 0) THEN
        write(msg,905) nmsh
        CALL rprint(msg)
      ENDIF

      return

  360 continue
*  To reach here, the number of mesh points created at some stage
*  of the refinement process was larger than the maximum permitted
*  value nmax.

      nmsh = nmold
      call dcopy(nmsh, xxold, 1, xx, 1)
      maxmsh = .true.


      return



  902 format(1h ,'im, jcomp, ermeas, normalized er',2i5,2(1pe11.3))
  903 format(1h ,'errmax',1pe11.3)
  904 format(1h ,'nmest, irefin',(10i5))
  905 format(1h ,'selcondmsh.  new nmsh =',i8)
  910 format(1h ,'ihcomp',(10i5))
      end

c ===================================================================================

      subroutine smpselcondmsh(ncomp, nmsh,
     *     nfxpnt, fixpnt,  nmax, xx,  irefin, intref, numadd,
     *     nmold, xxold, ddouble , maxmsh, r4, amg)

      implicit double precision (a-h,o-z)

      dimension  fixpnt(*)
      dimension  xx(*)
      dimension  irefin(*)
      dimension  xxold(*),  amg(*), r4(*)
      logical    ddouble , maxmsh, add

      intrinsic abs
      intrinsic max
      intrinsic int

*  blas: dcopy
*  double precision dlog

      logical use_c, comp_c
      common/algprs/ nminit, iprint, idum, use_c, comp_c

      parameter  ( zero = 0.0d+0, one = 1.0d+0,onep1 = 1.1d+0)
      parameter  ( erdcid = 5.0d+0 )
      character(len=150) msg




*  The routine smpselcondmsh performs selective mesh refinement, by adding
*  points to one or three interval(s) in the region indicated
*  by the integer intref (numadd gives the trial number of points to
*  added in each  interval)  and by adding and removing point
*  using  the monitor function based on the conditioning parameters.



      maxmsh = .false.

      ddouble = .false.
      nmold = nmsh
      ninter = nmsh - 1

*  Copy the current mesh into the xxold array.


         call dcopy(nmold, xx, 1, xxold, 1)

      ithres = 0
      thres = one

*  On input, the array amg represents the conditioning vector defined
*  over the components and mesh intervals (not mesh points).

      call moncondmsh(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,nptcond,r4,amg)



*  For each interval im,  the integer irefin(im) is calculated
*  based on an equidistrbution procedure involving the
*  threshold value thres.  If irefin(im) > 1, we add
*  points to interval im.  If irefin(im) = 1, we check whether
*  point im can be removed from the mesh.

*  nmest is a lower bound on the number of points in the new mesh.
*  We do not know in advance how many points will be removed,
*  so nmest is computed by assuming all eligible points are removed.


      nmest = nmsh

      do  im = 1,ninter
         irefin(im)=1
      enddo
      do 220 im = 1, ninter
         if  (r4(im) .ge. fatt_r1r3) then
              irefin(im) = nptcond
              nmest = nmest + nptcond - 1
         endif

  220 continue

      if(numadd .gt. 49) then
         numadd = 49
      elseif (numadd .lt. 4) then
         numadd = 4
      endif

        if ( intref .eq. 1) then
            irefin(1) = max(numadd,irefin(1))
            nmest = nmest + irefin(1) - 1
        elseif (intref .eq. ninter) then
            irefin(ninter) = max(numadd,irefin(ninter))
            nmest = nmest + irefin(ninter) - 1
        else
          if(numadd .gt. 9) then
             numadd = 9
           elseif (numadd .lt. 4) then
              numadd = 4
            endif
            irefin(intref-1) = max(numadd , irefin(intref-1))
            irefin(intref)   = max(numadd , irefin(intref))
            irefin(intref+1) = max(numadd , irefin(intref+1))
            nmest = nmest +
     *            irefin(intref-1)+irefin(intref)+irefin(intref+1) - 1
        end if



      if (nmest .gt. nmax) then

         go to 360

      endif


*  It appears that we can perform the desired selective mesh
*  refinement.

*  Now begin running through the mesh, adding and possibly deleting
*  points as indicated by the irefin array.

*  The integer new is a count of the number of intervals in
*  the tentative mesh being generated by the refinement strategy.

      new = 1

*  The first interval is treated as a special case, since xx(1)
*  always remains in the mesh, and cannot be a fixed point.

      rlen = xxold(2) - xx(1)
      slen = rlen

      if (new + irefin(1)  .gt. nmax) goto 360

      if (irefin(1).gt.1) then
         dx = rlen/irefin(1)
         do 230 j = 2, irefin(1)
            new = new + 1
            xx(new) = xx(1) + (j-1)*dx
  230    continue
      endif

*  The fixed points specified by the fixpnt array cannot be
*  removed from the mesh.  The value fxnext indicates the 'next'
*  fixed point.  When no further fixed points remain to be processed
*  (or if nfxpnt = 0), fxnext is set to a value strictly larger than
*  the last mesh point, so that no mesh point can equal fxnext.
*  This way we need to compare only the values of xxold(i)
*  and fxnext.

      ifxcnt = 1
      if (nfxpnt .eq. 0) then
         fxnext = onep1*abs(xxold(nmsh))
      else
         fxnext = fixpnt(ifxcnt)
      endif

*  jtkout is a counter of the number of consecutive points that
*  have been removed from the mesh.

      jtkout = 0
      do 330 im = 2, ninter

         rlold = rlen
         rlen = xxold(im+1) - xxold(im)

*  If xxold(im) is the next fixed point, it cannot be removed
*  and so we don't test its error estimates.

         if(xxold(im) .eq. fxnext)  then

            ifxcnt = ifxcnt + 1
            if(ifxcnt .gt. nfxpnt) then
               fxnext = onep1*abs(xxold(nmsh))
            else
               fxnext = fixpnt(ifxcnt)
            endif

         elseif (irefin(im).eq.1) then

*  If xxold(im) is not a fixed point and irefin(im) = 1, we may wish
*  to remove point im from the mesh.

*  If we are considering removing points and jtkout = 0, this
*  is the first point in a possible consecutive set to be removed,
*  and we initialize phihat, which represents a maximum of
*  certain estimates.
*  If jtkout is not zero, previous points contiguous to this
*  point have been removed, and phihat retains its previous value.

            slen = slen + rlen

            if ( jtkout .lt. 2
     *             .and. r4(im) .le. 5e-1*fatt_r3) then

*  Increment the counter of removed points.
*  'Remove' the mesh point xxold(im) by not including it.

               jtkout = jtkout+1
               go to 330
            endif
*        end of logic for irefin(im) = 1
         endif
          if (new + irefin(im) .gt. nmax) then

*  If the new mesh contains too many points, branch out of the
*  loop to try alternative strategies.

         go to 360
         end if
         jtkout = 0
         new = new + 1
         xx(new) = xxold(im)
         if (irefin(im) .gt. 1) then
            dx = rlen/irefin(im)
            do 300 j = 2, irefin(im)
              new = new + 1
              xx(new) = xxold(im) + (j-1)*dx
  300       continue
         endif
         slen = rlen

         if (new .gt. nmax) then

*  If the new mesh contains too many points, branch out of the
*  loop to try alternative strategies.

            go to 360

         elseif (new .gt. 3*ninter)  then

*  Here, the new mesh does not exceed the specified maximum,
*  but has more than 3 times as many intervals as the old mesh.
*  Try doubling the mesh if possible.
            if (iprint .eq. 1) THEN
              write(msg,*) 'smpselcondmsh'
              CALL rprint(msg)
            ENDIF

            nmsh = nmold
            call dcopy(nmsh, xxold, 1, xx, 1)
            call dblmsh (nmsh, nmax, xx, nmold, xxold, maxmsh)
            ddouble = .true.
            return

         endif



  330 continue

*  To end up here, we have processed the entire interval,
*  and have neither exceeded the specified maximum nor
*  exceeded three times the number of intervals in the old
*  mesh.  The last mesh point remains unchanged.

      new = new + 1
      if (new   .gt. nmax) goto 360
      xx(new) = xxold(nmsh)
      nmsh = new
      maxmsh = .false.
      if (iprint .ge. 0) THEN
        write(msg,905) nmsh
        CALL rprint(msg)
      ENDIF

      return

  360 continue
*  To reach here, the number of mesh points created at some stage
*  of the refinement process was larger than the maximum permitted
*  value nmax.

      nmsh = nmold
      call dcopy(nmsh, xxold, 1, xx, 1)
      maxmsh = .true.


      return



  905 format(1h ,'smpselcondmsh.  new nmsh =',i8)
  910 format(1h ,'ihcomp',(10i5))
      end

c ===================================================================================

      subroutine moncondmsh(nmsh,xx,r1,r2,r3,fatt_r1r3,fatt_r3,
     *                        nptcond,r4,amg)

      implicit double precision (a-h,o-z)

      dimension  xx(*)
      dimension  amg(*), r4(*)

      intrinsic abs
      intrinsic max
      intrinsic int


*  blas: dcopy
*  double precision dlog

      logical use_c, comp_c
      common/algprs/ nminit, iprint,idum,  use_c, comp_c
      common/monpar/ sfatt_alpha, sfatt_r3, sfatt_r1r3

      parameter  ( zero = 0.0d+0, one = 1.0d+0 )
      character(len=150) msg

* the function moncond compute the monitor function based on the
* conditioning parameters and the factor used to perform the mesh selection

       do i=1,nmsh-1
          r4(i) = (xx(i+1)-xx(i))*dabs(amg(i+1)- amg(i))
       end do


       r2 = r4(1)
       do i=2,nmsh-1
          r2 = r2 + r4(i)
       end do
       do i=1,nmsh-1
          r4(i)=r4(i)+(xx(i+1)-xx(i))*(r2/(xx(nmsh)-xx(1)))*sfatt_alpha
       end do
c new matlab monitor function
c       do i=1,nmsh-1
c          r4(i) = r4(i)+(r2/(xx(nmsh)-xx(1)))*sfatt_alpha
c       end do


       r1 = r4(1)
       do i=2,nmsh-1
          r1 = max(r1,r4(i))
       end do
       do i=1,nmsh-1
         r4(i) = r4(i)/r1
       end do
       r1 = one

        r2 = r4(1)
        r1m = r4(1)
       do i=2,nmsh-1
          r1m = min(r1m,r4(i))
          r2 = r2 + r4(i)
       end do

       r3 = r2/(nmsh-1)

       fatt_r3  = r3*sfatt_r3

       fatt_r1r3= max(r3,r1*sfatt_r1r3)
c vecchio 0.5 nuovo 0.65
       nptm = 0
       nptr = 0
       do i=1,nmsh-1
           if (r4(i) .ge. fatt_r1r3)  nptm = nptm + 1
           if (r4(i) .le. fatt_r3)  nptr = nptr + 1
        enddo

        if (nptm .le. 1) then
           nptcond =  14
        elseif (nptm .le. 2) then
           nptcond =  10
        elseif (nptm .le. 4) then
           nptcond =  8
        elseif (nptm .le. 8) then
            nptcond = 6
        elseif (nptm .le. nmsh/20 ) then
            nptcond = 4
        else
            nptcond = 2
        endif


       if (iprint .eq. 1) THEN
         write(msg,901)r1,r3,fatt_r1r3,nptcond,nptm,nptr
         CALL rprint(msg)
       ENDIF



  901 format(1h ,'moncondmsh.', (1pe11.3), 2(1pe11.3), 3i10)

      end

c ===================================================================================
C
        DOUBLE PRECISION FUNCTION ABDNRM(NBLOKS,NTOP,NBOT,NOVRLP,
     *                            NRWBLK,NCLBLK,TOP,A,BOT)

C******************************************************************
C
C       ABDNRM IS USED IN CONJUNCTION WITH DONEST TO COMPUTE THE
C       CONDITION NUMBER OF AN ALMOST BLOCK DIAGONAL MATRIX LIKE
C       THE ONES HANDLED BY COLROW. [SEE COMMENTS IN COLROW, CRDCMP]

C               *****  AUXILIARY PROGRAMS  *****

C       THE BLAS ARE REQUIRED BY ABDNRM.
C       SPECIFICALLY, DASUM IS USED.
C

C******************************************************************

        INTEGER NBLOKS,NTOP,NBOT,NOVRLP,NRWBLK,NCLBLK
        DOUBLE PRECISION TOP(NTOP,*),A(NRWBLK,NCLBLK,*),BOT(NBOT,*)
        INTEGER J,K
        DOUBLE PRECISION MAX,DASUM,TEMP
      character(len=150) msg
        TEMP = 0.0D0
C
C       FIRST, GO OVER THE COLUMNS OF TOP AND THE FIRST BLOCK:
C
        DO 10 J = 1,NOVRLP
           TEMP = MAX(TEMP, DASUM(NTOP,TOP(1,J),1) +
     *                      DASUM(NRWBLK,A(1,J,1),1))
 10           CONTINUE
        DO 40 K = 1,NBLOKS-1
C
C          IN EACH BLOCK:
C
C          FIRST, THE COLUMNS FROM THE KTH BLOCK ALONE:
C
           DO 20 J = NOVRLP+1,NRWBLK
              TEMP = MAX(TEMP, DASUM(NRWBLK,A(1,J,K),1))
 20                 CONTINUE
C
C          NOW, TH COLUMNS WHICH INTERSECT BOTH THE KTH AND
C          (K+1)ST BLOCKS.
C
           DO 30 J = NRWBLK+1,NCLBLK
              TEMP = MAX(TEMP, DASUM(NRWBLK,A(1,J,K),1) +
     *                         DASUM(NRWBLK,A(1,J-NRWBLK,K+1),1))
 30                 CONTINUE
 40                    CONTINUE
C
C       FINALLY, THE COLUMNS OF THE LAST BLOCK WHICH DO NOT OVERLAP
C       WITH ANYTHING.
C
        DO 50 J = NOVRLP+1,NRWBLK
           TEMP = MAX(TEMP, DASUM(NRWBLK,A(1,J,NBLOKS),1))
 50           CONTINUE
C
C       AND THOSE COLUMNS OVERLAPPING WITH BOTH THE LAST BLOCK AND THE
C       BOTTOM BLOCK.
C
        DO 60 J = NRWBLK+1,NCLBLK
           TEMP = MAX(TEMP, DASUM(NRWBLK,A(1,J,NBLOKS),1) +
     *                      DASUM(NBOT,BOT(1,J-NRWBLK),1))
 60           CONTINUE
        ABDNRM = TEMP
        RETURN
        END

      SUBROUTINE DONEST (N, V, X, ISGN, EST, KASE)
      INTEGER N, ISGN(N), KASE
      DOUBLE PRECISION V(N), X(N), EST

C
C     DONEST ESTIMATES THE 1-NORM OF A SQUARE, DOUBLE PRECISION MATRIX  A.
C     REVERSE COMMUNICATION IS USED FOR EVALUATING
C     MATRIX-VECTOR PRODUCTS.
C
C     ON ENTRY
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX.  N .GE. 1.
C
C        ISGN    INTEGER(N)
C                USED AS WORKSPACE.
C
C        KASE    INTEGER
C                = 0.
C
C     ON INTERMEDIATE RETURNS
C
C        KASE    = 1 OR 2.
C
C        X       DOUBLE PRECISION(N)
C                MUST BE OVERWRITTEN BY
C
C                     A*X,             IF KASE=1,
C                     TRANSPOSE(A)*X,  IF KASE=2,
C
C                AND DONEST MUST BE RE-CALLED, WITH ALL THE OTHER
C                PARAMETERS UNCHANGED.
C
C     ON FINAL RETURN
C
C        KASE    = 0.
C
C        EST     DOUBLE PRECISION
C                CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).
C
C        V       DOUBLE PRECISION(N)
C                = A*W,   WHERE  EST = NORM(V)/NORM(W)
C                         (W  IS NOT RETURNED).
C
C     THIS VERSION DATED MARCH 16, 1988.
C     NICK HIGHAM, UNIVERSITY OF MANCHESTER.
C
C     MODIFIED FOR DOUBLE PRECISION ON JUNE 11, 1996.
C
C     REFERENCE
C     N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
C     THE ONE-NORM OF A REAL OR COMPLEX MATRIX, WITH APPLICATIONS
C     TO CONDITION  ESTIMATION, NUMERICAL ANALYSIS REPORT NO. 135,
C     UNIVERSITY OF MANCHESTER, MANCHESTER M13 9PL, ENGLAND.
C
C     SUBROUTINES AND FUNCTIONS
C     BLAS     IDAMAX, DASUM, DCOPY
C     GENERIC  ABS, NINT, FLOAT, SIGN
C
        INTRINSIC FLOAT
C        DOUBLE PRECISION FLOAT

        INTRINSIC ABS
        DOUBLE PRECISION ABS

        INTRINSIC SIGN
        DOUBLE PRECISION SIGN


      DOUBLE PRECISION DASUM
      INTEGER IDAMAX

      INTEGER ITMAX
      PARAMETER (ITMAX = 5)
      DOUBLE PRECISION ZERO,ONE,TWO
      PARAMETER (ZERO = 0.0D0)
      PARAMETER (ONE = 1.0D0)
      PARAMETER (TWO = 2.0D0)
C
C     INTERNAL VARIABLES
      INTEGER I, ITER, J, JLAST, JUMP
      DOUBLE PRECISION ALTSGN, ESTOLD, TEMP
C
      SAVE
C
      IF (KASE .EQ. 0) THEN
         DO 10,I = 1,N
            X(I) = ONE/FLOAT(N)
 10             CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      ENDIF
C
      GOTO (100, 200, 300, 400, 500) JUMP
C
C     ................ ENTRY   (JUMP = 1)
C     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
C
 100   CONTINUE
      IF (N .EQ. 1) THEN
         V(1) = X(1)
         EST = ABS(V(1))
C        ... QUIT
         GOTO 510
      ENDIF
      EST = DASUM(N,X,1)
C
      DO 110,I = 1,N
         X(I) = SIGN(ONE,X(I))
         ISGN(I) = NINT(X(I))
 110      CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
C
C     ................ ENTRY   (JUMP = 2)
C     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
C
 200   CONTINUE
      J = IDAMAX(N,X,1)
      ITER = 2
C
C     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
 220   CONTINUE
      DO 230,I = 1,N
         X(I) = ZERO
 230      CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      RETURN
C
C     ................ ENTRY   (JUMP = 3)
C     X HAS BEEN OVERWRITTEN BY A*X.
C
 300   CONTINUE
      CALL DCOPY(N,X,1,V,1)
      ESTOLD = EST
      EST = DASUM(N,V,1)
      DO 310,I = 1,N
         IF ( NINT( SIGN(ONE,X(I)) ) .NE. ISGN(I) ) GOTO 320
 310      CONTINUE
C     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GOTO 410
C
 320   CONTINUE
C     TEST FOR CYCLING.
      IF (EST .LE. ESTOLD) GOTO 410
C
      DO 330,I = 1,N
         X(I) = SIGN(ONE,X(I))
         ISGN(I) = NINT(X(I))
 330      CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
C
C     ................ ENTRY   (JUMP = 4)
C     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
C
 400   CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF (   (  X(JLAST) .NE. ABS(X(J))  ) .AND.
     +       (ITER .LT. ITMAX)   ) THEN
         ITER = ITER + 1
         GOTO 220
      ENDIF
C
C     ITERATION COMPLETE.  FINAL STAGE.
C
 410   CONTINUE
      ALTSGN = ONE
      DO 420,I = 1,N
         X(I) = ALTSGN * (ONE + FLOAT(I-1)/FLOAT(N-1))
         ALTSGN = -ALTSGN
 420      CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
C
C     ................ ENTRY   (JUMP = 5)
C     X HAS BEEN OVERWRITTEN BY A*X.
C
 500   CONTINUE
      TEMP = TWO*DASUM(N,X,1)/FLOAT(3*N)
      IF (TEMP. GT. EST) THEN
         CALL DCOPY(N,X,1,V,1)
         EST = TEMP
      ENDIF
C
 510   KASE = 0
      RETURN
C
      END


c ===================================================================================
C
C     THE AUGUST 27 1992 VERSION OF COLROW IN WHICH X IS NO LONGER
C     REQUIRED, WITH THE SOLUTION BEING RETURNED IN B, THE RIGHT
C     HAND SIDE.  IN ADDITION, ALL VARIABLES ARE EXPLICITLY DECLARED.
C     A PARAMETER "JOB" IS INCLUDED, TO SPECIFY WHICH OF A.X = B OR
C     TRANSPOSE(A).X = B IS TO BE SOLVED.

      SUBROUTINE COLROW(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B,IFLAG,JOB)
C
C***************************************************************
C
C  THIS PROGRAM SOLVES ONE OF THE LINEAR SYSTEMS  A*X = B OR
C  TRANSPOSE(A)*X = B, WHERE  A IS AN ALMOST BLOCK DIAGONAL
C  MATRIX OF THE FORM
C
C               TOPBLK
C               ARRAY(1)
C                     ARRAY(2)
C                          .
C                             .
C                                .
C                                   .
C                                    ARRAY(NBLOKS)
C                                           BOTBLK
C
C  WHERE
C           TOPBLK IS  NRWTOP  BY NOVRLP
C           ARRAY(K), K=1,NBLOKS, ARE NRWBLK BY NRWBLK+NOVRLP
C           BOTBLK IS NRWBOT BY NOVRLP,
C  AND
C           NOVRLP = NRWTOP + NRWBOT
C  WITH
C           NOVRLP.LE.NRWBLK .
C
C  THE LINEAR SYSTEM IS OF ORDER  N = NBLOKS*NRWBLK + NOVRLP.
C
C  THE METHOD IMPLEMENTED IS BASED ON GAUSS ELIMINATION WITH
C  ALTERNATE ROW AND COLUMN ELIMINATION WITH PARTIAL PIVOTING,
C  WHICH PRODUCES A STABLE DECOMPOSITION OF THE MATRIX  A
C  WITHOUT INTRODUCING FILL-IN.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  TO OBTAIN A SINGLE PRECISION VERSION OF THIS PACKAGE, REMOVE
C  ALL DOUBLE PRECISION STATEMENTS.  THERE IS ONE SUCH STATEMENT
C  IN C O L R O W, THREE IN C R D C M P, AND TWO IN C R S O L V.
C  IN ADDITION, REFERENCES TO BUILT-IN FUNCTIONS DABS AND DMAX1
C  MUST BE REPLACED BY ABS AND AMAX1, RESPECTIVELY.  DABS OCCURS
C  NINE TIMES, IN C R D C M P.  DMAX1 OCCURS FOUR TIMES, IN
C  C R D C M P.  FINALLY, ZERO IS INITIALISED TO 0.D0 IN A
C  DATA STATEMENT IN C R D C M P.  THIS MUST BE REPLACED BY:
C               DATA ZERO/0.0/
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY ...
C
C               N      - INTEGER
C                         THE ORDER OF THE LINEAR SYSTEM,
C                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         THE FIRST BLOCK OF THE ALMOST BLOCK
C                         DIAGONAL MATRIX A
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
C                         BY NCLBLK BLOCK OF THE MATRIX A
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         THE LAST BLOCK OF THE MATRIX A
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         WORK SPACE
C
C                    B - DOUBLE PRECISION(N)
C                         THE RIGHT HAND SIDE VECTOR
C
C               JOB    - INTEGER, INDICATING:
C                      = 0: SOLVE A*X = B;
C                      NON-ZERO: SOLVE TRANSPOSE(A)*X = B.
C
C       *** ON RETURN  ...
C
C               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE
C                        DESIRED DECOMPOSITION OF THE MATRIX A
C                        (IF IFLAG = 0)
C
C                PIVOT - INTEGER(N)
C                         RECORDS THE PIVOTING INDICES DETER-
C                         MINED IN THE DECOMPOSITION
C
C                    B - DOUBLE PRECISION(N)
C                         THE SOLUTION VECTOR (IF IFLAG = 0)
C
C               IFLAG  - INTEGER
C                         =  1, IF INPUT PARAMETERS ARE INVALID
C                         = -1, IF MATRIX IS SINGULAR
C                         =  0, OTHERWISE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  AUXILIARY PROGRAMS  *****
C
C       CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
C    *     BOTBLK,NRWBOT,PIVOT,IFLAG)
C            - DECOMPOSES THE MATRIX  A  USING MODIFIED
C              ALTERNATE ROW AND COLUMN ELIMINATON WITH
C              PARTIAL PIVOTING, AND IS USED FOR THIS
C              PURPOSE IN C O L R O W.
C              THE ARGUMENTS ARE AS IN C O L R O W.
C
C       CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
C    *     BOTBLK,NRWBOT,PIVOT,B,JOB)
C            - SOLVES THE SYSTEM A*X = B ONCE A IS DECOMPOSED.
C              THE ARGUMENTS ARE ALL AS IN C O L R O W.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       THE SUBROUTINE  C O L R O W  AUTOMATICALLY SOLVES THE
C  INPUT SYSTEM WHEN IFLAG=0.  C O L R O W  IS CALLED ONLY ONCE
C  FOR A GIVEN SYSTEM. THE SOLUTION FOR A SEQUENCE OF P RIGHT
C  HAND SIDES CAN BE OBTAINED BY ONE CALL TO  C O L R O W  AND
C  P-1 CALLS TO CRSLVE ONLY. SINCE THE ARRAYS TOPBLK,ARRAY,
C  BOTBLK AND PIVOT CONTAIN THE DECOMPOSITION OF THE GIVEN
C  COEFFICIENT MATRIX AND PIVOTING INFORMATION ON RETURN FROM
C  C O L R O W , THEY MUST NOT BE ALTERED BETWEEN SUCCESSIVE
C  CALLS TO CRSLVE WITH THE SAME LEFT HAND SIDES. FOR THE
C  SAME REASON, IF THE USER WISHES TO SAVE THE COEFFICIENT
C  MATRIX, THE ARRAYS TOPBLK,ARRAY,BOTBLK MUST BE COPIED
C  BEFORE A CALL TO  C O L R O W .
C
C*************************************************************************
        IMPLICIT NONE
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,B
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,PIVOT(*),
     *          IFLAG,JOB, IDAMAX,i,IFAIL
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),B(*)

c        NOUT = 6

C       DO THE FACTORIZATION USING CRDCMP:
C
        CALL CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
     *          BOTBLK,NRWBOT,PIVOT,IFLAG)
        IF(IFLAG.NE.0) RETURN
C
c     *****************solving the linear system********************
        job=0

        CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,NBLOKS,
     *          BOTBLK,NRWBOT,PIVOT,B,JOB)
        RETURN
        END


      SUBROUTINE INVERSE(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,
     *          NBLOKS,BOTBLK,NRWBOT,PIVOT,INMAT)
C******************************************************************
C     INVERSE COMPUTES THE INVERSE OF A MATRIX
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        IMPLICIT NONE
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,WORK,
     *          INMAT
        INTEGER N,NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,PIVOT(*)
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),WORK(N),
     *          INMAT(N,N)
        INTEGER K,L,J

        do 300 k=1,N
           do 330 l=1,N
              WORK(l)=0.0d0
              if (l.eq.k) WORK(l)=1.0d0
 330       continue

          CALL CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *                 NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,
     *                 WORK,0)

          do 440 l=1,N
             INMAT(l,k)=WORK(l)
 440         continue
 300      continue

          RETURN
          END

c ===================================================================================

      SUBROUTINE CRDCMP(N,TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,NCLBLK,
     *   NBLOKS,BOTBLK,NRWBOT,PIVOT,IFLAG)
C
C***************************************************************
C
C  C R D C M P DECOMPOSES THE ALMOST BLOCK DIAGONAL MATRIX A
C  USING MODIFIED ALTERNATE ROW AND COLUMN ELIMINATION WITH
C  PARTIAL PIVOTING.  THE MATRIX  A  IS STORED IN THE ARRAYS
C  TOPBLK, ARRAY, AND BOTBLK.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY ...
C
C               N      - INTEGER
C                         THE ORDER OF THE LINEAR SYSTEM,
C                         GIVEN BY NBLOKS*NRWBLK + NOVRLP
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         THE FIRST BLOCK OF THE ALMOST BLOCK
C                         DIAGONAL MATRIX A TO BE DECOMPOSED
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         ARRAY(,,K) CONTAINS THE K-TH NRWBLK
C                         BY NCLBLK BLOCK OF THE MATRIX A
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         THE LAST BLOCK OF THE MATRIX A
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         WORK SPACE
C
C       *** ON RETURN  ...
C
C               TOPBLK,ARRAY,BOTBLK - ARRAYS CONTAINING THE
C                        DESIRED DECOMPOSITION OF THE MATRIX A
C                        (IF IFLAG = 0)
C
C                PIVOT - INTEGER(N)
C                         RECORDS THE PIVOTING INDICES DETER-
C                         MINED IN THE DECOMPOSITION
C
C               IFLAG  - INTEGER
C                         =  1, IF INPUT PARAMETERS ARE INVALID
C                         = -1, IF MATRIX IS SINGULAR
C                         =  0, OTHERWISE
C
C***************************************************************
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER PIVOT(*)
      DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),BOTBLK(NRWBOT,*)
      DATA ZERO / 0.0D+0 /
      character(len=150) msg
C
C***************************************************************
C
C          ****  DEFINE THE CONSDTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
      IFLAG = 0
      PIVMAX = ZERO
      NRWTP1 = NRWTOP+1
      NROWEL = NRWBLK-NRWTOP
      NRWEL1 = NROWEL+1
      NVRLP0 = NOVRLP-1
C
C***************************************************************
C
C          ****  CHECK VALIDITY OF THE INPUT PARAMETERS....
C
C               IF PARAMETERS ARE INVALID THEN TERMINATE AT 10;
C                                         ELSE CONTINUE AT 100.
C
C***************************************************************
C
      IF (N.NE.NBLOKS*NRWBLK+NOVRLP) GO TO 10
      IF (NOVRLP.NE.NRWTOP+NRWBOT) GO TO 10
      IF (NCLBLK.NE.NOVRLP+NRWBLK) GO TO 10
      IF (NOVRLP.GT.NRWBLK) GO TO 10
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE ACCEPTABLE - CONTINUE AT 100.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      GO TO 20
   10 CONTINUE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          PARAMETERS ARE INVALID.  SET IFLAG = 1, AND TERMINATE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IFLAG = 1
      RETURN
   20 CONTINUE
C
C***************************************************************
C
C               ****  FIRST, IN TOPBLK....
C
C***************************************************************
C
C          ***  APPLY NRWTOP COLUMN ELIMINATIONS WITH COLUMN
C                 PIVOTING ....
C
C***************************************************************
C
      DO 110 I = 1, NRWTOP
         IPLUS1 = I+1
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         IPVT = I
         COLMAX = DABS(TOPBLK(I,I))
         DO 30 J = IPLUS1, NOVRLP
            TEMPIV = DABS(TOPBLK(I,J))
            IF (TEMPIV.LE.COLMAX) GO TO 30
            IPVT = J
            COLMAX = TEMPIV
   30    CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         IF (PIVMAX+COLMAX.EQ.PIVMAX) GO TO 410
         PIVMAX = DMAX1(COLMAX,PIVMAX)
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         PIVOT(I) = IPVT
         IF (IPVT.EQ.I) GO TO 60
         DO 40 L = I, NRWTOP
            SWAP = TOPBLK(L,IPVT)
            TOPBLK(L,IPVT) = TOPBLK(L,I)
            TOPBLK(L,I) = SWAP
   40    CONTINUE
         DO 50 L = 1, NRWBLK
            SWAP = ARRAY(L,IPVT,1)
            ARRAY(L,IPVT,1) = ARRAY(L,I,1)
            ARRAY(L,I,1) = SWAP
   50    CONTINUE
   60    CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM COLUMN
C                       ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         COLPIV = TOPBLK(I,I)
         DO 100 J = IPLUS1, NOVRLP
            COLMLT = TOPBLK(I,J)/COLPIV
            TOPBLK(I,J) = COLMLT
            IF (IPLUS1.GT.NRWTOP) GO TO 80
            DO 70 L = IPLUS1, NRWTOP
               TOPBLK(L,J) = TOPBLK(L,J)-COLMLT*TOPBLK(L,I)
   70       CONTINUE
   80       CONTINUE
            DO 90 L = 1, NRWBLK
               ARRAY(L,J,1) = ARRAY(L,J,1)-COLMLT*ARRAY(L,I,1)
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
C
C***************************************************************
C
C          ****  IN EACH BLOCK ARRAY(,,K)....
C
C***************************************************************
C
      INCR = 0
      DO 320 K = 1, NBLOKS
         KPLUS1 = K+1
C
C          *****************************************************
C
C          ***  FIRST APPLY NRWBLK-NRWTOP ROW ELIMINATIONS WITH
C                       ROW PIVOTING....
C
C          *****************************************************
C
         DO 180 J = NRWTP1, NRWBLK
            JPLUS1 = J+1
            JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            IPVT = JMINN
            ROWMAX = DABS(ARRAY(JMINN,J,K))
            LOOP = JMINN+1
            DO 120 I = LOOP, NRWBLK
               TEMPIV = DABS(ARRAY(I,J,K))
               IF (TEMPIV.LE.ROWMAX) GO TO 120
               IPVT = I
               ROWMAX = TEMPIV
  120       CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            IF (PIVMAX+ROWMAX.EQ.PIVMAX) GO TO 410
            PIVMAX = DMAX1(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            INCRJ = INCR+J
            PIVOT(INCRJ) = INCR+IPVT+NRWTOP
            IF (IPVT.EQ.JMINN) GO TO 140
            DO 130 L = J, NCLBLK
               SWAP = ARRAY(IPVT,L,K)
               ARRAY(IPVT,L,K) = ARRAY(JMINN,L,K)
               ARRAY(JMINN,L,K) = SWAP
  130       CONTINUE
  140       CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            ROWPIV = ARRAY(JMINN,J,K)
            DO 150 I = LOOP, NRWBLK
               ARRAY(I,J,K) = ARRAY(I,J,K)/ROWPIV
  150       CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            DO 170 L = JPLUS1, NCLBLK
               ROWMLT = ARRAY(JMINN,L,K)
               DO 160 I = LOOP, NRWBLK
                  ARRAY(I,L,K) = ARRAY(I,L,K)-ROWMLT*ARRAY(I,J,K)
  160          CONTINUE
  170       CONTINUE
  180    CONTINUE
C
C          *****************************************************
C
C          ***  NOW APPLY NRWTOP COLUMN ELIMINATIONS WITH
C                      COLUMN PIVOTING....
C
C          *****************************************************
C
         DO 310 I = NRWEL1, NRWBLK
            IPLUSN = I+NRWTOP
            IPLUS1 = I+1
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE COLUMN PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            IPVT = IPLUSN
            COLMAX = DABS(ARRAY(I,IPVT,K))
            LOOP = IPLUSN+1
            DO 190 J = LOOP, NCLBLK
               TEMPIV = DABS(ARRAY(I,J,K))
               IF (TEMPIV.LE.COLMAX) GO TO 190
               IPVT = J
               COLMAX = TEMPIV
  190       CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            IF (PIVMAX+COLMAX.EQ.PIVMAX) GO TO 410
            PIVMAX = DMAX1(COLMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE COLUMNS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            INCRN = INCR+IPLUSN
            PIVOT(INCRN) = INCR+IPVT
            IRWBLK = IPLUSN-NRWBLK
            IF (IPVT.EQ.IPLUSN) GO TO 240
            DO 200 L = I, NRWBLK
               SWAP = ARRAY(L,IPVT,K)
               ARRAY(L,IPVT,K) = ARRAY(L,IPLUSN,K)
               ARRAY(L,IPLUSN,K) = SWAP
  200       CONTINUE
            IPVBLK = IPVT-NRWBLK
            IF (K.EQ.NBLOKS) GO TO 220
            DO 210 L = 1, NRWBLK
               SWAP = ARRAY(L,IPVBLK,KPLUS1)
               ARRAY(L,IPVBLK,KPLUS1) = ARRAY(L,IRWBLK,KPLUS1)
               ARRAY(L,IRWBLK,KPLUS1) = SWAP
  210       CONTINUE
            GO TO 240
  220       CONTINUE
            DO 230 L = 1, NRWBOT
               SWAP = BOTBLK(L,IPVBLK)
               BOTBLK(L,IPVBLK) = BOTBLK(L,IRWBLK)
               BOTBLK(L,IRWBLK) = SWAP
  230       CONTINUE
  240       CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS AND PERFORM COLUMN
C                       ELIMINATION
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            COLPIV = ARRAY(I,IPLUSN,K)
            DO 300 J = LOOP, NCLBLK
               COLMLT = ARRAY(I,J,K)/COLPIV
               ARRAY(I,J,K) = COLMLT
               IF (I.EQ.NRWBLK) GO TO 260
               DO 250 L = IPLUS1, NRWBLK
                  ARRAY(L,J,K) = ARRAY(L,J,K)-COLMLT*ARRAY(L,IPLUSN,K)
  250          CONTINUE
  260          CONTINUE
               JRWBLK = J-NRWBLK
               IF (K.EQ.NBLOKS) GO TO 280
               DO 270 L = 1, NRWBLK
                  ARRAY(L,JRWBLK,KPLUS1) = ARRAY(L,JRWBLK,KPLUS1)-COLMLT
     *               *ARRAY(L,IRWBLK,KPLUS1)
  270          CONTINUE
               GO TO 300
  280          CONTINUE
               DO 290 L = 1, NRWBOT
                  BOTBLK(L,JRWBLK) = BOTBLK(L,JRWBLK)-COLMLT*BOTBLK(L,
     *               IRWBLK)
  290          CONTINUE
  300       CONTINUE
  310    CONTINUE
         INCR = INCR+NRWBLK
  320 CONTINUE
C
C***************************************************************
C
C          ****  FINALLY, IN BOTBLK....
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C          ***  APPLY NRWBOT ROW ELIMINATIONS WITH ROW
C                  PIVOTING....
C
C               IF BOT HAS JUST ONE ROW GO TO 500
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IF (NRWBOT.EQ.1) GO TO 400
      DO 390 J = NRWTP1, NVRLP0
         JPLUS1 = J+1
         JMINN = J-NRWTOP
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               DETERMINE ROW PIVOT AND PIVOT INDEX
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         IPVT = JMINN
         ROWMAX = DABS(BOTBLK(JMINN,J))
         LOOP = JMINN+1
         DO 330 I = LOOP, NRWBOT
            TEMPIV = DABS(BOTBLK(I,J))
            IF (TEMPIV.LE.ROWMAX) GO TO 330
            IPVT = I
            ROWMAX = TEMPIV
  330    CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               TEST FOR SINGULARITY:
C
C                       IF SINGULAR THEN TERMINATE AT 1000;
C                                   ELSE CONTINUE.
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         IF (PIVMAX+ROWMAX.EQ.PIVMAX) GO TO 410
         PIVMAX = DMAX1(ROWMAX,PIVMAX)
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               IF NECESSARY INTERCHANGE ROWS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         INCRJ = INCR+J
         PIVOT(INCRJ) = INCR+IPVT+NRWTOP
         IF (IPVT.EQ.JMINN) GO TO 350
         DO 340 L = J, NOVRLP
            SWAP = BOTBLK(IPVT,L)
            BOTBLK(IPVT,L) = BOTBLK(JMINN,L)
            BOTBLK(JMINN,L) = SWAP
  340    CONTINUE
  350    CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               COMPUTE MULTIPLIERS
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         ROWPIV = BOTBLK(JMINN,J)
         DO 360 I = LOOP, NRWBOT
            BOTBLK(I,J) = BOTBLK(I,J)/ROWPIV
  360    CONTINUE
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               PERFORM ROW ELIMINATION WITH COLUMN INDEXING
C
C             CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
         DO 380 L = JPLUS1, NOVRLP
            ROWMLT = BOTBLK(JMINN,L)
            DO 370 I = LOOP, NRWBOT
               BOTBLK(I,L) = BOTBLK(I,L)-ROWMLT*BOTBLK(I,J)
  370       CONTINUE
  380    CONTINUE
  390 CONTINUE
  400 CONTINUE
C
C***************************************************************
C
C          DONE PROVIDED THE LAST ELEMENT IS NOT ZERO
C
C***************************************************************
C
      IF (PIVMAX+DABS(BOTBLK(NRWBOT,NOVRLP)).NE.PIVMAX) RETURN
C
C***************************************************************
C
C       ****  MATRIX IS SINGULAR - SET IFLAG = - 1.
C                                  TERMINATE AT 1000.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
  410 CONTINUE
      IFLAG = -1
      RETURN
      END

c ===================================================================================

        SUBROUTINE CRSLVE(TOPBLK,NRWTOP,NOVRLP,ARRAY,NRWBLK,
     *             NCLBLK,NBLOKS,BOTBLK,NRWBOT,PIVOT,B,JOB)
C
C***************************************************************
C
C  C R S L V E  SOLVES THE LINEAR SYSTEM
C                       A*X = B
C  USING THE DECOMPOSITION ALREADY GENERATED IN  C R D C M P.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               *****  PARAMETERS  *****
C
C       *** ON ENTRY  ...
C
C               TOPBLK - DOUBLE PRECISION(NRWTOP,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NOVRLP - INTEGER
C                         THE NUMBER OF COLUMNS IN WHICH SUCC-
C                         ESSIVE BLOCKS OVERLAP, WHERE
C                                NOVRLP = NRWTOP + NRWBOT
C
C               NRWTOP - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK TOPBLK
C
C               ARRAY  - DOUBLE PRECISION(NRWBLK,NCLBLK,NBLOKS)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBLK - INTEGER
C                         NUMBER OF ROWS IN K-TH BLOCK
C
C               NCLBLK - INTEGER
C                         NUMBER OF COLUMNS IN K-TH BLOCK
C
C               NBLOKS - INTEGER
C                         NUMBER OF NRWBLK BY NCLBLK BLOCKS IN
C                         THE MATRIX A
C
C               BOTBLK - DOUBLE PRECISION(NRWBOT,NOVRLP)
C                         OUTPUT FROM  C R D C M P
C
C               NRWBOT - INTEGER
C                         NUMBER OF ROWS IN THE BLOCK BOTBLK
C
C                PIVOT - INTEGER(N)
C                         THE PIVOT VECTOR FROM  C R D C M P
C
C                    B - DOUBLE PRECISION(N)
C                         THE RIGHT HAND SIDE VECTOR
C
C               JOB    - INTEGER, INDICATING:
C                      = 0: SOLVE A*X = B;
C                      NON-ZERO: SOLVE TRANSPOSE(A)*X = B.
C
C       *** ON RETURN  ...
C
C                    B - DOUBLE PRECISION(N)
C                         THE SOLUTION VECTOR
C
C***************************************************************
C
        IMPLICIT NONE
        DOUBLE PRECISION TOPBLK,ARRAY,BOTBLK,B
        DOUBLE PRECISION DOTPRD,BJ,XINCRJ,BINCRJ,SWAP,BI
        INTEGER NRWTOP,NOVRLP,NRWBLK,NCLBLK,NBLOKS,NRWBOT,PIVOT(*),
     *          JOB
        DIMENSION TOPBLK(NRWTOP,*),ARRAY(NRWBLK,NCLBLK,*),
     *          BOTBLK(NRWBOT,*),B(*)
        INTEGER NRWTP1,NRWBK1,NVRLP1,NRWBT1,NROWEL,NVRLP0,NBLKS1,
     *          NBKTOP,J,I,LOOP,INCR,INCRJ,INCRI,JPIVOT,JRWTOP,
     *          LL,L1,IPLUSN,INCRN,NRWTP0,NRWEL1,K,INCRTP,NRWBTL,
     *          IPVTN,NRWELL,IPVTI,L
C
C***************************************************************
C
C          ****  DEFINE THE CONSTANTS USED THROUGHOUT  ****
C
C***************************************************************
C
        NRWTP1 = NRWTOP+1
        NRWBK1 = NRWBLK+1
        NVRLP1 = NOVRLP+1
        NRWTP0 = NRWTOP-1
        NRWBT1 = NRWBOT+1
        NROWEL = NRWBLK-NRWTOP
        NRWEL1 = NROWEL+1
        NVRLP0 = NOVRLP-1
        NBLKS1 = NBLOKS+1
        NBKTOP = NRWBLK+NRWTOP
C
C       IF JOB IS NON-ZERO, TRANSFER TO THE SECTION DEALING WITH
C       TRANSPOSE(A)*X = B.
C
        IF ( JOB .NE. 0 ) GO TO 530
C
C***************************************************************
C
C               ****  FORWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST, IN TOPBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 130 J = 1,NRWTOP
           B(J) = B(J)/TOPBLK(J,J)
           IF(J.EQ.NRWTOP)GO TO 120
              BJ = -B(J)
              LOOP = J+1
              DO 110 I = LOOP,NRWTOP
                 B(I) = B(I)+TOPBLK(I,J)*BJ
 110                        CONTINUE
 120                                CONTINUE
 130                                     CONTINUE
C
C       ********************************************************
C
C          ***  IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCR = 0
        DO 280 K = 1,NBLOKS
           INCRTP = INCR+NRWTOP
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 220 J = 1,NRWTOP
              INCRJ = INCR+J
              XINCRJ = -B(INCRJ)
              DO 210 I = 1,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
 210                        CONTINUE
 220                                CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 240 J = NRWTP1,NRWBLK
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 225
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
 225                        CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 230 I = LOOP,NRWBLK
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*BINCRJ
 230                        CONTINUE
 240                                CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 270 J = NRWBK1,NBKTOP
              INCRJ = INCR+J
              JRWTOP = J -NRWTOP
              B(INCRJ) = B(INCRJ)/ARRAY(JRWTOP,J,K)
              IF(J.EQ.NBKTOP)GO TO 260
                 XINCRJ = -B(INCRJ)
                 LOOP = J-NRWTP0
                 DO 250 I = LOOP,NRWBLK
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
 250                              CONTINUE
 260                                         CONTINUE
 270                                                 CONTINUE
           INCR = INCR+NRWBLK
 280            CONTINUE
C
C       ********************************************************
C
C          ***  FINALLY, IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD MODIFICATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        INCRTP = INCR+NRWTOP
        DO 320 J = 1,NRWTOP
           INCRJ = INCR+J
           XINCRJ = -B(INCRJ)
           DO 310 I = 1,NRWBOT
              INCRI = INCRTP+I
              B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
 310                  CONTINUE
 320                       CONTINUE
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               FORWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        IF(NRWBOT.EQ.1)GO TO 350
           DO 340 J = NRWTP1,NVRLP0
              INCRJ = INCR+J
              JPIVOT = PIVOT(INCRJ)
              IF(JPIVOT.EQ.INCRJ)GO TO 325
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(JPIVOT)
                 B(JPIVOT) = SWAP
 325                        CONTINUE
              BINCRJ = -B(INCRJ)
              LOOP = J-NRWTP0
              DO 330 I = LOOP,NRWBOT
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*BINCRJ
 330                        CONTINUE
 340                                CONTINUE
 350                                     CONTINUE
C
C***************************************************************
C
C               ****  BACKWARD RECURSION  ****
C
C***************************************************************
C
C          ***  FIRST IN BOTBLK....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 430 LL = 1,NRWBOT
           J = NVRLP1-LL
           INCRJ = INCR+J
           NRWBTL = NRWBT1-LL
           B(INCRJ) = B(INCRJ)/BOTBLK(NRWBTL,J)
           IF(LL.EQ.NRWBOT)GO TO 420
              XINCRJ = -B(INCRJ)
              LOOP = NRWBOT-LL
              DO 410 I = 1,LOOP
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+BOTBLK(I,J)*XINCRJ
 410                        CONTINUE
 420                                CONTINUE
 430                                     CONTINUE
C
C       ********************************************************
C
C          ***  THEN IN EACH BLOCK ARRAY(,,K)....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 490 L = 1,NBLOKS
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           K = NBLKS1-L
           INCR = INCR-NRWBLK
           DO 450 L1 = NRWEL1,NRWBLK
              I = NRWBLK+NRWEL1-L1
              IPLUSN = I+NRWTOP
              LOOP = IPLUSN+1
              INCRN = INCR+IPLUSN
              DOTPRD = B(INCRN)
              DO 440 J = LOOP,NCLBLK
                 INCRJ = INCR+J
                 DOTPRD = DOTPRD-ARRAY(I,J,K)*B(INCRJ)
 440                        CONTINUE
              B(INCRN) = DOTPRD
              IPVTN = PIVOT(INCRN)
              IF(INCRN.EQ.IPVTN)GO TO 445
                 SWAP = B(INCRN)
                 B(INCRN) = B(IPVTN)
                 B(IPVTN) = SWAP
 445                        CONTINUE
 450                                CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD MODIFICATION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           INCRTP = INCR+NRWTOP
           DO 460 J = NRWBK1,NCLBLK
              INCRJ = INCR+J
              XINCRJ = -B(INCRJ)
              DO 455 I = 1,NROWEL
                 INCRI = INCRTP+I
                 B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
 455                        CONTINUE
 460                                CONTINUE
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD SOLUTION
C
C          CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
           DO 480 LL = 1,NROWEL
              J = NRWBK1-LL
              INCRJ = INCR+J
              NRWELL = NRWEL1-LL
              B(INCRJ) = B(INCRJ)/ARRAY(NRWELL,J,K)
              IF(LL.EQ.NROWEL)GO TO 470
                 XINCRJ = -B(INCRJ)
                 LOOP = NROWEL-LL
                 DO 465 I = 1,LOOP
                    INCRI = INCRTP+I
                    B(INCRI) = B(INCRI)+ARRAY(I,J,K)*XINCRJ
 465                              CONTINUE
 470                                         CONTINUE
 480                                                 CONTINUE
 490                                                      CONTINUE
C
C       ********************************************************
C
C          ***  IN TOPBLK FINISH WITH....
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C               BACKWARD ELIMINATION
C
C       CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
        DO 520 L = 1,NRWTOP
           I = NRWTP1-L
           LOOP = I+1
           DOTPRD = B(I)
           DO 510 J = LOOP,NOVRLP
              DOTPRD = DOTPRD-TOPBLK(I,J)*B(J)
 510                  CONTINUE
           B(I) = DOTPRD
           IPVTI = PIVOT(I)
           IF(I.EQ.IPVTI)GO TO 515
                 SWAP = B(I)
                 B(I) = B(IPVTI)
                 B(IPVTI) = SWAP
 515                     CONTINUE
 520                          CONTINUE
C
C       RETURN FROM THE SOLUTION OF A.X = B.
        RETURN
C
C       IF JOB IS NON-ZERO, SOLVE TRANSPOSE(A)*X = B:
C
 530       CONTINUE

C       FIRST, FORWARD ELIMINATION OF RHS USING TRANSPOSE(U).

        DO 540 I = 1,NRWTOP
           IPVTI = PIVOT(I)
           IF ( I .NE. IPVTI ) THEN
              SWAP = B(I)
              B(I) = B(IPVTI)
              B(IPVTI) = SWAP
           ENDIF
           BI = -B(I)
           LOOP = I+1
           DO 535 J = LOOP,NOVRLP
              B(J) = B(J) + BI*TOPBLK(I,J)
 535                CONTINUE
 540                   CONTINUE

C       IN EACH BLOCK, K = 1,..,NBLOKS:

        INCR = NRWTOP
        DO 590 K = 1,NBLOKS

C          FIRST, THE FORWARD SOLUTION.

           DO 550 J = 1,NROWEL
              INCRJ = INCR + J
              DO 545 I = 1,J-1
                 B(INCRJ) = B(INCRJ) - ARRAY(I,NRWTOP+J,K)*B(INCR+I)
 545                      CONTINUE
              B(INCRJ) = B(INCRJ)/ARRAY(J,NRWTOP+J,K)
 550                 CONTINUE

C           FORWARD MODIFICATION.

            DO 570 I = 1,NOVRLP
               INCRI = INCR + NROWEL + I
               LOOP = NRWBLK + I
               DO 560 J = 1,NROWEL
                  INCRJ = INCR + J
                  B(INCRI) = B(INCRI) - ARRAY(J,LOOP,K)*B(INCRJ)
 560                        CONTINUE
 570                               CONTINUE

C           NOW, FORWARD ELIMINATION OF RHS USING TRANSPOSE(U). THIS
C           CORRESPONDS TO THE LOOP 540 ABOVE.

            INCR = INCR + NROWEL
            DO 580 I = 1,NRWTOP
               INCRI = INCR + I
               IPVTI = PIVOT(INCRI)
               IF ( INCRI .NE. IPVTI ) THEN
                  SWAP = B(INCRI)
                  B(INCRI) = B(IPVTI)
                  B(IPVTI) = SWAP
               ENDIF
               LOOP = NROWEL + I
               BI = -B(INCRI)
               DO 575 J = I+1,NOVRLP
                  INCRJ = INCR+J
                  L = NRWBLK + J
                  B(INCRJ) = B(INCRJ) + BI*ARRAY(LOOP,L,K)
 575                        CONTINUE
 580                               CONTINUE
            INCR = INCR + NRWTOP
 590           CONTINUE

C       FINALLY, FINISH WITH NRWBOT SOLUTIONS:

        DO 600 J = 1,NRWBOT
           INCRJ = INCR + J
           DO 595 I = 1,J-1
              B(INCRJ) = B(INCRJ) - BOTBLK(I,J+NRWTOP)*B(INCR+I)
 595                CONTINUE
           B(INCRJ) = B(INCRJ)/BOTBLK(J,J+NRWTOP)
 600          CONTINUE


C       NOW, THE BACKWARD PASS:


C       FIRST, BACKWARD SOLUTION IN BOTBLK:

        INCRJ = INCR + NRWBOT
        DO 610 J = 1,NRWBOT-1
           INCRJ = INCRJ - 1
           DO 605 I = NRWBOT-J+1,NRWBOT
              INCRI = INCR + I
              B(INCRJ) = B(INCRJ) - BOTBLK(I,NOVRLP-J)*B(INCRI)
 605                CONTINUE

           IF ( INCRJ .NE. PIVOT(INCRJ) ) THEN
              SWAP = B(INCRJ)
              B(INCRJ) = B(PIVOT(INCRJ))
              B(PIVOT(INCRJ)) = SWAP
           ENDIF
 610          CONTINUE

C       NOW DO THE DEFERRED OPERATIONS IN BOTBLOK:

        DO 620 J = 1,NRWTOP
           INCRJ = INCR - J + 1
           DO 615 I = 1,NRWBOT
              INCRI = INCR + I
              B(INCRJ) = B(INCRJ) - BOTBLK(I,NRWTP1-J)*B(INCRI)
 615                CONTINUE
 620                   CONTINUE


C       NOW, IN EACH BLOCK, K = NBLOKS,..,1:
        DO 800 K = NBLOKS,1,-1

C          FIRST, THE BACKSUBSTITUIONS:

           DO 630 J = 1,NRWTOP
              INCRJ = INCR - J + 1
              LOOP = NBKTOP - J + 1
              DO 625 I = 1,J-1
                 INCRI = INCR - I + 1
                 B(INCRJ) = B(INCRJ) - ARRAY(NRWBLK-I+1,LOOP,K)*B(INCRI)
 625                      CONTINUE
              B(INCRJ) = B(INCRJ)/ARRAY(NRWBLK-J+1,LOOP,K)
 630                CONTINUE

C          THEN THE BACKWARD SOLUTION IN THE KTH BLOCK:

           DO 650 J = 1,NROWEL
              INCRJ = INCR - NRWTOP -J + 1
              DO 645 I = 1,J+NRWTOP-1
                 INCRI = INCRJ + I
                 B(INCRJ) = B(INCRJ) -
     *           ARRAY(NRWBLK-NRWTOP-J+1+I,NRWBLK-J+1,K)*B(INCRI)
 645                      CONTINUE
              IF ( INCRJ .NE. PIVOT(INCRJ) ) THEN
                 SWAP = B(INCRJ)
                 B(INCRJ) = B(PIVOT(INCRJ))
                 B(PIVOT(INCRJ)) = SWAP
              ENDIF
 650                CONTINUE

C          NOW, THE DEFERRED OPERATIONS ON B:

           INCR = INCR - NRWBLK
           DO 660 J = 1,NRWTOP
              INCRJ = INCR + J - NRWTOP
              DO 655 I = 1,NRWBLK
                 INCRI = INCR + I
                 B(INCRJ) = B(INCRJ) - ARRAY(I,J,K)*B(INCRI)
 655                     CONTINUE
 660                           CONTINUE
 800                              CONTINUE

C       FINALLY, THE LAST SET OF BACK-SUBSTITUTIONS IN TOPBLK:

        DO 900 J = 1,NRWTOP
           INCRJ = NRWTOP -J + 1
           DO 850 I = INCRJ+1,NRWTOP
              B(INCRJ) = B(INCRJ) - TOPBLK(I,INCRJ)*B(I)
 850                CONTINUE
           B(INCRJ) = B(INCRJ)/TOPBLK(INCRJ,INCRJ)
 900          CONTINUE
C
C       RETURN FROM THE SOLUTION OF A-TRANSPOSE.X = B

        RETURN
        END

c ===================================================================================

      subroutine lufac(n, ndim, a, ip, ier)
      implicit double precision (a-h,o-z)
      dimension a(ndim,n), ip(n)
      intrinsic abs
*  blas: daxpy, dscal, dswap. idamax

      parameter ( zero = 0.0d+0, one = 1.0d+0 )

*  The subroutine lufac is a very simple code to compute the
*  LU decomposition (with partial pivoting) of the n by n matrix a.

*  The LU factors are overwritten on a.  The integer array ip
*  reflects the pairwise interchanges performed.  note that ip(k)
*  therefore does not give the index in the original array of
*  the k-th pivot.

*  On exit, the error flag ier is zero when no zero pivots are
*  encountered.  Otherwise, ier is equal to the index of the
*  step at which a zero pivot occurred.

      ier = 0
      ip(n) = 0

*  Begin loop over columns 1 through n-1.  k is the current
*  column index.

      do 100 k = 1, n-1

*  Find the row index ipiv of the element of largest magnitude in
*  column k.

         ipiv = k-1 + idamax(n-k+1, a(k,k), 1)
         piv = a(ipiv,k)
         if (piv .eq. zero) then
            ier = k
            return
         endif
         ip(k) = ipiv

*  Perform interchanges if necessary.

         if (ipiv .ne. k) then
            call dswap(n-k+1, a(ipiv,k), ndim, a(k,k), ndim)
         endif

*  Save the (negative) multipliers in the subdiagonal elements of
*  column k.

         call dscal(n-k, (-one/piv), a(k+1,k), 1)

*  Update the remaining matrix.  Note that a(i,k) now contains
*  the negative multipliers.

         do 50 j = k+1, n
            call daxpy(n-k, a(k,j), a(k+1,k), 1, a(k+1,j), 1)
   50    continue

*  End of loop over columns.

  100 continue
      if (a(n,n).eq.zero) ier = n
      return
      end

      subroutine lusol (n, ndim, a, ip, b, x)
      implicit double precision (a-h, o-z)
      dimension a(ndim,n), ip(n), b(n), x(n)

*  blas:  daxpy, dcopy

*  The subroutine lusol is a simple-minded routine to solve a
*  linear system whose LU factors have been computed by lufac.
*  On entry, the matrix a should contain the LU factors, and
*  ip should contain the interchange array constructed by lufac.


*  Copy the right-hand side b into x.

      call dcopy (n, b, 1, x, 1)

*  Forward solution with l (unit lower-triangular factor), which
*  is stored in the strict lower triangle of a.

      do 20 k = 1, n-1
         ipiv = ip(k)
         if (ipiv .ne. k) then
            tem = x(ipiv)
            x(ipiv) = x(k)
            x(k) = tem
         endif
         call daxpy ( n-k, x(k), a(k+1,k), 1, x(k+1), 1 )
   20 continue

*  Backward solution with u (upper-triangular factor), which is stored
*  in the upper triangle of a.

      do 40 kb = n, 1, -1
         x(kb) = x(kb)/a(kb,kb)
         call daxpy(kb-1, (-x(kb)), a(1,kb), 1, x(1), 1)
   40 continue

      return
      end

c ===================================================================================

      subroutine dcopy ( n, x, incx, y, incy )
      implicit double precision (a-h,o-z)
      integer            n, incx, incy
      dimension          x( * ), y( * )

c  dcopy  performs the operation
c
c     y := x
c
c  nag fortran 77 version of the blas routine dcopy .
c  nag fortran 77 o( n ) basic linear algebra routine.
c
c  -- written on 26-november-1982.
c     sven hammarling, nag central office.

      integer            i     , ix    , iy

      if( n.lt.1 )return

      if( ( incx.eq.incy ).and.( incy.gt.0 ) )then
         do 10, iy = 1, 1 + ( n - 1 )*incy, incy
            y( iy ) = x( iy )
   10    continue
      else
         if( incx.ge.0 )then
            ix = 1
         else
            ix = 1 - ( n - 1 )*incx
         end if
         if( incy.gt.0 )then
            do 20, iy = 1, 1 + ( n - 1 )*incy, incy
               y( iy ) = x( ix )
               ix      = ix + incx
   20       continue
         else
            iy = 1 - ( n - 1 )*incy
            do 30, i = 1, n
               y( iy ) = x( ix )
               iy      = iy + incy
               ix      = ix + incx
   30       continue
         end if
      end if

      return

*     end of dcopy .

      end

c ===================================================================================

      subroutine daxpy ( n, alpha, x, incx, y, incy )
      implicit double precision (a-h,o-z)
      integer            n, incx, incy
      dimension   x( * ), y( * )

c  daxpy  performs the operation
c
c     y := alpha*x + y
c
c
c  modified nag fortran 77 version of the blas routine daxpy .
c
c  -- written on 3-september-1982.
c     sven hammarling, nag central office.

      integer            i     , ix    , iy
      parameter        ( zero  = 0.0+0 )

      if( n    .lt.1    )return
      if( alpha.eq.zero )return

      if( ( incx.eq.incy ).and.( incx.gt.0 ) )then
         do 10, ix = 1, 1 + ( n - 1 )*incx, incx
            y( ix ) = alpha*x( ix ) + y( ix )
   10    continue
      else
         if( incy.ge.0 )then
            iy = 1
         else
            iy = 1 - ( n - 1 )*incy
         end if
         if( incx.gt.0 )then
            do 20, ix = 1, 1 + ( n - 1 )*incx, incx
               y( iy ) = alpha*x( ix ) + y( iy )
               iy      = iy + incy
   20       continue
         else
            ix = 1 - ( n - 1 )*incx
            do 30, i = 1, n
               y( iy ) = alpha*x( ix ) + y( iy )
               ix      = ix + incx
               iy      = iy + incy
   30       continue
         end if
      end if
      return

*     end of daxpy .

      end

c ===================================================================================

      double precision function ddot  ( n, x, incx, y, incy )
      implicit double precision (a-h,o-z)
      integer           n, incx, incy
      dimension         x( * ), y( * )

c  ddot   returns the value
c
c     ddot   = x'y
c
c
c  modified nag fortran 77 version of the blas routine ddot  .
c
c  -- written on 21-september-1982.
c     sven hammarling, nag central office.

      integer             i     , ix    , iy
      parameter         ( zero  = 0.0d+0 )

      sum = zero
      if( n.ge.1 )then
         if( ( incx.eq.incy ).and.( incx.gt.0 ) )then
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx
               sum = sum + x( ix )*y( ix )
   10       continue
         else
            if( incy.ge.0 )then
               iy = 1
            else
               iy = 1 - ( n - 1 )*incy
            end if
            if( incx.gt.0 )then
               do 20, ix = 1, 1 + ( n - 1 )*incx, incx
                  sum = sum + x( ix )*y( iy )
                  iy  = iy  + incy
   20          continue
            else
               ix = 1 - ( n - 1 )*incx
               do 30, i = 1, n
                  sum = sum + x( ix )*y( iy )
                  ix  = ix  + incx
                  iy  = iy  + incy
   30          continue
            end if
         end if
      end if

      ddot   = sum
      return

*     end of ddot  .

      end

c ===================================================================================

      subroutine dscal ( n, alpha, x, incx )
      implicit double precision (a-h,o-z)
      integer          n, incx
      dimension        x( * )

c  dscal  performs the operation
c
c     x := alpha*x
c
c
c  modified nag fortran 77 version of the blas routine dscal .
c
c  -- written on 26-november-1982.
c     sven hammarling, nag central office.

      integer            ix
      parameter        ( one   = 1.0d+0, zero  = 0.0d+0 )

      if( n.ge.1 )then
         if( alpha.eq.zero )then
            do 10, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = zero
   10       continue
         else if( alpha.eq.( -one ) )then
            do 20, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = -x( ix )
   20       continue
         else if( alpha.ne.one )then
            do 30, ix = 1, 1 + ( n - 1 )*incx, incx
               x( ix ) = alpha*x( ix )
   30       continue
         end if
      end if

      return

*     end of dscal .

      end

c ===================================================================================

      subroutine dswap ( n, x, incx, y, incy )
      implicit double precision (a-h,o-z)
      integer      n, incx, incy
      dimension    x( * ), y( * )

c  dswap  performs the operations
c
c     temp := x,   x := y,   y := temp.
c
c
c  modified nag fortran 77 version of the blas routine dswap .
c
c  -- written on 26-november-1982.
c     sven hammarling, nag central office.

      integer            i     , ix    , iy

      if( n.lt.1 )return

      if( ( incx.eq.incy ).and.( incy.gt.0 ) )then
         do 10, iy = 1, 1 + ( n - 1 )*incy, incy
            temp    = x( iy )
            x( iy ) = y( iy )
            y( iy ) = temp
   10    continue
      else
         if( incx.ge.0 )then
            ix = 1
         else
            ix = 1 - ( n - 1 )*incx
         end if
         if( incy.gt.0 )then
            do 20, iy = 1, 1 + ( n - 1 )*incy, incy
               temp    = x( ix )
               x( ix ) = y( iy )
               y( iy ) = temp
               ix      = ix + incx
   20       continue
         else
            iy = 1 - ( n - 1 )*incy
            do 30, i = 1, n
               temp    = x( ix )
               x( ix ) = y( iy )
               y( iy ) = temp
               iy      = iy + incy
               ix      = ix + incx
   30       continue
         end if
      end if

      return

*     end of dswap .

      end

c ===================================================================================

      integer function idamax( n, x, incx )
      implicit double precision (a-h,o-z)
      integer         n, incx
      dimension       x( * )

c  idamax returns the smallest value of i such that
c
c     abs( x( i ) ) = max( abs( x( j ) ) )
c                      j
c
c  nag fortran 77 version of the blas routine idamax.
c  nag fortran 77 o( n ) basic linear algebra routine.
c
c  -- written on 31-may-1983.
c     sven hammarling, nag central office.

      intrinsic           abs
      integer             i     , imax  , ix

      if( n.lt.1 )then
         idamax = 0
         return
      end if

      imax = 1
      if( n.gt.1 )then
         xmax = abs( x( 1 ) )
         ix   = 1
         do 10, i = 2, n
            ix = ix + incx
            if( xmax.lt.abs( x( ix ) ) )then
               xmax = abs( x( ix ) )
               imax = i
            end if
   10    continue
      end if

      idamax = imax
      return

*     end of idamax.

      end

c ===================================================================================

      subroutine dload ( n, const, x, incx )
      implicit double precision (a-h,o-z)
      dimension  x(*)
c
c  dload  performs the operation
c
c     x = const*e,   e' = ( 1  1 ... 1 ).
c
c
c  nag fortran 77 o( n ) basic linear algebra routine.
c
c  -- written on 22-september-1983.
c     sven hammarling, nag central office.
c
      parameter        ( zero = 0.0d+0 )

      if( n.lt.1 )return

      if( const.ne.zero )then
         do 10, ix = 1, 1 + ( n - 1 )*incx, incx
            x( ix ) = const
   10    continue
      else
         do 20, ix = 1, 1 + ( n - 1 )*incx, incx
            x( ix ) = zero
   20    continue
      end if

      return

*     end of dload .
      end

c ===================================================================================

      subroutine maxpy ( nrow, ncol, alpha, xmat, nrowy, ymat )
      implicit double precision (a-h,o-z)
      dimension xmat(nrow, ncol), ymat(nrowy, ncol)

*  Subroutine maxpy takes as input the scalar alpha and two matrices,
*  xmat and ymat.  xmat has declared row dimension nrow, and
*  ymat has declared row dimension nrowy, but both are
*  conceptually nrow by ncol.
*  On output, (new ymat) is alpha*xmat+ (old ymat), by analogy
*  with the vector blas routine saxpy.

      do 100 j = 1, ncol
      do 100 i = 1, nrow
         ymat(i,j) = ymat(i,j) + alpha*xmat(i,j)
  100 continue
      return
      end

c ===================================================================================

      subroutine matcop( nrow1, nrow2, nrow, ncol, xmat1, xmat2 )
      implicit double precision (a-h,o-z)
      dimension xmat1(nrow1, ncol), xmat2(nrow2, ncol)

*  Given 2 matrices xmat1 and xmat2, where xmat1 has declared
*  row dimension nrow1, xmat2 has declared row dimension nrow2,
*  and both have column dimension ncol, the routine matcop copies
*  rows 1 through nrow, and columns 1 through ncol from xmat1 into
*  xmat2.


      if (nrow .le. 0 .or. ncol .le. 0) return
      do 100 j = 1, ncol
      do 100 i = 1, nrow
           xmat2(i,j) = xmat1(i,j)
  100 continue
      return
      end

c ===================================================================================

      subroutine mtload( nrow, ncol, const, nrowx, xmat )
      implicit double precision (a-h,o-z)
      dimension xmat(nrowx, ncol)

*  mtload sets elements 1 through nrow, 1 through ncol, of the
*  matrix xmat (whose declared row dimension is nrowx) to the
*  scalar value const.

      if (nrow .le. 0 .or. ncol .le. 0)  return
      do 100 j = 1, ncol
      do 100 i = 1, nrow
         xmat(i,j) = const
  100 continue
      return
      end

c ===================================================================================

      subroutine mssq  ( nrow, ncol, xmat, scale, sumsq )
      implicit double precision (a-h,o-z)
      dimension   xmat(nrow, *)

*  Given the nrow by ncol matrix xmat, mssq returns values
*  scale and sumsq such that
*     (scale**2) * sumsq = sum of squares of xmat(i,j),
*  where  scale = max  abs(xmat(i,j)).

*  mssq is a stripped-down matrix version of the blas routine sssq.

      intrinsic    abs
      parameter   ( one   = 1.0d+0, zero  = 0.0d+0 )

      scale = zero
      sumsq = one
      if( nrow.ge.1 .and. ncol .ge. 1) then
         do 10 i = 1, nrow
         do 10 j = 1, ncol
            if( xmat(i,j) .ne. zero )then
               absxij = abs(xmat(i,j))
               if( scale .lt. absxij ) then
                  sumsq = one + sumsq* (scale/absxij)**2
                  scale = absxij
               else
                  sumsq = sumsq + (absxij/scale)**2
               end if
            end if
   10    continue
      end if
      return
      end

c ===================================================================================

      subroutine dssq  ( n, x, incx, scale, sumsq )
      implicit double precision (a-h,o-z)
      integer            n, incx
      dimension   x( * )

*  Given the n-vector x, dssq returns values scale and sumsq such that
*     (scale**2) * sumsq = sum of squares of x(i),
*  where  scale = max  abs(x(i)).

      intrinsic          abs
      parameter        ( one   = 1.0d+0, zero  = 0.0d+0 )

      scale = zero
      sumsq = one
      if( n.ge.1 )then
         do 10, ix = 1, 1 + ( n - 1 )*incx, incx
            if( x( ix ).ne.zero )then
               absxi = abs( x( ix ) )
               if( scale.lt.absxi )then
                  sumsq = one   + sumsq*( scale/absxi )**2
                  scale = absxi
               else
                  sumsq = sumsq + ( absxi/scale )**2
               end if
            end if
   10    continue
      end if

      return

*     end of dssq  .

      end




