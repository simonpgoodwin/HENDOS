! ==============================================================
   MODULE nbstep
!
! HENDOS - HErmite N-body Dynamics Of Systems
! TIMESTEPS
! starting forces and timestep levels, and the choice of each
! particle's block timestep level after every step
!
   USE nbvars
   USE nbtree
   USE nbac
   USE nbforce
   IMPLICIT NONE
!
   CONTAINS
! =====================================================================
   SUBROUTINE firststep
!
! forces on everyone and starting levels
! dt = min(etas*|a|/|j|, etap*tpair) - the |a|/|j| term is skipped for
! particles with zero acceleration (e.g. a star at a symmetry centre)
!
   INTEGER :: i
   DOUBLE PRECISION :: amod,jmod,dtw
!
! everyone is active, forces from the initial positions
   nact=n
   DO i=1,n
     ilist(i)=i
   END DO
   rp=r
   vp=v
   tcur8=0
   aclast=.TRUE.
   IF (useac) CALL acsetrs0
   IF (usetree) CALL treebuild
   CALL forces
   f0=f0p
   f1=f1p
   f2=0.d0
   f3=0.d0
   f2e=0.d0
   f3e=0.d0
   tk=0
!
!$OMP PARALLEL DO IF(n >= nompj) PRIVATE(i,amod,jmod,dtw) SCHEDULE(STATIC)
   DO i=1,n
     amod=SQRT(f0(1,i)**2 + f0(2,i)**2 + f0(3,i)**2)
     jmod=SQRT(f1(1,i)**2 + f1(2,i)**2 + f1(3,i)**2)
     dtw=etap*SQRT(tpair2(i))
     IF (amod > 0.d0 .AND. jmod > 0.d0) dtw=MIN(dtw,etas*amod/jmod)
     CALL getlevel(dtw,lev(i))
     dtk(i)=ISHFT(1_8,kmax - lev(i))
     IF (useac) CALL acstart(i)
   END DO
!$OMP END PARALLEL DO
!
   END SUBROUTINE firststep
! =====================================================================
   SUBROUTINE newlevel(i)
!
! new timestep level for particle i after its step has been corrected
!
! dt_A = sqrt(eta*(|a||a2| + |a1|^2)/(|a1||a3| + |a2|^2))  Aarseth
!   with a2 and a3 at the END of the step: a2 = f2 + f3*dt, a3 = f3
! dt_P = etap*min_j(sqrt(r^3/(mi+mj)), r/v)               pairwise
!   catches approaching encounters before a and its derivatives react
!
! levels can go down (smaller step) by any amount, but up by only one
! level and only when the current time is commensurate with 2*dt
!
   INTEGER, INTENT(IN) :: i
   INTEGER :: knew
   DOUBLE PRECISION :: amod,jmod,smod,cmod,s(1:3),dt,dtw,denom
!
   dt=DBLE(dtk(i))*tick
   s(1:3)=f2(1:3,i) + f3(1:3,i)*dt
   amod=SQRT(f0(1,i)**2 + f0(2,i)**2 + f0(3,i)**2)
   jmod=SQRT(f1(1,i)**2 + f1(2,i)**2 + f1(3,i)**2)
   smod=SQRT(s(1)**2 + s(2)**2 + s(3)**2)
   cmod=SQRT(f3(1,i)**2 + f3(2,i)**2 + f3(3,i)**2)
!
   dtw=etap*SQRT(tpair2(i))
   denom=jmod*cmod + smod*smod
   IF (denom > 0.d0) dtw=MIN(dtw,SQRT(eta*(amod*smod + jmod*jmod)/denom))
!
   CALL getlevel(dtw,knew)
!
   IF (knew > lev(i)) THEN
     lev(i)=knew
   ELSE IF (knew < lev(i)) THEN
     IF (MOD(tk(i),2_8*dtk(i)) == 0) lev(i)=lev(i) - 1
   END IF
   dtk(i)=ISHFT(1_8,kmax - lev(i))
!
   END SUBROUTINE newlevel
! =====================================================================
   INTEGER FUNCTION commlev(k,t8)
!
! smallest level >= k whose step divides the integer time t8, so a body
! at t8 can be given it without breaking the block structure
!
   INTEGER, INTENT(IN) :: k
   INTEGER(KIND=8), INTENT(IN) :: t8
!
   commlev=MAX(0,MIN(k,kmax))
   DO WHILE (commlev < kmax)
     IF (MOD(t8,ISHFT(1_8,kmax - commlev)) == 0) EXIT
     commlev=commlev + 1
   END DO
!
   END FUNCTION commlev
! =====================================================================
   SUBROUTINE acstart(i)
!
! neighbour scheme: start the regular step of a new body i, whose
! irregular level lev(i) and time tk(i) have just been set - the first
! regular step is the same as the irregular step
!
   INTEGER, INTENT(IN) :: i
!
   levr(i)=lev(i)
   dtreg(i)=dtk(i)
   treg(i)=tk(i)
!
   END SUBROUTINE acstart
! =====================================================================
   SUBROUTINE acaccept(i)
!
! neighbour scheme: the step of star i has been accepted (the main
! program has set f0, f1 from f0p, f1p and tk(i) to the new time).
! Keep the new irregular force, and after a regular step the new
! regular force too; f0, f1 are then the exact totals from the full sum
! rather than the extrapolated ones the corrector used.
! Called before newlevel.
!
   INTEGER, INTENT(IN) :: i
!
   fi(1:3,i)=fin(1:3,i)
   fid(1:3,i)=fidn(1:3,i)
   IF (lreg(i)) THEN
     fr(1:3,i)=frn(1:3,i)
     frd(1:3,i)=frdn(1:3,i)
     f0(1:3,i)=fin(1:3,i) + frn(1:3,i)
     f1(1:3,i)=fidn(1:3,i) + frdn(1:3,i)
     treg(i)=tk(i)
   END IF
!
   END SUBROUTINE acaccept
! =====================================================================
   SUBROUTINE aclevel(i)
!
! neighbour scheme: after newlevel has chosen the irregular level of
! star i, choose its regular level if this was a regular step, and
! make sure the irregular step is no longer than the regular step.
! The regular level can go down (smaller step) by any amount, up by one
! level at a time, and its step must divide the current time so that
! every regular time is also an irregular time of the star.
!
   INTEGER, INTENT(IN) :: i
   INTEGER :: kr
!
   IF (lreg(i)) THEN
     CALL getlevel(dtwr(i),kr)
     IF (kr < levr(i)) kr=levr(i) - 1
     levr(i)=commlev(kr,tk(i))
     dtreg(i)=ISHFT(1_8,kmax - levr(i))
     lreg(i)=.FALSE.
   END IF
   IF (lev(i) < levr(i)) THEN
     lev(i)=levr(i)
     dtk(i)=ISHFT(1_8,kmax - lev(i))
   END IF
!
   END SUBROUTINE aclevel
! =====================================================================
   SUBROUTINE getlevel(dtw,k)
!
! smallest level k with dtmax/2**k <= dtw, limited to 0 <= k <= kmax
!
   DOUBLE PRECISION, INTENT(IN) :: dtw
   INTEGER, INTENT(OUT) :: k
!
   IF (dtw >= dtmax) THEN
     k=0
   ELSE
     k=CEILING(LOG(dtmax/dtw)/LOG(2.d0))
   END IF
   IF (k > kmax) THEN
! steps that want to be shorter than the shortest level are held there
! and counted (unperturbed binaries are dealt with in module nbbin)
     k=kmax
!$OMP ATOMIC
     nclamp=nclamp + 1
   END IF
!
   END SUBROUTINE getlevel
!
   END MODULE nbstep
