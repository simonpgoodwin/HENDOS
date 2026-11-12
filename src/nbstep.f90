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
