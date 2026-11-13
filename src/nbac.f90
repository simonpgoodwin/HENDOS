! ==============================================================
   MODULE nbac
!
! HENDOS - HErmite N-body Dynamics Of Systems
! AHMAD-COHEN NEIGHBOUR SCHEME (switch useac) for direct summation
! (Ahmad & Cohen 1973, J. Comp. Phys., 12, 389; Hermite form after
! Makino & Aarseth 1992, PASJ, 44, 141)
!
! The force on a star is split in two:
!   irregular force   from its neighbours (the stars inside its
!                     neighbour sphere, radius rs), which changes fast.
!                     The Galactic tide is included in it.
!   regular force     from all the other stars, which changes slowly
! Each star has two timesteps. On most steps (irregular steps) only the
! short sum over the neighbour list is done and the regular force is
! taken from its last value and derivative (linear extrapolation). On
! regular steps, which are 2, 4, 8... irregular steps apart, the full
! sum over all stars is done, a new neighbour list is made and the
! regular force is corrected.
!
! How it fits the Hermite scheme. The arrays f0, f1 are still the total
! force and its derivative, so the predictor and the corrector in the
! main program are unchanged:
!   irregular step   f0p = (neighbour sum + tide) + fr + frd*(t - treg)
!                    f1p = (its derivative) + frd
!                    the regular part is linear in time, so it drops
!                    out of the corrector's 2nd and 3rd derivatives,
!                    which are those of the irregular force alone
!   regular step     the same, with the OLD neighbour list, then the
!                    Hermite corrector is applied to the regular force
!                    over the whole regular step: its 2nd and 3rd
!                    derivatives r2, r3 give the extra corrections
!                      dreg  = r2 dtr^4/24 + r3 dtr^5/120   to the position
!                      dvreg = r2 dtr^3/6  + r3 dtr^4/24    to the velocity
!                    which the main program adds. Then the new list
!                    replaces the old one and the force is split again.
! The regular step is chosen from the regular force and its derivatives
! with the Aarseth formula and the parameter etar. It is a block step
! like the irregular one, never shorter than it, and regular times are
! always whole multiples of the regular step, so every regular time is
! also one of the star's irregular times.
!
! Neighbour list: stars inside rs, plus stars inside 3 rs that are
! approaching and will be inside rs within one regular step. rs is
! changed a little at every regular step to keep the number of
! neighbours near nnbopt. kfmax places are kept free in every list for
! stars that appear when a group is opened (see below).
!
! When bodies are created, removed or moved by the binary and merger
! routines the lists are kept correct by acmoved, acremoved and
! acadded; a new body has acok = .false. and gets a full force sum and
! a new list the first time forces is called for it.
! (Tested with a Python model of this scheme before it was written:
! a 150-star Plummer sphere over 4 N-body times gave |dE/E| = 3e-8 with
! etar = 0.005 and a quarter to a half of the full force sums.)
!
! OpenMP: for big blocks the active stars are shared over the threads;
! for small blocks the full sums and list searches are shared instead.
!
   USE nbvars
   IMPLICIT NONE
   DOUBLE PRECISION :: rs0=0.d0          ! first guess of the neighbour radius
!
   CONTAINS
! =====================================================================
   SUBROUTINE acalloc
!
! choose the list sizes if they were left automatic, and allocate.
! Warns if there are fewer than 100 stars: the scheme then costs as
! much as plain direct summation and is less accurate (in a Python
! model at N = 30 it was no faster and about 7 times less accurate).
!
   IF (nalloc < 100) THEN
     WRITE(6,'(a)') ''
     WRITE(6,'(a,i0,a)') ' WARNING: the neighbour scheme (useac) is on for only ', &
&      nalloc, ' stars.'
     WRITE(6,'(a)') '          Below about 100 stars it is no faster than plain '// &
&      'direct summation'
     WRITE(6,'(a)') '          and is less accurate: set useac = .false. unless '// &
&      'you are testing it.'
     WRITE(6,'(a)') ''
   END IF
   IF (nnbopt == 0) nnbopt=MAX(10,NINT((DBLE(nalloc)/8.d0)**0.75d0))
   nnbopt=MIN(nnbopt,nalloc - 1)
   IF (nnbmax == 0) nnbmax=3*nnbopt + kfmax + 10
   IF (nnbmax < nnbopt + kfmax + 2) THEN
     WRITE(6,'(a)') ' ERROR: nnbmax must be at least nnbopt + 8'
     STOP 1
   END IF
   ALLOCATE(fi(1:3,1:nalloc),fid(1:3,1:nalloc),fr(1:3,1:nalloc),frd(1:3,1:nalloc))
   ALLOCATE(fin(1:3,1:nalloc),fidn(1:3,1:nalloc),frn(1:3,1:nalloc),frdn(1:3,1:nalloc))
   ALLOCATE(dreg(1:3,1:nalloc),dvreg(1:3,1:nalloc),rs(1:nalloc),dtwr(1:nalloc))
   ALLOCATE(nnb(1:nalloc),levr(1:nalloc),nbl(1:nnbmax,1:nalloc))
   ALLOCATE(treg(1:nalloc),dtreg(1:nalloc),lreg(1:nalloc),acok(1:nalloc))
   acok=.FALSE.
   lreg=.FALSE.
   nnb=0
   rs=0.d0
   treg=0
   dtreg=0
   levr=0
!
   END SUBROUTINE acalloc
! =====================================================================
   SUBROUTINE acsetrs0
!
! first guess of the neighbour radius: the radius of a sphere holding
! nnbopt stars at the mean density inside the rms radius
!
   INTEGER :: i
   DOUBLE PRECISION :: s,mt
!
   s=0.d0
   mt=0.d0
   DO i=1,n
     s=s + m(i)*SUM(r(1:3,i)**2)
     mt=mt + m(i)
   END DO
   rs0=SQRT(s/mt)*(DBLE(nnbopt)/DBLE(n))**(1.d0/3.d0)
!
   END SUBROUTINE acsetrs0
! =====================================================================
   SUBROUTINE acforce(i,inpar)
!
! force and its derivative on active star i at the block time tcur8,
! into f0p, f1p, by an irregular step, a regular step or (for a body
! with acok false) a first full sum. inpar is true when this is called
! from inside an OpenMP parallel loop over the active stars.
! Thread safe: only star i's entries are written.
!
   INTEGER, INTENT(IN) :: i
   LOGICAL, INTENT(IN) :: inpar
   LOGICAL :: par,init
   INTEGER :: nl
   INTEGER, ALLOCATABLE :: newl(:)
   DOUBLE PRECISION :: an(1:3),jn(1:3),at(1:3),jt(1:3),ao(1:3),jo(1:3)
   DOUBLE PRECISION :: a(1:3),j(1:3),tn2,rn2,ta2,ra2,dtr,den
   DOUBLE PRECISION :: fro(1:3),frdo(1:3),a1(1:3),a2(1:3),r2(1:3),r3(1:3),r2e(1:3)
!
   par=(.NOT. inpar) .AND. n >= nompj
   init=.NOT. acok(i)
!
   IF (.NOT. init) THEN
     IF (tcur8 /= treg(i) + dtreg(i)) THEN
! ---------------- irregular step: neighbours only
       CALL nbsum(i,nnb(i),nbl(:,i),an,jn,tn2,rn2)
       dtr=DBLE(tcur8 - treg(i))*tick
       a=an + fr(1:3,i) + frd(1:3,i)*dtr
       j=jn + frd(1:3,i)
       CALL actide(i,a,j)
       fin(1:3,i)=a - fr(1:3,i) - frd(1:3,i)*dtr
       fidn(1:3,i)=j - frd(1:3,i)
       f0p(1:3,i)=a
       f1p(1:3,i)=j
       tpair2(i)=tn2
       rmin2(i)=rn2
       lreg(i)=.FALSE.
!$OMP ATOMIC
       nirregular=nirregular + 1
       RETURN
     END IF
   END IF
!
! ---------------- regular step or first sum: all the stars
   CALL fullsum(i,par,at,jt,ta2,ra2)
   tpair2(i)=ta2
   rmin2(i)=ra2
!$OMP ATOMIC
   nregular=nregular + 1
!
   IF (.NOT. init) THEN
! the end of the irregular step, with the old list
     dtr=DBLE(dtreg(i))*tick
     CALL nbsum(i,nnb(i),nbl(:,i),ao,jo,tn2,rn2)
     a=ao + fr(1:3,i) + frd(1:3,i)*dtr
     j=jo + frd(1:3,i)
     CALL actide(i,a,j)
     f0p(1:3,i)=a
     f1p(1:3,i)=j
! Hermite corrector of the regular force (old list) over the regular step
     fro=at - ao
     frdo=jt - jo
     a1=fr(1:3,i) - fro
     a2=frd(1:3,i) + frdo
     r2=2.d0*(-3.d0*a1 - (frd(1:3,i) + a2)*dtr)/dtr**2
     r3=6.d0*(2.d0*a1 + a2*dtr)/dtr**3
     dreg(1:3,i)=r2*dtr**4/24.d0 + r3*dtr**5/120.d0
     dvreg(1:3,i)=r2*dtr**3/6.d0 + r3*dtr**4/24.d0
! regular timestep wanted (Aarseth formula on the regular force)
     r2e=r2 + r3*dtr
     den=SQRT(SUM(frdo**2))*SQRT(SUM(r3**2)) + SUM(r2e**2)
     IF (den > 0.d0) THEN
       dtwr(i)=SQRT(etar*(SQRT(SUM(fro**2))*SQRT(SUM(r2e**2)) + SUM(frdo**2))/den)
     ELSE
       dtwr(i)=dtmax
     END IF
     lreg(i)=.TRUE.
! with corrector iterations the new list is only made on the last one
     IF (.NOT. aclast) RETURN
   END IF
!
! new neighbour list, and the force split with it
   ALLOCATE(newl(1:nnbmax))
   CALL aclist(i,par,init,nl,newl)
   CALL nbsum(i,nl,newl,an,jn,tn2,rn2)
   a=at
   j=jt
   CALL actide(i,a,j)
   frn(1:3,i)=at - an
   frdn(1:3,i)=jt - jn
   fin(1:3,i)=a - frn(1:3,i)
   fidn(1:3,i)=j - frdn(1:3,i)
   nnb(i)=nl
   nbl(1:nl,i)=newl(1:nl)
   DEALLOCATE(newl)
!
   IF (init) THEN
! a new body: everything is set now (levr, dtreg are set by acstart)
     fi(1:3,i)=fin(1:3,i)
     fid(1:3,i)=fidn(1:3,i)
     fr(1:3,i)=frn(1:3,i)
     frd(1:3,i)=frdn(1:3,i)
     f0p(1:3,i)=a
     f1p(1:3,i)=j
     treg(i)=tcur8
     lreg(i)=.FALSE.
     acok(i)=.TRUE.
   END IF
!
   END SUBROUTINE acforce
! =====================================================================
   SUBROUTINE pairsum(i,j,xi,vi,mi,a,jk,tmin2,rm2)
!
! add the acceleration and jerk on star i (at xi, vi, mass mi) from
! star j to a, jk, and update the shortest pairwise timescale squared
! and nearest distance squared. The same pair interaction as in forces
! (module nbforce) - the two must be kept identical.
!
   INTEGER, INTENT(IN) :: i,j
   DOUBLE PRECISION, INTENT(IN) :: xi(1:3),vi(1:3),mi
   DOUBLE PRECISION, INTENT(INOUT) :: a(1:3),jk(1:3),tmin2,rm2
   DOUBLE PRECISION :: dx,dy,dz,dvx,dvy,dvz,d2u,rad2,rinv,rinv2,mor,drdv,dv2
!
   IF (i == j) RETURN
   dx=rp(1,j) - xi(1)
   dy=rp(2,j) - xi(2)
   dz=rp(3,j) - xi(3)
   dvx=vp(1,j) - vi(1)
   dvy=vp(2,j) - vi(2)
   dvz=vp(3,j) - vi(3)
   d2u=dx*dx + dy*dy + dz*dz
   rm2=MIN(rm2,d2u)
   rad2=d2u + eps2
   rinv=1.d0/SQRT(rad2)
   rinv2=rinv*rinv
   mor=m(j)*rinv*rinv2
   drdv=3.d0*(dx*dvx + dy*dvy + dz*dvz)*rinv2
   a(1)=a(1) + mor*dx
   a(2)=a(2) + mor*dy
   a(3)=a(3) + mor*dz
   jk(1)=jk(1) + mor*(dvx - drdv*dx)
   jk(2)=jk(2) + mor*(dvy - drdv*dy)
   jk(3)=jk(3) + mor*(dvz - drdv*dz)
   dv2=dvx*dvx + dvy*dvy + dvz*dvz
   tmin2=MIN(tmin2, rad2*rad2*rinv/(mi + m(j)), rad2/(dv2 + TINY(1.d0)))
!
   END SUBROUTINE pairsum
! =====================================================================
   SUBROUTINE nbsum(i,nl,list,a,jk,tmin2,rm2)
!
! acceleration and jerk on star i from the nl stars in list
!
   INTEGER, INTENT(IN) :: i,nl,list(:)
   DOUBLE PRECISION, INTENT(OUT) :: a(1:3),jk(1:3),tmin2,rm2
   INTEGER :: k
   DOUBLE PRECISION :: xi(1:3),vi(1:3)
!
   xi=rp(1:3,i)
   vi=vp(1:3,i)
   a=0.d0
   jk=0.d0
   tmin2=HUGE(1.d0)
   rm2=HUGE(1.d0)
   DO k=1,nl
     CALL pairsum(i,list(k),xi,vi,m(i),a,jk,tmin2,rm2)
   END DO
!
   END SUBROUTINE nbsum
! =====================================================================
   SUBROUTINE fullsum(i,par,a,jk,tmin2,rm2)
!
! acceleration and jerk on star i from all the other stars; the sum is
! shared over the threads if par is true
!
   INTEGER, INTENT(IN) :: i
   LOGICAL, INTENT(IN) :: par
   DOUBLE PRECISION, INTENT(OUT) :: a(1:3),jk(1:3),tmin2,rm2
   INTEGER :: j
   DOUBLE PRECISION :: xi(1:3),vi(1:3),mi
!
   xi=rp(1:3,i)
   vi=vp(1:3,i)
   mi=m(i)
   a=0.d0
   jk=0.d0
   tmin2=HUGE(1.d0)
   rm2=HUGE(1.d0)
!$OMP PARALLEL DO IF(par) SCHEDULE(STATIC) PRIVATE(j) &
!$OMP& REDUCTION(+:a,jk) REDUCTION(MIN:tmin2,rm2)
   DO j=1,n
     CALL pairsum(i,j,xi,vi,mi,a,jk,tmin2,rm2)
   END DO
!$OMP END PARALLEL DO
!
   END SUBROUTINE fullsum
! =====================================================================
   SUBROUTINE actide(i,a,jk)
!
! add the Galactic tide to the total acceleration a and jerk jk of
! star i (the same as finish in forces): a = T.x + 2 omega (vy,-vx,0),
! jerk = T.v + 2 omega (ay,-ax,0) with the TOTAL acceleration
!
   INTEGER, INTENT(IN) :: i
   DOUBLE PRECISION, INTENT(INOUT) :: a(1:3),jk(1:3)
!
   IF (.NOT. ltide) RETURN
   a=a + MATMUL(ttid,rp(1:3,i))
   a(1)=a(1) + 2.d0*omega*vp(2,i)
   a(2)=a(2) - 2.d0*omega*vp(1,i)
   jk=jk + MATMUL(ttid,vp(1:3,i))
   jk(1)=jk(1) + 2.d0*omega*a(2)
   jk(2)=jk(2) - 2.d0*omega*a(1)
!
   END SUBROUTINE actide
! =====================================================================
   SUBROUTINE aclist(i,par,init,nl,newl)
!
! new neighbour list of star i: the stars inside rs(i), and those
! inside 3 rs that are approaching and will be inside rs within one
! regular step. If there are too many, rs is reduced and the search is
! done again; for a new body (init) rs is also increased if there are
! too few. Afterwards rs is changed by at most 25% towards the radius
! that would hold nnbopt neighbours, for the next regular step.
! The search is shared over the threads if par is true (the stars are
! cut into one piece per thread, as for the active list).
!
   INTEGER, INTENT(IN) :: i
   LOGICAL, INTENT(IN) :: par,init
   INTEGER, INTENT(OUT) :: nl,newl(:)
   INTEGER :: it,ic,j,k,nper,nlim,nfound,cnt(1:nthreads)
   INTEGER, ALLOCATABLE :: scr(:)
   DOUBLE PRECISION :: rs2,dtr,dx,dy,dz,dvx,dvy,dvz,d2,rv,fut,third
!
   third=1.d0/3.d0
   nlim=nnbmax - kfmax
   IF (rs(i) <= 0.d0) rs(i)=rs0
   dtr=0.d0
   IF (.NOT. init) dtr=DBLE(dtreg(i))*tick
   nper=(n + nthreads - 1)/nthreads
   ALLOCATE(scr(1:n))
!
   DO it=1,30
     rs2=rs(i)**2
!$OMP PARALLEL DO IF(par) SCHEDULE(STATIC,1) &
!$OMP& PRIVATE(ic,j,k,dx,dy,dz,dvx,dvy,dvz,d2,rv,fut)
     DO ic=1,nthreads
       k=(ic - 1)*nper
       DO j=(ic - 1)*nper + 1,MIN(n,ic*nper)
         IF (j == i) CYCLE
         dx=rp(1,j) - rp(1,i)
         dy=rp(2,j) - rp(2,i)
         dz=rp(3,j) - rp(3,i)
         d2=dx*dx + dy*dy + dz*dz
         IF (d2 >= rs2) THEN
           IF (init .OR. d2 >= 9.d0*rs2) CYCLE
           dvx=vp(1,j) - vp(1,i)
           dvy=vp(2,j) - vp(2,i)
           dvz=vp(3,j) - vp(3,i)
           rv=dx*dvx + dy*dvy + dz*dvz
           IF (rv >= 0.d0) CYCLE
           fut=(dx + dvx*dtr)**2 + (dy + dvy*dtr)**2 + (dz + dvz*dtr)**2
           IF (fut >= rs2) CYCLE
         END IF
         k=k + 1
         scr(k)=j
       END DO
       cnt(ic)=k - (ic - 1)*nper
     END DO
!$OMP END PARALLEL DO
     nfound=SUM(cnt)
     IF (it == 30) EXIT
     IF (nfound > nlim) THEN
       rs(i)=rs(i)*MAX(0.5d0,(DBLE(nnbopt)/DBLE(nfound))**third)
       CYCLE
     END IF
     IF (init .AND. nfound < nnbopt/2 .AND. nfound < n - 1) THEN
       rs(i)=rs(i)*MIN(2.d0,(DBLE(nnbopt)/DBLE(MAX(nfound,1)))**third)
       CYCLE
     END IF
     EXIT
   END DO
!
! join the pieces up (never more than nlim, whatever happened above)
   nl=0
   DO ic=1,nthreads
     DO k=1,cnt(ic)
       IF (nl >= nlim) EXIT
       nl=nl + 1
       newl(nl)=scr((ic - 1)*nper + k)
     END DO
   END DO
   DEALLOCATE(scr)
!
   rs(i)=rs(i)*MIN(1.25d0,MAX(0.8d0,(DBLE(nnbopt)/DBLE(MAX(nl,1)))**third))
!
   END SUBROUTINE aclist
! =====================================================================
   SUBROUTINE acmoverow(src,dst)
!
! copy the neighbour-scheme data of body src into slot dst (called by
! moverow in module nbbin)
!
   INTEGER, INTENT(IN) :: src,dst
!
   fi(1:3,dst)=fi(1:3,src)
   fid(1:3,dst)=fid(1:3,src)
   fr(1:3,dst)=fr(1:3,src)
   frd(1:3,dst)=frd(1:3,src)
   rs(dst)=rs(src)
   nnb(dst)=nnb(src)
   nbl(1:nnb(src),dst)=nbl(1:nnb(src),src)
   levr(dst)=levr(src)
   treg(dst)=treg(src)
   dtreg(dst)=dtreg(src)
   acok(dst)=acok(src)
   lreg(dst)=.FALSE.
!
   END SUBROUTINE acmoverow
! =====================================================================
   SUBROUTINE acmoved(src,dst)
!
! body src has been moved to slot dst: change its index in every list
!
   INTEGER, INTENT(IN) :: src,dst
   INTEGER :: i,k
!
!$OMP PARALLEL DO IF(n >= nompj) PRIVATE(i,k) SCHEDULE(STATIC)
   DO i=1,n
     DO k=1,nnb(i)
       IF (nbl(k,i) == src) nbl(k,i)=dst
     END DO
   END DO
!$OMP END PARALLEL DO
!
   END SUBROUTINE acmoved
! =====================================================================
   SUBROUTINE acremoved(j,irep)
!
! body j is being removed (it has joined a group or merged, and body
! irep now stands for it; irep = 0 if nothing does): take j out of
! every list. In a list that held j but not irep, irep takes j's place.
!
   INTEGER, INTENT(IN) :: j,irep
   INTEGER :: i,k,kj
   LOGICAL :: has
!
!$OMP PARALLEL DO IF(n >= nompj) PRIVATE(i,k,kj,has) SCHEDULE(STATIC)
   DO i=1,n
     kj=0
     has=.FALSE.
     DO k=1,nnb(i)
       IF (nbl(k,i) == j) kj=k
       IF (nbl(k,i) == irep) has=.TRUE.
     END DO
     IF (kj == 0) CYCLE
     IF (irep > 0 .AND. irep /= i .AND. .NOT. has) THEN
       nbl(kj,i)=irep
     ELSE
       nbl(kj,i)=nbl(nnb(i),i)
       nnb(i)=nnb(i) - 1
     END IF
   END DO
!$OMP END PARALLEL DO
!
   END SUBROUTINE acremoved
! =====================================================================
   SUBROUTINE acadded(jn,i0)
!
! body jn has appeared beside body i0 (a group was opened): add jn to
! every list that holds i0, so its neighbours feel it at once. Lists
! keep kfmax places free for this; if a list is full all the same, jn
! is left out until that star's next regular step.
!
   INTEGER, INTENT(IN) :: jn,i0
   INTEGER :: i,k
   LOGICAL :: has
!
!$OMP PARALLEL DO IF(n >= nompj) PRIVATE(i,k,has) SCHEDULE(STATIC)
   DO i=1,n
     IF (i == jn .OR. i == i0) CYCLE
     has=.FALSE.
     DO k=1,nnb(i)
       IF (nbl(k,i) == i0) has=.TRUE.
     END DO
     IF (has .AND. nnb(i) < nnbmax) THEN
       nnb(i)=nnb(i) + 1
       nbl(nnb(i),i)=jn
     END IF
   END DO
!$OMP END PARALLEL DO
!
   END SUBROUTINE acadded
!
   END MODULE nbac
