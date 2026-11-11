! ==============================================================
   MODULE nbbin
!
! HENDOS - HErmite N-body Dynamics Of Systems
! UNPERTURBED BINARIES (no regularisation)
!
! A bound pair that is effectively unperturbed is replaced by a single
! centre-of-mass (c.o.m.) body: the internal orbit is then an exact
! Kepler orbit, advanced analytically, and the pair no longer forces
! tiny timesteps on the integration. When it becomes perturbed again
! it is 'opened' back into its two stars, which are then integrated
! directly with the Hermite scheme like any other stars.
!
! relative perturbation of a binary of mass mb and apocentre rapo by
! all other bodies k at distance d_k from its c.o.m.:
!    gamma = 2 rapo^3 / mb * sum_k m_k / d_k^3
! (tidal acceleration across the orbit / internal acceleration)
!
! CLOSE  star i after its step if: its step < dtclose, its nearest
!        neighbour j is a single star on the same block, the pair is
!        bound, and gamma < gclose
! OPEN   a c.o.m. body after each of its steps if gamma > gopen
!        (gopen > gclose gives hysteresis, so binaries do not flicker)
! Binaries in the initial conditions are closed at the start.
!
! The bodies are kept contiguous in slots 1..n of the particle arrays
! (n = number of bodies, ntot = number of stars), so the force loops
! are unchanged: closing puts the c.o.m. in one slot and moves the last
! body into the other; opening uses the c.o.m. slot and slot n+1.
! Only one open and one close are done per block (they are rare).
!
! All of this is done serially between blocks, outside any OpenMP
! parallel region; the O(N) loops inside are parallelised. With the
! tree on, the perturbation sum and nearest-neighbour search use it.
!
! Binaries are invisible in the output: snapshots contain all ntot
! stars, with closed binaries split into their two stars.
! The particle types of the two stars are kept in btype while the
! binary is closed (the c.o.m. body has itype = 0) and restored when
! it is opened.
!
   USE nbvars
   USE nbtree
   USE nbforce
   USE nbstep
   IMPLICIT NONE
!
   CONTAINS
! =====================================================================
   SUBROUTINE kepler(mu,r0,v0,dt,r1,v1)
!
! advance a bound two-body relative orbit (G(m1+m2) = mu) by dt with the
! f and g functions and the eccentric anomaly difference (Danby 1988).
! Exact, and dt is reduced modulo the period, so any number of orbits
! can be done in one call.
!
   DOUBLE PRECISION, INTENT(IN) :: mu,r0(1:3),v0(1:3),dt
   DOUBLE PRECISION, INTENT(OUT) :: r1(1:3),v1(1:3)
   DOUBLE PRECISION :: r0n,r1n,alpha,a,nmo,sig0,c1,c2,dm,dtr,x,xn,lo,hi
   DOUBLE PRECISION :: f,g,fd,gd,ce,se,ff,fp,twopi
   INTEGER :: it
!
   twopi=6.283185307179586d0
   r0n=SQRT(SUM(r0*r0))
   alpha=2.d0/r0n - SUM(v0*v0)/mu          ! 1/a
   a=1.d0/alpha
   nmo=SQRT(mu*alpha**3)                   ! mean motion
   sig0=SUM(r0*v0)/SQRT(mu)
   c1=1.d0 - r0n/a
   c2=sig0/SQRT(a)
   dm=MODULO(nmo*dt,twopi)
   dtr=dm/nmo
! Kepler's equation in difference form, dm = x - c1 sin x + c2 (1 - cos x)
! the root is within 1 of dm; Newton's method kept inside the bracket
   lo=dm - 1.d0
   hi=dm + 1.d0
   x=dm
   DO it=1,60
     ff=x - c1*SIN(x) + c2*(1.d0 - COS(x)) - dm
     fp=1.d0 - c1*COS(x) + c2*SIN(x)
     IF (ff > 0.d0) THEN
       hi=x
     ELSE
       lo=x
     END IF
     xn=x - ff/fp
     IF (xn <= lo .OR. xn >= hi) xn=0.5d0*(lo + hi)
     IF (ABS(xn - x) < 1.d-15) THEN
       x=xn
       EXIT
     END IF
     x=xn
   END DO
   ce=COS(x)
   se=SIN(x)
   f=1.d0 - a/r0n*(1.d0 - ce)
   g=dtr + (se - x)/nmo
   r1=f*r0 + g*v0
   r1n=SQRT(SUM(r1*r1))
   fd=-SQRT(mu*a)*se/(r1n*r0n)
   gd=1.d0 - a/r1n*(1.d0 - ce)
   v1=fd*r0 + gd*v0
!
   END SUBROUTINE kepler
! =====================================================================
   SUBROUTINE bincomp(k,xc,vc,tnow,x1,u1,x2,u2)
!
! positions and velocities of the two stars of binary k at internal
! time tnow, given its c.o.m. position xc and velocity vc
!
   INTEGER, INTENT(IN) :: k
   DOUBLE PRECISION, INTENT(IN) :: xc(1:3),vc(1:3),tnow
   DOUBLE PRECISION, INTENT(OUT) :: x1(1:3),u1(1:3),x2(1:3),u2(1:3)
   DOUBLE PRECISION :: mb,rr(1:3),vr(1:3)
!
   mb=bmass(1,k) + bmass(2,k)
   CALL kepler(mb,brel(1:3,k),brel(4:6,k),tnow - bepoch(k),rr,vr)
   x1=xc + bmass(2,k)/mb*rr
   u1=vc + bmass(2,k)/mb*vr
   x2=xc - bmass(1,k)/mb*rr
   u2=vc - bmass(1,k)/mb*vr
!
   END SUBROUTINE bincomp
! =====================================================================
   SUBROUTINE pairorbit(mb,rr,vr,ebind,rapo)
!
! specific orbital energy and apocentre of a relative orbit
! (rapo is huge if unbound)
!
   DOUBLE PRECISION, INTENT(IN) :: mb,rr(1:3),vr(1:3)
   DOUBLE PRECISION, INTENT(OUT) :: ebind,rapo
   DOUBLE PRECISION :: rn,a,h(1:3),ecc2
!
   rn=SQRT(SUM(rr*rr))
   ebind=0.5d0*SUM(vr*vr) - mb/rn
   IF (ebind >= 0.d0) THEN
     rapo=HUGE(1.d0)
     RETURN
   END IF
   a=-0.5d0*mb/ebind
   h(1)=rr(2)*vr(3) - rr(3)*vr(2)
   h(2)=rr(3)*vr(1) - rr(1)*vr(3)
   h(3)=rr(1)*vr(2) - rr(2)*vr(1)
   ecc2=MAX(0.d0,1.d0 - SUM(h*h)/(mb*a))
   rapo=a*(1.d0 + SQRT(ecc2))
!
   END SUBROUTINE pairorbit
! =====================================================================
   DOUBLE PRECISION FUNCTION pericentre(mb,rr,vr)
!
! two-body pericentre distance of a relative orbit (bound or unbound),
! q = h^2 / (mb (1 + e))
!
   DOUBLE PRECISION, INTENT(IN) :: mb,rr(1:3),vr(1:3)
   DOUBLE PRECISION :: h(1:3),h2,ebind,ecc
!
   h(1)=rr(2)*vr(3) - rr(3)*vr(2)
   h(2)=rr(3)*vr(1) - rr(1)*vr(3)
   h(3)=rr(1)*vr(2) - rr(2)*vr(1)
   h2=SUM(h*h)
   ebind=0.5d0*SUM(vr*vr) - mb/SQRT(SUM(rr*rr))
   ecc=SQRT(MAX(0.d0,1.d0 + 2.d0*ebind*h2/mb**2))
   pericentre=h2/(mb*(1.d0 + ecc))
!
   END FUNCTION pericentre
! =====================================================================
   SUBROUTINE collision(i,j,rsep,q,tnow)
!
! two stars have come (or are about to come) within rcoll: stop the
! run with an alert. Collisions are not modelled yet.
!
   INTEGER, INTENT(IN) :: i,j
   DOUBLE PRECISION, INTENT(IN) :: rsep,q,tnow
   DOUBLE PRECISION :: rsun
!
   rsun=6.957d8/3.0856775814913673d16       ! solar radius in pc
   WRITE(6,'(a)') ''
   WRITE(6,'(a)') REPEAT('*',78)
   WRITE(6,'(a,f0.6,a)') '  STELLAR COLLISION at t = ', tnow*tunit, ' Myr'
   WRITE(6,'(a)') REPEAT('*',78)
   CALL collstar(i)
   CALL collstar(j)
   WRITE(6,'(a,es12.4,a)') '  separation now          ', rsep/rsun, ' Rsun'
   WRITE(6,'(a,es12.4,a)') '  two-body pericentre     ', q/rsun, ' Rsun'
   WRITE(6,'(a,es12.4,a)') '  collision radius rcoll  ', rcoll, ' Rsun'
   WRITE(6,'(a,3f14.6,a)') '  position                ', rp(1:3,i), ' pc'
   WRITE(6,'(a)') '  The run has been stopped: stellar collisions are not '// &
&   'modelled yet.'
   WRITE(6,'(a)') REPEAT('*',78)
   STOP 1
!
   CONTAINS
   SUBROUTINE collstar(l)
   INTEGER, INTENT(IN) :: l
   IF (ib(l) == 0) THEN
     WRITE(6,'(a,i12,a,f12.5,a,a)') '  star  id ', id(l), '   mass ', m(l), &
&      ' Msun   type ', TRIM(typename(itype(l)))
   ELSE
     WRITE(6,'(a,2i12,a,2f12.5,a,a,1x,a)') '  binary ids ', bident(1:2,ib(l)), &
&      '   masses ', bmass(1:2,ib(l)), ' Msun   types ', &
&      TRIM(typename(btype(1,ib(l)))), TRIM(typename(btype(2,ib(l))))
   END IF
   END SUBROUTINE collstar
!
   END SUBROUTINE collision
! =====================================================================
   SUBROUTINE collcheck(tnow)
!
! after each block: for each active body that came within 2 rcoll of
! another at its force evaluation, find that neighbour and test the
! pair: a collision if they are already within rcoll, or approaching
! with a two-body pericentre inside rcoll. (The flyby timestep limit
! means sampled separations near pericentre are within ~20% of the
! true pericentre, so 2 rcoll is a safe trigger.) Serial; rare.
!
   DOUBLE PRECISION, INTENT(IN) :: tnow
   INTEGER :: ia,i,j
   DOUBLE PRECISION :: rr(1:3),vr(1:3),rsep,q
!
   DO ia=1,nact
     i=ilist(ia)
     IF (rmin2(i) >= 4.d0*rcolli*rcolli) CYCLE
     CALL nearest(i,j)
     IF (j == 0) CYCLE
     rr=rp(1:3,i) - rp(1:3,j)
     vr=vp(1:3,i) - vp(1:3,j)
     rsep=SQRT(SUM(rr*rr))
     q=pericentre(m(i) + m(j),rr,vr)
     IF (rsep < rcolli .OR. (q < rcolli .AND. SUM(rr*vr) < 0.d0)) &
&      CALL collision(i,j,rsep,q,tnow)
   END DO
!
   END SUBROUTINE collcheck
! =====================================================================
   SUBROUTINE perturb(xc,vc,mb,rapo,i1,i2,gam,tp2)
!
! perturbation gamma on a (possible) binary of mass mb and apocentre
! rapo at xc from all bodies except i1 and i2, and the squared shortest
! pairwise free-fall/flyby time of the c.o.m. (for its timestep).
! Uses rp, vp: every body's position predicted to the current block time
!
   DOUBLE PRECISION, INTENT(IN) :: xc(1:3),vc(1:3),mb,rapo
   INTEGER, INTENT(IN) :: i1,i2
   DOUBLE PRECISION, INTENT(OUT) :: gam,tp2
   INTEGER :: j
   DOUBLE PRECISION :: sum3,d2,dv2,dx,dy,dz
!
   IF (usetree .AND. treevalid) THEN
     CALL treeperturb(xc,vc,mb,i1,i2,sum3,tp2)
     gam=2.d0*rapo**3/mb*sum3
     RETURN
   END IF
   sum3=0.d0
   tp2=HUGE(1.d0)
!$OMP PARALLEL DO IF(n >= nompj) PRIVATE(j,d2,dv2,dx,dy,dz) &
!$OMP& REDUCTION(+:sum3) REDUCTION(MIN:tp2) SCHEDULE(STATIC)
   DO j=1,n
     IF (j == i1 .OR. j == i2) CYCLE
     dx=rp(1,j) - xc(1)
     dy=rp(2,j) - xc(2)
     dz=rp(3,j) - xc(3)
     d2=dx*dx + dy*dy + dz*dz
     dv2=(vp(1,j) - vc(1))**2 + (vp(2,j) - vc(2))**2 + (vp(3,j) - vc(3))**2
     sum3=sum3 + m(j)/(d2*SQRT(d2))
     tp2=MIN(tp2, d2*SQRT(d2)/(mb + m(j)), d2/(dv2 + TINY(1.d0)))
   END DO
!$OMP END PARALLEL DO
   gam=2.d0*rapo**3/mb*sum3
!
   END SUBROUTINE perturb
! =====================================================================
   SUBROUTINE nearest(i,jn)
!
! index jn of the body nearest to body i (positions rp at the block time)
!
   INTEGER, INTENT(IN) :: i
   INTEGER, INTENT(OUT) :: jn
   INTEGER :: j,jbest
   DOUBLE PRECISION :: d2,dbest,gbest
!
   IF (usetree .AND. treevalid) THEN
     CALL treenearest(i,jn)
     RETURN
   END IF
   jn=0
   gbest=HUGE(1.d0)
!$OMP PARALLEL IF(n >= nompj) PRIVATE(j,d2,dbest,jbest)
   dbest=HUGE(1.d0)
   jbest=0
!$OMP DO SCHEDULE(STATIC)
   DO j=1,n
     IF (j == i) CYCLE
     d2=(rp(1,j) - rp(1,i))**2 + (rp(2,j) - rp(2,i))**2 + (rp(3,j) - rp(3,i))**2
     IF (d2 < dbest) THEN
       dbest=d2
       jbest=j
     END IF
   END DO
!$OMP END DO
!$OMP CRITICAL (nearestcrit)
   IF (dbest < gbest .OR. (dbest == gbest .AND. jbest < jn)) THEN
     gbest=dbest
     jn=jbest
   END IF
!$OMP END CRITICAL (nearestcrit)
!$OMP END PARALLEL
!
   END SUBROUTINE nearest
! =====================================================================
   SUBROUTINE moverow(src,dst)
!
! copy body src into slot dst
!
   INTEGER, INTENT(IN) :: src,dst
!
   r(1:3,dst)=r(1:3,src)
   v(1:3,dst)=v(1:3,src)
   rp(1:3,dst)=rp(1:3,src)
   vp(1:3,dst)=vp(1:3,src)
   f0(1:3,dst)=f0(1:3,src)
   f1(1:3,dst)=f1(1:3,src)
   f2(1:3,dst)=f2(1:3,src)
   f3(1:3,dst)=f3(1:3,src)
   f2e(1:3,dst)=f2e(1:3,src)
   f3e(1:3,dst)=f3e(1:3,src)
   m(dst)=m(src)
   id(dst)=id(src)
   itype(dst)=itype(src)
   ib(dst)=ib(src)
   lev(dst)=lev(src)
   tk(dst)=tk(src)
   dtk(dst)=dtk(src)
   tpair2(dst)=tpair2(src)
   rmin2(dst)=rmin2(src)
   phi(dst)=phi(src)
!
   END SUBROUTINE moverow
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
   SUBROUTINE closebin(i,j,t8,tnow)
!
! replace single stars i and j (both at integer time t8 = internal time
! tnow) by their c.o.m. body
!
   INTEGER, INTENT(IN) :: i,j
   INTEGER(KIND=8), INTENT(IN) :: t8
   DOUBLE PRECISION, INTENT(IN) :: tnow
   INTEGER :: k,lo,hi,knew
   DOUBLE PRECISION :: mi,mj,mb,xc(1:3),vc(1:3),gam,tp2,dtw,amod,jmod
   DOUBLE PRECISION :: ebind,rapo
!
   mi=m(i)
   mj=m(j)
   mb=mi + mj
   k=nbin + 1
   nbin=k
   bmass(1,k)=mi
   bmass(2,k)=mj
   bident(1,k)=id(i)
   bident(2,k)=id(j)
   btype(1,k)=itype(i)
   btype(2,k)=itype(j)
   brel(1:3,k)=r(1:3,i) - r(1:3,j)
   brel(4:6,k)=v(1:3,i) - v(1:3,j)
   bepoch(k)=tnow
!
   xc=(mi*r(1:3,i) + mj*r(1:3,j))/mb
   vc=(mi*v(1:3,i) + mj*v(1:3,j))/mb
   lo=MIN(i,j)
   hi=MAX(i,j)
! c.o.m. derivatives: the mass-weighted mean, in which the mutual
! forces of the pair cancel exactly
   f0(1:3,lo)=(mi*f0(1:3,i) + mj*f0(1:3,j))/mb
   f1(1:3,lo)=(mi*f1(1:3,i) + mj*f1(1:3,j))/mb
   f2e(1:3,lo)=(mi*f2e(1:3,i) + mj*f2e(1:3,j))/mb
   f3e(1:3,lo)=(mi*f3e(1:3,i) + mj*f3e(1:3,j))/mb
   r(1:3,lo)=xc
   v(1:3,lo)=vc
   rp(1:3,lo)=xc
   vp(1:3,lo)=vc
   m(lo)=mb
   id(lo)=0
   itype(lo)=0
   ib(lo)=k
   tk(lo)=t8
!
! timestep of the c.o.m.
   CALL pairorbit(mb,brel(1:3,k),brel(4:6,k),ebind,rapo)
   CALL perturb(xc,vc,mb,rapo,i,j,gam,tp2)
   tpair2(lo)=tp2
   dtw=etap*SQRT(tp2)
   amod=SQRT(SUM(f0(1:3,lo)**2))
   jmod=SQRT(SUM(f1(1:3,lo)**2))
   IF (amod > 0.d0 .AND. jmod > 0.d0) dtw=MIN(dtw,etas*amod/jmod)
   CALL getlevel(dtw,knew)
   lev(lo)=commlev(knew,t8)
   dtk(lo)=ISHFT(1_8,kmax - lev(lo))
!
! remove slot hi by moving the last body into it
   IF (hi /= n) CALL moverow(n,hi)
   n=n - 1
   nclosed=nclosed + 1
   treevalid=.FALSE.
!
   END SUBROUTINE closebin
! =====================================================================
   SUBROUTINE openbin(i,t8,tnow)
!
! split the c.o.m. body i (at integer time t8 = internal time tnow)
! back into its two stars, in slot i and slot n+1, and start them off
! with forces from all other bodies (whose rp, vp are at t8)
!
   INTEGER, INTENT(IN) :: i
   INTEGER(KIND=8), INTENT(IN) :: t8
   DOUBLE PRECISION, INTENT(IN) :: tnow
   INTEGER :: k,j,levc,knew,ia,l
   DOUBLE PRECISION :: x1(1:3),u1(1:3),x2(1:3),u2(1:3),dtw,amod,jmod
!
   k=ib(i)
   levc=lev(i)
   CALL bincomp(k,r(1:3,i),v(1:3,i),tnow,x1,u1,x2,u2)
!
   j=n + 1
   n=j
   r(1:3,i)=x1
   v(1:3,i)=u1
   m(i)=bmass(1,k)
   id(i)=bident(1,k)
   itype(i)=btype(1,k)
   ib(i)=0
   r(1:3,j)=x2
   v(1:3,j)=u2
   m(j)=bmass(2,k)
   id(j)=bident(2,k)
   itype(j)=btype(2,k)
   ib(j)=0
!
! forces on the two stars
   rp(1:3,i)=x1
   vp(1:3,i)=u1
   rp(1:3,j)=x2
   vp(1:3,j)=u2
   IF (usetree) CALL treebuild
   nact=2
   ilist(1)=i
   ilist(2)=j
   CALL forces
   DO ia=1,2
     l=ilist(ia)
     f0(1:3,l)=f0p(1:3,l)
     f1(1:3,l)=f1p(1:3,l)
     f2(1:3,l)=0.d0
     f3(1:3,l)=0.d0
     f2e(1:3,l)=0.d0
     f3e(1:3,l)=0.d0
     tk(l)=t8
     amod=SQRT(SUM(f0(1:3,l)**2))
     jmod=SQRT(SUM(f1(1:3,l)**2))
     dtw=etap*SQRT(tpair2(l))
     IF (amod > 0.d0 .AND. jmod > 0.d0) dtw=MIN(dtw,etas*amod/jmod)
     CALL getlevel(dtw,knew)
! at least the c.o.m. level, so the step divides t8
     lev(l)=MAX(knew,levc)
     dtk(l)=ISHFT(1_8,kmax - lev(l))
   END DO
!
! remove binary k from the table, moving the last binary into its place
   IF (k /= nbin) THEN
     bmass(1:2,k)=bmass(1:2,nbin)
     bident(1:2,k)=bident(1:2,nbin)
     btype(1:2,k)=btype(1:2,nbin)
     brel(1:6,k)=brel(1:6,nbin)
     bepoch(k)=bepoch(nbin)
     DO l=1,n
       IF (ib(l) == nbin) THEN
         ib(l)=k
         EXIT
       END IF
     END DO
   END IF
   nbin=nbin - 1
   nopened=nopened + 1
!
   END SUBROUTINE openbin
! =====================================================================
   LOGICAL FUNCTION closeable(i,j,tnowc)
!
! true if single stars i and j form a bound pair with gamma < gclose
! (internal time tnowc is only used for a collision alert)
!
   INTEGER, INTENT(IN) :: i,j
   DOUBLE PRECISION, INTENT(IN) :: tnowc
   DOUBLE PRECISION :: mb,rr(1:3),vr(1:3),xc(1:3),vc(1:3),ebind,rapo,gam,tp2
!
   closeable=.FALSE.
   IF (ib(i) /= 0 .OR. ib(j) /= 0) RETURN
   mb=m(i) + m(j)
   rr=r(1:3,i) - r(1:3,j)
   vr=v(1:3,i) - v(1:3,j)
   CALL pairorbit(mb,rr,vr,ebind,rapo)
   IF (ebind >= 0.d0) RETURN
! the orbit of a closed binary is fixed, so check its pericentre now
   IF (rcoll > 0.d0) THEN
     IF (pericentre(mb,rr,vr) < rcolli) CALL collision(i,j, &
&      SQRT(SUM(rr*rr)),pericentre(mb,rr,vr),tnowc)
   END IF
   xc=(m(i)*r(1:3,i) + m(j)*r(1:3,j))/mb
   vc=(m(i)*v(1:3,i) + m(j)*v(1:3,j))/mb
   CALL perturb(xc,vc,mb,rapo,i,j,gam,tp2)
   closeable=(gam < gclose)
!
   END FUNCTION closeable
! =====================================================================
   SUBROUTINE bincheck(t8,tnow)
!
! called after each block (all active bodies at integer time t8 =
! internal time tnow): open at most one perturbed binary and close at
! most one unperturbed pair among the bodies that have just stepped
!
   INTEGER(KIND=8), INTENT(IN) :: t8
   DOUBLE PRECISION, INTENT(IN) :: tnow
   INTEGER :: ia,i,j,k,nact0
   INTEGER, ALLOCATABLE :: alist(:)
   DOUBLE PRECISION :: mb,gam,tp2,ebind,rapo,rr(1:3),vr(1:3),x1(1:3),u1(1:3)
   DOUBLE PRECISION :: x2(1:3),u2(1:3)
!
! copy the active list (forces in openbin overwrites ilist)
   nact0=nact
   ALLOCATE(alist(1:nact0))
   alist=ilist(1:nact0)
!
! opening
   DO ia=1,nact0
     i=alist(ia)
     k=ib(i)
     IF (k == 0) CYCLE
     mb=m(i)
     CALL bincomp(k,r(1:3,i),v(1:3,i),tnow,x1,u1,x2,u2)
     rr=x1 - x2
     vr=u1 - u2
     CALL pairorbit(mb,rr,vr,ebind,rapo)
     CALL perturb(r(1:3,i),v(1:3,i),mb,rapo,i,i,gam,tp2)
     IF (gam > gopen) THEN
       CALL openbin(i,t8,tnow)
       EXIT
     END IF
   END DO
!
! closing - only stars with short steps are candidates
   DO ia=1,nact0
     i=alist(ia)
     IF (i > n) CYCLE
     IF (ib(i) /= 0 .OR. tk(i) /= t8) CYCLE
     IF (DBLE(dtk(i))*tick >= dtclosei) CYCLE
     CALL nearest(i,j)
     IF (j == 0) CYCLE
     IF (tk(j) /= t8) CYCLE
     IF (closeable(i,j,tnow)) THEN
       CALL closebin(i,j,t8,tnow)
       EXIT
     END IF
   END DO
   DEALLOCATE(alist)
!
   END SUBROUTINE bincheck
! =====================================================================
   SUBROUTINE binstart(tnow)
!
! close all unperturbed binaries in the initial conditions (everyone is
! synchronised at integer time 0 = internal time tnow). After a close,
! slot i may hold a different (moved) body, so i only advances when
! nothing was closed.
!
   DOUBLE PRECISION, INTENT(IN) :: tnow
   INTEGER :: i,j
!
   i=1
   DO WHILE (i <= n)
     IF (ib(i) == 0) THEN
       CALL nearest(i,j)
       IF (j > 0) THEN
         IF (closeable(i,j,tnow)) THEN
           CALL closebin(i,j,0_8,tnow)
           CYCLE
         END IF
       END IF
     END IF
     i=i + 1
   END DO
!
   END SUBROUTINE binstart
!
   END MODULE nbbin
