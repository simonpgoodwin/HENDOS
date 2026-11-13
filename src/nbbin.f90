! ==============================================================
   MODULE nbbin
!
! HENDOS - HErmite N-body Dynamics Of Systems
! BINARIES, FEW-BODY GROUPS, COLLISIONS AND MERGERS (no regularisation)
!
! GROUPS. A compact set of stars that is effectively unperturbed by the
! rest of the cluster is replaced by a single centre-of-mass (c.o.m.)
! body, and its internal motion is followed separately, in isolation.
! The stars then no longer force tiny timesteps on the main
! integration (where every step costs a prediction of all N bodies).
! There are three kinds of group (bmode):
!   1  binary          internal orbit is an exact Kepler orbit
!   2  few-body group  3 to kfmax stars in a close encounter, integrated
!                      by the small Hermite integrator in module nbfew
!                      (compare the 'multiples' module of AMUSE). The
!                      group and the cluster advance together in time:
!                      the group is brought up to date whenever its
!                      c.o.m. body steps and whenever output is written
!   3  stable triple   a hierarchical triple that passes the Mardling &
!                      Aarseth (2001) stability test is frozen as two
!                      nested Kepler orbits (no secular evolution)
! When a group becomes perturbed it is 'opened' back into its stars,
! which are then integrated directly like any other stars.
!
! relative perturbation of a group of mass mb and gsize s by all other
! bodies k at distance d_k from its c.o.m.:
!    gamma = 2 s^3 / mb * sum_k m_k / d_k^3
! (tidal acceleration across the group / internal acceleration);
! s is the apocentre of a binary, or the diameter of a larger group
!
! CLOSE   star i after its step if: its step < dtclose, its nearest
!         neighbour j is a single star on the same block, the pair is
!         bound, and gamma < gclose                       -> binary
! FEW     if that pair is bound but perturbed, the pair and the next
!         nearest body are tried as a group; so are a perturbed group
!         and its nearest body. The group is formed if its own
!         gamma < gfew, and if there are enough bodies in the main
!         integration for it to pay (see tryfew)          -> few-body
! FREEZE  a few-body group of three that is a stable hierarchy
!                                                         -> triple
! OPEN    a binary after one of its c.o.m. steps if gamma > gopen, any
!         other group if gamma > 10 gfew (unless it can take its
!         perturber in and still be isolated, see FEW)
! (gopen > gclose and the factor 10 give hysteresis, so groups do not
! flicker.) Binaries in the initial conditions are closed at the start.
! A star that leaves a few-body group makes the group large, so its
! gamma rises and it is opened; what is left is closed again as a
! binary if it is one.
!
! COLLISIONS AND MERGERS. Two stars that pass within rcoll of each
! other have collided. With mergers on they are merged into one star
! (mass and momentum conserved) and the run carries on; a merger is
! announced on the screen and written to <runid>_mergers.dat. The
! orbital energy of the pair is lost from the stellar dynamics: it is
! added up in emerged so that the energy check still works. With
! mergers off the run stops with an alert, as before. Collisions inside
! few-body groups are treated in the same way.
!
! The bodies are kept contiguous in slots 1..n of the particle arrays
! (n = number of bodies, ntot = number of stars), so the force loops
! are unchanged: closing puts the c.o.m. in the lowest slot and moves
! the last bodies into the others; opening uses the c.o.m. slot and
! slots n+1 onwards. Only one open and one close are done per block.
!
! All of this is done serially between blocks, outside any OpenMP
! parallel region; the O(N) loops inside are parallelised. With the
! tree on, the perturbation sum and nearest-neighbour search use it.
! With the neighbour scheme on (module nbac) the neighbour lists are
! kept correct through every change.
!
! Groups are invisible in the output: snapshots contain all ntot
! stars. The types of the stars are kept in btype while they are in a
! group (the c.o.m. body has itype = 0) and restored when it is opened.
!
   USE nbvars
   USE nbtree
   USE nbac
   USE nbforce
   USE nbstep
   USE nbfew
   IMPLICIT NONE
! a few-body group of nm stars is only taken out of the main
! integration if that has at least nfewfac*nm*(nm-1) bodies (see tryfew)
   INTEGER, PARAMETER :: nfewfac=2
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
   SUBROUTINE grpcomp(k,xc,vc,tnow,nm,xm,um)
!
! positions xm and velocities um of the nm stars of group k at internal
! time tnow, given its c.o.m. position xc and velocity vc.
! For a few-body group this brings the group's own integration up to
! tnow, so it must not be called from inside a parallel region.
!
   INTEGER, INTENT(IN) :: k
   DOUBLE PRECISION, INTENT(IN) :: xc(1:3),vc(1:3),tnow
   INTEGER, INTENT(OUT) :: nm
   DOUBLE PRECISION, INTENT(OUT) :: xm(1:3,1:kfmax),um(1:3,1:kfmax)
   INTEGER :: l
   DOUBLE PRECISION :: mb,mt,rr(1:3),vr(1:3),ro(1:3),vo(1:3),xi(1:3),vi(1:3)
!
   xm=0.d0
   um=0.d0
   SELECT CASE (bmode(k))
   CASE (1)
     CALL bincomp(k,xc,vc,tnow,xm(1:3,1),um(1:3,1),xm(1:3,2),um(1:3,2))
   CASE (3)
! inner orbit (stars 1, 2) and outer orbit (star 3 about their c.o.m.)
     mb=bmass(1,k) + bmass(2,k)
     mt=mb + bmass(3,k)
     CALL kepler(mb,brel(1:3,k),brel(4:6,k),tnow - bepoch(k),rr,vr)
     CALL kepler(mt,brel2(1:3,k),brel2(4:6,k),tnow - bepoch(k),ro,vo)
     xi=xc - bmass(3,k)/mt*ro
     vi=vc - bmass(3,k)/mt*vo
     xm(1:3,1)=xi + bmass(2,k)/mb*rr
     um(1:3,1)=vi + bmass(2,k)/mb*vr
     xm(1:3,2)=xi - bmass(1,k)/mb*rr
     um(1:3,2)=vi - bmass(1,k)/mb*vr
     xm(1:3,3)=xc + mb/mt*ro
     um(1:3,3)=vc + mb/mt*vo
   CASE DEFAULT
     CALL fewadvance(k,xc,vc,tnow)
     DO l=1,bnum(k)
       xm(1:3,l)=xc + bx(1:3,l,k)
       um(1:3,l)=vc + bv(1:3,l,k)
     END DO
   END SELECT
   nm=bnum(k)
!
   END SUBROUTINE grpcomp
! =====================================================================
   SUBROUTINE fewadvance(k,xc,vc,tnow)
!
! bring the internal integration of few-body group k up to internal
! time tnow (xc, vc are its c.o.m. position and velocity, only used to
! report a merger). If two of its stars collide on the way they are
! merged inside the group (or the run is stopped, if mergers are off)
! and the integration carries on with one star fewer.
!
   INTEGER, INTENT(IN) :: k
   DOUBLE PRECISION, INTENT(IN) :: xc(1:3),vc(1:3),tnow
   INTEGER :: nm,ic,jc,idn,tn
   INTEGER(KIND=8) :: nst
   DOUBLE PRECISION :: trem,tdone,rc2,m1,m2,mn,rr(1:3),vr(1:3),rsep
   DOUBLE PRECISION :: xn(1:3),vn(1:3)
   LOGICAL, SAVE :: warned=.FALSE.
!
   rc2=0.d0
   IF (rcoll > 0.d0) rc2=rcolli*rcolli
   DO
     trem=tnow - bepoch(k)
     nm=bnum(k)
     IF (trem <= 0.d0 .OR. nm < 2) EXIT
     CALL fewint(nm,bmass(1:nm,k),bx(1:3,1:nm,k),bv(1:3,1:nm,k),trem,etafew, &
&      eps2,rc2,tdone,ic,jc,nst)
     nfewsteps=nfewsteps + nst
     IF (ic == 0) THEN
       bepoch(k)=tnow
       EXIT
     END IF
     bepoch(k)=bepoch(k) + tdone
     IF (ic < 0) THEN
! very many steps: carry on, but say so once
       IF (.NOT. warned) WRITE(6,'(a,f0.4,a)') '  ** NOTE t = ', tnow*tunit, &
&        ' Myr: a few-body group needs more than 2e7 steps in one go'
       warned=.TRUE.
       CYCLE
     END IF
!
! stars ic and jc of the group have collided
     m1=bmass(ic,k)
     m2=bmass(jc,k)
     mn=m1 + m2
     rr=bx(1:3,ic,k) - bx(1:3,jc,k)
     vr=bv(1:3,ic,k) - bv(1:3,jc,k)
     rsep=SQRT(SUM(rr*rr))
     xn=(m1*bx(1:3,ic,k) + m2*bx(1:3,jc,k))/mn
     vn=(m1*bv(1:3,ic,k) + m2*bv(1:3,jc,k))/mn
     IF (.NOT. mergers) THEN
       WRITE(6,'(a)') ''
       WRITE(6,'(a)') REPEAT('*',78)
       WRITE(6,'(a,f0.6,a)') '  STELLAR COLLISION at t = ', bepoch(k)*tunit, &
&        ' Myr (inside a few-body group)'
       WRITE(6,'(a,2i12,a,2f12.5,a)') '  star ids ', bident(ic,k), bident(jc,k), &
&        '   masses ', m1, m2, ' Msun'
       WRITE(6,'(a,3f14.6,a)') '  position                ', xc + xn, ' pc'
       WRITE(6,'(a)') '  The run has been stopped: mergers are switched off.'
       WRITE(6,'(a)') REPEAT('*',78)
       STOP 1
     END IF
     emerged=emerged + 0.5d0*m1*m2/mn*SUM(vr*vr) - m1*m2/SQRT(rsep*rsep + eps2)
     idn=bident(ic,k)
     IF (m2 > m1) idn=bident(jc,k)
     tn=mergetype(btype(ic,k),btype(jc,k))
     CALL mergelog(bepoch(k),bident(ic,k),m1,btype(ic,k),bident(jc,k),m2, &
&      btype(jc,k),idn,mn,tn,xc + xn,vc + vn,rsep)
! the merged star takes slot ic (ic < jc); the last star moves into jc
     bmass(ic,k)=mn
     bident(ic,k)=idn
     btype(ic,k)=tn
     bx(1:3,ic,k)=xn
     bv(1:3,ic,k)=vn
     IF (jc /= nm) THEN
       bmass(jc,k)=bmass(nm,k)
       bident(jc,k)=bident(nm,k)
       btype(jc,k)=btype(nm,k)
       bx(1:3,jc,k)=bx(1:3,nm,k)
       bv(1:3,jc,k)=bv(1:3,nm,k)
     END IF
     bnum(k)=nm - 1
     ntot=ntot - 1
     nmerged=nmerged + 1
   END DO
!
   END SUBROUTINE fewadvance
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
! two stars have come (or are about to come) within rcoll and mergers
! are switched off: stop the run with an alert
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
   WRITE(6,'(a)') '  The run has been stopped: mergers are switched off.'
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
   INTEGER FUNCTION mergetype(t1,t2)
!
! particle type of the star made by merging stars of types t1 and t2:
! the more compact of the two wins (black hole, then neutron star,
! then white dwarf, then star, then planet)
!
! >>> CHANGE HERE for different rules (e.g. two neutron stars making a
!     black hole); when stellar evolution is added it will decide this
!
   INTEGER, INTENT(IN) :: t1,t2
!
   IF (t1 == tbhole .OR. t2 == tbhole) THEN
     mergetype=tbhole
   ELSE IF (t1 == tnstar .OR. t2 == tnstar) THEN
     mergetype=tnstar
   ELSE IF (t1 == twdwarf .OR. t2 == twdwarf) THEN
     mergetype=twdwarf
   ELSE IF (t1 == tstar .OR. t2 == tstar) THEN
     mergetype=tstar
   ELSE
     mergetype=MAX(t1,t2)
   END IF
!
   END FUNCTION mergetype
! =====================================================================
   SUBROUTINE mergelog(tnow,id1,m1,t1,id2,m2,t2,idn,mn,tn,xn,vn,rsep)
!
! announce a merger on the screen and add a line to the merger file
! <runid>_mergers.dat in the data directory (made, with its header,
! at the first merger). Units Myr, Msun, pc, km/s; separation in Rsun.
!
   DOUBLE PRECISION, INTENT(IN) :: tnow,m1,m2,mn,xn(1:3),vn(1:3),rsep
   INTEGER, INTENT(IN) :: id1,t1,id2,t2,idn,tn
   CHARACTER(len=512) :: fname
   LOGICAL :: ex
   DOUBLE PRECISION :: rsun
!
   rsun=6.957d8/3.0856775814913673d16       ! solar radius in pc
   WRITE(6,'(a,f0.6,a,i0,a,i0,a,f0.4,a,f0.4,a,i0,a,f0.4,a,a)') &
&    '  ** MERGER t = ', tnow*tunit, ' Myr: stars ', id1, ' and ', id2, ' (', &
&    m1, ' and ', m2, ' Msun) are now star ', idn, ', ', mn, ' Msun, type ', &
&    TRIM(typename(tn))
!
   fname=TRIM(datadir)//TRIM(runid)//'_mergers.dat'
   INQUIRE(file=TRIM(fname),exist=ex)
   IF (ex) THEN
     OPEN(8,file=TRIM(fname),status='old',position='append')
   ELSE
     OPEN(8,file=TRIM(fname),status='new')
     WRITE(8,'(a)') '# format                        = NBODY-MERGERS 1'
     WRITE(8,'(a,a)') '# runid                         = ', TRIM(runid)
     WRITE(8,'(a)') '# units                         = mass Msun, length pc, '// &
&      'velocity km/s, time Myr'
     WRITE(8,'(a)') '# columns                       = time id1 m1 type1 id2 '// &
&      'm2 type2 id m type x y z vx vy vz separation'
     WRITE(8,'(a)') '# column_units                  = Myr - Msun - - Msun - - '// &
&      'Msun - pc pc pc km/s km/s km/s Rsun'
     WRITE(8,'(a)') '# column_notes                  = stars id1 and id2 merged '// &
&      'to make star id (the id of the more massive); x..vz are the new star'
     WRITE(8,'(a)') '# end_header'
   END IF
   WRITE(8,'(es24.16,2(i10,es24.16,1x,a12),i10,es24.16,1x,a12,6es24.16,es16.8)') &
&    tnow*tunit, id1, m1, typename(t1), id2, m2, typename(t2), idn, mn, &
&    typename(tn), xn, vn*vunit, rsep/rsun
   CLOSE(8)
!
   END SUBROUTINE mergelog
! =====================================================================
   SUBROUTINE mergestars(i,j,t8,tnow)
!
! merge single stars i and j into one star at their centre of mass
! (block time t8 = internal time tnow; their predicted positions rp, vp
! are used, which for the active star are its new position). Mass and
! momentum are conserved; the star keeps the id of the more massive.
! The orbital energy of the pair is added to emerged.
!
   INTEGER, INTENT(IN) :: i,j
   INTEGER(KIND=8), INTENT(IN) :: t8
   DOUBLE PRECISION, INTENT(IN) :: tnow
   INTEGER :: lo,hi,idn,tn,one(1:1)
   DOUBLE PRECISION :: mi,mj,mn,rr(1:3),vr(1:3),rsep,xc(1:3),vc(1:3)
!
   mi=m(i)
   mj=m(j)
   mn=mi + mj
   rr=rp(1:3,i) - rp(1:3,j)
   vr=vp(1:3,i) - vp(1:3,j)
   rsep=SQRT(SUM(rr*rr))
   xc=(mi*rp(1:3,i) + mj*rp(1:3,j))/mn
   vc=(mi*vp(1:3,i) + mj*vp(1:3,j))/mn
   emerged=emerged + 0.5d0*mi*mj/mn*SUM(vr*vr) - mi*mj/SQRT(rsep*rsep + eps2)
   idn=id(i)
   IF (mj > mi) idn=id(j)
   tn=mergetype(itype(i),itype(j))
   CALL mergelog(tnow,id(i),mi,itype(i),id(j),mj,itype(j),idn,mn,tn,xc,vc,rsep)
!
   lo=MIN(i,j)
   hi=MAX(i,j)
   r(1:3,lo)=xc
   v(1:3,lo)=vc
   rp(1:3,lo)=xc
   vp(1:3,lo)=vc
   m(lo)=mn
   id(lo)=idn
   itype(lo)=tn
   ib(lo)=0
   tk(lo)=t8
! remove slot hi by moving the last body into it
   IF (useac) THEN
     CALL acremoved(hi,lo)
     acok(lo)=.FALSE.
   END IF
   IF (hi /= n) THEN
     CALL moverow(n,hi)
     IF (useac) CALL acmoved(n,hi)
   END IF
   n=n - 1
   ntot=ntot - 1
   nmerged=nmerged + 1
   treevalid=.FALSE.
! forces and first timestep of the new star
   one(1)=lo
   CALL startbody(1,one,t8,0)
!
   END SUBROUTINE mergestars
! =====================================================================
   SUBROUTINE collcheck(t8,tnow)
!
! after each block (integer time t8 = internal time tnow): for each
! active body that came within 2 rcoll of another at its force
! evaluation, find that neighbour and test the pair: a collision if
! they are already within rcoll, or approaching with a two-body
! pericentre inside rcoll. (The flyby timestep limit means sampled
! separations near pericentre are within ~20% of the true pericentre,
! so 2 rcoll is a safe trigger.) The stars are then merged, or the run
! is stopped if mergers are off. If one of the two is a group it is
! opened instead, and its stars are tested on their own next steps.
! Serial; rare. At most one merger per block.
!
   INTEGER(KIND=8), INTENT(IN) :: t8
   DOUBLE PRECISION, INTENT(IN) :: tnow
   INTEGER :: ia,i,j,nact0,nm,list(1:kfmax)
   INTEGER, ALLOCATABLE :: alist(:)
   DOUBLE PRECISION :: rr(1:3),vr(1:3),rsep,q
!
! copy the active list (starting a new star overwrites ilist)
   nact0=nact
   ALLOCATE(alist(1:nact0))
   alist=ilist(1:nact0)
   DO ia=1,nact0
     i=alist(ia)
     IF (rmin2(i) >= 4.d0*rcolli*rcolli) CYCLE
     CALL nearest(i,j)
     IF (j == 0) CYCLE
     rr=rp(1:3,i) - rp(1:3,j)
     vr=vp(1:3,i) - vp(1:3,j)
     rsep=SQRT(SUM(rr*rr))
     q=pericentre(m(i) + m(j),rr,vr)
     IF (rsep < rcolli .OR. (q < rcolli .AND. SUM(rr*vr) < 0.d0)) THEN
       IF (.NOT. mergers) CALL collision(i,j,rsep,q,tnow)
       IF (ib(i) /= 0) THEN
         IF (tk(i) == t8) CALL opengroup(i,t8,tnow,nm,list)
       ELSE IF (ib(j) /= 0) THEN
         IF (tk(j) == t8) CALL opengroup(j,t8,tnow,nm,list)
       ELSE
         CALL mergestars(i,j,t8,tnow)
       END IF
       EXIT
     END IF
   END DO
   DEALLOCATE(alist)
!
   END SUBROUTINE collcheck
! =====================================================================
   SUBROUTINE gperturb(xc,vc,mb,gsize,nex,ex,gam,tp2)
!
! perturbation gamma on a (possible) group of mass mb and gsize `gsize`
! at xc from all bodies except the nex bodies in ex, and the squared
! shortest pairwise free-fall/flyby time of its c.o.m. (for its
! timestep). Uses rp, vp: every body's position predicted to the
! current block time.
!
   DOUBLE PRECISION, INTENT(IN) :: xc(1:3),vc(1:3),mb,gsize
   INTEGER, INTENT(IN) :: nex,ex(1:nex)
   DOUBLE PRECISION, INTENT(OUT) :: gam,tp2
   INTEGER :: j
   DOUBLE PRECISION :: sum3,d2,dv2,dx,dy,dz
!
   IF (usetree .AND. treevalid .AND. nex <= 2) THEN
     CALL treeperturb(xc,vc,mb,ex(1),ex(nex),sum3,tp2)
     gam=2.d0*gsize**3/mb*sum3
     RETURN
   END IF
   sum3=0.d0
   tp2=HUGE(1.d0)
!$OMP PARALLEL DO IF(n >= nompj) PRIVATE(j,d2,dv2,dx,dy,dz) &
!$OMP& REDUCTION(+:sum3) REDUCTION(MIN:tp2) SCHEDULE(STATIC)
   DO j=1,n
     IF (ANY(ex(1:nex) == j)) CYCLE
     dx=rp(1,j) - xc(1)
     dy=rp(2,j) - xc(2)
     dz=rp(3,j) - xc(3)
     d2=dx*dx + dy*dy + dz*dz
     dv2=(vp(1,j) - vc(1))**2 + (vp(2,j) - vc(2))**2 + (vp(3,j) - vc(3))**2
     sum3=sum3 + m(j)/(d2*SQRT(d2))
     tp2=MIN(tp2, d2*SQRT(d2)/(mb + m(j)), d2/(dv2 + TINY(1.d0)))
   END DO
!$OMP END PARALLEL DO
   gam=2.d0*gsize**3/mb*sum3
!
   END SUBROUTINE gperturb
! =====================================================================
   SUBROUTINE nearest(i,jn)
!
! index jn of the body nearest to body i (positions rp at the block time)
!
   INTEGER, INTENT(IN) :: i
   INTEGER, INTENT(OUT) :: jn
!
   IF (usetree .AND. treevalid) THEN
     CALL treenearest(i,jn)
   ELSE
     CALL nearestx(i,0,jn)
   END IF
!
   END SUBROUTINE nearest
! =====================================================================
   SUBROUTINE nearestx(i,jex,jn)
!
! index jn of the body nearest to body i, not counting body jex
! (jex = 0: any body). Direct search; 0 if there is no such body.
!
   INTEGER, INTENT(IN) :: i,jex
   INTEGER, INTENT(OUT) :: jn
   INTEGER :: j,jbest
   DOUBLE PRECISION :: d2,dbest,gbest
!
   jn=0
   gbest=HUGE(1.d0)
!$OMP PARALLEL IF(n >= nompj) PRIVATE(j,d2,dbest,jbest)
   dbest=HUGE(1.d0)
   jbest=0
!$OMP DO SCHEDULE(STATIC)
   DO j=1,n
     IF (j == i .OR. j == jex) CYCLE
     d2=(rp(1,j) - rp(1,i))**2 + (rp(2,j) - rp(2,i))**2 + (rp(3,j) - rp(3,i))**2
     IF (d2 < dbest) THEN
       dbest=d2
       jbest=j
     END IF
   END DO
!$OMP END DO
!$OMP CRITICAL (nearestcrit)
   IF (jbest > 0) THEN
     IF (dbest < gbest .OR. (dbest == gbest .AND. jbest < jn)) THEN
       gbest=dbest
       jn=jbest
     END IF
   END IF
!$OMP END CRITICAL (nearestcrit)
!$OMP END PARALLEL
!
   END SUBROUTINE nearestx
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
   IF (useac) CALL acmoverow(src,dst)
!
   END SUBROUTINE moverow
! =====================================================================
   SUBROUTINE startbody(nl,list,t8,levmin)
!
! start off the nl new bodies in list at integer time t8: forces from
! all other bodies (whose rp, vp are at t8), then the first timestep
! level of each, which is at least levmin and divides t8
!
   INTEGER, INTENT(IN) :: nl,list(1:nl),levmin
   INTEGER(KIND=8), INTENT(IN) :: t8
   INTEGER :: ia,l,knew
   DOUBLE PRECISION :: dtw,amod,jmod
!
   tcur8=t8
   aclast=.TRUE.
   IF (usetree) CALL treebuild
   nact=nl
   ilist(1:nl)=list(1:nl)
   CALL forces
   DO ia=1,nl
     l=list(ia)
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
     lev(l)=commlev(MAX(knew,levmin),t8)
     dtk(l)=ISHFT(1_8,kmax - lev(l))
     IF (useac) CALL acstart(l)
   END DO
!
   END SUBROUTINE startbody
! =====================================================================
   SUBROUTINE closegroup(nm,list,mode,t8,tnow,rel,rel2)
!
! replace the nm single stars in list (all at integer time t8 =
! internal time tnow) by one c.o.m. body, as a group of kind mode
! (1 binary, 2 few-body, 3 stable triple - then rel, rel2 are its
! inner and outer orbits and list is in the order inner, inner, outer)
!
   INTEGER, INTENT(IN) :: nm,list(1:nm),mode
   INTEGER(KIND=8), INTENT(IN) :: t8
   DOUBLE PRECISION, INTENT(IN) :: tnow,rel(1:6),rel2(1:6)
   INTEGER :: k,l,i,lo,h,knew,rest(1:kfmax),nrest,one(1:1)
   DOUBLE PRECISION :: mb,xc(1:3),vc(1:3),gam,tp2,dtw,amod,jmod,ebind,gsize
   DOUBLE PRECISION :: fc(1:3,1:4)
!
   k=nbin + 1
   nbin=k
   bnum(k)=nm
   bmode(k)=mode
   bepoch(k)=tnow
   mb=0.d0
   xc=0.d0
   vc=0.d0
   fc=0.d0
   DO l=1,nm
     i=list(l)
     bmass(l,k)=m(i)
     bident(l,k)=id(i)
     btype(l,k)=itype(i)
     mb=mb + m(i)
     xc=xc + m(i)*r(1:3,i)
     vc=vc + m(i)*v(1:3,i)
! c.o.m. derivatives: the mass-weighted mean, in which the mutual
! forces of the stars cancel exactly
     fc(1:3,1)=fc(1:3,1) + m(i)*f0(1:3,i)
     fc(1:3,2)=fc(1:3,2) + m(i)*f1(1:3,i)
     fc(1:3,3)=fc(1:3,3) + m(i)*f2e(1:3,i)
     fc(1:3,4)=fc(1:3,4) + m(i)*f3e(1:3,i)
   END DO
   xc=xc/mb
   vc=vc/mb
   fc=fc/mb
   gsize=0.d0
   SELECT CASE (mode)
   CASE (1)
     brel(1:3,k)=r(1:3,list(1)) - r(1:3,list(2))
     brel(4:6,k)=v(1:3,list(1)) - v(1:3,list(2))
     CALL pairorbit(mb,brel(1:3,k),brel(4:6,k),ebind,gsize)
     nclosed=nclosed + 1
   CASE (3)
     brel(1:6,k)=rel
     brel2(1:6,k)=rel2
   CASE DEFAULT
     nfewformed=nfewformed + 1
   END SELECT
   IF (mode /= 1) THEN
     DO l=1,nm
       i=list(l)
       IF (mode == 2) THEN
         bx(1:3,l,k)=r(1:3,i) - xc
         bv(1:3,l,k)=v(1:3,i) - vc
       END IF
       gsize=MAX(gsize,2.d0*SQRT(SUM((r(1:3,i) - xc)**2)))
     END DO
   END IF
!
! timestep of the c.o.m. (before the stars are removed)
   CALL gperturb(xc,vc,mb,gsize,nm,list,gam,tp2)
!
   lo=MINVAL(list(1:nm))
   f0(1:3,lo)=fc(1:3,1)
   f1(1:3,lo)=fc(1:3,2)
   f2e(1:3,lo)=fc(1:3,3)
   f3e(1:3,lo)=fc(1:3,4)
   r(1:3,lo)=xc
   v(1:3,lo)=vc
   rp(1:3,lo)=xc
   vp(1:3,lo)=vc
   m(lo)=mb
   id(lo)=0
   itype(lo)=0
   ib(lo)=k
   tk(lo)=t8
   tpair2(lo)=tp2
   dtw=etap*SQRT(tp2)
   amod=SQRT(SUM(f0(1:3,lo)**2))
   jmod=SQRT(SUM(f1(1:3,lo)**2))
   IF (amod > 0.d0 .AND. jmod > 0.d0) dtw=MIN(dtw,etas*amod/jmod)
   CALL getlevel(dtw,knew)
   lev(lo)=commlev(knew,t8)
   dtk(lo)=ISHFT(1_8,kmax - lev(lo))
!
! remove the other slots, highest first, moving the last body into each
   nrest=0
   DO l=1,nm
     IF (list(l) /= lo) THEN
       nrest=nrest + 1
       rest(nrest)=list(l)
     END IF
   END DO
   DO WHILE (nrest > 0)
     l=MAXLOC(rest(1:nrest),DIM=1)
     h=rest(l)
     rest(l)=rest(nrest)
     nrest=nrest - 1
     IF (useac) CALL acremoved(h,lo)
     IF (h /= n) THEN
       CALL moverow(n,h)
       IF (useac) CALL acmoved(n,h)
     END IF
     n=n - 1
   END DO
   treevalid=.FALSE.
!
! with the neighbour scheme the c.o.m. body needs its own list and
! force split, so it is started with a full force sum
   IF (useac) THEN
     acok(lo)=.FALSE.
     one(1)=lo
     CALL startbody(1,one,t8,0)
   END IF
!
   END SUBROUTINE closegroup
! =====================================================================
   SUBROUTINE opengroup(i,t8,tnow,nm,list)
!
! split the c.o.m. body i (at integer time t8 = internal time tnow)
! back into its nm stars, in slot i and slots n+1 onwards, and start
! them off with forces from all other bodies (whose rp, vp are at t8).
! list returns the slots of the stars.
!
   INTEGER, INTENT(IN) :: i
   INTEGER(KIND=8), INTENT(IN) :: t8
   DOUBLE PRECISION, INTENT(IN) :: tnow
   INTEGER, INTENT(OUT) :: nm,list(1:kfmax)
   INTEGER :: k,l,j,levc,mode
   DOUBLE PRECISION :: xm(1:3,1:kfmax),um(1:3,1:kfmax)
!
   k=ib(i)
   levc=lev(i)
   mode=bmode(k)
   CALL grpcomp(k,r(1:3,i),v(1:3,i),tnow,nm,xm,um)
   list=0
   DO l=1,nm
     IF (l == 1) THEN
       j=i
     ELSE
       j=n + 1
       n=j
       IF (useac) THEN
         nnb(j)=0
         acok(j)=.FALSE.
         CALL acadded(j,i)
       END IF
     END IF
     r(1:3,j)=xm(1:3,l)
     v(1:3,j)=um(1:3,l)
     rp(1:3,j)=xm(1:3,l)
     vp(1:3,j)=um(1:3,l)
     m(j)=bmass(l,k)
     id(j)=bident(l,k)
     itype(j)=btype(l,k)
     ib(j)=0
     phi(j)=phi(i)
     rmin2(j)=HUGE(1.d0)
     list(l)=j
   END DO
   IF (useac) acok(i)=.FALSE.
!
! remove group k from the table, moving the last group into its place
   IF (k /= nbin) THEN
     bnum(k)=bnum(nbin)
     bmode(k)=bmode(nbin)
     bmass(1:kfmax,k)=bmass(1:kfmax,nbin)
     bident(1:kfmax,k)=bident(1:kfmax,nbin)
     btype(1:kfmax,k)=btype(1:kfmax,nbin)
     brel(1:6,k)=brel(1:6,nbin)
     brel2(1:6,k)=brel2(1:6,nbin)
     bepoch(k)=bepoch(nbin)
     IF (fewbody) THEN
       bx(1:3,1:kfmax,k)=bx(1:3,1:kfmax,nbin)
       bv(1:3,1:kfmax,k)=bv(1:3,1:kfmax,nbin)
     END IF
     DO l=1,n
       IF (ib(l) == nbin) THEN
         ib(l)=k
         EXIT
       END IF
     END DO
   END IF
   nbin=nbin - 1
   IF (mode == 1) THEN
     nopened=nopened + 1
   ELSE
     nfewopened=nfewopened + 1
   END IF
   treevalid=.FALSE.
!
! forces on the stars and their first timesteps (at least the c.o.m.
! level, so the step divides t8)
   CALL startbody(nm,list(1:nm),t8,levc)
!
   END SUBROUTINE opengroup
! =====================================================================
   INTEGER FUNCTION pairstate(i,j)
!
! what single stars i and j are to each other:
!   0  nothing to do (not both single, unbound, or perturbed)
!   1  a bound pair with gamma < gclose: can be closed as a binary
!   2  a bound pair whose pericentre is inside rcoll: a collision
!
   INTEGER, INTENT(IN) :: i,j
   INTEGER :: two(1:2)
   DOUBLE PRECISION :: mb,rr(1:3),vr(1:3),xc(1:3),vc(1:3),ebind,rapo,gam,tp2
!
   pairstate=0
   IF (ib(i) /= 0 .OR. ib(j) /= 0) RETURN
   mb=m(i) + m(j)
   rr=r(1:3,i) - r(1:3,j)
   vr=v(1:3,i) - v(1:3,j)
   CALL pairorbit(mb,rr,vr,ebind,rapo)
   IF (ebind >= 0.d0) RETURN
! the orbit of a closed binary is fixed, so check its pericentre now
   IF (rcoll > 0.d0) THEN
     IF (pericentre(mb,rr,vr) < rcolli) THEN
       pairstate=2
       RETURN
     END IF
   END IF
   xc=(m(i)*r(1:3,i) + m(j)*r(1:3,j))/mb
   vc=(m(i)*v(1:3,i) + m(j)*v(1:3,j))/mb
   two(1)=i
   two(2)=j
   CALL gperturb(xc,vc,mb,rapo,2,two,gam,tp2)
   IF (gam < gclose) pairstate=1
!
   END FUNCTION pairstate
! =====================================================================
   SUBROUTINE paircollide(i,j,t8,tnow)
!
! the bound pair i, j has a pericentre inside rcoll: merge the stars,
! or stop the run if mergers are off
!
   INTEGER, INTENT(IN) :: i,j
   INTEGER(KIND=8), INTENT(IN) :: t8
   DOUBLE PRECISION, INTENT(IN) :: tnow
   DOUBLE PRECISION :: rr(1:3),vr(1:3)
!
   rr=r(1:3,i) - r(1:3,j)
   vr=v(1:3,i) - v(1:3,j)
   IF (.NOT. mergers) CALL collision(i,j,SQRT(SUM(rr*rr)), &
&    pericentre(m(i) + m(j),rr,vr),tnow)
   CALL mergestars(i,j,t8,tnow)
!
   END SUBROUTINE paircollide
! =====================================================================
   SUBROUTINE grpgamma(i,tnow,gam)
!
! perturbation gamma of the group whose c.o.m. body is i (at internal
! time tnow) from all other bodies. A few-body group is brought up to
! date first.
!
   INTEGER, INTENT(IN) :: i
   DOUBLE PRECISION, INTENT(IN) :: tnow
   DOUBLE PRECISION, INTENT(OUT) :: gam
   INTEGER :: k,nm,l,one(1:1)
   DOUBLE PRECISION :: xm(1:3,1:kfmax),um(1:3,1:kfmax),ebind,gsize,tp2
!
   k=ib(i)
   CALL grpcomp(k,r(1:3,i),v(1:3,i),tnow,nm,xm,um)
   IF (bmode(k) == 1) THEN
     CALL pairorbit(m(i),xm(1:3,1) - xm(1:3,2),um(1:3,1) - um(1:3,2),ebind,gsize)
   ELSE
     gsize=0.d0
     DO l=1,nm
       gsize=MAX(gsize,2.d0*SQRT(SUM((xm(1:3,l) - r(1:3,i))**2)))
     END DO
   END IF
   one(1)=i
   CALL gperturb(r(1:3,i),v(1:3,i),m(i),gsize,1,one,gam,tp2)
!
   END SUBROUTINE grpgamma
! =====================================================================
   SUBROUTINE tryfreeze(k)
!
! few-body group k has three stars: if it is a stable hierarchy (and
! the inner binary does not touch), freeze it as two nested Kepler
! orbits (mode 3). Its stars are put in the order inner, inner, outer.
! The group must be up to date (bx, bv at time bepoch).
!
   INTEGER, INTENT(IN) :: k
   LOGICAL :: stable
   INTEGER :: i1,i2,i3,idx(1:3),ti(1:3),tt(1:3)
   DOUBLE PRECISION :: rel(1:6),rel2(1:6),qin,qout,tm(1:3)
!
   CALL hierarchy(bmass(1:3,k),bx(1:3,1:3,k),bv(1:3,1:3,k),stable,i1,i2,i3, &
&    rel,rel2,qin,qout)
   IF (.NOT. stable) RETURN
   IF (rcoll > 0.d0 .AND. qin < rcolli) RETURN
   idx=(/i1,i2,i3/)
   tm=bmass(idx,k)
   ti=bident(idx,k)
   tt=btype(idx,k)
   bmass(1:3,k)=tm
   bident(1:3,k)=ti
   btype(1:3,k)=tt
   brel(1:6,k)=rel
   brel2(1:6,k)=rel2
   bmode(k)=3
   nfrozen=nfrozen + 1
!
   END SUBROUTINE tryfreeze
! =====================================================================
   SUBROUTINE tryfew(nb,bod,t8,tnow,formed)
!
! try to make a few-body group out of the nb bodies in bod (single
! stars or groups, all at integer time t8 = internal time tnow). It is
! made if it has 3 to kfmax stars, its perturbation by the rest of
! the cluster is below gfew, and it would pay (below). Any groups among
! the bodies are opened first, then all the stars are closed as one
! group. A group of three that is a stable hierarchy is frozen as two
! Kepler orbits straight away.
!
! WHEN A GROUP PAYS. Taking nm stars out saves work only if the main
! integration is big enough. There, each step of a close pair costs
! about 2.5 n pair interactions (two force sums over the n bodies and
! the prediction of all of them). In the few-body integrator a step
! costs about 1.3 nm (nm-1) (two evaluations of every pair), and it
! takes 2 to 4 times as many steps per orbit as the main integrator
! (measured on binaries with eccentricity 0.9 to 0, at the standard
! eta and etafew; it is also about 1000 times more accurate per orbit).
! The costs are equal when
!      n = (1.1 to 2.3) nm (nm-1)
! so a group is only made if n >= nfewfac nm (nm-1), nfewfac = 2:
!      stars in the group        3     4     5     6
!      bodies needed            12    24    40    60
! The break-even depends on the size of the group, not on a fixed
! number of stars, and it is far below that of the neighbour scheme
! (150), because the saving per step is n/nm^2, not n/neighbours.
! Accuracy does not depend on n at all: it is set by gfew and etafew.
! A stable triple is always worth it, whatever n is, because frozen as
! Kepler orbits it costs nothing: so three stars that are a stable
! hierarchy are taken out even if n is below the number above.
!
! >>> CHANGE HERE (nfewfac, above CONTAINS) once the break-even has
!     been measured with the compiled code: benchmark/bench_kira.py
!     --skip-kira --modes nofew direct
!
   INTEGER, INTENT(IN) :: nb,bod(1:nb)
   INTEGER(KIND=8), INTENT(IN) :: t8
   DOUBLE PRECISION, INTENT(IN) :: tnow
   LOGICAL, INTENT(OUT) :: formed
   INTEGER :: b,i,k,nm,nk,l,nl,list(1:kfmax),lk(1:kfmax),i1,i2,i3
   DOUBLE PRECISION :: xg(1:3,1:kfmax),ug(1:3,1:kfmax),xm(1:3,1:kfmax)
   DOUBLE PRECISION :: um(1:3,1:kfmax),mm(1:kfmax)
   DOUBLE PRECISION :: mt,xc(1:3),vc(1:3),gsize,gam,tp2,dum(1:6)
   DOUBLE PRECISION :: rel(1:6),rel2(1:6),qin,qout
   LOGICAL :: stable
!
   formed=.FALSE.
! the stars the group would have
   nm=0
   mt=0.d0
   xc=0.d0
   vc=0.d0
   DO b=1,nb
     i=bod(b)
     k=ib(i)
     nk=1
     IF (k /= 0) nk=bnum(k)
     IF (nm + nk > kfmax) RETURN
     IF (k == 0) THEN
       xm(1:3,nm + 1)=r(1:3,i)
       um(1:3,nm + 1)=v(1:3,i)
       mm(nm + 1)=m(i)
     ELSE
       CALL grpcomp(k,r(1:3,i),v(1:3,i),tnow,nk,xg,ug)
       IF (nm + nk > kfmax) RETURN
       xm(1:3,nm + 1:nm + nk)=xg(1:3,1:nk)
       um(1:3,nm + 1:nm + nk)=ug(1:3,1:nk)
       mm(nm + 1:nm + nk)=bmass(1:nk,k)
     END IF
     nm=nm + nk
     mt=mt + m(i)
     xc=xc + m(i)*r(1:3,i)
     vc=vc + m(i)*v(1:3,i)
   END DO
   IF (nm < 3) RETURN
   xc=xc/mt
   vc=vc/mt
   gsize=0.d0
   DO l=1,nm
     gsize=MAX(gsize,2.d0*SQRT(SUM((xm(1:3,l) - xc)**2)))
   END DO
   CALL gperturb(xc,vc,mt,gsize,nb,bod,gam,tp2)
   IF (gam >= gfew) RETURN
!
! would it pay? Not if the main integration is too small - unless it
! is a stable triple, which will cost nothing once it is frozen
   IF (n < nfewfac*nm*(nm - 1)) THEN
     IF (nm /= 3) RETURN
     CALL hierarchy(mm(1:3),xm(1:3,1:3),um(1:3,1:3),stable,i1,i2,i3,rel,rel2,qin,qout)
     IF (.NOT. stable) RETURN
     IF (rcoll > 0.d0 .AND. qin < rcolli) RETURN
   END IF
!
! make it: open the groups, then close all the stars together
   nl=0
   DO b=1,nb
     i=bod(b)
     IF (ib(i) /= 0) THEN
       CALL opengroup(i,t8,tnow,nk,lk)
       list(nl + 1:nl + nk)=lk(1:nk)
       nl=nl + nk
     ELSE
       nl=nl + 1
       list(nl)=i
     END IF
   END DO
   IF (nl < 3) RETURN
   dum=0.d0
   CALL closegroup(nl,list(1:nl),2,t8,tnow,dum,dum)
   formed=.TRUE.
! a stable triple is frozen at once (the new group is the last one)
   IF (nl == 3) CALL tryfreeze(nbin)
!
   END SUBROUTINE tryfew
! =====================================================================
   SUBROUTINE bincheck(t8,tnow)
!
! called after each block (all active bodies at integer time t8 =
! internal time tnow). Among the bodies that have just stepped:
!  - groups: freeze a few-body group that has become a stable triple;
!    open one group that has become perturbed (or take its perturber
!    into a few-body group, if that group would be isolated)
!  - single stars with short steps: close one unperturbed pair, merge a
!    pair that is about to collide, or make a few-body group of a
!    perturbed pair and the next nearest body
!
   INTEGER(KIND=8), INTENT(IN) :: t8
   DOUBLE PRECISION, INTENT(IN) :: tnow
   INTEGER :: ia,i,j,l,k,nact0,nm,list(1:kfmax),two(1:2),three(1:3)
   INTEGER, ALLOCATABLE :: alist(:)
   DOUBLE PRECISION :: gam,gmax,dum(1:6)
   LOGICAL :: done
!
! copy the active list (starting new bodies overwrites ilist)
   nact0=nact
   ALLOCATE(alist(1:nact0))
   alist=ilist(1:nact0)
   dum=0.d0
!
! groups
   DO ia=1,nact0
     i=alist(ia)
     IF (i > n) CYCLE
     IF (ib(i) == 0 .OR. tk(i) /= t8) CYCLE
     CALL grpgamma(i,tnow,gam)
     k=ib(i)
     gmax=gopen
     IF (bmode(k) /= 1) gmax=10.d0*gfew
     IF (bmode(k) == 2) THEN
       IF (bnum(k) < 3) THEN
! mergers have left fewer than three stars
         CALL opengroup(i,t8,tnow,nm,list)
         EXIT
       END IF
       IF (bnum(k) == 3 .AND. gam <= gmax) CALL tryfreeze(k)
     END IF
     IF (gam > gmax) THEN
       done=.FALSE.
       IF (fewbody) THEN
         CALL nearest(i,j)
         IF (j > 0) THEN
           IF (tk(j) == t8) THEN
             two(1)=i
             two(2)=j
             CALL tryfew(2,two,t8,tnow,done)
           END IF
         END IF
       END IF
       IF (.NOT. done) CALL opengroup(i,t8,tnow,nm,list)
       EXIT
     END IF
   END DO
!
! single stars - only stars with short steps are candidates
   DO ia=1,nact0
     i=alist(ia)
     IF (i > n) CYCLE
     IF (ib(i) /= 0 .OR. tk(i) /= t8) CYCLE
     IF (DBLE(dtk(i))*tick >= dtclosei) CYCLE
     CALL nearest(i,j)
     IF (j == 0) CYCLE
     IF (tk(j) /= t8) CYCLE
     two(1)=i
     two(2)=j
     done=.FALSE.
     IF (ib(j) /= 0) THEN
! the nearest body is a group
       IF (fewbody) CALL tryfew(2,two,t8,tnow,done)
     ELSE
       SELECT CASE (pairstate(i,j))
       CASE (1)
         CALL closegroup(2,two,1,t8,tnow,dum,dum)
         done=.TRUE.
       CASE (2)
         CALL paircollide(i,j,t8,tnow)
         done=.TRUE.
       CASE DEFAULT
         IF (fewbody) THEN
           CALL nearestx(i,j,l)
           IF (l > 0) THEN
             IF (tk(l) == t8) THEN
               three(1)=i
               three(2)=j
               three(3)=l
               CALL tryfew(3,three,t8,tnow,done)
             END IF
           END IF
         END IF
       END SELECT
     END IF
     IF (done) EXIT
   END DO
   DEALLOCATE(alist)
!
   END SUBROUTINE bincheck
! =====================================================================
   SUBROUTINE binstart(tnow)
!
! at the start (everyone is synchronised at integer time 0 = internal
! time tnow): close all unperturbed binaries in the initial conditions,
! merge pairs that would collide, and make few-body groups of
! perturbed pairs with their nearest body where those are isolated
! (stable triples among them are frozen at their first step). After a
! change, slot i may hold a different body, so i only advances when
! nothing was done.
!
   DOUBLE PRECISION, INTENT(IN) :: tnow
   INTEGER :: i,j,l,two(1:2),three(1:3)
   DOUBLE PRECISION :: dum(1:6)
   LOGICAL :: done
!
   dum=0.d0
   i=1
   DO WHILE (i <= n)
     done=.FALSE.
     IF (ib(i) == 0) THEN
       CALL nearest(i,j)
       IF (j > 0) THEN
         IF (ib(j) == 0) THEN
           two(1)=i
           two(2)=j
           SELECT CASE (pairstate(i,j))
           CASE (1)
             CALL closegroup(2,two,1,0_8,tnow,dum,dum)
             done=.TRUE.
           CASE (2)
             CALL paircollide(i,j,0_8,tnow)
             done=.TRUE.
           CASE DEFAULT
! only bound pairs: a perturbed binary and its perturber
             IF (fewbody) THEN
               IF (0.5d0*SUM((v(1:3,i) - v(1:3,j))**2) < (m(i) + m(j))/ &
&                  SQRT(SUM((r(1:3,i) - r(1:3,j))**2))) THEN
                 CALL nearestx(i,j,l)
                 IF (l > 0) THEN
                   three(1)=i
                   three(2)=j
                   three(3)=l
                   CALL tryfew(3,three,0_8,tnow,done)
                 END IF
               END IF
             END IF
           END SELECT
         END IF
       END IF
     END IF
     IF (.NOT. done) i=i + 1
   END DO
!
   END SUBROUTINE binstart
!
   END MODULE nbbin
