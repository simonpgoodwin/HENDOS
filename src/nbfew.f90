! ==============================================================
   MODULE nbfew
!
! HENDOS - HErmite N-body Dynamics Of Systems
! FEW-BODY TOOLS
! A small integrator for an isolated group of a few stars, and the
! stability test for a hierarchical triple. These routines know nothing
! about the rest of the code: everything comes in and goes out through
! their arguments, in internal units (G = 1). They are used by the
! group routines in module nbbin.
!
!   fewint      integrates the group for a given time
!   fewacc      accelerations and jerks inside the group
!   elements    semi-major axis, eccentricity and angular momentum of
!               a two-body orbit
!   hierarchy   is a triple a stable hierarchy (inner binary + distant
!               third star)?
!
   IMPLICIT NONE
!
   CONTAINS
! =====================================================================
   SUBROUTINE fewacc(nm,mm,x,v,e2,a,jk,tmin2,r2min,ic,jc)
!
! accelerations a and jerks jk of the nm stars of a group from each
! other, the square of the shortest pairwise free-fall or flyby time
! tmin2, and the closest pair ic, jc with its squared separation r2min
!
   INTEGER, INTENT(IN) :: nm
   DOUBLE PRECISION, INTENT(IN) :: mm(1:nm),x(1:3,1:nm),v(1:3,1:nm),e2
   DOUBLE PRECISION, INTENT(OUT) :: a(1:3,1:nm),jk(1:3,1:nm),tmin2,r2min
   INTEGER, INTENT(OUT) :: ic,jc
   INTEGER :: i,j
   DOUBLE PRECISION :: d(1:3),dv(1:3),d2,rad2,rinv,rinv3,rv,w(1:3)
!
   a=0.d0
   jk=0.d0
   tmin2=HUGE(1.d0)
   r2min=HUGE(1.d0)
   ic=0
   jc=0
   DO i=1,nm - 1
     DO j=i + 1,nm
       d=x(1:3,j) - x(1:3,i)
       dv=v(1:3,j) - v(1:3,i)
       d2=SUM(d*d)
       IF (d2 < r2min) THEN
         r2min=d2
         ic=i
         jc=j
       END IF
       rad2=d2 + e2
       rinv=1.d0/SQRT(rad2)
       rinv3=rinv/rad2
       rv=3.d0*SUM(d*dv)/rad2
       w=(dv - rv*d)*rinv3
       a(1:3,i)=a(1:3,i) + mm(j)*d*rinv3
       a(1:3,j)=a(1:3,j) - mm(i)*d*rinv3
       jk(1:3,i)=jk(1:3,i) + mm(j)*w
       jk(1:3,j)=jk(1:3,j) - mm(i)*w
       tmin2=MIN(tmin2, rad2*rad2*rinv/(mm(i) + mm(j)), rad2/(SUM(dv*dv) + TINY(1.d0)))
     END DO
   END DO
!
   END SUBROUTINE fewacc
! =====================================================================
   SUBROUTINE fewint(nm,mm,x,v,tint,etaf,e2,rc2,tdone,ic,jc,nst)
!
! integrate the isolated group of nm stars (masses mm, positions x and
! velocities v relative to the group's centre of mass) for a time tint.
! 4th-order Hermite scheme (predict, evaluate, correct, evaluate) with
! one timestep shared by all the stars,
!    dt = etaf * (shortest pairwise free-fall or flyby time)
! so the step follows the closest pair through its pericentre.
!
! It stops early
!   - if two stars come within sqrt(rc2) of each other (rc2 > 0):
!     ic, jc are then the pair and tdone the time reached
!   - after nstmax steps: ic = -1 and tdone is the time reached
! otherwise ic = 0 and tdone = tint. nst returns the number of steps.
!
   INTEGER, INTENT(IN) :: nm
   DOUBLE PRECISION, INTENT(IN) :: mm(1:nm),tint,etaf,e2,rc2
   DOUBLE PRECISION, INTENT(INOUT) :: x(1:3,1:nm),v(1:3,1:nm)
   DOUBLE PRECISION, INTENT(OUT) :: tdone
   INTEGER, INTENT(OUT) :: ic,jc
   INTEGER(KIND=8), INTENT(OUT) :: nst
   INTEGER(KIND=8), PARAMETER :: nstmax=20000000_8
   INTEGER :: i,k1,k2
   LOGICAL :: last
   DOUBLE PRECISION :: a0(1:3,1:nm),j0(1:3,1:nm),a1(1:3,1:nm),j1(1:3,1:nm)
   DOUBLE PRECISION :: xp(1:3,1:nm),vp(1:3,1:nm),v1(1:3,1:nm)
   DOUBLE PRECISION :: t,dt,tmin2,r2min,cx(1:3),cv(1:3),mt
!
   t=0.d0
   nst=0
   tdone=0.d0
   CALL fewacc(nm,mm,x,v,e2,a0,j0,tmin2,r2min,k1,k2)
   DO
     IF (rc2 > 0.d0 .AND. r2min < rc2) THEN
       ic=k1
       jc=k2
       tdone=t
       RETURN
     END IF
     IF (nst >= nstmax) THEN
       ic=-1
       jc=0
       tdone=t
       RETURN
     END IF
     IF (t >= tint) EXIT
     dt=etaf*SQRT(tmin2)
     last=.FALSE.
     IF (t + dt >= tint) THEN
       dt=tint - t
       last=.TRUE.
     END IF
! predict
     xp=x + v*dt + a0*dt**2/2.d0 + j0*dt**3/6.d0
     vp=v + a0*dt + j0*dt**2/2.d0
! evaluate and correct
     CALL fewacc(nm,mm,xp,vp,e2,a1,j1,tmin2,r2min,k1,k2)
     v1=v + (a0 + a1)*dt/2.d0 + (j0 - j1)*dt**2/12.d0
     x=x + (v + v1)*dt/2.d0 + (a0 - a1)*dt**2/12.d0
     v=v1
! evaluate at the corrected positions: the start of the next step
     CALL fewacc(nm,mm,x,v,e2,a0,j0,tmin2,r2min,k1,k2)
     nst=nst + 1
     IF (last) THEN
       t=tint
     ELSE
       t=t + dt
     END IF
   END DO
   ic=0
   jc=0
   tdone=tint
!
! keep the centre of mass exactly at the origin
   mt=SUM(mm)
   cx=0.d0
   cv=0.d0
   DO i=1,nm
     cx=cx + mm(i)*x(1:3,i)
     cv=cv + mm(i)*v(1:3,i)
   END DO
   DO i=1,nm
     x(1:3,i)=x(1:3,i) - cx/mt
     v(1:3,i)=v(1:3,i) - cv/mt
   END DO
!
   END SUBROUTINE fewint
! =====================================================================
   SUBROUTINE elements(mu,rr,vr,a,ecc,h)
!
! semi-major axis a (negative if unbound), eccentricity ecc and
! specific angular momentum vector h of the relative orbit rr, vr of
! two bodies with G*(total mass) = mu
!
   DOUBLE PRECISION, INTENT(IN) :: mu,rr(1:3),vr(1:3)
   DOUBLE PRECISION, INTENT(OUT) :: a,ecc,h(1:3)
   DOUBLE PRECISION :: e
!
   e=0.5d0*SUM(vr*vr) - mu/SQRT(SUM(rr*rr))
   a=-0.5d0*mu/e
   h(1)=rr(2)*vr(3) - rr(3)*vr(2)
   h(2)=rr(3)*vr(1) - rr(1)*vr(3)
   h(3)=rr(1)*vr(2) - rr(2)*vr(1)
   ecc=SQRT(MAX(0.d0,1.d0 + 2.d0*e*SUM(h*h)/mu**2))
!
   END SUBROUTINE elements
! =====================================================================
   SUBROUTINE hierarchy(mm,x,v,stable,i1,i2,i3,rel,rel2,qin,qout)
!
! is the triple (masses mm, positions x, velocities v) a stable
! hierarchy? The inner binary is the bound pair with the smallest
! semi-major axis; the third star must be bound to it, and the
! pericentre of the outer orbit must satisfy the stability criterion of
! Mardling & Aarseth (2001, MNRAS, 321, 398),
!   q_out/a_in > 2.8 [ (1 + m3/(m1+m2)) (1 + e_out)/sqrt(1 - e_out) ]^(2/5)
!                    * (1 - 0.3 i/pi)
! with i the angle between the two orbital planes. Returns
!   i1, i2, i3   the inner pair and the outer star
!   rel          inner orbit: r, v of star i1 relative to star i2
!   rel2         outer orbit: r, v of star i3 relative to the inner c.o.m.
!   qin, qout    pericentre distances of the inner and outer orbits
!
! >>> CHANGE HERE for a different stability criterion
!
   DOUBLE PRECISION, INTENT(IN) :: mm(1:3),x(1:3,1:3),v(1:3,1:3)
   LOGICAL, INTENT(OUT) :: stable
   INTEGER, INTENT(OUT) :: i1,i2,i3
   DOUBLE PRECISION, INTENT(OUT) :: rel(1:6),rel2(1:6),qin,qout
   INTEGER :: p,ia,ib,ic
   INTEGER, PARAMETER :: pa(3)=(/1,1,2/),pb(3)=(/2,3,3/),pc(3)=(/3,2,1/)
   DOUBLE PRECISION :: a,ecc,h(1:3),abest,ain,ein,hin(1:3),aout,eout,hout(1:3)
   DOUBLE PRECISION :: mb,mt,xc(1:3),vc(1:3),cosi,inc,crit,pi
!
   pi=3.141592653589793d0
   stable=.FALSE.
   i1=0
   i2=0
   i3=0
   rel=0.d0
   rel2=0.d0
   qin=0.d0
   qout=0.d0
   abest=HUGE(1.d0)
   DO p=1,3
     ia=pa(p)
     ib=pb(p)
     CALL elements(mm(ia) + mm(ib),x(1:3,ia) - x(1:3,ib),v(1:3,ia) - v(1:3,ib),a,ecc,h)
     IF (a > 0.d0 .AND. a < abest) THEN
       abest=a
       i1=ia
       i2=ib
       i3=pc(p)
     END IF
   END DO
   IF (i1 == 0) RETURN
!
   mb=mm(i1) + mm(i2)
   mt=mb + mm(i3)
   rel(1:3)=x(1:3,i1) - x(1:3,i2)
   rel(4:6)=v(1:3,i1) - v(1:3,i2)
   CALL elements(mb,rel(1:3),rel(4:6),ain,ein,hin)
   xc=(mm(i1)*x(1:3,i1) + mm(i2)*x(1:3,i2))/mb
   vc=(mm(i1)*v(1:3,i1) + mm(i2)*v(1:3,i2))/mb
   rel2(1:3)=x(1:3,i3) - xc
   rel2(4:6)=v(1:3,i3) - vc
   CALL elements(mt,rel2(1:3),rel2(4:6),aout,eout,hout)
   IF (aout <= 0.d0 .OR. eout >= 1.d0) RETURN
   qin=ain*(1.d0 - ein)
   qout=aout*(1.d0 - eout)
   cosi=SUM(hin*hout)/SQRT(SUM(hin*hin)*SUM(hout*hout) + TINY(1.d0))
   inc=ACOS(MAX(-1.d0,MIN(1.d0,cosi)))
   crit=2.8d0*((1.d0 + mm(i3)/mb)*(1.d0 + eout)/SQRT(1.d0 - eout))**0.4d0 &
&       *(1.d0 - 0.3d0*inc/pi)
   stable=(qout > crit*ain)
!
   END SUBROUTINE hierarchy
!
   END MODULE nbfew
