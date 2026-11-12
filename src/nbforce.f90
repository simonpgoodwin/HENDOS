! ==============================================================
   MODULE nbforce
!
! HENDOS - HErmite N-body Dynamics Of Systems
! GRAVITY
! accelerations and jerks on the active particles (direct summation or
! the tree, plus the Galactic tide), and the energies and potentials
!
   USE nbvars
   USE nbtree
   IMPLICIT NONE
!
   CONTAINS
! =====================================================================
   SUBROUTINE forces
!
! accelerations f0p and jerks f1p on the active particles ilist(1:nact)
! from all particles at the predicted positions rp, vp.
! Also returns tpair2(i) = min over j of the square of the pairwise
! free-fall time r^3/(mi+mj) and flyby time r^2/v^2, and rmin2(i), the
! squared (unsoftened) distance to the nearest body, for the collision
! check.
!
! With usetree the Barnes-Hut tree is used (module nbtree), otherwise an
! O(nact*N) direct sum, parallelised with OpenMP in one of two ways:
!  - big blocks (nact >= nompmin): threads share out the active stars,
!    each thread does whole i-sums (no communication)
!  - small blocks (most blocks with individual timesteps): the active
!    stars are done one at a time and threads share out the j-sum,
!    combined with reductions (if n >= nompj, otherwise serial)
! Newton's 3rd law is NOT used as it causes races between threads.
! The two copies of the pair interaction below must be kept identical.
!
   INTEGER :: ia,i,j
   DOUBLE PRECISION :: ax,ay,az,jx,jy,jz,dx,dy,dz,dvx,dvy,dvz,xi,yi,zi
   DOUBLE PRECISION :: vxi,vyi,vzi,mi,rad2,rinv,rinv2,mor,drdv,dv2,tmin2
   DOUBLE PRECISION :: d2u,rm2
!
   IF (usetree) THEN
! ---------------- Barnes-Hut tree: each walk is independent, so the
! active stars are shared over the threads (dynamic: walks differ)
!$OMP PARALLEL DO IF(nact >= nompmin) SCHEDULE(DYNAMIC,4) &
!$OMP& PRIVATE(ia,i,ax,ay,az,jx,jy,jz,tmin2,rm2)
     DO ia=1,nact
       i=ilist(ia)
       CALL treeforce(i,ax,ay,az,jx,jy,jz,tmin2,rm2)
       CALL finish(i,rp(1,i),rp(2,i),rp(3,i),vp(1,i),vp(2,i),vp(3,i), &
&        ax,ay,az,jx,jy,jz,tmin2,rm2)
     END DO
!$OMP END PARALLEL DO
!
   ELSE IF (nact >= nompmin) THEN
! ---------------- big block: parallel over active stars
!$OMP PARALLEL DO SCHEDULE(STATIC) &
!$OMP& PRIVATE(ia,i,j,ax,ay,az,jx,jy,jz,dx,dy,dz,dvx,dvy,dvz,xi,yi,zi, &
!$OMP&         vxi,vyi,vzi,mi,rad2,rinv,rinv2,mor,drdv,dv2,tmin2,d2u,rm2)
     DO ia=1,nact
       i=ilist(ia)
       xi=rp(1,i)
       yi=rp(2,i)
       zi=rp(3,i)
       vxi=vp(1,i)
       vyi=vp(2,i)
       vzi=vp(3,i)
       mi=m(i)
       ax=0.d0
       ay=0.d0
       az=0.d0
       jx=0.d0
       jy=0.d0
       jz=0.d0
       tmin2=HUGE(1.d0)
       rm2=HUGE(1.d0)
       DO j=1,n
         IF (i == j) CYCLE
! separations are j - i, so the force on i is +m(j)*dr/r^3
         dx=rp(1,j) - xi
         dy=rp(2,j) - yi
         dz=rp(3,j) - zi
         dvx=vp(1,j) - vxi
         dvy=vp(2,j) - vyi
         dvz=vp(3,j) - vzi
         d2u=dx*dx + dy*dy + dz*dz
         rm2=MIN(rm2,d2u)
         rad2=d2u + eps2
         rinv=1.d0/SQRT(rad2)
         rinv2=rinv*rinv
         mor=m(j)*rinv*rinv2
         drdv=3.d0*(dx*dvx + dy*dvy + dz*dvz)*rinv2
         ax=ax + mor*dx
         ay=ay + mor*dy
         az=az + mor*dz
         jx=jx + mor*(dvx - drdv*dx)
         jy=jy + mor*(dvy - drdv*dy)
         jz=jz + mor*(dvz - drdv*dz)
! pairwise timescales (squared): free-fall r^3/M and flyby r^2/v^2
         dv2=dvx*dvx + dvy*dvy + dvz*dvz
         tmin2=MIN(tmin2, rad2*rad2*rinv/(mi + m(j)), rad2/(dv2 + TINY(1.d0)))
       END DO
       CALL finish(i,xi,yi,zi,vxi,vyi,vzi,ax,ay,az,jx,jy,jz,tmin2,rm2)
     END DO
!$OMP END PARALLEL DO
!
   ELSE
! ---------------- small block: one active star at a time, parallel j-sum
     DO ia=1,nact
       i=ilist(ia)
       xi=rp(1,i)
       yi=rp(2,i)
       zi=rp(3,i)
       vxi=vp(1,i)
       vyi=vp(2,i)
       vzi=vp(3,i)
       mi=m(i)
       ax=0.d0
       ay=0.d0
       az=0.d0
       jx=0.d0
       jy=0.d0
       jz=0.d0
       tmin2=HUGE(1.d0)
       rm2=HUGE(1.d0)
!$OMP PARALLEL DO IF(n >= nompj) SCHEDULE(STATIC) &
!$OMP& PRIVATE(j,dx,dy,dz,dvx,dvy,dvz,rad2,rinv,rinv2,mor,drdv,dv2,d2u) &
!$OMP& REDUCTION(+:ax,ay,az,jx,jy,jz) REDUCTION(MIN:tmin2,rm2)
       DO j=1,n
         IF (i == j) CYCLE
         dx=rp(1,j) - xi
         dy=rp(2,j) - yi
         dz=rp(3,j) - zi
         dvx=vp(1,j) - vxi
         dvy=vp(2,j) - vyi
         dvz=vp(3,j) - vzi
         d2u=dx*dx + dy*dy + dz*dz
         rm2=MIN(rm2,d2u)
         rad2=d2u + eps2
         rinv=1.d0/SQRT(rad2)
         rinv2=rinv*rinv
         mor=m(j)*rinv*rinv2
         drdv=3.d0*(dx*dvx + dy*dvy + dz*dvz)*rinv2
         ax=ax + mor*dx
         ay=ay + mor*dy
         az=az + mor*dz
         jx=jx + mor*(dvx - drdv*dx)
         jy=jy + mor*(dvy - drdv*dy)
         jz=jz + mor*(dvz - drdv*dz)
         dv2=dvx*dvx + dvy*dvy + dvz*dvz
         tmin2=MIN(tmin2, rad2*rad2*rinv/(mi + m(j)), rad2/(dv2 + TINY(1.d0)))
       END DO
!$OMP END PARALLEL DO
       CALL finish(i,xi,yi,zi,vxi,vyi,vzi,ax,ay,az,jx,jy,jz,tmin2,rm2)
     END DO
   END IF
!
   CONTAINS
!
   SUBROUTINE finish(i,xi,yi,zi,vxi,vyi,vzi,ax,ay,az,jx,jy,jz,tmin2,rm2)
! add the Galactic tide and store the results for star i.
! Everything is passed as arguments (NOT host association) so each
! OpenMP thread works on its own private copies.
! tide (rotating frame): a = T.x + 2 omega (vy, -vx, 0)
! jerk = T.v + 2 omega (ay, -ax, 0) using the TOTAL acceleration,
! so all of a must be finished before the jerk terms are added
   INTEGER, INTENT(IN) :: i
   DOUBLE PRECISION, INTENT(IN) :: xi,yi,zi,vxi,vyi,vzi,tmin2,rm2
   DOUBLE PRECISION, INTENT(INOUT) :: ax,ay,az,jx,jy,jz
   IF (ltide) THEN
     ax=ax + ttid(1,1)*xi + ttid(1,2)*yi + ttid(1,3)*zi + 2.d0*omega*vyi
     ay=ay + ttid(2,1)*xi + ttid(2,2)*yi + ttid(2,3)*zi - 2.d0*omega*vxi
     az=az + ttid(3,1)*xi + ttid(3,2)*yi + ttid(3,3)*zi
     jx=jx + ttid(1,1)*vxi + ttid(1,2)*vyi + ttid(1,3)*vzi + 2.d0*omega*ay
     jy=jy + ttid(2,1)*vxi + ttid(2,2)*vyi + ttid(2,3)*vzi - 2.d0*omega*ax
     jz=jz + ttid(3,1)*vxi + ttid(3,2)*vyi + ttid(3,3)*vzi
   END IF
   f0p(1,i)=ax
   f0p(2,i)=ay
   f0p(3,i)=az
   f1p(1,i)=jx
   f1p(2,i)=jy
   f1p(3,i)=jz
   tpair2(i)=tmin2
   rmin2(i)=rm2
   END SUBROUTINE finish
!
   END SUBROUTINE forces
! =====================================================================
   SUBROUTINE energy(ke,pot,ptid)
!
! kinetic, gravitational potential (negative) and tidal potential
! energy in internal units, softened consistently. Only meaningful when
! all particles are synchronised. Multiply by vunit**2 for Msun (km/s)^2
! With a tide, ke + pot + ptid is the Jacobi energy, which is conserved
! (the Coriolis force does no work). Tidal potential = -x.T.x/2
! Also sets phi(i), the gravitational potential per unit mass of each
! star (the full i-j loop is twice the work of i<j, but parallelises
! with no races and gives phi for the snapshots)
!
   DOUBLE PRECISION, INTENT(OUT) :: ke,pot,ptid
   INTEGER :: i,j
   DOUBLE PRECISION :: poti,dx,dy,dz,xi,yi,zi,tx(1:3)
!
   ke=0.d0
   ptid=0.d0
!$OMP PARALLEL DO IF(n >= nompj) PRIVATE(i,tx) REDUCTION(+:ke,ptid) &
!$OMP& SCHEDULE(STATIC)
   DO i=1,n
     ke=ke + 0.5d0*m(i)*(v(1,i)**2 + v(2,i)**2 + v(3,i)**2)
     IF (ltide) THEN
       tx(1:3)=MATMUL(ttid,r(1:3,i))
       ptid=ptid - 0.5d0*m(i)*(r(1,i)*tx(1) + r(2,i)*tx(2) + r(3,i)*tx(3))
     END IF
   END DO
!$OMP END PARALLEL DO
!
!$OMP PARALLEL DO PRIVATE(i,j,poti,dx,dy,dz,xi,yi,zi) SCHEDULE(STATIC)
   DO i=1,n
     xi=r(1,i)
     yi=r(2,i)
     zi=r(3,i)
     poti=0.d0
     DO j=1,n
       IF (j == i) CYCLE
       dx=r(1,j) - xi
       dy=r(2,j) - yi
       dz=r(3,j) - zi
       poti=poti - m(j)/SQRT(dx*dx + dy*dy + dz*dz + eps2)
     END DO
     phi(i)=poti
   END DO
!$OMP END PARALLEL DO
   pot=0.5d0*SUM(m(1:n)*phi(1:n))
!
   END SUBROUTINE energy
!
   END MODULE nbforce
