! ==============================================================
   MODULE nbtree
!
! HENDOS - HErmite N-body Dynamics Of Systems
! BARNES-HUT TREE (switch usetree) for large associations
!
! Octree over the bodies' predicted positions rp, with up to nleaf
! bodies per leaf. Each node has monopole, quadrupole and octupole
! moments about its centre of mass (raw moments S2_ij = sum m x_i x_j,
! S3_ijk = sum m x_i x_j x_k, x relative to the node c.o.m.), its c.o.m.
! velocity, and bmax, the radius about the c.o.m. containing its bodies.
!
! Opening criterion: a node is used as a whole if bmax/d < theta, d the
! distance from the field star to the node c.o.m.; otherwise it is
! opened. Leaves that are opened are summed directly (exactly as in the
! direct code, including softening).
! Acceleration from an accepted node (r from node c.o.m. to the star):
!   a = -M r/r^3
!       + (3 S2.r - s r)/r^5 - 5 (3 q - r^2 s) r / (2 r^7)
!       + (15 C - 6 r (T.r) - 3 r^2 T)/(2 r^7) - 7 (15 c - 9 r^2 T.r) r/(6 r^9)
!   s = tr S2, q = r.S2.r, C_i = S3_ijk r_j r_k, c = C.r, T_k = S3_iik
! Jerk from an accepted node: monopole (with the node c.o.m. velocity)
! plus the quadrupole term from the internal motion of the node's
! bodies, i.e. the quadrupole formula with S2 replaced by
! dS2/dt = sum m (x u + u x), u velocities relative to the node c.o.m.
! (tested: roughly halves the jerk error of the monopole alone).
!
! Nodes are numbered in depth-first order, so the first child of an
! internal node is node+1, and tnext(node) is the next node once the
! node's subtree is finished - the tree is walked without a stack.
! The bodies of each node are a contiguous range of tperm.
!
! With block timesteps the topology is rebuilt (treebuild) once, on
! average, every body has stepped, and after any binary is closed or
! opened; in between only the moments and bmax are recomputed from the
! new predicted positions (treemoments, O(N)), which keeps the opening
! criterion valid as bodies move.
!
! The Galactic tide is added to every star after the force sum, so it is
! included whether the tree is used or not. Energies are always direct.
!
   USE nbvars
   IMPLICIT NONE
   INTEGER :: nnode,nnodemax
   INTEGER, ALLOCATABLE :: tperm(:),tfirst(:),tcount(:),tnext(:)
   LOGICAL, ALLOCATABLE :: tleaf(:)
   DOUBLE PRECISION, ALLOCATABLE :: tm(:),tcom(:,:),tcv(:,:),tbmax(:)
   DOUBLE PRECISION, ALLOCATABLE :: ts2(:,:,:),ts3(:,:),ttr(:),tt3(:,:)
   DOUBLE PRECISION, ALLOCATABLE :: tsd(:,:,:),tsdtr(:)
   INTEGER, PARAMETER :: maxdepth=40
!
   CONTAINS
! =====================================================================
   SUBROUTINE treealloc
!
! node arrays: the depth-first build uses at most ~2N/nleaf + N chains;
! 3N + 1000 is generous, and overflow stops with a message
!
   nnodemax=3*ntot + 1000
   ALLOCATE(tperm(1:ntot),tfirst(1:nnodemax),tcount(1:nnodemax))
   ALLOCATE(tnext(1:nnodemax),tleaf(1:nnodemax))
   ALLOCATE(tm(1:nnodemax),tcom(1:3,1:nnodemax),tcv(1:3,1:nnodemax))
   ALLOCATE(tbmax(1:nnodemax),ts2(1:3,1:3,1:nnodemax),ts3(1:10,1:nnodemax))
   ALLOCATE(ttr(1:nnodemax),tt3(1:3,1:nnodemax))
   ALLOCATE(tsd(1:3,1:3,1:nnodemax),tsdtr(1:nnodemax))
!
   END SUBROUTINE treealloc
! =====================================================================
   SUBROUTINE treebuild
!
! build the tree topology for bodies 1..n at positions rp, then the
! moments
!
   INTEGER :: i
   DOUBLE PRECISION :: xmin(1:3),xmax(1:3),centre(1:3),half
!
   DO i=1,n
     tperm(i)=i
   END DO
   xmin=MINVAL(rp(1:3,1:n),DIM=2)
   xmax=MAXVAL(rp(1:3,1:n),DIM=2)
   centre=0.5d0*(xmin + xmax)
   half=0.5d0*MAXVAL(xmax - xmin)*1.0001d0 + TINY(1.d0)
   nnode=1
   CALL buildnode(1,1,n,centre,half,0)
   treevalid=.TRUE.
   ntreebuild=ntreebuild + 1
   CALL treemoments
!
   END SUBROUTINE treebuild
! =====================================================================
   RECURSIVE SUBROUTINE buildnode(node,ifirst,ilast,centre,half,depth)
!
! node holds bodies tperm(ifirst:ilast) in the cube centre +- half.
! Split into octants in place; children are created depth first.
!
   INTEGER, INTENT(IN) :: node,ifirst,ilast,depth
   DOUBLE PRECISION, INTENT(IN) :: centre(1:3),half
   INTEGER :: cnt(0:7),start(0:7),pos(0:7),oct,k,ip,c
   INTEGER, ALLOCATABLE :: tmp(:)
   DOUBLE PRECISION :: cc(1:3)
!
   tfirst(node)=ifirst
   tcount(node)=ilast - ifirst + 1
   IF (tcount(node) <= nleaf .OR. depth >= maxdepth) THEN
     tleaf(node)=.TRUE.
     tnext(node)=nnode + 1
     RETURN
   END IF
   tleaf(node)=.FALSE.
!
! sort the node's bodies by octant (counting sort)
   cnt=0
   DO k=ifirst,ilast
     oct=octant(rp(1:3,tperm(k)),centre)
     cnt(oct)=cnt(oct) + 1
   END DO
   start(0)=ifirst
   DO oct=1,7
     start(oct)=start(oct-1) + cnt(oct-1)
   END DO
   pos=start
   ALLOCATE(tmp(ifirst:ilast))
   DO k=ifirst,ilast
     ip=tperm(k)
     oct=octant(rp(1:3,ip),centre)
     tmp(pos(oct))=ip
     pos(oct)=pos(oct) + 1
   END DO
   tperm(ifirst:ilast)=tmp
   DEALLOCATE(tmp)
!
! children (non-empty octants only), depth first
   DO oct=0,7
     IF (cnt(oct) == 0) CYCLE
     nnode=nnode + 1
     IF (nnode > nnodemax) THEN
       WRITE(6,'(a)') ' ERROR: tree node overflow - try a larger nleaf'
       STOP 1
     END IF
     c=nnode
     cc(1)=centre(1) + SIGN(0.5d0*half,DBLE(IAND(oct,1)) - 0.5d0)
     cc(2)=centre(2) + SIGN(0.5d0*half,DBLE(IAND(oct,2)) - 0.5d0)
     cc(3)=centre(3) + SIGN(0.5d0*half,DBLE(IAND(oct,4)) - 0.5d0)
     CALL buildnode(c,start(oct),start(oct) + cnt(oct) - 1,cc,0.5d0*half,depth + 1)
   END DO
   tnext(node)=nnode + 1
!
   END SUBROUTINE buildnode
! =====================================================================
   INTEGER FUNCTION octant(x,centre)
   DOUBLE PRECISION, INTENT(IN) :: x(1:3),centre(1:3)
   octant=0
   IF (x(1) >= centre(1)) octant=octant + 1
   IF (x(2) >= centre(2)) octant=octant + 2
   IF (x(3) >= centre(3)) octant=octant + 4
   END FUNCTION octant
! =====================================================================
   SUBROUTINE treemoments
!
! moments, c.o.m. velocity and bmax of every node from the current rp,
! vp: leaves directly from their bodies (in parallel), internal nodes
! from their children with parallel-axis shifts, deepest first
! (children always have larger node numbers than their parent)
!
   INTEGER :: node,k,ip,c,ii,jj
   DOUBLE PRECISION :: x(1:3),u(1:3),mt,d(1:3),w(1:3),b
!
!$OMP PARALLEL DO IF(n >= nompj) SCHEDULE(DYNAMIC,64) PRIVATE(node,k,ip,x,u,mt,b,ii,jj)
   DO node=1,nnode
     IF (.NOT. tleaf(node)) CYCLE
     mt=0.d0
     tcom(1:3,node)=0.d0
     tcv(1:3,node)=0.d0
     DO k=tfirst(node),tfirst(node) + tcount(node) - 1
       ip=tperm(k)
       mt=mt + m(ip)
       tcom(1:3,node)=tcom(1:3,node) + m(ip)*rp(1:3,ip)
       tcv(1:3,node)=tcv(1:3,node) + m(ip)*vp(1:3,ip)
     END DO
     tm(node)=mt
     tcom(1:3,node)=tcom(1:3,node)/mt
     tcv(1:3,node)=tcv(1:3,node)/mt
     ts2(1:3,1:3,node)=0.d0
     ts3(1:10,node)=0.d0
     tsd(1:3,1:3,node)=0.d0
     b=0.d0
     DO k=tfirst(node),tfirst(node) + tcount(node) - 1
       ip=tperm(k)
       x=rp(1:3,ip) - tcom(1:3,node)
       u=vp(1:3,ip) - tcv(1:3,node)
       CALL addpoint(m(ip),x,ts2(1:3,1:3,node),ts3(1:10,node))
       DO jj=1,3
         DO ii=1,3
           tsd(ii,jj,node)=tsd(ii,jj,node) + m(ip)*(x(ii)*u(jj) + u(ii)*x(jj))
         END DO
       END DO
       b=MAX(b,SUM(x*x))
     END DO
     tbmax(node)=SQRT(b)
   END DO
!$OMP END PARALLEL DO
!
   DO node=nnode,1,-1
     IF (tleaf(node)) CYCLE
     mt=0.d0
     tcom(1:3,node)=0.d0
     tcv(1:3,node)=0.d0
     c=node + 1
     DO WHILE (c < tnext(node))
       mt=mt + tm(c)
       tcom(1:3,node)=tcom(1:3,node) + tm(c)*tcom(1:3,c)
       tcv(1:3,node)=tcv(1:3,node) + tm(c)*tcv(1:3,c)
       c=tnext(c)
     END DO
     tm(node)=mt
     tcom(1:3,node)=tcom(1:3,node)/mt
     tcv(1:3,node)=tcv(1:3,node)/mt
     ts2(1:3,1:3,node)=0.d0
     ts3(1:10,node)=0.d0
     tsd(1:3,1:3,node)=0.d0
     b=0.d0
     c=node + 1
     DO WHILE (c < tnext(node))
       d=tcom(1:3,c) - tcom(1:3,node)
       w=tcv(1:3,c) - tcv(1:3,node)
       CALL shiftadd(tm(c),d,ts2(1:3,1:3,c),ts3(1:10,c), &
&        ts2(1:3,1:3,node),ts3(1:10,node))
       DO jj=1,3
         DO ii=1,3
           tsd(ii,jj,node)=tsd(ii,jj,node) + tsd(ii,jj,c) + &
&            tm(c)*(d(ii)*w(jj) + w(ii)*d(jj))
         END DO
       END DO
       b=MAX(b,SQRT(SUM(d*d)) + tbmax(c))
       c=tnext(c)
     END DO
     tbmax(node)=b
   END DO
!
! derived quantities used in every evaluation: trace of S2 and T_k = S3_iik
   DO node=1,nnode
     ttr(node)=ts2(1,1,node) + ts2(2,2,node) + ts2(3,3,node)
     tsdtr(node)=tsd(1,1,node) + tsd(2,2,node) + tsd(3,3,node)
     tt3(1,node)=ts3(1,node) + ts3(4,node) + ts3(6,node)
     tt3(2,node)=ts3(2,node) + ts3(7,node) + ts3(9,node)
     tt3(3,node)=ts3(3,node) + ts3(8,node) + ts3(10,node)
   END DO
!
   END SUBROUTINE treemoments
! =====================================================================
   SUBROUTINE addpoint(mp,x,s2,s3)
!
! add a point mass mp at x (relative to the c.o.m.) to the raw moments.
! s3 packed: xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
!
   DOUBLE PRECISION, INTENT(IN) :: mp,x(1:3)
   DOUBLE PRECISION, INTENT(INOUT) :: s2(1:3,1:3),s3(1:10)
   INTEGER :: i,j
!
   DO j=1,3
     DO i=1,3
       s2(i,j)=s2(i,j) + mp*x(i)*x(j)
     END DO
   END DO
   s3(1)=s3(1) + mp*x(1)*x(1)*x(1)
   s3(2)=s3(2) + mp*x(1)*x(1)*x(2)
   s3(3)=s3(3) + mp*x(1)*x(1)*x(3)
   s3(4)=s3(4) + mp*x(1)*x(2)*x(2)
   s3(5)=s3(5) + mp*x(1)*x(2)*x(3)
   s3(6)=s3(6) + mp*x(1)*x(3)*x(3)
   s3(7)=s3(7) + mp*x(2)*x(2)*x(2)
   s3(8)=s3(8) + mp*x(2)*x(2)*x(3)
   s3(9)=s3(9) + mp*x(2)*x(3)*x(3)
   s3(10)=s3(10) + mp*x(3)*x(3)*x(3)
!
   END SUBROUTINE addpoint
! =====================================================================
   SUBROUTINE shiftadd(mc,d,s2c,s3c,s2,s3)
!
! add a child's moments (about its own c.o.m., which is at d from the
! parent c.o.m.) to the parent's: the parallel-axis theorem
!   S2_ij  += S2c_ij + mc d_i d_j
!   S3_ijk += S3c_ijk + d_i S2c_jk + d_j S2c_ik + d_k S2c_ij + mc d_i d_j d_k
!
   DOUBLE PRECISION, INTENT(IN) :: mc,d(1:3),s2c(1:3,1:3),s3c(1:10)
   DOUBLE PRECISION, INTENT(INOUT) :: s2(1:3,1:3),s3(1:10)
   INTEGER, PARAMETER :: ia(10)=(/1,1,1,1,1,1,2,2,2,3/)
   INTEGER, PARAMETER :: ja(10)=(/1,1,1,2,2,3,2,2,3,3/)
   INTEGER, PARAMETER :: ka(10)=(/1,2,3,2,3,3,2,3,3,3/)
   INTEGER :: l,i,j,k
!
   s2=s2 + s2c
   DO j=1,3
     DO i=1,3
       s2(i,j)=s2(i,j) + mc*d(i)*d(j)
     END DO
   END DO
   DO l=1,10
     i=ia(l)
     j=ja(l)
     k=ka(l)
     s3(l)=s3(l) + s3c(l) + d(i)*s2c(j,k) + d(j)*s2c(i,k) + d(k)*s2c(i,j) &
&         + mc*d(i)*d(j)*d(k)
   END DO
!
   END SUBROUTINE shiftadd
! =====================================================================
   SUBROUTINE treeforce(i,ax,ay,az,jx,jy,jz,tmin2,rm2)
!
! acceleration, jerk, pairwise timescale and nearest distance for body i
! at rp(i), vp(i) from the tree. Thread safe (all results in arguments).
!
   INTEGER, INTENT(IN) :: i
   DOUBLE PRECISION, INTENT(OUT) :: ax,ay,az,jx,jy,jz,tmin2,rm2
   INTEGER :: node,k,j
   DOUBLE PRECISION :: xi(1:3),vi(1:3),mi,rr(1:3),dv(1:3),d2,d,rinv,rinv2,mor
   DOUBLE PRECISION :: drdv,rad2,d2u,dv2,s,q,sr(1:3),cv(1:3),c,tr,r5,r7,r9,a(1:3)
   DOUBLE PRECISION :: theta2
!
   xi=rp(1:3,i)
   vi=vp(1:3,i)
   mi=m(i)
   ax=0.d0
   ay=0.d0
   az=0.d0
   jx=0.d0
   jy=0.d0
   jz=0.d0
   tmin2=HUGE(1.d0)
   rm2=HUGE(1.d0)
   theta2=theta*theta
   node=1
   DO WHILE (node <= nnode)
     rr=xi - tcom(1:3,node)
     d2=SUM(rr*rr)
     IF (tbmax(node)**2 < theta2*d2) THEN
! ---- accepted node: multipoles for a, monopole for the jerk
       d=SQRT(d2)
       rinv=1.d0/d
       rinv2=rinv*rinv
       r5=rinv2*rinv2*rinv
       r7=r5*rinv2
       r9=r7*rinv2
       a=-tm(node)*rr*rinv2*rinv
       s=ttr(node)
       sr=MATMUL(ts2(1:3,1:3,node),rr)
       q=SUM(rr*sr)
       a=a + (3.d0*sr - s*rr)*r5 - 2.5d0*(3.d0*q - d2*s)*rr*r7
       CALL octc(ts3(1:10,node),rr,cv)
       c=SUM(cv*rr)
       tr=SUM(tt3(1:3,node)*rr)
       a=a + 0.5d0*(15.d0*cv - 6.d0*rr*tr - 3.d0*d2*tt3(1:3,node))*r7 &
&          - 7.d0/6.d0*(15.d0*c - 9.d0*d2*tr)*rr*r9
       ax=ax + a(1)
       ay=ay + a(2)
       az=az + a(3)
       dv=vi - tcv(1:3,node)
       mor=tm(node)*rinv2*rinv
       drdv=3.d0*SUM(rr*dv)*rinv2
       jx=jx - mor*(dv(1) - drdv*rr(1))
       jy=jy - mor*(dv(2) - drdv*rr(2))
       jz=jz - mor*(dv(3) - drdv*rr(3))
! internal motion of the node: quadrupole formula with dS2/dt
       s=tsdtr(node)
       sr=MATMUL(tsd(1:3,1:3,node),rr)
       q=SUM(rr*sr)
       a=(3.d0*sr - s*rr)*r5 - 2.5d0*(3.d0*q - d2*s)*rr*r7
       jx=jx + a(1)
       jy=jy + a(2)
       jz=jz + a(3)
       dv2=SUM(dv*dv)
       tmin2=MIN(tmin2, d2*d/(mi + tm(node)), d2/(dv2 + TINY(1.d0)))
       rm2=MIN(rm2,(d - tbmax(node))**2)
       node=tnext(node)
     ELSE IF (tleaf(node)) THEN
! ---- opened leaf: direct sum, exactly as in the direct code
       DO k=tfirst(node),tfirst(node) + tcount(node) - 1
         j=tperm(k)
         IF (j == i) CYCLE
         rr=rp(1:3,j) - xi
         dv=vp(1:3,j) - vi
         d2u=SUM(rr*rr)
         rm2=MIN(rm2,d2u)
         rad2=d2u + eps2
         rinv=1.d0/SQRT(rad2)
         rinv2=rinv*rinv
         mor=m(j)*rinv*rinv2
         drdv=3.d0*SUM(rr*dv)*rinv2
         ax=ax + mor*rr(1)
         ay=ay + mor*rr(2)
         az=az + mor*rr(3)
         jx=jx + mor*(dv(1) - drdv*rr(1))
         jy=jy + mor*(dv(2) - drdv*rr(2))
         jz=jz + mor*(dv(3) - drdv*rr(3))
         dv2=SUM(dv*dv)
         tmin2=MIN(tmin2, rad2*rad2*rinv/(mi + m(j)), rad2/(dv2 + TINY(1.d0)))
       END DO
       node=tnext(node)
     ELSE
! ---- open the node: go to its first child
       node=node + 1
     END IF
   END DO
!
   END SUBROUTINE treeforce
! =====================================================================
   SUBROUTINE octc(s3,r,cv)
!
! C_i = S3_ijk r_j r_k from the packed octupole moments
!
   DOUBLE PRECISION, INTENT(IN) :: s3(1:10),r(1:3)
   DOUBLE PRECISION, INTENT(OUT) :: cv(1:3)
   DOUBLE PRECISION :: xx,yy,zz,xy,xz,yz
!
   xx=r(1)*r(1)
   yy=r(2)*r(2)
   zz=r(3)*r(3)
   xy=r(1)*r(2)
   xz=r(1)*r(3)
   yz=r(2)*r(3)
   cv(1)=s3(1)*xx + 2.d0*s3(2)*xy + 2.d0*s3(3)*xz + s3(4)*yy + 2.d0*s3(5)*yz + s3(6)*zz
   cv(2)=s3(2)*xx + 2.d0*s3(4)*xy + 2.d0*s3(5)*xz + s3(7)*yy + 2.d0*s3(8)*yz + s3(9)*zz
   cv(3)=s3(3)*xx + 2.d0*s3(5)*xy + 2.d0*s3(6)*xz + s3(8)*yy + 2.d0*s3(9)*yz + s3(10)*zz
!
   END SUBROUTINE octc
! =====================================================================
   SUBROUTINE treeperturb(xc,vc,mb,i1,i2,sum3,tp2)
!
! sum over bodies (except i1, i2) of m/d^3 at xc, and the squared
! shortest pairwise free-fall/flyby time, using monopoles of accepted
! nodes - for the binary perturbation gamma
!
   DOUBLE PRECISION, INTENT(IN) :: xc(1:3),vc(1:3),mb
   INTEGER, INTENT(IN) :: i1,i2
   DOUBLE PRECISION, INTENT(OUT) :: sum3,tp2
   INTEGER :: node,k,j
   DOUBLE PRECISION :: rr(1:3),d2,dv2,theta2
!
   sum3=0.d0
   tp2=HUGE(1.d0)
   theta2=theta*theta
   node=1
   DO WHILE (node <= nnode)
     rr=tcom(1:3,node) - xc
     d2=SUM(rr*rr)
     IF (tbmax(node)**2 < theta2*d2) THEN
       dv2=SUM((tcv(1:3,node) - vc)**2)
       sum3=sum3 + tm(node)/(d2*SQRT(d2))
       tp2=MIN(tp2, d2*SQRT(d2)/(mb + tm(node)), d2/(dv2 + TINY(1.d0)))
       node=tnext(node)
     ELSE IF (tleaf(node)) THEN
       DO k=tfirst(node),tfirst(node) + tcount(node) - 1
         j=tperm(k)
         IF (j == i1 .OR. j == i2) CYCLE
         rr=rp(1:3,j) - xc
         d2=SUM(rr*rr)
         dv2=SUM((vp(1:3,j) - vc)**2)
         sum3=sum3 + m(j)/(d2*SQRT(d2))
         tp2=MIN(tp2, d2*SQRT(d2)/(mb + m(j)), d2/(dv2 + TINY(1.d0)))
       END DO
       node=tnext(node)
     ELSE
       node=node + 1
     END IF
   END DO
!
   END SUBROUTINE treeperturb
! =====================================================================
   SUBROUTINE treenearest(i,jn)
!
! nearest body to body i, pruning nodes whose bounding sphere is
! further away than the best found so far
!
   INTEGER, INTENT(IN) :: i
   INTEGER, INTENT(OUT) :: jn
   INTEGER :: node,k,j
   DOUBLE PRECISION :: best,d,d2
!
   jn=0
   best=HUGE(1.d0)
   node=1
   DO WHILE (node <= nnode)
     d=SQRT(SUM((rp(1:3,i) - tcom(1:3,node))**2)) - tbmax(node)
     IF (d > 0.d0 .AND. d*d >= best) THEN
       node=tnext(node)
     ELSE IF (tleaf(node)) THEN
       DO k=tfirst(node),tfirst(node) + tcount(node) - 1
         j=tperm(k)
         IF (j == i) CYCLE
         d2=SUM((rp(1:3,j) - rp(1:3,i))**2)
         IF (d2 < best) THEN
           best=d2
           jn=j
         END IF
       END DO
       node=tnext(node)
     ELSE
       node=node + 1
     END IF
   END DO
!
   END SUBROUTINE treenearest
!
   END MODULE nbtree
