! ==============================================================
! HENDOS - HErmite N-body Dynamics Of Systems
! BOX FRACTAL GENERATORS (Goodwin & Whitworth 2004, A&A, 413, 929)
!
!   makefractal     box fractal with inherited or random velocities
!   larsonfractal   box fractal with velocities that follow Larson's
!                   law (it calls makefractal for the positions)
!
! These are the Sheffield routines makefractal and makefractal2,
! rearranged: makefractal makes the positions and all the velocity
! fields that are built from the fractal itself, and larsonfractal adds
! the correlated (Larson) velocities of makefractal2. They are called
! by the program in makefractal.f90, use ran3 and gasdev
! (icrandom.f90) and work in single precision.
!
! What came from where:
!   positions   coresize > 0   as makefractal2
!               coresize = 0   as the original makefractal
!   velocities  makefractal    vtype 'f','u'   as makefractal2
!                              vtype 's'       as the original makefractal
!               larsonfractal  as makefractal2 with vtype 'c'
! Lines marked '! HENDOS:' differ from the routines as supplied.
! ==============================================================
!
   SUBROUTINE makefractal(nstar,rstar,vstar,fdim,ndiv,idum, radius, vscale, coresize, vtype)
!
! makes fractal of nstar particles of dimension fdim, using ndiv
! subunits forcing the number of cells if force=.true.
! Simon Goodwin & Ant Whitworth (2004, A&A, 413, 929)
!
! fdim (real) = fractal dimension
! ndiv (int) = number of divisions of box along each side (>=2), so
!              each box has ndiv**3 children. 2 is the usual value.
!              HENDOS: now an argument (it was set to 2 inside)
! nstar (int) = number of stars in *final* distribution
! rstar (3,n) = positions to return
! vstar (3,n) = velocities to return
! force (logical) = make box dimension such that integer no. children
! radius (real) = radius of the sphere (half the side of the first box)
! vscale (real) = scale of the velocities
! coresize (real) = > 0: boxes are divided until the children of the
!                   next level would be closer together than this, and
!                   each position is given a Gaussian offset of this
!                   size at every level
!                   = 0: boxes are divided until there are more than
!                   3.5*nstar of them, and each position is given a
!                   uniform offset of up to a quarter of the width of
!                   its box at every level
! vtype (character) = velocity field:
!     'f'  children inherit their parent's velocity plus a random part
!          proportional to the size of their box
!     's'  children inherit their parent's velocity plus a random part
!          that is the same size at every level (the original makefractal)
!     'u'  uncorrelated: Gaussian of dispersion vscale for every star
!     'n'  none: all velocities zero (used by larsonfractal, which makes
!          its own)
!     (for velocities that follow Larson's law use larsonfractal)
! idum = random number seed
   IMPLICIT NONE
   INTEGER :: ndiv,nstar,idum,isnext,iforce
   INTEGER :: is,nsmax,islast,flagabort,k
   INTEGER :: k1,k2,k3,i,isfirst,nsnow,nalloc
   REAL :: fdim,dum,delta,dr(3),prob,rad,cov(3)
   REAL :: how,how2,ran3,com(3), radius, vscale
   REAL :: gasdev, coresize
   LOGICAL :: force,cutsphere,inherit
   INTEGER :: n,isurv,level
   REAL :: rstar(1:3,nstar),vstar(1:3,nstar)
   INTEGER, DIMENSION(:), ALLOCATABLE :: forced
   REAL, DIMENSION(:,:), ALLOCATABLE :: r,v
   CHARACTER*1 :: vtype
!
! set initial parameters/settings
! HENDOS: the starting level (lstart, nstart) of makefractal2 is gone -
! it was set and printed but never used; the fractal always starts
! from one box.
!
! enforce fractal dimension
   force=.TRUE.
!
! probability of a child surviving
   prob=REAL(ndiv)**(fdim-3.)
!
! do you want to cut a sphere out?  It'll be a cube if you don't.
   cutsphere=.true.
!
! HENDOS: only the velocity fields this routine makes
   IF (vtype/='f' .AND. vtype/='s' .AND. vtype/='u' .AND. vtype/='n') &
&    STOP 'makefractal: vtype must be f, s, u or n (use larsonfractal for Larson velocities)'
!
! HENDOS: do the children start from their parent's velocity?
   inherit=(vtype=='f' .OR. vtype=='s')
!
! HENDOS: at least nstar*50, enough for the last level whatever ndiv
! is, and no division by zero if coresize=0
   nalloc=MAX(nstar*50,INT(3.5*nstar + 1.)*(ndiv**3 + 1) + 8)
   IF (coresize>0.) nalloc=MAX(nalloc,INT((radius/coresize)**fdim*50.))
   ALLOCATE(r(1:3,1:nalloc))
   ALLOCATE(v(1:3,1:nalloc))
   ALLOCATE(forced(1:ndiv**3))
!
! ===================================================================
!
! list characteristics
   WRITE(6,*) 'fractal dimension', fdim
   WRITE(6,*) 'ndiv', ndiv
   WRITE(6,*) 'cluster radius', radius
   WRITE(6,*) 'core size', coresize
   WRITE(6,*) 'velocity type ', vtype
   WRITE(6,*) 'number of stars', nstar
   WRITE(6,*) 'survival chance', prob
   WRITE(6,*) 'forcing?', force
   if (cutsphere) then
     write(6,*) 'Will cut out a sphere'
   else
     write(6,*) 'Will leave as a cube'
   end if
!
! warning
! HENDOS: was REAL(start)**fdim in makefractal2 - start is not a
! variable. The number of children of each box is ndiv**fdim.
   how=REAL(ndiv)**fdim
   WRITE(6,*) '# of children expected', how
   how2=how - INT(how)
   IF (how2>0.2 .AND. how2<0.8) THEN
     WRITE(6,*) 'WARNING:'
     IF (how2<0.5) THEN
       dum=LOG(REAL(INT(how)))/LOG(REAL(ndiv))
       WRITE(6,*) 'could look like fdim=', dum
     ELSE
       dum=LOG(REAL(INT(how)+1))/LOG(REAL(ndiv))
       WRITE(6,*) 'could look like fdim=', dum
     END IF
   END IF
!
! fractal dimension this will probably be
   IF (force) THEN
     iforce=INT(how + 0.5)
     WRITE(6,*) 'forcing to survivors=', iforce
     dum=LOG(REAL(iforce))/LOG(REAL(ndiv))
     WRITE(6,*) 'which means fdim=', dum
   END IF
!
! build fractal cube
   r=0.
   v=0.
   nsmax=1
   islast=1
   flagabort=0
   delta=radius/REAL(ndiv)
   level=0
   DO
     level=level + 1
! HENDOS: stop rather than run off the end of the arrays
     IF (nsmax*(iforce + 1) + 1 > nalloc) STOP 'too many boxes for the arrays - use a larger coresize'
     isfirst=nsmax + 1
     isnext=isfirst
! loop over all parents
     DO is=1,nsmax
!
! if forced decide on survivors now
       IF (force) THEN
         forced=0
         DO n=1,iforce
4          isurv=INT(ran3(idum)*REAL(ndiv**3)) + 1
           IF (forced(isurv)==1) GOTO 4
           forced(isurv)=1
         END DO
       END IF
! loop over all ndiv**3 possible children
       n=0
       DO k1=1,ndiv
         dr(1)=REAL(2*k1-1-ndiv)*delta
         DO k2=1,ndiv
           dr(2)=REAL(2*k2-1-ndiv)*delta
           DO k3=1,ndiv
             n=n+1
             dr(3)=REAL(2*k3-1-ndiv)*delta
             r(1:3,isnext)=r(1:3,is) + dr(1:3)
             IF (inherit) THEN
! keep parent's velocity
               v(1:3,isnext)=v(1:3,is)
             END IF
! do children survive?
! forced or not forced here
             IF (force) THEN
               IF (forced(n)==1) THEN
                 islast=isnext
                 isnext=isnext + 1
               END IF
             ELSE
               IF (ran3(idum)<prob) THEN
                 islast=isnext
                 isnext=isnext + 1
               END IF
             END IF
           END DO
         END DO
       END DO
     END DO
! are there children (we hope so)
     IF (islast/=nsmax) THEN
! delete parents
       DO is=isfirst,islast
         r(1:3,is-nsmax)=r(1:3,is)
         v(1:3,is-nsmax)=v(1:3,is)
       END DO
! new total number
       nsmax=islast - nsmax
       DO is=nsmax+1,islast
         r(1:3,is)=0.
         v(1:3,is)=0.
       END DO
! resize delta
       delta=delta/REAL(ndiv)
     ELSE
       flagabort=1
       EXIT
     END IF
     IF (flagabort==1) STOP 'no children'
! add noise to positions
! HENDOS: with a core size the noise in makefractal2 was
! (2.*gasdev(idum) - 1.) * coresize / 2.  The '- 1.' (left from a
! uniform deviate) shifted every star by -coresize/2 at every level;
! it is removed, the spread is the same. With no core size the noise
! is that of the original makefractal.
     DO is=1,nsmax
       DO k=1,3
         IF (coresize>0.) THEN
           r(k,is)=r(k,is) + gasdev(idum) * coresize
         ELSE
           r(k,is)=r(k,is) + (2.*ran3(idum) - 1.)*delta
         END IF
       END DO
     END DO
! add noise to velocities ('u' is set once the stars are chosen)
     IF (vtype=='f') THEN
! proportional to the box size (makefractal2)
       DO is=1,nsmax
         DO k=1,3
           v(k,is)=v(k,is) + delta * gasdev(idum) * vscale
         END DO
       END DO
     ELSE IF (vtype=='s') THEN
! the same at every level (original makefractal)
       DO is=1,nsmax
         DO k=1,3
           v(k,is)=v(k,is) + 2.**(fdim-1.) * gasdev(idum) * vscale
         END DO
       END DO
     END IF
     WRITE(6,*) 'done', delta, nsmax
     IF (coresize>0.) THEN
! Michiel alternative: exit if distance at next step would be lower than core size
       IF (delta * REAL(ndiv) < coresize) EXIT
     ELSE
! exit if enough particles
       IF (nsmax>3.5*nstar) EXIT
     END IF
   END DO
   WRITE(6,*) 'Made initial cubic distribution'
!
! cut out sphere
! HENDOS: only if cutsphere is set (it was always cut), and the stars
! inside are moved down the arrays in order, positions and velocities
! together. The original swap loops lost the last star, and the one in
! makefractal moved the positions without the velocities.
   IF (cutsphere) THEN
     nsnow=0
     DO is=1,nsmax
       rad=r(1,is)**2 + r(2,is)**2 + r(3,is)**2
       IF (rad<=radius ** 2) THEN
         nsnow=nsnow + 1
         r(1:3,nsnow)=r(1:3,is)
         v(1:3,nsnow)=v(1:3,is)
       END IF
     END DO
     nsmax=nsnow
     WRITE(6,*) 'Sphere cut out'
   END IF
! HENDOS: message says what to do
   IF (nsmax<nstar) STOP 'nsmax<nstar - use a smaller coresize or change random number seed'
!
! we now have a fractal of nsmax stars
! need to randomly remove stars until we get nstar
! HENDOS: test at the top, so nothing is removed if nsmax=nstar already
   DO
     IF (nsmax==nstar) EXIT
! star to remove
     is=INT(ran3(idum)*REAL(nsmax)) + 1
     r(1:3,is)=r(1:3,nsmax)
     v(1:3,is)=v(1:3,nsmax)
     nsmax=nsmax - 1
   END DO
   WRITE(6,*) 'Cut to size', nsmax
!
! now put r in rstar array
   DO is=1,nstar
     rstar(1:3,is)=r(1:3,is)
     vstar(1:3,is)=v(1:3,is)
   END DO
!
   IF (vtype == 'u') THEN
     DO is=1,nstar
       DO k=1,3
         vstar(k,is)=gasdev(idum) * vscale
       END DO
     END DO
   END IF


!
! set coms to zero
   com=0.
   cov=0.
   DO i=1,nstar
     com(1:3)=com(1:3) + rstar(1:3,i)
     cov(1:3)=cov(1:3) + vstar(1:3,i)
   END DO
   com=com/REAL(nstar)
   cov=cov/REAL(nstar)
   WRITE(6,*) 'COM', com(1), com(2), com(3)
   WRITE(6,*) 'COV', cov(1), cov(2), cov(3)
   DO i=1,nstar
     rstar(1:3,i)=rstar(1:3,i) - com(1:3)
     vstar(1:3,i)=vstar(1:3,i) - cov(1:3)
   END DO
!
   END SUBROUTINE makefractal
!
! ==============================================================
!
   SUBROUTINE larsonfractal(nstar,rstar,vstar,fdim,ndiv,idum, radius, vscale, coresize)
!
! makes a box fractal of nstar particles whose velocities follow
! Larson's law: the velocity difference between two stars grows as the
! square root of their separation.
!
! The positions are made by makefractal (same fdim, ndiv, radius and
! coresize - see there). The velocities are then drawn together as
! Gaussian random numbers whose covariance falls linearly with
! separation,  <v_i v_j> = vscale**2 * (D - d_ij/2)  for each component,
! with D ten times the largest separation, which gives
! <(v_i - v_j)**2> = vscale**2 * d_ij. This is the vtype='c' option of
! makefractal2. The method is written up in larson.tex.
!
! It needs an nstar x nstar matrix: memory 8*nstar**2 bytes, and a
! time that grows as nstar**3.
!
! nstar (int) = number of stars
! rstar (3,n) = positions to return
! vstar (3,n) = velocities to return
! fdim (real) = fractal dimension
! ndiv (int) = number of divisions of box along each side (usually 2)
! idum = random number seed
! radius (real) = radius of the sphere
! vscale (real) = scale of the velocities
! coresize (real) = core size (0 = none), as in makefractal
   IMPLICIT NONE
   INTEGER :: nstar,ndiv,idum
   INTEGER :: is,js,ks,i
   REAL :: fdim,radius,vscale,coresize,gasdev,cov(3)
   REAL :: max_distance, distance, larson_random
   REAL :: rstar(1:3,nstar),vstar(1:3,nstar)
! HENDOS: allocatable (they were automatic arrays, which overflow the
! stack for large nstar)
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: larson_correlate
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: larson_draw
!
! the fractal itself, with no velocities
   CALL makefractal(nstar,rstar,vstar,fdim,ndiv,idum,radius,vscale,coresize,'n')
!
   ! find the maximum distance between two stars to set the absolute scaling of the velocity dispersion
   WRITE(6, *) 'setting Larson velocities'
   ALLOCATE(larson_correlate(1:nstar,1:nstar))
   ALLOCATE(larson_draw(1:nstar))
   max_distance=0
   DO is=1,nstar
     DO js=1,nstar
       distance = SQRT((rstar(1, is) - rstar(1, js)) ** 2. + (rstar(2, is) - rstar(2, js)) ** 2.&
&                        + (rstar(3, is) - rstar(3, js)) ** 2.)
       IF (distance > max_distance) THEN
         max_distance = distance
       END IF
     END DO
   END DO
   max_distance = max_distance * 10.

   ! build the correlation matrix
   DO is=1,nstar
     DO js=1,nstar
       IF (is==js) THEN
         larson_correlate(is,js) = vscale ** 2. * max_distance
       ELSE
         distance = SQRT((rstar(1, is) - rstar(1, js)) ** 2. + (rstar(2, is) - rstar(2, js)) ** 2. &
              &+ (rstar(3, is) - rstar(3, js)) ** 2.)
         !larson_correlate(is,js) = vscale ** 2. * SQRT(max_distance ** 2. -&
         !                     &(max_distance * distance))
         larson_correlate(is,js) = vscale ** 2. * (max_distance - distance / 2.)
         larson_correlate(js,is) = larson_correlate(is, js)
       END IF
     END DO
   END DO

   ! Cholesky decomposition
   do is = 1, nstar
      larson_draw(is) = larson_correlate(is, is)
      larson_correlate(is, is) = 1

      do js = is+1, nstar
          larson_correlate(js, is) = larson_correlate(js, is) / larson_draw(is)
          larson_correlate(is, js) = 0.
      end do

      do js = is+1, nstar
          do ks = js, nstar
              larson_correlate(ks, js) = larson_correlate(ks, js) - larson_correlate(ks, is)&
                      &* larson_correlate(js, is) * larson_draw(is)
          end do
      end do
   end do

   ! set velocities to zero
   do is = 1, nstar
     do js=1,3
       vstar(js, is) = 0
     end do
   end do

   ! draw velocities
   ! HENDOS: MAX(..,0) guards the square root against rounding error
   do is = 1, nstar
     do js = 1, 3
       larson_random = gasdev(idum) * SQRT(MAX(larson_draw(is),0.d0))
       do ks = 1, nstar
         vstar(js, ks) = vstar(js, ks) + REAL(larson_random * larson_correlate(ks, is))
         !write(6, *) is, ks, larson_correlate(is, ks)
       enddo
     enddo
   enddo
   DEALLOCATE(larson_correlate,larson_draw)
!
! set cov to zero
   cov=0.
   DO i=1,nstar
     cov(1:3)=cov(1:3) + vstar(1:3,i)
   END DO
   cov=cov/REAL(nstar)
   WRITE(6,*) 'COV', cov(1), cov(2), cov(3)
   DO i=1,nstar
     vstar(1:3,i)=vstar(1:3,i) - cov(1:3)
   END DO
!
   END SUBROUTINE larsonfractal
!
