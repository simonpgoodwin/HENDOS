! ==============================================================
! HENDOS - HErmite N-body Dynamics Of Systems
! YOUR OWN INITIAL CONDITIONS
!
! This is the file to edit to make initial conditions your own way
! without touching the rest of the generator. It has three routines,
! each called by the program makeics (makeics.f90) at one point:
!
!   userstructure   your own positions and velocities
!                   used when  distribution = 'user'  in ics.nml
!   usermasses      your own masses
!                   used when  imf = 'user'  in ics.nml
!   userfinal       any last changes, to anything
!                   ALWAYS called, after everything else and just
!                   before the file is written. As supplied it does
!                   nothing.
!
! As supplied, userstructure makes a uniform sphere and usermasses a
! Salpeter mass function: they work, and are there to be replaced.
!
! After editing, type  make  in HENDOS/run.
!
! Random numbers:   x=DBLE(ran3(idum))     uniform between 0 and 1
!                   x=DBLE(gasdev(idum))   Gaussian, mean 0, variance 1
! Always pass idum on and never set it yourself: it carries the
! sequence started from the seed in ics.nml, so that one seed always
! gives the same stars.
!
! See src/ics/ICS_README.md for more, including worked examples.
! ==============================================================
!
   SUBROUTINE userstructure(nstar,rstar,vstar,idum)
!
! YOUR OWN STRUCTURE - called in STEP 1 when distribution = 'user'
!
! Return a position and a velocity for each of the nstar stars.
!
! nstar (int) = number of stars                               (given)
! rstar (3,n) = positions to return, in units in which the length that
!               radius (in ics.nml) refers to is 1. makeics multiplies
!               them by radius to get pc.
! vstar (3,n) = velocities to return. Any units if the velocities are
!               to be scaled to the virial ratio (the standard); km/s
!               if setvirial = .false. in ics.nml.
! idum = random number state                                  (pass on)
!
! The masses are not known here: every star counts as mass 1. They are
! given in STEP 2, and makeics then moves everything to the centre of
! mass and scales the velocities.
!
! AS SUPPLIED: a uniform sphere of radius 1 with an isotropic Gaussian
! velocity for each star.
!
   IMPLICIT NONE
   INTEGER :: nstar,idum,i
   DOUBLE PRECISION :: rstar(1:3,nstar),vstar(1:3,nstar)
   DOUBLE PRECISION :: rad,ct,st,phi,pi
   REAL :: ran3,gasdev
!
   pi=3.141592653589793d0
!
   DO i=1,nstar
! radius: uniform in volume, so the cube root of a uniform number
     rad=DBLE(ran3(idum))**(1.d0/3.d0)
! a random direction
     ct=2.d0*DBLE(ran3(idum)) - 1.d0
     st=SQRT(1.d0 - ct*ct)
     phi=2.d0*pi*DBLE(ran3(idum))
     rstar(1,i)=rad*st*COS(phi)
     rstar(2,i)=rad*st*SIN(phi)
     rstar(3,i)=rad*ct
! velocity: Gaussian in each direction
     vstar(1,i)=DBLE(gasdev(idum))
     vstar(2,i)=DBLE(gasdev(idum))
     vstar(3,i)=DBLE(gasdev(idum))
   END DO
!
   END SUBROUTINE userstructure
!
! ==============================================================
!
   SUBROUTINE usermasses(nstar,m,rstar,idum)
!
! YOUR OWN MASSES - called in STEP 2 when imf = 'user'
!
! Return a mass in Msun for each of the nstar stars.
!
! nstar (int) = number of stars                               (given)
! m (n) = masses to return (Msun)
! rstar (3,n) = positions, in model units (not yet pc, but the shape
!               is right), in case the mass is to depend on where the
!               star is. Do not change them here.             (given)
! idum = random number state                                  (pass on)
!
! AS SUPPLIED: a Salpeter power law, dN/dm ~ m**(-2.35), between mlow
! and mhigh, set just below.
!
   IMPLICIT NONE
   INTEGER :: nstar,idum,i
   DOUBLE PRECISION :: m(nstar),rstar(1:3,nstar)
   DOUBLE PRECISION :: slope,mlow,mhigh,x,p
   REAL :: ran3
!
   slope=2.35d0
   mlow=0.5d0
   mhigh=50.d0
!
   p=1.d0 - slope
   DO i=1,nstar
     x=DBLE(ran3(idum))
     m(i)=(mlow**p + x*(mhigh**p - mlow**p))**(1.d0/p)
   END DO
!
   END SUBROUTINE usermasses
!
! ==============================================================
!
   SUBROUTINE userfinal(nstar,m,rstar,vstar,ptype,idum)
!
! YOUR OWN LAST CHANGES - always called, in STEP 6
!
! Everything is finished and in its final units when this is called:
! the stars are in their centre-of-mass frame and their velocities have
! been scaled. Whatever you do here is what is written to the file.
! Nothing is rescaled afterwards - if you change masses or velocities
! the virial ratio changes with them (the summary on the screen shows
! what you ended up with), and if you want the centre of mass back at
! the origin you must put it there.
!
! nstar (int) = number of stars                               (given)
! m (n) = masses (Msun)
! rstar (3,n) = positions (pc)
! vstar (3,n) = velocities (km/s)
! ptype (n) = particle types: 'star', 'white_dwarf', 'neutron_star',
!             'black_hole' or 'planet'
! idum = random number state                                  (pass on)
!
! AS SUPPLIED: does nothing. Some things it could do are below - take
! out the leading '!' of the lines of one to use it.
!
   IMPLICIT NONE
   INTEGER :: nstar,idum
   DOUBLE PRECISION :: m(nstar),rstar(1:3,nstar),vstar(1:3,nstar)
   CHARACTER(len=12) :: ptype(nstar)
   INTEGER :: i,imax,imin
   DOUBLE PRECISION :: r2,r2min,omega,tmp
   REAL :: ran3,gasdev
!
! ---- EXAMPLE 1: make every star above 8 Msun a black hole -----------
!   DO i=1,nstar
!     IF (m(i) > 8.d0) ptype(i)='black_hole'
!   END DO
!
! ---- EXAMPLE 2: put the most massive star at the centre -------------
! (swap its mass with that of the star nearest the origin)
!   imax=1
!   imin=1
!   r2min=HUGE(1.d0)
!   DO i=1,nstar
!     IF (m(i) > m(imax)) imax=i
!     r2=rstar(1,i)**2 + rstar(2,i)**2 + rstar(3,i)**2
!     IF (r2 < r2min) THEN
!       r2min=r2
!       imin=i
!     END IF
!   END DO
!   tmp=m(imax)
!   m(imax)=m(imin)
!   m(imin)=tmp
!
! ---- EXAMPLE 3: add solid-body rotation about the z axis ------------
! (omega in km/s per pc)
!   omega=0.5d0
!   DO i=1,nstar
!     vstar(1,i)=vstar(1,i) - omega*rstar(2,i)
!     vstar(2,i)=vstar(2,i) + omega*rstar(1,i)
!   END DO
!
   END SUBROUTINE userfinal
!
