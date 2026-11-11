! ==============================================================
! HENDOS - HErmite N-body Dynamics Of Systems
! PLUMMER SPHERE
!
!   makeplummer    positions and velocities of a Plummer sphere
!
! Called by the program makeics (makeics.f90). Uses ran3 (icrandom.f90).
! ==============================================================
!
   SUBROUTINE makeplummer(nstar,rstar,vstar,idum)
!
! positions and velocities of nstar stars in a Plummer sphere
! (Aarseth, Henon & Wielen 1974, A&A, 37, 183), cut at 10 Plummer radii.
!
! The sphere is returned in model units: G = 1, total mass = 1 and
! Plummer radius = 1, for stars of equal mass. Nothing here knows the
! masses or the physical size - makeics gives the stars their masses,
! scales the positions to the half-mass radius asked for and scales the
! velocities to the virial ratio.
!
! nstar (int) = number of stars
! rstar (3,n) = positions to return (Plummer radius = 1)
! vstar (3,n) = velocities to return (model units)
! idum = random number seed (already started: see ran3)
!
   IMPLICIT NONE
   INTEGER :: nstar,idum,i
   DOUBLE PRECISION :: rstar(1:3,nstar),vstar(1:3,nstar)
   DOUBLE PRECISION :: x1,x2,x3,rad,ve,ct,st,phi,pi
   REAL :: ran3
!
   pi=3.141592653589793d0
!
   DO i=1,nstar
! radius from the cumulative mass, cut at r=10
     DO
       x1=DBLE(ran3(idum))
       IF (x1 <= 0.d0) CYCLE
       rad=1.d0/SQRT(x1**(-2.d0/3.d0) - 1.d0)
       IF (rad < 10.d0) EXIT
     END DO
     x2=DBLE(ran3(idum))
     x3=DBLE(ran3(idum))
     ct=2.d0*x2 - 1.d0
     st=SQRT(1.d0 - ct*ct)
     phi=2.d0*pi*x3
     rstar(1,i)=rad*st*COS(phi)
     rstar(2,i)=rad*st*SIN(phi)
     rstar(3,i)=rad*ct
! speed by rejection from g(q) = q^2 (1-q^2)^3.5
     DO
       x1=DBLE(ran3(idum))
       x2=DBLE(ran3(idum))
       IF (0.1d0*x2 < x1*x1*(1.d0 - x1*x1)**3.5d0) EXIT
     END DO
     ve=SQRT(2.d0)*(1.d0 + rad*rad)**(-0.25d0)
     x2=DBLE(ran3(idum))
     x3=DBLE(ran3(idum))
     ct=2.d0*x2 - 1.d0
     st=SQRT(1.d0 - ct*ct)
     phi=2.d0*pi*x3
     vstar(1,i)=x1*ve*st*COS(phi)
     vstar(2,i)=x1*ve*st*SIN(phi)
     vstar(3,i)=x1*ve*ct
   END DO
!
   END SUBROUTINE makeplummer
!
