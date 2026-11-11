! ==============================================================
! HENDOS - HErmite N-body Dynamics Of Systems
! STELLAR MASSES
!
!   maschberger    masses drawn from the Maschberger (2013) IMF
!
! Called by the program makeics (makeics.f90). Uses ran3 (icrandom.f90).
! ==============================================================
!
   SUBROUTINE maschberger(nstar,m,idum,alpha,beta,mu,mlo,mup)
!
! masses of nstar single stars from the L3 initial mass function of
! Maschberger (2013, MNRAS, 429, 1725),
!    p(m) ~ (m/mu)**(-alpha) * (1 + (m/mu)**(1-alpha))**(-beta)
! between mlo and mup. It is a Salpeter-like power law of slope alpha
! at high masses that turns over smoothly near mu, and it can be
! sampled exactly from one random number pl3 with
!    G(m) = (1 + (m/mu)**(1-alpha))**(1-beta)
!    m = mu * ((pl3*(G(mup) - G(mlo)) + G(mlo))**(1/(1-beta)) - 1)**(1/(1-alpha))
! This is the Sheffield sampling code, in double precision.
!
! nstar (int) = number of stars
! m (n) = masses to return (same units as mu, mlo, mup: Msun)
! idum = random number seed (already started: see ran3)
! alpha = high-mass slope (2.3)
! beta = low-mass shape parameter (1.4)
! mu = scale mass (0.2 Msun)
! mlo, mup = lowest and highest mass
!
   IMPLICIT NONE
   INTEGER :: nstar,idum,i
   DOUBLE PRECISION :: m(nstar),alpha,beta,mu,mlo,mup
   DOUBLE PRECISION :: gmup,gmlo,pl3,mmm
   REAL :: ran3
!
   gmup=(1.d0 + (mup/mu)**(1.d0-alpha))**(1.d0-beta)
   gmlo=(1.d0 + (mlo/mu)**(1.d0-alpha))**(1.d0-beta)
!
   DO i=1,nstar
     pl3=DBLE(ran3(idum))
     mmm=pl3*(gmup - gmlo) + gmlo
     mmm=mmm**(1.d0/(1.d0-beta))
     mmm=mmm - 1.d0
     mmm=mmm**(1.d0/(1.d0-alpha))
     m(i)=mmm*mu
! rounding must not take a mass outside the limits
     m(i)=MIN(MAX(m(i),mlo),mup)
   END DO
!
   END SUBROUTINE maschberger
!
