! ==============================================================
! HENDOS - HErmite N-body Dynamics Of Systems
! RANDOM NUMBERS for the initial conditions generators
!
!   ran3(idum)    uniform deviate between 0 and 1 (never exactly 1)
!   gasdev(idum)  Gaussian deviate, zero mean and unit variance,
!                 made from ran3
!
! Every generator in HENDOS takes its random numbers from ran3, so one
! seed always gives the same initial conditions on any machine.
! ==============================================================
      REAL FUNCTION ran3(idum)
!
!  Uniform random number by Knuth's subtractive method (Knuth, The
!  Art of Computer Programming vol. 2, section 3.2.2) - the method of
!  the routine ran3 in Numerical Recipes, with the same calling
!  convention. Written for HENDOS from the published algorithm.
!
!  Call it first with a NEGATIVE idum to start (or restart) the
!  sequence; idum is then set to 1 and should be left alone.
!
!  The state is a table of 55 integers below mbig. Each call replaces
!  one entry by its difference from the entry 31 places further on.
!
      IMPLICIT NONE
      INTEGER :: idum
      INTEGER, PARAMETER :: mbig=1000000000, mseed=161803398
      INTEGER, SAVE :: ma(1:55), inext, inextp
      LOGICAL, SAVE :: started=.FALSE.
      INTEGER :: i,ii,k,mj,mk
!
      IF (idum < 0 .OR. .NOT. started) THEN
          started = .TRUE.
!  fill the table from the seed, in a scrambled order
          mj = MOD(ABS(mseed - ABS(idum)),mbig)
          ma(55) = mj
          mk = 1
          DO i = 1,54
              ii = MOD(21*i,55)
              ma(ii) = mk
              mk = mj - mk
              IF (mk < 0) mk = mk + mbig
              mj = ma(ii)
          END DO
!  run the table through four times to mix it
          DO k = 1,4
              DO i = 1,55
                  ma(i) = ma(i) - ma(1 + MOD(i + 30,55))
                  IF (ma(i) < 0) ma(i) = ma(i) + mbig
              END DO
          END DO
          inext = 0
          inextp = 31
          idum = 1
      END IF
!
!  next number: step both pointers round the table and subtract
      inext = inext + 1
      IF (inext == 56) inext = 1
      inextp = inextp + 1
      IF (inextp == 56) inextp = 1
      mj = ma(inext) - ma(inextp)
      IF (mj < 0) mj = mj + mbig
      ma(inext) = mj
      ran3 = REAL(mj)/REAL(mbig)
!  single precision can round the largest values up to 1
      IF (ran3 >= 1.) ran3 = NEAREST(1.,-1.)
!
      RETURN
      END FUNCTION ran3
! ==============================================================
      REAL FUNCTION gasdev(idum)
!  Gaussian random number, zero mean and unit variance, from two
!  uniform deviates by the polar Box-Muller method. Each pair of
!  uniform deviates gives two Gaussian ones: the second is kept
!  for the next call.
      IMPLICIT NONE
      INTEGER :: idum
      REAL :: ran3,v1,v2,rsq,fac
      REAL, SAVE :: spare
      LOGICAL, SAVE :: have=.FALSE.
!
      IF (have) THEN
          gasdev = spare
          have = .FALSE.
          RETURN
      END IF
      DO
          v1 = 2.*ran3(idum) - 1.
          v2 = 2.*ran3(idum) - 1.
          rsq = v1*v1 + v2*v2
          IF (rsq < 1. .AND. rsq > 0.) EXIT
      END DO
      fac = SQRT(-2.*LOG(rsq)/rsq)
      spare = v1*fac
      have = .TRUE.
      gasdev = v2*fac
!
      RETURN
      END FUNCTION gasdev
