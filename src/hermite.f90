! ==============================================================
   PROGRAM hermite
!
! HENDOS - HErmite N-body Dynamics Of Systems
! (Sheffield N-body code for stellar dynamics)
!
! 4th-order Hermite integrator (Makino & Aarseth 1992, PASJ, 44, 141)
! individual block timesteps in powers of two (Makino 1991, PASJ, 43, 859)
! timestep = min(Aarseth criterion, pairwise free-fall/flyby criterion)
! predictor raised in order with each star's a2, a3 from its last step
! optional P(EC)^n iteration of the corrector
! optional Barnes-Hut tree with octupole moments           (module nbtree)
! optional Galactic tide (Hill's approximation, rotating frame)
! optional Ahmad-Cohen neighbour scheme for direct sums   (module nbac)
! unperturbed binaries treated as one body with an exact Kepler orbit,
!   isolated few-body encounters integrated separately, stable triples
!   as nested Kepler orbits; all opened when perturbed (no
!   regularisation)                               (modules nbbin, nbfew)
! stellar collisions: the stars are merged (or the run stops)
! optional binary restart files for long runs       (module nbrestart)
! all the expensive loops are parallelised with OpenMP (compile with
!   -fopenmp; the number of threads is set with OMP_NUM_THREADS)
!
! This file is the main program only: set-up, then the loop over output
! intervals and blocks (predict, forces, correct, new timesteps).
! Everything else is in the modules - see src/README.md for what each
! routine does and where it is.
!
! block times are held as 8-byte integers in units of tick=dtmax/2**kmax,
! so they are exact. All particles are synchronised every dtmax (= tout),
! at which point the integer clock is reset to zero.
!
! usage:   ./hermite params.nml  [ | tee data/<runid>.log ]
! all run parameters are in the parameter file (default params.nml)
!
! DIRECTORIES: source in HENDOS/src, built and run in HENDOS/run (which
! holds the Makefile, executables and parameter files); snapshot and
! analysis files are written to the data directory (parameter datadir,
! normally HENDOS/run/data)
!
! OUTPUT (all in Msun, pc, km/s, Myr) - file names start with runid
!   screen     standard or verbose (parameter screen), every tout
!   snapshots  <runid>_snap_<time>Myr.dat every tsnap   (module nbout)
!   analysis   <runid>_ana_<name>.dat every tana      (module nbana)
!   mergers    <runid>_mergers.dat, one line per merger  (module nbbin)
!   restart    <runid>_restart.bin every trestart    (module nbrestart)
!
! INITIAL CONDITIONS: a snapshot-format file (NBODY-SNAPSHOT), e.g.
! from makeics or a previous run (restarts from its time), or a
! legacy file with lines  m x y z vx vy vz  (e.g. from bigfractals)
!
! PARTICLE TYPES: every particle has a type (star, white_dwarf,
! neutron_star, black_hole, planet - the table is in module nbvars),
! read from the initial conditions and written to the snapshots.
! Particles in files without types are stars.
!
   USE nbvars
   USE nbsetup
   USE nbtree
   USE nbac
   USE nbforce
   USE nbstep
   USE nbbin
   USE nbout
   USE nbrestart
   USE nbana
   IMPLICIT NONE
   INTEGER :: i,ia,it,iout,iout0,nout,nsnap,nana,nrest,nsince,ic,k,nper
   INTEGER, ALLOCATABLE :: iscr(:),cact(:),coff(:)
   INTEGER(KIND=8) :: tblock,tend8,nblock,nforce
   DOUBLE PRECISION :: dt,dtp,time
   DOUBLE PRECISION :: a1(1:3),a2(1:3),dt2,dt3,dt4,dt5
!
! ============================================================
! parameter file name from the command line
   IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
     CALL GET_COMMAND_ARGUMENT(1,pfilename)
   ELSE
     pfilename='params.nml'
   END IF
!
! units, then run parameters (physical units)
   CALL setunits
   CALL readparams(pfilename)
   CALL checkdatadir
   eps2=eps*eps                       ! length unit is pc already
!
! ============================================================
! initial conditions (sets n, time0 and allocates arrays), or the
! complete state of an earlier run from its restart file
   iout0=0
   nblock=0
   nforce=0
   IF (restart) THEN
     CALL readrestart(time,iout0,nblock,nforce)
   ELSE
     CALL readic
     IF (useac) CALL acalloc
   END IF
!
! number of output intervals, snapshot, analysis and restart frequency
   nout=NINT((tend - time0)/tout)
   IF (nout < 1) THEN
     WRITE(6,'(a,f0.4,a)') ' ERROR: tend must be at least one tout after '// &
&      'the start time of the initial conditions (', time0, ' Myr)'
     STOP
   END IF
   IF (restart .AND. iout0 >= nout) THEN
     WRITE(6,'(a,f0.4,a)') ' ERROR: the restart file is already at ', time, &
&      ' Myr - increase tend to carry the run on'
     STOP
   END IF
   nsnap=0
   IF (tsnap > 0.d0) nsnap=MAX(1,NINT(tsnap/tout))
   nana=0
   IF (tana > 0.d0) nana=MAX(1,NINT(tana/tout))
   nrest=0
   IF (trestart > 0.d0) nrest=MAX(1,NINT(trestart/tout))
!
! a new run must not overwrite the files of an earlier run
   IF (time0 == 0.d0 .AND. .NOT. restart) CALL checkrunid(nsnap,nana)
   IF (usetree) CALL treealloc
!
! work space for making the list of active particles in parallel
   ALLOCATE(iscr(1:nalloc),cact(1:nthreads),coff(1:nthreads))
!
! block timestep set-up
   dtmax=tout/tunit
   tick=dtmax/2.d0**kmax
   tend8=ISHFT(1_8,kmax)              ! end of one output interval in ticks
!
   IF (restart) THEN
     IF (nana > 0) CALL ana_init
   ELSE
! initial forces and timestep levels, then start-of-run information
     nclamp=0
     CALL firststep
     IF (rcoll > 0.d0) CALL collcheck(0_8,time0/tunit)
     IF (closebinaries) CALL binstart(time0/tunit)
     CALL startinfo(nout)
!
     time=time0
     IF (nsnap > 0) THEN
! (not when restarting from that same snapshot)
       IF (TRIM(snapname(time)) /= TRIM(icfile)) CALL snapshot(time,nforce)
     END IF
     IF (nana > 0) THEN
       CALL ana_init
       IF (time0 == 0.d0) CALL ana_do(time)
     END IF
   END IF
   nsince=0
!
! =====================================================================
! main loop - over output intervals, then over blocks
   DO iout=iout0 + 1,nout
     tbase=(time0 + DBLE(iout - 1)*tout)/tunit
! the tree is always rebuilt at the start of an output interval, so
! that a restarted run does exactly what the original would have done
     nsince=ntreebuild
!
     DO
!
! find next block time and the active particles
       tblock=HUGE(1_8)
!$OMP PARALLEL DO IF(n >= nompj) REDUCTION(MIN:tblock) SCHEDULE(STATIC)
       DO i=1,n
         tblock=MIN(tblock,tk(i) + dtk(i))
       END DO
!$OMP END PARALLEL DO
! the list of active particles, in order of index. The particles are cut
! into one piece per thread; each thread lists the active ones of its
! piece in the scratch list iscr, then the pieces are joined up in ilist
       nper=(n + nthreads - 1)/nthreads
!$OMP PARALLEL DO IF(n >= nompj) PRIVATE(ic,i,k) SCHEDULE(STATIC,1)
       DO ic=1,nthreads
         k=(ic - 1)*nper
         DO i=(ic - 1)*nper + 1,MIN(n,ic*nper)
           IF (tk(i) + dtk(i) == tblock) THEN
             k=k + 1
             iscr(k)=i
           END IF
         END DO
         cact(ic)=k - (ic - 1)*nper
       END DO
!$OMP END PARALLEL DO
       coff(1)=0
       DO ic=2,nthreads
         coff(ic)=coff(ic-1) + cact(ic-1)
       END DO
       nact=coff(nthreads) + cact(nthreads)
!$OMP PARALLEL DO IF(nact >= nompmin) PRIVATE(ic,k) SCHEDULE(STATIC,1)
       DO ic=1,nthreads
         DO k=1,cact(ic)
           ilist(coff(ic) + k)=iscr((ic - 1)*nper + k)
         END DO
       END DO
!$OMP END PARALLEL DO
!
! predictor - predict ALL particles to tblock for the force calculation.
! Uses F0, F1 and also f2e, f3e (2nd and 3rd derivatives from each
! star's previous step). The extra terms are essential with the tide:
! the Coriolis force depends on velocity, and a 2nd-order velocity
! prediction would make the whole scheme only 3rd order
!$OMP PARALLEL DO IF(n >= nompj) PRIVATE(i,dtp) SCHEDULE(STATIC)
       DO i=1,n
         dtp=DBLE(tblock - tk(i))*tick
         rp(1:3,i)=r(1:3,i) + v(1:3,i)*dtp + f0(1:3,i)*dtp**2/2.d0 + &
&                  f1(1:3,i)*dtp**3/6.d0 + f2e(1:3,i)*dtp**4/24.d0 + &
&                  f3e(1:3,i)*dtp**5/120.d0
         vp(1:3,i)=v(1:3,i) + f0(1:3,i)*dtp + f1(1:3,i)*dtp**2/2.d0 + &
&                  f2e(1:3,i)*dtp**3/6.d0 + f3e(1:3,i)*dtp**4/24.d0
       END DO
!$OMP END PARALLEL DO
!
! tree: rebuild it every ntreebuild block steps (or after a binary was
! closed/opened); on the steps in between just restock it - recompute
! the moments from the new predicted positions. nsince counts the block
! steps since the last rebuild
       IF (usetree) THEN
         IF (.NOT. treevalid .OR. nsince >= ntreebuild) THEN
           CALL treebuild
           nsince=0
         ELSE
           CALL treemoments
         END IF
         nsince=nsince + 1
       END IF
!
! base for the corrector of the active stars: the 3rd-order prediction
! from F0 and F1 only (the corrector adds the f2, f3 terms itself)
!$OMP PARALLEL DO IF(nact >= nompmin) PRIVATE(ia,i,dt) SCHEDULE(STATIC)
       DO ia=1,nact
         i=ilist(ia)
         dt=DBLE(dtk(i))*tick
         rc(1:3,i)=r(1:3,i) + v(1:3,i)*dt + f0(1:3,i)*dt**2/2.d0 + &
&                  f1(1:3,i)*dt**3/6.d0
         vc(1:3,i)=v(1:3,i) + f0(1:3,i)*dt + f1(1:3,i)*dt**2/2.d0
       END DO
!$OMP END PARALLEL DO
!
! evaluate and correct the active particles (repeated niter times)
       tcur8=tblock
       DO it=1,niter
         aclast=(it == niter)
!
! forces on active particles from all (predicted) particles
         CALL forces
!
! corrector terms f2 and f3 (2nd and 3rd derivatives at start of step)
!$OMP PARALLEL DO IF(nact >= nompmin) SCHEDULE(STATIC) &
!$OMP& PRIVATE(ia,i,dt,dt2,dt3,dt4,dt5,a1,a2)
         DO ia=1,nact
           i=ilist(ia)
           dt=DBLE(dtk(i))*tick
           dt2=dt*dt
           dt3=dt2*dt
           dt4=dt3*dt
           dt5=dt4*dt
           a1(1:3)=f0(1:3,i) - f0p(1:3,i)
           a2(1:3)=f1(1:3,i) + f1p(1:3,i)
           f2(1:3,i)=2.d0*(-3.d0*a1(1:3) - (f1(1:3,i) + a2(1:3))*dt)/dt2
           f3(1:3,i)=6.d0*(2.d0*a1(1:3) + a2(1:3)*dt)/dt3
!
! add correction (the corrected values are used by any further iteration)
           rp(1:3,i)=rc(1:3,i) + f2(1:3,i)*dt4/24.d0 + f3(1:3,i)*dt5/120.d0
           vp(1:3,i)=vc(1:3,i) + f2(1:3,i)*dt3/6.d0  + f3(1:3,i)*dt4/24.d0
! neighbour scheme: on a regular step, the correction from the regular
! force over the whole regular step as well
           IF (useac) THEN
             IF (lreg(i)) THEN
               rp(1:3,i)=rp(1:3,i) + dreg(1:3,i)
               vp(1:3,i)=vp(1:3,i) + dvreg(1:3,i)
             END IF
           END IF
         END DO
!$OMP END PARALLEL DO
       END DO
!
! accept step and get new levels
!$OMP PARALLEL DO IF(nact >= nompmin) PRIVATE(ia,i,dt) SCHEDULE(STATIC)
       DO ia=1,nact
         i=ilist(ia)
         r(1:3,i)=rp(1:3,i)
         v(1:3,i)=vp(1:3,i)
         f0(1:3,i)=f0p(1:3,i)
         f1(1:3,i)=f1p(1:3,i)
         dt=DBLE(dtk(i))*tick
         f2e(1:3,i)=f2(1:3,i) + f3(1:3,i)*dt
         f3e(1:3,i)=f3(1:3,i)
         tk(i)=tblock
         IF (useac) CALL acaccept(i)
         CALL newlevel(i)
         IF (useac) CALL aclevel(i)
       END DO
!$OMP END PARALLEL DO
       nblock=nblock + 1
       nforce=nforce + nact*niter
!
! stellar collisions (merge the stars or stop the run), then open,
! close and freeze binaries and few-body groups (serial, rare)
       IF (rcoll > 0.d0) CALL collcheck(tblock,tbase + DBLE(tblock)*tick)
       IF (closebinaries) CALL bincheck(tblock,tbase + DBLE(tblock)*tick)
!
! end of output interval - everyone is synchronised here
       IF (tblock == tend8) EXIT
     END DO
!
! reset integer clock (regular times too: they are all at the end of
! the interval, or keep their place relative to it)
     tk=0
     IF (useac) treg(1:n)=treg(1:n) - tend8
     time=time0 + DBLE(iout)*tout     ! Myr
!
! output ----------------------------------------------
! screen first: it also updates the diagnostics used by the others
     CALL statusout(iout,nout,time,nblock,nforce)
     IF (nsnap > 0) THEN
       IF (MOD(iout,nsnap) == 0) CALL snapshot(time,nforce)
     END IF
     IF (nana > 0) THEN
       IF (MOD(iout,nana) == 0) CALL ana_do(time)
     END IF
     IF (nrest > 0) THEN
       IF (MOD(iout,nrest) == 0 .OR. iout == nout) &
&        CALL writerestart(time,iout,nblock,nforce)
     END IF
!
   END DO
!
   IF (nana > 0) CALL ana_end
   CALL endinfo(time,nforce)
!
   END PROGRAM hermite
