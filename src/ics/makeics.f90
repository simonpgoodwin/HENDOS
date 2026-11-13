   PROGRAM makeics
!
! HENDOS - HErmite N-body Dynamics Of Systems
! THE INITIAL CONDITIONS GENERATOR
!
! usage:   ./makeics ics.nml
! With nensemble > 1 in &ics it makes an ENSEMBLE: nensemble sets of
! initial conditions that differ only in the random number seed (seed,
! seed+1, seed+2, ...), written to files with the member number in
! their names: outfile = 'data/plum1k.dat' gives data/plum1k_e001.dat,
! data/plum1k_e002.dat, ... Member k is exactly what a single run with
! seed + k - 1 would make. hermite runs the whole ensemble if it is
! given the same nensemble (see run/PARAMS_README.md).
! The input is in the parameter file (default ics.nml), described in
! run/PARAMS_README.md. How to change this code is described below
! and, at more length, in src/ics/ICS_README.md.
!
! =====================================================================
! WHAT IT DOES
! =====================================================================
! The main program (below) is seven steps, one CALL each:
!
!   STEP 1  STRUCTURE   positions and velocities of n stars, in model
!                       units, every mass unity       (makestructure)
!   STEP 2  MASSES      a mass and a type for each star  (makemasses)
!   STEP 3  SIZE        positions scaled to pc           (in the main program)
!   STEP 4  FRAME       moved to the centre of mass      (tocom)
!   STEP 5  VELOCITIES  put into km/s: scaled to the virial ratio
!                       (virialise), or for a Larson fractal set to
!                       the size-linewidth relation      (larsonise)
!   STEP 6  YOUR OWN    anything you like, last of all   (userfinal)
!   STEP 7  OUTPUT      summary and file                 (report, writeics)
!
! =====================================================================
! THE STARS
! =====================================================================
! Everything about the stars is in four arrays, used by every step:
!     m(i)         mass of star i
!     r(1:3,i)     position x, y, z
!     v(1:3,i)     velocity vx, vy, vz
!     ptype(i)     particle type: 'star', 'white_dwarf', 'neutron_star',
!                  'black_hole' or 'planet'
! for i = 1 to n. Their units change as the steps go by:
!                    m          r                  v
!   after STEP 1     1          model units        arbitrary
!   after STEP 2     Msun       model units        arbitrary
!   after STEP 3     Msun       pc                 arbitrary
!   after STEP 5     Msun       pc                 km/s
! 'model units' are whatever the structure routine works in (Plummer
! radius = 1 for a Plummer sphere, total radius = 1 for a fractal);
! STEP 1 also sets rscale, the number that turns them into pc.
!
! =====================================================================
! HOW TO CHANGE THINGS
! =====================================================================
! Without touching this file - edit src/ics/usersubs.f90, then make:
!   your own structure      set distribution = 'user' in ics.nml and
!                           write it in  userstructure
!   your own masses         set imf = 'user' and write them in  usermasses
!   anything else at all    write it in  userfinal, which is called
!                           after everything else and has the last word
!                           (e.g. change types, move or add to
!                           velocities, swap stars round)
!
! In this file - every place is marked with  ! >>> CHANGE HERE
!   a new named structure   add a CASE in  makestructure
!   a new named mass function   add a CASE in  makemasses
!   a new fractal velocity field   add a CASE in  makestructure
!   particle types          set ptype in  makemasses
!   what radius means       rscale in  makestructure
!   how velocities are set  STEP 5 in the main program, and  virialise
!   leave out or reorder a step   the main program: each step is one
!                           CALL, and the comment on each routine says
!                           what it needs and what it changes
!   a new input parameter   three places, each marked  ! >>> NEW PARAMETER:
!                           declare it, add it to its namelist and give
!                           it a standard value in  readinput, and put
!                           it in run/ics.nml
! A new routine in a new file must be added to ICSRC in run/Makefile.
!
! Routines called (all in src/ics):
!   makeplummer                          plummersubs.f90
!   makefractal, larsonfractal           fractalsubs.f90
!   maschberger                          masssubs.f90
!   userstructure, usermasses, userfinal usersubs.f90
!   ran3, gasdev                         icrandom.f90
! Random numbers: x=ran3(idum) is uniform between 0 and 1, x=gasdev(idum)
! is Gaussian. Always pass idum, and never change it yourself.
! =====================================================================
!
   IMPLICIT NONE
!
! ---- input parameters (read from the parameter file by readinput) ----
! >>> NEW PARAMETER: declare it here
! &ics
   CHARACTER(len=256) :: outfile       ! file to write
   INTEGER :: n                        ! number of stars
   INTEGER :: seed                     ! random number seed
   INTEGER :: nensemble                ! number of sets to make (1 = just one)
   CHARACTER(len=16) :: distribution   ! 'plummer', 'fractal' or 'user'
   DOUBLE PRECISION :: radius          ! pc: half-mass (plummer), total (fractal)
   DOUBLE PRECISION :: qvir            ! virial ratio
   LOGICAL :: setvirial                ! scale the velocities to qvir?
! &masses
   CHARACTER(len=16) :: imf            ! 'maschberger', 'equal' or 'user'
   DOUBLE PRECISION :: mlo,mup         ! mass limits (Msun)
   DOUBLE PRECISION :: alpha,beta,mu   ! Maschberger IMF parameters
   DOUBLE PRECISION :: mequal          ! mass of every star if imf='equal'
! &fractal
   DOUBLE PRECISION :: fdim            ! fractal dimension
   INTEGER :: ndiv                     ! boxes per side at each division
   CHARACTER(len=16) :: velocities     ! 'standard','inherit','random','larson'
   DOUBLE PRECISION :: coresize        ! core size (pc), 0 = none
   DOUBLE PRECISION :: vlarson         ! Larson coefficient (km/s at 1 pc)
!
! ---- the stars ----
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: m        ! masses
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: r,v    ! positions, velocities
   CHARACTER(len=12), ALLOCATABLE, DIMENSION(:) :: ptype   ! particle types
!
! ---- other ----
   INTEGER :: idum                     ! random number state for ran3, gasdev
   DOUBLE PRECISION :: rscale          ! pc per model length (set in STEP 1)
   LOGICAL :: larson                   ! is this a Larson fractal?
   DOUBLE PRECISION :: mtot            ! total mass (Msun)
   DOUBLE PRECISION :: ekin,epot       ! energies (set by energies)
   DOUBLE PRECISION :: qmade           ! virial ratio of what was made
   DOUBLE PRECISION :: gcons           ! G in pc (km/s)^2 / Msun
   CHARACTER(len=256) :: pfile         ! parameter file name
   INTEGER :: imem                     ! member of the ensemble being made
   INTEGER :: seedmem                  ! its random number seed
   CHARACTER(len=256) :: outname       ! its file name
   INTEGER :: izero                    ! 0, to restart gasdev (see below)
   REAL :: gdum
   REAL, EXTERNAL :: gasdev
!
   gcons=4.300917270619619d-3
!
! parameter file name from the command line, then read it
   IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
     CALL GET_COMMAND_ARGUMENT(1,pfile)
   ELSE
     pfile='ics.nml'
   END IF
   CALL readinput
!
   ALLOCATE(m(1:n))
   ALLOCATE(r(1:3,1:n))
   ALLOCATE(v(1:3,1:n))
   ALLOCATE(ptype(1:n))
!
! =====================================================================
! one pass through the seven steps for each member of the ensemble
! (just one pass if nensemble = 1). Only the seed and the file name
! differ between members.
! =====================================================================
   DO imem=1,nensemble
   seedmem=seed + imem - 1
   outname=outfile
   IF (nensemble > 1) outname=ensname(outfile,imem)
!
! start the random numbers (ran3 needs a negative number to start)
   idum=-ABS(seedmem)
   IF (idum == 0) idum=-1
! (and gasdev must forget the number it kept from the member before,
! so that member k is exactly what a single run with its seed makes)
   izero=0
   gdum=gasdev(izero)
!
! =====================================================================
! STEP 1 - STRUCTURE
! positions and velocities in model units, every mass set to unity.
! Also sets rscale and larson.
! =====================================================================
   CALL makestructure
!
! =====================================================================
! STEP 2 - MASSES
! each star's mass (Msun) and type
! =====================================================================
   CALL makemasses
   mtot=SUM(m)
!
! =====================================================================
! STEP 3 - SIZE
! positions from model units to pc
! =====================================================================
   r=r*rscale
!
! =====================================================================
! STEP 4 - FRAME
! centre of mass at the origin and at rest
! =====================================================================
   CALL tocom
!
! =====================================================================
! STEP 5 - VELOCITIES
! from arbitrary units to km/s
! >>> CHANGE HERE to set the size of the velocities in another way
! =====================================================================
   IF (larson) THEN
! Larson fractal: the observed size-linewidth relation, qvir not used
     CALL larsonise
   ELSE IF (setvirial) THEN
! everything else: one factor, to give exactly the virial ratio qvir
     CALL virialise
   END IF
! (with setvirial = .false. the velocities are left exactly as STEP 1
!  made them, and are taken to be in km/s already)
!
! =====================================================================
! STEP 6 - YOUR OWN CHANGES
! userfinal is in usersubs.f90. As supplied it does nothing.
! Everything is in its final units here: Msun, pc, km/s.
! =====================================================================
   CALL userfinal(n,m,r,v,ptype,idum)
   mtot=SUM(m)
!
! =====================================================================
! STEP 7 - OUTPUT
! =====================================================================
   CALL report
   CALL writeics
   END DO
!
   IF (nensemble > 1) WRITE(6,'(a,i0,a,a,a,a)') '  ensemble of ', nensemble, &
&    ' written: ', TRIM(ensname(outfile,1)), ' to ', TRIM(ensname(outfile,nensemble))
!
   CONTAINS
! =====================================================================
   FUNCTION ensname(name,k) RESULT(str)
!
! the name of member k of an ensemble: _e and the member number (3
! digits, more if needed) put before the file extension,
!   data/plum1k.dat -> data/plum1k_e007.dat
! hermite uses the same rule (ensname in src/nbsetup.f90) - keep the
! two the same.
!
   CHARACTER(len=*), INTENT(IN) :: name
   INTEGER, INTENT(IN) :: k
   CHARACTER(len=256) :: str
   CHARACTER(len=16) :: num
   INTEGER :: idot,islash
!
   IF (k < 1000) THEN
     WRITE(num,'(a,i3.3)') '_e', k
   ELSE
     WRITE(num,'(a,i0)') '_e', k
   END IF
   idot=INDEX(TRIM(name),'.',BACK=.TRUE.)
   islash=INDEX(TRIM(name),'/',BACK=.TRUE.)
   IF (idot > islash + 1) THEN
     str=name(1:idot - 1)//TRIM(num)//name(idot:LEN_TRIM(name))
   ELSE
     str=TRIM(name)//TRIM(num)
   END IF
!
   END FUNCTION ensname
! =====================================================================
! =====================================================================
   SUBROUTINE makestructure
!
! STEP 1. Positions r and velocities v of the n stars.
!   needs:    the input parameters, idum
!   changes:  r (model units; r*rscale is in pc)
!             v (arbitrary units, STEP 5 sets their size)
!             m (all set to 1)
!             rscale (the number that takes r to pc)
!             larson (true for a Larson fractal)
!
! >>> CHANGE HERE to add a new structure:
!     1. write a routine that returns positions and velocities (in its
!        own file, like makeplummer in plummersubs.f90)
!     2. add a CASE below that calls it and sets rscale
!     3. add its name to the error message in CASE DEFAULT
!     4. add the file to ICSRC in run/Makefile
!     Nothing else needs to know: STEPS 2 to 7 work for any structure.
!
   REAL, ALLOCATABLE, DIMENSION(:,:) :: rs,vs
   DOUBLE PRECISION :: mcut,rhfac
!
   larson=.FALSE.
!
   SELECT CASE (TRIM(distribution))
!
! ---- Plummer sphere -------------------------------------------------
   CASE ('plummer')
! made with Plummer radius = 1, cut at 10 Plummer radii
     CALL makeplummer(n,r,v,idum)
! radius is the HALF-MASS radius. For a Plummer sphere cut at 10 it is
! rhfac = 1.288 Plummer radii (1.305 uncut): mcut is the mass fraction
! inside the cut, and rhfac the radius holding half of that
! >>> CHANGE HERE if radius should mean something else, e.g. for the
!     Plummer radius itself use  rscale=radius
     mcut=1.d3/(1.d0 + 1.d2)**1.5d0
     rhfac=1.d0/SQRT((0.5d0*mcut)**(-2.d0/3.d0) - 1.d0)
     rscale=radius/rhfac
!
! ---- box fractal ----------------------------------------------------
   CASE ('fractal')
! the fractal routines work in single precision, in a sphere of radius
! 1, so the core size is passed as a fraction of the radius
     ALLOCATE(rs(1:3,1:n))
     ALLOCATE(vs(1:3,1:n))
! >>> CHANGE HERE to add a new fractal velocity field: add a CASE
     SELECT CASE (TRIM(velocities))
     CASE ('standard')
       CALL makefractal(n,rs,vs,REAL(fdim),ndiv,idum,1.0,1.0,REAL(coresize/radius),'s')
     CASE ('inherit')
       CALL makefractal(n,rs,vs,REAL(fdim),ndiv,idum,1.0,1.0,REAL(coresize/radius),'f')
     CASE ('random')
       CALL makefractal(n,rs,vs,REAL(fdim),ndiv,idum,1.0,1.0,REAL(coresize/radius),'u')
     CASE ('larson')
       CALL larsonfractal(n,rs,vs,REAL(fdim),ndiv,idum,1.0,1.0,REAL(coresize/radius))
       larson=.TRUE.
     CASE DEFAULT
       CALL bad('velocities must be ''standard'', ''inherit'', ''random'' or ''larson''')
     END SELECT
     r=DBLE(rs)
     v=DBLE(vs)
     DEALLOCATE(rs,vs)
! radius is the TOTAL radius
     rscale=radius
!
! ---- your own (usersubs.f90) ----------------------------------------
   CASE ('user')
     CALL userstructure(n,r,v,idum)
! userstructure returns positions with the length radius refers to = 1
     rscale=radius
!
   CASE DEFAULT
     CALL bad('distribution must be ''plummer'', ''fractal'' or ''user''')
   END SELECT
!
! every structure comes back with all masses unity
   m=1.d0
!
   END SUBROUTINE makestructure
! =====================================================================
   SUBROUTINE makemasses
!
! STEP 2. The mass and type of each star.
!   needs:    n, idum, the &masses parameters (r is there too, in
!             model units, for masses that depend on position)
!   changes:  m (Msun), ptype
! All stars are single. Binaries and multiples are not made yet.
!
! >>> CHANGE HERE to add a new mass function:
!     1. write a routine that returns n masses (in its own file, like
!        maschberger in masssubs.f90)
!     2. add a CASE below that calls it
!     3. add its name to the error message in CASE DEFAULT
!     4. add the file to ICSRC in run/Makefile
!
   SELECT CASE (TRIM(imf))
!
! ---- Maschberger (2013) IMF -----------------------------------------
   CASE ('maschberger')
     CALL maschberger(n,m,idum,alpha,beta,mu,mlo,mup)
!
! ---- all the same mass ----------------------------------------------
   CASE ('equal')
     m=mequal
!
! ---- your own (usersubs.f90) ----------------------------------------
   CASE ('user')
     CALL usermasses(n,m,r,idum)
!
   CASE DEFAULT
     CALL bad('imf must be ''maschberger'', ''equal'' or ''user''')
   END SELECT
!
! particle types
! >>> CHANGE HERE to give stars other types, e.g.
!        WHERE (m > 20.d0) ptype='black_hole'
!     The names must be ones the N-body code knows (typename in
!     src/nbvars.f90).
   ptype='star'
!
   END SUBROUTINE makemasses
! =====================================================================
   SUBROUTINE tocom
!
! STEP 4. Move to the centre of mass and the centre of velocity.
!   needs:    m, r, v, mtot
!   changes:  r, v
!
   INTEGER :: i
   DOUBLE PRECISION :: com(1:3),cov(1:3)
!
   com=0.d0
   cov=0.d0
   DO i=1,n
     com(1:3)=com(1:3) + m(i)*r(1:3,i)
     cov(1:3)=cov(1:3) + m(i)*v(1:3,i)
   END DO
   com=com/mtot
   cov=cov/mtot
   DO i=1,n
     r(1:3,i)=r(1:3,i) - com(1:3)
     v(1:3,i)=v(1:3,i) - cov(1:3)
   END DO
!
   END SUBROUTINE tocom
! =====================================================================
   SUBROUTINE virialise
!
! STEP 5, the usual case. Multiply all velocities by one factor so that
! kinetic energy / potential energy = qvir exactly, for the masses
! (Msun) and positions (pc) the stars now have.
!   needs:    m (Msun), r (pc), v (any units), qvir
!   changes:  v (km/s)
!
! >>> CHANGE HERE to scale the velocities to something else: fac is
!     the one number every velocity is multiplied by
!
   DOUBLE PRECISION :: fac
!
   CALL energies
   IF (ekin <= 0.d0) THEN
     IF (qvir > 0.d0) CALL bad('all velocities are zero, cannot scale them to qvir')
   ELSE
     fac=SQRT(qvir*epot/ekin)
     v=v*fac
   END IF
!
   END SUBROUTINE virialise
! =====================================================================
   SUBROUTINE larsonise
!
! STEP 5, Larson fractals only. Put the velocities on the observed
! size-linewidth relation, in km/s. qvir is NOT used - the virial ratio
! is whatever these velocities give for the masses and size.
!   needs:    v (as made by larsonfractal), rscale, vlarson
!   changes:  v (km/s)
!
! larsonfractal was called with radius = 1 and vscale = 1, so the mean
! square difference of one velocity component between two stars is
! equal to their separation in units of the radius. Multiplying by
! vlarson*sqrt(radius in pc) makes it
!      <(v_i - v_j)^2> = vlarson^2 * (separation in pc)     (km/s)^2
! per component. A region filled uniformly with stars to radius R pc
! then has 1D velocity dispersion 0.72*vlarson*sqrt(R) km/s, so the
! standard vlarson = 1 is the relation of Solomon et al. (1987, ApJ,
! 319, 730), sigma = 0.72 (R/pc)^0.5 km/s. See larson.tex.
!
   v=v*vlarson*SQRT(rscale)
!
   END SUBROUTINE larsonise
! =====================================================================
   SUBROUTINE energies
!
! kinetic energy ekin and potential energy epot (positive) of the stars
! as they are now. In Msun (km/s)^2 if m is in Msun, r in pc and v in
! km/s. Does not change the stars.
!
   INTEGER :: i,j
   DOUBLE PRECISION :: epoti,dr(1:3)
!
   ekin=0.d0
   DO i=1,n
     ekin=ekin + 0.5d0*m(i)*(v(1,i)**2 + v(2,i)**2 + v(3,i)**2)
   END DO
   epot=0.d0
   DO i=1,n-1
     epoti=0.d0
     DO j=i+1,n
       dr(1:3)=r(1:3,i) - r(1:3,j)
       epoti=epoti + m(j)/SQRT(dr(1)**2 + dr(2)**2 + dr(3)**2)
     END DO
     epot=epot + epoti*m(i)
   END DO
   epot=gcons*epot
!
   END SUBROUTINE energies
! =====================================================================
   SUBROUTINE report
!
! STEP 7. Summary on the screen of what was actually made (measured
! from the stars, after any changes in userfinal). Sets qmade.
! Does not change the stars.
!
   INTEGER :: i,it
   DOUBLE PRECISION :: rlo,rhi,rmid,rmax,minside,sig,com(1:3)
!
! virial ratio of what was made
   CALL energies
   qmade=0.d0
   IF (epot > 0.d0) qmade=ekin/epot
! centre of mass (userfinal may have moved it)
   com=0.d0
   DO i=1,n
     com(1:3)=com(1:3) + m(i)*r(1:3,i)
   END DO
   com=com/mtot
! half-mass radius about the centre of mass, by bisection
   rmax=0.d0
   DO i=1,n
     rmax=MAX(rmax,SQRT(SUM((r(1:3,i) - com(1:3))**2)))
   END DO
   rlo=0.d0
   rhi=rmax
   DO it=1,60
     rmid=0.5d0*(rlo + rhi)
     minside=0.d0
     DO i=1,n
       IF (SUM((r(1:3,i) - com(1:3))**2) <= rmid*rmid) minside=minside + m(i)
     END DO
     IF (minside >= 0.5d0*mtot) THEN
       rhi=rmid
     ELSE
       rlo=rmid
     END IF
   END DO
! mass-weighted 3D velocity dispersion
   sig=SQRT(2.d0*ekin/mtot)
!
   WRITE(6,'(a)') ' ------------------------------------------------------------'
   WRITE(6,'(a)') '  HENDOS initial conditions'
   WRITE(6,'(a,a)')         '  distribution            ', TRIM(distribution)
   IF (TRIM(distribution) == 'fractal') &
&  WRITE(6,'(a,a)')         '  fractal velocities      ', TRIM(velocities)
   WRITE(6,'(a,i12)')       '  stars                   ', n
   WRITE(6,'(a,a)')         '  masses                  ', TRIM(imf)
   WRITE(6,'(a,f12.3,a)')   '  total mass              ', mtot, ' Msun'
   WRITE(6,'(a,3f12.4,a)')  '  mass min, mean, max     ', MINVAL(m), mtot/DBLE(n), &
&                           MAXVAL(m), ' Msun'
   WRITE(6,'(a,f12.4,a)')   '  radius asked for        ', radius, ' pc'
   WRITE(6,'(a,f12.4,a)')   '  half-mass radius made   ', rhi, ' pc'
   WRITE(6,'(a,f12.4,a)')   '  furthest star           ', rmax, ' pc'
   IF (larson) THEN
     WRITE(6,'(a,f12.4,a)') '  Larson coefficient      ', vlarson, ' km/s at 1 pc'
     WRITE(6,'(a,f12.4,a)') '  virial ratio made       ', qmade, &
&                           ' (from the Larson velocities; qvir is not used)'
   ELSE IF (setvirial) THEN
     WRITE(6,'(a,f12.4,a,f8.4,a)') '  virial ratio made       ', qmade, &
&                           ' (asked for', qvir, ')'
   ELSE
     WRITE(6,'(a,f12.4,a)') '  virial ratio made       ', qmade, &
&                           ' (velocities not scaled: setvirial is off)'
   END IF
   WRITE(6,'(a,f12.4,a)')   '  3D velocity dispersion  ', sig, ' km/s'
   IF (nensemble > 1) WRITE(6,'(a,i12,a,i0,a,i0,a)') '  ensemble member         ', &
&    imem, ' of ', nensemble, '  (seed ', seedmem, ')'
   WRITE(6,'(a,a)')         '  written to              ', TRIM(outname)
   WRITE(6,'(a)') ' ------------------------------------------------------------'
!
   END SUBROUTINE report
! =====================================================================
   SUBROUTINE writeics
!
! STEP 7. Write the file in the snapshot format (NBODY-SNAPSHOT 2) read
! by hermite: a header of '# key = value' lines, then one line per star
!     id m x y z vx vy vz type
! Must come after report (which sets qmade). Does not change the stars.
!
! >>> CHANGE HERE to record something else in the header: add a line
!     like the others. hermite ignores header lines it does not know.
!     Do not change the columns - the N-body code reads them.
!
   INTEGER :: i,ios,dv(1:8)
   CHARACTER(len=19) :: date
!
   CALL DATE_AND_TIME(values=dv)
   WRITE(date,'(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') dv(1), '-', dv(2), &
&    '-', dv(3), ' ', dv(5), ':', dv(6), ':', dv(7)
   OPEN(3,file=TRIM(outname),status='replace',iostat=ios)
   IF (ios /= 0) CALL bad('cannot write '//TRIM(outname)//' - does the directory exist?')
   WRITE(3,'(a)') '# format                        = NBODY-SNAPSHOT 2'
   WRITE(3,'(a)') '# creator                       = HENDOS (makeics.f90)'
   WRITE(3,'(a,a)') '# date                          = ', date
   WRITE(3,'(a,a)') '# parameter_file                = ', TRIM(pfile)
   WRITE(3,'(a,a)') '# distribution                  = ', TRIM(distribution)
   WRITE(3,'(a,es24.16)') '# radius_pc                     = ', radius
   IF (TRIM(distribution) == 'plummer') THEN
     WRITE(3,'(a)') '# radius_is                     = half-mass radius'
     WRITE(3,'(a,es24.16)') '# plummer_radius_pc             = ', rscale
     WRITE(3,'(a)') '# plummer_cut                   = 10 Plummer radii'
   ELSE IF (TRIM(distribution) == 'fractal') THEN
     WRITE(3,'(a)') '# radius_is                     = total radius'
     WRITE(3,'(a,es24.16)') '# fractal_dimension             = ', fdim
     WRITE(3,'(a,i0)') '# fractal_ndiv                  = ', ndiv
     WRITE(3,'(a,a)') '# fractal_velocities            = ', TRIM(velocities)
     WRITE(3,'(a,es24.16)') '# fractal_coresize_pc           = ', coresize
     IF (larson) WRITE(3,'(a,es24.16)') '# larson_coefficient_km_s       = ', vlarson
   END IF
   WRITE(3,'(a,a)') '# masses                        = ', TRIM(imf)
   IF (TRIM(imf) == 'maschberger') THEN
     WRITE(3,'(a,5es24.16)') '# imf_alpha_beta_mu_mlo_mup     = ', alpha, beta, mu, mlo, mup
   ELSE IF (TRIM(imf) == 'equal') THEN
     WRITE(3,'(a,es24.16)') '# mass_of_each_star_Msun        = ', mequal
   END IF
   WRITE(3,'(a,i0)') '# random_seed                   = ', seedmem
! an ensemble: which member this is. The members differ only in the seed
   IF (nensemble > 1) THEN
     WRITE(3,'(a,a)') '# ensemble                      = ', TRIM(outfile)
     WRITE(3,'(a,i0)') '# ensemble_member               = ', imem
     WRITE(3,'(a,i0)') '# ensemble_size                 = ', nensemble
     WRITE(3,'(a,i0)') '# ensemble_first_seed           = ', seed
   END IF
   IF (setvirial .AND. .NOT. larson) &
&    WRITE(3,'(a,es24.16)') '# virial_ratio_asked_for        = ', qvir
   WRITE(3,'(a,es24.16)') '# virial_ratio                  = ', qmade
   WRITE(3,'(a,es24.16)') '# time_Myr                      = ', 0.d0
   WRITE(3,'(a,i0)') '# n                             = ', n
   WRITE(3,'(a,es24.16)') '# total_mass_Msun               = ', mtot
   WRITE(3,'(a)') '# units                         = mass Msun, length pc, '// &
&    'velocity km/s, time Myr'
   WRITE(3,'(a)') '# frame                         = centre of mass'
   WRITE(3,'(a)') '# columns                       = id m x y z vx vy vz type'
   WRITE(3,'(a)') '# column_units                  = - Msun pc pc pc km/s km/s km/s -'
   WRITE(3,'(a)') '# end_header'
   DO i=1,n
     WRITE(3,'(i10,7es24.16,1x,a)') i, m(i), r(1:3,i), v(1:3,i), TRIM(ptype(i))
   END DO
   CLOSE(3)
!
   END SUBROUTINE writeics
! =====================================================================
   SUBROUTINE readinput
!
! Set the standard value of every parameter, read the parameter file
! and check it. &ics must be there; &masses and &fractal are optional.
! The names of distribution, imf and velocities are checked where they
! are used (makestructure, makemasses), so a new one is added there only.
!
   INTEGER :: ios
   CHARACTER(len=256) :: msg
! >>> NEW PARAMETER: add it to the namelist of its section
   NAMELIST /ics/ outfile, n, seed, distribution, radius, qvir, setvirial, &
&    nensemble
   NAMELIST /masses/ imf, mlo, mup, alpha, beta, mu, mequal
   NAMELIST /fractal/ fdim, ndiv, velocities, coresize, vlarson
!
! standard values
! >>> NEW PARAMETER: give it a standard value here
   outfile=''
   n=1000
   seed=12345
   nensemble=1
   distribution='plummer'
   radius=1.d0
   qvir=0.5d0
   setvirial=.TRUE.
   imf='maschberger'
   mlo=0.1d0
   mup=2.d0
   alpha=2.3d0
   beta=1.4d0
   mu=0.2d0
   mequal=1.d0
   fdim=1.6d0
   ndiv=2
   velocities='standard'
   coresize=0.d0
   vlarson=1.d0
!
   OPEN(2,file=TRIM(pfile),status='old',action='read',iostat=ios)
   IF (ios /= 0) CALL bad('cannot open parameter file '//TRIM(pfile))
   READ(2,nml=ics,iostat=ios,iomsg=msg)
   IF (ios < 0) CALL bad('no &ics section found in '//TRIM(pfile))
   CALL nmlerror('ics',ios,msg)
   REWIND(2)
   READ(2,nml=masses,iostat=ios,iomsg=msg)
   CALL nmlerror('masses',ios,msg)
   REWIND(2)
   READ(2,nml=fractal,iostat=ios,iomsg=msg)
   CALL nmlerror('fractal',ios,msg)
   CLOSE(2)
!
   outfile=ADJUSTL(outfile)
   distribution=ADJUSTL(distribution)
   imf=ADJUSTL(imf)
   velocities=ADJUSTL(velocities)
!
! checks
   IF (LEN_TRIM(outfile) == 0) CALL bad('outfile must be set in &ics')
   IF (INDEX(TRIM(outfile),' ') > 0) CALL bad('outfile name must not contain spaces')
   IF (n < 2) CALL bad('n must be at least 2')
   IF (nensemble < 1) CALL bad('nensemble must be at least 1')
   IF (radius <= 0.d0) CALL bad('radius must be > 0')
   IF (qvir < 0.d0) CALL bad('qvir must be >= 0')
   IF (TRIM(imf) == 'maschberger') THEN
     IF (mlo <= 0.d0 .OR. mup <= mlo) CALL bad('need 0 < mlo < mup')
     IF (mu <= 0.d0) CALL bad('mu must be > 0')
     IF (alpha <= 1.d0 .OR. beta <= 1.d0) CALL bad('alpha and beta must be > 1')
   END IF
   IF (TRIM(imf) == 'equal' .AND. mequal <= 0.d0) CALL bad('mequal must be > 0')
   IF (TRIM(distribution) == 'fractal') THEN
     IF (fdim <= 0.d0 .OR. fdim > 3.d0) CALL bad('fdim must be between 0 and 3')
     IF (ndiv < 2 .OR. ndiv > 10) CALL bad('ndiv must be between 2 and 10')
     IF (coresize < 0.d0 .OR. coresize >= radius) &
&      CALL bad('coresize must be 0, or between 0 and radius')
     IF (vlarson <= 0.d0) CALL bad('vlarson must be > 0')
   END IF
!
   END SUBROUTINE readinput
! =====================================================================
   SUBROUTINE nmlerror(name,ios,msg)
! stop with the compiler's message if a section of the parameter file
! is wrong (ios > 0); a missing section (ios < 0) keeps its standard values
   CHARACTER(len=*), INTENT(IN) :: name,msg
   INTEGER, INTENT(IN) :: ios
   IF (ios > 0) THEN
     WRITE(6,'(a)') ' ERROR in &'//name//' section of '//TRIM(pfile)//':'
     WRITE(6,'(a,a)') '   ', TRIM(msg)
     STOP
   END IF
   END SUBROUTINE nmlerror
! =====================================================================
   SUBROUTINE bad(text)
! stop with an error message
   CHARACTER(len=*), INTENT(IN) :: text
   WRITE(6,'(a,a)') ' ERROR: ', text
   STOP
   END SUBROUTINE bad
!
   END PROGRAM makeics
