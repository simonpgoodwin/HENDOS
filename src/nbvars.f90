   MODULE nbvars
! HENDOS - HErmite N-body Dynamics Of Systems
! module for main variables
! internally G=1, mass unit = 1 Msun, length unit = 1 pc
! all input and output is in Msun, pc, km/s and Myr - conversions
! happen only on reading and writing
    IMPLICIT NONE
! particle types -----------------------------------------------
! every particle has a type, held as an integer code itype(i) which
! indexes typename; the names are what appear in all files.
! TO ADD A NEW TYPE: increase ntypes and add its name to typename (no
! spaces, at most 12 characters, lower case, padded to 12 with blanks),
! and a named code below if the code needs to refer to it. Nothing
! else needs changing: reading, snapshots and screen output all use
! this table.
    INTEGER, PARAMETER :: ntypes=5
    CHARACTER(len=12), PARAMETER :: typename(1:ntypes)=(/ 'star        ', &
&     'white_dwarf ', 'neutron_star', 'black_hole  ', 'planet      ' /)
    INTEGER, PARAMETER :: tstar=1,twdwarf=2,tnstar=3,tbhole=4,tplanet=5
! parameters read from the parameter file (see readparams) ----
!   run set-up, physical units
    CHARACTER(len=32) :: runid       ! run name, used in output file names
    CHARACTER(len=256) :: icfile     ! initial conditions file
    DOUBLE PRECISION :: tend         ! end time (Myr)
    DOUBLE PRECISION :: tout         ! output interval and largest step (Myr)
    DOUBLE PRECISION :: eps          ! softening length (pc)
!   integrator
    DOUBLE PRECISION :: eta          ! Aarseth timestep parameter
    DOUBLE PRECISION :: etap         ! pairwise free-fall/flyby timestep parameter
    DOUBLE PRECISION :: etas         ! starting timestep parameter
    INTEGER :: niter                 ! number of corrector iterations P(EC)^niter
    INTEGER :: kmax                  ! deepest allowed timestep level, dtmin=dtmax/2**kmax
!   Galactic tide (input units km/s/kpc)
    CHARACTER(len=16) :: tidal       ! 'none', 'solar' or 'tensor'
    DOUBLE PRECISION :: oorta,oortb  ! Oort constants A, B (km/s/kpc)
    DOUBLE PRECISION :: rho0         ! local total mass density (Msun/pc^3)
    DOUBLE PRECISION :: tij(1:3,1:3) ! user tidal tensor, rotating frame ((km/s/kpc)^2)
    DOUBLE PRECISION :: omegat       ! user frame rotation rate about z (km/s/kpc)
!   binaries (no regularisation - see module nbbin)
    LOGICAL :: closebinaries         ! treat unperturbed binaries as one body
    DOUBLE PRECISION :: dtclose      ! close candidates: stars with step < dtclose (yr)
    DOUBLE PRECISION :: gclose       ! close if perturbation gamma < gclose
    DOUBLE PRECISION :: gopen        ! open if perturbation gamma > gopen
!   Barnes-Hut tree (see module nbtree)
    LOGICAL :: usetree               ! tree forces instead of direct summation
    DOUBLE PRECISION :: theta        ! opening angle
    INTEGER :: nleaf                 ! max bodies per leaf
    INTEGER :: ntreebuild            ! rebuild the tree every ntreebuild block steps
!   collisions
    DOUBLE PRECISION :: rcoll        ! stop if two stars come within rcoll (Rsun), 0 = off
!   computing
    INTEGER :: nompmin               ! min. active stars to share them over threads
    INTEGER :: nompj                 ! min. N to share the j-sum (and O(N) loops)
!   output
    CHARACTER(len=256) :: datadir    ! directory (path) for snapshot and analysis files
                                     !   (held with a trailing /, blank = run directory)
    CHARACTER(len=16) :: screen      ! screen output 'standard' or 'verbose'
    DOUBLE PRECISION :: tsnap        ! snapshot interval (Myr), 0 = none
    DOUBLE PRECISION :: dewarn       ! warn if |dE/E| exceeds this
!   analysis switches
    DOUBLE PRECISION :: tana        ! analysis interval (Myr), 0 = none
    LOGICAL :: lagrangian            ! Lagrangian radii
! derived and internal ----------------------------------------
    LOGICAL :: lverbose              ! verbose screen output
    INTEGER :: ntot                  ! number of stars (n = number of bodies)
    DOUBLE PRECISION :: tbase        ! internal time at start of output interval
    DOUBLE PRECISION :: dtclosei     ! dtclose in internal units
    DOUBLE PRECISION :: rcolli       ! rcoll in pc
    LOGICAL :: treevalid=.FALSE.     ! tree topology matches the bodies
    INTEGER(KIND=8) :: ntreebuilt=0  ! number of tree builds so far
    INTEGER :: nbin                  ! binaries currently treated as one body
    INTEGER(KIND=8) :: nclosed,nopened ! binaries closed and opened so far
    CHARACTER(len=256) :: pfilename  ! parameter file name
    CHARACTER(len=64) :: nmlmissing  ! parameter file sections not found
    CHARACTER(len=32) :: icformat    ! format of the initial conditions file
    DOUBLE PRECISION :: time0        ! start time (Myr), from the IC file
    DOUBLE PRECISION :: comshift(1:3),covshift(1:3) ! c.o.m. shift of ICs (pc, km/s)
    INTEGER :: n                     ! number of particles
    INTEGER :: nact                  ! number of active particles in current block
    INTEGER(KIND=8) :: nclamp        ! no. of times a step was limited by kmax
    DOUBLE PRECISION :: eps2         ! softening length squared (0 = unsoftened)
    DOUBLE PRECISION :: dtmax        ! largest step (level 0) = tout, internal units
    DOUBLE PRECISION :: tick         ! integer time unit = dtmax/2**kmax
! Galactic tide in internal units - acceleration on star at x, v is
!   a = ttid.x + 2*omega*(vy, -vx, 0)
! in a frame rotating at omega about z, x pointing away from the
! Galactic centre, y in the direction of Galactic rotation
    LOGICAL :: ltide                 ! tide on/off
    DOUBLE PRECISION :: ttid(1:3,1:3) ! effective tidal tensor (incl. centrifugal)
    DOUBLE PRECISION :: omega        ! frame rotation rate
! unit conversions - set in setunits
    DOUBLE PRECISION :: vunit        ! internal velocity unit in km/s
    DOUBLE PRECISION :: tunit        ! internal time unit in Myr
! particle data
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: r,v,rp,vp,rc,vc
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: f0,f1,f0p,f1p,f2,f3
! 2nd and 3rd derivatives of acceleration at the end of each star's last
! step - used to raise the order of the predictor
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: f2e,f3e
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: m,tpair2
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: phi   ! grav. potential per unit mass
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rmin2 ! closest neighbour distance^2 (last step)
    INTEGER, ALLOCATABLE, DIMENSION(:) :: id             ! permanent star identifier
    INTEGER, ALLOCATABLE, DIMENSION(:) :: itype          ! particle type code (0 = closed binary)
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ib             ! binary index of body, 0 = single
! closed binaries: masses, ids, types, relative r, v (star 1 - star 2) at epoch
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: bmass,brel
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bident,btype
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: bepoch
! block timestep data - times are integers in units of tick
!   lev(i)  : timestep level, dt(i) = dtmax/2**lev(i)
!   tk(i)   : time of last force evaluation of i
!   dtk(i)  : step of i = 2**(kmax-lev(i))
    INTEGER, ALLOCATABLE, DIMENSION(:) :: lev,ilist
    INTEGER(KIND=8), ALLOCATABLE, DIMENSION(:) :: tk,dtk
!
   CONTAINS
! =====================================================================
   INTEGER FUNCTION typecode(name)
!
! integer code of the particle type called name (any case, leading
! blanks ignored), or 0 if there is no such type
!
   CHARACTER(len=*), INTENT(IN) :: name
   CHARACTER(len=12) :: low
   INTEGER :: k,ic
!
   typecode=0
   IF (LEN_TRIM(ADJUSTL(name)) > 12) RETURN
   low=ADJUSTL(name)
   DO k=1,LEN_TRIM(low)
     ic=IACHAR(low(k:k))
     IF (ic >= 65 .AND. ic <= 90) low(k:k)=ACHAR(ic + 32)
   END DO
   DO k=1,ntypes
     IF (low == typename(k)) typecode=k
   END DO
!
   END FUNCTION typecode
!
   END MODULE nbvars
