! ==============================================================
   MODULE nbsetup
!
! HENDOS - HErmite N-body Dynamics Of Systems
! SET-UP AND INPUT
! units, the parameter file, the data directory, array allocation and
! reading the initial conditions - everything done once, before the
! integration starts
!
   USE nbvars
   IMPLICIT NONE
!
   CONTAINS
! =====================================================================
   SUBROUTINE setunits
!
! internal units G=1, M=1 Msun, R=1 pc
! velocity unit = sqrt(G Msun/pc) = 0.06558 km/s
! time unit = pc/velocity unit = 14.91 Myr
!
   DOUBLE PRECISION :: gmsun,pc,myr
!
   gmsun=1.32712440018d20       ! G*Msun in m^3 s^-2 (IAU)
   pc=3.0856775814913673d16     ! m
   myr=3.15576d13               ! s (Julian)
!
   vunit=SQRT(gmsun/pc)/1.d3
   tunit=pc/(vunit*1.d3)/myr
!
   END SUBROUTINE setunits
! =====================================================================
   SUBROUTINE readparams(pfile)
!
! read the parameter file - Fortran namelist groups:
!   &run         run set-up (physical units)   - must be present
!   &output      screen, snapshots             - optional
!   &analysis    in-run analysis               - optional
!   &tides       Galactic tidal field          - optional
!   &tree        Barnes-Hut tree               - optional
!   &integrator  timestep/integrator settings  - optional
!   &binaries    unperturbed binary treatment  - optional
!   &collisions  stellar collision check       - optional
!   &computing   performance settings          - optional
! groups can be in any order; any parameter left out keeps the default
! set below. Values are checked here; they are echoed to the screen by
! startinfo (verbose mode).
!
! tides: in a frame rotating with the cluster's circular orbit (Hill's
! approximation) the extra acceleration is a = T.x - 2 Omega x v.
! tidal='solar' builds T from the Oort constants and local density:
!   Omega = A - B,  T = diag(4A(A-B), 0, -(4 pi G rho0 + 2(A^2-B^2)))
! tidal='tensor' takes T (tij) and Omega (omegat) directly.
!
   CHARACTER(len=*), INTENT(IN) :: pfile
   INTEGER :: ios,k,l
   LOGICAL :: bad
   DOUBLE PRECISION :: fconv,pi
   CHARACTER(len=256) :: msg
   NAMELIST /run/ runid, icfile, tend, tout, eps
   NAMELIST /output/ datadir, screen, tsnap, dewarn
   NAMELIST /analysis/ tana, lagrangian
   NAMELIST /tides/ tidal, oorta, oortb, rho0, tij, omegat
   NAMELIST /integrator/ eta, etap, etas, niter, kmax
   NAMELIST /binaries/ closebinaries, dtclose, gclose, gopen
   NAMELIST /collisions/ rcoll
   NAMELIST /tree/ usetree, theta, nleaf
   NAMELIST /computing/ nompmin, nompj
!
! defaults --------------------------------------------------
   runid=''
   icfile=''
   tend=1.d0
   tout=0.1d0
   eps=0.d0
   datadir='data'
   screen='standard'
   tsnap=0.d0
   dewarn=1.d-4
   tana=0.d0
   lagrangian=.TRUE.
   tidal='none'
   oorta=15.3d0         ! Bovy 2017, MNRAS, 468, L63
   oortb=-11.9d0
   rho0=0.097d0         ! McKee, Parravano & Hollenbach 2015, ApJ, 814, 13
   tij=0.d0
   omegat=0.d0
   eta=0.02d0
   etap=0.2d0
   etas=0.01d0
   niter=1
   kmax=40
   closebinaries=.TRUE.
   rcoll=10.d0
   usetree=.FALSE.
   theta=0.5d0
   nleaf=8
   dtclose=100.d0
   gclose=1.d-6
   gopen=1.d-5
   nompmin=32
   nompj=1024
   nmlmissing=''
!
   OPEN(2,file=TRIM(pfile),status='old',action='read',iostat=ios)
   IF (ios /= 0) THEN
     WRITE(6,'(a,a)') ' ERROR: cannot open parameter file ', TRIM(pfile)
     STOP
   END IF
!
! &run - required
   READ(2,nml=run,iostat=ios,iomsg=msg)
   IF (ios < 0) THEN
     WRITE(6,'(a,a)') ' ERROR: no &run section found in ', TRIM(pfile)
     STOP
   END IF
   CALL nmlerror('run')
!
! optional sections
   REWIND(2)
   READ(2,nml=output,iostat=ios,iomsg=msg)
   CALL nmlerror('output')
   REWIND(2)
   READ(2,nml=analysis,iostat=ios,iomsg=msg)
   CALL nmlerror('analysis')
   REWIND(2)
   READ(2,nml=tides,iostat=ios,iomsg=msg)
   CALL nmlerror('tides')
   REWIND(2)
   READ(2,nml=integrator,iostat=ios,iomsg=msg)
   CALL nmlerror('integrator')
   REWIND(2)
   READ(2,nml=binaries,iostat=ios,iomsg=msg)
   CALL nmlerror('binaries')
   REWIND(2)
   READ(2,nml=collisions,iostat=ios,iomsg=msg)
   CALL nmlerror('collisions')
   REWIND(2)
   READ(2,nml=tree,iostat=ios,iomsg=msg)
   CALL nmlerror('tree')
   REWIND(2)
   READ(2,nml=computing,iostat=ios,iomsg=msg)
   CALL nmlerror('computing')
   CLOSE(2)
!
   runid=ADJUSTL(runid)
   icfile=ADJUSTL(icfile)
   screen=ADJUSTL(screen)
   tidal=ADJUSTL(tidal)
! data directory: any path, relative to the run directory or absolute
! (a leading ~/ is replaced by the home directory); held with one
! trailing /, blank = the run directory
   datadir=ADJUSTL(datadir)
   IF (datadir(1:2) == '~/') THEN
     CALL GET_ENVIRONMENT_VARIABLE('HOME',msg)
     datadir=TRIM(msg)//datadir(2:)
   END IF
   l=LEN_TRIM(datadir)
   DO WHILE (l > 0)
     IF (datadir(l:l) /= '/') EXIT
     datadir(l:l)=' '
     l=l - 1
   END DO
   IF (TRIM(datadir) == '.') datadir=''
   IF (LEN_TRIM(datadir) > 0) datadir=TRIM(datadir)//'/'
!
! checks ----------------------------------------------------
   bad=.FALSE.
   IF (INDEX(TRIM(datadir),' ') > 0) THEN
     WRITE(6,'(a)') ' ERROR: datadir name must not contain spaces'
     bad=.TRUE.
   END IF
! runid: required, used in all file names, so letters, digits - _ . only
   IF (LEN_TRIM(runid) == 0) THEN
     WRITE(6,'(a)') ' ERROR: runid must be set in &run - it is used in '// &
&      'all output file names'
     bad=.TRUE.
   ELSE IF (.NOT. safename(runid)) THEN
     WRITE(6,'(a)') ' ERROR: runid '''//TRIM(runid)//''' may only contain '// &
&      'letters, digits, - _ and . (no spaces)'
     bad=.TRUE.
   END IF
   IF (LEN_TRIM(icfile) == 0) THEN
     WRITE(6,'(a)') ' ERROR: icfile must be set in &run'
     bad=.TRUE.
   ELSE IF (INDEX(TRIM(icfile),' ') > 0) THEN
     WRITE(6,'(a)') ' ERROR: icfile name must not contain spaces'
     bad=.TRUE.
   END IF
   IF (tout <= 0.d0) THEN
     WRITE(6,'(a)') ' ERROR: tout must be > 0'
     bad=.TRUE.
   END IF
   IF (eps < 0.d0) THEN
     WRITE(6,'(a)') ' ERROR: eps must be >= 0'
     bad=.TRUE.
   END IF
   IF (TRIM(screen) == 'verbose') THEN
     lverbose=.TRUE.
   ELSE IF (TRIM(screen) == 'standard') THEN
     lverbose=.FALSE.
   ELSE
     WRITE(6,'(a)') ' ERROR: screen must be ''standard'' or ''verbose'''
     bad=.TRUE.
   END IF
   IF (tsnap < 0.d0 .OR. tana < 0.d0) THEN
     WRITE(6,'(a)') ' ERROR: tsnap and tana must be >= 0 (0 = off)'
     bad=.TRUE.
   END IF
   IF (eta <= 0.d0 .OR. etap <= 0.d0 .OR. etas <= 0.d0) THEN
     WRITE(6,'(a)') ' ERROR: eta, etap and etas must be > 0'
     bad=.TRUE.
   END IF
   IF (niter < 1) THEN
     WRITE(6,'(a)') ' ERROR: niter must be >= 1'
     bad=.TRUE.
   END IF
   IF (kmax < 1 .OR. kmax > 60) THEN
     WRITE(6,'(a)') ' ERROR: kmax must be between 1 and 60'
     bad=.TRUE.
   END IF
!
! tides - convert to internal units, 1 km/s/kpc = fconv internal
   fconv=1.d0/(vunit*1.d3)
   pi=3.141592653589793d0
   ttid=0.d0
   omega=0.d0
   ltide=.FALSE.
   SELECT CASE (TRIM(tidal))
   CASE ('none')
     ltide=.FALSE.
   CASE ('solar')
     ltide=.TRUE.
     IF (rho0 < 0.d0) THEN
       WRITE(6,'(a)') ' ERROR: rho0 must be >= 0'
       bad=.TRUE.
     END IF
     omega=(oorta - oortb)*fconv
     ttid(1,1)=4.d0*oorta*(oorta - oortb)*fconv**2
     ttid(2,2)=0.d0
     ttid(3,3)=-(4.d0*pi*rho0 + 2.d0*(oorta**2 - oortb**2)*fconv**2)
   CASE ('tensor')
     ltide=.TRUE.
     omega=omegat*fconv
! a tidal tensor is symmetric - use the symmetric part, warn if needed
     DO k=1,3
       DO l=1,3
         ttid(k,l)=0.5d0*(tij(k,l) + tij(l,k))*fconv**2
       END DO
     END DO
     IF (MAXVAL(ABS(tij - TRANSPOSE(tij))) > 1.d-10*MAXVAL(ABS(tij))) &
&      WRITE(6,'(a)') ' WARNING: tij is not symmetric - using its symmetric part'
   CASE DEFAULT
     WRITE(6,'(a)') ' ERROR: tidal must be ''none'', ''solar'' or ''tensor'''
     bad=.TRUE.
   END SELECT
!
   IF (closebinaries) THEN
     IF (gclose <= 0.d0 .OR. gopen <= gclose) THEN
       WRITE(6,'(a)') ' ERROR: need 0 < gclose < gopen'
       bad=.TRUE.
     END IF
     IF (dtclose <= 0.d0) THEN
       WRITE(6,'(a)') ' ERROR: dtclose must be > 0'
       bad=.TRUE.
     END IF
! the Kepler orbits of closed binaries are unsoftened
     IF (eps > 0.d0) THEN
       WRITE(6,'(a)') ' WARNING: binary treatment switched off because '// &
&        'eps > 0 (Kepler orbits are unsoftened)'
       closebinaries=.FALSE.
     END IF
   END IF
   dtclosei=dtclose*1.d-6/tunit
   IF (usetree) THEN
     IF (theta <= 0.d0 .OR. theta >= 1.d0) THEN
       WRITE(6,'(a)') ' ERROR: theta must be between 0 and 1'
       bad=.TRUE.
     END IF
     IF (nleaf < 1) THEN
       WRITE(6,'(a)') ' ERROR: nleaf must be >= 1'
       bad=.TRUE.
     END IF
   END IF
   IF (rcoll < 0.d0) THEN
     WRITE(6,'(a)') ' ERROR: rcoll must be >= 0 (0 = no collision check)'
     bad=.TRUE.
   END IF
! solar radius 6.957e8 m (IAU nominal) in pc
   rcolli=rcoll*6.957d8/3.0856775814913673d16
!
   IF (bad) STOP
!
! warnings for intervals that are not whole multiples of tout
   IF (tsnap > 0.d0) THEN
     IF (ABS(tsnap/tout - NINT(tsnap/tout)) > 1.d-6 .OR. tsnap < tout) &
&      WRITE(6,'(a,f0.4,a)') ' WARNING: tsnap is not a multiple of tout, '// &
&      'snapshots every ', DBLE(MAX(1,NINT(tsnap/tout)))*tout, ' Myr'
   END IF
   IF (tana > 0.d0) THEN
     IF (ABS(tana/tout - NINT(tana/tout)) > 1.d-6 .OR. tana < tout) &
&      WRITE(6,'(a,f0.4,a)') ' WARNING: tana is not a multiple of tout, '// &
&      'analysis every ', DBLE(MAX(1,NINT(tana/tout)))*tout, ' Myr'
   END IF
   IF (tsnap > 0.d0 .AND. tout < 1.d-4) &
&    WRITE(6,'(a)') ' WARNING: tout < 1e-4 Myr - snapshot file names '// &
&    '(4 decimal places in Myr) may not be unique'
!
   CONTAINS
!
   SUBROUTINE nmlerror(name)
! stop with the compiler's message on a bad section, note missing ones
   CHARACTER(len=*), INTENT(IN) :: name
   IF (ios < 0) THEN
     nmlmissing=TRIM(nmlmissing)//' &'//name
   ELSE IF (ios > 0) THEN
     WRITE(6,'(a)') ' ERROR in &'//name//' section of '//TRIM(pfile)//':'
     WRITE(6,'(a,a)') '   ', TRIM(msg)
     STOP
   END IF
   END SUBROUTINE nmlerror
!
   LOGICAL FUNCTION safename(str)
! true if str contains only letters, digits, - _ and .
   CHARACTER(len=*), INTENT(IN) :: str
   INTEGER :: i
   CHARACTER :: c
   safename=.TRUE.
   DO i=1,LEN_TRIM(str)
     c=str(i:i)
     IF (.NOT. ((c >= 'a' .AND. c <= 'z') .OR. (c >= 'A' .AND. c <= 'Z') &
&      .OR. (c >= '0' .AND. c <= '9') .OR. c == '-' .OR. c == '_' &
&      .OR. c == '.')) safename=.FALSE.
   END DO
   END FUNCTION safename
!
   END SUBROUTINE readparams
! =====================================================================
   SUBROUTINE checkdatadir
!
! make the data directory if it does not exist, and stop now (rather
! than at the first snapshot) if files cannot be written there
!
   INTEGER :: ios
!
   IF (LEN_TRIM(datadir) == 0) RETURN
   CALL EXECUTE_COMMAND_LINE('mkdir -p '//TRIM(datadir))
   OPEN(9,file=TRIM(datadir)//'.hendos_write_test',status='replace',iostat=ios)
   IF (ios /= 0) THEN
     WRITE(6,'(a)') ' ERROR: cannot write to the data directory '//TRIM(datadir)
     WRITE(6,'(a)') '        Check datadir in the &output section of the '// &
&     'parameter file.'
     STOP
   END IF
   CLOSE(9,status='delete')
!
   END SUBROUTINE checkdatadir
! =====================================================================
   SUBROUTINE allocarrays
!
!
   ALLOCATE(r(1:3,1:n))
   ALLOCATE(v(1:3,1:n))
   ALLOCATE(rp(1:3,1:n))
   ALLOCATE(vp(1:3,1:n))
   ALLOCATE(rc(1:3,1:n))
   ALLOCATE(vc(1:3,1:n))
   ALLOCATE(f0(1:3,1:n))
   ALLOCATE(f1(1:3,1:n))
   ALLOCATE(f0p(1:3,1:n))
   ALLOCATE(f1p(1:3,1:n))
   ALLOCATE(f2(1:3,1:n))
   ALLOCATE(f3(1:3,1:n))
   ALLOCATE(f2e(1:3,1:n))
   ALLOCATE(f3e(1:3,1:n))
   ALLOCATE(m(1:n))
   ALLOCATE(tpair2(1:n))
   ALLOCATE(phi(1:n))
   ALLOCATE(id(1:n))
   ALLOCATE(itype(1:n))
   ALLOCATE(ib(1:n))
   ALLOCATE(rmin2(1:n))
   ALLOCATE(bmass(1:2,1:n/2+1))
   ALLOCATE(bident(1:2,1:n/2+1))
   ALLOCATE(btype(1:2,1:n/2+1))
   ALLOCATE(brel(1:6,1:n/2+1))
   ALLOCATE(bepoch(1:n/2+1))
   ALLOCATE(lev(1:n))
   ALLOCATE(ilist(1:n))
   ALLOCATE(tk(1:n))
   ALLOCATE(dtk(1:n))
!
   END SUBROUTINE allocarrays
! =====================================================================
   SUBROUTINE readic
!
! read the initial conditions file icfile. Two formats are recognised:
!
! NBODY-SNAPSHOT (written by makeics and as snapshots)
!   header lines '# key = value' ending with '# end_header', then one
!   line per star
!     version 2:  id m x y z vx vy vz type [anything else is ignored]
!     version 1:  id m x y z vx vy vz      [anything else is ignored]
!   type is the name of the particle type (star, white_dwarf, ... see
!   typename in nbvars); version 1 files have no types, so every
!   particle is a star. The version is the number after NBODY-SNAPSHOT
!   in the 'format' header line.
!   the start time is taken from 'time_Myr' in the header, so a run can
!   be restarted from any snapshot
!
! legacy (e.g. bigfractals): lines  m x y z vx vy vz, optional # header
!   lines at the top, ids are set to 1..n, the start time to 0 and
!   every particle is a star
!
! units Msun, pc, km/s. Velocities are converted to internal units and
! everything is moved to the centre-of-mass frame.
!
   INTEGER :: i,k,ios,nhead,nhdr,ncol,ieq,ip,ifmt
   LOGICAL :: snapfmt
   DOUBLE PRECISION :: dum(1:8),com(1:3),cov(1:3),mtot
   CHARACTER(len=1024) :: line
   CHARACTER(len=64) :: key
   CHARACTER(len=32) :: tname
!
   OPEN(3,file=TRIM(icfile),status='old',action='read',iostat=ios)
   IF (ios /= 0) THEN
     WRITE(6,'(a,a)') ' ERROR: cannot open initial conditions file ', TRIM(icfile)
     STOP
   END IF
!
! read the header lines (all lines at the top starting with #)
   nhead=0
   nhdr=-1
   time0=0.d0
   snapfmt=.FALSE.
   ifmt=1
   DO
     READ(3,'(a)',iostat=ios) line
     IF (ios /= 0) EXIT
     line=ADJUSTL(line)
     IF (line(1:1) /= '#') EXIT
     nhead=nhead + 1
     ieq=INDEX(line,'=')
     IF (ieq > 2) THEN
       key=ADJUSTL(line(2:ieq-1))
       SELECT CASE (TRIM(key))
       CASE ('format')
         ip=INDEX(line(ieq+1:),'NBODY-SNAPSHOT')
         IF (ip > 0) THEN
           snapfmt=.TRUE.
! format version: the number after NBODY-SNAPSHOT (1 if there is none)
           READ(line(ieq+ip+14:),*,iostat=ios) ifmt
           IF (ios /= 0) ifmt=1
         END IF
       CASE ('time_Myr')
         READ(line(ieq+1:),*,iostat=ios) time0
       CASE ('n')
         READ(line(ieq+1:),*,iostat=ios) nhdr
       END SELECT
     END IF
   END DO
   IF (snapfmt) THEN
     ncol=8
     WRITE(icformat,'(a,i0)') 'NBODY-SNAPSHOT ', ifmt
   ELSE
     ncol=7
     icformat='legacy (m x y z vx vy vz)'
     time0=0.d0
   END IF
!
! count stars
   REWIND(3)
   DO i=1,nhead
     READ(3,*)
   END DO
   n=0
   DO
     READ(3,*,iostat=ios) dum(1:ncol)
     IF (ios /= 0) EXIT
     n=n + 1
   END DO
   IF (n < 2) THEN
     WRITE(6,'(a,a)') ' ERROR: fewer than 2 stars read from ', TRIM(icfile)
     STOP
   END IF
   IF (nhdr >= 0 .AND. nhdr /= n) WRITE(6,'(a,i0,a,i0,a)') &
&    ' WARNING: header of initial conditions says n = ', nhdr, &
&    ' but ', n, ' stars were read'
!
! read them
   CALL allocarrays
   REWIND(3)
   DO i=1,nhead
     READ(3,*)
   END DO
   DO i=1,n
     IF (snapfmt .AND. ifmt >= 2) THEN
! with types: each line is read as text first (skipping blank lines),
! so a missing type is caught rather than taken from the next line
       DO
         READ(3,'(a)',iostat=ios) line
         IF (ios /= 0) EXIT
         IF (LEN_TRIM(line) > 0) EXIT
       END DO
       IF (ios == 0) READ(line,*,iostat=ios) id(i), m(i), r(1:3,i), v(1:3,i), tname
       IF (ios /= 0) THEN
         WRITE(6,'(a,i0,a)') ' ERROR: cannot read particle ', i, ' of '// &
&          TRIM(icfile)//' - each line needs  id m x y z vx vy vz type'
         STOP
       END IF
       itype(i)=typecode(tname)
       IF (itype(i) == 0) THEN
         WRITE(6,'(a,i0,a)') ' ERROR: particle ', i, ' of '//TRIM(icfile)// &
&          ' has unknown type '''//TRIM(tname)//''''
         WRITE(6,'(a,99(1x,a))') '        known types:', &
&          (TRIM(typename(k)),k=1,ntypes)
         STOP
       END IF
     ELSE IF (snapfmt) THEN
       READ(3,*) id(i), m(i), r(1:3,i), v(1:3,i)
       itype(i)=tstar
     ELSE
       READ(3,*) m(i), r(1:3,i), v(1:3,i)
       id(i)=i
       itype(i)=tstar
     END IF
   END DO
   CLOSE(3)
!
! km/s to internal units
   v=v/vunit
!
! move to centre of mass frame
   com=0.d0
   cov=0.d0
   mtot=0.d0
   DO i=1,n
     com(1:3)=com(1:3) + m(i)*r(1:3,i)
     cov(1:3)=cov(1:3) + m(i)*v(1:3,i)
     mtot=mtot + m(i)
   END DO
   com=com/mtot
   cov=cov/mtot
   DO i=1,n
     r(1:3,i)=r(1:3,i) - com(1:3)
     v(1:3,i)=v(1:3,i) - cov(1:3)
   END DO
   comshift=com
   covshift=cov*vunit
!
! all bodies start as single stars
   ntot=n
   ib=0
   nbin=0
   nclosed=0
   nopened=0
!
   END SUBROUTINE readic
!
   END MODULE nbsetup
