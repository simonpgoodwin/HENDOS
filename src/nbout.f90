! ==============================================================
   MODULE nbout
!
! HENDOS - HErmite N-body Dynamics Of Systems
! all output: screen (standard or verbose) and snapshot files
! everything written is in physical units - Msun, pc, km/s, Myr, yr
! binaries treated as one body are invisible here: energies include
! their internal energy, and snapshots contain all ntot stars
!
   USE nbvars
   USE nbforce
   USE nbbin
   IMPLICIT NONE
! diagnostics of the current (synchronised) state - set by diagnose
   DOUBLE PRECISION :: ekin,epot,etid,etot  ! energies, internal units
   DOUBLE PRECISION :: etot0                ! initial energy
   DOUBLE PRECISION :: eprev                ! energy at previous output
   DOUBLE PRECISION :: demax                ! largest |dE/E| so far
   DOUBLE PRECISION :: rlag(1:3)            ! 10, 50, 90% Lagrangian radii (pc)
   DOUBLE PRECISION :: cm(1:3),cv(1:3)      ! centre of mass (pc, internal)
   DOUBLE PRECISION :: sig3d                ! 3D velocity dispersion (internal)
! counters at the previous output, for per-interval numbers
   INTEGER(KIND=8) :: nblock0,nforce0,nclamp0,nclosed0,nopened0
   DOUBLE PRECISION :: wlast
   INTEGER(KIND=8) :: clock0,crate
   CHARACTER(len=78), PARAMETER :: dline=REPEAT('=',78)
   CHARACTER(len=78), PARAMETER :: sline=REPEAT('-',78)
!
   CONTAINS
! =====================================================================
   FUNCTION hms(sec) RESULT(str)
!
! seconds as h:mm:ss, right adjusted
!
   DOUBLE PRECISION, INTENT(IN) :: sec
   CHARACTER(len=12) :: str
   INTEGER(KIND=8) :: is,h,mi,s
!
   is=NINT(MAX(sec,0.d0),KIND=8)
   h=is/3600_8
   mi=MOD(is,3600_8)/60_8
   s=MOD(is,60_8)
   WRITE(str,'(i0,a,i2.2,a,i2.2)') h, ':', mi, ':', s
   str=ADJUSTR(str)
!
   END FUNCTION hms
! =====================================================================
   FUNCTION walltime() RESULT(w)
!
! wall-clock seconds since the start of the integration
!
   DOUBLE PRECISION :: w
   INTEGER(KIND=8) :: c
!
   CALL SYSTEM_CLOCK(c)
   w=DBLE(c - clock0)/DBLE(crate)
!
   END FUNCTION walltime
! =====================================================================
   FUNCTION typecounts() RESULT(str)
!
! the number of particles of each type present, as text, e.g.
! 'star 990, black_hole 10' (stars inside closed binaries are counted)
!
   CHARACTER(len=256) :: str
   CHARACTER(len=16) :: num
   INTEGER :: cnt(1:ntypes),i,k
!
   cnt=0
   DO i=1,n
     k=ib(i)
     IF (k == 0) THEN
       cnt(itype(i))=cnt(itype(i)) + 1
     ELSE
       cnt(btype(1,k))=cnt(btype(1,k)) + 1
       cnt(btype(2,k))=cnt(btype(2,k)) + 1
     END IF
   END DO
   str=''
   DO k=1,ntypes
     IF (cnt(k) == 0) CYCLE
     WRITE(num,'(i0)') cnt(k)
     IF (LEN_TRIM(str) > 0) str=TRIM(str)//','
     str=TRIM(str)//' '//TRIM(typename(k))//' '//TRIM(num)
   END DO
   str=ADJUSTL(str)
!
   END FUNCTION typecounts
! =====================================================================
   SUBROUTINE indexsort(nn,a,idx)
!
! heapsort index: a(idx(1)) <= a(idx(2)) <= ... (a is not changed)
!
   INTEGER, INTENT(IN) :: nn
   DOUBLE PRECISION, INTENT(IN) :: a(1:nn)
   INTEGER, INTENT(OUT) :: idx(1:nn)
   INTEGER :: i,j,l,ir,indxt
   DOUBLE PRECISION :: q
!
   DO i=1,nn
     idx(i)=i
   END DO
   IF (nn < 2) RETURN
   l=nn/2 + 1
   ir=nn
   DO
     IF (l > 1) THEN
       l=l - 1
       indxt=idx(l)
       q=a(indxt)
     ELSE
       indxt=idx(ir)
       q=a(indxt)
       idx(ir)=idx(1)
       ir=ir - 1
       IF (ir == 1) THEN
         idx(1)=indxt
         EXIT
       END IF
     END IF
     i=l
     j=l + l
     DO WHILE (j <= ir)
       IF (j < ir) THEN
         IF (a(idx(j)) < a(idx(j+1))) j=j + 1
       END IF
       IF (q < a(idx(j))) THEN
         idx(i)=idx(j)
         i=j
         j=j + j
       ELSE
         j=ir + 1
       END IF
     END DO
     idx(i)=indxt
   END DO
!
   END SUBROUTINE indexsort
! =====================================================================
   SUBROUTINE diagnose(time)
!
! energies (and potential per star, phi), centre of mass, velocity
! dispersion and 10/50/90% Lagrangian radii about the centre of mass
!
   DOUBLE PRECISION, INTENT(IN) :: time
   INTEGER :: i,k
   DOUBLE PRECISION :: mt,frac(1:3),mb,rr(1:3),vr(1:3)
!
   frac=(/0.1d0, 0.5d0, 0.9d0/)
   CALL energy(ekin,epot,etid)
! add the internal energy of binaries treated as one body, at the
! current phase of their orbits
   DO i=1,n
     k=ib(i)
     IF (k == 0) CYCLE
     mb=bmass(1,k) + bmass(2,k)
     CALL kepler(mb,brel(1:3,k),brel(4:6,k),time/tunit - bepoch(k),rr,vr)
     ekin=ekin + 0.5d0*bmass(1,k)*bmass(2,k)/mb*SUM(vr*vr)
     epot=epot - bmass(1,k)*bmass(2,k)/SQRT(SUM(rr*rr))
   END DO
   etot=ekin + epot + etid
!
   mt=SUM(m(1:n))
   cm=0.d0
   cv=0.d0
   DO i=1,n
     cm(1:3)=cm(1:3) + m(i)*r(1:3,i)
     cv(1:3)=cv(1:3) + m(i)*v(1:3,i)
   END DO
   cm=cm/mt
   cv=cv/mt
   sig3d=0.d0
   DO i=1,n
     sig3d=sig3d + m(i)*SUM((v(1:3,i) - cv(1:3))**2)
   END DO
   sig3d=SQRT(sig3d/mt)
!
   CALL lagradii(3,frac,rlag)
!
   END SUBROUTINE diagnose
! =====================================================================
   SUBROUTINE lagradii(nf,frac,rl)
!
! Lagrangian radii (pc) about the centre of mass cm for the mass
! fractions frac(1:nf), which must be in increasing order
!
   INTEGER, INTENT(IN) :: nf
   DOUBLE PRECISION, INTENT(IN) :: frac(1:nf)
   DOUBLE PRECISION, INTENT(OUT) :: rl(1:nf)
   INTEGER :: i,k
   INTEGER, ALLOCATABLE :: idx(:)
   DOUBLE PRECISION, ALLOCATABLE :: r2(:)
   DOUBLE PRECISION :: mt,mc
!
   ALLOCATE(r2(1:n))
   ALLOCATE(idx(1:n))
   DO i=1,n
     r2(i)=SUM((r(1:3,i) - cm(1:3))**2)
   END DO
   CALL indexsort(n,r2,idx)
   mt=SUM(m(1:n))
   rl=SQRT(r2(idx(n)))
   mc=0.d0
   k=1
   DO i=1,n
     mc=mc + m(idx(i))
     DO WHILE (k <= nf)
       IF (mc < frac(k)*mt) EXIT
       rl(k)=SQRT(r2(idx(i)))
       k=k + 1
     END DO
     IF (k > nf) EXIT
   END DO
   DEALLOCATE(r2,idx)
!
   END SUBROUTINE lagradii
! =====================================================================
   SUBROUTINE levelhist(hist)
!
! number of stars on each timestep level
!
   INTEGER, INTENT(OUT) :: hist(0:60)
   INTEGER :: i
!
   hist=0
   DO i=1,n
     hist(lev(i))=hist(lev(i)) + 1
   END DO
!
   END SUBROUTINE levelhist
! =====================================================================
   SUBROUTINE writelevels
!
! table of occupied timestep levels (verbose mode)
!
   INTEGER :: hist(0:60),k
!
   CALL levelhist(hist)
   WRITE(6,'(a)') '   timestep levels       level      step (yr)       stars'
   DO k=0,kmax
     IF (hist(k) > 0) WRITE(6,'(22x,i6,es15.4,i12)') k, &
&                       tout*1.d6/2.d0**k, hist(k)
   END DO
!
   END SUBROUTINE writelevels
! =====================================================================
   SUBROUTINE startinfo(nout)
!
! screen output at the start of the run
!
   INTEGER, INTENT(IN) :: nout
   DOUBLE PRECISION :: mt,fconv,rjac,tcr,trh,mmin,mmax
   INTEGER :: i,k
!
   CALL diagnose(time0)
   etot0=etot
   eprev=etot
   demax=0.d0
   mt=SUM(m(1:n))
   fconv=1.d0/(vunit*1.d3)
!
! standard - always written
   WRITE(6,'(a)') dline
   WRITE(6,'(a,a)') '  HENDOS - HErmite N-body Dynamics Of Systems     run: ', TRIM(runid)
   WRITE(6,'(a)') dline
   WRITE(6,'(a,a)')          '  initial conditions     ', TRIM(icfile)
   WRITE(6,'(a,i14)')        '  particles              ', ntot
   WRITE(6,'(a,a)')          '  particle types         ', TRIM(typecounts())
   WRITE(6,'(a,f14.2,a)')    '  total mass             ', mt, ' Msun'
   WRITE(6,'(a,f14.4,a)')    '  half-mass radius       ', rlag(2), ' pc'
   WRITE(6,'(a,f14.4)')      '  virial ratio           ', -ekin/epot
   WRITE(6,'(a,f14.4,a,f12.4,a)') '  run from               ', time0, &
&                            ' Myr  to', time0 + DBLE(nout)*tout, ' Myr'
   IF (LEN_TRIM(datadir) > 0) THEN
     WRITE(6,'(a,a)')        '  output files in        ', TRIM(datadir)
   ELSE
     WRITE(6,'(a)')          '  output files in        ./'
   END IF
   WRITE(6,'(a,f14.4,a)')    '  output every           ', tout, ' Myr'
   IF (tsnap > 0.d0) THEN
     WRITE(6,'(a,f14.4,a)')  '  snapshots every        ', tsnap, ' Myr'
   ELSE
     WRITE(6,'(a)')          '  snapshots                       none'
   END IF
   IF (ltide) THEN
     WRITE(6,'(a,a14)')      '  Galactic tide          ', TRIM(tidal)
     IF (ttid(1,1) > 0.d0) THEN
       rjac=(mt/ttid(1,1))**(1.d0/3.d0)
       WRITE(6,'(a,f14.4,a)') '  Jacobi radius          ', rjac, ' pc'
     END IF
   ELSE
     WRITE(6,'(a)')          '  Galactic tide                   none'
   END IF
   IF (usetree) THEN
     WRITE(6,'(a,f14.2,a)')  '  forces                 ', theta, &
&                            '  (tree opening angle; energies direct)'
   ELSE
     WRITE(6,'(a)')          '  forces                        direct'
   END IF
   WRITE(6,'(a,a14)')        '  screen output          ', TRIM(screen)
!
! verbose - parameters and initial conditions in detail
   IF (lverbose) THEN
     WRITE(6,'(a)') sline
     WRITE(6,'(a,a)') '  PARAMETERS from ', TRIM(pfilename)
     IF (LEN_TRIM(nmlmissing) > 0) WRITE(6,'(a,a)') &
&      '  sections not in file, standard values used: ', TRIM(nmlmissing)
     WRITE(6,'(a,es14.4,a)') '    eps        ', eps, ' pc'
     WRITE(6,'(a,f14.4)')    '    eta        ', eta
     WRITE(6,'(a,f14.4)')    '    etap       ', etap
     WRITE(6,'(a,f14.4)')    '    etas       ', etas
     WRITE(6,'(a,i14)')      '    niter      ', niter
     WRITE(6,'(a,i14,a,es10.3,a)') '    kmax       ', kmax, &
&                            '   (shortest step', tout*1.d6/2.d0**kmax, ' yr)'
     WRITE(6,'(a,i14)')      '    nompmin    ', nompmin
     WRITE(6,'(a,i14)')      '    nompj      ', nompj
     WRITE(6,'(a,es14.4)')   '    dewarn     ', dewarn
     WRITE(6,'(a,l14)')      '    closebinaries', closebinaries
     WRITE(6,'(a,f14.4,a)')  '    rcoll      ', rcoll, ' Rsun (0 = no check)'
     WRITE(6,'(a,l14)')      '    usetree    ', usetree
     IF (usetree) THEN
       WRITE(6,'(a,f14.4)')  '    theta      ', theta
       WRITE(6,'(a,i14)')    '    nleaf      ', nleaf
       WRITE(6,'(a,i14)')    '    ntreebuild ', ntreebuild
     END IF
     IF (closebinaries) THEN
       WRITE(6,'(a,es14.4,a)') '    dtclose    ', dtclose, ' yr'
       WRITE(6,'(a,es14.4)')   '    gclose     ', gclose
       WRITE(6,'(a,es14.4)')   '    gopen      ', gopen
     END IF
     IF (ltide) THEN
       IF (TRIM(tidal) == 'solar') THEN
         WRITE(6,'(a,f14.4,a)') '    oorta      ', oorta, ' km/s/kpc'
         WRITE(6,'(a,f14.4,a)') '    oortb      ', oortb, ' km/s/kpc'
         WRITE(6,'(a,f14.4,a)') '    rho0       ', rho0, ' Msun/pc^3'
       END IF
       WRITE(6,'(a,f14.4,a)')   '    Omega      ', omega/fconv, ' km/s/kpc'
       WRITE(6,'(a,3f12.2,a)')  '    T row 1    ', ttid(1,1:3)/fconv**2, ' (km/s/kpc)^2'
       WRITE(6,'(a,3f12.2,a)')  '    T row 2    ', ttid(2,1:3)/fconv**2, ' (km/s/kpc)^2'
       WRITE(6,'(a,3f12.2,a)')  '    T row 3    ', ttid(3,1:3)/fconv**2, ' (km/s/kpc)^2'
     END IF
!
     mmin=HUGE(1.d0)
     mmax=0.d0
     DO i=1,n
       k=ib(i)
       IF (k == 0) THEN
         mmin=MIN(mmin,m(i))
         mmax=MAX(mmax,m(i))
       ELSE
         mmin=MIN(mmin,bmass(1,k),bmass(2,k))
         mmax=MAX(mmax,bmass(1,k),bmass(2,k))
       END IF
     END DO
     WRITE(6,'(a)') sline
     WRITE(6,'(a)') '  INITIAL CONDITIONS'
     WRITE(6,'(a,a)')         '    file format             ', TRIM(icformat)
     WRITE(6,'(a,3es12.4,a)') '    centre of mass shift    ', comshift, ' pc'
     WRITE(6,'(a,3es12.4,a)') '    c.o.m. velocity shift   ', covshift, ' km/s'
     WRITE(6,'(a,3es12.4)')   '    mass min, mean, max     ', mmin, &
&                             mt/DBLE(ntot), mmax
     WRITE(6,'(a,3f12.4,a)')  '    Lagrangian 10,50,90%    ', rlag, ' pc'
     WRITE(6,'(a,f12.4,a)')   '    3D velocity dispersion  ', sig3d*vunit, ' km/s'
     WRITE(6,'(a,es12.4,a)')  '    kinetic energy          ', ekin*vunit**2, ' Msun km2/s2'
     WRITE(6,'(a,es12.4,a)')  '    potential energy        ', epot*vunit**2, ' Msun km2/s2'
     IF (ltide) &
&    WRITE(6,'(a,es12.4,a)')  '    tidal energy            ', etid*vunit**2, ' Msun km2/s2'
     WRITE(6,'(a,es12.4,a)')  '    total energy            ', etot*vunit**2, ' Msun km2/s2'
! crossing time 2 r_h / sigma_3D and Spitzer half-mass relaxation time
! t_rh = 0.138 N^1/2 r_h^3/2 / ((G M/N)^1/2 ln(0.11 N))
     tcr=2.d0*rlag(2)/sig3d*tunit
     WRITE(6,'(a,f12.4,a)')   '    crossing time 2r_h/sig  ', tcr, ' Myr'
     IF (ntot > 10) THEN
       trh=0.138d0*SQRT(DBLE(ntot))*rlag(2)**1.5d0/SQRT(mt/DBLE(ntot))/ &
&          LOG(0.11d0*DBLE(ntot))*tunit
       WRITE(6,'(a,f12.4,a)') '    half-mass relaxation    ', trh, ' Myr'
     END IF
     IF (closebinaries) WRITE(6,'(a,i12,a,i0,a)') '    binaries as one body    ', &
&      nbin, '   (', n, ' bodies)'
     CALL writelevels
   END IF
   WRITE(6,'(a)') dline
!
! start the clocks
   CALL SYSTEM_CLOCK(clock0,crate)
   wlast=0.d0
   nblock0=0
   nforce0=0
   nclamp0=nclamp
   nclosed0=nclosed
   nopened0=nopened
!
   END SUBROUTINE startinfo
! =====================================================================
   SUBROUTINE statusout(iout,nout,time,nblock,nforce)
!
! screen output at each output time
! standard: one line;  verbose: a block of diagnostics
! warnings are written in both modes
!
   INTEGER, INTENT(IN) :: iout,nout
   DOUBLE PRECISION, INTENT(IN) :: time
   INTEGER(KIND=8), INTENT(IN) :: nblock,nforce
   DOUBLE PRECISION :: de,dei,w,wi,eta_s
   INTEGER(KIND=8) :: nbi,nfi
!
   CALL diagnose(time)
   de=(etot - etot0)/ABS(etot0)
   dei=(etot - eprev)/ABS(etot0)
   demax=MAX(demax,ABS(de))
   w=walltime()
   wi=w - wlast
   eta_s=w*DBLE(nout - iout)/DBLE(iout)
   nbi=nblock - nblock0
   nfi=nforce - nforce0
!
   IF (.NOT. lverbose) THEN
! standard - a header every 50 lines, then one line per output
     IF (MOD(iout - 1,50) == 0) WRITE(6,'(a)') &
&      '   time(Myr)     done       dE/E       Q  r_h (pc)   wall time   time left'
     WRITE(6,'(f12.4,f8.1,a,es11.2,f8.3,f10.4,2a12)') time, &
&      100.d0*DBLE(iout)/DBLE(nout), '%', de, -ekin/epot, rlag(2), &
&      hms(w), hms(eta_s)
   ELSE
! verbose
     WRITE(6,'(a)') ''
     WRITE(6,'(a,f0.4,a,i0,a,i0,a)') ' ---- t = ', time, ' Myr    (output ', &
&      iout, ' of ', nout, ') ----------------------------------------'
     WRITE(6,'(a,es12.4,a,es12.4)') '   energy (Msun km2/s2)     kinetic ', &
&      ekin*vunit**2, '   potential ', epot*vunit**2
     IF (ltide) THEN
       WRITE(6,'(a,es12.4,a,es12.4)') '                            tidal   ', &
&        etid*vunit**2, '   total     ', etot*vunit**2
     ELSE
       WRITE(6,'(a,es12.4)') '                            total   ', etot*vunit**2
     END IF
     WRITE(6,'(a,es11.3,a,es11.3,a,es10.2)') '   energy error   since start', de, &
&      '   this interval', dei, '   max', demax
     WRITE(6,'(a,f10.4,a,f10.4,a)') '   virial ratio            ', -ekin/epot, &
&      '      sigma_3D', sig3d*vunit, ' km/s'
     WRITE(6,'(a,3f10.4)') '   Lagrangian radii (pc)   10/50/90%  ', rlag
     WRITE(6,'(a,3es11.3,a)') '   centre of mass          ', cm, ' pc'
     WRITE(6,'(a,3es11.3,a)') '   c.o.m. velocity         ', cv*vunit, ' km/s'
     WRITE(6,'(a,i12,a,i14)') '   this interval   blocks  ', nbi, &
&      '   star steps', nfi
     IF (nbi > 0) WRITE(6,'(a,f12.1)') '                   mean stars per block', &
&      DBLE(nfi)/DBLE(nbi)
     CALL writelevels
     WRITE(6,'(a,i0)') '   steps held at kmax this interval ', nclamp - nclamp0
     IF (usetree) WRITE(6,'(a,i12,a,i0)') '   tree nodes              ', nnode, &
&      '   builds so far ', ntreebuilt
     IF (closebinaries) WRITE(6,'(a,i0,a,i0,a,i0,a)') '   binaries as one body    ', &
&      nbin, '   (closed ', nclosed - nclosed0, ', opened ', nopened - nopened0, &
&      ' this interval)'
     WRITE(6,'(a,a,a,a,a,a)') '   wall time   interval', hms(wi), '   total', &
&      hms(w), '   left', hms(eta_s)
     IF (wi > 0.d0) WRITE(6,'(a,es11.3,a)') '   speed      ', DBLE(nfi)/wi, &
&      ' star steps per second'
   END IF
!
! warnings - both modes
   IF (ABS(de) > dewarn) WRITE(6,'(a,f0.4,a,es9.2,a,es9.2)') &
&    '  ** WARNING t = ', time, ' Myr: |dE/E| = ', ABS(de), ' is above dewarn =', dewarn
   IF (nclamp > nclamp0) WRITE(6,'(a,f0.4,a,i0,a)') '  ** NOTE t = ', time, &
&    ' Myr: ', nclamp - nclamp0, ' steps were held at the shortest timestep (kmax)'
!
   eprev=etot
   wlast=w
   nblock0=nblock
   nforce0=nforce
   nclamp0=nclamp
   nclosed0=nclosed
   nopened0=nopened
!
   END SUBROUTINE statusout
! =====================================================================
   SUBROUTINE endinfo(time,nforce)
!
! summary at the end of the run - both modes
!
   DOUBLE PRECISION, INTENT(IN) :: time
   INTEGER(KIND=8), INTENT(IN) :: nforce
   DOUBLE PRECISION :: w
!
   w=walltime()
   WRITE(6,'(a)') dline
   WRITE(6,'(a,f14.4,a)')  '  finished at            ', time, ' Myr'
   WRITE(6,'(a,es14.3)')   '  final energy error     ', (etot - etot0)/ABS(etot0)
   WRITE(6,'(a,es14.3)')   '  largest energy error   ', demax
   WRITE(6,'(a,a14)')      '  wall time              ', ADJUSTL(hms(w))
   WRITE(6,'(a,i14)')      '  star steps             ', nforce
   IF (w > 0.d0) WRITE(6,'(a,es14.3,a)') '  speed                  ', &
&    DBLE(nforce)/w, ' star steps per second'
   IF (time > time0) WRITE(6,'(a,es14.3)') '  steps per star per Myr ', &
&    DBLE(nforce)/DBLE(ntot)/(time - time0)
   IF (lverbose .AND. closebinaries) WRITE(6,'(a,i0,a,i0,a,i0)') &
&    '  binaries closed ', nclosed, ', opened ', nopened, ', still closed ', nbin
   WRITE(6,'(a)') dline
!
   END SUBROUTINE endinfo
! =====================================================================
   SUBROUTINE hdrs(u,key,val)
! header line '# key = value' on unit u, text value
   INTEGER, INTENT(IN) :: u
   CHARACTER(len=*), INTENT(IN) :: key,val
   WRITE(u,'(a,a,a,a)') '# ', key//REPEAT(' ',MAX(0,30-LEN(key))), '= ', TRIM(val)
   END SUBROUTINE hdrs
!
   SUBROUTINE hdrd(u,key,val)
! header line with one real value
   INTEGER, INTENT(IN) :: u
   CHARACTER(len=*), INTENT(IN) :: key
   DOUBLE PRECISION, INTENT(IN) :: val
   CHARACTER(len=32) :: str
   WRITE(str,'(es24.16)') val
   CALL hdrs(u,key,ADJUSTL(str))
   END SUBROUTINE hdrd
!
   SUBROUTINE hdrv(u,key,val,nv)
! header line with nv (up to 9) real values
   INTEGER, INTENT(IN) :: u,nv
   CHARACTER(len=*), INTENT(IN) :: key
   DOUBLE PRECISION, INTENT(IN) :: val(1:nv)
   CHARACTER(len=512) :: str
   WRITE(str,'(9es24.16)') val(1:nv)
   CALL hdrs(u,key,ADJUSTL(str))
   END SUBROUTINE hdrv
!
   SUBROUTINE hdri(u,key,val)
! header line with an integer value
   INTEGER, INTENT(IN) :: u
   CHARACTER(len=*), INTENT(IN) :: key
   INTEGER(KIND=8), INTENT(IN) :: val
   CHARACTER(len=24) :: str
   WRITE(str,'(i0)') val
   CALL hdrs(u,key,str)
   END SUBROUTINE hdri
!
   SUBROUTINE hdrcommon(u)
! metadata common to every output file: run, units and frame
   INTEGER, INTENT(IN) :: u
   CALL hdrs(u,'creator','HENDOS (hermite.f90)')
   CALL hdrs(u,'date',datestamp())
   CALL hdrs(u,'runid',runid)
   CALL hdrs(u,'parameter_file',pfilename)
   CALL hdrs(u,'initial_conditions',icfile)
   CALL hdri(u,'n',INT(ntot,KIND=8))
   CALL hdrd(u,'total_mass_Msun',SUM(m(1:n)))
   CALL hdrs(u,'units','mass Msun, length pc, velocity km/s, time Myr, '// &
&    'energy Msun (km/s)^2')
   IF (ltide) THEN
     CALL hdrs(u,'frame','rotating with circular Galactic orbit: x away '// &
&      'from Galactic centre, y along rotation')
   ELSE
     CALL hdrs(u,'frame','centre of mass of the initial conditions, non-rotating')
   END IF
   END SUBROUTINE hdrcommon
! =====================================================================
   FUNCTION datestamp() RESULT(str)
   CHARACTER(len=19) :: str
   INTEGER :: dv(1:8)
   CALL DATE_AND_TIME(values=dv)
   WRITE(str,'(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') dv(1), '-', dv(2), &
&    '-', dv(3), ' ', dv(5), ':', dv(6), ':', dv(7)
   END FUNCTION datestamp
! =====================================================================
   FUNCTION snapname(time) RESULT(fname)
!
! snapshot file name  <runid>_snap_<time>Myr.dat  (in the data
! directory, which is included in the name returned) with the time zero
! padded to a fixed width (6 digits . 4 decimals), so that ls lists a
! run's snapshots in time order, e.g. test01_snap_000012.2500Myr.dat
!
   DOUBLE PRECISION, INTENT(IN) :: time
   CHARACTER(len=512) :: fname
   CHARACTER(len=16) :: tstr
   INTEGER(KIND=8) :: it
!
   it=NINT(time*1.d4,KIND=8)
   WRITE(tstr,'(i6.6,a,i4.4)') it/10000_8, '.', MOD(it,10000_8)
   fname=TRIM(datadir)//TRIM(runid)//'_snap_'//TRIM(tstr)//'Myr.dat'
!
   END FUNCTION snapname
! =====================================================================
   SUBROUTINE snapshot(time,nforce)
!
! write snapshot <runid>_snap_<time>Myr.dat (see snapname)
! must be called just after diagnose (statusout or startinfo), so the
! energies, phi and Lagrangian radii are for the current state.
!
! FORMAT (NBODY-SNAPSHOT 2) - plain text
!   header lines '# key = value', ending with '# end_header'
!   then one line per star:  id m x y z vx vy vz type phi dt
!   type is the name of the particle type (see typename in nbvars)
! The initial conditions generators write the same format (first 9
! columns only), and any snapshot can be used as initial conditions.
! (NBODY-SNAPSHOT 1 was the same without the type column.)
!
   DOUBLE PRECISION, INTENT(IN) :: time
   INTEGER(KIND=8), INTENT(IN) :: nforce
   INTEGER :: i,k,l,ns
   INTEGER, ALLOCATABLE :: sid(:),sord(:),styp(:)
   DOUBLE PRECISION, ALLOCATABLE :: sdat(:,:),skey(:)
   CHARACTER(len=512) :: fname
   DOUBLE PRECISION :: fconv,tv(1:9),phit,tx(1:3),dtyr,sep
   DOUBLE PRECISION :: x1(1:3),u1(1:3),x2(1:3),u2(1:3)
!
   fname=snapname(time)
   fconv=1.d0/(vunit*1.d3)
!
   OPEN(4,file=TRIM(fname),status='replace')
   CALL hdrs(4,'format','NBODY-SNAPSHOT 2')
   CALL hdrcommon(4)
   CALL hdrs(4,'particle_types',typecounts())
   CALL hdrd(4,'time_Myr',time)
   CALL hdrs(4,'tidal',tidal)
   IF (ltide) THEN
     CALL hdrd(4,'tidal_omega_km_s_kpc',omega/fconv)
     tv=RESHAPE(TRANSPOSE(ttid),(/9/))/fconv**2
     CALL hdrv(4,'tidal_tensor_km2_s2_kpc2',tv,9)
     IF (TRIM(tidal) == 'solar') THEN
       CALL hdrd(4,'oort_A_km_s_kpc',oorta)
       CALL hdrd(4,'oort_B_km_s_kpc',oortb)
       CALL hdrd(4,'rho0_Msun_pc3',rho0)
     END IF
   END IF
   CALL hdrd(4,'softening_pc',eps)
   CALL hdrd(4,'eta',eta)
   CALL hdrd(4,'etap',etap)
   CALL hdri(4,'niter',INT(niter,KIND=8))
   CALL hdri(4,'kmax',INT(kmax,KIND=8))
   CALL hdrd(4,'dtmax_yr',tout*1.d6)
   CALL hdrd(4,'energy_kinetic',ekin*vunit**2)
   CALL hdrd(4,'energy_potential',epot*vunit**2)
   CALL hdrd(4,'energy_tidal',etid*vunit**2)
   CALL hdrd(4,'energy_total',etot*vunit**2)
   CALL hdrd(4,'energy_initial',etot0*vunit**2)
   CALL hdrd(4,'energy_error',(etot - etot0)/ABS(etot0))
   CALL hdrd(4,'energy_error_max',demax)
   CALL hdrd(4,'virial_ratio',-ekin/epot)
   CALL hdrv(4,'lagrangian_radii_10_50_90_pc',rlag,3)
   CALL hdrv(4,'centre_of_mass_pc',cm,3)
   CALL hdrv(4,'centre_of_mass_velocity_km_s',cv*vunit,3)
   CALL hdrd(4,'wall_time_s',walltime())
   CALL hdri(4,'star_steps',nforce)
   CALL hdrs(4,'columns','id m x y z vx vy vz type phi dt')
   CALL hdrs(4,'column_units','- Msun pc pc pc km/s km/s km/s - (km/s)^2 yr')
   CALL hdrs(4,'column_notes','type = particle type; phi = gravitational + '// &
&    'tidal potential per unit mass; dt = current timestep; rows in id order')
   WRITE(4,'(a)') '# end_header'
!
! all ntot stars: binaries treated as one body are split into their two
! stars (phi of each includes its companion), then written in id order
   ALLOCATE(sid(1:ntot),sdat(1:9,1:ntot),sord(1:ntot),skey(1:ntot))
   ALLOCATE(styp(1:ntot))
   ns=0
   DO i=1,n
     phit=phi(i)
     IF (ltide) THEN
       tx(1:3)=MATMUL(ttid,r(1:3,i))
       phit=phit - 0.5d0*(r(1,i)*tx(1) + r(2,i)*tx(2) + r(3,i)*tx(3))
     END IF
     dtyr=DBLE(dtk(i))*tick*tunit*1.d6
     k=ib(i)
     IF (k == 0) THEN
       ns=ns + 1
       sid(ns)=id(i)
       styp(ns)=itype(i)
       sdat(1:9,ns)=(/m(i), r(1:3,i), v(1:3,i)*vunit, phit*vunit**2, dtyr/)
     ELSE
       CALL bincomp(k,r(1:3,i),v(1:3,i),time/tunit,x1,u1,x2,u2)
       sep=SQRT(SUM((x1 - x2)**2))
       ns=ns + 1
       sid(ns)=bident(1,k)
       styp(ns)=btype(1,k)
       sdat(1:9,ns)=(/bmass(1,k), x1, u1*vunit, (phit - bmass(2,k)/sep)*vunit**2, dtyr/)
       ns=ns + 1
       sid(ns)=bident(2,k)
       styp(ns)=btype(2,k)
       sdat(1:9,ns)=(/bmass(2,k), x2, u2*vunit, (phit - bmass(1,k)/sep)*vunit**2, dtyr/)
     END IF
   END DO
   skey(1:ns)=DBLE(sid(1:ns))
   CALL indexsort(ns,skey(1:ns),sord(1:ns))
   DO l=1,ns
     i=sord(l)
     WRITE(4,'(i10,7es24.16,1x,a12,2es16.8)') sid(i), sdat(1:7,i), &
&      typename(styp(i)), sdat(8:9,i)
   END DO
   DEALLOCATE(sid,sdat,sord,skey,styp)
   CLOSE(4)
   IF (lverbose) WRITE(6,'(a,a)') '   snapshot written: ', TRIM(fname)
!
   END SUBROUTINE snapshot
!
   END MODULE nbout
