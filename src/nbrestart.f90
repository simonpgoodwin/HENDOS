! ==============================================================
   MODULE nbrestart
!
! HENDOS - HErmite N-body Dynamics Of Systems
! RESTART FILES (parameters trestart and restart) - for long runs
!
! A restart file holds the complete state of a run at an output time,
! in binary at full precision: every body with its force derivatives
! and timestep levels, the groups (binaries, few-body groups, triples),
! the neighbour lists and split forces, the counters and the energy the
! run started with. A run continued from it carries on as the original
! would have. (A snapshot can also be used to start a run, but that is
! a new run: text precision, new first timesteps, energy zero reset.)
!
!   writing   every trestart Myr (a whole number of output intervals)
!             and at the end, to <runid>_restart.bin in the data
!             directory. Each file replaces the one before; it is
!             written under a temporary name first, so a crash while
!             writing never destroys the previous file.
!   reading   set  restart = .true.  in &run and run hermite again with
!             the same parameter file. tend may be increased to extend
!             a finished run. The physics switches (tree, neighbour
!             scheme, few-body) must be the same as when it was written.
!
! FILE FORMAT - made to be added to. The file is a list of sections.
! Each section starts with a record holding its name (16 characters)
! and the number of records that follow in it. A reader that does not
! know a section skips that many records, so a new section (stellar
! evolution, say) can be added without breaking anything, and sections
! are only written when their physics is switched on.
!
! >>> TO ADD A SECTION: write it in writerestart with secthead, and add
!     a CASE for its name in readrestart
!
   USE nbvars
   USE nbsetup
   USE nbac
   USE nbout
   IMPLICIT NONE
   INTEGER, PARAMETER :: rstversion=1
!
   CONTAINS
! =====================================================================
   FUNCTION restartname() RESULT(fname)
   CHARACTER(len=512) :: fname
   fname=TRIM(datadir)//TRIM(runid)//'_restart.bin'
   END FUNCTION restartname
! =====================================================================
   SUBROUTINE secthead(u,name,nrec)
! start a section: its name and the number of records that follow
   INTEGER, INTENT(IN) :: u,nrec
   CHARACTER(len=*), INTENT(IN) :: name
   CHARACTER(len=16) :: tag
   tag=name
   WRITE(u) tag, nrec
   END SUBROUTINE secthead
! =====================================================================
   SUBROUTINE writerestart(time,iout,nblock,nforce)
!
! write the restart file. Must be called at an output time, just after
! statusout (everyone is synchronised and the integer clock is zero).
!
   DOUBLE PRECISION, INTENT(IN) :: time
   INTEGER, INTENT(IN) :: iout
   INTEGER(KIND=8), INTENT(IN) :: nblock,nforce
   CHARACTER(len=512) :: fname
   INTEGER :: nnbm
!
   fname=restartname()
   OPEN(9,file=TRIM(fname)//'.tmp',form='unformatted',status='replace')
!
   nnbm=0
   IF (useac) nnbm=nnbmax
   CALL secthead(9,'HEADER',1)
   WRITE(9) rstversion, nalloc, n, ntot, nbin, kfmax, nnbm
!
   CALL secthead(9,'RUN',1)
   WRITE(9) time, time0, iout, nblock, nforce, comshift, covshift, icformat
!
   CALL secthead(9,'COUNTERS',1)
   WRITE(9) nclamp, nclosed, nopened, nfewformed, nfewopened, nfrozen, &
&    nfewsteps, nmerged, emerged, ntreebuilt, nregular, nirregular
!
   CALL secthead(9,'ENERGY',1)
   WRITE(9) etot0, eprev, demax, nmergedsnap
!
   CALL secthead(9,'PARTICLES',17)
   WRITE(9) r(1:3,1:n)
   WRITE(9) v(1:3,1:n)
   WRITE(9) f0(1:3,1:n)
   WRITE(9) f1(1:3,1:n)
   WRITE(9) f2(1:3,1:n)
   WRITE(9) f3(1:3,1:n)
   WRITE(9) f2e(1:3,1:n)
   WRITE(9) f3e(1:3,1:n)
   WRITE(9) m(1:n)
   WRITE(9) tpair2(1:n)
   WRITE(9) rmin2(1:n)
   WRITE(9) phi(1:n)
   WRITE(9) id(1:n)
   WRITE(9) itype(1:n)
   WRITE(9) ib(1:n)
   WRITE(9) lev(1:n)
   WRITE(9) dtk(1:n)
!
   IF (nbin > 0) THEN
     CALL secthead(9,'GROUPS',8)
     WRITE(9) bnum(1:nbin)
     WRITE(9) bmode(1:nbin)
     WRITE(9) bmass(1:kfmax,1:nbin)
     WRITE(9) bident(1:kfmax,1:nbin)
     WRITE(9) btype(1:kfmax,1:nbin)
     WRITE(9) brel(1:6,1:nbin)
     WRITE(9) brel2(1:6,1:nbin)
     WRITE(9) bepoch(1:nbin)
     IF (fewbody) THEN
       CALL secthead(9,'FEWBODY',2)
       WRITE(9) bx(1:3,1:kfmax,1:nbin)
       WRITE(9) bv(1:3,1:kfmax,1:nbin)
     END IF
   END IF
!
   IF (useac) THEN
     CALL secthead(9,'NEIGHBOURS',11)
     WRITE(9) rs0, nnbopt
     WRITE(9) fi(1:3,1:n)
     WRITE(9) fid(1:3,1:n)
     WRITE(9) fr(1:3,1:n)
     WRITE(9) frd(1:3,1:n)
     WRITE(9) rs(1:n)
     WRITE(9) nnb(1:n)
     WRITE(9) nbl(1:nnbmax,1:n)
     WRITE(9) levr(1:n)
     WRITE(9) dtreg(1:n)
     WRITE(9) acok(1:n)
   END IF
!
! >>> NEW SECTIONS GO HERE
!
   CALL secthead(9,'END',0)
   CLOSE(9)
   CALL EXECUTE_COMMAND_LINE('mv -f '//TRIM(fname)//'.tmp '//TRIM(fname))
   IF (lverbose) WRITE(6,'(a,a)') '   restart file written: ', TRIM(fname)
!
   END SUBROUTINE writerestart
! =====================================================================
   SUBROUTINE readrestart(time,iout,nblock,nforce)
!
! read the restart file, allocate the arrays and set everything up to
! carry on from the time it was written
!
   DOUBLE PRECISION, INTENT(OUT) :: time
   INTEGER, INTENT(OUT) :: iout
   INTEGER(KIND=8), INTENT(OUT) :: nblock,nforce
   CHARACTER(len=512) :: fname
   CHARACTER(len=16) :: tag
   INTEGER :: ios,nrec,k,iver,nn,nt,nb,kf,nnbm,nopt
   LOGICAL :: gothead,gotpart,gotac
!
   fname=restartname()
   OPEN(9,file=TRIM(fname),form='unformatted',status='old',action='read',iostat=ios)
   IF (ios /= 0) THEN
     WRITE(6,'(a,a)') ' ERROR: restart = .true. but there is no restart file ', &
&      TRIM(fname)
     STOP 1
   END IF
   gothead=.FALSE.
   gotpart=.FALSE.
   gotac=.FALSE.
   nn=0
   nb=0
!
   DO
     READ(9,iostat=ios) tag, nrec
     IF (ios /= 0) THEN
       WRITE(6,'(a,a)') ' ERROR: the restart file is incomplete: ', TRIM(fname)
       STOP 1
     END IF
     IF (TRIM(tag) /= 'HEADER' .AND. TRIM(tag) /= 'END' .AND. .NOT. gothead) THEN
       WRITE(6,'(a)') ' ERROR: the restart file does not start with its header'
       STOP 1
     END IF
     SELECT CASE (TRIM(tag))
     CASE ('END')
       EXIT
     CASE ('HEADER')
       READ(9) iver, nalloc, nn, nt, nb, kf, nnbm
       IF (iver > rstversion .OR. kf /= kfmax) THEN
         WRITE(6,'(a)') ' ERROR: the restart file was written by a different '// &
&          'version of HENDOS'
         STOP 1
       END IF
       IF ((nnbm > 0) .NEQV. useac) THEN
         WRITE(6,'(a)') ' ERROR: useac must be the same as when the restart '// &
&          'file was written'
         STOP 1
       END IF
! the arrays are allocated for the number of stars the run began with
       n=nalloc
       CALL allocarrays
       IF (useac) THEN
         nnbmax=nnbm
         nnbopt=MAX(nnbopt,1)
         CALL acalloc
       END IF
       n=nn
       ntot=nt
       nbin=nb
       gothead=.TRUE.
     CASE ('RUN')
       READ(9) time, time0, iout, nblock, nforce, comshift, covshift, icformat
     CASE ('COUNTERS')
       READ(9) nclamp, nclosed, nopened, nfewformed, nfewopened, nfrozen, &
&        nfewsteps, nmerged, emerged, ntreebuilt, nregular, nirregular
     CASE ('ENERGY')
       READ(9) etot0, eprev, demax, nmergedsnap
     CASE ('PARTICLES')
       READ(9) r(1:3,1:n)
       READ(9) v(1:3,1:n)
       READ(9) f0(1:3,1:n)
       READ(9) f1(1:3,1:n)
       READ(9) f2(1:3,1:n)
       READ(9) f3(1:3,1:n)
       READ(9) f2e(1:3,1:n)
       READ(9) f3e(1:3,1:n)
       READ(9) m(1:n)
       READ(9) tpair2(1:n)
       READ(9) rmin2(1:n)
       READ(9) phi(1:n)
       READ(9) id(1:n)
       READ(9) itype(1:n)
       READ(9) ib(1:n)
       READ(9) lev(1:n)
       READ(9) dtk(1:n)
       gotpart=.TRUE.
     CASE ('GROUPS')
       READ(9) bnum(1:nbin)
       READ(9) bmode(1:nbin)
       READ(9) bmass(1:kfmax,1:nbin)
       READ(9) bident(1:kfmax,1:nbin)
       READ(9) btype(1:kfmax,1:nbin)
       READ(9) brel(1:6,1:nbin)
       READ(9) brel2(1:6,1:nbin)
       READ(9) bepoch(1:nbin)
     CASE ('FEWBODY')
       IF (.NOT. fewbody) THEN
         WRITE(6,'(a)') ' ERROR: fewbody must be the same as when the '// &
&          'restart file was written'
         STOP 1
       END IF
       READ(9) bx(1:3,1:kfmax,1:nbin)
       READ(9) bv(1:3,1:kfmax,1:nbin)
     CASE ('NEIGHBOURS')
       READ(9) rs0, nopt
       nnbopt=nopt
       READ(9) fi(1:3,1:n)
       READ(9) fid(1:3,1:n)
       READ(9) fr(1:3,1:n)
       READ(9) frd(1:3,1:n)
       READ(9) rs(1:n)
       READ(9) nnb(1:n)
       READ(9) nbl(1:nnbmax,1:n)
       READ(9) levr(1:n)
       READ(9) dtreg(1:n)
       READ(9) acok(1:n)
       gotac=.TRUE.
     CASE DEFAULT
! a section this version does not know: skip it
       WRITE(6,'(a,a,a)') ' NOTE: section ', TRIM(tag), ' of the restart file '// &
&        'is not known to this version of HENDOS and was skipped'
       DO k=1,nrec
         READ(9)
       END DO
     END SELECT
   END DO
   CLOSE(9)
   IF (.NOT. gotpart .OR. (useac .AND. .NOT. gotac)) THEN
     WRITE(6,'(a,a)') ' ERROR: the restart file is incomplete: ', TRIM(fname)
     STOP 1
   END IF
!
! everyone is synchronised at the time of the file, integer clock zero
   rp(1:3,1:n)=r(1:3,1:n)
   vp(1:3,1:n)=v(1:3,1:n)
   tk(1:n)=0
   IF (useac) THEN
     treg(1:n)=0
     lreg(1:n)=.FALSE.
   END IF
   treevalid=.FALSE.
!
! screen, clocks and the counters of the output intervals
   WRITE(6,'(a)') dline
   WRITE(6,'(a,a)') '  HENDOS - HErmite N-body Dynamics Of Systems     run: ', TRIM(runid)
   WRITE(6,'(a)') dline
   WRITE(6,'(a,a)')        '  restarted from         ', TRIM(fname)
   WRITE(6,'(a,f14.4,a)')  '  time                   ', time, ' Myr'
   WRITE(6,'(a,f14.4,a)')  '  end time               ', tend, ' Myr'
   WRITE(6,'(a,i14)')      '  particles              ', ntot
   WRITE(6,'(a,i14)')      '  OpenMP threads         ', nthreads
   WRITE(6,'(a,es14.3)')   '  energy error so far    ', demax
   WRITE(6,'(a)') dline
   CALL SYSTEM_CLOCK(clock0,crate)
   wlast=0.d0
   ioutstart=iout
   nblock0=nblock
   nforce0=nforce
   nclamp0=nclamp
   nclosed0=nclosed
   nopened0=nopened
   nmerged0=nmerged
   nfewformed0=nfewformed
   nfewopened0=nfewopened
   nregular0=nregular
   nirregular0=nirregular
!
   END SUBROUTINE readrestart
!
   END MODULE nbrestart
