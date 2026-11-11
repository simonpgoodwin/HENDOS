! ==============================================================
   MODULE nbana
!
! HENDOS - HErmite N-body Dynamics Of Systems
! IN-RUN ANALYSIS
!
! Each analysis has an on/off switch in the &analysis section of the
! parameter file and writes its own file  <runid>_ana_<name>.dat
! (in the data directory) with a metadata header (same style as the
! snapshots) followed by one line per analysis time. Analyses are done
! every tana Myr, which is a multiple of tout, so all stars are
! synchronised.
! On a restart (initial conditions time > 0) existing files are
! appended to rather than overwritten.
!
! ana_do is called just after statusout, so the diagnostics in nbout
! (centre of mass, 10/50/90% Lagrangian radii, sigma_3D) are already
! up to date and can be used without recomputing them.
! Energies (direct-N, expensive) are only computed for the screen
! output and snapshots, not by analyses.
!
! checkrunid (also here, as it must know every output file name) stops
! a new run from overwriting the files of an earlier one.
!
! TO ADD A NEW ANALYSIS <name>:
!   1. add a LOGICAL switch <name> to the 'analysis switches' in nbvars,
!      to the &analysis namelist and defaults in readparams, and to
!      params.nml
!   2. give it a unit number below
!   3. open it in ana_init with  CALL ana_open(unit,'<name>',columns,units)
!   4. write SUBROUTINE ana_<name>(time) and call it from ana_do
!   5. close it in ana_end
!   6. add its file to the check in checkrunid (bottom of this file)
!
   USE nbvars
   USE nbout
   IMPLICIT NONE
! unit numbers
   INTEGER, PARAMETER :: ulagr=11      ! Lagrangian radii
! mass fractions for the Lagrangian radii analysis
   INTEGER, PARAMETER :: nlagr=9
   DOUBLE PRECISION, PARAMETER :: flagr(1:nlagr)=(/0.01d0, 0.02d0, 0.05d0, &
&    0.1d0, 0.2d0, 0.3d0, 0.5d0, 0.7d0, 0.9d0/)
!
   CONTAINS
! =====================================================================
   FUNCTION ananame(name) RESULT(fname)
! analysis file name  <runid>_ana_<name>.dat  (in the data directory,
! which is included in the name returned)
   CHARACTER(len=*), INTENT(IN) :: name
   CHARACTER(len=512) :: fname
   fname=TRIM(datadir)//TRIM(runid)//'_ana_'//name//'.dat'
   END FUNCTION ananame
! =====================================================================
   SUBROUTINE ana_open(u,name,columns,units)
!
! open an analysis file and, if it is new, write its header
!
   INTEGER, INTENT(IN) :: u
   CHARACTER(len=*), INTENT(IN) :: name,columns,units
   CHARACTER(len=512) :: fname
   LOGICAL :: ex
!
   fname=ananame(name)
   INQUIRE(file=TRIM(fname),exist=ex)
   IF (time0 > 0.d0 .AND. ex) THEN
! restart - carry on with the existing file
     OPEN(u,file=TRIM(fname),status='old',position='append')
   ELSE
     OPEN(u,file=TRIM(fname),status='replace')
     CALL hdrs(u,'format','NBODY-ANALYSIS 1')
     CALL hdrs(u,'analysis',name)
     CALL hdrcommon(u)
     CALL hdrd(u,'analysis_interval_Myr',tana)
     CALL hdrs(u,'columns',columns)
     CALL hdrs(u,'column_units',units)
     WRITE(u,'(a)') '# end_header'
   END IF
   IF (lverbose) WRITE(6,'(a,a)') '   analysis file: ', TRIM(fname)
!
   END SUBROUTINE ana_open
! =====================================================================
   SUBROUTINE ana_init
!
! open the files for the analyses that are switched on
!
   IF (lagrangian) CALL ana_open(ulagr,'lagrangian', &
&    'time r1% r2% r5% r10% r20% r30% r50% r70% r90% sigma3d', &
&    'Myr pc pc pc pc pc pc pc pc pc km/s')
!
   END SUBROUTINE ana_init
! =====================================================================
   SUBROUTINE ana_do(time)
!
! do all the analyses that are switched on
!
   DOUBLE PRECISION, INTENT(IN) :: time
!
   IF (lagrangian) CALL ana_lagrangian(time)
!
   END SUBROUTINE ana_do
! =====================================================================
   SUBROUTINE ana_end
!
   IF (lagrangian) CLOSE(ulagr)
!
   END SUBROUTINE ana_end
! =====================================================================
   SUBROUTINE ana_lagrangian(time)
!
! Lagrangian radii about the centre of mass and 3D velocity dispersion
!
   DOUBLE PRECISION, INTENT(IN) :: time
   DOUBLE PRECISION :: rl(1:nlagr)
!
   CALL lagradii(nlagr,flagr,rl)
   WRITE(ulagr,'(f16.6,9es15.6,es15.6)') time, rl, sig3d*vunit
   FLUSH(ulagr)
!
   END SUBROUTINE ana_lagrangian
! =====================================================================
   SUBROUTINE checkrunid(nsnap,nana)
!
! for a new run (starting at t=0) stop if output files with this runid
! already exist, so an earlier run is never overwritten by accident
!
   INTEGER, INTENT(IN) :: nsnap,nana
   LOGICAL :: ex,clash
   CHARACTER(len=512) :: fname
!
   clash=.FALSE.
   IF (nsnap > 0) THEN
     fname=snapname(0.d0)
     INQUIRE(file=TRIM(fname),exist=ex)
     IF (ex .AND. TRIM(fname) /= TRIM(icfile)) clash=.TRUE.
   END IF
   IF (nana > 0) THEN
     IF (lagrangian) THEN
       INQUIRE(file=TRIM(ananame('lagrangian')),exist=ex)
       IF (ex) clash=.TRUE.
     END IF
   END IF
   IF (clash) THEN
     WRITE(6,'(a)') ' ERROR: output files for runid '''//TRIM(runid)// &
&     ''' already exist here.'
     WRITE(6,'(a)') '        Choose a new runid in the parameter file, or '// &
&     'move or delete the old files.'
     STOP
   END IF
!
   END SUBROUTINE checkrunid
!
   END MODULE nbana
