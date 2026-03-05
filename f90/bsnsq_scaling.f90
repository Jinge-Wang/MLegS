PROGRAM BSNSQ_SCALING
!=======================================================================
! [USAGE]:
! Scaling-test driver for STEP_BOUSSI_ETDAB_CN.
!
! Initialises the simulation from the Q-vortex configuration defined in
! read.input, runs a warm-up Richardson startup step through
! PT_SOLVER%INITIALIZE (timer OFF), then executes exactly 100 timed
! steps of STEP_BOUSSI_ETDAB_CN with the wall-clock profiler enabled.
!
! Timed regions reported by TIMER_REPORT:
!   STEP_ETDAB_CN   - full STEP_BOUSSI_ETDAB_CN per step
!   NONLIN          - complete CALC_BOUSSI_NONLIN per step
!     NONLIN_COMP   - gradient + vector products (pure computation)
!     NONLIN_TOFP   - TOFP space conversions (FFF/PFF -> PPP)
!     NONLIN_PROJ   - PROJECT + rebuild to PFF (projection)
!   LINEAR          - CALC_BOUSSI_ETDAB_SUB (linear ETD update)
!
! The MPI_BARRIER inside each TIMER_START / TIMER_END call guarantees
! all ranks are synchronised; only rank-0 accumulates and prints times.
!
! Intended use (SLURM sbatch, see run_scaling.sh):
!   OMP_NUM_THREADS=1 srun -n <NPROCS> ./bin/bsnsq_scaling_exec
! Run for decreasing NPROCS to obtain a strong-scaling profile.
!=======================================================================
   USE OMP_LIB
   USE MPI
   USE MOD_MISC
   USE MOD_BANDMAT
   USE MOD_EIG
   USE MOD_FD
   USE MOD_LIN_LEGENDRE
   USE MOD_SCALAR3
   USE MOD_FFT
   USE MOD_LAYOUT
   USE MOD_LEGOPS
   USE MOD_BOUSSINESQ
   USE MOD_MARCH
   USE MOD_DIAGNOSTICS
   USE MOD_INIT
   IMPLICIT NONE

   INTEGER,  PARAMETER :: NSTEPS_TIMED = 100   ! timed steps per run

   INTEGER         :: IT, III, NPROCS
   TYPE(SCALAR)    :: PSI_TOT, CHI_TOT, B_PER

   ! ------------------------------------------------------------------
   ! 1.  Environment and grid setup (reads read.input, calls MPI_INIT)
   ! ------------------------------------------------------------------
   CALL SETUP_ENVIRONMENT('noecho')
   CALL SETUP_GRID(FILES%SAVEDIR)

   IF (MPI_RANK .EQ. 0) THEN
      CALL MPI_COMM_SIZE(MPI_COMM_IVP, NPROCS, IERR)
      WRITE(*,'(A)') ''
      WRITE(*,'(A,I6,A)') '>>> BSNSQ_SCALING: starting with ', NPROCS, ' MPI procs'
      WRITE(*,'(A,I0,A)') '    Timed steps: ', NSTEPS_TIMED
      WRITE(*,'(A)') ''
   END IF

   ! ------------------------------------------------------------------
   ! 2.  Initial conditions: background Q-vortex, zero density pert.
   ! ------------------------------------------------------------------
   CALL QVORTEX(PSI_TOT, CHI_TOT)
   CALL ALLOCATE(B_PER); B_PER%E = 0.D0

   ! Adjust monitoring mode indices (negative -> wrap-around convention)
   DO III = 1, SIZE(MONITORDATA%MK, 1)
      IF (MONITORDATA%MK(III,1) .LT. 0) &
         MONITORDATA%MK(III,:) = -MONITORDATA%MK(III,:)
      IF (MONITORDATA%MK(III,2) .LT. 0) &
         MONITORDATA%MK(III,2) = 2*NXCHOP - 1 + MONITORDATA%MK(III,2)
   END DO

   ! ------------------------------------------------------------------
   ! 3.  Warm-up / Richardson startup — timer stays OFF
   !     PT_SOLVER%INITIALIZE performs one FE startup step and sets up
   !     the ETD operators; we do not want those costs in the report.
   ! ------------------------------------------------------------------
   TIMER_ENABLED = .FALSE.
   CALL PT_SOLVER%INITIALIZE(PSI_TOT, CHI_TOT, B_PER)

   ! Suppress file I/O during the timed loop
   PT_SOLVER%SAVE_VARIABLE = .FALSE.
   FILES%N = 1

   ! ------------------------------------------------------------------
   ! 4.  Enable timer and zero accumulators
   ! ------------------------------------------------------------------
   TIMER_ENABLED = .TRUE.
   CALL TIMER_RESET()

   ! ------------------------------------------------------------------
   ! 5.  Timed loop — call STEP_BOUSSI_ETDAB_CN directly (bypasses the
   !     diagnostic / output overhead inside PT_SOLVER%TIME_STEPPING)
   ! ------------------------------------------------------------------
   DO IT = 1, NSTEPS_TIMED
      CALL STEP_BOUSSI_ETDAB_CN(PT_SOLVER, PSI_TOT, CHI_TOT, B_PER)
   END DO

   ! ------------------------------------------------------------------
   ! 6.  Print timer report
   ! ------------------------------------------------------------------
   CALL TIMER_REPORT()

   ! ------------------------------------------------------------------
   ! 7.  Finalise (deallocates everything, calls MPI_FINALIZE)
   ! ------------------------------------------------------------------
   CALL PT_SOLVER%FINALIZE(PSI_TOT, CHI_TOT, B_PER)

END PROGRAM BSNSQ_SCALING
