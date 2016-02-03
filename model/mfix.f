!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: MFIX                                                    !
!  Author: M. Syamlal                                 Date: 29-JAN-92  !
!                                                                      !
!  Purpose: The main module in the MFIX program                        !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
!
!> \mainpage Multiphase Flow with Interphase eXchanges
!!
!! MFIX is a general-purpose computer code developed at the National
!! Energy Technology Laboratory, NETL, for describing the hydrodynamics,
!! heat transfer, and chemical reactions in fluid-solid systems.
!!
!! It has been used for describing bubbling and circulating fluidized
!! beds and spouted beds. MFiX calculations give transient data on the
!! three-dimensional distribution of pressure, velocity, temperature,
!! and species mass fractions. MFiX code is based on a generally
!! accepted set of multiphase flow equations. The code is used as a
!! "test-stand" for testing and developing multiphase flow constitutive
!!  equations.
!!
!! \section Notice
!! Neither the United States Government nor any agency thereof, nor any
!! of their employees, makes any warranty, expressed or implied, or
!! assumes any legal liability or responsibility for the accuracy,
!! completeness, or usefulness of any information, apparatus, product,
!! or process disclosed or represents that its use would not infringe
!! privately owned rights.
!!
!! * MFIX is provided without any user support for applications in the
!!   user's immediate organization. It should not be redistributed in
!!   whole or in part.
!!
!! * The use of MFIX is to be acknowledged in any published paper based
!!   on computations using this software by citing the MFIX theory
!!   manual. Some of the submodels are being developed by researchers
!!   outside of NETL. The use of such submodels is to be acknowledged
!!   by citing the appropriate papers of the developers of the submodels.
!!
!! * The authors would appreciate receiving any reports of bugs or other
!!   difficulties with the software, enhancements to the software, and
!!   accounts of practical applications of this software.
!!
!! \section Disclaimer
!! This report was prepared as an account of work sponsored by an agency
!! of the United States Government. Neither the United States Government
!! nor any agency thereof, nor any of their employees, makes any
!! warranty, express or implied, or assumes any legal liability or
!! responsibility for the accuracy, completeness, or usefulness of any
!! information, apparatus, product, or process disclosed, or represents
!! that its use would not infringe privately owned rights. Reference
!! herein to any specific commercial product, process, or service by
!! trade name, trademark, manufacturer, or otherwise does not
!! necessarily constitute or imply its endorsement, recommendation, or
!! favoring by the United States Government or any agency thereof. The
!! views and opinions of authors expressed herein do not necessarily
!! state or reflect those of the United States Government or any
!! agency thereof.

PROGRAM MFIX

   !-----------------------------------------------
   ! Modules
   !-----------------------------------------------
   USE ADJUST_DT, ONLY: ADJUSTDT
   USE INTERACTIVE, ONLY: CHECK_INTERACT_ITER
   USE LEQSOL, ONLY: SOLVER_STATISTICS, REPORT_SOLVER_STATS, MAX_NIT
   USE MAIN, ONLY: SETUP, START, END, NIT_TOTAL, REALLY_FINISH, CMD_LINE_ARGS_COUNT, CMD_LINE_ARGS, ADD_COMMAND_LINE_ARGUMENT, IER
   USE PARAM1, ONLY: UNDEFINED, UNDEFINED_I
   USE RUN, ONLY: NSTEP, AUTO_RESTART, AUTOMATIC_RESTART, ITER_RESTART, TIME, TSTOP, DT, INTERACTIVE_MODE, INTERACTIVE_NITS
   USE STEP, ONLY: PRE_STEP, PRE_ITERATE, DO_ITERATION, POST_ITERATE, POST_STEP, LOG_CONVERGED, LOG_DIVERGED

   IMPLICIT NONE

   !-----------------------------------------------
   ! Local variables
   !-----------------------------------------------
   ! flag indicating convergence status with MUSTIT = 0,1,2 implying
   ! complete convergence, non-covergence and divergence respectively
   INTEGER :: MUSTIT

   INTEGER :: II
   CHARACTER(LEN=80) :: tmp
   INTEGER :: NIT
   LOGICAL :: CALL_POST_ITERATE

   DO II=1, COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(II,tmp)
      CALL ADD_COMMAND_LINE_ARGUMENT(tmp)
   ENDDO

   CALL SETUP

   ! AUTO RESTART LOOP
   DO
      CALL START
      REALLY_FINISH = .FALSE.
      DO
         CALL PRE_STEP
         IF (REALLY_FINISH) EXIT
         ! Advance the solution in time by iteratively solving the equations
         DO
            CALL PRE_ITERATE(NIT,MUSTIT)
            ! Begin iterations
            !-----------------------------------------------------------------
            CALL_POST_ITERATE = .TRUE.
            DO
               MUSTIT = 0
               NIT = NIT + 1
               CALL DO_ITERATION(NIT,MUSTIT)

               !  If not converged continue iterations; else exit subroutine.
1000           CONTINUE

               ! Display residuals
               CALL DISPLAY_RESID (NIT)

               ! Determine course of simulation: converge, non-converge, diverge?
               IF (MUSTIT == 0) THEN
                  IF (DT==UNDEFINED .AND. NIT==1) CYCLE   !Iterations converged
                  CALL LOG_CONVERGED(NIT)
                  IER = 0
                  CALL_POST_ITERATE = .FALSE.
                  EXIT   ! for if mustit =0 (converged)
               ELSEIF (MUSTIT==2 .AND. DT/=UNDEFINED) THEN
                  CALL LOG_DIVERGED(NIT)
                  IER = 1
                  CALL_POST_ITERATE = .FALSE.
                  EXIT  ! for if mustit =2 (diverged)
               ENDIF

               ! not converged (mustit = 1, !=0,2 )
               IF(INTERACTIVE_MODE .AND. INTERACTIVE_NITS/=UNDEFINED_I) THEN
                  CALL CHECK_INTERACT_ITER(MUSTIT)
                  IF(MUSTIT == 1) THEN
                     CYCLE
                  ELSE
                     GOTO 1000
                  ENDIF
               ELSEIF(NIT < MAX_NIT) THEN
                  MUSTIT = 0
                  CYCLE
               ENDIF ! continue iterate
            ENDDO
            ! ----------------------------------------------------------------<<<

            IF (CALL_POST_ITERATE) THEN
               CALL POST_ITERATE(NIT)
            ENDIF

            ! SOF: MFIX will not go the next time step if MAX_NIT is reached,
            ! instead it will decrease the time step. (IER changed from 0 to 1)
            IER = 1
            IF (ADJUSTDT(IER,NIT)) EXIT
         ENDDO
         CALL POST_STEP(NIT)
         IF (REALLY_FINISH) EXIT
      ENDDO

      IF(SOLVER_STATISTICS) CALL REPORT_SOLVER_STATS(NIT_TOTAL, NSTEP)
      IF(AUTO_RESTART.AND.AUTOMATIC_RESTART.AND.ITER_RESTART.LE.10) THEN
         CYCLE
      ELSE
         EXIT
      ENDIF
   ENDDO

   CALL END

   STOP

END PROGRAM MFIX
