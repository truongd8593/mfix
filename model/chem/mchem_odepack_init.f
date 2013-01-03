!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: MCHEM_ODEPACK_INIT                                     C
!     Purpose: controlling values for ODEAPCK(reference to ODEPACK manual)C
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE MCHEM_ODEPACK_INIT
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE mchem
      USE run
      IMPLICIT NONE
!
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------  
!
!====================================================================
!      CALL DLSODA(FEX, NEQ, ODSPEC, T1, ODE_dt, ITOL, RTOL, ATOL, &
!      ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT)
!1)  NEQ, ODSPEC, T1, ODE_dt, ITASK, ISTATE, RWORK, LRW, IWORK, LIW are
!    set in react.f for CALL_DI=.true. and usrfg.f for CALL_ISAT=.true.
!    
!2)  For simpility, usrs just need to provide the following in this 
!    subroutines
!        ITOL according to the ATOL is a scalar or array
!        RTOL, ATOL tolerance parameters
!         EWT(I) = RTOL*ABS(Y(i)) + ATOL if ITOL = 1
!         EWT(I) = RTOL*ABS(Y(i)) + ATOL(i) if ITOL = 2
!        IOPT  0 indicate no  optional input
!        JT    Jacobian type indicator
!     * If the ITOL and number of species change, the defination in 
!       mchem_mod.f should be changed.
!     * If the usrs want to change calling to ODEPACK, please refer to
!       ODEPACK manual.
!====================================================================
! 
      INTEGER         NL
!
!     The part that user must modify
!     user define the tolerances for ODEPACK according ITOL

      IF(.NOT. SOLID_RO_V) THEN !constant particle density
         IF(.NOT. (allocated(ATOL))) allocate(ATOL(12))
         ITOL = 2
         RTOL = 1.0D-3
         DO NL = 1, 12
            ATOL(NL) = 1.0D-5
         END DO

         ATOL(2) = 1.0d-2
         ATOL(9) = 1.0d-2      
         JT     =1
         IOPT   = 0
      ENDIF
!QX
      IF(SOLID_RO_V) THEN !variable particle density
         IF(.NOT. (allocated(ATOL))) allocate(ATOL(25))
         ITOL = 2
         RTOL = 1.0D-3
         
         DO NL = 1, 25
            ATOL(NL) = 1.0D-5
         END DO

         JT     =2
         
         IOPT   = 1
      ENDIF

      RETURN
      END SUBROUTINE MCHEM_ODEPACK_INIT
