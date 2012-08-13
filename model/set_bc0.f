!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SET_BC0                                                 C
!  Purpose: This module does the initial setting of all boundary       C
!           conditions. The user specifications of the boundary        C
!           conditions are checked for veracity in check_data_07       C
!                                                                      C
!  Author: M. Syamlal                                 Date: 29-JAN-92  C
!  Reviewer: P. Nicoletti, W. Rogers, S. Venkatesan   Date: 29-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Check whether the velocity components have the correct     C
!           sign at MASS_INFLOW cells                                  C
!  Author: M. Syamlal                                 Date: 08-JUL-92  C
!  Reviewer: M. Syamlal                               Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: BC_DEFINED, BC_TYPE, BC_DT_0, TIME, BC_Jet_g0,C
!                        BC_K_b, BC_K_t, BC_J_s, BC_J_n, BC_I_w,       C
!                        BC_I_e, BC_PLANE, BC_EP_g, BC_P_g, BC_T_g,    C
!                        BC_T_s,  BC_U_g, BC_V_g, BC_W_g,              C
!                        MMAX, BC_ROP_s, BC_U_s, BC_V_s, BC_W_s        C
!  Variables modified: BC_TIME, BC_V_g, I, J, K, IJK, EP_g, P_g, T_g,  C
!                      T_s, U_g, V_g, W_g, ROP_s, U_s, V_s, W_s,       C
!                      M                                               C
!                                                                      C
!  Local variables: L, IJK1, IJK2, IJK3                                C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_BC0 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE compar     
      USE mpi_utility 
      USE physprop
      USE constant
      USE bc
      USE fldvar
      USE indices
      USE run
      USE funits 
      USE scales 
      USE scalars
      USE boundfunijk 
      USE toleranc
      USE sendrecv
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Local index for boundary condition 
      INTEGER ::  L 
! indices 
      INTEGER :: I, J, K, IJK, M, N 
! Local index for setting U velocity b.c. 
      INTEGER :: IJK1 
! Local index for setting V velocity b.c. 
      INTEGER :: IJK2 
! Local index for setting W velocity b.c. 
      INTEGER :: IJK3 
! Dummy variable for gas pressure 
      DOUBLE PRECISION :: PJ 
! number densities for use in GHD theory only
      DOUBLE PRECISION :: nM, nTOT
! Temporary variable to gather fluid_at
      INTEGER, DIMENSION(:), ALLOCATABLE :: FLAG_G
! Logical function to identify a fluid cell in global coordiantes
      LOGICAL :: FLUID_AT_G
!-----------------------------------------------
! External functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: EOSG, CALC_MW 
!----------------------------------------------
! Include statement functions
!----------------------------------------------
      INCLUDE 'sc_p_g1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'sc_p_g2.inc'
!----------------------------------------------

! Global function to determine FLUID_AT_G
      FLUID_AT_G(IJK)    = FLAG_G(IJK) .EQ. 1

!  Setup for cyclic boundary conditions

! Perform the following computations only on the root processor
! and broadcasted everywhere!
      IF(myPE.EQ.root) THEN

! allocate temp var
         ALLOCATE(FLAG_G(IJKMAX3))
         CALL GATHER(FLAG, FLAG_G, root)
         IJK1 = FUNIJK_GL(IMAX1/2 + 1,JMAX1,KMAX1/2 + 1)

! Exact implementation as in the serial code. In the serial version 
! CYCLE has to be replaced by EXIT to have the same meaning as in 
! the original version
         DO IJK = IJK1, ijkmax3 
            IF (FLUID_AT_G(IJK)) EXIT  
         ENDDO 
      
         PJ = UNDEFINED 
         DO L = 1, DIMENSION_BC 
            IF (BC_DEFINED(L) .AND. BC_TYPE(L)=='P_OUTFLOW') PJ = BC_P_G(L) 
         ENDDO 
      
         IF (CYCLIC .AND. PJ == UNDEFINED) THEN 
            IJK_P_G = IJK 
         ELSE 
            IF (PJ == UNDEFINED .AND. RO_G0 .NE. UNDEFINED) THEN ! added incompressible
               IJK_P_G = IJK 
            ELSE 
               IJK_P_G = UNDEFINED_I 
            ENDIF 
         ENDIF 

      ELSE   ! else branch of if(myPE.eq.root)
         ALLOCATE(FLAG_G(1))
         CALL GATHER(FLAG, FLAG_G, root)
      ENDIF   ! end if/else (myPE.eq.root)

! Deallocate storage
      DEALLOCATE(FLAG_G)

! Broadcast the values
      CALL BCAST(PJ)
      CALL BCAST(IJK_P_G)


! Set the boundary conditions - Defining the field variables at the 
! boundaries according to the user specifications. These are not the
! real values in the wall cells, only initial guesses      
      DO L = 1, DIMENSION_BC 
         IF (BC_DEFINED(L)) THEN 


! setting wall (nsw, fsw, psw) boundary conditions
!----------------------------------------------------------------->>>
            IF (BC_TYPE(L)=='FREE_SLIP_WALL' .OR. &
                BC_TYPE(L)=='NO_SLIP_WALL' .OR. &
                BC_TYPE(L)=='PAR_SLIP_WALL') THEN

               DO K = BC_K_B(L), BC_K_T(L) 
               DO J = BC_J_S(L), BC_J_N(L) 
               DO I = BC_I_W(L), BC_I_E(L) 

                  IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                  IJK = BOUND_FUNIJK(I,J,K) 

                  IF (WALL_AT(IJK)) THEN
! this conditional is probably unnecessary since this section already
! falls under the 'wall' if/else branch for the bc_type check...
                     IF(BC_Tw_g(L) /= UNDEFINED)&
                        T_g(IJK) = BC_Tw_g(L) 
                     IF(SMAX > 0) & 
                        WHERE(BC_Tw_s(L,:SMAX) /= UNDEFINED)&
                          T_s(IJK,:SMAX) = BC_Tw_s(L,:SMAX) 
                     IF(MMAX > 0) & 
                        WHERE(BC_Thetaw_m(L,:MMAX) /= UNDEFINED)&
                          Theta_m(IJK,:MMAX) = BC_Thetaw_m(L,:MMAX) 
                     IF(NMAX(0) > 0) & 
                        WHERE (BC_Xw_G(L,:NMAX(0)) /= UNDEFINED) X_G(IJK,:&
                          NMAX(0)) = BC_Xw_G(L,:NMAX(0))
                     DO M = 1, SMAX 
                        IF(NMAX(M) > 0) & 
                        WHERE (BC_Xw_s(L, M, :NMAX(M)) /= UNDEFINED) X_s(IJK,M,:&
                          NMAX(M)) = BC_Xw_s(L, M, :NMAX(M))
                     ENDDO
                     IF(NScalar > 0) & 
                        WHERE (BC_ScalarW(L,:NScalar) /= UNDEFINED)&
                          Scalar(IJK,:NScalar) = BC_ScalarW(L,:NScalar) 
                  ENDIF   ! end if (wall_at(ijk))

               ENDDO   ! end do (i=bc_i_w(l),bc_i_e(l))
               ENDDO   ! end do (j=bc_j_s(l),bc_j_n(l))
               ENDDO   ! end do (k=bc_k_b(l),bc_k_t(l))
! end setting (FSW, NSW or PSW)               
!-----------------------------------------------------------------<<<

            ELSE   ! else branch of if(FSW, NSW or PSW)


! setting all other 'non-wall' boundary conditions 
!----------------------------------------------------------------->>>


! initializing for time dependent mass inflow
               BC_JET_G(L) = UNDEFINED 
               IF (BC_DT_0(L) /= UNDEFINED) THEN 
                  BC_TIME(L) = TIME + BC_DT_0(L) 
                  BC_OUT_N(L) = 0 
                  BC_MOUT_G(L) = ZERO 
                  BC_VOUT_G(L) = ZERO 
                  IF (SMAX > 0) THEN 
                     BC_MOUT_S(L,:SMAX) = ZERO 
                     BC_VOUT_S(L,:SMAX) = ZERO 
                  ENDIF 
                  BC_JET_G(L) = BC_JET_G0(L) 
                  IF (BC_JET_G(L) /= UNDEFINED) THEN 
                     SELECT CASE (TRIM(BC_PLANE(L)))  
                     CASE ('W')  
                        BC_U_G(L) = BC_JET_G(L) 
                     CASE ('E')  
                        BC_U_G(L) = BC_JET_G(L) 
                     CASE ('S')  
                        BC_V_G(L) = BC_JET_G(L) 
                     CASE ('N')  
                        BC_V_G(L) = BC_JET_G(L) 
                     CASE ('B')  
                        BC_W_G(L) = BC_JET_G(L) 
                     CASE ('T')  
                        BC_W_G(L) = BC_JET_G(L) 
                     END SELECT 
                  ENDIF 
               ELSE 
                  BC_TIME(L) = UNDEFINED 
               ENDIF 


               DO K = BC_K_B(L), BC_K_T(L) 
               DO J = BC_J_S(L), BC_J_N(L) 
               DO I = BC_I_W(L), BC_I_E(L)

                  IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                  IJK = BOUND_FUNIJK(I,J,K) 

! this conditional may still be necessary to avoid cyclic boundaries.
! cyclic boundaries would fall under this branch but are considered 
! 'wall_at'
                  IF (.NOT.WALL_AT(IJK)) THEN 

! setting p_outflow_at or outflow_at boundary conditions
! ---------------------------------------------------------------->>>
                     IF (P_OUTFLOW_AT(IJK) .OR. OUTFLOW_AT(IJK)) THEN

! Why check if the bc quantity for the field variable is defined before
! performing the assignment?: 
! For a new run the field variable will be undefined in these locations
! so any check is irrelevant. However, for a restart run the field
! variable in the boundary cell may have an existing value based on the
! flow field. In this case, the check acts to maintain the existing value
! unless given a user defined BC value. This is unlike the PI, MI and MO 
! treatment where the values of field variables in boundary cells are
! strictly overwritten regardless if the BC is defined. In addition, 
! check_data_07 basically ensures such quantities will be defined for 
! MI or PI if they are required.

                        P_STAR(IJK) = ZERO 
                        P_G(IJK) = SCALE(BC_P_G(L)) 
                        IF (BC_EP_G(L) /= UNDEFINED) EP_G(IJK) = BC_EP_G(L) 
! unless gas flow is entering back into the domain and ep_g in the
! boundary cell is defined, t_g will be set to its adjacent fluid cell
! value in set_bc1. so why do this here? if during the initial start-up,
! the gas has backflow and ep_g is defined then t_g may be undefined in
! the boundary.. otherwise this seems unnecessary?
                        IF (BC_T_G(L) /= UNDEFINED) THEN
                          T_G(IJK) = BC_T_G(L) 
                        ELSE
                          T_G(IJK) = TMIN
                        ENDIF
 
                        IF (NMAX(0) > 0) THEN 
                          WHERE (BC_X_G(L,:NMAX(0)) /= UNDEFINED) &
                            X_G(IJK,:NMAX(0)) = BC_X_G(L,:NMAX(0)) 
                        ENDIF 
  
                        IF (NScalar > 0) THEN 
                          WHERE (BC_Scalar(L,:NScalar) /= UNDEFINED)&
                            Scalar(IJK,:NScalar) = BC_Scalar(L,:NScalar)
                        ENDIF  
   
                        IF (K_Epsilon) THEN 
                          IF (BC_K_Turb_G(L) /= UNDEFINED)&
                            K_Turb_G(IJK) = BC_K_Turb_G(L) 
                          IF (BC_E_Turb_G(L) /= UNDEFINED)&
                            E_Turb_G(IJK) = BC_E_Turb_G(L)
                        ENDIF
 
                        DO M = 1, SMAX 
                          IF (BC_ROP_S(L,M) /= UNDEFINED) &
                            ROP_S(IJK,M) = BC_ROP_S(L,M) 
                          IF(BC_T_S(L,M)/=UNDEFINED) &
                            T_S(IJK,M)=BC_T_S(L,M) 
                          IF (BC_THETA_M(L,M) /= UNDEFINED) &
                            THETA_M(IJK,M) = BC_THETA_M(L,M) 

                          IF (NMAX(M) > 0) THEN 
                            WHERE (BC_X_S(L,M,:NMAX(M)) /= UNDEFINED) &
                              X_S(IJK,M,:NMAX(M)) = BC_X_S(L,M,:NMAX(M))
                          ENDIF 
                        ENDDO 

! for GHD theory to compute mixture BC of velocity and density
                        IF(TRIM(KT_TYPE) == 'GHD') THEN
                          ROP_S(IJK,MMAX) = ZERO
                          nTOT = ZERO
                          THETA_M(IJK,MMAX) = ZERO
                          DO M = 1, SMAX 
                            IF (BC_ROP_S(L,M) /= UNDEFINED) THEN
                              ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) + BC_ROP_S(L,M) 
                              nM = BC_ROP_S(L,M)*6d0/(PI*D_p(IJK,M)**3*RO_S(M))
                              nTOT = nTOT + nM
                              IF (BC_THETA_M(L,M) /= UNDEFINED) THETA_M(IJK,MMAX) = &
                                 THETA_M(IJK,MMAX) + nM*BC_THETA_M(L,M) 
                            ENDIF
                          ENDDO 
                          IF(ROP_S(IJK,MMAX) > ZERO) THETA_M(IJK,MMAX) = THETA_M(IJK,MMAX) / nTOT
                        ENDIF
! end setting P_outflow_at or outflow_at boundary conditions
! ----------------------------------------------------------------<<<


                     ELSE   ! else branch of if(P_outflow_at .or. outflow_at)
                            
! setting for p_inflow, mass_inflow or mass_outflow boundary
! conditions
! ---------------------------------------------------------------->>>

! Unlike PO or O boundary condition branch above, in this branch no
! checks are made to determine if the BC is defined before it is
! assigned to the field variable. However, check_data_07 generally
! ensures such BC quantities will be defined if they are needed for 
! PI and MI boundaries.

                        P_STAR(IJK) = ZERO 
                        P_G(IJK) = SCALE(BC_P_G(L)) 
                        EP_G(IJK) = BC_EP_G(L)
                        T_G(IJK) = BC_T_G(L) 

                        IF (NMAX(0) > 0) THEN 
                          X_G(IJK,:NMAX(0)) = BC_X_G(L,:NMAX(0)) 
                        ENDIF 

                        IF (NScalar > 0) THEN 
                          Scalar(IJK,:NScalar) = BC_Scalar(L,:NScalar) 
                        ENDIF  

                        IF (K_Epsilon) THEN 
                          K_Turb_G(IJK) = BC_K_Turb_G(L) 
                          E_Turb_G(IJK) = BC_E_Turb_G(L)
                        ENDIF

                        DO M = 1, SMAX 
                          ROP_S(IJK,M) = BC_ROP_S(L,M) 
                          T_S(IJK,M) = BC_T_S(L,M) 
                          THETA_M(IJK,M) = BC_THETA_M(L,M) 
                          IF (NMAX(M) > 0) THEN 
                             X_S(IJK,M,:NMAX(M)) = BC_X_S(L,M,:NMAX(M)) 
                          ENDIF 
                        ENDDO 

! Changed to make consistent approach                       
                        IJK1 = IJK 
                        IJK2 = IJK 
                        IJK3 = IJK 
                        SELECT CASE (TRIM(BC_PLANE(L)))  
                          CASE ('W')  
                            IJK1 = BOUND_FUNIJK(IM1(I),J,K) 
                          CASE ('S')  
                            IJK2 = BOUND_FUNIJK(I,JM1(J),K) 
                          CASE ('B')  
                            IJK3 = BOUND_FUNIJK(I,J,KM1(K)) 
                        END SELECT 

! When the boundary plane is W, S, or B the velocity components
! need to be set for the boundary cell and adjacent fluid cell.
                        U_G(IJK) = BC_U_G(L) 
                        V_G(IJK) = BC_V_G(L) 
                        W_G(IJK) = BC_W_G(L) 
                        U_G(IJK1) = BC_U_G(L) 
                        V_G(IJK2) = BC_V_G(L) 
                        W_G(IJK3) = BC_W_G(L) 
 
                        IF (MMAX > 0) THEN 
                           U_S(IJK,:SMAX) = BC_U_S(L,:SMAX) 
                           V_S(IJK,:SMAX) = BC_V_S(L,:SMAX) 
                           W_S(IJK,:SMAX) = BC_W_S(L,:SMAX) 
                           U_S(IJK1,:SMAX) = BC_U_S(L,:SMAX) 
                           V_S(IJK2,:SMAX) = BC_V_S(L,:SMAX) 
                           W_S(IJK3,:SMAX) = BC_W_S(L,:SMAX) 
                        ENDIF 

! for GHD theory to compute mixture BC of velocity and density
                        IF(TRIM(KT_TYPE) == 'GHD') THEN
                           ROP_S(IJK,MMAX) = ZERO
                           nTOT = ZERO 
                           THETA_M(IJK,MMAX) = ZERO 
                           U_S(IJK,MMAX) =  ZERO 
                           V_S(IJK,MMAX) =  ZERO 
                           W_S(IJK,MMAX) =  ZERO 
                           U_S(IJK1,MMAX) =  ZERO 
                           V_S(IJK2,MMAX) =  ZERO 
                           W_S(IJK3,MMAX) =  ZERO 
                           DO M = 1, SMAX 
                              ROP_S(IJK,MMAX) = ROP_S(IJK,MMAX) + BC_ROP_S(L,M) 
                              nM = BC_ROP_S(L,M)*6d0/(PI*D_p(IJK,M)**3*RO_S(M))
                              nTOT = nTOT + nM
                              THETA_M(IJK,MMAX) = THETA_M(IJK,MMAX) + nM*BC_THETA_M(L,M) 
                              U_S(IJK,MMAX) = U_S(IJK,MMAX) + BC_ROP_S(L,M)*BC_U_S(L,M) 
                              V_S(IJK,MMAX) = V_S(IJK,MMAX) + BC_ROP_S(L,M)*BC_V_S(L,M) 
                              W_S(IJK,MMAX) = W_S(IJK,MMAX) + BC_ROP_S(L,M)*BC_W_S(L,M) 
                              SELECT CASE (TRIM(BC_PLANE(L)))    
                                CASE ('W')  
                                  U_S(IJK1,MMAX) = U_S(IJK1,MMAX) + BC_ROP_S(L,M)*BC_U_S(L,M)     
                                CASE ('S')  
                                  V_S(IJK2,MMAX) = V_S(IJK2,MMAX) + BC_ROP_S(L,M)*BC_V_S(L,M)     
                                CASE ('B')  
                                  W_S(IJK3,MMAX) = W_S(IJK3,MMAX) + BC_ROP_S(L,M)*BC_W_S(L,M) 
                              END SELECT 
                           ENDDO 
                           IF(ROP_S(IJK,MMAX) > ZERO) THEN
                              THETA_M(IJK,MMAX) = THETA_M(IJK,MMAX) / nTOT 
                              U_S(IJK,MMAX) = U_S(IJK,MMAX) / ROP_S(IJK,MMAX)
                              V_S(IJK,MMAX) = V_S(IJK,MMAX) / ROP_S(IJK,MMAX)
                              W_S(IJK,MMAX) = W_S(IJK,MMAX) / ROP_S(IJK,MMAX) 
                              SELECT CASE (TRIM(BC_PLANE(L)))   
                                CASE ('W')   
                                  U_S(IJK1,MMAX) = U_S(IJK1,MMAX) / ROP_S(IJK,MMAX)   
                                CASE ('S')   
                                  V_S(IJK2,MMAX) = V_S(IJK2,MMAX) / ROP_S(IJK,MMAX)  
                                CASE ('B')   
                                  W_S(IJK3,MMAX) = W_S(IJK3,MMAX) / ROP_S(IJK,MMAX) 
                              END SELECT 
                           ENDIF
                        ENDIF  ! end if (trim(kt_type) ='ghd')


                        IF (.NOT.MASS_OUTFLOW_AT(IJK)) THEN 
! this is addressed in later call to set_mw_mix_g, so why do it?
                           IF (MW_AVG == UNDEFINED) MW_MIX_G(IJK) = &
                             CALC_MW(X_G,DIMENSION_3,IJK,NMAX(0),MW_G)
! this is addressed in later call to set_ro_g, so why do it?
                          IF (RO_G0 == UNDEFINED) RO_G(IJK) = &
                             EOSG(MW_MIX_G(IJK),P_G(IJK),T_G(IJK)) 
! this is also addressed in later call to set_ro_g, so why do it?
                           ROP_G(IJK) = EP_G(IJK)*RO_G(IJK) 
                        ENDIF 


                     ENDIF     ! end if/else branch of if(P_outflow_at or outflow_at)
! end setting for p_inflow, mass_inflow or mass_outflow boundary
! conditions
! ----------------------------------------------------------------<<<                     
                  ENDIF       ! end if (not wall_at)
     
               ENDDO   ! do i
               ENDDO   ! do j
               ENDDO   ! do k

            ENDIF      ! end if/else branch if(FSW, NSW or PSW)
! end setting all boundary conditions             
!-----------------------------------------------------------------<<<

         ENDIF         ! if (bc_defined)
      ENDDO            ! do dimension_bc

! FIX AEOLUS Sofiane's bug fix to make T_g nonzero in k=0,1 ghost layers
! when k-decomposition employed
      call send_recv(T_G,2)
      call send_recv(P_G,2)
      call send_recv(X_G,2)

      RETURN  
      END SUBROUTINE SET_BC0 


