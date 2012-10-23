!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP                                           C
!  Purpose: Calculate the indicated physical properties that vary      C
!           with time if directed to do so by the corresponding flag   C
!                                                                      C
!  Author: M. Syamlal                                 Date: 17-JUL-92  C
!  Reviewer: P. Nicoletti                             Date: 11-DEC-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Mods for MFIX 2.0 (old name CALC_PHYSPROP)                 C
!  Author: M. Syamlal                                 Date: 23-APR-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number: 2                                                  C
!  Purpose: allow SI                                                   C
!  Author: S. Dartevelle                              Date: 01-Jul-02  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!    Perry, R.H., and C.H. Chilton, Chemical Engineer's Handbook, 5th  C
!      edition, McGraw-Hill Kogakusha, Tokyo, 1973.                    C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!  Local variables: None                                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PHYSICAL_PROP(DENSITY, PSIZE, SP_HEAT, IER) 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE physprop
      USE geometry
      USE indices
      USE run
      USE toleranc 
      USE constant
      USE scalars
      USE compar 
      USE funits 
      USE usr         
      USE mpi_utility
      USE discretelement
      USE cutcell 
      IMPLICIT NONE
!-----------------------------------------------
! Global parameters
!-----------------------------------------------
! gas constant
      DOUBLE PRECISION, PARAMETER :: RGAS = 1.987207D0  !cal/mol.K
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
! Flags to tell whether to calculate or not
      LOGICAL, INTENT(IN) :: DENSITY(0:DIMENSION_M), &
                             PSIZE(0:DIMENSION_M),&
                             SP_HEAT(0:DIMENSION_M) 
! Error index
      INTEGER, INTENT(INOUT) :: IER
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, IMJK, IJMK, IJKM, I
! solids phase index and species index      
      INTEGER :: M, N
! average molecular weight
      DOUBLE PRECISION :: MW
      LOGICAL :: ALL_IS_ERROR
! alias for gas temperature
      DOUBLE PRECISION :: TGX
!-----------------------------------------------
! External functions
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: EOSG, calc_CpoR,&
                                    calc_ICpoR
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'usrnlst.inc'
      INCLUDE 'cp_fun1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------


!!$omp  parallel do private(IJK, MW, N)  
      DO IJK = IJKSTART3, IJKEND3 
         IF (.NOT.WALL_AT(IJK)) THEN 
            TGX = T_g(IJK)


! Fluid Density
! ---------------------------------------------------------------->>>
            IF (DENSITY(0)) THEN 
               IF (MW_AVG == UNDEFINED) THEN 

! Average molecular weight: Xg1/Mw1 + Xg2/Mw2 + Xg3/Mw3 + ....
!                 MW = CALC_MW(X_g, DIMENSION_3, IJK, NMAX(0), MW_g)

                  IF(.NOT.database_read) call read_database0(IER)
! calculating the average molecular weight of the fluid
                  MW = ZERO
                  N = 1
                  IF (NMAX(0) > 0) THEN
                     MW = MW + SUM(X_G(IJK,:NMAX(0))/MW_G(:NMAX(0)))
                     N = NMAX(0) + 1
                  ENDIF
                  MW = ONE/MAX(MW,OMW_MAX)
                  MW_MIX_G(IJK) = MW
! calculating the fluid density and bulk density
                  RO_G(IJK) = EOSG(MW,P_G(IJK),T_G(IJK))
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
               ELSE
                  RO_G(IJK) = EOSG(MW_AVG,P_G(IJK),T_G(IJK))
                  ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
               ENDIF
       
               IF(RO_G(IJK) < ZERO) THEN
                  IF(CARTESIAN_GRID) THEN 
                     WRITE(*,1001) I_OF(IJK), J_OF(IJK), K_OF(IJK), &
                        RO_G(IJK), P_G(IJK), T_G(IJK), CUT_CELL_AT(IJK), & 
                        SMALL_CELL_AT(IJK), xg_e(I_OF(IJK)), &
                        yg_n(J_of(ijk)), zg_t(k_of(ijk))
                  ELSE
                     WRITE(*,1000) I_OF(IJK), J_OF(IJK), K_OF(IJK), &
                        RO_G(IJK), P_G(IJK), T_G(IJK)
                  ENDIF
                  Neg_RHO_G = .TRUE. !this will reduce dt instead of exiting
               ENDIF
            ENDIF
! ----------------------------------------------------------------<<<



! Fluid specific heat
! ---------------------------------------------------------------->>>
! Constant pressure specific heat of air in cal/g.K
! 1 Cal = 4.183925 J
            IF (SP_HEAT(0) .AND. C_PG0==UNDEFINED) THEN
               IF(.NOT.database_read) call read_database0(IER)

! calculating an average specific heat of the fluid               
!             IF(C(23) == ONE) THEN
                C_PG(IJK) = ZERO
                DO N = 1, NMAX(0)
                   C_PG(IJK) = C_PG(IJK) + X_g(IJK,N) *&
                      calc_CpoR(T_G(IJK),Thigh_g(N),Tlow_g(N),Tcom_g(N),&
                      Ahigh_g(1,N),Alow_g(1,N)) * RGAS / MW_g(N) 
                ENDDO
!             ELSE
!                C_pg(IJK) =  X_g(IJK,CH4)*CPCH4(TGX)  &
!                   + X_g(IJK,CO2)*CPCO2(TGX)  &
!                   + X_g(IJK,O2)*CPO2(TGX) + X_g(IJK,H2O)*CPH2O(TGX) &
!                   + X_g(IJK,N2)*CPN2(TGX) 
!             ENDIF
      
              IF (UNITS == 'SI') C_PG(IJK) = 4183.925D0*C_PG(IJK)    !in J/kg K
            ENDIF   ! end if (sp_heat(0) .and. c_pg0=undefined)

         ENDIF   ! end if (.not.wall_at(ijk))
      ENDDO   ! end do (ijk=ijkstart3,ijkend3)
! ----------------------------------------------------------------<<<


! In case of negative density force exit from the physical property
! calculation routine and reduce the time step
      CALL GLOBAL_ALL_OR(Neg_RHO_G, ALL_IS_ERROR)
      IF(ALL_IS_ERROR) THEN
         Neg_RHO_G = .TRUE.
         RETURN
      ENDIF


! if using a discrete element model then quantities related to the
! continuum solids phase do not need to be calculated. note that the
! code could be changed so that the logical flags directing this
! routine are appropriately set for the discrete element case
      IF(.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID) THEN

!!$omp  parallel do private(IJK)           
         DO IJK = IJKSTART3, IJKEND3 
            DO M=1,MMAX             
               IF (.NOT.WALL_AT(IJK)) THEN

! rong : diameter
! ---------------------------------------------------------------->>>
                  IF (CALL_DQMOM) THEN            
                     IF(PSIZE(M)) THEN 
                        N=phase4scalar(M) ! can N /= M ? (sof)
                        IF(EP_s(IJK,N)>small_number) D_p(IJK,M)= Scalar(IJK,N)
                     ENDIF   ! end if (psize(m))
                  ENDIF   ! end if (call_dqmom)
! ----------------------------------------------------------------<<<

! Solids specific heat
! ---------------------------------------------------------------->>>
! Specific heat of solids (Coal = 0.3 cal/g.K)
!    Perry & Chilton(1973) -- Table 3-201 on page 3-136
! Specific heat of solids (Ash =  0.310713 cal/g.K)
!    Dobran et al., 1991
                  IF (SP_HEAT(M) .AND. C_PS0==UNDEFINED) THEN
                     IF(.NOT.database_read) call read_database0(IER)
                     C_PS(IJK, M) = ZERO
                     DO N = 1, NMAX(M)
                        C_PS(IJK, M) = C_PS(IJK,M) + X_s(IJK,M,N) * &
                           calc_CpoR(T_S(IJK,M),Thigh_s(M,N),&
                           Tlow_s(M,N),Tcom_s(M,N),Ahigh_s(1,M,N),&
                           Alow_s(1,M,N)) * RGAS / MW_s(M,N) 
                     ENDDO
!                     C_PS(IJK,M) = 0.310713d0
                     IF (UNITS == 'SI') C_PS(IJK,M) = 4183.925D0*C_PS(IJK,M)    !in J/kg K
                  ENDIF   ! end if (sp_heat(m) .and. c_ps0==undefined)
! ----------------------------------------------------------------<<<

               ENDIF   ! end if (.not.wall_at(ijk))
               
            ENDDO   ! end do (m=1,mmax)
         ENDDO   ! end do (ijk=ijkstart3,ijkend3)

      ENDIF   ! end if (.not.discrete_element)
              
      RETURN 

 1000 FORMAT(1X,'Message from: PHYSICAL_PROP',/& 
            'WARNING: Gas density negative in this cell: ', /&
            'I = ',I4,2X,' J = ',I4,2X,' K = ',I4, /&
            'Values of variables: ','RO_g = ', G12.5, 2X, &
            'P_g = ', G12.5, 2X, 'T_g = ', G12.5, /&
            'Suggestion: Lower UR_FAC(1) in mfix.dat')

 1001 FORMAT(1X,'Message from: PHYSICAL_PROP',/& 
            'WARNING: Gas density negative in this cell: ', /&
            'I = ',I4,2X,' J = ',I4,2X,' K = ',I4, /&
            'Values of variables: ','RO_g = ', G12.5, 2X, &
            'P_g = ', G12.5, 2X, 'T_g = ', G12.5, /&
            'CUT CELL, SMALL CELL ?', 2(2x,L1), /&
            'East, North, and Top coodinate = ', 3(2x, g17.8), / & 
            'Suggestion: Lower UR_FAC(1) in mfix.dat')

      END SUBROUTINE PHYSICAL_PROP 
