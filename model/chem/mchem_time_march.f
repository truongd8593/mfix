!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: MCHEM_TIME_MARCH                                       C
!     Purpose: Called in time_march.f to do rxns calcs                    C
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE MCHEM_TIME_MARCH
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE run
      USE physprop
      USE fldvar
      USE rxns      
      USE mpi_utility   
      USE toleranc
      USE mchem
      IMPLICIT NONE
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
      INTEGER          IJK
!
!                      ODEs solved in ISAT or DI                             
      DOUBLE PRECISION Spec0(NSpec)
! 
!                      Mass transfer between phases
      DOUBLE PRECISION R_tmp(0:MMAX,0:MMAX)
!
!                      Time step for ODEs
      DOUBLE PRECISION ODE_dt
!
!                      Loop indices
      INTEGER          NL, NM, NLM, NODE  
!QX
      DOUBLE PRECISION sum_ep_s
      INTEGER          iii, count_istate_n1
!
      INCLUDE 'ep_s1.inc' 
      INCLUDE 'function.inc'
      INCLUDE 'ep_s2.inc'
!QX
      IF(.NOT. SOLID_RO_V) THEN ! constant particle density
         !
!     Time step for Integration
!
       IF (CALL_DI .AND. (ISATdt .EQ. UNDEFINED)) THEN     
          ODE_dt = DT
       ELSE 
          ODE_dt = ISATdt
       END IF
!QX
       RESTART_REACTION = .false.
!QX end
!      
!     Time marching
!
      DO            
         IF (TIME .LE. TIME_isat) EXIT

             DO IJK = ijkstart3, ijkend3 
                 
                 IF (FLUID_AT(IJK)) THEN

!
!     Assigning the values of ODEs for integeration
!     The ODEs have the format:
!         [ Ro_g, T_g, X_g ]
!         [ EP_s, T_s, X_s, d_p for phase 1]
!         [ EP_s, T_s, X_s, d_p for phase 2]
!         ......
!         [ EP_s, T_s, X_s, d_p for phase m]
!                   Gas phase 
!  
                    NODE = 0
                    Spec0(1) = RO_G(IJK)
                    Spec0(2) = T_G(IJK)

                    DO NL = 1, NMAX(0)
                       Spec0(NL+2) = X_G(IJK,NL)
                    END DO
                    NODE = 2 + NMAX(0)
!
!                   Solid phase 
!  
                     DO NL = 1, MMAX
                       Spec0(NODE+1) = EP_S(IJK,NL)  
                       Spec0(NODE+2) = T_S(IJK,NL)
                       NODE = NODE + 2
                       DO NM = 1, NMAX(NL)
                          Spec0(NODE+NM) = X_S(IJK,NL,NM)
                       END DO
                       NODE = NODE + NMAX(NL) + 1
                       Spec0(NODE) = D_P(IJK,NL) 
                    END DO   

!no rxn in freeboard
                  IF ((1.0D0-EP_g(IJK)) .LE. ZERO_EP_S) THEN
                      Sum_R_g(IJK) = 0.0d0
                      Sum_R_s(IJK,1) = 0.0d0
                      R_tmp(0,1) = 0.0d0
                  ELSE
!
!    CALL ISAT or DI calculations
!
                    DO NL = 1, NSpec
                       IF (Spec0(NL) .LE. 1.0D-20) THEN
                          Spec0(NL) = 0.0d0
                       END IF
                    END DO
!
                    CALL REACT(Spec0, R_tmp, IJK, NSpec, ODE_dt)
!
                    DO NL = 1, NSpec
                       IF (Spec0(NL) .LE. 1.0D-20) THEN
                          Spec0(NL) = 0.0d0
                       END IF
                    END DO
!
!    Tranfer data back to mfix
! 
                    NODE = 0
                    T_G(IJK) = Spec0(2)  
                    DO NL = 1, NMAX(0)
                       X_G(IJK,NL) = Spec0(NL+2)
                    END DO 
                    NODE = 2 + NMAX(0)
                    DO NL = 1, MMAX
                       ROP_S(IJK,NL) = Spec0(NODE+1) * RO_S(NL)
                       T_S(IJK, NL) = Spec0(NODE+2)
                       NODE = NODE + 2
                       DO NM = 1, NMAX(NL)
                         X_S(IJK,NL,NM) = Spec0(NODE+NM)
                       END DO
                       NODE = NODE + NMAX(NL) + 1
                       D_P(IJK,NL) = Spec0(NODE)
                    END DO   
                       
                    EP_G(IJK) = ONE - SUM(ROP_S(IJK, 1:MMAX)/RO_S(1:MMAX))

                 End if
!
!     Interface mass transfer
!
                    DO NL = 0, MMAX 
                       DO NM = NL + 1, MMAX 
                          NLM = NL + 1 + (NM - 1)*NM/2 
                          IF (R_TMP(NL,NM) /= UNDEFINED) THEN 
                              R_PHASE(IJK,NLM) = R_TMP(NL,NM) 
                          ELSE IF (R_TMP(NM,NL) /= UNDEFINED) THEN 
                              R_PHASE(IJK,NLM) = -R_TMP(NM,NL) 
                          ENDIF 
                       END DO 
                    END DO 
!
!            
                END IF 
            END DO


         TIME_isat = TIME_isat + ODE_dt
!         TIME_isat = TIME_isat + ISATdt
      END DO

      ENDIF !end of constant particle density
!!
      IF(SOLID_RO_V) THEN !variable particle density

!     Time step for Integration
!
       IF (CALL_DI .AND. (ISATdt .EQ. UNDEFINED)) THEN     
          ODE_dt = DT
       ELSE 
          ODE_dt = ISATdt
       END IF

!QX
       mchem_istate = .false.
       RESTART_REACTION = .false.
       count_istate_n1 = 0
!      
!     Time marching
!
      DO            
         IF (TIME .LE. TIME_isat) EXIT

             DO IJK = ijkstart3, ijkend3 
                 
                 IF (FLUID_AT(IJK)) THEN

!t_s < 300 K , not reaction calcultion
         IF(T_s(IJK,2) .lt. 300 .and. T_g(IJK) .lt. 300) then
            goto 1111
         endif
!!
!
!     Assigning the values of ODEs for integeration
!     The ODEs have the format:
!         [ Ro_g, T_g, X_g ]
!         [ EP_s, T_s, X_s, d_p for phase 1]
!         [ EP_s, T_s, X_s, d_p for phase 2]
!         ......
!         [ EP_s, T_s, X_s, d_p for phase m]
!                   Gas phase 
!  
                    NODE = 0
                    Spec0(1) = RO_G(IJK)
                    Spec0(2) = T_G(IJK)

                    DO NL = 1, NMAX(0)
                       Spec0(NL+2) = X_G(IJK,NL)
                    END DO
                    NODE = 2 + NMAX(0)
!
!                   Solid phase 
!QX
                     DO NL = 1, MMAX
                        SPec0(NODE+1) = RO_SV(IJK,NL)
                        Spec0(NODE+2) = EP_S(IJK,NL)  
                        Spec0(NODE+3) = T_S(IJK,NL)
                        NODE = NODE + 3
                        Spec0(NODE+1) = D_P(IJK,NL) 
                        NODE = NODE + 1

                       DO NM = 1, NMAX(NL)
                          Spec0(NODE+NM) = X_S(IJK,NL,NM)
                       END DO

                       NODE = NODE + NMAX(NL)
                    END DO   
 
!no rxn in freeboard
                  IF ((SUM(X_G(IJK,1:8)) .LE. 0.d0) .and. (EP_S(IJK,2) .LE. ZERO_EP_S)) THEN
                      Sum_R_g(IJK) = 0.0d0
                      Sum_R_s(IJK,1) = 0.0d0
                      Sum_R_s(IJK,2) = 0.0d0
                      R_tmp(0,1) = 0.0d0
                      R_tmp(0,2) = 0.0d0
                      R_tmp(1,2) = 0.0d0
                  ELSE
!
!    CALL ISAT or DI calculations
!
                    DO NL = 1, NSpec
                       IF (Spec0(NL) .LE. 1.0D-20) THEN
                          Spec0(NL) = 0.0d0
                       END IF
                    END DO
!QX
                    DO NL = 0, MMAX 
                       DO NM = NL + 1, MMAX 
                          NLM = NL + 1 + (NM - 1)*NM/2 
                          R_TMP(NL,NM) = ZERO
                          R_TMP(NM,NL)  = ZERO
                       END DO 
                    END DO 

                    CALL REACT(Spec0, R_tmp, IJK, NSpec, ODE_dt)

!QX  reduce time step
                    if(mchem_istate .and. count_istate_n1 .le. 20) then
                       write(*,*) 'redece time step'
                       RESTART_REACTION = .true.
                       DT = 0.9d0 * DT
                       count_istate_n1 = count_istate_n1 + 1
                       return
                    endif
!end
                    DO NL = 1, NSpec
                       IF (Spec0(NL) .LE. 1.0D-20) THEN
                          Spec0(NL) = 0.0d0
                       END IF
                    END DO
!
!
!    Tranfer data back to mfix
! 
                    NODE = 0
                    T_G(IJK) = Spec0(2)  

                    DO NL = 1, NMAX(0)
                       X_G(IJK,NL) = Spec0(NL+2)
                       X_g(IJK,NL) = MIN(ONE,MAX(ZERO,X_g(IJK,NL)))
                    END DO
                    NODE = 2 + NMAX(0)
!QX
                    NODE = NODE + 4 + NMAX(1) ! pass by solid phase one
                    IF(NMAX(1) .ne. 1) then
                       print*,'ERROR in NMAX(1)'
                       STOP
                       ENDIF

                    DO NL = 2, MMAX
                       RO_SV(IJK,NL) = Spec0(NODE+1)
                       ROP_S(IJK,NL) = Spec0(NODE+2) * RO_SV(IJK,NL)
!end
                       T_S(IJK, NL) = Spec0(NODE+3)
                       NODE = NODE + 3

                       D_P(IJK,NL) = Spec0(NODE+1)
                       NODE = NODE + 1

                       DO NM = 1, NMAX(NL)
                         X_S(IJK,NL,NM) = Spec0(NODE+NM)
                         X_S(IJK,NL,NM) = MIN(ONE,MAX(ZERO,X_S(IJK,NL,NM)))
                       END DO

                       NODE = NODE + NMAX(NL)
                    END DO   
!QX                       
                    sum_ep_s = 0.d0
                    do iii=1,MMAX
                       if(ROP_s(IJK,iii) .gt. 0.d0) then
                          sum_ep_s = sum_ep_s + ROP_s(IJK,iii)/RO_SV(IJK,iii)
                       endif
                    enddo

                    EP_G(IJK) = ONE - sum_ep_s                    
                    ROP_G(IJK) = Spec0(1) * EP_G(IJK)
!end
!
                 End if
!
!     Interface mass transfer
!
                    DO NL = 0, MMAX 
                       DO NM = NL + 1, MMAX 
                          NLM = NL + 1 + (NM - 1)*NM/2 
                          IF (R_TMP(NL,NM) /= UNDEFINED) THEN 
                             R_PHASE(IJK,NLM) = R_TMP(NL,NM) 
                          ELSE IF (R_TMP(NM,NL) /= UNDEFINED) THEN 
                             R_PHASE(IJK,NLM) = -R_TMP(NM,NL) 
                          ENDIF
                       END DO
                       
                    END DO
!

1111                continue

                END IF 
            END DO

         TIME_isat = TIME_isat + ODE_dt
!         TIME_isat = TIME_isat + ISATdt
      END DO
     ENDIF !end variable particle density
!
      RETURN  
      END SUBROUTINE MCHEM_TIME_MARCH 
