!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_MI                                                 !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_MASS_INLET

      USE compar
      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry

      IMPLICIT NONE

!-----------------------------------------------
! Local variables 
!-----------------------------------------------

      INTEGER I, J, K         ! Loop Counter
      INTEGER IP, NP          ! Counter, Index
      INTEGER LS              ! Loop Counter start value

      REAL RAND_x, RAND_y, RAND_z
!-----------------------------------------------

      LS = 1 

      DO IP = 1, PI_COUNT

!       Check to see if MAX_PIS has been exceeded, if so, STOP
        IF(PIS .GT. MAX_PIS) THEN
          WRITE (UNIT_LOG, 1000)
          CALL MFIX_EXIT(myPE)
        ENDIF 

!       Find the first free space in the particle existance array
        DO NP = LS, MAX_PIS
          IF(.NOT.PEA(NP)) EXIT
        ENDDO

        PEA(NP) = .TRUE.  !Set the flage in the particle existance array
      
        PRINT*,'NEW PARTICLE INDEX: NP',NP      

        PIS = PIS + 1 !Increment the particle in system value by one

!       Set the initial position values
        IF(DES_BC_X_w /= DES_BC_X_e)THEN
          CALL RANDOM_NUMBER(RAND_x)
          DES_POS_OLD(NP,1) = (RAND_x*(DES_BC_X_e-DES_BC_X_w))+DES_BC_X_w
        ELSE
          DES_POS_OLD(NP,1) = DES_BC_X_w
        ENDIF

        PRINT*,'DES_POS_OLD(NP,1):',DES_POS_OLD(NP,1)

        IF(DES_BC_Y_s /= DES_BC_Y_n)THEN
          CALL RANDOM_NUMBER(RAND_y)
          DES_POS_OLD(NP,2) = (RAND_y*(DES_BC_Y_n-DES_BC_Y_s))+DES_BC_Y_s
        ELSE
          DES_POS_OLD(NP,2) = DES_BC_Y_s
        ENDIF

        PRINT*,'DES_POS_OLD(NP,2):',DES_POS_OLD(NP,2)

        IF(DIMN == 3) THEN
          IF(DES_BC_Z_b /= DES_BC_Z_t)THEN
            CALL RANDOM_NUMBER(RAND_z)
            DES_POS_OLD(NP,3) = (RAND_z*(DES_BC_Z_t-DES_BC_Z_b))+DES_BC_Z_b
          ELSE
            DES_POS_OLD(NP,3) = DES_BC_Z_b
          ENDIF
        ENDIF


        DES_POS_NEW(NP,:) = DES_POS_OLD(NP,:)
           
!       Set the initial velocity values
        DES_VEL_OLD(NP,:) = 0
        DES_VEL_NEW(NP,:) = DES_VEL_OLD(NP,:)

!       Set the initial angular velocity values
        OMEGA_OLD(NP,:) = 0

!       Set the particle radius value
        DES_RADIUS(NP) = DES_RADIUS(1)  !<------------------------------------ FIX LATER

!       Set the particle radius value
        RO_Sol(NP) = RO_Sol(1)  !<-------------------------------------- FIX LATER

!       Calculate the new particle's Volume, Mass, OMOI
        PVOL(NP) = PI * DES_RADIUS(NP)**3
        PMASS(NP) = PVOL(NP) * RO_Sol(NP)
        OMOI(NP) = 5 / (2 * PMASS(NP) * DES_RADIUS(NP)**2) 

!       Deterine the fluid cell IJK contianing the new particle
        DO I = IMIN1, IMAX3
          IF((DES_POS_NEW(NP,1) .GE. XE(I-1)) .AND. &
             (DES_POS_NEW(NP,1) .LT. XE(I))) THEN
            PIJK(NP,1) = I
            EXIT
          ENDIF
        ENDDO
        DO J = JMIN1, JMAX3
          IF((DES_POS_NEW(NP,2) .GE. YN(J-1)) .AND. &
             (DES_POS_NEW(NP,2) .LT. YN(J))) THEN
            PIJK(NP,2) = J
            EXIT
          ENDIF
        ENDDO
        IF(DIMN.EQ.2) THEN
          PIJK(NP,3)  = 1
        ELSE
          DO K = KMIN1, KMAX3
            IF((DES_POS_NEW(NP,3) .GT. ZT(K-1)) .AND. &
               (DES_POS_NEW(NP,3) .LE. ZT(K))) THEN 
              PIJK(NP,3) = K
              EXIT
            ENDIF
          ENDDO
        ENDIF


      ENDDO

      CALL NEIGHBOUR

      PRINT*,'From: DES_MASS_INLET - End'
      PRINT*,'Have not died yet!'


 1000 FORMAT(/1X,70('*')//&
         ' From: DES_MASS_INLET -',/&
         ' Message: Maximum number of particles in the system MAX_PIS',/&
         ' has been exceeded.  Increase the value in mfix.dat',//&
         1X,70('*')/)
       
      RETURN
      END SUBROUTINE DES_MASS_INLET
