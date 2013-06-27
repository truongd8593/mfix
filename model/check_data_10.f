!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_DATA_10                                          C
!  Purpose: Check point source specifications.                         C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_DATA_10

      use param
      use run, only: SPECIES_EQ
      use physprop, only: MMAX
      use physprop, only: NMAX

      use run
      use rxns
      use ps
      use compar
      use geometry
      use mpi_utility

      implicit none

      INTEGER :: IJK, I, J, K, M, N

      INTEGER PSV

      CHARACTER*64 eMsg

      DOUBLE PRECISION lSum

      INTEGER :: iErr

      INCLUDE 'function.inc'

      POINT_SOURCE = .FALSE.
      PS_DEFINED = .FALSE.

      iErr = 0

! DETERMINE WHICH BOUNDARY CONDITION INDICES HAVE VALUES
      PSV_LP: do PSV = 1, DIMENSION_PS

         if(PS_X_W(PSV) /= UNDEFINED) PS_DEFINED(PSV) = .TRUE. 
         if(PS_X_E(PSV) /= UNDEFINED) PS_DEFINED(PSV) = .TRUE. 
         if(PS_Y_S(PSV) /= UNDEFINED) PS_DEFINED(PSV) = .TRUE. 
         if(PS_Y_N(PSV) /= UNDEFINED) PS_DEFINED(PSV) = .TRUE. 
         if(PS_Z_B(PSV) /= UNDEFINED) PS_DEFINED(PSV) = .TRUE. 
         if(PS_Z_T(PSV) /= UNDEFINED) PS_DEFINED(PSV) = .TRUE. 
         if(PS_I_W(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE. 
         if(PS_I_E(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE. 
         if(PS_J_S(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE. 
         if(PS_J_N(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE. 
         if(PS_K_B(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE. 
         if(PS_K_T(PSV) /= UNDEFINED_I) PS_DEFINED(PSV) = .TRUE. 

         if(.NOT.PS_DEFINED(PSV)) cycle PSV_LP

! Flag that one or more point sources has been detected.
         POINT_SOURCE = .TRUE.

! West Face.
         if(PS_I_W(PSV)==UNDEFINED_I) then
            if(PS_X_W(PSV)==UNDEFINED) then
               if(myPE == PE_IO) then
                  write(*,1001) PSV, 'PS_X_W and PS_I_W'
                  write(*,1000)
                  write(UNIT_LOG,1001) PSV, 'PS_X_W and PS_I_W'
                  write(UNIT_LOG,1000)
               endif
               iErr = 1001
            else
               CALL CALC_CELL (XMIN, PS_X_W(PSV), DX, IMAX, I)
               PS_I_W(PSV) = I
            endif
         else
            PS_X_W(PSV) = 0.0d0
            do i = 2, PS_I_W(PSV)-1
               PS_X_W(PSV) = PS_X_W(PSV) + DX(I)
            enddo
         endif

! East Face.
         if(PS_I_E(PSV)==UNDEFINED_I) then
            if(PS_X_E(PSV)==UNDEFINED) then
               if(myPE == PE_IO) then
                  write(*,1001) PSV, 'PS_X_E and PS_I_E'
                  write(*,1000)
                  write(UNIT_LOG,1001) PSV, 'PS_X_E and PS_I_E'
                  write(UNIT_LOG,1000)
               endif
               iErr = 1001
            else
               CALL CALC_CELL (XMIN, PS_X_E(PSV), DX, IMAX, I)
               PS_I_E(PSV) = I
            endif
         else
            PS_X_E(PSV) = 0.0d0
            do i = 2, PS_I_W(PSV)
               PS_X_E(PSV) = PS_X_E(PSV) + DX(I)
            enddo
         endif

! South Face.
         if(PS_J_S(PSV)==UNDEFINED_I) then
            if(PS_Y_S(PSV)==UNDEFINED) then
               if(myPE == PE_IO) then
                  write(*,1001) PSV, 'PS_Y_S and PS_J_S'
                  write(*,1000)
                  write(UNIT_LOG,1001) PSV, 'PS_Y_S and PS_J_E'
                  write(UNIT_LOG,1000)
               endif
               iErr = 1001
            else
               CALL CALC_CELL (ZERO, PS_Y_S(PSV), DY, JMAX, J)
               PS_J_S(PSV) = J
            endif
         else
            PS_Y_S(PSV) = 0.0d0
            do j = 2, PS_J_S(PSV)-1
               PS_Y_S(PSV) = PS_Y_S(PSV) + DY(J)
            enddo
         endif


! North Face.
         if(PS_J_N(PSV)==UNDEFINED_I) then
            if(PS_Y_N(PSV)==UNDEFINED) then
               if(myPE == PE_IO) then
                  write(*,1001) PSV, 'PS_Y_N and PS_J_N'
                  write(*,1000)
                  write(UNIT_LOG,1001) PSV, 'PS_Y_N and PS_J_N'
                  write(UNIT_LOG,1000)
               endif
               iErr = 1001
            else
               CALL CALC_CELL (ZERO, PS_Y_N(PSV), DY, JMAX, J)
               PS_J_N(PSV) = J
            endif
         else
            PS_Y_N(PSV) = 0.0d0
            do j = 2, PS_J_N(PSV)
               PS_Y_N(PSV) = PS_Y_N(PSV) + DY(J)
            enddo
         endif


         if(NO_K) THEN
            PS_Z_B(PSV) = ZERO;     PS_K_B(PSV) = 1
            PS_Z_T(PSV) = ZLENGTH;  PS_K_T(PSV) = 1
         else
! Bottom Face.
            if(PS_K_B(PSV)==UNDEFINED_I) then
               if(PS_Z_B(PSV)==UNDEFINED) then
                  if(myPE == PE_IO) then
                     write(*,1001) PSV, 'PS_Z_B and PS_Z_B'
                     write(*,1000)
                     write(UNIT_LOG,1001) PSV, 'PS_K_B and PS_K_B'
                     write(UNIT_LOG,1000)
                  endif
                  iErr = 1001
               else
                  CALL CALC_CELL (ZERO, PS_Z_B(PSV), DZ, KMAX, K)
                  PS_K_B(PSV) = K
               endif
            else
               PS_Z_B(PSV) = 0.0d0
               do k = 2, PS_K_B(PSV)-1
                  PS_Z_B(PSV) = PS_Z_B(PSV) + DZ(K)
               enddo
            endif

! Top Face.
            if(PS_K_T(PSV)==UNDEFINED_I) then
               if(PS_Z_T(PSV)==UNDEFINED) then
                  if(myPE == PE_IO) then
                     write(*,1001) PSV, 'PS_Z_T and PS_Z_T'
                     write(*,1000)
                     write(UNIT_LOG,1001) PSV, 'PS_K_T and PS_K_T'
                     write(UNIT_LOG,1000)
                  endif
                  iErr = 1001
               else
                  CALL CALC_CELL (ZERO, PS_Z_T(PSV), DZ, KMAX, K)
                  PS_K_T(PSV) = K
               endif
            else
               PS_Z_T(PSV) = 0.0d0
               do k = 2, PS_K_T(PSV)
                  PS_Z_T(PSV) = PS_Z_T(PSV) + DZ(K)
               enddo
            endif
         endif

! If there was missing bounds for the point source, the user was
! notified. Skip the following checks as they can result in a floating
! point exception. However, continue reviewing the rest of the point
! sources so that as errors as possible are reported.
         if(iErr /= 0) cycle PSV_LP

! Gas Phase:
!-----------------------------------------------------------------------
! Verify that the mass flow rate and velocity data is given.
         CALL CHECK_PS_FLOW(PSV, 0, PS_MASSFLOW_G(PSV), PS_U_g(PSV),   &
            PS_V_g(PSV), PS_W_g(PSV), iErr)

! Verify that the sum of species mass fractions add to one.
         CALL CHECK_PS_SPECIES(PSV, 0, PS_MASSFLOW_G(PSV), &
            DIM_N_g, PS_X_g(PSV,:),  iErr)

         CALL CHECK_PS_ENERGY(PSV, 0, PS_MASSFLOW_g(PSV), PS_T_g(PSV), &
            iErr)

! Solids Phases:
!-----------------------------------------------------------------------
         M_LP: do M=1,MMAX
! Verify that the mass flow rate and velocity data is given.
            CALL CHECK_PS_FLOW(PSV, M, PS_MASSFLOW_S(PSV,M),           &
               PS_U_s(PSV,M), PS_V_s(PSV,M), PS_W_s(PSV,M), iErr)

! Verify that the sum of species mass fractions add to one.
            CALL CHECK_PS_SPECIES(PSV, M, PS_MASSFLOW_S(PSV,M),        &
               DIM_N_s, PS_X_s(PSV,M,:),  iErr)

! Verify that energy data is provided.
            CALL CHECK_PS_ENERGY(PSV, M, PS_MASSFLOW_s(PSV,M),         &
               PS_T_s(PSV,M), iErr)

         enddo M_LP

      enddo PSV_LP

  500 if(iErr /= 0) then
         if(myPE == PE_IO) write(*,9999)
         CALL mfix_exit(myPE)
      endif

      return

 1000 FORMAT(/' Please correct the mfix.dat file.',/1X,70('*')/) 

 1001 FORMAT(/1X,70('*')/' From: CHECK_DATA_10 -',                     &
         ' Invalid specifications for point source ',I3,'.',/          &
         ' Error: 1001: Neither ',A,' are defined!')

 1002 FORMAT(/1X,70('*')/' From: CHECK_DATA_10',/' Error 1002:',       &
         ' Invalid specifications for point source ',I3,'!',/1x,A,     &
         ' has a negative mass flow rate (PS_MASSFLOW_',A1,' < 0.0).',/&
         ' Point sources can only add mass to a system.')

 1003 FORMAT(/1X,70('*')/' From: CHECK_DATA_10',/' Error 1003:',       &
         ' Invalid specifications for point source ',I3,'!',/          &
         ' The velocity of ',A,' is undefined:'/' PS_U_',A1,&
         ' = PS_V_',A1,' = PS_W_',A1,' = UNDEFINED.')

 9999 FORMAT(2/,' One or more fatal errors found in check_data_10.',&
         ' Calling MFIX_EXIT.',/)

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_PS_FLOW                                          C
!  Purpose: Check point source specifications.                         C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_PS_FLOW(lPSV, lM, mFlow, lU, lV, lW, lErr)

      INTEGER, intent(in) :: lPSV  ! PS Index
      INTEGER, intent(in) :: lM    ! Phase index

      DOUBLE PRECISION, intent(inout) :: mFlow       ! PS_MASSFLOW_x
      DOUBLE PRECISION, intent(inout) :: lU, lV, lW  ! PS_U/V/W_x

      INTEGER, intent(inout) :: lErr

      CHARACTER*1  :: cM
      CHARACTER*16 :: pMsg


      if(lM == 0) then
         pMsg = 'Gas phase'
         cM = 'g'
      else
         pMsg = ''; write(pMsg,"('Solids phase',I2)") lM
         cM = 's'
      endif


! Mass flow is undefined --> Velocity must also be undefined.
!```````````````````````````````````````````````````````````````````````
      if(mFlow == UNDEFINED) then
         if(lU /= UNDEFINED .OR. &
            lV /= UNDEFINED .OR. &
            lW /= UNDEFINED) then
            if(myPE == PE_IO) then
               write(*,1004) lPSV, trim(pMsg)
               write(UNIT_LOG,1004) lPSV, trim(pMsg)
            endif
            lErr = 1004; goto 500
         else
! Zero out the entries.
            mFlow = ZERO
            lU = ZERO
            lV = ZERO
            lW = ZERO
            return
         endif

! Mass flow is zero --> Velocity must also be zero.
!```````````````````````````````````````````````````````````````````````
      elseif(mFlow == ZERO) then

         if(abs(lU+lV+lW) > small_number) then
            if(myPE == PE_IO) then
               write(*,1005) lPSV, trim(pMsg)
               write(UNIT_LOG,1005) lPSV, trim(pMsg)
            endif
            lErr = 1005; goto 500
         endif


! Mass flow is negative --> ERROR
!```````````````````````````````````````````````````````````````````````
      elseif(mFlow < ZERO) then
         if(myPE == PE_IO) then
            write(*,1002) lPSV, trim(pMsg), cM
            write(UNIT_LOG,1002) lPSV, trim(pMsg), cM
         endif
         lErr = 1002; goto 500

! Mass flow is specified:
! Velocity does not have to be defined (no momentum source). If the
! components are UNDEFINED, zero them out.
!```````````````````````````````````````````````````````````````````````
      else
         if(lU == UNDEFINED) lU = ZERO
         if(lV == UNDEFINED) lV = ZERO
         if(lW == UNDEFINED) lW = ZERO

!         if( (lU + lV + lW) < small_number) then
!            if(myPE == PE_IO) then
!               write(*,1003) PSV, trim(pMsg)
!               write(UNIT_LOG,1003) PSV, trim(pMsg)
!            endif
!            lErr = 1003; goto 500
!         endif
      endif

      return  ! No input errors.


! Finish writing error message and return.
  500 if(myPE == PE_IO) then
         write(*,1000)
         write(UNIT_LOG,1000)
      endif

      return     


 1000 FORMAT(/' Please correct the mfix.dat file.',/1X,70('*')/) 

 1002 FORMAT(/1X,70('*')/' From: CHECK_DATA_10 -',                     &
         ' Invalid specifications for point source ',I3,'.',/          &
         ' Error 1002: ',A,' has a negative mass flow rate',/          &
         ' (PS_MASSFLOW_',A1,' < 0.0). Point sources can only add',    &
         ' mass to a system.')

 1003 FORMAT(/1X,70('*')/' From: CHECK_DATA_10 -',                     &
         ' Invalid specifications for point source ',I3,'.',/          &
         ' Error 1003: ',A,' velocity is undefined.')

 1004 FORMAT(/1X,70('*')/' From: CHECK_DATA_10 -',                     &
         ' Invalid specifications for point source ',I3,'.',/          &
         ' Error 1004: ',A,' mass flux undefined but velocity is given.')


 1005 FORMAT(/1X,70('*')/' From: CHECK_DATA_10 -',                     &
         ' Invalid specifications for point source ',I3,'.',/          &
         ' Error 1005: ',A,' mass flux is zero but velocity is given.')


      END SUBROUTINE CHECK_PS_FLOW




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_PS_SPECIES                                       C
!  Purpose: Check point source specifications.                         C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_PS_SPECIES(lPSV, lM, mFlow, lDIM, lX, lErr)


      INTEGER, intent(in) :: lPSV  ! PS Index
      INTEGER, intent(in) :: lM    ! Phase index
      INTEGER, intent(in) :: lDIM  ! Dimension of species array

      DOUBLE PRECISION, intent(inout) :: mFlow       ! PS_MASSFLOW_x
      DOUBLE PRECISION, intent(inout) :: lX(lDIM)    ! PS_X_x

      INTEGER, intent(inout) :: lErr

      LOGICAL, external :: COMPARE


      CHARACTER*1  :: cM
      CHARACTER*16 :: pMsg


! Without mass flow, species composition isn't important.
      if(mFlow == ZERO) then
         lX = ZERO
         lX(1) = ONE
         return
      endif

! Clear out mass fractions that are undefined.
      do N=1,NMAX(lM)
         if(lX(N) == UNDEFINED) lX(N) = ZERO
      enddo
! Make sure data isn't given for an unknown species.
      do N=NMAX(lM)+1,lDIM
         if(lX(N) == UNDEFINED) then
            lX(N) = ZERO
         else
            lErr = 1103; goto 500
         endif
      enddo

! Verify that the mass fractions are non-negative.
      do N=1,NMAX(lM)
         if(lX(N) < ZERO) then
            lErr = 1101; goto 500
         endif
      enddo

! At this point, the species mass fractions - even if the species -
! equations are not solved - should sum to either zero or one.
      lSum = sum(lX)


! If the species equations are solved, the sum must equal one.
      if(SPECIES_EQ(lM)) then
! Verify that the mass fractions sum to one.
         if(.NOT.COMPARE(ONE,lSum)) then
            lErr = 1102; goto 500
         endif

! Verify that the mass fractions sum to one or zero.
      else
         if(.NOT.COMPARE(ONE,lSum) .AND. .NOT.COMPARE(ZERO,lSum)) then
            lErr = 1104; goto 500
         elseif(COMPARE(ZERO,lSum))then
            lX(1) = ONE
         endif
      endif

      return  ! No errors.


! Error Message Handeling:
!-----------------------------------------------------------------------
  500 if(myPE == PE_IO) then

! Set error message entires.
         if(lM == 0) then
            pMsg = 'Gas phase';                                 cM = 'g'
         else
            pMsg = ''; write(pMsg,"('Solids phase',I2)") lM;    cM = 's'
         endif

         write(*,1100) lPSV         ! Message header.
         write(UNIT_LOG,1100) lPSV

         SELECT case (lErr)         ! Error specific messages.
         case (1101)  
            write(*,1101) trim(pMsg)
            write(UNIT_LOG,1101) trim(pMsg)

         case (1102)
            write(*,1102) trim(pMsg)
            write(UNIT_LOG,1102) trim(pMsg)

         case (1103)
            write(*,1103) trim(pMsg), cM
            write(UNIT_LOG,1103) trim(pMsg), cM

         case DEFAULT
            write(*,1198) trim(pMsg)
            write(UNIT_LOG,1198) trim(pMsg)

         end SELECT

         write(*,1199)              ! Message footer.
         write(UNIT_LOG,1199)

      endif

      return     

! Error Messages -------------------------------------------------------

 1101 FORMAT(' Error 1101: ',A,' has a negative species mass fraction.')
 1102 FORMAT(' Error 1102: ',A,' mass fractions do not sum to one.')
 1103 FORMAT(' Error 1103: ',A,' has one or more mass fractions', &
         ' specified ',/' for an unlisted species. Check NMAX_',A1)
 1104 FORMAT(' Error 1104: ',A,'has species mass fractions that do',   &
         ' not sum',/' to either one or zero.')

! Generic Error Header/Footer/Default Messages -------------------------

 1100 FORMAT(/1X,70('*')/' From: CHECK_DATA_10 -',                     &
         ' Invalid specifications for point source ',I3,'.')
 1198 FORMAT(' Error 1198: ',A,' caught an unclassified error.')
 1199 FORMAT(/' Please correct the mfix.dat file.',/1X,70('*')/) 


      END SUBROUTINE CHECK_PS_SPECIES



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_PS_ENERGY                                        C
!  Purpose: Check point source specifications.                         C
!                                                                      C
!  Author: J. Musser                                  Date: 10-JUN-13  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CHECK_PS_ENERGY(lPSV, lM, mFlow, lT, lErr)

      INTEGER, intent(in) :: lPSV  ! PS Index
      INTEGER, intent(in) :: lM    ! Phase index

      DOUBLE PRECISION, intent(inout) :: mFlow   ! PS_MASSFLOW_x
      DOUBLE PRECISION, intent(inout) :: lT      ! PS_T_x

      INTEGER, intent(inout) :: lErr

      CHARACTER*16 :: pMsg

      if(.NOT.ENERGY_EQ) return
      if(mFlow == ZERO) return

! Set error message flags.
      if(lM == 0) then
         pMsg = 'Gas phase'
      else
         pMsg = ''; write(pMsg,"('Solids phase',I2)") lM
      endif

! Verify that a temperature is provided.
      if(lT == UNDEFINED) then
         if(myPE == PE_IO) then
            write(*,1201) lPSV, trim(pMsg)
            write(UNIT_LOG,1201) lPSV, trim(pMsg)
         endif
         lErr = 1201; goto 500

! Verify that a given temperature is physical.
      elseif(lT <= ZERO) then
         if(myPE == PE_IO) then
            write(*,1202) lPSV, trim(pMsg)
            write(UNIT_LOG,1202) lPSV, trim(pMsg)
         endif
         lErr = 1202; goto 500
      endif

      return  ! No errors.

! Finish writing error message and return.
  500 if(myPE == PE_IO) then
         write(*,1200)
         write(UNIT_LOG,1200)
      endif

      return     


 1200 FORMAT(/' Please correct the mfix.dat file.',/1X,70('*')/) 

 1201 FORMAT(/1X,70('*')/' From: CHECK_DATA_10 -',                     &
         ' Invalid specifications for point source ',I3,'.',/          &
         ' Error 1201: ',A,' temperature is undefined.')

 1202 FORMAT(/1X,70('*')/' From: CHECK_DATA_10 -',                     &
         ' Invalid specifications for point source ',I3,'.',/          &
         ' Error 1202: ',A,' temperature is unphysical.')

      END SUBROUTINE CHECK_PS_ENERGY


      END SUBROUTINE CHECK_DATA_10
