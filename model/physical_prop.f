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
      SUBROUTINE PHYSICAL_PROP(IER, LEVEL)

      use compar
      use funits 
      use geometry
      use indices
      use mpi_utility
      use param1 
      use physprop

      use coeff, only: DENSITY  ! Density
      use coeff, only: SP_HEAT  ! Specific heat 
      use coeff, only: PSIZE    ! Particle diameter

      implicit none

! Dummy arguments
!-----------------------------------------------------------------------
! Global error Flag.
      INTEGER, intent(inout) :: IER
      INTEGER, intent(in) :: LEVEL

! Local variables
!-----------------------------------------------------------------------
! Arrays for storing errors:
! 100 - Negative gas phase density
! 101 - Invalid temperature in calc_CpoR
! 10x - Unclassified
      INTEGER :: Err_l(0:numPEs-1)  ! local
      INTEGER :: Err_g(0:numPEs-1)  ! global

! Initialize error flags.
      Err_l = 0

! Calculate density only. This is invoked several times within iterate,
! making it the most frequently called.
      if(LEVEL == 0) then
         if(DENSITY(0)) CALL PHYSICAL_PROP_ROg
         if(any(DENSITY(1:DIMENSION_M))) CALL PHYSICAL_PROP_ROs

! Calculate everything except density. This is called at the start of
! each iteration.
      elseif(LEVEL == 1) then
         if(SP_HEAT(0)) CALL PHYSICAL_PROP_CPg
         if(any(SP_HEAT(1:DIMENSION_M))) CALL PHYSICAL_PROP_CPs
         if(any(PSIZE(1:DIMENSION_M))) CALL PHYSICAL_PROP_Dp


! Calculate everything. This is invoked via calc_coeff_all as part of
! the initialization (before starting the time march) and at the start
! of each step step thereafter.
      elseif(LEVEL == 2) then
         if(DENSITY(0)) CALL PHYSICAL_PROP_ROg
         if(SP_HEAT(0)) CALL PHYSICAL_PROP_CPg
         if(any(DENSITY(1:DIMENSION_M))) CALL PHYSICAL_PROP_ROs
         if(any(SP_HEAT(1:DIMENSION_M))) CALL PHYSICAL_PROP_CPs
         if(any(PSIZE(1:DIMENSION_M))) CALL PHYSICAL_PROP_Dp
      endif


! In case of negative density force exit from the physical property
! calculation routine and reduce the time step
      CALL global_all_sum(Err_l, Err_g)
      IER = maxval(Err_g)
      if(IER == 0) return


! Error handeling. - Local.
!-----------------------------------------------------------------------
! An invalid temperature was found by calc_CpoR. This is a fatal run-
! time error and forces a call to MFIX_EXIT.
      IF(IER == 101) then
         if(myPE == PE_IO) then
            write(*,2000)
            write(UNIT_LOG,2000)
         endif
         CALL MFIX_EXIT(myPE)
      ENDIF

      return

 2000 FORMAT(/1X,70('*')/' From: PHYSICAL_PROP',/' Fatal Error 2000:', &
         ' calc_CpoR reporetd an invalid temperature.',/,'See Cp.log', &
         ' for details. Calling MFIX_EXIT.',/1X,70('*')/) 

      contains

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_ROg                                       C
!  Purpose: Calculate the gas phase density.                           C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_ROg


! Global variables:
!-----------------------------------------------------------------------
! Gas phase species mass fractions.
      use fldvar, only: X_g
! Gas phase temperature.
      use fldvar, only: T_g
! Gas phase density (compressible).
      use fldvar, only: RO_g
! Gas phase pressure.
      use fldvar, only: P_g
! Gas phase volume fraction.
      use fldvar, only: EP_g
! Gas phase material density.
      use fldvar, only: ROP_g
! Maximum value for molecular weight (divided by one)
      use toleranc, only: OMW_MAX
! Run time flag for generating negative gas density log files
      use run, only: REPORT_NEG_DENSITY

      implicit none

! Local Variables:
!-----------------------------------------------------------------------
! Average molecular weight
      DOUBLE PRECISION :: MW

! Loop indicies
      INTEGER :: IJK   ! Computational cell
      INTEGER :: M     ! Solids phase
      INTEGER :: N     ! Species index

! Equation of State - GAS
      DOUBLE PRECISION, EXTERNAL :: EOSG

! Flag to write log header
      LOGICAL wHeader

      include 'function.inc'

! Initialize:
      wHeader = .TRUE.

! Average molecular weight: Xg1/Mw1 + Xg2/Mw2 + Xg3/Mw3 + ....
      IF(.NOT.database_read) call read_database0(IER)

      IJK_LP: DO IJK = IJKSTART3, IJKEND3 
         IF(WALL_AT(IJK)) cycle IJK_LP
         IF (MW_AVG == UNDEFINED) THEN 
! Calculating the average molecular weight of the fluid.
            MW = SUM(X_G(IJK,:NMAX(0))/MW_G(:NMAX(0)))
            MW = ONE/MAX(MW,OMW_MAX)
            MW_MIX_G(IJK) = MW
! Calculate the fluid density and bulk density
            RO_G(IJK) = EOSG(MW,P_G(IJK),T_G(IJK))
            ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
         ELSE
            RO_G(IJK) = EOSG(MW_AVG,P_G(IJK),T_G(IJK))
            ROP_G(IJK) = RO_G(IJK)*EP_G(IJK)
         ENDIF

         IF(RO_G(IJK) < ZERO) THEN
            Err_l(myPE) = 100
            IF(REPORT_NEG_DENSITY)CALL ROgErr_LOG(IJK, wHeader)
         ENDIF
      ENDDO IJK_LP


      RETURN
      END SUBROUTINE PHYSICAL_PROP_ROg


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_ROs                                       C
!  Purpose: Calculate solids phase (variable) density.                 C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_ROs

      implicit none

      RETURN
      END SUBROUTINE PHYSICAL_PROP_ROs



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_CPg                                       C
!  Purpose: Calculate the gas phase constant pressure specific heat.   C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
! Notes:                                                               C
!  > Unit conversion: 1 cal = 4.183925 J                               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_CPg

! Global Variables:
!-----------------------------------------------------------------------
! Universal gas constant in cal/mol.K
      use constant, only: RGAS => GAS_CONST_cal
! Maximum value for molecular weight (divided by one)
      use toleranc, only: OMW_MAX
! Gas phase species mass fractions.
      use fldvar, only: X_g
! Gas phase temperature.
      use fldvar, only: T_g
! Units: CGS/SI
      use run, only: UNITS

      implicit none

! Local Variables:
!-----------------------------------------------------------------------
! Species specific heat.
      DOUBLE PRECISION :: lCp
! Error flag returned from calc_CpoR
      INTEGER :: lCP_Err
! Average molecular weight
      DOUBLE PRECISION :: MW
! Loop indicies
      INTEGER :: IJK   ! Computational cell
      INTEGER :: M     ! Solids phase
      INTEGER :: N     ! Species index

! Function to evaluate Cp polynomial.
      DOUBLE PRECISION, EXTERNAL :: calc_CpoR

      include 'function.inc'

!-----------------------------------------------------------------------

! Ensure that the database was read. This *should* have been caught by
! check_data_05 but this call remains to prevent an accident.
      IF(.NOT.database_read) CALL read_database0(IER)

      lCP_Err = 0

      IJK_LP: DO IJK = IJKSTART3, IJKEND3 
         IF(WALL_AT(IJK)) CYCLE IJK_LP
! Calculating an average specific heat for the fluid.
         C_PG(IJK) = ZERO
         DO N = 1, NMAX(0)
            lCp = calc_CpoR(T_G(IJK), 0, N, lCP_Err)
            C_PG(IJK) = C_PG(IJK) + X_g(IJK,N) * lCp * RGAS / MW_g(N) 
         ENDDO
      ENDDO IJK_LP

! The database calculation always returns cal/g.K thus the following
! conversion is needed if using SI units. ** Vector operation
      IF(UNITS == 'SI') C_PG = 4.183925d3 * C_PG  !in J/kg K

      RETURN
      END SUBROUTINE PHYSICAL_PROP_CPg




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_CPs                                       C
!  Purpose: Calculate solids phase constant pressure specific heat.    C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_CPs

! Universal gas constant in cal/mol.K
      use constant, only: RGAS => GAS_CONST_cal
! Units: CGS/SI
      use run, only: UNITS
      use toleranc, only: OMW_MAX
      use fldvar, only: T_s
      use fldvar, only: X_s

      implicit none

      DOUBLE PRECISION :: lCp

! Loop indicies
      INTEGER :: IJK   ! Computational cell
      INTEGER :: M     ! Solids phase
      INTEGER :: N     ! Species index

! Local error flag indicating that the Cp is out of range.
      INTEGER :: lCP_Err

! Function to evaluate Cp polynomial.
      DOUBLE PRECISION, EXTERNAL :: calc_CpoR

      include 'function.inc'

! Ensure that the database was read. This *should* have been caught by
! check_data_05 but this call remains to prevent an accident.
      IF(.NOT.database_read) CALL read_database0(IER)

      lCP_Err = 0

      M_LP: DO M=1, MMAX
         IJK_LP: DO IJK = IJKSTART3, IJKEND3 
            IF(WALL_AT(IJK)) CYCLE IJK_LP
! Calculating an average specific heat for the fluid.
            C_PS(IJK, M) = ZERO

            DO N = 1, NMAX(M)
               lCp = calc_CpoR(T_S(IJK,M), M, N, lCP_Err)
               C_PS(IJK,M) = C_PS(IJK,M) + X_s(IJK,M,N) * &
                  (lCp * RGAS / MW_s(M,N))
            ENDDO

         ENDDO IJK_LP
      ENDDO M_LP

! The database calculation always returns cal/g.K thus the following
! conversion is needed if using SI units. ** Vector operation
      IF(UNITS == 'SI') C_PS = 4.183925d3 * C_PS  !in J/kg K

      END SUBROUTINE PHYSICAL_PROP_CPs




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: PHYSICAL_PROP_CPs                                       C
!  Purpose: Calculate solids phase constant pressure specific heat.    C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE PHYSICAL_PROP_Dp

      use run, only: CALL_DQMOM

      use fldvar, only: scalar
      use scalars, only: phase4scalar

      use fldvar, only: ROP_s
      use fldvar, only: RO_s
      use fldvar, only: D_p

      implicit none

! Loop indicies
      INTEGER :: IJK   ! Computational cell
      INTEGER :: M     ! Solids phase

! Map from true index to map.
      INTEGER :: lM

      include 'ep_s1.inc'
      include 'function.inc'
      include 'ep_s2.inc'


      IF(.NOT.CALL_DQMOM) return

      M_LP: DO M=1, MMAX

         lM = phase4scalar(M) ! Map from scalar eq to solids phase

         IF(.NOT.PSIZE(M)) CYCLE M_LP
         IJK_LP: DO IJK = IJKSTART3, IJKEND3 
            IF(WALL_AT(IJK)) CYCLE IJK_LP

            IF(EP_s(IJK,lM) > small_number) D_p(IJK,M)= Scalar(IJK,lM)

         ENDDO IJK_LP
      ENDDO M_LP

      return
      END SUBROUTINE PHYSICAL_PROP_Dp



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: NegROg_LOG                                              C
!  Purpose: Record information about the location and conditions that  C
!           resulted in a negative gas phase density.                  C
!                                                                      C
!  Author: J. Musser                                  Date: 28-JUN-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE ROgErr_LOG(IJK, tHeader)

! Simulation time
      use run, only: TIME
! Gas phase temperature.
      use fldvar, only: T_g
! Gas phase density (compressible).
      use fldvar, only: RO_g
! Gas phase pressure.
      use fldvar, only: P_g
      use cutcell 

      INTEGER, intent(in) :: IJK
      LOGICAL, intent(inout) :: tHeader

      LOGICAL :: lExists
      CHARACTER*32 :: lFile
      INTEGER, parameter :: lUnit = 4868
      LOGICAL, save :: fHeader = .TRUE.


      lFile = ''; 
      if(numPEs > 1) then
         write(lFile,"('ROgErr_',I4.4,'.log')") myPE
      else
         write(lFile,"('ROgErr.log')")
      endif
      inquire(file=trim(lFile),exist=lExists)
      if(lExists) then
         open(lUnit,file=trim(adjustl(lFile)),                         &
            status='old', position='append')
      else
         open(lUnit,file=trim(adjustl(lFile)), status='new')
      endif

      if(fHeader) then
         write(lUnit,1000)
         fHeader = .FALSE.
      endif

      if(tHeader) then
         write(lUnit,"(/2x,'Simulation time: ',g11.5)") TIME
         tHeader = .FALSE.
      endif

      write(lUnit,1001) IJK, I_OF(IJK), J_OF(IJK), K_OF(IJK)
      write(lUnit,"(6x,A,1X,g11.5,$)") 'RO_g:', RO_g(IJK)
      write(lUnit,"(2x,A,1X,g11.5,$)") 'P_g:', P_g(IJK)
      write(lUnit,"(2x,A,1X,g11.5)") 'T_g:', T_g(IJK)
      if(CARTESIAN_GRID) then
         write(lUnit,"(6x,A,1X,L1,$)") 'Cut Cell:', CUT_CELL_AT(IJK)
         write(lUnit,"(6x,A,1X,L1)") 'Small Cell:', SMALL_CELL_AT(IJK)
         write(lUnit,"(6x,'Coordinates (E/N/T): ',1X,3(2x, g17.8))") &
            xg_e(I_OF(IJK)), yg_n(J_of(ijk)), zg_t(k_of(ijk))
      endif

      close(lUnit)

      RETURN

 1000 FORMAT(2X,'One or more cells have reported a negative gas',      &
         ' density (RO_g(IJK)). If',/2x,'this is a persistent issue,', &
         ' lower UR_FAC(1) in mfix.dat.')

 1001 FORMAT(/4X,'IJK: ',I8,7X,'I: ',I4,'  J: ',I4,'  K: ',I4)

      END SUBROUTINE ROgErr_LOG

      END SUBROUTINE PHYSICAL_PROP
