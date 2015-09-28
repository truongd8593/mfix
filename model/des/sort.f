! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

module sort

  implicit none
  public :: SORT_PARTICLES, FIND_STATE_BOUNDS
  private :: Partition, Swap

contains

  subroutine FIND_STATE_BOUNDS
    USE discretelement
    implicit none
    integer :: state, ii
    state = 1
    do ii=1, MAX_PIP
       if (state .eq. PARTICLE_STATE(ii)) then
          STATE_BOUNDS(2,state) = ii
       else
          state = state + 1
          STATE_BOUNDS(1,state) = ii
          STATE_BOUNDS(2,state) = ii
       endif
    end do
  end subroutine FIND_STATE_BOUNDS

  recursive subroutine SORT_PARTICLES(first,last)
    use discretelement
    implicit none
    integer, intent(in) :: first, last
    integer :: iq

    if(last-first > 0) then
       call Partition(iq, first, last)
       call SORT_PARTICLES(first,iq-1)
       call SORT_PARTICLES(iq,last)
    endif
  end subroutine SORT_PARTICLES

  subroutine Partition(marker, first, last)
    use discretelement
    implicit none
    integer, intent(out) :: marker
    integer, intent(in) :: first, last
    integer :: ii, jj
    real :: pp      ! pivot value
    pp = PARTICLE_STATE(first)
    ii = first - 1
    jj = last  + 1

    do
       jj = jj-1
       do
          if (PARTICLE_STATE(jj) <= pp) exit
          jj = jj-1
       end do
       ii = ii+1
       do
          if (PARTICLE_STATE(ii) >= pp) exit
          ii = ii+1
       end do
       if (ii < jj) then
          ! exchange particles ii and jj
          call swap(ii,jj)
       elseif (ii == jj) then
          marker = ii+1
          return
       else
          marker = ii
          return
       endif
    end do

  end subroutine Partition

  subroutine swap(ii, jj)

    USE des_rxns
    USE des_thermo
    USE mfix_pic
    USE discretelement
    USE particle_filter
    USE run

    implicit none

    integer, intent(in) :: ii, jj

    call real_swap(des_radius)
    call real_swap(RO_Sol)
    call real_swap(PVOL)
    call real_swap(PMASS)
    call real_swap(OMOI)
    call real_swap2(DES_POS_NEW)
    call real_swap2(DES_VEL_NEW)
    call real_swap2(OMEGA_NEW)
    call real_swap2(PPOS)
    call byte_swap(PARTICLE_STATE)
    call integer_swap(iglobal_id)
    call integer_swap2_reverse(pijk)
    call integer_swap(dg_pijk)
    call integer_swap(dg_pijkprv)
    call logical_swap(ighost_updated)
    call real_swap2(FC)
    call real_swap2(TOW)
    call real_swap(F_GP)
    call integer_swap2(WALL_COLLISION_FACET_ID)
    call real_swap3(WALL_COLLISION_PFT)
    call real_swap2(DRAG_FC)

    call integer_swap2(NEIGHBOR_INDEX)
    call integer_swap2(NEIGHBOR_INDEX_OLD)

    IF(PARTICLE_ORIENTATION) call real_swap2(ORIENTATION)

    IF(FILTER_SIZE > 0) THEN
       call integer_swap2(FILTER_CELL)
       call real_swap2(FILTER_WEIGHT)
    ENDIF

    IF(MPPIC) THEN
       call real_swap(DES_STAT_WT)
       call real_swap2_reverse(PS_GRAD)
       call real_swap2(AVGSOLVEL_P)
       call real_swap(EPG_P)
    ENDIF

    IF(USE_COHESION) THEN
       call real_swap(PostCohesive)
    ENDIF

    IF (DO_OLD) THEN
       call real_swap2(DES_POS_OLD)
       call real_swap2(DES_VEL_OLD)
       call real_swap2(DES_ACC_OLD)
       call real_swap2(OMEGA_OLD)
       call real_swap2(ROT_ACC_OLD)
    ENDIF

    IF(ENERGY_EQ)THEN
       call real_swap(DES_T_s_OLD)
       call real_swap(DES_T_s_NEW)
       call real_swap(DES_C_PS)
       call real_swap2_reverse(DES_X_s)
       call real_swap(Q_Source)

       IF (INTG_ADAMS_BASHFORTH) &
            call real_swap(Q_Source0)
    ENDIF

    IF(ANY_SPECIES_EQ)THEN
       call real_swap2_reverse( DES_R_sp )
       call real_swap2_reverse( DES_R_sc )

       IF (INTG_ADAMS_BASHFORTH) THEN
          call real_swap( dMdt_OLD )
          call real_swap2_reverse( dXdt_OLD )
       ENDIF

       call real_swap( Qint )
    ENDIF

    IF(DES_USR_VAR_SIZE > 0) &
         call real_swap2(DES_USR_VAR)

  contains

    SUBROUTINE BYTE_SWAP(AA)
      IMPLICIT NONE

      INTEGER(KIND=1), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: AA
      INTEGER(KIND=1) :: byte_tmp

      byte_tmp = AA(ii)
      AA(ii) = AA(jj)
      AA(jj) = byte_tmp

    END SUBROUTINE BYTE_SWAP

    SUBROUTINE INTEGER_SWAP(AA)
      IMPLICIT NONE

      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: AA
      INTEGER :: temp

      temp = AA(ii)
      AA(ii) = AA(jj)
      AA(jj) = temp

    END SUBROUTINE INTEGER_SWAP

    SUBROUTINE INTEGER_SWAP2_reverse(AA)
      IMPLICIT NONE

      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: AA
      INTEGER temp(size(AA,2))

      temp(:) = AA(ii,:)
      AA(ii,:) = AA(jj,:)
      AA(jj,:) = temp(:)

    END SUBROUTINE INTEGER_SWAP2_reverse

    SUBROUTINE INTEGER_SWAP2(AA)
      IMPLICIT NONE

      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: AA
      INTEGER temp(size(AA,1))

      temp(:) = AA(:,ii)
      AA(:,ii) = AA(:,jj)
      AA(:,jj) = temp(:)

    END SUBROUTINE INTEGER_SWAP2

    SUBROUTINE LOGICAL_SWAP(AA)
      IMPLICIT NONE

      LOGICAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: AA
      LOGICAL :: temp

      temp = AA(ii)
      AA(ii) = AA(jj)
      AA(jj) = temp

    END SUBROUTINE LOGICAL_SWAP

    SUBROUTINE LOGICAL_SWAP2(AA)
      IMPLICIT NONE

      LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: AA
      LOGICAL temp(size(AA,1))

      temp(:) = AA(:,ii)
      AA(:,ii) = AA(:,jj)
      AA(:,jj) = temp(:)

    END SUBROUTINE LOGICAL_SWAP2

    SUBROUTINE REAL_SWAP(AA)
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: AA
      DOUBLE PRECISION :: temp

      temp = AA(ii)
      AA(ii) = AA(jj)
      AA(jj) = temp

    END SUBROUTINE REAL_SWAP

    SUBROUTINE REAL_SWAP2(AA)
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: AA
      DOUBLE PRECISION temp(size(AA,1))

      temp(:) = AA(:,ii)
      AA(:,ii) = AA(:,jj)
      AA(:,jj) = temp(:)

    END SUBROUTINE REAL_SWAP2

    SUBROUTINE REAL_SWAP2_reverse(AA)
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: AA
      DOUBLE PRECISION temp(size(AA,2))

      temp(:) = AA(ii,:)
      AA(ii,:) = AA(jj,:)
      AA(jj,:) = temp(:)

    END SUBROUTINE REAL_SWAP2_REVERSE

    SUBROUTINE REAL_SWAP3(AA)
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: AA
      DOUBLE PRECISION temp(size(AA,1),size(AA,2))

      temp(:,:) = AA(:,:,ii)
      AA(:,:,ii) = AA(:,:,jj)
      AA(:,:,jj) = temp(:,:)

    END SUBROUTINE REAL_SWAP3

    SUBROUTINE LOGICAL_SWAP2_REVERSE(AA)
      IMPLICIT NONE

      LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: AA
      LOGICAL temp(size(AA,2))

      temp(:) = AA(ii,:)
      AA(ii,:) = AA(jj,:)
      AA(jj,:) = temp(:)

    END SUBROUTINE LOGICAL_SWAP2_REVERSE

  end subroutine swap

end module sort
