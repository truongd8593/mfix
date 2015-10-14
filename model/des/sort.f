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
    return
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

  recursive subroutine SORT_PARTICLES(first,last,revert)
    use discretelement
    implicit none
    integer, intent(in) :: first, last
    logical, intent(in) :: revert
    integer :: iq

    if(last-first > 0) then
       call Partition(iq, first, last,revert)
       call SORT_PARTICLES(first,iq-1,revert)
       call SORT_PARTICLES(iq,last,revert)
    endif
  end subroutine SORT_PARTICLES

  subroutine Partition(marker, first, last, revert)
    use discretelement
    implicit none
    integer, intent(out) :: marker
    integer, intent(in) :: first, last
    logical, intent(in) :: revert
    integer :: ii, jj
    real :: pp      ! pivot value
    if (revert) then
       pp = orig_index(first)
    else
       pp = PARTICLE_STATE(first)
    end if
    ii = first - 1
    jj = last  + 1

    do
       jj = jj-1
       do
          if (revert) then
             if (orig_index(jj) <= pp) exit
          else
             if (PARTICLE_STATE(jj) <= pp) exit
          end if
          jj = jj-1
       end do
       ii = ii+1
       do
          if (revert) then
             if (orig_index(ii) >= pp) exit
          else
             if (PARTICLE_STATE(ii) >= pp) exit
          endif
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
    use functions
    use geometry

    implicit none

    integer, intent(in) :: ii, jj
    integer :: ier, i, j, k, l, pi, nn

ier = 0
    l = ii
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

         IF(I > IEND1 .OR. I < ISTART1) IER = 1
         IF(J > JEND1 .OR. J < JSTART1) IER = 1
         IF(DO_K .AND. (K > KEND1 .OR. K < KSTART1)) IER = 1

      IF(is_normal(L) .and. IER .ne. 0) THEN
         print *,"III=",istart1,i,iend1
         print *,"JJJ=",jstart1,j,jend1
         print *,"KKK=",kstart1,k,kend1

         print *,"L = ",l,particle_state(L)
         print *,"mype=",mype
         stop 444666
      ENDIF

ier = 0
    l = jj
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

         IF(I > IEND1 .OR. I < ISTART1) IER = 1
         IF(J > JEND1 .OR. J < JSTART1) IER = 1
         IF(DO_K .AND. (K > KEND1 .OR. K < KSTART1)) IER = 1

      IF(is_normal(L) .and. IER .ne. 0) THEN
         print *,"III=",istart1,i,iend1
         print *,"JJJ=",jstart1,j,jend1
         print *,"KKK=",kstart1,k,kend1

         print *,"L = ",l,particle_state(L)
         print *,"mype=",mype
         stop 666444
      ENDIF

    if (PIJK(ii,4)>0) then
    do pi = 1, PINC(PIJK(ii,4))
       if (ii .eq. PIC(PIJK(ii,4))%p(pi)) then
          PIC(PIJK(ii,4))%p(pi) = jj
          !if (ii.eq.1) print *,"PUTTING ",jj," IN ",PIJK(ii,4)
          exit
       endif
    end do
 end if

    if (PIJK(jj,4)>0) then
    do pi = 1, PINC(PIJK(jj,4))
       if (jj .eq. PIC(PIJK(jj,4))%p(pi)) then
          PIC(PIJK(jj,4))%p(pi) = ii
          !if (jj.eq.1) print *,"PUTTING ",ii," IN ",PIJK(jj,4)
          exit
       endif
    end do
 end if

    if (dg_pijk(ii)>0) then
    do pi = 1, dg_pic(dg_pijk(ii))%isize
       if (ii .eq. dg_pic(dg_pijk(ii))%p(pi)) then
          dg_pic(dg_pijk(ii))%p(pi) = jj
          !if (ii.eq.1) print *,"putting ",jj," in ",dg_pijk(ii)
          exit
       endif
    end do
 end if

    if (dg_pijk(jj)>0) then
    do pi = 1, dg_pic(dg_pijk(jj))%isize
       if (jj .eq. dg_pic(dg_pijk(jj))%p(pi)) then
          dg_pic(dg_pijk(jj))%p(pi) = ii
          !if (jj.eq.1) print *,"putting ",ii," in ",dg_pijk(jj)
          exit
       endif
    end do
 end if

    call real_swap(des_radius)
    call real_swap(RO_Sol)
    call real_swap(PVOL)
    call real_swap(PMASS)
    call real_swap(OMOI)
    call real_swap2(DES_POS_NEW)
    call real_swap2(DES_VEL_NEW)
    call real_swap2(OMEGA_NEW)
    call real_swap2(PPOS)
    !call byte_swap(PARTICLE_STATE)
    call integer_swap(PARTICLE_STATE)
    call integer_swap(iglobal_id)
    call integer_swap2_reverse(pijk)
    call integer_swap(dg_pijk)
    call integer_swap(orig_index)
    call integer_swap(dg_pijkprv)
    call logical_swap(ighost_updated)
    call real_swap2(FC)
    call real_swap2(TOW)
    call real_swap(F_GP)
    call integer_swap2(WALL_COLLISION_FACET_ID)
    call real_swap3(WALL_COLLISION_PFT)
    call real_swap2(DRAG_FC)

    do nn=1, neigh_num
       if (neighbors(nn).eq.ii) neighbors(nn) = jj
       if (neighbors(nn).eq.jj) neighbors(nn) = ii
    enddo

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

ier = 0
    l = ii
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

         IF(I > IEND1 .OR. I < ISTART1) IER = 1
         IF(J > JEND1 .OR. J < JSTART1) IER = 1
         IF(DO_K .AND. (K > KEND1 .OR. K < KSTART1)) IER = 1

      IF(is_normal(L) .and. IER .ne. 0) THEN
         print *,"III=",istart1,i,iend1
         print *,"JJJ=",jstart1,j,jend1
         print *,"KKK=",kstart1,k,kend1

         print *,"ii,jj = ",ii,jj,particle_state(ii),particle_state(jj)
         print *,"mype=",mype
         print *,"do_k=",do_k
         stop 262626
      ENDIF

ier = 0
    l = jj
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)

         IF(I > IEND1 .OR. I < ISTART1) IER = 1
         IF(J > JEND1 .OR. J < JSTART1) IER = 1
         IF(DO_K .AND. (K > KEND1 .OR. K < KSTART1)) IER = 1

      IF(is_normal(L) .and. IER .ne. 0) THEN
         print *,"III=",istart1,i,iend1
         print *,"JJJ=",jstart1,j,jend1
         print *,"KKK=",kstart1,k,kend1

         print *,"ii,jj = ",ii,jj,particle_state(ii),particle_state(jj)
         print *,"mype=",mype
         stop 1717171
      ENDIF

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
