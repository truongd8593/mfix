!------------------------------------------------------------------------
! Module           : desmpi_wrapper
! Purpose          : Contains wrapper class for mpi communications- send,recv
!                    and scatter, gather
! Author           : Pradeep.G
!------------------------------------------------------------------------
      module desmpi_wrapper
      use parallel_mpi
      use mpi_utility
      use compar

      interface des_mpi_irecv
      module procedure des_mpi_irecv_db
      end interface
      interface des_mpi_isend
      module procedure des_mpi_isend_db
      end interface
      interface des_mpi_scatterv
      module procedure des_mpi_scatterv_i,des_mpi_scatterv_db
      end interface
      interface des_mpi_gatherv
      module procedure des_mpi_gatherv_i,des_mpi_gatherv_db
      end interface
      contains

!------------------------------------------------------------------------
! Subroutine       : des_mpi_barrier
! Purpose          : Wrapper class for barrier
!
!------------------------------------------------------------------------
      subroutine des_mpi_barrier
      implicit none
! local variables
      character(len=80), parameter :: name = 'des_mpi_barrier'
      return
      end subroutine

!------------------------------------------------------------------------
! Subroutine       : des_mpi_irecv_db
! Purpose          : Wrapper class for mpi_irecv
!------------------------------------------------------------------------
      subroutine des_mpi_irecv_db(precvbuf,precvcnt,ptoproc,ptag,preq,perr)
      implicit none
! dummy variables
      double precision, dimension(:) :: precvbuf
      integer :: precvcnt,ptoproc,ptag,preq,perr
      return
      end subroutine

!------------------------------------------------------------------------
! Subroutine       : des_mpi_isend_db
! Purpose          : Wrapper class for mpi_isend
!------------------------------------------------------------------------
      subroutine des_mpi_isend_db(psendbuf,psendcnt,ptoproc,ptag,preq,perr)
      implicit none
! dummy variables
      double precision, dimension(:) :: psendbuf
      integer :: psendcnt,ptoproc,ptag,preq,perr
      return
      end subroutine

!------------------------------------------------------------------------
! Subroutine       : des_mpi_wait
! Purpose          : Wrapper class for mpi_wait
!------------------------------------------------------------------------
      subroutine des_mpi_wait(preq,perr)
      implicit none
! dummy variables
      integer :: preq,perr
      return
      end subroutine


!------------------------------------------------------------------------
! Subroutine       : des_mpi_scatterv_db
! Purpose          : Wrapper class for mpi_scatterv
!------------------------------------------------------------------------
      subroutine des_mpi_scatterv_db(prootbuf,pscattercnts,pdispls, &
                           pprocbuf,precvcnt,proot,perr )
      implicit none
! dummy variables
      double precision, dimension(:):: prootbuf,pprocbuf
      integer, dimension (:) :: pdispls,pscattercnts
      integer :: precvcnt,proot,perr

      pprocbuf = prootbuf
      return
      end subroutine
!------------------------------------------------------------------------
! Subroutine       : des_mpi_scatterv_i
! Purpose          : Wrapper class for mpi_scatterv
!------------------------------------------------------------------------
      subroutine des_mpi_scatterv_i(prootbuf,pscattercnts,pdispls, &
                           pprocbuf, precvcnt,proot,perr )
      implicit none
! dummy variables
      integer, dimension(:) :: prootbuf,pprocbuf
      integer, dimension (:) :: pdispls,pscattercnts
      integer :: precvcnt,proot,perr

      pprocbuf = prootbuf
      return
      end subroutine
!------------------------------------------------------------------------
! Subroutine       : des_mpi_gatherv_db
! Purpose          : Wrapper class for mpi_gatherv
!------------------------------------------------------------------------
      subroutine des_mpi_gatherv_db(psendbuf,psendcnt,precvbuf, &
                                    precvcnts, pdispls,proot,perr )
      implicit none
! dummy variables
      double precision, dimension(:) :: psendbuf,precvbuf
      integer, dimension (:) :: pdispls,precvcnts
      integer :: psendcnt,proot,perr
      precvbuf = psendbuf
      return
      end subroutine
!------------------------------------------------------------------------
! Subroutine       : des_mpi_gatherv_i
! Purpose          : Wrapper class for mpi_gatherv
!------------------------------------------------------------------------
      subroutine des_mpi_gatherv_i(psendbuf,psendcnt,precvbuf, &
                                    precvcnts, pdispls,proot,perr )
      implicit none
! dummy variables
      integer, dimension(:) :: psendbuf,precvbuf
      integer, dimension (:) :: pdispls,precvcnts
      integer :: psendcnt,proot,perr

      precvbuf = psendbuf
      return
      end subroutine
!------------------------------------------------------------------------
! Subroutine       : des_mpi_stop
! Purpose          : Wrapper class for mpi_abort
!------------------------------------------------------------------------
      subroutine des_mpi_stop
      implicit none
      stop
      end subroutine

      end module
