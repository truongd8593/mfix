      subroutine mfix_exit (myid)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
      use compar
      use mpi_utility
      implicit none
      integer, optional, intent(in) :: myid

      INTEGER :: mylid

      if (.not. present(myid)) then
         mylid = myPE
      else
         mylid = myid
      endif

      call exitMPI(mylid)

      stop  
      end subroutine mfix_exit 
