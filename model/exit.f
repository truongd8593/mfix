      subroutine mfix_exit (myid)
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
      USE funits
      use compar
      use mpi_utility
      implicit none
      integer, optional, intent(in) :: myid
      logical op

      INTEGER :: mylid

      inquire(unit=unit_log,exist=op)
      if(op) close(unit_log)

      inquire(unit=unit_out,exist=op)
      if(op) close(unit_out)

      if (.not. present(myid)) then
         mylid = myPE
      else
         mylid = myid
      endif

      call exitMPI(mylid)

      stop  
      end subroutine mfix_exit 
