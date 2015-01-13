!**********************************************************************!
!                                                                      !
!**********************************************************************!
      module des_data_unpack

      use desmpi

      interface unpack_dbuf
         module procedure unpack_db0, unpack_db1
      end interface unpack_dbuf

      contains


!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      subroutine unpack_db0(lbuf,idata,pface)
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(inout) :: idata

      idata = drecvbuf(lbuf,pface) 
      lbuf = lbuf + 1

      return
      end subroutine unpack_db0

!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      subroutine unpack_db1(lbuf,idata,pface)
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(inout) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      idata = drecvbuf(lbuf:lbuf+lsize-1,pface) 
      lbuf = lbuf + lsize

      return
      end subroutine unpack_db1

      end module des_data_unpack

