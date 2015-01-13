!**********************************************************************!
!                                                                      !
!**********************************************************************!
      module des_data_pack

      use desmpi

      interface pack_dbuf
         module procedure pack_db0, pack_db1
      end interface pack_dbuf

      contains


!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      subroutine pack_db0(lbuf,idata,pface)
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(in) :: idata

      dsendbuf(lbuf,pface) = idata
      lbuf = lbuf + 1

      return
      end subroutine pack_db0

!----------------------------------------------------------------------!
!                                                                      !
!----------------------------------------------------------------------!
      subroutine pack_db1(lbuf,idata,pface)
      integer, intent(inout) :: lbuf
      integer, intent(in) :: pface
      double precision, intent(in) :: idata(:)

      integer :: lsize

      lsize = size(idata)

      dsendbuf(lbuf:lbuf+lsize-1,pface) = idata
      lbuf = lbuf + lsize

      return
      end subroutine pack_db1

      end module des_data_pack

