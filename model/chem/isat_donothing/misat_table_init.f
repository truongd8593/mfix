!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: MISAT_TABLE_INIT                                       C
!     Purpose: controlling values for ISAT (reference to ISAT manual)     C
!                                                                         C
!     Author: Nan Xie                                   Date: 02-Aug-04   C
!     Reviewer:                                         Date:             C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE MISAT_TABLE_INIT
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE mchem
      IMPLICIT NONE
!
!-----------------------------------------------
!     G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!     L o c a l   V a r i a b l e s
!-----------------------------------------------   
!
!====================================================================
!      CALL isatab(idtab,mode,nx,spec0,nf,nh,nhd,usrfg,iusr,rusr,info, &
!     &     rinfo,spect,ga,ha,stats)
!     1)nx, spec0, nf, nh, nhd, usrfg, iusr, rusr, spect, ga, ha, stats
!       are set in react.f 
!     2)  For simpility, usrs just need to provide the following in this 
!        subroutines
!       idtab, mode, info, rinfo
!==================================================================== 
!     unique identifier of the table (idtab >=0 )
!


      RETURN
      END SUBROUTINE MISAT_TABLE_INIT
