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
      idtab = 0
!
      mode = 0
!
!     ISAT controlling parameter
!     input array for info()
! -----------------------------
!    set scale if info(1) = 1
      info(1) = 1
!   the table is created from isat_#_P.tab if info(3) = 1
      info(3) = 0
!    check point the table
!    for mpi 1 for node 0 only; 2 for all nodes
      info(4) = 0
!    the number of trees to be used
      info(5) = 16
!    replace if table becomes full if info(10)=1
      info(10) = 1
!    to generate ISATAB performence output on the file isat_#_P.op
!    1 for node 0; 2 for all nodes
      info(11) = 2
!    return statistics in array stats, as well as performing normal 
!    ISATAB operations
      info(12) = 1
!    perform DI if info(30) = 1
      info(30) = 1
! MPI control
!    0 for log file for each processor;1 for processor 0 only
      info(42) = 0

!     
!     input array for rinfo()
!--------------------------------
!     error tolerance
       rinfo(1) = 1.d-2
!     the maximum storage (in megabytes) allowed for the ISATAB table
       rinfo(9) = 1000
!     output increment for .op file
!      rinfo(10) = 1.02
!++++++++++++++++ change++++++++++++++
!     The scale factors for x(50+nx),  must be positive
       rinfo(51) = 1.d-4
       rinfo(52) = 882.d0
       rinfo(53) = 0.15d0
       rinfo(54) = 1.d-6
       rinfo(55) = 0.05d0
       rinfo(56) = 0.05d0
       rinfo(57) = 1.0d0
       rinfo(58) = 1.d0
       rinfo(59) = 882.d0
       rinfo(60) = 0.5d-3
       rinfo(61) = 1.0d0
       rinfo(62) = 1.0d-2
!     The scale factors for fa(50+nx+nf),  must be positive
       rinfo(63) = 1.d-4
       rinfo(64) = 882.d0
       rinfo(65) = 0.15d0
       rinfo(66) = 1.d-6
       rinfo(67) = 0.05d0
       rinfo(68) = 0.05d0
       rinfo(69) = 1.0d0
       rinfo(70) = 1.d0
       rinfo(71) = 882.d0
       rinfo(72) = 0.5d-3
       rinfo(73) = 1.0d0
       rinfo(74) = 1.0d-2
!
!

      RETURN
      END SUBROUTINE MISAT_TABLE_INIT
