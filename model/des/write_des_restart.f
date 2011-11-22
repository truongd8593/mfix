!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_DES_RESTART                                      C
!  Purpose: Writing DES data for restart
!                                                                      C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 26-Jul-06  C
!  Reviewer:Sreekanth Pannala                         Date: 31-Oct-06  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE WRITE_DES_RESTART

      USE param1      
      USE discretelement
      USE run
      USE des_bc

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER BCV     ! Loop counter for no. of DES_BCMI
      LOGICAL ASSOC   ! If a derived data type is associated

!-----------------------------------------------


      OPEN (UNIT=901,FILE=TRIM(RUN_NAME)//'_DES.RES',FORM='Unformatted',STATUS='unknown')

      REWIND (901)

      WRITE (901) PARTICLES
      WRITE (901) VTP_FINDEX
      WRITE (901) TECPLOT_FINDEX      
      WRITE (901) DTSOLID

      WRITE (901) DES_POS_NEW
      WRITE (901) DES_VEL_NEW
      WRITE (901) OMEGA_NEW

      WRITE (901) DES_RADIUS
      WRITE (901) RO_Sol

      WRITE (901) NEIGHBOURS

! Needed for particle contact history in tangential direction
      WRITE (901) PFT
      WRITE (901) PN
      WRITE (901) PV

! J. Musser DES boundary condition data
      WRITE (901) PIP
      WRITE (901) PEA

! These arrays are allocated only if inlet exists
      IF(DES_MI)THEN
         WRITE (901) DES_MI_TIME
         WRITE (901) MI_FACTOR
         WRITE (901) MI_WINDOW

         DO BCV =1, DES_BCMI

            WRITE (901) PARTICLE_PLCMNT(BCV)

            IF(PARTICLE_PLCMNT(BCV) == 'ORDR')THEN

               ASSOC = ASSOCIATED(MI_ORDER(BCV)%VALUE)
               WRITE (901) ASSOC
               IF(ASSOC)THEN
                  WRITE (901) SIZE(MI_ORDER(BCV)%VALUE)
                  WRITE (901) MI_ORDER(BCV)%VALUE
               ENDIF

               ASSOC = ASSOCIATED(I_OF_MI(BCV)%VALUE)
               WRITE (901) ASSOC
               IF(ASSOC)THEN
                  WRITE (901) SIZE(I_OF_MI(BCV)%VALUE)
                  WRITE (901) I_OF_MI(BCV)%VALUE
               ENDIF

               ASSOC = ASSOCIATED(J_OF_MI(BCV)%VALUE)
               WRITE (901) ASSOC
               IF(ASSOC)THEN
                  WRITE (901) SIZE(J_OF_MI(BCV)%VALUE)
                  WRITE (901) J_OF_MI(BCV)%VALUE
               ENDIF

            ENDIF  
         ENDDO
      ENDIF
 
      CLOSE(901)

      END SUBROUTINE WRITE_DES_RESTART 

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
!                                       
!  module name: des_write_restart        
!  purpose: writing des data for restart  
!                         
!  Author : Pradeep G      
!  Purpose : Reads either single restart file or multiple restart files     
!            (based on bdist_io) flag 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^c
      subroutine des_write_restart

      use param1      
      use compar
      use discretelement
      use run
      use des_bc
      use desmpi 
      use machine 
      use cdist 
      use mpi_utility
      implicit none
!-----------------------------------------------
! local variables
!-----------------------------------------------
      integer lboundnum 
      integer li,lj,lres_unit,lnext_rec,ltot_pip,ltor_dimn,lsize
      integer lproc,lparcnt,lglocnt,lgathercnts(0:numpes-1)
      logical lassoc
      character(30) lfilename

!-----------------------------------------------
! remove *********************
      ltor_dimn = 1 + (dimn-2)*2 
! open the restart file  
      if (bdist_io) then 
         lres_unit = 901 
         write(lfilename,'(A,I4.4,A)') trim(run_name)//'_DES',mype,'.res'
         open (unit=lres_unit,file=lfilename,form='unformatted',status='unknown', &
               access='direct',recl=open_n1)
      else 
         if(mype.eq.pe_io) then 
            lres_unit = 901
            write(lfilename,'(A)') trim(run_name)//'_DES.res'
            open (unit=lres_unit,file=lfilename,form='unformatted',status='unknown', &
                  access='direct',recl=open_n1)
         end if 
      end if 

      if (bdist_io) then 
         write(lres_unit,rec=1) pip,ighost_cnt,vtp_findex,tecplot_findex,dtsolid 
         lnext_rec = 6
         do li = 1,dimn
            call des_writepar(lres_unit,des_pos_new(:,li),pip,lnext_rec)
         end do 
         call des_writepar(lres_unit,iglobal_id,pip,lnext_rec)
         do li = 2,4
            call des_writepar(lres_unit,pea(:,li),pip,lnext_rec)
         end do 
         do li = 1,dimn
            call des_writepar(lres_unit,des_vel_new(:,li),pip,lnext_rec)
         end do 
         do li = 1,ltor_dimn
            call des_writepar(lres_unit,omega_new(:,li),pip,lnext_rec)
         end do 
         call des_writepar(lres_unit,des_radius,pip,lnext_rec)
         call des_writepar(lres_unit,ro_sol,pip,lnext_rec)
         
         IF(MPPIC) call des_writepar(lres_unit,des_stat_wt(:),pip, lnext_rec)

         call des_writepar(lres_unit,neighbours(:,1),pip,lnext_rec)
         call des_writepar(lres_unit,pn(:,1),pip,lnext_rec)
         do li = 2,maxneighbors
            call des_writepar(lres_unit,neighbours(:,li),pip,lnext_rec)
            call des_writepar(lres_unit,pn(:,li),pip,lnext_rec)
            call des_writepar(lres_unit,pv(:,li),pip,lnext_rec)
            do lj = 1,dimn  
               call des_writepar(lres_unit,pft(:,li,lj),pip,lnext_rec)
            end do 
         end do 
! Write details for mass inlet and outlet parameters
         if(des_mi)then
            write (lres_unit,rec=lnext_rec) des_mi_time;lnext_rec=lnext_rec+1
            do lboundnum =1, des_bcmi
               write (lres_unit,rec=lnext_rec) mi_factor(lboundnum),mi_window(lboundnum),particle_plcmnt(lboundnum)
               lnext_rec = lnext_rec + 1 
               if(particle_plcmnt(lboundnum) == 'ORDR')then
                  lassoc = associated(mi_order(lboundnum)%value)
                  write (lres_unit,rec=lnext_rec)lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     lsize =size(mi_order(lboundnum)%value)
                     write (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     call out_bin_512i(lres_unit,mi_order(lboundnum)%value,lsize,lnext_rec)
                  endif
                  lassoc = associated(i_of_mi(lboundnum)%value)
                  write (lres_unit,rec=lnext_rec) lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     lsize =size(i_of_mi(lboundnum)%value)
                     write (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     call out_bin_512i(lres_unit,i_of_mi(lboundnum)%value,lsize,lnext_rec)
                  endif
                  lassoc = associated(j_of_mi(lboundnum)%value)
                  write (lres_unit,rec=lnext_rec) lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     lsize =size(j_of_mi(lboundnum)%value)
                     write (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     call out_bin_512i(lres_unit,j_of_mi(lboundnum)%value,lsize,lnext_rec)
                  endif
               endif  
            enddo
         endif

       else 
! for serial IO set the parameters for gather  
         lglocnt = 10
         lparcnt = pip - ighost_cnt
         call global_sum(lparcnt,lglocnt) 
         allocate (dprocbuf(lparcnt),drootbuf(lglocnt),iprocbuf(lparcnt),irootbuf(lglocnt))
         igath_sendcnt = lparcnt 
         lgathercnts = 0
         lgathercnts(mype) = lparcnt
         call global_sum(lgathercnts,igathercnts)
         idispls(0) = 0 
         do lproc = 1,numpes-1 
            idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)  
         end do 

         if(mype.eq.pe_io) write(lres_unit,rec=1) lglocnt,vtp_findex,tecplot_findex,dtsolid 
         lnext_rec = 6
         do li = 1,dimn
            call des_gatherwrite(lres_unit,des_pos_new(:,li),lglocnt,lnext_rec)
         end do 
         call des_gatherwrite(lres_unit,iglobal_id,lglocnt,lnext_rec)
         do li = 2,4
            call des_gatherwrite(lres_unit,pea(:,li),lglocnt,lnext_rec)
         end do 
         do li = 1,dimn
            call des_gatherwrite(lres_unit,des_vel_new(:,li),lglocnt,lnext_rec)
         end do 
         do li = 1,ltor_dimn
            call des_gatherwrite(lres_unit,omega_new(:,li),lglocnt,lnext_rec)
         end do 
         call des_gatherwrite(lres_unit,des_radius,lglocnt,lnext_rec)
         call des_gatherwrite(lres_unit,ro_sol,lglocnt,lnext_rec)

         IF(MPPIC) call des_gatherwrite(lres_unit,des_stat_wt,lglocnt,lnext_rec)

! in case of neighbours use global id instead of local id  
         call des_gatherwrite(lres_unit,neighbours(:,1),lglocnt,lnext_rec)
         call des_gatherwrite(lres_unit,pn(:,1),lglocnt,lnext_rec)
         do li = 2,maxneighbors
            call des_gatherwrite(lres_unit,neighbours(:,li),lglocnt,lnext_rec,ploc2glb=.true.)
            call des_gatherwrite(lres_unit,pn(:,li),lglocnt,lnext_rec,ploc2glb=.true.)
            call des_gatherwrite(lres_unit,pv(:,li),lglocnt,lnext_rec)
            do lj = 1,dimn  
               call des_gatherwrite(lres_unit,pft(:,li,lj),lglocnt,lnext_rec)
            end do 
         end do 
         deallocate (dprocbuf,drootbuf,iprocbuf,irootbuf)
!set neighbour indices from global number to local number

! Write details for mass inlet and outlet parameters
         if(des_mi .and. mype .eq. pe_io)then
            write (lres_unit,rec=lnext_rec) des_mi_time;lnext_rec=lnext_rec+1
            do lboundnum =1, des_bcmi
               write (lres_unit,rec=lnext_rec) mi_factor(lboundnum),mi_window(lboundnum),particle_plcmnt(lboundnum)
               lnext_rec = lnext_rec + 1 
               if(particle_plcmnt(lboundnum) == 'ORDR')then
                  lassoc = associated(mi_order(lboundnum)%value)
                  write (lres_unit,rec=lnext_rec)lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     lsize =size(mi_order(lboundnum)%value)
                     write (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     call out_bin_512i(lres_unit,mi_order(lboundnum)%value,lsize,lnext_rec)
                  endif
                  lassoc = associated(i_of_mi(lboundnum)%value)
                  write (lres_unit,rec=lnext_rec) lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     lsize =size(i_of_mi(lboundnum)%value)
                     write (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     call out_bin_512i(lres_unit,i_of_mi(lboundnum)%value,lsize,lnext_rec)
                  endif
                  lassoc = associated(j_of_mi(lboundnum)%value)
                  write (lres_unit,rec=lnext_rec) lassoc;lnext_rec=lnext_rec+1
                  if(lassoc)then
                     lsize =size(j_of_mi(lboundnum)%value)
                     write (lres_unit,rec=lnext_rec) lsize;lnext_rec=lnext_rec+1
                     call out_bin_512i(lres_unit,j_of_mi(lboundnum)%value,lsize,lnext_rec)
                  endif
               endif  
            enddo
         endif
          
      end if 
      if(bdist_io .or.mype .eq. pe_io) close(lres_unit)
       
      end subroutine des_write_restart 

      
