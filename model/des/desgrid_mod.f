!------------------------------------------------------------------------
! Module           : desgrid 
! Purpose          : Defines the desgrid and sets the indices; also sets the
!                    communication between the desgrid for ghost cell exchange 
! Author           : Pradeep G
! Comments         : The variables are named similar to the fluid grid
!                    more descriptions about naming could be found in 
!                    compar_mod under dmp_modules folder
!------------------------------------------------------------------------
      module desgrid
      use param1
      use funits
      use geometry
      use compar
      use discretelement
      use constant
      use desmpi_wrapper

! variables related to global block structure 
      integer  :: dg_imin1,dg_imin2,dg_jmin1,dg_jmin2,dg_kmin1,dg_kmin2, &
                  dg_imax1,dg_imax2,dg_jmax1,dg_jmax2,dg_kmax1,dg_kmax2 

! variables contain information of local indices of all processors
     integer,dimension (:),allocatable ::  &
                  dg_istart1_all,dg_istart2_all,dg_iend1_all,dg_iend2_all, &
                  dg_jstart1_all,dg_jstart2_all,dg_jend1_all,dg_jend2_all, &
                  dg_kstart1_all,dg_kstart2_all,dg_kend1_all,dg_kend2_all, &
                  dg_isize_all,dg_jsize_all,dg_ksize_all

! variables relate to processor specific details
      integer  :: dg_istart1,dg_istart2,dg_jstart1,dg_jstart2,dg_kstart1,dg_kstart2, &
                  dg_iend1,dg_iend2,dg_jend1,dg_jend2,dg_kend1,dg_kend2
      integer  :: dg_ijkstart2,dg_ijkend2,dg_ijksize2

! variables related to processor domain size  
      double precision :: dg_xstart,dg_xend,dg_ystart,dg_yend,dg_zstart,dg_zend, &
                          dg_dxinv,dg_dyinv,dg_dzinv 

! Variables related to cyclic boundary  used in desgrid_functions  
      integer,dimension(:,:),allocatable :: dg_cycoffset,icycoffset

! particle in cell related variable 
      type iap2
         integer :: isize
         integer, dimension(:), pointer:: p
      end type iap2
      type(iap2), dimension(:), allocatable:: dg_pic  
      integer, dimension(:), allocatable :: dg_pijk,dg_pijkprv 

! constants required for functions computing local and global ijkvalues   
      integer dg_c1_gl,dg_c2_gl,dg_c3_gl,dg_c1_lo,dg_c2_lo,dg_c3_lo
      integer, dimension(:), allocatable :: dg_c1_all,dg_c2_all,dg_c3_all     

      contains 

!------------------------------------------------------------------------
! Subroutine       : desgrid_init 
! Purpose          : sets indices for desgrid and defines constants 
!                    required for dg_funijk,dg_funijk_gl
!                    communication between the desgrid for ghost cell exchange 
!------------------------------------------------------------------------
      subroutine desgrid_init()
      implicit none 

! local varibles 
      double precision :: ldomlen 
      double precision :: ltempdx,ltempdy,ltempdz 
      integer :: lijkproc,liproc,ljproc,lkproc,lis,lie,ljs,lje,lks,lke 
      integer :: lijk
      include 'des/desgrid_functions.inc' 

!set indices for all processors 
     allocate (dg_istart1_all(0:numpes-1),dg_istart2_all(0:numpes-1),dg_iend1_all(0:numpes-1),dg_iend2_all(0:numpes-1), &
               dg_jstart1_all(0:numpes-1),dg_jstart2_all(0:numpes-1),dg_jend1_all(0:numpes-1),dg_jend2_all(0:numpes-1), &
               dg_kstart1_all(0:numpes-1),dg_kstart2_all(0:numpes-1),dg_kend1_all(0:numpes-1),dg_kend2_all(0:numpes-1), &
               dg_isize_all(0:nodesi-1),dg_jsize_all(0:nodesj-1),dg_ksize_all(0:nodesk-1))
     dg_istart1_all=0;dg_istart2_all=0;dg_iend1_all=0;dg_iend2_all=0
     dg_jstart1_all=0;dg_jstart2_all=0;dg_jend1_all=0;dg_jend2_all=0
     dg_kstart1_all=0;dg_kstart2_all=0;dg_kend1_all=0;dg_kend2_all=0
     dg_isize_all=0;dg_jsize_all=0;dg_ksize_all=0

       
!set grid size based on user input desgridsearch_<ijk>max 
      ltempdx = xlength/desgridsearch_imax 
      ltempdy = ylength/desgridsearch_jmax 
      if (.not. no_k) then
         ltempdz = zlength/desgridsearch_kmax 
      end if  
      dg_ksize_all(:) = 1

      lijkproc = 0 
      do lkproc=0,nodesk-1
         do ljproc=0,nodesj-1
            do liproc=0,nodesi-1
               dg_isize_all(liproc) = (xe(iend1_all(lijkproc))-xe(istart1_all(lijkproc)-1))/ltempdx
               dg_jsize_all(ljproc) = (yn(jend1_all(lijkproc))-yn(jstart1_all(lijkproc)-1))/ltempdy
               if(.not. no_k) dg_ksize_all(lkproc) = (zt(kend1_all(lijkproc))-zt(kstart1_all(lijkproc)-1))/ltempdz
               dg_istart1_all(lijkproc) = sum(dg_isize_all(0:liproc-1)) + 2  
               dg_jstart1_all(lijkproc) = sum(dg_jsize_all(0:ljproc-1)) + 2 
               dg_kstart1_all(lijkproc) = sum(dg_ksize_all(0:lkproc-1)) + 2 
               dg_iend1_all(lijkproc) = dg_isize_all(liproc)+dg_istart1_all(lijkproc)-1
               dg_jend1_all(lijkproc) = dg_jsize_all(ljproc)+dg_jstart1_all(lijkproc)-1
               dg_kend1_all(lijkproc) = dg_ksize_all(lkproc)+dg_kstart1_all(lijkproc)-1
               dg_istart2_all(lijkproc) = dg_istart1_all(lijkproc)-1
               dg_jstart2_all(lijkproc) = dg_jstart1_all(lijkproc)-1
               dg_kstart2_all(lijkproc) = dg_kstart1_all(lijkproc)-1
               dg_iend2_all(lijkproc) =  dg_iend1_all(lijkproc)+1
               dg_jend2_all(lijkproc) = dg_jend1_all(lijkproc)+1
               dg_kend2_all(lijkproc) = dg_kend1_all(lijkproc)+1
               lijkproc = lijkproc + 1
            end do 
         end do 
      end do 
      if(no_k) then
         dg_kstart1_all(:) = 1
         dg_kend1_all(:) = 1
         dg_kstart2_all(:) = 1
         dg_kend2_all(:) = 1
      end if 
      
! set indices for global block 
      dg_imin2 = 1
      dg_imin1 = 2
      dg_imax1 = dg_imin1+sum(dg_isize_all(0:nodesi-1))-1
      dg_imax2 = dg_imax1+1 
      dg_jmin2 = 1
      dg_jmin1 = 2
      dg_jmax1 = dg_jmin1+sum(dg_jsize_all(0:nodesj-1))-1
      dg_jmax2 = dg_jmax1+1 
      if (no_k) then 
         dg_kmin2 = 1
         dg_kmin1 = 1
         dg_kmax1 = 1
         dg_kmax2 = 1
      else 
         dg_kmin2 = 1
         dg_kmin1 = 2
         dg_kmax1 = dg_kmin1+sum(dg_ksize_all(0:nodesk-1))-1
         dg_kmax2 = dg_kmax1+1 
      end if 

! set offset indices for periodic boundaries
      lijkproc = mype 
      liproc = iofproc(lijkproc)
      ljproc = jofproc(lijkproc)
      lkproc = kofproc(lijkproc)
      allocate(dg_cycoffset(2*dimn,3),icycoffset(2*dimn,3))
      dg_cycoffset(:,:) = 0; icycoffset(:,:) = 0 
      if (des_periodic_walls_x) then  
         if(liproc.eq.0) then 
            dg_cycoffset(2,1)= (dg_imax2-dg_imin1)
            icycoffset(2,1)= (imax2-imin1)
         end if 
         if(liproc.eq.nodesi-1) then 
            dg_cycoffset(1,1)=-(dg_imax2-dg_imin1)
            icycoffset(1,1)= -(imax2-imin1)
         end if 
      end if  
      if (des_periodic_walls_y) then  
         if(ljproc.eq.0) then 
            dg_cycoffset(4,2)= (dg_jmax2-dg_jmin1)
            icycoffset(4,2)= (jmax2-jmin1)
         end if 
         if(ljproc.eq.nodesj-1) then 
            dg_cycoffset(3,2)=-(dg_jmax2-dg_jmin1)
            icycoffset(3,2)= -(jmax2-jmin1)
         end if 
      end if  
      if (des_periodic_walls_z) then  
         if(lkproc.eq.0) then 
            dg_cycoffset(6,3)=(dg_kmax2-dg_kmin1)
            icycoffset(6,3)= (kmax2-kmin1)
         end if 
         if(lkproc.eq.nodesk-1) then 
            dg_cycoffset(5,3)=-(dg_kmax2-dg_kmin1)
            icycoffset(5,3)= -(kmax2-kmin1)
         end if 
      end if  
        

! set the indices and variables for the current processor
      dg_istart2 = dg_istart2_all(mype)
      dg_istart1 = dg_istart1_all(mype)
      dg_iend1 = dg_iend1_all(mype)
      dg_iend2 = dg_iend2_all(mype)
      dg_jstart2 = dg_jstart2_all(mype)
      dg_jstart1 = dg_jstart1_all(mype)
      dg_jend1 = dg_jend1_all(mype)
      dg_jend2 = dg_jend2_all(mype)
      dg_kstart2 = dg_kstart2_all(mype)
      dg_kstart1 = dg_kstart1_all(mype)
      dg_kend1 = dg_kend1_all(mype)
      dg_kend2 = dg_kend2_all(mype)

! set constants required for dg_funijk dg_funijk_gl for all processors 
      allocate(dg_c1_all(0:numpes-1),dg_c2_all(0:numpes-1),dg_c3_all(0:numpes-1))     

      dg_c1_all=0;dg_c2_all=0;dg_c3_all=0     
      do lijkproc = 0,numpes-1
         dg_c1_all(lijkproc) = (dg_jend2_all(lijkproc)-dg_jstart2_all(lijkproc)+1)
         dg_c2_all(lijkproc) = dg_c1_all(lijkproc)*(dg_iend2_all(lijkproc)-dg_istart2_all(lijkproc)+1)
         dg_c3_all(lijkproc) = -dg_c1_all(lijkproc)*dg_istart2_all(lijkproc) &
                                  -dg_c2_all(lijkproc)*dg_kstart2_all(lijkproc) &
                                  -dg_jstart2_all(lijkproc)+1
      end do 
! set global constants
      dg_c1_gl = (dg_jmax2-dg_jmin2+1)
      dg_c2_gl = dg_c1_gl*(dg_imax2-dg_imin2+1)
      dg_c3_gl = -dg_c1_gl*imin2-dg_c2_gl*kmin2-dg_jmin2+1

!set local constants 
      dg_c1_lo = dg_c1_all(mype)
      dg_c2_lo = dg_c2_all(mype) 
      dg_c3_lo = dg_c3_all(mype)

      dg_ijksize2 = (dg_iend2-dg_istart2+1)* & 
                    (dg_jend2-dg_jstart2+1)* & 
                    (dg_kend2-dg_kstart2+1)  
      dg_ijkstart2 = dg_funijk(dg_istart2,dg_jstart2,dg_kstart2)
      dg_ijkend2 = dg_funijk(dg_iend2,dg_jend2,dg_kend2)

! Confirmation checks 
      if (dg_ijkstart2.ne.1) then 
         if (dmp_log) write(unit_log,'(A)') "Error in dg_funijk: dg_ijkstart2 is not correct"  
         call des_mpi_stop
      end if 
      if (dg_ijkend2.ne.dg_ijksize2) then 
         if (dmp_log) write(unit_log,'(A)') "Error in dg_funijk: dg_ijkend2 is not correct"  
         call des_mpi_stop
      end if 
       

! set the domain length and dx,dy and dz values used in particle_in_cell 
! to pin the particles 
      dg_xstart = xe(istart2)
      dg_xend = xe(iend1)
      dg_ystart = yn(jstart2)
      dg_yend = yn(jend1)
      dg_dxinv = (dg_iend1-dg_istart1+1)/(dg_xend-dg_xstart)
      dg_dyinv = (dg_jend1-dg_jstart1+1)/(dg_yend-dg_ystart)
      if(no_k) then 
         dg_zstart = 1
         dg_zend = 1
         dg_dzinv = 1
      else 
         dg_zstart = zt(kstart2)
         dg_zend = zt(kend1)
         dg_dzinv = (dg_kend1-dg_kstart1+1)/(dg_zend-dg_zstart)
      end if 

! allocate the desgridsearch related variables  
      allocate(dg_pic(dg_ijksize2)); dg_pic(:)%isize = 0
      do lijk = 1,dg_ijksize2
         nullify(dg_pic(lijk)%p)
      end do 
      allocate(dg_pijk(max_pip),dg_pijkprv(max_pip)) ;dg_pijk = 0;dg_pijkprv=0;
!pradeep **************REMOVE
!      call des_dbggrid 
      end subroutine  

!------------------------------------------------------------------------
! Subroutine       : desgrid_check 
! Purpose          : checks the input parameter, if no user inout then 
!                    sets the grid size based on maximum diameter 
!------------------------------------------------------------------------
      subroutine desgrid_check() 

      implicit none 
    
! local variables  
      double precision :: max_diam,tmp_factor,dl_tmp  

      max_diam = 2.0d0*max_radius
      tmp_factor = 3.0d0*(max_diam)
      if (desgridsearch_imax == undefined_i) then
         dl_tmp = xlength/tmp_factor
         desgridsearch_imax = int(dl_tmp)
         if (desgridsearch_imax <= 0) desgridsearch_imax = 1
         if(dmp_log) write(unit_log,'(3x,a,i8)') &
            'desgridsearch_imax was set to ', desgridsearch_imax
      else
         dl_tmp = xlength/dble(desgridsearch_imax)
         if (dl_tmp < max_diam) then
            if(dmp_log) write(unit_log,1037) 'x', 'x', 'i', 'i'                    
            call mfix_exit(mype)
         endif
      endif
      if (desgridsearch_jmax == undefined_i) then
         dl_tmp = ylength/tmp_factor
         desgridsearch_jmax = int(dl_tmp)
         if (desgridsearch_jmax <= 0) desgridsearch_jmax = 1
         if(dmp_log) write(unit_log,'(3x,a,i8)') &
            'desgridsearch_jmax was set to ', desgridsearch_jmax
      else
         dl_tmp = ylength/dble(desgridsearch_jmax)
         if (dl_tmp < max_diam) then
            if(dmp_log) write(unit_log,1037) 'y', 'y', 'j', 'j'
            call mfix_exit(mype)
         endif
      endif
      if (dimn .eq. 2) then
         if (desgridsearch_kmax == undefined_i) then
            desgridsearch_kmax = 1
         elseif(desgridsearch_kmax /= 1) then
            desgridsearch_kmax = 1            
            if(dmp_log) write(unit_log,'(3x,a,i8)') &
               'desgridsearch_kmax was set to ', desgridsearch_kmax
         endif            
      else
         if (desgridsearch_kmax == undefined_i) then
             dl_tmp = zlength/tmp_factor
             desgridsearch_kmax = int(dl_tmp)
             if (desgridsearch_kmax <= 0) desgridsearch_kmax = 1
             if(dmp_log) write(unit_log,'(3x,a,i8)') &
               'desgridsearch_kmax was set to ', desgridsearch_kmax
         else
             dl_tmp = zlength/dble(desgridsearch_kmax)
             if (dl_tmp < max_diam) then
                if(dmp_log) write(unit_log,1037) 'z', 'z', 'k', 'k'
                call mfix_exit(mype)
            endif
         endif
      endif   ! end if/else dimn == 2

 1037 FORMAT(/1X,70('*')//' From: DESGRID_CHECK',/' Message: ',&
          'The neighbor search grid is too fine in the ',A, &
          '-direction',/10X,'with a particle diameter > ',A, &
          'length/desgridsearch_',A,'max. This will',/10X,'create ',&
          'problems for the search method and detecting neighbors',/10X,&
          'Decrease desgridsearch_',A,'max in mfix.dat to coarsen ',&
          'grid.',/1X,70('*')/)
      end subroutine  

!------------------------------------------------------------------------
! Subroutine       : desgrid_pic
! Purpose          : it updates the particle in cell information 
!                    and also located the particle
!                    moving across processor boundary
! Parameters       : plocate - Locate the particles
! Comments         : plocate should be set to true only when particle has to 
!                    be located; if it is false only PIC information will be 
!                    updated by this routine 
!------------------------------------------------------------------------
      subroutine desgrid_pic(plocate)

      implicit none 
! dummy variables
      logical plocate
! local variables
      integer, dimension(dg_ijksize2) :: lpic,lindx 
      integer li,lj,lk,lijk,lcurpar,lparcount,lcurpic
      logical, save :: first_pass = .true.
      integer :: lallocstat,lallocerr
      include 'des/desgrid_functions.inc' 

! locate the particles including ghost cells
      lparcount = 1 
      lpic(:) = 0 
      if (plocate) then 
         dg_pijkprv(:)= dg_pijk(:)
         do lcurpar = 1,max_pip 
            if(lparcount.gt.pip) exit 
            if(.not.pea(lcurpar,1)) cycle 
            lparcount = lparcount + 1
            li = min(dg_iend2,max(dg_istart2,iofpos(des_pos_new(lcurpar,1))))
            lj = min(dg_jend2,max(dg_jstart2,jofpos(des_pos_new(lcurpar,2))))
            if(no_k) then 
               lk = 1 
            else 
               lk = min(dg_kend2,max(dg_kstart2,kofpos(des_pos_new(lcurpar,3))))
            end if 
            dg_pijk(lcurpar) = dg_funijk(li,lj,lk)
            lijk = dg_pijk(lcurpar)
            lpic(lijk) = lpic(lijk) + 1
         end do 
      else 
         do lcurpar = 1,max_pip 
            if(lparcount.gt.pip) exit 
            if(.not.pea(lcurpar,1)) cycle 
            lparcount = lparcount + 1
            lijk = dg_pijk(lcurpar)
            lpic(lijk) = lpic(lijk) + 1
         end do 
      end if 

      if (first_pass) then
         dg_pijkprv(:) = dg_pijk(:)
         first_pass = .false.
      end if 
      
! redfine the array of dg_pic 
      do lijk = dg_ijkstart2,dg_ijkend2
         lcurpic = lpic(lijk)
         if(dg_pic(lijk)%isize.ne.lcurpic) then 
            if(dg_pic(lijk)%isize.gt.0) deallocate(dg_pic(lijk)%p)
            if(lcurpic.gt.0) allocate(dg_pic(lijk)%p(lcurpic))
            dg_pic(lijk)%isize = lcurpic 
         end if 
      end do 

! assign the particle info in pic array 
      lindx(:) = 1
      lparcount = 1
      do lcurpar = 1, max_pip
         if(lparcount.gt.pip) exit
         if(.not.pea(lcurpar,1)) cycle
         lparcount = lparcount + 1
         lijk = dg_pijk(lcurpar) 
         dg_pic(lijk)%p(lindx(lijk)) = lcurpar
         lindx(lijk) = lindx(lijk) +  1
      enddo

      end subroutine   

!------------------------------------------------------------------------
! subroutine       : desgrid_neigh_build ()  
! Purpose          : This particles build the neigh list for the particles 
!                    currently active in the system  
!------------------------------------------------------------------------
      subroutine desgrid_neigh_build () 

      implicit none 
      include 'des/desgrid_functions.inc'
 
! local variables 
      integer lcurpar,lparcount,lmaxneigh,lkoffset
      integer lijk,lic,ljc,lkc,li,lj,lk,ltotpic,lpicloc,lneigh,lneighcnt
      double precision lsearch_rad,ldist,ldistvec(dimn)

! loop through neighbours and build the contact particles list for particles 
! present in the system
      lkoffset = dimn-2 
      lparcount=1 
      do lcurpar =1,max_pip
         if (lparcount.gt.pip) exit
         if (.not. pea(lcurpar,1)) cycle 
         lparcount = lparcount +1 
         if (pea(lcurpar,4)) cycle 
         lneighcnt = 0
         lijk = dg_pijk(lcurpar)
         lic = dg_iof_lo(lijk)
         ljc = dg_jof_lo(lijk)
         lkc = dg_kof_lo(lijk)
         do lk = lkc-lkoffset,lkc+lkoffset 
         do lj = ljc-1,ljc+1 
         do li = lic-1,lic+1 
            lijk = dg_funijk(li,lj,lk) 
            ltotpic =dg_pic(lijk)%isize
            do lpicloc = 1,ltotpic  
               lneigh = dg_pic(lijk)%p(lpicloc)
               if (lneigh.eq.lcurpar) cycle
               lsearch_rad = factor_RLM*(des_radius(lcurpar)+des_radius(lneigh))
               ldistvec = des_pos_new(lcurpar,:)-des_pos_new(lneigh,:)
               ldist = sqrt(dot_product(ldistvec,ldistvec))
               if (ldist.gt.lsearch_rad) cycle 
               lneighcnt = lneighcnt + 1
               if(lneighcnt .gt. MN) then 
                  write(*,'("Particle ",I7," has more than allowed maximum (",I3,") neighbours")') iglobal_id(lcurpar),MN
                  call des_mpi_stop
               end if   
               neighbours(lcurpar,1) = lneighcnt 
               neighbours(lcurpar,lneighcnt+1) = lneigh
            end do 
         end do 
         end do 
         end do 
      end do    
      end subroutine 

!------------------------------------------------------------------------
! subroutine       : des_dbggrid 
! Purpose          : For printing the indices used for desgrid  
!------------------------------------------------------------------------
      subroutine des_dbggrid()
 
      implicit none  
! local varaiables 
      integer lproc,liproc,ljproc,lkproc
      character (30) filename      
      include 'des/desgrid_functions.inc'
      write(filename,'("dbg_desgridn",I4.4,".dat")') mype 
      open(44,file=filename)
      do lproc = 0,numpes-1 
         write(44,*) "Information for Proc =", lproc 
         liproc= iofproc(lproc)
         ljproc= jofproc(lproc)
         lkproc= kofproc(lproc)
         write(44,*) "i,j,k location of proc",liproc,ljproc,lkproc
         write(44,*) "i,j,k size of proc", dg_isize_all(liproc),dg_jsize_all(ljproc),dg_ksize_all(lkproc)
         write(44,*) "-------------------------------------------------"
         write(44,*) "indices        start     end"
         write(44,*) "-------------------------------------------------"
         write(44,*) "for i1:      ",dg_istart1_all(lproc),dg_iend1_all(lproc)
         write(44,*) "for i2:      ",dg_istart2_all(lproc),dg_iend2_all(lproc)
         write(44,*) "for j1:      ",dg_jstart1_all(lproc),dg_jend1_all(lproc)
         write(44,*) "for j2:      ",dg_jstart2_all(lproc),dg_jend2_all(lproc)
         write(44,*) "for k1:      ",dg_kstart1_all(lproc),dg_kend1_all(lproc)
         write(44,*) "for k2:      ",dg_kstart2_all(lproc),dg_kend2_all(lproc)
      end do  
      write(44,*) "-------------------------------------------------"
      write(44,*) "Local Start and end"
      write(44,*) "-------------------------------------------------"
      write(44,*) "for i1:      ",dg_istart1,dg_iend1
      write(44,*) "for i2:      ",dg_istart2,dg_iend2
      write(44,*) "for j1:      ",dg_jstart1,dg_jend1
      write(44,*) "for j2:      ",dg_jstart2,dg_jend2
      write(44,*) "for k1:      ",dg_kstart1,dg_kend1
      write(44,*) "for k2:      ",dg_kstart2,dg_kend2
      write(44,*) "-------------------------------------------------"
      write(44,*) "global Start and end"
      write(44,*) "-------------------------------------------------"
      write(44,*) "for i1:      ",dg_imin1,dg_imax1
      write(44,*) "for i2:      ",dg_imin2,dg_imax2
      write(44,*) "for j1:      ",dg_jmin1,dg_jmax1
      write(44,*) "for j2:      ",dg_jmin2,dg_jmax2
      write(44,*) "for k1:      ",dg_kmin1,dg_kmax1
      write(44,*) "for k2:      ",dg_kmin2,dg_kmax2
      close(44) 
      end subroutine 

      end module 
