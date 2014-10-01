!------------------------------------------------------------------------
! Module           : desgrid
! Purpose          : Defines the desgrid and sets the indices; also sets the
!                    communication between the desgrid for ghost cell exchange
! Author           : Pradeep G
! Comments         : The variables are named similar to the fluid grid
!                    more descriptions about naming could be found in
!                    compar_mod under dmp_modules folder
!
! Contains following subroutines:
!    desgrid_init, desgrid_check, desgrid_pic,
!    desgrid_neigh_build, des_dbggrid
!------------------------------------------------------------------------
      module desgrid

!-----------------------------------------------
! Modules
!-----------------------------------------------
      use param1
      use funits
      use geometry
      use compar
      use discretelement
      use constant
      use desmpi_wrapper
!-----------------------------------------------

! variables related to global block structure
!---------------------------------------------------------------------//
      integer  :: dg_imin1, dg_imax1, &
                  dg_imin2, dg_imax2

      integer  :: dg_jmin1, dg_jmax1, &
                  dg_jmin2, dg_jmax2

      integer  :: dg_kmin1, dg_kmax1, &
                  dg_kmin2, dg_kmax2

! variables contain information of local indices of all processors
!---------------------------------------------------------------------//
      integer,dimension (:),allocatable ::      &
                  dg_istart1_all, dg_iend1_all, &
                  dg_istart2_all, dg_iend2_all, &
                  dg_isize_all

      integer,dimension (:),allocatable ::      &
                  dg_jstart1_all, dg_jend1_all, &
                  dg_jstart2_all, dg_jend2_all, &
                  dg_jsize_all

      integer,dimension (:),allocatable ::      &
                  dg_kstart1_all, dg_kend1_all, &
                  dg_kstart2_all, dg_kend2_all, &
                  dg_ksize_all

! variables relate to processor specific details
!---------------------------------------------------------------------//
      integer  :: dg_istart,  dg_iend,  &
                  dg_istart1, dg_iend1, &
                  dg_istart2, dg_iend2

      integer  :: dg_jstart,  dg_jend,  &
                  dg_jstart1, dg_jend1, &
                  dg_jstart2, dg_jend2

      integer  :: dg_kstart,  dg_kend,  &
                  dg_kstart1, dg_kend1, &
                  dg_kstart2, dg_kend2

      integer  :: dg_ijkstart2, dg_ijkend2, dg_ijksize2


! variables related to processor domain size
!---------------------------------------------------------------------//
      double precision :: dg_xstart, dg_xend, dg_dxinv
      double precision :: dg_ystart, dg_yend, dg_dyinv
      double precision :: dg_zstart, dg_zend, dg_dzinv

! Variables related to cyclic boundary  used in desgrid_functions
      integer,dimension(:,:),allocatable :: dg_cycoffset, icycoffset

! particle in cell related variable
      type iap2
         integer :: isize
         integer, dimension(:), pointer:: p
      end type iap2
      type(iap2), dimension(:), allocatable:: dg_pic
      integer, dimension(:), allocatable :: dg_pijk,dg_pijkprv

! constants required for functions computing local and global ijkvalues
      integer dg_c1_gl, dg_c2_gl, dg_c3_gl  ! global
      integer dg_c1_lo, dg_c2_lo, dg_c3_lo  ! local

      integer, dimension(:), allocatable :: dg_c1_all,dg_c2_all,dg_c3_all


      contains

! following functions are to get proc number from I,J,K position
        integer function procijk(fi,fj,fk)
          implicit none
          integer fi,fj,fk
          procijk =fi+fj*nodesi+fk*nodesi*nodesj
        end function procijk

        integer function iofproc(fijk)
          implicit none
          integer fijk
          iofproc = mod(mod(fijk,nodesi*nodesj),nodesi)
        end function iofproc

        integer function jofproc(fijk)
          implicit none
          integer fijk
          jofproc = mod(fijk - iofproc(fijk),nodesi*nodesj)/nodesi
        end function jofproc

        integer function kofproc(fijk)
          implicit none
          integer fijk
          kofproc = (fijk-iofproc(fijk)-jofproc(fijk)*nodesi)/(nodesi*nodesj)
        end function kofproc

! functions for desgrid indices
        integer function dg_funijk(fi,fj,fk)
          implicit none
          integer fi,fj,fk
          dg_funijk = fj+fi*dg_c1_lo+fk*dg_c2_lo+dg_c3_lo
        end function dg_funijk

        integer function dg_funijk_gl(fi,fj,fk)
          implicit none
          integer fi,fj,fk
          dg_funijk_gl = fj+fi*dg_c1_gl+fk*dg_c2_gl+dg_c3_gl
        end function dg_funijk_gl

        integer function dg_funijk_proc(fi,fj,fk,fproc)
          implicit none
          integer fi,fj,fk,fproc
          dg_funijk_proc = fj+fi*dg_c1_all(fproc)+fk*dg_c2_all(fproc)+dg_c3_all(fproc)
        end function dg_funijk_proc

! to find surrounding cells
        integer function dg_funim(fijk)
          implicit none
          integer fijk
          dg_funim = fijk - dg_c1_lo
        end function dg_funim

        integer function dg_funip(fijk)
          implicit none
          integer fijk
          dg_funip = fijk + dg_c1_lo
        end function dg_funip

        integer function dg_funjm(fijk)
          implicit none
          integer fijk
          dg_funjm = fijk - 1
        end function dg_funjm

        integer function dg_funjp(fijk)
          implicit none
          integer fijk
          dg_funjp = fijk + 1
        end function dg_funjp

        integer function dg_funkm(fijk)
          implicit none
          integer fijk
          dg_funkm = fijk - dg_c2_lo
        end function dg_funkm

        integer function dg_funkp(fijk)
          implicit none
          integer fijk
          dg_funkp = fijk + dg_c2_lo
        end function dg_funkp

! following are used to find the i,j,k values from global ijk and local ijk
        integer function dg_jof_gl(fijk)
          implicit none
          integer fijk
          dg_jof_gl = mod(mod(fijk-1,dg_c2_gl),dg_c1_gl)+dg_jmin2
        end function dg_jof_gl

        integer function dg_iof_gl(fijk)
          implicit none
          integer fijk
          dg_iof_gl = (mod(fijk-dg_jof_gl(fijk)+dg_jmin2-1,dg_c2_gl))/dg_c1_gl + dg_imin2
        end function dg_iof_gl

        integer function dg_kof_gl(fijk)
          implicit none
          integer fijk
          dg_kof_gl = (fijk-dg_c3_gl-dg_iof_gl(fijk)*dg_c1_gl-dg_jof_gl(fijk))/dg_c2_gl
        end function dg_kof_gl

        integer function dg_jof_lo(fijk)
          implicit none
          integer fijk
          dg_jof_lo = mod(mod(fijk-1,dg_c2_lo),dg_c1_lo)+dg_jstart2
        end function dg_jof_lo

        integer function dg_iof_lo(fijk)
          implicit none
          integer fijk
          dg_iof_lo = (mod(fijk-dg_jof_lo(fijk)+dg_jstart2-1,dg_c2_lo))/dg_c1_lo + dg_istart2
        end function dg_iof_lo

        integer function dg_kof_lo(fijk)
          implicit none
          integer fijk
          dg_kof_lo = (fijk-dg_c3_lo-dg_iof_lo(fijk)*dg_c1_lo-dg_jof_lo(fijk))/dg_c2_lo
        end function dg_kof_lo

! converting ijk from current proc to another
        integer function dg_ijkconv(fijk,fface,fto_proc)
          implicit none
          integer fijk,fto_proc,fface
          dg_ijkconv = dg_funijk_proc(dg_iof_lo(fijk)+dg_cycoffset(fface,1), &
               dg_jof_lo(fijk)+dg_cycoffset(fface,2), &
               dg_kof_lo(fijk)+dg_cycoffset(fface,3), fto_proc)
        end function dg_ijkconv

! location i,j,k from position
        integer function iofpos(fpos)
          implicit none
          double precision fpos
          iofpos = floor((fpos-dg_xstart)*dg_dxinv) + dg_istart1
        end function iofpos

        integer function jofpos(fpos)
          implicit none
          double precision fpos
          jofpos = floor((fpos-dg_ystart)*dg_dyinv) + dg_jstart1
        end function jofpos

        integer function kofpos(fpos)
          implicit none
          double precision fpos
          kofpos = floor((fpos-dg_zstart)*dg_dzinv) + dg_kstart1
        end function kofpos

!------------------------------------------------------------------------
! Subroutine       : desgrid_init
! Purpose          : sets indices for desgrid and defines constants
!                    required for dg_funijk,dg_funijk_gl
!                    communication between the desgrid for ghost cell exchange
!------------------------------------------------------------------------
      subroutine desgrid_init()
      implicit none
!-----------------------------------------------
! Local varibles
!-----------------------------------------------
      double precision :: ldomlen
      double precision :: ltempdx,ltempdy,ltempdz
      integer :: lijkproc,liproc,ljproc,lkproc,lis,lie,ljs,lje,lks,lke
      integer :: lijk
!-----------------------------------------------

! set indices for all processors
     allocate (dg_istart1_all(0:numpes-1), dg_iend1_all(0:numpes-1))
     allocate (dg_istart2_all(0:numpes-1), dg_iend2_all(0:numpes-1))
     allocate (dg_isize_all(0:nodesi-1))

     allocate (dg_jstart1_all(0:numpes-1), dg_jend1_all(0:numpes-1))
     allocate (dg_jstart2_all(0:numpes-1), dg_jend2_all(0:numpes-1))
     allocate (dg_jsize_all(0:nodesj-1))

     allocate (dg_kstart1_all(0:numpes-1), dg_kend1_all(0:numpes-1))
     allocate (dg_kstart2_all(0:numpes-1), dg_kend2_all(0:numpes-1))
     allocate (dg_ksize_all(0:nodesk-1))


     dg_istart1_all=0; dg_iend1_all=0
     dg_istart2_all=0; dg_iend2_all=0
     dg_isize_all=0

     dg_jstart1_all=0; dg_jend1_all=0
     dg_jstart2_all=0; dg_jend2_all=0
     dg_jsize_all=0

     dg_kstart1_all=0; dg_kend1_all=0
     dg_kstart2_all=0; dg_kend2_all=0
     dg_ksize_all=0


! set grid size based on user input desgridsearch_<ijk>max
      ltempdx = xlength/desgridsearch_imax
      ltempdy = ylength/desgridsearch_jmax
      if(do_k) ltempdz = zlength/desgridsearch_kmax
      dg_ksize_all(:) = 1

      lijkproc = 0
      do lkproc=0,nodesk-1
         do ljproc=0,nodesj-1
            do liproc=0,nodesi-1
               dg_isize_all(liproc) = (xe(iend1_all(lijkproc))-xe(istart1_all(lijkproc)-1))/ltempdx
               dg_jsize_all(ljproc) = (yn(jend1_all(lijkproc))-yn(jstart1_all(lijkproc)-1))/ltempdy
               if(do_k) dg_ksize_all(lkproc) = (zt(kend1_all(lijkproc))-zt(kstart1_all(lijkproc)-1))/ltempdz
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

!       Defining new set of varaibles to define upper and lower bound of the
! indices to include actual physical boundaries of the problem.
      dg_istart = dg_istart1;   dg_iend = dg_iend1
      dg_jstart = dg_jstart1;   dg_jend = dg_jend1
      dg_kstart = dg_kstart1;   dg_kend = dg_kend1

      if(dg_istart .eq. dg_imin1) dg_istart = dg_imin2
      if(dg_iend   .eq. dg_imax1) dg_iend   = dg_imax2

      if(dg_jstart .eq. dg_jmin1) dg_jstart = dg_jmin2
      if(dg_jend   .eq. dg_jmax1) dg_jend   = dg_jmax2

      if(dg_kstart .eq. dg_kmin1) dg_kstart = dg_kmin2
      if(dg_kend   .eq. dg_kmax1) dg_kend   = dg_kmax2

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

! set local constants
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
         if (dmp_log) write(unit_log,'(A)')&
            "Error in dg_funijk: dg_ijkstart2 is not correct"
         call des_mpi_stop
      end if
      if (dg_ijkend2.ne.dg_ijksize2) then
         if (dmp_log) write(unit_log,'(A)') &
            "Error in dg_funijk: dg_ijkend2 is not correct"
         call des_mpi_stop
      end if


! set the domain length and dx,dy and dz values used in particle_in_cell
! to bin the particles
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
      allocate(dg_pic(dg_ijksize2))
      dg_pic(:)%isize = 0
      do lijk = 1,dg_ijksize2
         nullify(dg_pic(lijk)%p)
      end do
      allocate(dg_pijk(max_pip),dg_pijkprv(max_pip))
      dg_pijk = 0
      dg_pijkprv=0

!      call des_dbggrid
      end subroutine desgrid_init




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
!-----------------------------------------------
! Dummy arguments
!-----------------------------------------------
      logical, INTENT(IN) :: plocate
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer, dimension(dg_ijksize2) :: lpic,lindx
      integer li,lj,lk,lijk,lcurpar,lparcount,lcurpic
      logical, save :: first_pass = .true.
      integer :: lallocstat,lallocerr
!-----------------------------------------------

! locate the particles including ghost cells
      lparcount = 1
      lpic(:) = 0
      if (plocate) then
         dg_pijkprv(:)= dg_pijk(:)
         do lcurpar = 1,max_pip
            if(lparcount.gt.pip) exit
            if(.not.pea(lcurpar,1)) cycle
            lparcount = lparcount + 1
            li = min(dg_iend2,max(dg_istart2,iofpos(des_pos_new(1,lcurpar))))
            lj = min(dg_jend2,max(dg_jstart2,jofpos(des_pos_new(2,lcurpar))))
            if(no_k) then
               lk = 1
            else
               lk = min(dg_kend2,max(dg_kstart2,kofpos(des_pos_new(3,lcurpar))))
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
!$omp parallel do default(shared)                               &
!$omp private(lijk,lcurpic) schedule (guided,50)
      do lijk = dg_ijkstart2,dg_ijkend2
         lcurpic = lpic(lijk)
         if(dg_pic(lijk)%isize.ne.lcurpic) then
            if(dg_pic(lijk)%isize.gt.0) deallocate(dg_pic(lijk)%p)
            if(lcurpic.gt.0) allocate(dg_pic(lijk)%p(lcurpic))
            dg_pic(lijk)%isize = lcurpic
         end if
      end do
!$omp end parallel do


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

!      open (unit=100,file='desgrid.txt',status='unknown')
!      write(100,*)lpic
!      close(100)

      end subroutine desgrid_pic


!------------------------------------------------------------------------
! subroutine       : desgrid_neigh_build ()
! Purpose          : This particles build the neigh list for the particles
!                    currently active in the system
!------------------------------------------------------------------------
      subroutine desgrid_neigh_build ()

!-----------------------------------------------
! Modules
!-----------------------------------------------
      Use des_thermo
      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer lcurpar,lmaxneigh,lkoffset
      integer lijk,lic,ljc,lkc,li,lj,lk,ltotpic,lpicloc,lneigh,lneighcnt
      double precision lsearch_rad,ldistsquared
      double precision :: ldistvec(3)
      double precision :: lcurpar_pos(3)
      double precision :: lcur_off
      integer il_off,iu_off,jl_off,ju_off,kl_off,ku_off
!-----------------------------------------------

! loop through neighbours and build the contact particles list for particles
! present in the system
      lkoffset = dimn-2

      do lcurpar =1,max_pip
         if (.not. pea(lcurpar,1)) cycle
         if (pea(lcurpar,2)) cycle
         if (pea(lcurpar,4)) cycle
         lneighcnt = 0
         lijk = dg_pijk(lcurpar)
         lic = dg_iof_lo(lijk)
         ljc = dg_jof_lo(lijk)
         lkc = dg_kof_lo(lijk)

                 il_off = 1
                 iu_off = 1
                 jl_off = 1
                 ju_off = 1
                 kl_off = 1
                 ku_off = 1

                 lcurpar_pos(:) = des_pos_new(:,lcurpar)
!   The desgrid size should not be less than 2*dia*rlm_factor
                 lcur_off = (lcurpar_pos(1)-dg_xstart)*dg_dxinv - &
                            floor((lcurpar_pos(1)-dg_xstart)*dg_dxinv)
                 if(lcur_off .ge. 0.5) then
                        il_off = 0
                 else
                        iu_off = 0
                 endif

                 lcur_off = (lcurpar_pos(2)-dg_ystart)*dg_dyinv - &
                            floor((lcurpar_pos(2)-dg_ystart)*dg_dyinv)
                 if(lcur_off .ge. 0.5) then
                        jl_off = 0
                 else
                        ju_off = 0
                 endif

                 if(NO_K)then
                        kl_off = 0
                        ku_off = 0
                 else
                         lcur_off = (lcurpar_pos(3)-dg_zstart)*dg_dzinv - &
                                        floor((lcurpar_pos(3)-dg_zstart)*dg_dzinv)
                         if(lcur_off .ge. 0.5) then
                                kl_off = 0
                         else
                                ku_off = 0
                         endif
                 endif

         do lk = lkc-kl_off,lkc+ku_off
         do lj = ljc-jl_off,ljc+ju_off
         do li = lic-il_off,lic+iu_off
            lijk = dg_funijk(li,lj,lk)
            ltotpic =dg_pic(lijk)%isize
            do lpicloc = 1,ltotpic
               lneigh = dg_pic(lijk)%p(lpicloc)
               if (lneigh <= lcurpar) cycle
               lsearch_rad = factor_RLM*(des_radius(lcurpar)+des_radius(lneigh))
               ldistvec = lcurpar_pos(:)-des_pos_new(:,lneigh)
               ldistsquared = dot_product(ldistvec,ldistvec)
               if (ldistsquared.gt.lsearch_rad*lsearch_rad) cycle
               lneighcnt = lneighcnt + 1
               if(lneighcnt .gt. MN) then
                  write(*,'(A,I7,A,I3,A)') 'Particle ', &
                     IGLOBAL_ID(LCURPAR), &
                     'has more than allowed maximum ', MN, &
                     'neighbours'
                  call des_mpi_stop
               end if
               neighbours(lcurpar,1) = lneighcnt
               neighbours(lcurpar,lneighcnt+1) = lneigh
               if (pea(lcurpar,1) .and. .not.pea(lcurpar,4) .and. pea(lneigh,1)) THEN
                  call collision_add(lcurpar, lneigh)
               endif
            end do
         end do
         end do
         end do
      end do

      end subroutine desgrid_neigh_build

!------------------------------------------------------------------------
! subroutine       : des_dbggrid
! Purpose          : For printing the indices used for desgrid
!------------------------------------------------------------------------
      subroutine des_dbggrid()
      implicit none
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      integer lproc,liproc,ljproc,lkproc
      character (30) filename
!-----------------------------------------------

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
      end subroutine des_dbggrid

      end module
