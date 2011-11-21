! dummy variables       
      integer fi,fj,fk,fijk,fto_proc,fproc,fface
      double precision fpos  
! statement functions 
! 
      integer procijk          ! IJK value of processor from given i,j,k position 
      integer iofproc          ! i location of processor for given ijk value
      integer jofproc          ! j location of processor for given ijk value
      integer kofproc          ! k location of processor for given ijk value
      integer dg_funijk        ! ijk of desgrid local
      integer dg_funijk_gl     ! ijk of desgrid global 
      integer dg_funijk_proc   ! ijk of desgrid local and for a processor
      integer dg_ijkconv       ! Convert ijk from one proc to another
      integer dg_iof_gl        ! i of global ijk 
      integer dg_jof_gl        ! j of global ijk 
      integer dg_kof_gl        ! k of glbal ijk  
      integer dg_iof_lo        ! i of local ijk
      integer dg_jof_lo        ! j of local ijk 
      integer dg_kof_lo        ! k of local ijk 
      integer dg_funim         ! iminus ijk of  given ijk
      integer dg_funip         ! iplus ijk of given ijk
      integer dg_funjm         ! jminus ijk of  given ijk
      integer dg_funjp         ! jplus ijk of given ijk
      integer dg_funkm         ! kminus ijk of  given ijk
      integer dg_funkp         ! kplus ijk of given ijk
      integer iofpos           ! i of desgrid based on the given x location
      integer jofpos           ! j of desgrid based on the given y location
      integer kofpos           ! k of desgrid based on the given z location

! following statement functions are introduced to get proc number based on 
! there position in I,J,K
      procijk(fi,fj,fk) =fi+fj*nodesi+fk*nodesi*nodesj 
      iofproc(fijk) = mod(mod(fijk,nodesi*nodesj),nodesi)  
      jofproc(fijk) = mod(fijk - iofproc(fijk),nodesi*nodesj)/nodesi
      kofproc(fijk) = (fijk-iofproc(fijk)-jofproc(fijk)*nodesi)/(nodesi*nodesj)

! functions for desgrid indices 
      dg_funijk(fi,fj,fk) = fj+fi*dg_c1_lo+fk*dg_c2_lo+dg_c3_lo
      dg_funijk_gl(fi,fj,fk) = fj+fi*dg_c1_gl+fk*dg_c2_gl+dg_c3_gl 
      dg_funijk_proc(fi,fj,fk,fproc) = fj+fi*dg_c1_all(fproc)+fk*dg_c2_all(fproc)+dg_c3_all(fproc)

! to find surrounding cells 
      dg_funim(fijk) = fijk - dg_c1_lo 
      dg_funip(fijk) = fijk + dg_c1_lo 
      dg_funjm(fijk) = fijk - 1
      dg_funjp(fijk) = fijk + 1
      dg_funkm(fijk) = fijk - dg_c2_lo 
      dg_funkp(fijk) = fijk + dg_c2_lo 

! following are used to find the i,j,k values from global ijk and local ijk
      dg_jof_gl(fijk) = mod(mod(fijk-1,dg_c2_gl),dg_c1_gl)+dg_jmin2
      dg_iof_gl(fijk) = (mod(fijk-dg_jof_gl(fijk)+dg_jmin2-1,dg_c2_gl))/dg_c1_gl + dg_imin2
      dg_kof_gl(fijk) = (fijk-dg_c3_gl-dg_iof_gl(fijk)*dg_c1_gl-dg_jof_gl(fijk))/dg_c2_gl

      dg_jof_lo(fijk) = mod(mod(fijk-1,dg_c2_lo),dg_c1_lo)+dg_jstart2
      dg_iof_lo(fijk) = (mod(fijk-dg_jof_lo(fijk)+dg_jstart2-1,dg_c2_lo))/dg_c1_lo + dg_istart2
      dg_kof_lo(fijk) = (fijk-dg_c3_lo-dg_iof_lo(fijk)*dg_c1_lo-dg_jof_lo(fijk))/dg_c2_lo


! converting ijk from current proc to another  
      dg_ijkconv(fijk,fface,fto_proc) = dg_funijk_proc(dg_iof_lo(fijk)+dg_cycoffset(fface,1), &
                                                           dg_jof_lo(fijk)+dg_cycoffset(fface,2), &
                                                           dg_kof_lo(fijk)+dg_cycoffset(fface,3), fto_proc) 

! location i,j,k from position 
      iofpos(fpos) = floor((fpos-dg_xstart)*dg_dxinv) + dg_istart1   
      jofpos(fpos) = floor((fpos-dg_ystart)*dg_dyinv) + dg_jstart1   
      kofpos(fpos) = floor((fpos-dg_zstart)*dg_dzinv) + dg_kstart1   

