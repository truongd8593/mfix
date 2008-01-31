!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:  CELL_NEAR_WALL                                        C
!  Purpose: FIND THE CELLS THAT IN THE NEIGHBORHOOD OF WALLS OR INTERNAL C
!           SURFACES                                                   C
!                                                                      C
!  Author: Tingwen Li                                 Date: 17-Jan-08  C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

!  This subroutine check the cells in the domain and mark them with different label 
!  to distinguish is the cell is in the neighborhood of walls or internal surfaces.
!  It is excuted once only. 
!  variables: 
!  c_near_w(ijk,1:7)
!           (ijk,7)= -2: wall or flow cells  
!                    -1: fluid cells               
!                    0: fluid cell is near the wall, face contact
!                    1: fluid cell is near the wall, edge contact    
!                    2: fluid cell is near the wall, node contact
!                    3: fluid cell is near the internal surface, face contact
!                    4: fluid cell is near the interanl surface, edge contact
!           (ijk,1~6): cell walls for 0,3
!           (ijk,1~3): cell nodes and edge for 1,2,4 
! Attention!!!
! For edge and node contact, only 1 contact is considered  
! Here cells with pflow and mflow are treated as wall to prevent particles from exitting the domain
! If you want particles to exit the domain freely, disable flow_at(i,j,k) check 
! Any question, please contact litw@interchange.ubc.ca

      SUBROUTINE CELL_NEAR_WALL

      USE discretelement
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE run
      USE geometry
      USE matrix
      USE indices
      USE physprop
      USE drag
      USE constant
      USE compar
      USE sendrecv

      USE toleranc         
      USE is               
      USE bc               
      
      IMPLICIT NONE
      
      integer:: i,j,k,ijk,edgenum, nodenum

      INCLUDE 'function.inc'

      do ijk=IJKSTART3, IJKEND3
         c_near_w(ijk,7)=-2     ! init
         c_near_w(ijk,1:6)=0    ! init
         edgenum=0              ! init
         nodenum=0              ! init
!------------------------------------
         if(fluid_at(ijk))then
            c_near_w(ijk,7)=-1  ! fluid cell   
            I =I_OF(ijk)
            J =J_OF(ijk)
            K =K_OF(ijk) 
            
!     check if any face contact to wall or pflow, mflow bcs
            if(wall_at(funijk(i-1,j,k)) .or. flow_at(funijk(i-1,j,k)))then
               c_near_w(ijk,7)=0
               c_near_w(ijk,1)=1
            end if
            if(wall_at(funijk(i+1,j,k)) .or. flow_at(funijk(i+1,j,k)))then
               c_near_w(ijk,7)=0
               c_near_w(ijk,2)=1
            end if
            if(wall_at(funijk(i,j-1,k)) .or. flow_at(funijk(i,j-1,k)))then
               c_near_w(ijk,7)=0
               c_near_w(ijk,3)=1
            end if
            if(wall_at(funijk(i,j+1,k)) .or. flow_at(funijk(i,j+1,k)))then
               c_near_w(ijk,7)=0
               c_near_w(ijk,4)=1
            end if
            IF(DIMN.EQ.3) THEN
               if(wall_at(funijk(i,j,k-1)) .or. flow_at(funijk(i,j,k-1)))then
                  c_near_w(ijk,7)=0
                  c_near_w(ijk,5)=1
               end if
               if(wall_at(funijk(i,j,k+1)) .or. flow_at(funijk(i,j,k+1)))then
                  c_near_w(ijk,7)=0
                  c_near_w(ijk,6)=1
               end if
            END IF

!     check if edge contact to wall
!     here the possibility that one fluid cell edge contact with wall cell 
!     and node contact with another wall cell is not considered since one cell
!     is not enough to capture the flow field of continue phase. For three 
!     dimensional case, it can be very complex, have to be double checked to make sure
!     everything is included
            if(c_near_w(ijk,7).eq.-1)then ! edge check
               if(wall_at(funijk(i-1,j-1,k)) .or. flow_at(funijk(i-1,j-1,k)))then
                  c_near_w(ijk,7)=1
                  c_near_w(ijk,1)=-1
                  c_near_w(ijk,2)=-1
                  c_near_w(ijk,3)=0
                  edgenum=edgenum+1
               end if	
               if(wall_at(funijk(i+1,j+1,k)) .or. flow_at(funijk(i+1,j+1,k)))then
                  c_near_w(ijk,7)=1
                  c_near_w(ijk,1)=0
                  c_near_w(ijk,2)=0
                  c_near_w(ijk,3)=0      	   	
                  edgenum=edgenum+1
               end if	
               if(wall_at(funijk(i-1,j+1,k)) .or. flow_at(funijk(i-1,j+1,k)))then
                  c_near_w(ijk,7)=1
                  c_near_w(ijk,1)=-1
                  c_near_w(ijk,2)=0
                  c_near_w(ijk,3)=0      	   	
                  edgenum=edgenum+1
               end if	
               if(wall_at(funijk(i+1,j-1,k)) .or. flow_at(funijk(i+1,j-1,k)))then
                  c_near_w(ijk,7)=1
                  c_near_w(ijk,1)=0
                  c_near_w(ijk,2)=-1
                  c_near_w(ijk,3)=0      	   	
                  edgenum=edgenum+1
               end if	
               IF(DIMN.EQ.3) THEN
                  if(wall_at(funijk(i-1,j,k-1)) .or. flow_at(funijk(i-1,j,k-1)))then
                     c_near_w(ijk,7)=1
                     c_near_w(ijk,1)=-1
                     c_near_w(ijk,2)=0      	   	
                     c_near_w(ijk,3)=-1
                     edgenum=edgenum+1      	   	
                  end if	
                  if(wall_at(funijk(i+1,j,k+1)) .or. flow_at(funijk(i+1,j,k+1)))then
                     c_near_w(ijk,7)=1
                     c_near_w(ijk,1)=0
                     c_near_w(ijk,2)=0      	   	
                     c_near_w(ijk,3)=0
                     edgenum=edgenum+1
                  end if	
                  if(wall_at(funijk(i-1,j,k+1)) .or. flow_at(funijk(i-1,j,k+1)))then
                     c_near_w(ijk,7)=1
                     c_near_w(ijk,1)=-1
                     c_near_w(ijk,2)=0      	   	
                     c_near_w(ijk,3)=0
                     edgenum=edgenum+1
                  end if	
                  if(wall_at(funijk(i+1,j,k-1)) .or. flow_at(funijk(i+1,j,k-1)))then
                     c_near_w(ijk,7)=1
                     c_near_w(ijk,1)=0
                     c_near_w(ijk,2)=0      	   	
                     c_near_w(ijk,3)=-1
                     edgenum=edgenum+1
                  end if
                  if(wall_at(funijk(i,j-1,k-1)) .or. flow_at(funijk(i,j-1,k-1)))then
                     c_near_w(ijk,7)=1
                     c_near_w(ijk,1)=0      	   	
                     c_near_w(ijk,2)=-1
                     c_near_w(ijk,3)=-1
                     edgenum=edgenum+1      	   	
                  end if	
                  if(wall_at(funijk(i,j+1,k+1)) .or. flow_at(funijk(i,j+1,k+1)))then
                     c_near_w(ijk,7)=1
                     c_near_w(ijk,1)=0      	   	      	   	
                     c_near_w(ijk,2)=0
                     c_near_w(ijk,3)=0
                     edgenum=edgenum+1      	   	
                  end if	
                  if(wall_at(funijk(i,j+1,k-1)) .or. flow_at(funijk(i,j+1,k-1)))then
                     c_near_w(ijk,7)=1
                     c_near_w(ijk,1)=0      	   	      	   	
                     c_near_w(ijk,2)=0
                     c_near_w(ijk,3)=-1
                     edgenum=edgenum+1
                  end if	
                  if(wall_at(funijk(i,j-1,k+1)) .or. flow_at(funijk(i,j-1,k+1)))then
                     c_near_w(ijk,7)=1
                     c_near_w(ijk,1)=0     
                     c_near_w(ijk,2)=-1
                     c_near_w(ijk,3)=0
                     edgenum=edgenum+1      	   	
                  end if
               END IF     
               if(edgenum.gt.1)then
                  write(*,*)'Error from cell_near_wall',edgenum, ' edges contact with wall'
                  write(*,*)'Currently, only one edge contacts to wall cell is permitted!'
                  stop
               end if 	   
            end if              ! c_near_w(ijk,7) edge check

!     check if edge contact to wall through node, only for 3d case               	 
            IF(DIMN.EQ.3 .and. c_near_w(ijk,7).eq.-1) THEN ! node check
               if(wall_at(funijk(i-1,j-1,k-1)) .or. flow_at(funijk(i-1,j-1,k-1)))then
                  c_near_w(ijk,7)=2
                  c_near_w(ijk,1)=-1
                  c_near_w(ijk,2)=-1
                  c_near_w(ijk,3)=-1
                  nodenum=nodenum+1      	   	
               end if	
               if(wall_at(funijk(i+1,j-1,k-1)) .or. flow_at(funijk(i+1,j-1,k-1)))then
                  c_near_w(ijk,7)=2
                  c_near_w(ijk,1)=0
                  c_near_w(ijk,2)=-1
                  c_near_w(ijk,3)=-1
                  nodenum=nodenum+1      	   	
               end if	
               if(wall_at(funijk(i-1,j+1,k-1)) .or. flow_at(funijk(i-1,j+1,k-1)))then
                  c_near_w(ijk,7)=2
                  c_near_w(ijk,1)=-1
                  c_near_w(ijk,2)=0
                  c_near_w(ijk,3)=-1
                  nodenum=nodenum+1      	   	
               end if	
               if(wall_at(funijk(i+1,j+1,k-1)) .or. flow_at(funijk(i+1,j+1,k-1)))then
                  c_near_w(ijk,7)=2
                  c_near_w(ijk,1)=0
                  c_near_w(ijk,2)=0
                  c_near_w(ijk,3)=-1
                  nodenum=nodenum+1      	   	
               end if	
               if(wall_at(funijk(i-1,j-1,k+1)) .or. flow_at(funijk(i-1,j-1,k+1)))then
                  c_near_w(ijk,7)=2
                  c_near_w(ijk,1)=-1
                  c_near_w(ijk,2)=-1
                  c_near_w(ijk,3)=0
                  nodenum=nodenum+1      	   	
               end if	
               if(wall_at(funijk(i+1,j-1,k+1)) .or. flow_at(funijk(i+1,j-1,k+1)))then
                  c_near_w(ijk,7)=2
                  c_near_w(ijk,1)=0
                  c_near_w(ijk,2)=-1
                  c_near_w(ijk,3)=0
                  nodenum=nodenum+1      	   	
               end if	
               if(wall_at(funijk(i-1,j+1,k+1)) .or. flow_at(funijk(i-1,j+1,k+1)))then
                  c_near_w(ijk,7)=2
                  c_near_w(ijk,1)=-1
                  c_near_w(ijk,2)=0
                  c_near_w(ijk,3)=0
                  nodenum=nodenum+1      	   	
               end if	
               if(wall_at(funijk(i+1,j+1,k+1)) .or. flow_at(funijk(i+1,j+1,k+1)))then
                  c_near_w(ijk,7)=2
                  c_near_w(ijk,1)=0
                  c_near_w(ijk,2)=0
                  c_near_w(ijk,3)=0
                  nodenum=nodenum+1      	   	
               end if
               if(nodenum.gt.1)then
                  write(*,*)'Error from cell_near_wall',nodenum, ' nodes contact with wall'
                  write(*,*)'Currently, only one node contacts to wall cell is permitted!'
                  stop
               end if 	
            END IF  
            
!     check if any face contact to internal surface       
            if(c_near_w(ijk,7).eq.-1) then ! internal surface check
               if(sip_at_e(funijk(i-1,j,k)))then
                  c_near_w(ijk,7)=3
                  c_near_w(ijk,1)=1
               end if	         
               if(sip_at_e(ijk))then
                  c_near_w(ijk,7)=3
                  c_near_w(ijk,2)=1
               end if	
               if(sip_at_n(funijk(i,j-1,k)))then
                  c_near_w(ijk,7)=3
                  c_near_w(ijk,3)=1    
               end if
               if(sip_at_n(ijk))then
                  c_near_w(ijk,7)=3
                  c_near_w(ijk,4)=1    
               end if
               IF(DIMN.EQ.3) THEN
                  if(sip_at_t(funijk(i,j,k-1)))then
                     c_near_w(ijk,7)=3
                     c_near_w(ijk,5)=1
                  end if
                  if(sip_at_t(ijk))then
                     c_near_w(ijk,7)=3
                     c_near_w(ijk,6)=1
                  end if
               END IF  
            end if              !c_near_w(ijk,7).eq.-1

!     check if any edge contact to internal surface
            if(c_near_w(ijk,7).eq.-1) then 
               if(sip_at_e(funijk(i,j+1,k)) .OR. sip_at_n(funijk(i+1,j,k)))then
                  c_near_w(ijk,7)=4
                  c_near_w(ijk,1)=0
                  c_near_w(ijk,2)=0
                  c_near_w(ijk,3)=0
               end if	
               if(sip_at_e(funijk(i,j-1,k)) .OR. sip_at_n(funijk(i+1,j-1,k)))then
                  c_near_w(ijk,7)=4
                  c_near_w(ijk,1)=0
                  c_near_w(ijk,2)=-1
                  c_near_w(ijk,3)=0
               end if	
               if(sip_at_e(funijk(i-1,j-1,k)) .OR. sip_at_n(funijk(i-1,j-1,k)))then
                  c_near_w(ijk,7)=4
                  c_near_w(ijk,1)=-1
                  c_near_w(ijk,2)=-1
                  c_near_w(ijk,3)=0
               end if	
               if(sip_at_e(funijk(i-1,j+1,k)) .OR. sip_at_n(funijk(i-1,j,k)))then
                  c_near_w(ijk,7)=4
                  c_near_w(ijk,1)=-1
                  c_near_w(ijk,2)=0
                  c_near_w(ijk,3)=0
               end if	

               IF(DIMN.EQ.3) THEN
                  if(sip_at_e(funijk(i,j+1,k)) .OR. sip_at_n(funijk(i,j,k+1)))then
                     c_near_w(ijk,7)=4
                     c_near_w(ijk,1)=0
                     c_near_w(ijk,2)=0
                     c_near_w(ijk,3)=0
                  end if	
                  if(sip_at_e(funijk(i,j-1,k)) .OR. sip_at_n(funijk(i,j-1,k+1)))then
                     c_near_w(ijk,7)=4
                     c_near_w(ijk,1)=0
                     c_near_w(ijk,2)=-1
                     c_near_w(ijk,3)=0
                  end if	
                  if(sip_at_e(funijk(i,j-1,k-1)) .OR. sip_at_n(funijk(i,j-1,k-1)))then
                     c_near_w(ijk,7)=4
                     c_near_w(ijk,1)=0
                     c_near_w(ijk,2)=-1
                     c_near_w(ijk,3)=-1
                  end if	
                  if(sip_at_e(funijk(i,j+1,k-1)) .OR. sip_at_n(funijk(i,j,k-1)))then
                     c_near_w(ijk,7)=4
                     c_near_w(ijk,1)=0
                     c_near_w(ijk,2)=0
                     c_near_w(ijk,3)=-1
                  end if
                  
                  if(sip_at_e(funijk(i,j,k+1)) .OR. sip_at_n(funijk(i+1,j,k)))then
                     c_near_w(ijk,7)=4
                     c_near_w(ijk,1)=0
                     c_near_w(ijk,2)=0
                     c_near_w(ijk,3)=0
                  end if	
                  if(sip_at_e(funijk(i,j,k-1)) .OR. sip_at_n(funijk(i+1,j,k-1)))then
                     c_near_w(ijk,7)=4
                     c_near_w(ijk,1)=0
                     c_near_w(ijk,2)=0
                     c_near_w(ijk,3)=-1
                  end if	
                  if(sip_at_e(funijk(i-1,j,k-1)) .OR. sip_at_n(funijk(i-1,j,k-1)))then
                     c_near_w(ijk,7)=4
                     c_near_w(ijk,1)=-1
                     c_near_w(ijk,2)=0
                     c_near_w(ijk,3)=-1
                  end if	
                  if(sip_at_e(funijk(i-1,j,k+1)) .OR. sip_at_n(funijk(i-1,j,k)))then
                     c_near_w(ijk,7)=4
                     c_near_w(ijk,1)=-1
                     c_near_w(ijk,2)=0
                     c_near_w(ijk,3)=0
                  end if	         
               END IF  
            end if              !c_near_w(ijk,7).eq.-1       
            
         end if                 !fluid_at
!------------------------------------            	
!      debugging
!      if(c_near_w(ijk,7).ge.0)then
!      write(*,104)i,j,c_near_w(ijk,7),c_near_w(ijk,1:4)
!      end if
!104   format(7(i6,1x))      
      end do! ijk
    RETURN
  END SUBROUTINE CELL_NEAR_WALL
