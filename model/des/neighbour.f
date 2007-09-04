!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: NEIGHBOUR                                              C
!  Purpose: DES - Neighbors search; N-Square, Quadtree(2D)/Octree(3D)  C
!           Now also Cell linked search Aug 07                         C
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!  Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!  Comments: Added cell linked-list search                             C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE NEIGHBOUR

      USE param1
      USE discretelement
      IMPLICIT NONE

      INTEGER I, II, LL, CO, NI, TEMP
      

      DO I = 1, PARTICLES
         PPOS(I,:) = DES_POS_NEW(I,:)
         DO II = 1, MAXNEIGHBORS
            NEIGHBOURS(I,II) = -1
         END DO
         NEIGHBOURS(I,1) = 0
      END DO

      IF (DO_NSQUARE) THEN
         CALL NSQUARE
      ELSE
         N2CT = ZERO
      END IF

      IF (DO_QUADTREE) THEN
         CALL QUADTREE
      ELSE
         QUADCT = ZERO
      END IF

      IF (DO_OCTREE) THEN
         CALL OCTREE
      ELSE
         OCTCT = ZERO
      END IF
      IF (DO_GRID_BASED_SEARCH) THEN 
         CALL GRID_BASED_NEIGHBOR_SEARCH
      ENDIF 
      
      RETURN
      END SUBROUTINE NEIGHBOUR

      
      SUBROUTINE GRID_BASED_NEIGHBOR_SEARCH
        
        USE param1
        USE discretelement
        USE geometry
        
      IMPLICIT NONE

      INTEGER I, II, LL, CO, NI, TEMP, JJ, KK , J, K
      INTEGER KM1, KP1, IM1, IP1, JM1, JP1, PNO, NPG, PC(3), IP2, NLIM
      
      DOUBLE PRECISION  DIST(DIMN), DISTMAG, R_LM
      
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
      MAX_RADIUS = DES_RADIUS(1)
      PRINT*, 'IN CELL LINKED LIST SEARCH'
      DO LL = 1, PARTICLES
         
         NEIGHBOURS(LL,1) = NEIGHBOURS(LL,1) + 1
         NLIM  = NEIGHBOURS(LL,1) + 1
         IF(NLIM.GT.MAXNEIGHBORS) THEN 
            WRITE(*,*) 'NLIM =', NLIM,' > MAXNEIGHBORS =', MAXNEIGHBORS, ' FOR PARTICLE ', LL 
            STOP 
         end IF
         
         NEIGHBOURS(LL,NLIM) = LL
         
         PC(:) =  PIJK(LL,1:3)
         
         II = PC(1)
         JJ = PC(2)
         KK = PC(3)
         IP1 = II
         IM1 = II
         JP1 = JJ
         JM1 = JJ
         KM1 = KK
         KP1 = KK
         
         IF((XE(II) - DES_POS_NEW(LL,1)).LT.(DES_RADIUS(LL)+MAX_RADIUS))  IP1 = MIN(IMAX1, II+1)
         IF((DES_POS_NEW(LL,1) - XE(II-1)).LT.(DES_RADIUS(LL)+MAX_RADIUS))  IM1 = MAX(IMIN1, II-1)
         IF((YN(JJ) - DES_POS_NEW(LL,2)).LT.(DES_RADIUS(LL)+MAX_RADIUS))  JP1 = MIN(JMAX1, JJ+1)
         IF((DES_POS_NEW(LL,2) - YN(JJ-1)).LT.(DES_RADIUS(LL)+MAX_RADIUS))  JM1 = MAX(JMIN1, JJ-1)
         IF(DIMN.EQ.3) THEN 
            IF((ZT(KK) - DES_POS_NEW(LL,3)).LT.(DES_RADIUS(LL)+MAX_RADIUS))  KP1 = MIN(KMAX1, KK+1)
            IF((DES_POS_NEW(LL,3) - ZT(KK-1)).LT.(DES_RADIUS(LL)+MAX_RADIUS))  KM1 = MAX(KMIN1, KK-1)
         ENDIF
         
         DO K = KM1, KP1
            DO J = JM1, JP1
               DO I = IM1, IP1
                  
                  If (ASSOCIATED(PIC(I,J,K)%p)) then
                     NPG = SIZE(PIC(I,J,K)%p)
                  Else
                     NPG = 0
                  Endif
                  
                  Do IP2 = 1,NPG
                     PNO = PIC(I,J,K)%p(ip2)
                     if(PNO.GT.LL) then 
                        R_LM = DES_RADIUS(LL) + DES_RADIUS(PNO)
                        R_LM = FACTOR_RLM*R_LM
                        DIST(:) = DES_POS_NEW(PNO,:) - DES_POS_NEW(LL,:)
                        DISTMAG = SQRT(DES_DOTPRDCT(DIST,DIST))
                        IF(R_LM - DISTMAG.gt.SMALL_NUMBER) THEN
                           
                           NEIGHBOURS(LL,1) = NEIGHBOURS(LL,1) + 1
                           NLIM  = NEIGHBOURS(LL,1) + 1
                           IF(NLIM.GT.MAXNEIGHBORS) THEN 
                              WRITE(*,*) 'NLIM =', NLIM,' > MAXNEIGHBORS =', MAXNEIGHBORS, ' FOR PARTICLE ', LL 
                              WRITE(*,*) 'EITHER REDUCE THE R_LM FACTOR OR INCREASE MN IN MFIX.DAT'
                              STOP 
                           end IF
                           NEIGHBOURS(LL,NLIM) = PNO

                           NEIGHBOURS(PNO,1) = NEIGHBOURS(PNO,1) + 1
                           NLIM  = NEIGHBOURS(PNO,1) + 1
                           IF(NLIM.GT.MAXNEIGHBORS) THEN 
                              WRITE(*,*) 'NLIM =', NLIM,' > MAXNEIGHBORS =', MAXNEIGHBORS, ' FOR PARTICLE ', PNO
                              WRITE(*,*) 'EITHER REDUCE THE R_LM FACTOR OR INCREASE MN IN MFIX.DAT'
                              STOP 
                           end IF
                           NEIGHBOURS(PNO,NLIM) = LL
                           
                        end IF !contact condition
                     end if !PNO.GT.LL
                  end Do !IP2
               end DO
            end DO
         end DO
      end DO
    END SUBROUTINE GRID_BASED_NEIGHBOR_SEARCH
    
