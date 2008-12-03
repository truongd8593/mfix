!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GRID_BASED_NEIGHBOR_SEARCH                             C
!>  Purpose: Cell linked search 
!                                                                      C
!  Author: Rahul Garg                                 Date: 01-Aug-07  C
!  Reviewer: Sreekanth Pannala                        Date: 24-OCT-08  C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
      
      SUBROUTINE GRID_BASED_NEIGHBOR_SEARCH
      
      USE param1
      USE discretelement
      USE geometry
      
      IMPLICIT NONE
      LOGICAL PER_COND, ALREADY_NEIGHBOURS
      INTEGER I, II, LL, CO, NI, TEMP, JJ, KK , J, K, NEIGH_L
      INTEGER KM1, KP1, IM1, IP1, JM1, JP1, PNO, NPG, PC(3), IP2, NLIM
      
      DOUBLE PRECISION  DIST(DIMN), DISTMAG, R_LM, LX, LY, LZ, XPER_FAC, YPER_FAC, ZPER_FAC,  CORD_PNO(DIMN)
      
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
      PRINT*, 'IN CELL LINKED LIST SEARCH'

      LX = XE(IMAX1) - XE(IMIN2)
      LY = YN(JMAX1) - YN(JMIN2)
      LZ = ZT(KMAX1) - ZT(KMIN2)
      
!WRITE(*,*) ' DOMAIN LENGTH = ', LX, LY, LZ, INTX_PER
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

         IP1 = II+1
         IM1 = II-1
         JP1 = JJ+1
         JM1 = JJ-1
         
         IF(DIMN.EQ.3) THEN 
            KP1 = KK+1
            KM1 = KK-1
         end IF


!     !$         IF((XE(II) - DES_POS_NEW(LL,1)).LT.(DES_RADIUS(LL)+MAX_RADIUS))  IP1 = (II+1)
!     !$         IF((DES_POS_NEW(LL,1) - XE(II-1)).LT.(DES_RADIUS(LL)+MAX_RADIUS))  IM1 = (II-1)
!     !$         IF((YN(JJ) - DES_POS_NEW(LL,2)).LT.(DES_RADIUS(LL)+MAX_RADIUS))  JP1 =  (JJ+1)
!     !$         IF((DES_POS_NEW(LL,2) - YN(JJ-1)).LT.(DES_RADIUS(LL)+MAX_RADIUS))  JM1 = (JJ-1)
!     !$         IF(DIMN.EQ.3) THEN 
!     !$            IF((ZT(KK) - DES_POS_NEW(LL,3)).LT.(DES_RADIUS(LL)+MAX_RADIUS))  KP1 = (KK+1)
!     !$            IF((DES_POS_NEW(LL,3) - ZT(KK-1)).LT.(DES_RADIUS(LL)+MAX_RADIUS))  KM1 = (KK-1)
!     !$         ENDIF
!     !$         IF(LL.EQ.1) PRINT*, 'CHEC',PC(1), PIJK(2,1), IP1, (XE(II) - DES_POS_NEW(LL,1)),(DES_RADIUS(LL)+MAX_RADIUS), DES_RADIUS(LL), MAX_RADIUS
         DO KK = KM1, KP1
            DO JJ = JM1, JP1
               DO II = IM1, IP1
                  
                  I = II
                  J = JJ
                  K = KK
                  
                  XPER_FAC = 0
                  YPER_FAC = 0
                  ZPER_FAC = 0
                  PER_COND = .FALSE.
                  IF(II.GT.IMAX1) THEN 
                     IF(INTX_PER) THEN 
                        I = IMIN1
                        XPER_FAC = one
                        PER_COND = .true.
!WRITE(*,*) 'cond true EAST',I,J,K,SIZE(PIC(I,J,K)%p)
                     ELSE
                        I = IMAX1
                     ENDIF
                  ENDIF
                  
                  IF(II.LT.IMIN1) THEN 
                     IF(INTX_PER) THEN 
                        I = IMAX1
                        XPER_FAC = -one
                        PER_COND = .true.
!WRITE(*,*) 'cond true WEST', I,J,K,SIZE(PIC(I,J,K)%p)
                     ELSE 
                        I = IMIN1
                     ENDIF
                  ENDIF
                  
                  IF(JJ.GT.JMAX1) THEN 
                     IF(INTY_PER) THEN 
                        J = JMIN1
                        YPER_FAC = one
                        PER_COND = .true.
!WRITE(*,*) 'cond true NORTH', I,J,K,SIZE(PIC(I,J,K)%p)
                     ELSE
                        J = JMAX1
                     ENDIF
                  ENDIF
                  
                  IF(JJ.LT.JMIN1) THEN 
                     IF(INTY_PER) THEN 
                        J = JMAX1
                        YPER_FAC = -one
                        PER_COND = .true.
!WRITE(*,*) 'cond true SOUTH', I,J,K,SIZE(PIC(I,J,K)%p)
                     ELSE
                        J = JMIN1
                     ENDIF
                  ENDIF
                  
                  IF(DIMN.EQ.3) THEN 
                     IF(KK.GT.KMAX1) THEN 
                        IF(INTZ_PER) THEN
                           K = KMIN1
                           ZPER_FAC = one
                        ELSE
                           K = KMAX1
                        ENDIF
                     ENDIF
                     
                     IF(KK.LT.KMIN1) THEN 
                        IF(INTZ_PER) THEN
                           K = KMAX1
                           ZPER_FAC = -one
                        ELSE
                           K = KMIN1
                        ENDIF
                     ENDIF
                  ENDIF
                  
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
                        CORD_PNO(1) = DES_POS_NEW(PNO,1) + XPER_FAC*LX
                        CORD_PNO(2) = DES_POS_NEW(PNO,2) + YPER_FAC*LY
                        IF(DIMN.EQ.3) THEN 
                           CORD_PNO(3) = DES_POS_NEW(PNO,3) + ZPER_FAC*LZ
                        ENDIF 
                        
                        
                        DIST(:) = CORD_PNO(:) - DES_POS_NEW(LL,:)
                        DISTMAG = SQRT(DES_DOTPRDCT(DIST,DIST))
                        
                        ALREADY_NEIGHBOURS = .FALSE.
                        
!IF(LL.EQ.1.AND.PNO.EQ.2) THEN 
!   WRITE(*,*)'CORD=', CORD_PNO(1)
!   WRITE(*,*)'DISTMAG', DISTMAG, R_LM-DISTMAG
!ENDIF
                        
                        DO NEIGH_L = 2, NEIGHBOURS(LL,1)+1
                           IF(PNO.EQ. NEIGHBOURS(LL,NEIGH_L)) ALREADY_NEIGHBOURS=.true.
                        ENDDO
                        
                        IF(R_LM - DISTMAG.gt.SMALL_NUMBER.AND.(.NOT.ALREADY_NEIGHBOURS)) THEN 
                           IF(PER_COND) THEN 
!WRITE(*,*) 'pC = ', pc
!WRITE(*,*) 'II, JJ = ', II, JJ
!WRITE(*,*) 'I, J = ', I, J
!WRITE(*,*) 'XYPER_FAC ', XPER_FAC, YPER_FAC
!WRITE(*,*) 'DES_VEL_NEW = ', DES_POS_NEW(PNO,:)
!WRITE(*,*) 'MODIFIED POSITION = ', CORD_PNO(:)
                           ENDIF
                           
                           
                           NEIGHBOURS(LL,1) = NEIGHBOURS(LL,1) + 1
                           NLIM  = NEIGHBOURS(LL,1) + 1
                           IF(NLIM.GT.MAXNEIGHBORS) THEN 
                              WRITE(*,*) 'NLIM =', NLIM,' > MAXNEIGHBORS =', MAXNEIGHBORS, ' FOR PARTICLE LL', LL 
                              WRITE(*,*) 'EITHER REDUCE THE R_LM FACTOR OR INCREASE MN IN MFIX.DAT'
                              PRINT*,'POSL = ',DES_POS_NEW(LL,:)
                              DO NEIGH_L = 2, NEIGHBOURS(LL,1)+1

                                 DIST(:) = DES_POS_NEW(NEIGHBOURS(LL,NEIGH_L),:) - DES_POS_NEW(LL,:)
                                 DISTMAG = SQRT(DES_DOTPRDCT(DIST,DIST))
                                 PRINT*,'LL =',NEIGHBOURS(LL,NEIGH_L), DES_POS_NEW(NEIGHBOURS(LL,NEIGH_L),:)
                                 PRINT*,DISTMAG, FACTOR_RLM*(DES_RADIUS(LL) + DES_RADIUS(NEIGHBOURS(LL,NEIGH_L))), &
                                 DES_RADIUS(LL),  DES_RADIUS(NEIGHBOURS(LL,NEIGH_L)), FACTOR_RLM
                              ENDDO
                              STOP 
                           end IF
                           NEIGHBOURS(LL,NLIM) = PNO

                           NEIGHBOURS(PNO,1) = NEIGHBOURS(PNO,1) + 1
                           NLIM  = NEIGHBOURS(PNO,1) + 1
                           IF(NLIM.GT.MAXNEIGHBORS) THEN 
                              WRITE(*,*) 'NLIM =', NLIM,' > MAXNEIGHBORS =', MAXNEIGHBORS, ' FOR PARTICLE PNO', PNO
                              WRITE(*,*) 'EITHER REDUCE THE R_LM FACTOR OR INCREASE MN IN MFIX.DAT'
                              STOP 
                           end IF
                           NEIGHBOURS(PNO,NLIM) = LL
                           
                        end IF  !contact condition
                     end if     !PNO.GT.LL
                  end Do        !IP2
               end DO
            end DO
         end DO
      end DO
      END SUBROUTINE GRID_BASED_NEIGHBOR_SEARCH
    
