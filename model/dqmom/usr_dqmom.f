!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR_DQMOM                                              C
!  Purpose: This routine is called from the time loop and is           C
!           user-definable tosolve teh population equation             C
!                                                                      C
!  Author: rong fan                                   Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE USR_DQMOM
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
      USE param 
      USE param1 
      USE run
      USE physprop
      USE geometry
      USE fldvar
      USE output
      USE indices
      USE rxns
      USE constant
      Use ambm
      USE compar
      USE scalars          
      Use usr
      IMPLICIT NONE
!-----------------------------------------------
!
!  Include files defining common blocks here
!
!

!  Define local variables here
!
       INTEGER IJK
!             phase index
       double precision t1,t2
!             beginning time and stop time
       double precision eps
!             error       
       double precision h1,hmin
       integer nok, nbad,I,J,n,K
       double precision YY,XX
       double precision max1, min1, max_min
!
!  Include files defining statement functions here
        INCLUDE 'ep_s1.inc'
        INCLUDE 'ep_s2.inc'
        INCLUDE 'function.inc'
   IF(CALL_DQMOM) THEN


!
!  Insert user-defined code here
        if (time<= 1E-15) then
          t1 =time
          t2= time
        else
          t1= time-dt 
          t2= time
        end if
        
          eps=1.0E-3
          h1=1.0E-4
          hmin=0
          nok=0
          nbad=0
          n=2*Nscalar
!           write(*,*) 'before ode'
!           write(*,*) t1,t2
!            YY=-0.5*DY(1)              
!           DO J = JMIN1,JMAX1
!            YY=YY+0.5*(DY(J-1)+DY(J))
!          XX=-0.5*DX(1)  
!             DO I = IMIN1,IMAX1
!           XX=XX+0.5*(DX(I-1)+DX(I))
!              IJK = FUNIJK(I,J,1)
!               write(*,'(2(F8.3),2(G20.15))') XX,YY,Scalar(IJK,1),Scalar(IJK,2)
!            END DO
!         END DO
 
          

           DO IJK = ijkstart3, ijkend3 
                 
            IF (FLUID_AT(IJK)) THEN
               DO I=1,Nscalar
                ystart(I)=ROP_s(IJK,I)/RO_S(IJK,I)
               ENDDO
                 

               DO I=Nscalar+1,2*Nscalar
               ystart(I)=Scalar(IJK,I-Nscalar)
               ENDDO
   
               IJK_INDEX=IJK

               max1=ystart(1)
 
               DO K=2,Nscalar 
                  max1=MAX(max1,ystart(K))
               ENDDO

               min1=ystart(1)
   
               DO K=2,Nscalar 
                  min1=MIN(min1,ystart(K))
               ENDDO 
              
    
               IF(max1>1.0e-3) THEN 
               call odeint(ystart,n+1,t1,t2,eps,h1,&
               hmin,nok,nbad)               
               ENDIF
               

               DO I=1,Nscalar
               ROP_s(IJK,I) = ystart(I)*RO_S(IJK,I)
               Scalar(IJK,I) = ystart(I+Nscalar)
               ENDDO

   
                      
             end if
            END DO
            endif
 
         
      RETURN  
      END SUBROUTINE USR_DQMOM 
