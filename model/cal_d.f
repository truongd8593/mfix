
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE CAL_D calculates shear velocity
!

      SUBROUTINE CAL_D(V_sh)

      USE param 
      USE param1 
      USE parallel 
      USE matrix 
      USE scales 
      USE constant
      USE physprop
      USE fldvar
      USE visc_s
      USE rxns
      USE toleranc 
      USE geometry
      USE indices
      USE is
      USE tau_s
      USE bc
      USE vshear

      DOUBLE PRECISION V_sh,dis
      DOUBLE PRECISION xdist(IMAX2,JMAX2)
      DOUBLE PRECISION xdist3(IMAX2,JMAX2,KMAX2),cnter3(IMAX2,JMAX2,KMAX2)

      INTEGER I1,J1,I,J
      INCLUDE 'b_force1.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s2.inc'
      INCLUDE 'b_force2.inc'


      IF (NO_K) THEN

!calculate distances
      Do  J1= 1, JMAX2
       	xdist(1,J1)=1d0/(ODX(1))
	DO  I1 = 2, IMAX2
	xdist(I1,J1)=1d0/(ODX(I1))+xdist((I1-1),J1)
        END DO
      END DO 

     Do  IJK= 1, IJKMAX2
          I = I_OF(IJK)
          J = J_OF(IJK)	
	 
	dis=xdist(I,J)

!shear velocity alligned u momentum cells

	VSHE(IJK)=V_sh-2d0*(1d0-(dis-(1d0/ODX(I))/xlength))*V_sh

!shear velocity alligned with scalar, v, w momentum cells

	VSH(IJK)=V_sh-2d0*(1d0-(dis-(1d0/ODX(I))/xlength))*V_sh&
	-2d0*(1d0/(2d0*ODX(I)*xlength))*V_sh
	
      END DO	


      ELSE

!calculate distances
      Do K1=1,KMAX2      
      Do  J1= 1, JMAX2
       
	xdist3(1,J1,K1)=1d0/(ODX(1))
	DO  I1 = 2, IMAX2
	xdist3(I1,J1,K1)=1d0/(ODX(I1))+xdist3((I1-1),J1,K1)
        END DO
      END DO 
      END DO


      Do  IJK= 1, IJKMAX2
          I = I_OF(IJK)
          J = J_OF(IJK)	
	  K = K_OF(IJK)
	
	  dis=xdist3(I,J,K)

!shear velocity alligned u momentum cells

 	VSHE(IJK)=V_sh-2d0*(1d0-(dis-(1d0/ODX(I))/xlength))*V_sh

!shear velocity alligned with scalar, v, w momentum cells

	VSH(IJK)=V_sh-2d0*(1d0-(dis-(1d0/ODX(I))/xlength))*V_sh&
	-2d0*(1d0/(2d0*ODX(I)*xlength))*V_sh
	
      END DO	

      END IF
      RETURN
      END SUBROUTINE CAL_D
