!include user defined modules in this file.  Remember to allocate user-defined
!adjustable arrays in usr0.f


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: devol.inc                                              C
!  Purpose: Common block containing data relating to devolatilization  C
!                                                                      C
!  Author: S. Venkatesan                              Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References: MGAS code, Wen, et al. (1982)       C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C


MODULE usr


      Use param
      Use param1
      
      
!                      Type of Coal
      CHARACTER*20     COALNAM
 
!                      other parameters
      DOUBLE PRECISION :: FTC=0.9, FTH=0.08, FTO=0.02, FTN=0.0, FTS=0.0,&
!                       DOCO=0.14, DOCO2=0.13, DOH2O=0.73,&
                       DOCO=0.00, DOCO2=0.0, DOH2O=1.0,&
		       DHH2=0.46, DHCH4=0.54, DHC2H6=0.0, DHC2H4=0.0,&
                       DHC3H8=0.0, DHC6H6=0.0,&
		       COCO=1.0, COCO2=0.0, COH2O=0.0,&
		       CHH2=0.0, CHCH4=1.0, CHC2H6=0.0, CHC2H4=0.0, CHC3H8=0.0,&
		       CHC6H6=0.0,&
                       TOLFC
 
      DOUBLE PRECISION HEATD, HEATC, FVC, FVH, FVO, FVN, FVS,&
                       DAFC, EP_A, f_EP_A, AK2, AE2, AK5, AE5, WG3,&
                       AKM, AEM, AKD, AED, AKC, AEC
!
!                      Total number of gas or solids species;
!                      defined in PHYSICAL_PROP.INC
!
!                      Devolatilization coefficients
      DOUBLE PRECISION BETAD (14), ALPHAD
!
!                      Stoichiometric coefficients in tar combustion rxn
!                      Tar + f3_1 O2 --> f3_3 CO2 + f3_6 H2O
      DOUBLE PRECISION F3_1, F3_3, F3_6
!
!                      Heat of tar combustion reaction
      DOUBLE PRECISION HEATF3
!
!                      Cracking coefficients
      DOUBLE PRECISION BETAC (14), ALPHAC
!
!                      VMSTAR (Gregory and Littlejohn devol. constant)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VMSTAR
!
!                      Sherwood number
      DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: N_sh
!
      
      
END MODULE usr
