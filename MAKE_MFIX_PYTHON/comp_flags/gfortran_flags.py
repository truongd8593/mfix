# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:47:12 2013

@author: rgarg
"""
import sys 
libs_d = ""
libs  = ""

OBJ_EXT = ""
FORTRAN_EXT = ""
FORTRAN_CMD=""
omp = ""


compile_d1= ""
compile_d2 = ""
compile_d3 = ""

compile1 = ""
compile2 = ""
compile3 = ""

link = ""
link_d = ""


#Set the below appropriately. For pathscale, module names are capitalized.
#so set the below to true there
makefile_upper_case = False
#the above is needed only needed for gfortran because gfortran uses a different 
#option to direct object and module files. 

def get_comp_flags():

    global libs_d, libs
    global OBJ_EXT, FORTRAN_EXT, omp, FORTRAN_CMD
    global compile_d1, compile_d2, compile_d3
    global compile1, compile2, compile3
    global link, link_d
    
########################################################################################
########################################################################################
#*****************************USER INPUT BEGINNING**************************************  

    libs_d = " ${DPO}odepack.a ${DPO}blas90.a ${DPO}dgtsv90.a"
    libs = "${DPO}odepack.a ${DPO}blas90.a   ${DPO}dgtsv90.a"
    
    OBJ_EXT     ="o"
    FORTRAN_EXT ="f"
    omp         ="-mp"
    FORTRAN_CMD ="gfortran"
    compile_d1   = "-c -I. -fconvert='big-endian' -ffree-form -ffree-line-length-0 -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow -I${mpi_include} -I${DPO} -g "
    compile_d2   = "-c -I. -fconvert='big-endian' -ffree-form -ffree-line-length-0 -I${mpi_include} -I${DPO} -g "
    compile_d3   = "-c -I. -fconvert='big-endian' -I${mpi_include} -I${DPO} -g "
    link_d       = "-g"
    
    compile1="-c -I. -fconvert='big-endian' -ffree-form -ffree-line-length-0 -m64 -Ofast -flto -mtune=corei7-avx -march=corei7-avx -masm=intel -funroll-loops -I${mpi_include} -I${DPO} "
    
    compile2    = "$compile1"
    compile3="-c -I. -fconvert='big-endian' -I${mpi_include} -I${DPO} -O3 "
    link        = " "

########################################################################################
########################################################################################
#******************************USER INPUT END*******************************************


if __name__ == "__main__":
    get_comp_flags()
    print compile_d1
    
    print compile_d2
    
