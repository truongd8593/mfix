# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 10:16:52 2013

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

    FORTRAN_CMD  = "ifort"
    OBJ_EXT      = "o"
    FORTRAN_EXT  = "f"
    
    libs_d   = " ${DPO}odepack.a ${DPO}blas90.a ${DPO}dgtsv90.a " 
    libs     = "${DPO}odepack.a ${DPO}blas90.a   ${DPO}dgtsv90.a "
    
    omp      = "-openmp -assume cc_omp"

    compile_common_d = " -convert big_endian -assume byterecl -O0 -fpe0  -CB -traceback" + " "    
    compile_d1       = "-g -c -I. -FR" + compile_common_d
    compile_d2       = "-g -c -I. -FR" + compile_common_d 
    compile_d3       = "-g -c -I." + compile_common_d
    
    link_d =" "
    link   = " "
    
    compile_common  = " -c -I. -w -w95 -i_dynamic -ip -convert big_endian -assume byterecl" + " "
    optimize_common = " -O2 "
    
    compile1        = compile_common + " -FR " + optimize_common
    compile2        = compile_common + " -FR " + optimize_common
    compile3        = compile_common + optimize_common

    
########################################################################################
########################################################################################
#******************************USER INPUT END*******************************************


if __name__ == "__main__":
    get_comp_flags()
    print compile_d1
    
    print compile_d2
    
