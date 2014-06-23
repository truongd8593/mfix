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
FORTRAN_CMD =""
omp = ""


compile_d1 = ""
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
    libs   = "${DPO}odepack.a ${DPO}blas90.a   ${DPO}dgtsv90.a"

    OBJ_EXT          = "o"
    FORTRAN_EXT      = "f"
    FORTRAN_CMD      = "pgi"
    omp              = "-mp"
    compile_common   = "-c -I. -Mnosave -Mfreeform -Mrecursive -Mreentrant -byteswapio" + " "
    
    compile1         =" -O -Mdalign" + compile_common    
    compile2         ="-O1 -Mdalign" + compile_common
    compile3         ="-O1 -Mdalign -c -I. -Mnosave -Mfixed -Mrecursive -Mreentrant -byteswapio"
    link             =" "
   
    compile_common_d = "-c -I. -Mnosave -byteswapio" + " "    
    debug_common     = "-g -Mbounds -Mchkptr -Mchkfpstk -Mchkstk -Ktrap=fp "   
    compile_d1       = "-Mfreeform" + compile_common_d   + debug_common
    compile_d2       = "-Mfreeform" + compile_common_d   + debug_common 
    compile_d3       = compile_common_d   + debug_common 
    link_d           = " "   

########################################################################################
########################################################################################
#******************************USER INPUT END*******************************************

if __name__ == "__main__":
    get_comp_flags()
    print compile_d1
    
    print compile_d2
    
